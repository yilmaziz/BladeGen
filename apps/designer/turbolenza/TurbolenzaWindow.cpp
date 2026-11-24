// The TURBOlenza window: tree, pages, the solve on a worker thread, post-processing windows,
// saving, and the offscreen self-test (see TurbolenzaWindow.hpp).
#include "TurbolenzaWindow.hpp"

#include "ParameterPanel.hpp"
#include "PlotWidget.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDateTime>
#include <QDir>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QStackedWidget>
#include <QStatusBar>
#include <QTimer>
#include <QTreeWidget>
#include <QVBoxLayout>

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <numbers>

namespace Designer {

using namespace BladeGen;
namespace T = BladeGen::Turbolenza;

namespace {

enum Page : int { kBlade = 0, kCascade, kMesh, kSolver, kRun, kPost };
constexpr double kDeg = 180.0 / std::numbers::pi;
const QStringList kExportNames = {"surface_mach.csv", "surface_pressure.csv", "convergence_history.csv", "fields.csv", "run_summary.json"};

}  // namespace

// ── worker thread -> GUI ─────────────────────────────────────────────────────

struct TurbolenzaWindow::Channel
{
    std::mutex mu;
    std::vector<T::PassRecord> records;  ///< not yet taken by the GUI
    std::vector<std::string> lines;
    std::atomic<bool> cancel{false};
    std::atomic<bool> done{false};
    std::optional<T::Case> c;
    std::optional<T::RunResult> result;
    std::string error;  ///< build_case failed
};

/// The run's observer (v1 monitor.ConsoleProgress + the GUI channel): every line goes to the
/// program's terminal at once (flushed) and to the window's log; observation only.
class TurbolenzaWindow::Observer : public T::RunObserver
{
public:
    explicit Observer(Channel& ch) : ch_(ch) {}
    void say(const std::string& line)
    {
        std::cout << line << std::endl;
        std::lock_guard<std::mutex> g(ch_.mu);
        ch_.lines.push_back(line);
    }
    void start(const T::SolveInfo& i) override
    {
        auto tol = [](double v) {
            char b[32];
            std::snprintf(b, sizeof b, "%.1e", v);
            return std::string(b);
        };
        char b[512];
        std::snprintf(b, sizeof b, "%s solve: %d x %d nodes, at most %d passes; tolerances rho %s, W %s, Kutta %s; mdot' %.6g (%s); "
                      "first beta2 %.4f deg; %d thread(s), %s",
                      i.mode.c_str(), i.n_m, i.n_eta, i.max_passes, tol(i.density_tolerance).c_str(), tol(i.velocity_tolerance).c_str(),
                      i.kutta_tolerance ? tol(*i.kutta_tolerance).c_str() : "-", i.mass_flow_per_depth, i.mass_flow_source.c_str(),
                      i.exit_angle_start * kDeg, i.threads, i.linear_solver.c_str());
        say(b);
    }
    void update(const T::PassRecord& r) override
    {
        const std::string line = T::format_record(r);
        std::cout << line << std::endl;
        std::lock_guard<std::mutex> g(ch_.mu);
        ch_.lines.push_back(line);
        ch_.records.push_back(r);
    }
    void finish(bool converged, const std::string& error) override { say(converged ? std::string("converged") : "stopped: " + error); }

private:
    Channel& ch_;
};

/// A post-processing window (top level); onClose runs when the user closes it.
class TurbolenzaWindow::PostWindow : public QWidget
{
public:
    PostWindow(const QString& title, QWidget* parent) : QWidget(parent, Qt::Window)
    {
        setWindowTitle(title);
        box = new QVBoxLayout(this);
        resize(1000, 680);
    }
    QVBoxLayout* box = nullptr;
    std::function<void()> onClose;

protected:
    void closeEvent(QCloseEvent* e) override
    {
        if (onClose) onClose();
        e->accept();
    }
};

// ── window ───────────────────────────────────────────────────────────────────

TurbolenzaWindow::TurbolenzaWindow(Session& s, QWidget* parent) : QMainWindow(parent, Qt::Window), s_(s)
{
    setWindowTitle("TURBOlenza");
    resize(1500, 950);
    auto* split = new QSplitter(Qt::Horizontal, this);
    tree_ = new QTreeWidget(split);
    tree_->setHeaderHidden(true);
    tree_->setMinimumWidth(190);
    stack_ = new QStackedWidget(split);
    stack_->addWidget(bladePage());
    stack_->addWidget(cascadePage());
    stack_->addWidget(meshPage());
    stack_->addWidget(solverPage());
    stack_->addWidget(runPage());
    stack_->addWidget(postPage());
    split->addWidget(tree_);
    split->addWidget(stack_);
    split->setStretchFactor(1, 1);
    split->setSizes({220, 1280});
    setCentralWidget(split);
    status_ = new QLabel(this);
    statusBar()->addWidget(status_, 1);
    buildTree();
    buildMenus();
    poll_ = new QTimer(this);
    poll_->setInterval(100);
    connect(poll_, &QTimer::timeout, this, [this] { pollRun(); });
    takeLayer(std::min(2, static_cast<int>(s_.design.layers.size()) - 1));
    showPage(kBlade);
}

TurbolenzaWindow::~TurbolenzaWindow()
{
    if (mesh_watch_) mesh_watch_->waitForFinished();
    if (ch_) ch_->cancel = true;
    if (worker_.joinable()) worker_.join();
}

void TurbolenzaWindow::closeEvent(QCloseEvent* e)
{
    if (running()) {
        const auto b = QMessageBox::question(this, "TURBOlenza", "A solve is running. Stop it and close?");
        if (b != QMessageBox::Yes) {
            e->ignore();
            return;
        }
        stopRun();
        worker_.join();  // the solver stops after its current pass
        pollRun();
    }
    for (auto* w : {res_win_, contour_win_, loading_win_})
        if (w) w->hide();
    e->accept();
}

void TurbolenzaWindow::buildTree()
{
    auto* pre = new QTreeWidgetItem(tree_, {"Pre-processing"});
    auto add = [&](QTreeWidgetItem* parent, const QString& name, int page, int which) {
        auto* it = parent ? new QTreeWidgetItem(parent, {name}) : new QTreeWidgetItem(tree_, {name});
        it->setData(0, Qt::UserRole, page);
        it->setData(0, Qt::UserRole + 1, which);
        items_.push_back(it);
        return it;
    };
    pre->setData(0, Qt::UserRole, kBlade);
    add(pre, "Blade", kBlade, 0);
    add(pre, "Cascade", kCascade, 0);
    add(pre, "Mesh", kMesh, 0);
    add(pre, "Solver", kSolver, 0);
    add(nullptr, "Run", kRun, 0);
    auto* post = add(nullptr, "Post-processing", kPost, 0);
    add(post, "Residuals", kPost, 1);
    add(post, "Contours", kPost, 2);
    add(post, "Blade loading", kPost, 3);
    tree_->expandAll();
    connect(tree_, &QTreeWidget::itemClicked, this, [this](QTreeWidgetItem* it, int) {
        showPage(it->data(0, Qt::UserRole).toInt());
        const int which = it->data(0, Qt::UserRole + 1).toInt();
        if (which == 1) showResiduals();
        else if (which == 2) showContours();
        else if (which == 3) showLoading();
    });
}

void TurbolenzaWindow::buildMenus()
{
    auto* file = menuBar()->addMenu("&File");
    file->addAction("&Save results...", this, [this] { saveDialog(); });
    file->addSeparator();
    file->addAction("&Close", this, [this] { close(); });
    auto* run = menuBar()->addMenu("&Run");
    run->addAction("&Generate mesh", this, [this] {
        showPage(kMesh);
        generateMesh();
    });
    run->addAction("&Run", this, [this] { startRun(); });
    run->addAction("&Stop", this, [this] { stopRun(); });
    auto* post = menuBar()->addMenu("&Post-processing");
    post->addAction("&Residuals", this, [this] { showResiduals(); });
    post->addAction("&Contours", this, [this] { showContours(); });
    post->addAction("&Blade loading", this, [this] { showLoading(); });
}

void TurbolenzaWindow::showPage(int index)
{
    stack_->setCurrentIndex(index);
    for (auto* it : items_)
        if (it->data(0, Qt::UserRole).toInt() == index && it->data(0, Qt::UserRole + 1).toInt() == 0) {
            tree_->blockSignals(true);
            tree_->setCurrentItem(it);
            tree_->blockSignals(false);
        }
    if (index == kCascade) refreshCascade();
}

void TurbolenzaWindow::takeLayer(int k)
{
    if (k < 0 || k >= static_cast<int>(s_.design.layers.size())) return;
    source_layer_ = k;
    layer_ = s_.design.layers[static_cast<std::size_t>(k)];  // a copy: the 3-D design is not touched
    syncing_ = true;
    layer_combo_->setCurrentIndex(k);
    syncing_ = false;
    syncBlade();
    setCascadeDefaults(layer_.stagger_deg, layer_.solidity);
    section_.reset();
    section_error_.clear();
    cascade_.reset();
    mesh_case_.reset();
    refreshProfile();
    status_->setText(QString("blade: a copy of layer %1 of the design").arg(QString::fromStdString(layer_.name)));
}

const T::BladeSection* TurbolenzaWindow::section()
{
    if (!section_ && section_error_.isEmpty()) {
        QApplication::setOverrideCursor(Qt::WaitCursor);
        try {
            section_ = T::section_from_layer(layer_, ev_, &fit_);
        } catch (const std::exception& e) {
            section_error_ = e.what();
        }
        QApplication::restoreOverrideCursor();
        refreshProfile();
    }
    return section_ ? &*section_ : nullptr;
}

// ── run ──────────────────────────────────────────────────────────────────────

void TurbolenzaWindow::startRun()
{
    if (running()) return;
    const auto* sec = section();
    if (!sec) {
        run_status_->setText("the section cannot be built: " + section_error_);
        showPage(kRun);
        return;
    }
    run_settings_ = settings();
    records_.clear();
    result_.reset();
    run_case_.reset();
    log_->clear();
    if (res_win_ && res_win_->isVisible()) res_win_->close();  // the plot comes back to the Run page
    ch_ = std::make_unique<Channel>();
    const double tk = run_settings_.solver_mode == T::SolverMode::CoupledKutta ? run_settings_.kutta_tolerance
                                                                               : std::numeric_limits<double>::quiet_NaN();
    Tl::drawResiduals(residuals_, {}, run_settings_.density_tolerance, run_settings_.velocity_tolerance.value_or(run_settings_.density_tolerance), tk,
                      "running");
    worker_ = std::thread([ch = ch_.get(), st = run_settings_, b = *sec] {
        Observer obs(*ch);
        try {
            T::Case c = T::build_case(st, b);
            obs.say(T::format_case_summary(c));
            obs.say(T::format_grid_summary(c));
            {
                std::lock_guard<std::mutex> g(ch->mu);
                ch->c = c;
            }
            T::RunResult r = T::run_case(c, &obs, &ch->cancel);
            obs.say(T::format_run_report(r));
            std::lock_guard<std::mutex> g(ch->mu);
            ch->result = std::move(r);
        } catch (const std::exception& e) {
            obs.say(std::string("the case cannot be run: ") + e.what());
            std::lock_guard<std::mutex> g(ch->mu);
            ch->error = e.what();
        }
        ch->done = true;
    });
    run_->setEnabled(false);
    run2_->setEnabled(false);
    stop_->setEnabled(true);
    stop2_->setEnabled(true);
    run_status_->setText(QString("running: %1, %2 thread(s)").arg(QString::fromStdString(run_settings_.case_name)).arg(run_settings_.threads));
    status_->setText("solving...");
    last_plot_ms_ = 0;
    poll_->start();
    showPage(kRun);
}

void TurbolenzaWindow::stopRun()
{
    if (ch_ && running()) {
        ch_->cancel = true;
        run_status_->setText("stopping after the current pass...");
    }
}

void TurbolenzaWindow::pollRun()
{
    if (!ch_) return;
    std::vector<std::string> lines;
    std::vector<T::PassRecord> recs;
    {
        std::lock_guard<std::mutex> g(ch_->mu);
        lines.swap(ch_->lines);
        recs.swap(ch_->records);
    }
    if (!lines.empty()) {
        QString t;
        for (const auto& l : lines) t += QString::fromStdString(l).trimmed() + "\n";
        t.chop(1);
        log_->appendPlainText(t);
    }
    records_.insert(records_.end(), recs.begin(), recs.end());
    const qint64 now = QDateTime::currentMSecsSinceEpoch();
    if (!recs.empty() && now - last_plot_ms_ >= 300) {  // at most ~3 redraws a second
        last_plot_ms_ = now;
        updateLivePlot(false);
    }
    if (ch_->done && worker_.joinable()) {
        worker_.join();
        runFinished();
    }
}

void TurbolenzaWindow::updateLivePlot(bool final)
{
    const auto& st = run_settings_;
    const double tk = st.solver_mode == T::SolverMode::CoupledKutta ? st.kutta_tolerance : std::numeric_limits<double>::quiet_NaN();
    QString status = QString("%1 passes").arg(records_.size());
    if (final && result_)
        status = result_->converged ? QString("converged in %1 passes").arg(result_->iterations)
                                    : QString("NOT CONVERGED after %1 passes: %2").arg(result_->iterations).arg(QString::fromStdString(result_->error));
    Tl::drawResiduals(residuals_, records_, st.density_tolerance, st.velocity_tolerance.value_or(st.density_tolerance), tk, status);
    if (final && result_) residuals_->setBadge({status, result_->converged ? 0 : 1});
}

void TurbolenzaWindow::runFinished()
{
    poll_->stop();
    pollRun();  // the last lines
    {
        std::lock_guard<std::mutex> g(ch_->mu);
        if (ch_->c) run_case_ = std::move(ch_->c);
        if (ch_->result) result_ = std::move(ch_->result);
        if (!ch_->error.empty()) run_status_->setText("the case cannot be run: " + QString::fromStdString(ch_->error));
    }
    run_->setEnabled(true);
    run2_->setEnabled(true);
    stop_->setEnabled(false);
    stop2_->setEnabled(false);
    if (!result_) {
        status_->setText("the case cannot be run");
        updateLivePlot(true);
        return;
    }
    const auto& r = *result_;
    run_status_->setText(r.converged ? QString("converged in %1 passes, %2 s").arg(r.iterations).arg(r.runtime.wall_time_s, 0, 'f', 2)
                                     : QString("NOT CONVERGED after %1 passes: %2").arg(r.iterations).arg(QString::fromStdString(r.error)));
    status_->setText(run_status_->text());
    report_->setText(QString::fromStdString(T::format_run_report(r)));
    updateLivePlot(true);  // the one history, drawn once more with its final status
    if (r.solution) openPostWindows();
}

// ── post-processing ──────────────────────────────────────────────────────────

void TurbolenzaWindow::openPostWindows()
{
    showResiduals();
    showContours();
    showLoading();
}

void TurbolenzaWindow::showResiduals()
{
    if (!res_win_) {
        res_win_ = new PostWindow("TURBOlenza: Residuals", this);
        res_win_->onClose = [this] {  // the plot goes back to the Run page
            run_plot_note_->hide();
            run_plot_box_->addWidget(residuals_);
            residuals_->show();
        };
    }
    if (residuals_->parentWidget() != res_win_) {
        res_win_->box->addWidget(residuals_);
        run_plot_note_->show();
    }
    res_win_->show();
    res_win_->raise();
}

void TurbolenzaWindow::showContours()
{
    if (!contour_win_) {
        contour_win_ = new PostWindow("TURBOlenza: Contours", this);
        auto* row = new QHBoxLayout;
        contour_field_ = new QComboBox(contour_win_);
        contour_field_->addItems({"static pressure", "relative Mach number"});
        bands_ = new QSpinBox(contour_win_);
        bands_->setRange(4, 48);
        bands_->setValue(16);
        isolines_ = new QCheckBox("iso-lines", contour_win_);
        isolines_->setChecked(true);
        row->addWidget(new QLabel("field", contour_win_));
        row->addWidget(contour_field_);
        row->addWidget(new QLabel("bands", contour_win_));
        row->addWidget(bands_);
        row->addWidget(isolines_);
        row->addStretch(1);
        contour_win_->box->addLayout(row);
        contour_plot_ = new PlotWidget(contour_win_);
        contour_win_->box->addWidget(contour_plot_, 1);
        connect(contour_field_, &QComboBox::currentIndexChanged, this, [this] { redrawContours(); });
        connect(bands_, &QSpinBox::valueChanged, this, [this] { redrawContours(); });
        connect(isolines_, &QCheckBox::toggled, this, [this] { redrawContours(); });
    }
    redrawContours();
    contour_win_->show();
    contour_win_->raise();
}

void TurbolenzaWindow::redrawContours()
{
    if (!contour_plot_) return;
    if (!result_ || !result_->solution) {
        contour_plot_->setBadge({"no solution yet: run a case", 1});
        contour_plot_->update();
        return;
    }
    const bool mach = contour_field_->currentIndex() == 1;
    const auto& r = *result_;
    const QString title = QString("%1: β_out %2°%3")
                              .arg(mach ? "Relative Mach number" : "Static pressure [Pa]")
                              .arg(r.solution->exit_angle * kDeg, 0, 'f', 3)
                              .arg(r.converged ? "" : "  [NOT CONVERGED]");
    Tl::drawContours(contour_plot_, *r.solution, mach ? Tl::ContourField::Mach : Tl::ContourField::Pressure, bands_->value(),
                     isolines_->isChecked(), title);
    contour_plot_->setBadge({});
}

void TurbolenzaWindow::showLoading()
{
    if (!loading_win_) {
        loading_win_ = new PostWindow("TURBOlenza: Blade loading", this);
        auto* row = new QHBoxLayout;
        loading_q_ = new QComboBox(loading_win_);
        loading_q_->addItems({"Mach number", "static pressure p", "Dunavant S = (Pt1 - p) / q1"});
        row->addWidget(new QLabel("quantity", loading_win_));
        row->addWidget(loading_q_);
        row->addStretch(1);
        loading_win_->box->addLayout(row);
        loading_plot_ = new PlotWidget(loading_win_);
        loading_win_->box->addWidget(loading_plot_, 1);
        connect(loading_q_, &QComboBox::currentIndexChanged, this, [this] { redrawLoading(); });
    }
    redrawLoading();
    loading_win_->show();
    loading_win_->raise();
}

void TurbolenzaWindow::redrawLoading()
{
    if (!loading_plot_) return;
    if (!result_ || !result_->solution) {
        loading_plot_->setBadge({"no solution yet: run a case", 1});
        loading_plot_->update();
        return;
    }
    const int q = loading_q_->currentIndex();
    const auto& r = *result_;
    const QString what = q == 0 ? "Blade-surface Mach number" : q == 1 ? "Blade-surface static pressure" : "Blade-surface pressure coefficient";
    Tl::drawLoading(loading_plot_, *r.solution, q == 0 ? Tl::LoadingQuantity::Mach : q == 1 ? Tl::LoadingQuantity::Pressure : Tl::LoadingQuantity::DunavantS,
                    QString("%1: β_out %2°%3").arg(what).arg(r.solution->exit_angle * kDeg, 0, 'f', 3).arg(r.converged ? "" : "  [NOT CONVERGED]"));
    loading_plot_->setBadge({});
}

QStringList TurbolenzaWindow::saveResults(const QString& dir, bool overwrite)
{
    QStringList out;
    if (!result_ || !run_case_) throw std::runtime_error("no run to save");
    // export_run returns the names of the files it wrote, relative to dir
    for (const auto& f : T::export_run(*result_, *run_case_, dir.toStdString(), overwrite))
        out << QDir(dir).filePath(QString::fromStdString(f));
    // the plots (each window laid out at its size even when it is hidden)
    showResiduals();
    showContours();
    showLoading();
    auto png = [&](PlotWidget* w, const QString& name) {
        const QString path = QDir(dir).filePath(name);
        if (!overwrite && QFileInfo::exists(path)) throw std::runtime_error(("file exists: " + path).toStdString());
        if (!w->image().save(path)) throw std::runtime_error(("cannot write " + path).toStdString());
        out << path;
    };
    QApplication::processEvents();
    png(residuals_, "residuals.png");
    const int f0 = contour_field_->currentIndex(), q0 = loading_q_->currentIndex();
    contour_field_->setCurrentIndex(0);
    png(contour_plot_, "contour_pressure.png");
    contour_field_->setCurrentIndex(1);
    png(contour_plot_, "contour_mach.png");
    contour_field_->setCurrentIndex(f0);
    const char* names[3] = {"loading_mach.png", "loading_pressure.png", "loading_dunavant_S.png"};
    for (int q = 0; q < 3; ++q) {
        loading_q_->setCurrentIndex(q);
        png(loading_plot_, names[q]);
    }
    loading_q_->setCurrentIndex(q0);
    return out;
}

void TurbolenzaWindow::saveDialog()
{
    if (!result_ || !run_case_) {
        QMessageBox::information(this, "TURBOlenza", "There is no run to save yet.");
        return;
    }
    const QString dir = QFileDialog::getExistingDirectory(this, "Save the TURBOlenza results into");
    if (dir.isEmpty()) return;
    bool exists = false;
    for (const auto& n : kExportNames) exists = exists || QFileInfo::exists(QDir(dir).filePath(n));
    bool overwrite = false;
    if (exists) {
        if (QMessageBox::question(this, "TURBOlenza", "The directory already holds results. Overwrite them?") != QMessageBox::Yes) return;
        overwrite = true;
    }
    try {
        const auto files = saveResults(dir, overwrite);
        status_->setText(QString("%1 files written to %2").arg(files.size()).arg(dir));
    } catch (const std::exception& e) {
        QMessageBox::warning(this, "TURBOlenza", QString("The results were not saved: %1").arg(e.what()));
    }
}

// ── self-test ────────────────────────────────────────────────────────────────

int TurbolenzaWindow::selfTest(const QString& dir, int passes)
{
    QDir().mkpath(dir);
    auto fail = [](const char* what, int code) {
        std::fprintf(stderr, "turbolenza self-test: %s\n", what);
        return code;
    };
    show();
    QApplication::processEvents();
    // 1. the Mid layer of the open design, copied
    takeLayer(2);
    if (layer_.name != s_.design.layers[2].name) return fail("the Mid layer was not taken", 21);
    const auto* sec = section();
    if (!sec) return fail(("the Mid layer gives no section: " + section_error_.toStdString()).c_str(), 22);
    std::printf("turbolenza self-test: section %s, SS %d / PS %d control points, fit residual %.2e / %.2e c\n", sec->designation.c_str(),
                fit_.suction.control_points, fit_.pressure.control_points, fit_.suction.max_residual, fit_.pressure.max_residual);
    QApplication::processEvents();
    stack_->widget(kBlade)->grab().save(dir + "/tl_blade.png");
    // 2. the cascade at the layer's stagger and solidity, the inlet at zero incidence, M1 0.5
    showPage(kCascade);
    if (!cascade_) return fail("the cascade was not built", 23);
    beta1_->setValue(std::round(10.0 * cascade_->kappa_le * kDeg) / 10.0);
    mach_given_->setChecked(true);
    mach1_->setValue(0.5);
    refreshCascade();
    QApplication::processEvents();
    stack_->widget(kCascade)->grab().save(dir + "/tl_cascade.png");
    // 3. a very_coarse mesh
    showPage(kMesh);
    level_->setCurrentText("very_coarse");
    if (!generateMesh(true)) return fail(("the very_coarse mesh failed: " + mesh_error_->text().toStdString()).c_str(), 24);
    QApplication::processEvents();
    stack_->widget(kMesh)->grab().save(dir + "/tl_mesh.png");
    // 4. a short solve
    showPage(kSolver);
    max_rho_->setValue(passes);
    max_k_->setValue(passes);
    QApplication::processEvents();
    stack_->widget(kSolver)->grab().save(dir + "/tl_solver.png");
    startRun();
    if (!running()) return fail("the solve did not start", 25);
    QElapsedTimer clock;
    clock.start();
    while (running() || (ch_ && !ch_->done)) {
        QApplication::processEvents(QEventLoop::AllEvents, 50);
        if (clock.elapsed() > 600000) {
            stopRun();
            return fail("the solve took more than 10 minutes", 26);
        }
    }
    QApplication::processEvents();
    if (!result_) return fail("the solve returned no result", 27);
    if (records_.empty()) return fail("no pass was reported", 28);
    if (records_.size() != result_->records.size()) return fail("the live history and the run's records differ", 29);
    std::printf("turbolenza self-test: %zu passes reported (%s)\n", records_.size(),
                result_->converged ? "converged" : result_->error.c_str());
    grab().save(dir + "/tl_run.png");
    if (!result_->solution) return fail("the run has no solution to post-process", 30);
    // 5. post-processing windows, saved results
    for (auto* w : {res_win_, contour_win_, loading_win_})
        if (!w || !w->isVisible()) return fail("a post-processing window did not open", 31);
    QStringList files;
    try {
        files = saveResults(dir, true);
    } catch (const std::exception& e) {
        return fail(e.what(), 32);
    }
    for (const auto& f : files)
        if (!QFileInfo(f).exists() || QFileInfo(f).size() == 0) return fail(("missing or empty: " + f.toStdString()).c_str(), 33);
    // the one residual widget: it is in the Residuals window and back on the Run page when that closes
    if (residuals_->parentWidget() != res_win_) return fail("the residual plot is not in the Residuals window", 34);
    res_win_->close();
    QApplication::processEvents();
    if (residuals_->window() != this) return fail("the residual plot did not return to the Run page", 35);
    std::printf("turbolenza self-test OK: %lld files in %s\n", static_cast<long long>(files.size()), qPrintable(dir));
    return 0;
}

}  // namespace Designer
