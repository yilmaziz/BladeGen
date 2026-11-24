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
#include <QTableWidget>
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

enum Page : int { kBlade = 0, kCascade, kMesh, kSolver, kRun, kPost, kStudy };
constexpr double kDeg = 180.0 / std::numbers::pi;
const QStringList kExportNames = {"surface_mach.csv", "surface_pressure.csv", "convergence_history.csv", "fields.csv", "run_summary.json",
                                  "bl_suction.csv", "bl_pressure.csv", "bl_history.csv"};

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

struct TurbolenzaWindow::StudyChannel
{
    std::mutex mu;
    std::vector<std::string> lines;
    std::vector<std::pair<std::string, T::PassRecord>> records;  ///< (level, pass), not yet taken by the GUI
    std::vector<T::MeshStudyEntry> finished;                    ///< levels done, not yet taken
    std::atomic<bool> cancel{false};
    std::atomic<bool> done{false};
    std::optional<T::MeshStudyResult> result;
    std::string error;  ///< the study itself failed
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
        const std::string bl = i.bl_iteration ? "[BL iteration " + std::to_string(*i.bl_iteration) + "] " : std::string();
        std::snprintf(b, sizeof b, "%s%s solve: %d x %d nodes, at most %d passes; tolerances rho %s, W %s, Kutta %s; mdot' %.6g (%s); "
                      "first beta2 %.4f deg; %d thread(s), %s",
                      bl.c_str(), i.mode.c_str(), i.n_m, i.n_eta, i.max_passes, tol(i.density_tolerance).c_str(), tol(i.velocity_tolerance).c_str(),
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
    stack_->addWidget(studyPage());
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
    connect(poll_, &QTimer::timeout, this, [this] {
        pollRun();
        pollStudy();
    });
    takeLayer(std::min(2, static_cast<int>(s_.design.layers.size()) - 1));
    showPage(kBlade);
}

TurbolenzaWindow::~TurbolenzaWindow()
{
    if (mesh_watch_) mesh_watch_->waitForFinished();
    if (ch_) ch_->cancel = true;
    if (sch_) sch_->cancel = true;
    if (worker_.joinable()) worker_.join();
    if (study_worker_.joinable()) study_worker_.join();
}

void TurbolenzaWindow::closeEvent(QCloseEvent* e)
{
    if (running() || studyRunning()) {
        const auto b = QMessageBox::question(this, "TURBOlenza", "A solve is running. Stop it and close?");
        if (b != QMessageBox::Yes) {
            e->ignore();
            return;
        }
        stopRun();
        stopStudy();
        if (worker_.joinable()) worker_.join();  // the solver stops after its current pass
        if (study_worker_.joinable()) study_worker_.join();
        pollRun();
        pollStudy();
    }
    for (auto* w : {res_win_, contour_win_, loading_win_, bl_win_, study_beta_win_, study_load_win_})
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
    auto* solver = add(pre, "Solver", kSolver, 0);
    add(solver, "Mesh independence", kStudy, 0);
    add(nullptr, "Run", kRun, 0);
    auto* post = add(nullptr, "Post-processing", kPost, 0);
    add(post, "Residuals", kPost, 1);
    add(post, "Contours", kPost, 2);
    add(post, "Blade loading", kPost, 3);
    add(post, "Boundary layer", kPost, 6);
    add(post, "Mesh study: beta_out", kPost, 4);
    add(post, "Mesh study: blade loading", kPost, 5);
    tree_->expandAll();
    connect(tree_, &QTreeWidget::itemClicked, this, [this](QTreeWidgetItem* it, int) {
        showPage(it->data(0, Qt::UserRole).toInt());
        const int which = it->data(0, Qt::UserRole + 1).toInt();
        if (which == 1) showResiduals();
        else if (which == 2) showContours();
        else if (which == 3) showLoading();
        else if (which == 4) showStudyBetaOut();
        else if (which == 5) showStudyLoading();
        else if (which == 6) showBoundaryLayer();
    });
}

void TurbolenzaWindow::buildMenus()
{
    auto* file = menuBar()->addMenu("&File");
    file->addAction("&Save results...", this, [this] { saveDialog(); });
    file->addAction("Save &mesh study...", this, [this] { saveStudyDialog(); });
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
    post->addAction("Boundary &layer", this, [this] { showBoundaryLayer(); });
    post->addAction("Mesh study: &beta_out", this, [this] { showStudyBetaOut(); });
    post->addAction("Mesh study: blade &loading", this, [this] { showStudyLoading(); });
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
    stream_.reset();
    stream_error_.clear();
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
    if (running() || studyRunning()) return;
    QString input_error;
    const auto input = caseInput(input_error);
    if (!input) {
        run_status_->setText(input_error);
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
    worker_ = std::thread([ch = ch_.get(), st = run_settings_, in = *input] {
        Observer obs(*ch);
        try {
            T::Case c = in.stream ? T::build_case(st, *in.stream) : T::build_case(st, in.section);
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
    updateRunButtons();
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
    updateRunButtons();
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

void TurbolenzaWindow::updateRunButtons()
{
    const bool run = running(), study = studyRunning();
    for (auto* b : {run_, run2_, study_run_}) b->setEnabled(!run && !study);
    stop_->setEnabled(run);
    stop2_->setEnabled(run);
    study_stop_->setEnabled(study);
    study_save_->setEnabled(!study && study_ && !study_->levels.empty());
}

// ── mesh study ───────────────────────────────────────────────────────────────

void TurbolenzaWindow::startStudy()
{
    if (running() || studyRunning()) return;
    std::vector<std::string> levels;
    for (std::size_t k = 0; k < study_check_.size() && k < T::kGridLevels.size(); ++k)
        if (study_check_[k]->isChecked()) levels.emplace_back(T::kGridLevels[k]);
    showPage(kStudy);
    if (levels.empty()) {
        study_status_->setText("check at least one grid level");
        return;
    }
    QString input_error;
    const auto input = caseInput(input_error);
    if (!input) {
        study_status_->setText(input_error);
        return;
    }
    study_settings_ = settings();
    study_ = T::MeshStudyResult{};
    study_records_.clear();
    study_level_.clear();
    study_table_->setRowCount(0);
    log_->clear();
    sch_ = std::make_unique<StudyChannel>();
    redrawStudyResiduals("running");
    study_worker_ = std::thread([ch = sch_.get(), st = study_settings_, in = *input, levels] {
        auto say = [ch](const std::string& line) {
            std::cout << line << std::endl;
            std::lock_guard<std::mutex> g(ch->mu);
            ch->lines.push_back(line);
        };
        const int n = static_cast<int>(levels.size());
        T::MeshStudyCallbacks cb;
        cb.level_started = [&](const std::string& level, int k, int n_m, int n_eta) {
            char t[256];
            std::snprintf(t, sizeof t, "[%s] mesh study level %d of %d: %d x %d nodes", level.c_str(), k + 1, n, n_m, n_eta);
            say(t);
        };
        cb.pass = [&](const std::string& level, const T::PassRecord& r) {
            const std::string line = "[" + level + "] " + T::format_record(r);
            std::cout << line << std::endl;
            std::lock_guard<std::mutex> g(ch->mu);
            ch->lines.push_back(line);
            ch->records.emplace_back(level, r);
        };
        cb.level_finished = [&](const T::MeshStudyEntry& e) {
            char t[1024];
            if (e.converged)
                std::snprintf(t, sizeof t, "[%s] converged in %d passes: beta_out %.4f deg, max Mach %.4f, %.2f s", e.level.c_str(), e.passes,
                              e.beta_out_deg, e.mach_max, e.wall_time_s);
            else
                std::snprintf(t, sizeof t, "[%s] NOT CONVERGED after %d passes (%.2f s): %s", e.level.c_str(), e.passes, e.wall_time_s,
                              e.error.c_str());
            say(t);
            std::lock_guard<std::mutex> g(ch->mu);
            ch->finished.push_back(e);
        };
        try {
            T::MeshStudyResult r =
                T::run_mesh_study(st, in.section, levels, cb, &ch->cancel, in.stream ? &*in.stream : nullptr);
            say(r.cancelled ? "mesh study stopped: the remaining levels were skipped" : "mesh study finished");
            std::lock_guard<std::mutex> g(ch->mu);
            ch->result = std::move(r);
        } catch (const std::exception& e) {
            say(std::string("the mesh study failed: ") + e.what());
            std::lock_guard<std::mutex> g(ch->mu);
            ch->error = e.what();
        }
        ch->done = true;
    });
    updateRunButtons();
    study_status_->setText(QString("mesh study running: %1 level(s), %2 thread(s)").arg(levels.size()).arg(study_settings_.threads));
    status_->setText("mesh study running...");
    study_plot_ms_ = 0;
    poll_->start();
}

void TurbolenzaWindow::stopStudy()
{
    if (sch_ && studyRunning()) {
        sch_->cancel = true;
        study_status_->setText("stopping after the current pass...");
    }
}

void TurbolenzaWindow::pollStudy()
{
    if (!sch_) return;
    std::vector<std::string> lines;
    std::vector<std::pair<std::string, T::PassRecord>> recs;
    std::vector<T::MeshStudyEntry> done;
    {
        std::lock_guard<std::mutex> g(sch_->mu);
        lines.swap(sch_->lines);
        recs.swap(sch_->records);
        done.swap(sch_->finished);
    }
    if (!lines.empty()) {
        QString t;
        for (const auto& l : lines) t += QString::fromStdString(l).trimmed() + "\n";
        t.chop(1);
        log_->appendPlainText(t);
    }
    for (auto& [level, r] : recs) {
        if (level != study_level_) {  // a new level: its own live history
            study_level_ = level;
            study_records_.clear();
        }
        study_records_.push_back(r);
    }
    const qint64 now = QDateTime::currentMSecsSinceEpoch();
    if (!recs.empty() && now - study_plot_ms_ >= 300) {  // at most ~3 redraws a second
        study_plot_ms_ = now;
        redrawStudyResiduals(QString("[%1] %2 passes").arg(QString::fromStdString(study_level_)).arg(study_records_.size()));
    }
    for (auto& e : done) {
        addStudyRow(e);
        if (study_) study_->levels.push_back(std::move(e));
    }
    if (sch_->done && study_worker_.joinable()) {
        study_worker_.join();
        studyFinished();
    }
}

void TurbolenzaWindow::studyFinished()
{
    poll_->stop();
    pollStudy();  // the last lines and levels
    std::string error;
    {
        std::lock_guard<std::mutex> g(sch_->mu);
        if (sch_->result) study_ = std::move(sch_->result);  // the same levels, in the order run
        error = sch_->error;
    }
    updateRunButtons();
    int converged = 0;
    const int n = study_ ? static_cast<int>(study_->levels.size()) : 0;
    if (study_)
        for (const auto& e : study_->levels) converged += e.converged ? 1 : 0;
    if (!error.empty())
        study_status_->setText("the mesh study failed: " + QString::fromStdString(error));
    else
        study_status_->setText(QString("mesh study %1: %2 of %3 level(s) converged")
                                   .arg(study_ && study_->cancelled ? "stopped" : "finished")
                                   .arg(converged)
                                   .arg(n));
    status_->setText(study_status_->text());
    redrawStudyResiduals(QString("[%1] %2 passes").arg(QString::fromStdString(study_level_)).arg(study_records_.size()));
    if (n > 0) {
        showStudyBetaOut();
        showStudyLoading();
    }
}

void TurbolenzaWindow::addStudyRow(const T::MeshStudyEntry& e)
{
    const int row = study_table_->rowCount();
    study_table_->insertRow(row);
    auto f = [](double v, int d) { return std::isfinite(v) ? QString::number(v, 'f', d) : QString("-"); };
    auto cell = [&](int col, const QString& t) {
        auto* it = new QTableWidgetItem(t);
        study_table_->setItem(row, col, it);
        return it;
    };
    cell(0, QString::fromStdString(e.level));
    cell(1, e.n_m > 0 ? QString("%1 x %2").arg(e.n_m).arg(e.n_eta) : QString("-"));
    cell(2, f(e.dm_over_chord, 5));
    cell(3, QString::number(e.passes));
    auto* st = cell(4, e.converged ? QString("converged") : QString::fromStdString(e.error));
    if (!e.converged) {
        st->setForeground(QColor(200, 0, 0));
        st->setToolTip(QString::fromStdString(e.error));
    }
    cell(5, f(e.beta_out_deg, 4));
    cell(6, f(e.mach_max, 4));
    cell(7, f(e.wall_time_s, 2));
}

void TurbolenzaWindow::redrawStudyResiduals(const QString& status)
{
    const auto& st = study_settings_;
    const double tk = st.solver_mode == T::SolverMode::CoupledKutta ? st.kutta_tolerance : std::numeric_limits<double>::quiet_NaN();
    Tl::drawResiduals(study_residuals_, study_records_, st.density_tolerance, st.velocity_tolerance.value_or(st.density_tolerance), tk, status);
}

void TurbolenzaWindow::showStudyBetaOut()
{
    if (!study_beta_win_) {
        study_beta_win_ = new PostWindow("TURBOlenza: Mesh study, beta_out", this);
        auto* row = new QHBoxLayout;
        study_x_ = new QComboBox(study_beta_win_);
        study_x_->addItems({"dm / c", "total nodes"});
        row->addWidget(new QLabel("x axis", study_beta_win_));
        row->addWidget(study_x_);
        row->addStretch(1);
        study_beta_win_->box->addLayout(row);
        study_beta_plot_ = new PlotWidget(study_beta_win_);
        study_beta_win_->box->addWidget(study_beta_plot_, 1);
        connect(study_x_, &QComboBox::currentIndexChanged, this, [this] { redrawStudyBetaOut(); });
    }
    redrawStudyBetaOut();
    study_beta_win_->show();
    study_beta_win_->raise();
}

void TurbolenzaWindow::redrawStudyBetaOut()
{
    if (!study_beta_plot_) return;
    if (!study_ || study_->levels.empty()) {
        study_beta_plot_->setBadge({"no mesh study yet", 1});
        study_beta_plot_->update();
        return;
    }
    Tl::drawStudyBetaOut(study_beta_plot_, study_->levels, study_x_->currentIndex() == 1);
}

void TurbolenzaWindow::showStudyLoading()
{
    if (!study_load_win_) {
        study_load_win_ = new PostWindow("TURBOlenza: Mesh study, blade loading", this);
        auto* row = new QHBoxLayout;
        study_q_ = new QComboBox(study_load_win_);
        study_q_->addItems({"Mach number", "static pressure p", "Dunavant S = (Pt1 - p) / q1"});
        row->addWidget(new QLabel("quantity", study_load_win_));
        row->addWidget(study_q_);
        row->addStretch(1);
        study_load_win_->box->addLayout(row);
        study_load_plot_ = new PlotWidget(study_load_win_);
        study_load_win_->box->addWidget(study_load_plot_, 1);
        connect(study_q_, &QComboBox::currentIndexChanged, this, [this] { redrawStudyLoading(); });
    }
    redrawStudyLoading();
    study_load_win_->show();
    study_load_win_->raise();
}

void TurbolenzaWindow::redrawStudyLoading()
{
    if (!study_load_plot_) return;
    if (!study_ || study_->levels.empty()) {
        study_load_plot_->setBadge({"no mesh study yet", 1});
        study_load_plot_->update();
        return;
    }
    const int q = study_q_->currentIndex();
    Tl::drawStudyLoading(study_load_plot_, study_->levels,
                         q == 0 ? Tl::LoadingQuantity::Mach : q == 1 ? Tl::LoadingQuantity::Pressure : Tl::LoadingQuantity::DunavantS);
}

QStringList TurbolenzaWindow::saveStudy(const QString& dir, bool overwrite)
{
    QStringList out;
    if (!study_ || study_->levels.empty()) throw std::runtime_error("no mesh study to save");
    // write_mesh_study returns the files it wrote, relative to dir
    for (const auto& f : T::write_mesh_study(*study_, dir.toStdString(), overwrite)) out << QDir(dir).filePath(QString::fromStdString(f));
    showStudyBetaOut();
    showStudyLoading();
    auto png = [&](PlotWidget* w, const QString& name) {
        const QString path = QDir(dir).filePath(name);
        if (!overwrite && QFileInfo::exists(path)) throw std::runtime_error(("file exists: " + path).toStdString());
        if (!w->image().save(path)) throw std::runtime_error(("cannot write " + path).toStdString());
        out << path;
    };
    QApplication::processEvents();
    const int x0 = study_x_->currentIndex(), q0 = study_q_->currentIndex();
    study_x_->setCurrentIndex(0);
    png(study_beta_plot_, "mesh_study_beta_out.png");
    study_x_->setCurrentIndex(x0);
    const char* names[3] = {"mesh_study_loading_mach.png", "mesh_study_loading_pressure.png", "mesh_study_loading_dunavant_S.png"};
    for (int q = 0; q < 3; ++q) {
        study_q_->setCurrentIndex(q);
        png(study_load_plot_, names[q]);
    }
    study_q_->setCurrentIndex(q0);
    return out;
}

void TurbolenzaWindow::saveStudyDialog()
{
    if (!study_ || study_->levels.empty() || studyRunning()) {
        QMessageBox::information(this, "TURBOlenza", "There is no finished mesh study to save yet.");
        return;
    }
    const QString dir = QFileDialog::getExistingDirectory(this, "Save the TURBOlenza mesh study into");
    if (dir.isEmpty()) return;
    bool exists = QFileInfo::exists(QDir(dir).filePath("mesh_study.csv"));
    for (const auto& e : study_->levels) exists = exists || QFileInfo::exists(QDir(dir).filePath(QString::fromStdString(e.level)));
    bool overwrite = false;
    if (exists) {
        if (QMessageBox::question(this, "TURBOlenza", "The directory already holds a mesh study. Overwrite it?") != QMessageBox::Yes) return;
        overwrite = true;
    }
    try {
        const auto files = saveStudy(dir, overwrite);
        status_->setText(QString("%1 files written to %2").arg(files.size()).arg(dir));
    } catch (const std::exception& e) {
        QMessageBox::warning(this, "TURBOlenza", QString("The mesh study was not saved: %1").arg(e.what()));
    }
}

// ── post-processing ──────────────────────────────────────────────────────────

void TurbolenzaWindow::openPostWindows()
{
    showResiduals();
    showContours();
    showLoading();
    if (result_ && result_->viscous) showBoundaryLayer();
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
        loading_inviscid_ = new QCheckBox("also the inviscid loading (BL iteration 0)", loading_win_);
        loading_inviscid_->setChecked(true);
        loading_inviscid_->setToolTip("a viscous run: the reported loading is the coupled one; this adds the inviscid "
                                      "solution's (dashed)");
        row->insertWidget(row->count() - 1, loading_inviscid_);
        loading_plot_ = new PlotWidget(loading_win_);
        loading_win_->box->addWidget(loading_plot_, 1);
        connect(loading_q_, &QComboBox::currentIndexChanged, this, [this] { redrawLoading(); });
        connect(loading_inviscid_, &QCheckBox::toggled, this, [this] { redrawLoading(); });
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
    const bool coupled = r.viscous && r.viscous->coupling == "two_way";
    loading_inviscid_->setEnabled(coupled);
    QString title = QString("%1: β_out %2°%3").arg(what).arg(r.solution->exit_angle * kDeg, 0, 'f', 3).arg(r.converged ? "" : "  [NOT CONVERGED]");
    if (coupled) title += QString("  (coupled; inviscid β_out %1°)").arg(r.viscous->inviscid_exit_angle * kDeg, 0, 'f', 3);
    Tl::drawLoading(loading_plot_, *r.solution, q == 0 ? Tl::LoadingQuantity::Mach : q == 1 ? Tl::LoadingQuantity::Pressure : Tl::LoadingQuantity::DunavantS,
                    title, coupled && loading_inviscid_->isChecked() ? &*r.viscous : nullptr);
    loading_plot_->setBadge({});
}

bool TurbolenzaWindow::referenceCase() const
{
    const T::CaseSettings& st = run_case_ ? run_case_->settings : run_settings_;
    return std::abs(st.inlet_flow_angle_deg - 45.0) < 0.05 && st.inlet_mach && std::abs(*st.inlet_mach - 0.65) < 0.005 &&
           std::abs(st.solidity - 1.0) < 0.005 && std::abs(st.stagger_deg - 33.0) < 0.05;
}

void TurbolenzaWindow::showBoundaryLayer()
{
    if (!bl_win_) {
        bl_win_ = new PostWindow("TURBOlenza: Boundary layer (Aungier 2003 Sec. 5.6)", this);
        auto* row = new QHBoxLayout;
        bl_q_ = new QComboBox(bl_win_);
        bl_q_->addItems({"shape factor H", "θ / c", "δ* / c", "skin friction c_f", "Re_θ"});
        bl_overlay_ = new QCheckBox("Aungier 2003 Fig. 5-15 (H)", bl_win_);
        bl_overlay_->setToolTip("the digitized H of Aungier's Fig. 5-15 (NACA 65(12A10)-10, β1 45°, M1 0.65), shown by "
                                "default for that case");
        row->addWidget(new QLabel("quantity", bl_win_));
        row->addWidget(bl_q_);
        row->addWidget(bl_overlay_);
        row->addStretch(1);
        bl_win_->box->addLayout(row);
        auto* split = new QSplitter(Qt::Vertical, bl_win_);
        bl_plot_ = new PlotWidget(split);
        auto* hist = new QSplitter(Qt::Horizontal, split);
        bl_change_plot_ = new PlotWidget(hist);
        bl_beta_plot_ = new PlotWidget(hist);
        split->addWidget(bl_plot_);
        split->addWidget(hist);
        split->setStretchFactor(0, 3);
        split->setStretchFactor(1, 2);
        bl_win_->box->addWidget(split, 1);
        connect(bl_q_, &QComboBox::currentIndexChanged, this, [this] { redrawBoundaryLayer(); });
        connect(bl_overlay_, &QCheckBox::toggled, this, [this] { redrawBoundaryLayer(); });
    }
    bl_overlay_->blockSignals(true);
    bl_overlay_->setChecked(referenceCase());
    bl_overlay_->blockSignals(false);
    redrawBoundaryLayer();
    bl_win_->show();
    bl_win_->raise();
}

void TurbolenzaWindow::redrawBoundaryLayer()
{
    if (!bl_plot_) return;
    if (!result_ || !result_->viscous) {
        for (auto* w : {bl_plot_, bl_change_plot_, bl_beta_plot_}) {
            w->setBadge({"no boundary layer: run a case with the boundary layer on (Solver page)", 1});
            w->update();
        }
        return;
    }
    const auto& r = *result_;
    const auto& v = *r.viscous;
    const int q = bl_q_->currentIndex();
    const T::BlIteration* last = v.history.empty() ? nullptr : &v.history.back();
    auto xc = [](double x) { return std::isfinite(x) ? QString::number(x, 'f', 3) : QString("none"); };
    // history.back() is the iteration whose layers are shown (layers_iteration): the values match the curves
    QString title = QString("Boundary layer, %1 (%2, %3), BL iteration %10: transition SS %4 / PS %5, separation SS %6 / PS %7, "
                            "loss ω %8%9")
                        .arg(QString::fromStdString(v.coupling))
                        .arg(QString::fromStdString(v.laminar_cf_form))
                        .arg(QString::fromStdString(v.wake_model))
                        .arg(xc(last ? last->transition_x_c_ss : std::nan("")))
                        .arg(xc(last ? last->transition_x_c_ps : std::nan("")))
                        .arg(xc(last ? last->separation_x_c_ss : std::nan("")))
                        .arg(xc(last ? last->separation_x_c_ps : std::nan("")))
                        .arg(v.loss_omega, 0, 'f', 4)
                        .arg(r.converged ? "" : "  [NOT CONVERGED]")
                        .arg(v.layers_iteration);
    const Tl::BlQuantity which = q == 0   ? Tl::BlQuantity::H
                                 : q == 1 ? Tl::BlQuantity::ThetaOverC
                                 : q == 2 ? Tl::BlQuantity::DeltaStarOverC
                                 : q == 3 ? Tl::BlQuantity::Cf
                                          : Tl::BlQuantity::ReTheta;
    Tl::drawBoundaryLayer(bl_plot_, v, run_case_ ? run_case_->settings.chord : 1.0, which, bl_overlay_->isChecked(), title);
    Tl::drawBlChange(bl_change_plot_, v, run_case_ ? run_case_->settings.bl_tolerance : std::nan(""));
    Tl::drawBlBetaOut(bl_beta_plot_, v);
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
    if (result_->viscous) {
        showBoundaryLayer();
        QApplication::processEvents();
        const int b0 = bl_q_->currentIndex();
        const char* bl_names[5] = {"bl_H.png", "bl_theta.png", "bl_delta_star.png", "bl_cf.png", "bl_re_theta.png"};
        for (int q = 0; q < 5; ++q) {
            bl_q_->setCurrentIndex(q);
            png(bl_plot_, bl_names[q]);
        }
        bl_q_->setCurrentIndex(b0);
        png(bl_change_plot_, "bl_history_change.png");
        png(bl_beta_plot_, "bl_history_beta_out.png");
    }
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
    // 6. the boundary layer (one way) on the same case, solved to convergence: its window and files
    showPage(kSolver);
    viscous_->setChecked(true);
    bl_coupling_->setCurrentIndex(0);
    max_rho_->setValue(5000);
    max_k_->setValue(5000);
    startRun();
    if (!running()) return fail("the viscous solve did not start", 36);
    clock.restart();
    while (running() || (ch_ && !ch_->done)) {
        QApplication::processEvents(QEventLoop::AllEvents, 50);
        if (clock.elapsed() > 600000) {
            stopRun();
            return fail("the viscous solve took more than 10 minutes", 37);
        }
    }
    QApplication::processEvents();
    if (!result_ || !result_->viscous || result_->viscous->suction.x_c.empty())
        return fail(("the viscous run has no boundary layer: " + (result_ ? result_->error : std::string("no result"))).c_str(), 38);
    if (!bl_win_ || !bl_win_->isVisible()) return fail("the Boundary layer window did not open", 39);
    std::printf("turbolenza self-test: one-way boundary layer, transition SS %.3f / PS %.3f, loss %.4f\n",
                result_->viscous->suction.transition_x_c, result_->viscous->pressure_side.transition_x_c,
                result_->viscous->loss_omega);
    const QString vdir = QDir(dir).filePath("viscous");
    QDir().mkpath(vdir);
    QStringList vfiles;
    try {
        vfiles = saveResults(vdir, true);
    } catch (const std::exception& e) {
        return fail(e.what(), 40);
    }
    for (const char* name : {"bl_suction.csv", "bl_pressure.csv", "bl_history.csv", "bl_H.png", "bl_history_change.png"})
        if (!QFileInfo(QDir(vdir).filePath(name)).exists()) return fail((std::string("missing: ") + name).c_str(), 41);
    viscous_->setChecked(false);
    // 6. not in v1: the Mid layer on its stream surface, a rotating row (rpm for beta1 about kappa_LE), its mesh
    showPage(kCascade);
    surface_->setCurrentIndex(1);
    row_->setCurrentIndex(1);
    mach_abs_->setValue(0.3);
    alpha1_->setValue(0.0);
    const auto* blade = streamBlade();
    if (!blade) return fail(("the layer gives no stream-surface blade: " + stream_error_.toStdString()).c_str(), 42);
    {
        const T::Cascade c = T::make_cascade(*blade);
        const double r_in = blade->chord * blade->sheet->radius(-ext_up_->value());
        const T::VelocityTriangle t0 = T::velocity_triangle(1000.0, 0.0, 0.3, pt_->value(), tt_->value(),
                                                            T::PerfectGas(gamma_->value(), gas_r_->value()), r_in);
        rpm_->setValue(std::round(t0.C_m * std::tan(c.kappa_le) / r_in * 60.0 / (2.0 * std::numbers::pi)));
    }
    refreshCascade();
    QApplication::processEvents();
    stack_->widget(kCascade)->grab().save(dir + "/tl_cascade_rotor.png");
    showPage(kMesh);
    if (!generateMesh(true)) return fail(("the rotor mesh failed: " + mesh_error_->text().toStdString()).c_str(), 43);
    if (!mesh_case_ || !mesh_case_->triangle) return fail("the rotor case has no velocity triangle", 44);
    std::printf("turbolenza self-test: rotor %s at %.0f rpm, beta1 %.4f deg, M1_rel %.4f, %d x %d nodes\n",
                blade->source.c_str(), rpm_->value(), mesh_case_->triangle->beta * kDeg, mesh_case_->triangle->mach_relative,
                mesh_case_->grid.n_m(), mesh_case_->grid.n_eta());
    QApplication::processEvents();
    stack_->widget(kMesh)->grab().save(dir + "/tl_mesh_rotor.png");
    surface_->setCurrentIndex(0);
    row_->setCurrentIndex(0);
    std::printf("turbolenza self-test OK: %lld + %lld files in %s\n", static_cast<long long>(files.size()),
                static_cast<long long>(vfiles.size()), qPrintable(dir));
    return 0;
}

}  // namespace Designer
