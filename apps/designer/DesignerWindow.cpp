#include "DesignerWindow.hpp"

#include "ParameterPanel.hpp"
#include "design/DesignSolid.hpp"

#include <QApplication>
#include <QCloseEvent>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDockWidget>
#include <QEvent>
#include <QFileDialog>
#include <QFormLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QMenuBar>
#include <QMessageBox>
#include <QScrollArea>
#include <QSpinBox>
#include <QStatusBar>
#include <QToolBar>
#include <QTreeWidget>
#include <QtConcurrent/QtConcurrentRun>

#include <algorithm>
#include <cstdio>

namespace Designer {

using namespace BladeGen;

DesignerWindow::DesignerWindow()
{
    setWindowTitle("BladeGen Designer");
    mdi_ = new QMdiArea(this);
    mdi_->setViewMode(QMdiArea::SubWindowView);
    mdi_->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    mdi_->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    setCentralWidget(mdi_);
    mdi_->viewport()->installEventFilter(this);  // the equal-size grid follows the viewport size

    panel_ = new ParameterPanel(s_, this);
    auto* scroll = new QScrollArea(this);
    scroll->setWidget(panel_);
    scroll->setWidgetResizable(true);
    scroll->setMinimumWidth(430);
    auto* dock = new QDockWidget("Design parameters", this);
    dock->setObjectName("parameters");
    dock->setWidget(scroll);
    addDockWidget(Qt::LeftDockWidgetArea, dock);

    issues_ = new QTreeWidget(this);
    issues_->setHeaderLabels({"layer / issue", "severity", "where"});
    issues_->header()->setSectionResizeMode(0, QHeaderView::Stretch);
    auto* idock = new QDockWidget("Geometry checks", this);
    idock->setObjectName("issues");
    idock->setWidget(issues_);
    addDockWidget(Qt::BottomDockWidgetArea, idock);
    connect(issues_, &QTreeWidget::itemClicked, this, [this](QTreeWidgetItem* it) {
        const int k = it->data(0, Qt::UserRole).toInt();
        if (k >= 0 && k < static_cast<int>(s_.design.layers.size()) && k != s_.active_layer) {
            s_.active_layer = k;
            panel_->sync();
            refreshViews();
        }
    });

    status_ = new QLabel(this);
    statusBar()->addPermanentWidget(status_, 1);
    fit_watch_ = new QFutureWatcher<QString>(this);
    connect(fit_watch_, &QFutureWatcher<QString>::finished, this, [this] {
        const QString msg = fit_watch_->result();
        panel_->setFitStatus(msg, false);
        panel_->sync();
        submit();
    });

    connect(panel_, &ParameterPanel::designChanged, this, [this] {
        dirty_ = true;
        updateTitle();
        submit();
    });
    connect(panel_, &ParameterPanel::activeLayerChanged, this, [this] { refreshViews(); });
    connect(panel_, &ParameterPanel::fitRequested, this, &DesignerWindow::startFit);
    connect(&engine_, &Engine::finished, this, &DesignerWindow::onResult, Qt::QueuedConnection);

    buildMenus();
    defaultLayout();
    resize(1700, 1000);
    submit();
}

DesignerWindow::~DesignerWindow()
{
    if (fit_watch_->isRunning()) fit_watch_->waitForFinished();
}

void DesignerWindow::closeEvent(QCloseEvent* e)
{
    if (dirty_) {
        const auto b = QMessageBox::question(this, "BladeGen Designer", "Save the changes to the design?",
                                             QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
        if (b == QMessageBox::Cancel || (b == QMessageBox::Save && !save(false))) {
            e->ignore();
            return;
        }
    }
    e->accept();
}

void DesignerWindow::buildMenus()
{
    auto* file = menuBar()->addMenu("&File");
    file->addAction("&New design (Rotor 67 custom layers)", QKeySequence::New, this, [this] { newDesign(false); });
    file->addAction("New design (NACA 65 layers)", this, [this] { newDesign(true); });
    file->addAction("&Open design...", QKeySequence::Open, this, &DesignerWindow::openDialog);
    file->addAction("&Save", QKeySequence::Save, this, [this] { save(false); });
    file->addAction("Save &as...", QKeySequence::SaveAs, this, [this] { save(true); });
    file->addSeparator();
    auto* ex = file->addMenu("&Export");
    ex->addAction("Blade solid (STEP / IGES / BREP, OpenCASCADE)...", this, &DesignerWindow::exportStep);
    ex->addAction("Blade, hub and casing mesh (OBJ)...", this, &DesignerWindow::exportObj);
    ex->addAction("Layer coordinates (CSV)...", this, &DesignerWindow::exportCsv);
    ex->addAction("Window image (PNG)...", this, &DesignerWindow::exportImage);
    file->addSeparator();
    file->addAction("E&xit", QKeySequence::Quit, this, &QWidget::close);

    auto* view = menuBar()->addMenu("&View");
    auto* nv = view->addMenu("&New view");
    for (const auto& grp : contentGroups()) {
        for (Content c : grp) nv->addAction(contentName(c), this, [this, c] { addView(c); });
        nv->addSeparator();
    }
    view->addAction("6 views (3 × 2)", this, [this] { gridLayout(6); });
    view->addAction("4 views (2 × 2)", this, [this] { gridLayout(4); });
    view->addAction("2 views (side by side)", this, [this] { gridLayout(2); });
    view->addAction("&Tile", this, [this] {
        grid_n_ = 0;
        mdi_->tileSubWindows();
    });
    view->addAction("&Cascade", this, [this] {
        grid_n_ = 0;
        mdi_->cascadeSubWindows();
    });
    view->addAction("&Default layout", this, &DesignerWindow::defaultLayout);
    view->addAction("Zoom fit all", this, [this] { refreshViews(true); });

    auto* help = menuBar()->addMenu("&Help");
    help->addAction("&About", this, [this] {
        QMessageBox::about(this, "BladeGen Designer",
                           QString("<b>BladeGen Designer %1</b><p>Parametric C2 blades in the (m', θ) plane, mapped conformally "
                                   "onto the stream surfaces of the meridional channel (Voß et al., Computers and Fluids 292, 2025).</p>"
                                   "<p>Right-click a view to choose its content, layer and plot options.</p>")
                               .arg(BLADEGEN_VERSION));
    });

    auto* tb = addToolBar("main");
    tb->setObjectName("toolbar");
    tb->addAction("New", this, [this] { newDesign(false); });
    tb->addAction("Open", this, &DesignerWindow::openDialog);
    tb->addAction("Save", this, [this] { save(false); });
    tb->addSeparator();
    tb->addAction("New view", this, [this] { addView(Content::BladeToBlade); });
    tb->addAction("3D Solid", this, [this] { addView(Content::Solid3D); })->setToolTip("new view: the OpenCASCADE blade solid with hub and casing");
    tb->addSeparator();
    for (const int n : {6, 4, 2})
        tb->addAction(QString("%1 views").arg(n), this, [this, n] { gridLayout(n); })
            ->setToolTip(n == 6 ? "six views in equal cells (3 × 2)" : n == 4 ? "four views in equal cells (2 × 2)" : "two views side by side");
    tb->addAction("Tile", this, [this] {
        grid_n_ = 0;
        mdi_->tileSubWindows();
    });
    tb->addSeparator();
    tb->addAction("Export STEP", this, &DesignerWindow::exportStep);
}

ContentView* DesignerWindow::addView(Content c)
{
    auto* v = new ContentView(s_, c);
    auto* sub = mdi_->addSubWindow(v);
    sub->setAttribute(Qt::WA_DeleteOnClose);
    sub->setWindowTitle(v->windowTitle());
    connect(v, &ContentView::titleChanged, sub, &QMdiSubWindow::setWindowTitle);
    connect(v, &ContentView::newViewWanted, this, [this](Content cc) { addView(cc); });
    connect(v, &ContentView::meshPointsWanted, this, &DesignerWindow::meshPointsDialog);
    connect(v, &ContentView::designEdited, this, [this](bool done) {
        dirty_ = true;
        updateTitle();
        if (done) panel_->sync();  // the knobs follow the dragged control point on release
        submit();
    });
    sub->resize(560, 400);
    sub->show();
    return v;
}

QList<ContentView*> DesignerWindow::views() const
{
    QList<ContentView*> out;
    for (auto* w : mdi_->subWindowList())
        if (auto* v = qobject_cast<ContentView*>(w->widget())) out << v;
    return out;
}

void DesignerWindow::defaultLayout()
{
    mdi_->closeAllSubWindows();
    for (Content c : {Content::MeridionalBeta, Content::View3D, Content::BladeToBlade, Content::ProfileCurvature,
                      Content::BetaAll, Content::ThicknessAll})
        addView(c);
    gridLayout(6);
}

void DesignerWindow::gridLayout(int n)
{
    n = std::clamp(n, 1, 6);
    auto subs = mdi_->subWindowList(QMdiArea::CreationOrder);
    if (auto* act = mdi_->activeSubWindow(); act && subs.indexOf(act) >= n) {  // keep the active view
        subs.removeOne(act);
        subs.insert(n - 1, act);
    }
    for (int i = n; i < subs.size(); ++i) subs[i]->close();
    // fill up with the default contents not shown yet
    QList<Content> shown;
    for (int i = 0; i < std::min(n, static_cast<int>(subs.size())); ++i)
        if (auto* v = qobject_cast<ContentView*>(subs[i]->widget())) shown << v->content();
    for (Content c : {Content::MeridionalBeta, Content::View3D, Content::BladeToBlade, Content::ProfileCurvature, Content::BetaAll,
                      Content::ThicknessAll}) {
        if (static_cast<int>(shown.size()) >= n) break;
        if (shown.contains(c)) continue;
        addView(c);
        shown << c;
    }
    grid_n_ = n;
    arrangeGrid();
}

void DesignerWindow::arrangeGrid()
{
    if (grid_n_ <= 0) return;
    const auto subs = mdi_->subWindowList(QMdiArea::CreationOrder);
    const int cols = grid_n_ == 6 ? 3 : grid_n_ == 1 ? 1 : 2, rows = (grid_n_ + cols - 1) / cols;
    const QSize a = mdi_->viewport()->size();
    int k = 0;
    // cells in creation order
    for (auto* s : subs) {
        if (k >= grid_n_) break;
        if (s->isHidden()) continue;  // closed, not deleted yet
        const int r = k / cols, c = k % cols;
        const int x0 = a.width() * c / cols, x1 = a.width() * (c + 1) / cols;
        const int y0 = a.height() * r / rows, y1 = a.height() * (r + 1) / rows;
        s->showNormal();
        s->setGeometry(x0, y0, x1 - x0, y1 - y0);
        ++k;
    }
}

bool DesignerWindow::eventFilter(QObject* o, QEvent* e)
{
    if (o == mdi_->viewport() && e->type() == QEvent::Resize && grid_n_ > 0) arrangeGrid();
    return QMainWindow::eventFilter(o, e);
}

void DesignerWindow::submit()
{
    status_->setText("evaluating ...");
    engine_.submit(s_.design);
}

void DesignerWindow::onResult()
{
    s_.result = engine_.latest();
    refreshViews();
    updateIssues();
    if (!s_.result) return;
    QStringList maps;
    for (std::size_t k = 0; k < s_.result->layers.size() && k < s_.design.layers.size(); ++k)
        if (s_.result->layers[k].error.empty())
            maps << QString("%1: %2").arg(QString::fromStdString(s_.design.layers[k].name),
                                          s_.result->layers[k].surface.cylindrical() ? "cylindrical wrap" : "surface of revolution");
    status_->setText(QString("%1 · geometry %2 ms · %3 blade rows × %4 points · %5")
                         .arg(s_.result->ok() ? "OK" : "ERRORS")
                         .arg(1000.0 * s_.result->seconds, 0, 'f', 1)
                         .arg(s_.result->blade.nj)
                         .arg(s_.result->blade.ni)
                         .arg(maps.join(" · ")));
}

void DesignerWindow::refreshViews(bool refit)
{
    for (auto* v : views()) v->refresh(refit);
}

void DesignerWindow::updateIssues()
{
    issues_->clear();
    const auto R = s_.result;
    if (!R) return;
    for (const auto& e : R->errors) {
        auto* it = new QTreeWidgetItem(issues_, {QString::fromStdString(e), "ERROR", "design"});
        it->setForeground(0, QColor(200, 0, 0));
        it->setData(0, Qt::UserRole, -1);
    }
    for (std::size_t k = 0; k < R->layers.size() && k < s_.design.layers.size(); ++k) {
        const auto& L = R->layers[k];
        const QString name = QString::fromStdString(s_.design.layers[k].name);
        auto* top = new QTreeWidgetItem(issues_);
        top->setData(0, Qt::UserRole, static_cast<int>(k));
        if (!L.error.empty()) {
            top->setText(0, name + ": " + QString::fromStdString(L.error));
            top->setText(1, "ERROR");
            top->setForeground(0, QColor(200, 0, 0));
            continue;
        }
        const auto& q = L.section->quality;
        const int ne = q.count(Profile::Severity::Error), nw = q.count(Profile::Severity::Warning);
        top->setText(0, QString("%1 — %2").arg(name, QString::fromStdString(L.section->designation)));
        top->setText(1, ne ? QString("%1 error(s)").arg(ne) : nw ? QString("%1 warning(s)").arg(nw) : "OK");
        top->setForeground(1, ne ? QColor(200, 0, 0) : nw ? QColor(210, 120, 0) : QColor(0, 140, 0));
        for (const auto& is : q.issues) {
            auto* it = new QTreeWidgetItem(top, {QString::fromStdString(is.message), QString::fromStdString(Profile::to_string(is.severity)),
                                                 QString("%1 s = %2").arg(is.side == 0 ? "SS" : is.side == 1 ? "PS" : "profile").arg(is.s, 0, 'f', 4)});
            it->setData(0, Qt::UserRole, static_cast<int>(k));
            const QColor c = is.severity == Profile::Severity::Error ? QColor(200, 0, 0) : is.severity == Profile::Severity::Warning ? QColor(210, 120, 0) : QColor(90, 90, 90);
            it->setForeground(0, c);
            it->setForeground(1, c);
        }
        if (ne || nw) top->setExpanded(true);
    }
}

void DesignerWindow::settle()
{
    engine_.wait();
    QApplication::processEvents();
    onResult();
    QApplication::processEvents();
}

void DesignerWindow::updateTitle()
{
    const QString f = s_.path.isEmpty() ? "untitled" : QFileInfo(s_.path).fileName();
    setWindowTitle(QString("BladeGen Designer %1 — %2%3").arg(BLADEGEN_VERSION, f, dirty_ ? " *" : ""));
}

void DesignerWindow::newDesign(bool naca)
{
    s_.design = naca ? Design::BladeDesign::naca65_defaults() : Design::BladeDesign::defaults();
    s_.path.clear();
    s_.active_layer = 2;
    dirty_ = false;
    panel_->sync();
    updateTitle();
    submit();
    refreshViews(true);
}

bool DesignerWindow::openDesign(const QString& path, bool quiet)
{
    try {
        s_.design = Design::read_design(path.toStdString());
        s_.path = path;
        s_.active_layer = std::min(2, static_cast<int>(s_.design.layers.size()) - 1);
        dirty_ = false;
        panel_->sync();
        updateTitle();
        submit();
        refreshViews(true);
        return true;
    } catch (const std::exception& e) {
        if (quiet) std::fprintf(stderr, "bladegen_designer: %s: %s\n", qPrintable(path), e.what());
        else QMessageBox::warning(this, "Open design", QString("%1:\n%2").arg(path, e.what()));
        return false;
    }
}

void DesignerWindow::openDialog()
{
    const QString f = QFileDialog::getOpenFileName(this, "Open design", {}, "BladeGen design (*.json)");
    if (!f.isEmpty()) openDesign(f);
}

bool DesignerWindow::save(bool as)
{
    QString f = s_.path;
    if (as || f.isEmpty()) f = QFileDialog::getSaveFileName(this, "Save design", s_.path.isEmpty() ? "design.json" : s_.path, "BladeGen design (*.json)");
    if (f.isEmpty()) return false;
    try {
        Design::write_design(s_.design, f.toStdString());
        s_.path = f;
        dirty_ = false;
        updateTitle();
        return true;
    } catch (const std::exception& e) {
        QMessageBox::warning(this, "Save design", e.what());
        return false;
    }
}

void DesignerWindow::exportStep()
{
    const QString f = QFileDialog::getSaveFileName(this, "Export blade solid", "blade.step", "STEP (*.step *.stp);;IGES (*.igs);;BREP (*.brep);;STL (*.stl)");
    if (f.isEmpty()) return;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    try {
        const auto shapes = Design::build_design_solid(s_.design, engine_.evaluator());
        Occt::write_shape(shapes.blade, f.toStdString(), Occt::Unit::Millimetre);
        QApplication::restoreOverrideCursor();
        QMessageBox::information(this, "Export", QString("Wrote %1\nsolid valid: %2, volume %3 mm³, skin interpolation error %4 mm")
                                                     .arg(f)
                                                     .arg(shapes.blade_valid ? "yes" : "NO")
                                                     .arg(shapes.blade_volume * 1e9, 0, 'g', 6)
                                                     .arg(shapes.skin.max_error * 1e3, 0, 'g', 3));
    } catch (const std::exception& e) {
        QApplication::restoreOverrideCursor();
        QMessageBox::warning(this, "Export", e.what());
    }
}

void DesignerWindow::exportObj()
{
    if (!s_.result) return;
    const QString f = QFileDialog::getSaveFileName(this, "Export mesh", "blade_mesh.obj", "Wavefront OBJ (*.obj)");
    if (f.isEmpty()) return;
    try {
        Design::write_obj(*s_.result, f.toStdString());
    } catch (const std::exception& e) {
        QMessageBox::warning(this, "Export", e.what());
    }
}

void DesignerWindow::exportCsv()
{
    if (!s_.result) return;
    const QString f = QFileDialog::getSaveFileName(this, "Export layer coordinates", "layers.csv", "CSV (*.csv)");
    if (f.isEmpty()) return;
    try {
        Design::write_layer_csv(s_.design, *s_.result, f.toStdString());
    } catch (const std::exception& e) {
        QMessageBox::warning(this, "Export", e.what());
    }
}

void DesignerWindow::exportImage()
{
    const QString f = QFileDialog::getSaveFileName(this, "Save window image", "designer.png", "PNG (*.png)");
    if (!f.isEmpty()) grab().save(f);
}

void DesignerWindow::meshPointsDialog()
{
    QDialog dlg(this);
    dlg.setWindowTitle("Number of points");
    auto* f = new QFormLayout(&dlg);
    auto* rows = new QSpinBox(&dlg);
    rows->setRange(2, 400);
    rows->setValue(s_.design.mesh.n_rows);
    auto* body = new QSpinBox(&dlg);
    body->setRange(8, 1000);
    body->setValue(s_.design.mesh.n_body);
    auto* edge = new QSpinBox(&dlg);
    edge->setRange(4, 400);
    edge->setValue(s_.design.mesh.n_edge);
    f->addRow("spanwise rows", rows);
    f->addRow("points per side body", body);
    f->addRow("points per half edge", edge);
    auto* bb = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    f->addRow(bb);
    connect(bb, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(bb, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);
    if (dlg.exec() != QDialog::Accepted) return;
    s_.design.mesh.n_rows = rows->value();
    s_.design.mesh.n_body = body->value();
    s_.design.mesh.n_edge = edge->value();
    dirty_ = true;
    panel_->sync();
    submit();
}

void DesignerWindow::startFit()
{
    if (fit_watch_->isRunning()) return;
    const int k = s_.active_layer;
    panel_->setFitStatus("fitting the custom parameters to the NACA 65 section (a few seconds) ...", true);
    auto layer = s_.design.layers[static_cast<std::size_t>(k)];
    fit_watch_->setFuture(QtConcurrent::run([this, k, layer]() mutable -> QString {
        try {
            Design::fit_layer(layer, engine_.evaluator());
            QMetaObject::invokeMethod(this, [this, k, layer] {
                if (k < static_cast<int>(s_.design.layers.size())) {
                    auto& dst = s_.design.layers[static_cast<std::size_t>(k)];
                    dst.custom = layer.custom;
                    dst.fitted_to = layer.fitted_to;
                    dst.fit_rms = layer.fit_rms;
                    dst.method = Design::Method::Naca65Custom;
                    dirty_ = true;
                }
            }, Qt::BlockingQueuedConnection);
            return QString("fitted: RMS %1 % chord").arg(100.0 * layer.fit_rms, 0, 'f', 4);
        } catch (const std::exception& e) {
            return QString("fit failed: %1").arg(e.what());
        }
    }));
}

}  // namespace Designer
