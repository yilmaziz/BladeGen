// The pre-processing pages of the TURBOlenza window (see TurbolenzaWindow.hpp).
#include "TurbolenzaWindow.hpp"

#include "ParameterPanel.hpp"
#include "PlotWidget.hpp"

#include "turbolenza/Cascade.hpp"
#include "turbolenza/Grid.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QDoubleValidator>
#include <QFont>
#include <QFormLayout>
#include <QGroupBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QStackedWidget>
#include <QStandardItemModel>
#include <QTabWidget>
#include <QTableWidget>
#include <QTimer>
#include <QtConcurrent/QtConcurrentRun>
#include <QVBoxLayout>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <map>
#include <numbers>

namespace Designer {

using namespace BladeGen;
namespace T = BladeGen::Turbolenza;

namespace {

constexpr double kDeg = 180.0 / std::numbers::pi;

QDoubleSpinBox* dspin(QWidget* parent, double lo, double hi, int decimals, double step, const QString& suffix = {})
{
    auto* s = new QDoubleSpinBox(parent);
    s->setRange(lo, hi);
    s->setDecimals(decimals);
    s->setSingleStep(step);
    s->setKeyboardTracking(false);
    if (!suffix.isEmpty()) s->setSuffix(" " + suffix);
    s->setMinimumWidth(120);
    return s;
}

QLineEdit* sci(QWidget* parent, double v)
{
    auto* e = new QLineEdit(QString::number(v, 'g', 6), parent);
    auto* val = new QDoubleValidator(0.0, 1.0, 12, e);
    val->setNotation(QDoubleValidator::ScientificNotation);
    e->setValidator(val);
    e->setMaximumWidth(120);
    return e;
}

double sciValue(const QLineEdit* e)
{
    bool ok = false;
    const double v = e->text().toDouble(&ok);
    return ok ? v : std::numeric_limits<double>::quiet_NaN();  // refused by build_case, never repaired
}

QWidget* scrolled(QWidget* inner)
{
    auto* sa = new QScrollArea;
    sa->setWidget(inner);
    sa->setWidgetResizable(true);
    sa->setFrameShape(QFrame::NoFrame);
    return sa;
}

QLabel* help(const QString& t, QWidget* parent)
{
    auto* l = new QLabel(t, parent);
    l->setWordWrap(true);
    l->setStyleSheet("color: #444;");
    return l;
}

}  // namespace

// ── Blade ────────────────────────────────────────────────────────────────────

QWidget* TurbolenzaWindow::bladePage()
{
    auto* page = new QSplitter(Qt::Horizontal);
    auto* left = new QWidget;
    auto* v = new QVBoxLayout(left);
    v->addWidget(help("The blade section TURBOlenza solves: a COPY of one layer of the open BladeGen design. Edits here "
                      "do not change the 3-D design. The section is chord-normalized (LE (0, 0), TE (1, 0), suction "
                      "side +y); each side becomes the solver's clamped quintic B-spline.", left));
    auto* row = new QHBoxLayout;
    row->addWidget(new QLabel("Layer", left));
    layer_combo_ = new QComboBox(left);
    for (const auto& l : s_.design.layers) layer_combo_->addItem(QString::fromStdString(l.name));
    layer_combo_->setCurrentIndex(std::min(2, layer_combo_->count() - 1));
    row->addWidget(layer_combo_, 1);
    auto* take = new QPushButton("Take from layer", left);
    take->setToolTip("copy this layer of the open design (profile, stagger, 2-D solidity) into TURBOlenza");
    row->addWidget(take);
    v->addLayout(row);
    connect(take, &QPushButton::clicked, this, [this] { takeLayer(layer_combo_->currentIndex()); });

    auto* mrow = new QHBoxLayout;
    mrow->addWidget(new QLabel("construction:", left));
    method_ = new QComboBox(left);
    method_->addItems({"NACA 65", "Custom", "NACA 65 + Custom"});
    mrow->addWidget(method_, 1);
    v->addLayout(mrow);
    connect(method_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (syncing_ || i < 0) return;
        layer_.method = static_cast<Design::Method>(i);
        bladeChanged();
    });

    naca_box_ = new QGroupBox("NACA 65-(Cl0 A10) t/c", left);
    auto* nl = new QVBoxLayout(naca_box_);
    cl0_ = new Knob("Cl0 (A10 mean line)", 0.0, 2.5, 3, naca_box_);
    tc_ = new Knob("t_max / c", 0.02, 0.20, 4, naca_box_);
    nl->addWidget(cl0_);
    nl->addWidget(tc_);
    auto* vr = new QHBoxLayout;
    vr->addWidget(new QLabel("variant:", naca_box_));
    variant_ = new QComboBox(naca_box_);
    variant_->addItems({"kovach_sandercock", "emery_base"});
    vr->addWidget(variant_, 1);
    nl->addLayout(vr);
    fit_button_ = new QPushButton("Fit custom parameters to this NACA 65 section", naca_box_);
    nl->addWidget(fit_button_);
    v->addWidget(naca_box_);
    connect(cl0_, &Knob::changed, this, [this](double x) {
        layer_.naca.cl0 = x;
        bladeChanged();
    });
    connect(tc_, &Knob::changed, this, [this](double x) {
        layer_.naca.t_max_over_c = x;
        bladeChanged();
    });
    connect(variant_, &QComboBox::currentTextChanged, this, [this](const QString& t) {
        if (syncing_) return;
        layer_.naca.variant = Naca::variant_from_string(t.toStdString());
        bladeChanged();
    });
    connect(fit_button_, &QPushButton::clicked, this, [this] {
        QApplication::setOverrideCursor(Qt::WaitCursor);
        try {
            Design::fit_layer(layer_, ev_);
        } catch (const std::exception& e) {
            status_->setText(QString("fit failed: %1").arg(e.what()));
        }
        QApplication::restoreOverrideCursor();
        syncBlade();
        bladeChanged();
    });

    custom_box_ = new QGroupBox("Custom profile parameters", left);
    auto* cl = new QVBoxLayout(custom_box_);
    auto* tabs = new QTabWidget(custom_box_);
    std::map<std::string, QVBoxLayout*> fam;
    const std::map<std::string, QString> short_name = {{"Leading edge", "LE"}, {"Trailing edge", "TE"}, {"Suction side", "SS"},
                                                       {"Thickness distribution", "Thickness"}, {"Pressure side", "PS"}};
    for (const auto& info : Profile::param_table()) {
        if (!fam.count(info.family)) {
            auto* p = new QWidget(tabs);
            auto* pl = new QVBoxLayout(p);
            fam[info.family] = pl;
            const auto it = short_name.find(info.family);
            tabs->addTab(p, it != short_name.end() ? it->second : QString::fromStdString(info.family));
        }
        auto* k = new Knob(QString::fromStdString(info.label), info.lo, info.hi, 4, custom_box_, QString::fromStdString(info.name));
        fam[info.family]->addWidget(k);
        custom_.emplace_back(info.name, k);
        connect(k, &Knob::changed, this, [this, name = info.name](double x) {
            Profile::set_param(layer_.custom, name, x);
            bladeChanged();
        });
    }
    if (fam.count("Pressure side"))
        for (int i = 0; i < 3; ++i) {
            auto* k = new Knob(QString("free PS point %1: m'").arg(i + 1), 0.02, 0.98, 4, custom_box_, "ps.mprime_ctrl");
            fam["Pressure side"]->addWidget(k);
            ps_ctrl_.push_back(k);
            connect(k, &Knob::changed, this, [this, i](double x) {
                auto& c = layer_.custom.ps.mprime_ctrl;
                if (static_cast<std::size_t>(i) < c.size()) c[static_cast<std::size_t>(i)] = x;
                bladeChanged();
            });
        }
    auto* conn = new QWidget(tabs);
    auto* cn = new QVBoxLayout(conn);
    connector_ = new QComboBox(conn);
    connector_->addItems({"Fair", "C2Parametric", "G2Legacy"});
    cn->addWidget(new QLabel("edge connectors", conn));
    cn->addWidget(connector_);
    enforce_c1_ = new QCheckBox("thickness: derive shape_te so |dt/dξ| matches at the peak", conn);
    cn->addWidget(enforce_c1_);
    cn->addStretch(1);
    tabs->addTab(conn, "Connector");
    for (auto& [f, lay] : fam) lay->addStretch(1);
    connect(connector_, &QComboBox::currentTextChanged, this, [this](const QString& t) {
        if (syncing_) return;
        layer_.custom.connector = Profile::connector_mode_from_string(t.toStdString());
        bladeChanged();
    });
    connect(enforce_c1_, &QCheckBox::toggled, this, [this](bool on) {
        if (syncing_) return;
        layer_.custom.thickness.enforce_c1 = on;
        bladeChanged();
    });
    cl->addWidget(tabs);
    v->addWidget(custom_box_);
    blade_info_ = new QLabel(left);
    blade_info_->setWordWrap(true);
    blade_info_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    v->addWidget(blade_info_);
    v->addStretch(1);
    page->addWidget(scrolled(left));
    profile_plot_ = new PlotWidget(page);
    page->addWidget(profile_plot_);
    page->setStretchFactor(0, 2);
    page->setStretchFactor(1, 3);
    profile_timer_ = new QTimer(this);
    profile_timer_->setSingleShot(true);
    profile_timer_->setInterval(120);
    connect(profile_timer_, &QTimer::timeout, this, [this] { refreshProfile(); });
    return page;
}

void TurbolenzaWindow::syncBlade()
{
    syncing_ = true;
    method_->setCurrentIndex(static_cast<int>(layer_.method));
    cl0_->setValue(layer_.naca.cl0);
    tc_->setValue(layer_.naca.t_max_over_c);
    variant_->setCurrentText(QString::fromStdString(Naca::to_string(layer_.naca.variant)));
    for (auto& [name, k] : custom_) k->setValue(Profile::get_param(layer_.custom, name));
    for (std::size_t i = 0; i < ps_ctrl_.size(); ++i) {
        ps_ctrl_[i]->setEnabled(i < layer_.custom.ps.mprime_ctrl.size());
        if (i < layer_.custom.ps.mprime_ctrl.size()) ps_ctrl_[i]->setValue(layer_.custom.ps.mprime_ctrl[i]);
    }
    connector_->setCurrentText(QString::fromStdString(Profile::to_string(layer_.custom.connector)));
    enforce_c1_->setChecked(layer_.custom.thickness.enforce_c1);
    updateBladeEnabled();
    syncing_ = false;
}

void TurbolenzaWindow::updateBladeEnabled()
{
    naca_box_->setEnabled(layer_.method != Design::Method::Custom);
    fit_button_->setEnabled(layer_.method == Design::Method::Naca65Custom);
    custom_box_->setEnabled(layer_.uses_custom());
}

void TurbolenzaWindow::bladeChanged()
{
    if (syncing_) return;
    section_.reset();
    section_error_.clear();
    cascade_.reset();
    mesh_case_.reset();
    updateBladeEnabled();
    profile_timer_->start();
    cascade_timer_->start();
}

void TurbolenzaWindow::refreshProfile()
{
    try {
        const auto sec = ev_.section(layer_, Design::MeshOptions{});
        const auto& cam = sec->camber;
        Tl::drawProfile(profile_plot_, *sec, QString("Layer %1 (copy): %2").arg(QString::fromStdString(layer_.name), QString::fromStdString(sec->designation)));
        profile_plot_->setBadge({});
        QString t = QString("<b>%1</b><br>LE radius %2 c, TE radius %3 c, t_max %4 c at x/c %5<br>camber line angles: LE %6°, TE %7°")
                        .arg(QString::fromStdString(sec->designation))
                        .arg(sec->le_radius, 0, 'g', 5)
                        .arg(sec->te_radius, 0, 'g', 5)
                        .arg(cam.max_thickness, 0, 'f', 5)
                        .arg(cam.max_thickness_x, 0, 'f', 4)
                        .arg(cam.angle_deg.front(), 0, 'f', 3)
                        .arg(-cam.angle_deg.back(), 0, 'f', 3);
        if (sec->equivalent_arc_deg)
            t += QString("<br>equivalent circular arc (Aungier Eq. 4-6, TURBOlenza's χ): χ_LE = χ_TE = %1°").arg((*sec->equivalent_arc_deg)[0], 0, 'f', 4);
        else
            t += "<br>TURBOlenza's χ_LE, χ_TE: the camber line angles above";
        if (section_)
            t += QString("<br>solver curves: SS %1 / PS %2 control points, max fit residual %3 / %4 c")
                     .arg(fit_.suction.control_points)
                     .arg(fit_.pressure.control_points)
                     .arg(fit_.suction.max_residual, 0, 'e', 2)
                     .arg(fit_.pressure.max_residual, 0, 'e', 2);
        blade_info_->setText(t);
    } catch (const std::exception& e) {
        profile_plot_->setSeries({});
        profile_plot_->setBadge({QString("the profile cannot be built: %1").arg(e.what()), 2});
        profile_plot_->commit(true);
        blade_info_->setText(QString("<span style='color:#b00'>%1</span>").arg(e.what()));
    }
}

// ── Cascade ──────────────────────────────────────────────────────────────────

QWidget* TurbolenzaWindow::cascadePage()
{
    const T::CaseSettings d;
    auto* page = new QSplitter(Qt::Horizontal);
    auto* left = new QWidget;
    auto* v = new QVBoxLayout(left);
    v->addWidget(help("The cascade of the blade section and the far-inlet state. The pitch, the angle of attack and the "
                      "blade angles are derived, never entered.", left));
    auto* g = new QGroupBox("Cascade", left);
    auto* f = new QFormLayout(g);
    chord_ = dspin(g, 1e-4, 100.0, 5, 0.01, "m");
    solidity_ = dspin(g, 0.05, 10.0, 4, 0.01);
    stagger_ = dspin(g, -85.0, 85.0, 3, 0.5, "deg");
    f->addRow("chord c (true, the solver's c_ref)", chord_);
    f->addRow("solidity σ = c / s", solidity_);
    f->addRow("stagger γ (from axial)", stagger_);
    v->addWidget(g);
    auto* fl = new QGroupBox("Inlet flow and gas", left);
    auto* ff = new QFormLayout(fl);
    beta1_ = dspin(fl, -85.0, 85.0, 3, 0.5, "deg");
    mach_given_ = new QCheckBox("inlet Mach number M1", fl);
    mach1_ = dspin(fl, 0.0, 0.999, 4, 0.01);
    mdot_given_ = new QCheckBox("mass flow override ṁ' [kg/(m s)]", fl);
    mdot_ = dspin(fl, 0.0, 1e7, 6, 1.0);
    pt_ = dspin(fl, 1.0, 1e8, 2, 100.0, "Pa");
    tt_ = dspin(fl, 1.0, 5000.0, 3, 1.0, "K");
    gamma_ = dspin(fl, 1.01, 2.0, 4, 0.01);
    gas_r_ = dspin(fl, 1.0, 1e5, 3, 1.0, "J/(kg K)");
    omega_ = dspin(fl, -1e5, 1e5, 3, 1.0, "rad/s");
    omega_->setToolTip("only 0 (a stationary cascade) is solved: a rotating row is refused");
    ff->addRow("inlet flow angle β1 (from axial, + toward +y)", beta1_);
    ff->addRow(mach_given_, mach1_);
    ff->addRow(mdot_given_, mdot_);
    ff->addRow("inlet total pressure Pt", pt_);
    ff->addRow("inlet total temperature Tt", tt_);
    ff->addRow("ratio of specific heats γ", gamma_);
    ff->addRow("gas constant R", gas_r_);
    ff->addRow("rotation ω", omega_);
    v->addWidget(fl);
    auto* rg = new QGroupBox("Reynolds number (sets Pt; the tunnel holds Pt and Tt)", left);
    auto* rf = new QFormLayout(rg);
    re_ref_ = dspin(rg, 1e3, 1e10, 0, 10000.0);
    re_ref_->setValue(550000.0);
    m_ref_ = dspin(rg, 0.01, 0.99, 4, 0.01);
    m_ref_->setValue(0.30);
    auto* re_button = new QPushButton("Set P0 from a reference Reynolds number", rg);
    re_button->setToolTip("Re = rho1 W1 c / mu(T1) at the reference Mach number (chord-based, the inlet state, Sutherland at "
                          "the inlet static temperature, Tt as given): solved for Pt");
    re_info_ = new QLabel(rg);
    re_info_->setWordWrap(true);
    re_info_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    rf->addRow("reference Re (chord)", re_ref_);
    rf->addRow("at the Mach number", m_ref_);
    rf->addRow(re_button);
    rf->addRow(re_info_);
    v->addWidget(rg);
    connect(re_button, &QPushButton::clicked, this, [this] { setP0FromReynolds(); });
    derived_ = new QLabel(left);
    derived_->setWordWrap(true);
    derived_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    v->addWidget(derived_);
    v->addStretch(1);
    // defaults: CaseSettings (v1 run_full_solver.py)
    chord_->setValue(d.chord);
    solidity_->setValue(d.solidity);
    stagger_->setValue(d.stagger_deg);
    beta1_->setValue(d.inlet_flow_angle_deg);
    mach_given_->setChecked(d.inlet_mach.has_value());
    mach1_->setValue(d.inlet_mach.value_or(0.5));
    mdot_given_->setChecked(d.mass_flow_per_depth_override.has_value());
    mdot_->setValue(d.mass_flow_per_depth_override.value_or(0.0));
    pt_->setValue(d.inlet_total_pressure);
    tt_->setValue(d.inlet_total_temperature);
    gamma_->setValue(d.gamma);
    gas_r_->setValue(d.gas_constant);
    omega_->setValue(d.omega);
    mach1_->setEnabled(mach_given_->isChecked());
    mdot_->setEnabled(mdot_given_->isChecked());
    page->addWidget(scrolled(left));
    cascade_plot_ = new PlotWidget(page);
    page->addWidget(cascade_plot_);
    page->setStretchFactor(0, 2);
    page->setStretchFactor(1, 3);
    cascade_timer_ = new QTimer(this);
    cascade_timer_->setSingleShot(true);
    cascade_timer_->setInterval(150);
    connect(cascade_timer_, &QTimer::timeout, this, [this] {
        if (stack_->currentIndex() == 1) refreshCascade();
    });
    auto stale = [this] {
        cascade_.reset();
        mesh_case_.reset();
        cascade_timer_->start();
    };
    for (auto* s : {chord_, solidity_, stagger_, beta1_, mach1_, mdot_, pt_, tt_, gamma_, gas_r_, omega_})
        connect(s, &QDoubleSpinBox::valueChanged, this, stale);
    connect(mach_given_, &QCheckBox::toggled, this, [this, stale](bool on) {
        mach1_->setEnabled(on);
        stale();
    });
    connect(mdot_given_, &QCheckBox::toggled, this, [this, stale](bool on) {
        mdot_->setEnabled(on);
        stale();
    });
    return page;
}

void TurbolenzaWindow::setP0FromReynolds()
{
    try {
        const T::SutherlandLaw air;
        const double P0 = T::total_pressure_for_reynolds(re_ref_->value(), m_ref_->value(), tt_->value(), chord_->value(),
                                                         gamma_->value(), gas_r_->value(), air);
        pt_->setValue(P0);
        QString t = QString("Pt = %1 Pa gives Re = %2 at M = %3").arg(P0, 0, 'f', 2).arg(re_ref_->value(), 0, 'g', 7).arg(m_ref_->value());
        if (mach_given_->isChecked())
            t += QString("; at the case's M1 = %1: Re = %2")
                     .arg(mach1_->value())
                     .arg(T::chord_reynolds_number(mach1_->value(), pt_->value(), tt_->value(), chord_->value(), gamma_->value(),
                                                   gas_r_->value(), air),
                          0, 'g', 7);
        re_info_->setText(t);
    } catch (const std::exception& e) {
        re_info_->setText(QString("<span style='color:#b00'>%1</span>").arg(e.what()));
    }
}

void TurbolenzaWindow::setCascadeDefaults(double stagger_deg, double solidity)
{
    stagger_->setValue(stagger_deg);
    solidity_->setValue(solidity);
}

void TurbolenzaWindow::refreshCascade()
{
    const auto* sec = section();
    const double chord = chord_->value(), sol = solidity_->value(), b1 = beta1_->value(), gam = stagger_->value();
    QString t = QString("<b>derived</b><br>pitch s = c / σ = %1 m<br>angle of attack β1 - γ = %2°").arg(chord / sol, 0, 'g', 7).arg(b1 - gam, 0, 'f', 4);
    if (!sec) {
        derived_->setText(t + QString("<br><span style='color:#b00'>the section cannot be built: %1</span>").arg(section_error_));
        return;
    }
    try {
        if (!cascade_) cascade_ = T::make_cascade(*sec, sol, gam);
        const auto& c = *cascade_;
        t += QString("<br>χ_LE %1°, χ_TE %2° (camber %3°)<br>κ_LE = γ + χ_LE = %4°, κ_TE = γ - χ_TE = %5°"
                     "<br>incidence β1 - κ_LE = %6°<br>throat o/s = %7, o/c = %8")
                 .arg(c.chi_le * kDeg, 0, 'f', 4)
                 .arg(c.chi_te * kDeg, 0, 'f', 4)
                 .arg(c.camber_angle * kDeg, 0, 'f', 4)
                 .arg(c.kappa_le * kDeg, 0, 'f', 4)
                 .arg(c.kappa_te * kDeg, 0, 'f', 4)
                 .arg(b1 - c.kappa_le * kDeg, 0, 'f', 4)
                 .arg(c.throat_over_pitch, 0, 'f', 5)
                 .arg(c.throat_over_chord, 0, 'f', 5);
        derived_->setText(t);
        Tl::drawCascade(cascade_plot_, c, b1);
        cascade_plot_->setBadge({});
    } catch (const std::exception& e) {
        derived_->setText(t + QString("<br><span style='color:#b00'>%1</span>").arg(e.what()));
        cascade_plot_->setBadge({QString("the cascade cannot be built: %1").arg(e.what()), 2});
        cascade_plot_->update();
    }
}

// ── Mesh ─────────────────────────────────────────────────────────────────────

QWidget* TurbolenzaWindow::meshPage()
{
    const T::CaseSettings d;
    auto* page = new QSplitter(Qt::Horizontal);
    auto* left = new QWidget;
    auto* v = new QVBoxLayout(left);
    v->addWidget(help("The (m, η) grid of one passage (Aungier 2003 Secs. 5.2, 5.3): constant Δm and Δη; the first and "
                      "last blade stations lie inside the passage, one node on each surface (Fig. 5-5). A standard "
                      "level is chosen from the geometry; an explicit n_passage must be admissible.", left));
    auto* g = new QGroupBox("Grid", left);
    auto* f = new QFormLayout(g);
    level_ = new QComboBox(g);
    for (const char* l : T::kGridLevels) level_->addItem(l);
    level_->addItem("explicit (n_passage x n_eta)");
    level_->setCurrentText(QString::fromStdString(d.grid_level.value_or("medium")));
    n_passage_ = new QSpinBox(g);
    n_passage_->setRange(2, 5000);
    n_passage_->setValue(45);
    n_eta_ = new QSpinBox(g);
    n_eta_->setRange(3, 2000);
    n_eta_->setValue(21);
    ext_up_ = dspin(g, 0.01, 50.0, 3, 0.1, "chords");
    ext_down_ = dspin(g, 0.01, 50.0, 3, 0.1, "chords");
    ext_up_->setValue(d.upstream_extension_chords);
    ext_down_->setValue(d.downstream_extension_chords);
    edge_nodes_ = new QComboBox(g);
    edge_nodes_->addItem("inside the passage (Aungier, default)");
    edge_nodes_->addItem("on the LE and TE (touching, Fig. 5-5, experiment)");
    edge_nodes_->setCurrentIndex(d.edge_nodes == T::EdgeNodes::Touching ? 1 : 0);
    edge_nodes_->setToolTip("Inside: the first and last blade stations lie inside the passage, one node on each surface "
                            "(Aungier 2003 Fig. 5-5, recommended). Touching: they pass through the nose tip and the tail "
                            "end, both wall nodes on that point, the side lines ending there; Aungier reports that this "
                            "layout can cause local numerical instability.");
    f->addRow("grid level", level_);
    f->addRow("First / last blade nodes", edge_nodes_);
    f->addRow("n_passage (cells LE -> TE)", n_passage_);
    f->addRow("n_eta (nodes across)", n_eta_);
    f->addRow("upstream extension", ext_up_);
    f->addRow("downstream extension", ext_down_);
    v->addWidget(g);
    auto* gen = new QPushButton("Generate mesh", left);
    gen_button_ = gen;
    v->addWidget(gen);
    mesh_error_ = new QLabel(left);
    mesh_error_->setWordWrap(true);
    mesh_error_->setStyleSheet("color: #b00;");
    mesh_error_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    v->addWidget(mesh_error_);
    mesh_info_ = new QLabel(left);
    mesh_info_->setWordWrap(true);
    mesh_info_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    v->addWidget(mesh_info_);
    v->addStretch(1);
    auto explicitGrid = [this] {
        const bool ex = level_->currentIndex() == static_cast<int>(T::kGridLevels.size());
        n_passage_->setEnabled(ex);
        n_eta_->setEnabled(ex);
        mesh_case_.reset();
    };
    connect(level_, &QComboBox::currentIndexChanged, this, explicitGrid);
    connect(edge_nodes_, &QComboBox::currentIndexChanged, this, [this] { mesh_case_.reset(); });
    explicitGrid();
    for (auto* s : {n_passage_, n_eta_}) connect(s, &QSpinBox::valueChanged, this, [this] { mesh_case_.reset(); });
    for (auto* s : {ext_up_, ext_down_}) connect(s, &QDoubleSpinBox::valueChanged, this, [this] { mesh_case_.reset(); });
    connect(gen, &QPushButton::clicked, this, [this] { generateMesh(); });
    mesh_watch_ = new QFutureWatcher<MeshJob>(this);
    connect(mesh_watch_, &QFutureWatcher<MeshJob>::finished, this, [this] { meshFinished(); });
    page->addWidget(scrolled(left));
    grid_plot_ = new PlotWidget(page);
    page->addWidget(grid_plot_);
    page->setStretchFactor(0, 2);
    page->setStretchFactor(1, 3);
    return page;
}

bool TurbolenzaWindow::generateMesh(bool wait)
{
    if (mesh_watch_->isRunning()) return false;
    mesh_error_->clear();
    mesh_info_->clear();
    mesh_case_.reset();
    const auto* sec = section();
    if (!sec) {
        mesh_error_->setText("the section cannot be built: " + section_error_);
        return false;
    }
    // build_case (cascade + grid; a standard level searches its candidates, ~0.35 s each) runs off
    // the GUI thread
    gen_button_->setEnabled(false);
    mesh_info_->setText("generating the mesh...");
    mesh_watch_->setFuture(QtConcurrent::run([st = settings(), b = *sec] {
        MeshJob j;
        try {
            j.c = T::build_case(st, b);
        } catch (const T::GridGenerationError& e) {
            j.err = QString("GridGenerationError: %1").arg(e.what());
            for (const auto& [why, n] : e.reasons) j.err += QString("\n  %1: %2 candidate(s)").arg(QString::fromStdString(why)).arg(n);
        } catch (const T::GridTopologyError& e) {
            j.err = QString("GridTopologyError: %1").arg(e.what());
        } catch (const T::GridError& e) {
            j.err = QString("GridError: %1").arg(e.what());
        } catch (const T::CaseError& e) {
            j.err = QString("CaseError: %1").arg(e.what());
        } catch (const T::GeometryError& e) {
            j.err = QString("GeometryError: %1").arg(e.what());
        } catch (const std::exception& e) {
            j.err = QString("error: %1").arg(e.what());
        }
        return j;
    }));
    if (!wait) return true;
    while (mesh_watch_->isRunning() || !gen_button_->isEnabled()) QApplication::processEvents(QEventLoop::AllEvents, 50);
    return mesh_case_.has_value();
}

void TurbolenzaWindow::meshFinished()
{
    MeshJob j = mesh_watch_->result();
    gen_button_->setEnabled(true);
    mesh_info_->clear();
    if (!j.c) {
        mesh_error_->setText(j.err);
        grid_plot_->setSeries({});
        grid_plot_->setFills({});
        grid_plot_->setBadge({j.err, 2});
        grid_plot_->commit(true);
        return;
    }
    mesh_case_ = std::move(j.c);
    const auto& c = *mesh_case_;
    const auto& gr = c.grid;
    Tl::drawGrid(grid_plot_, gr);
    grid_plot_->setBadge({});
    // the edge stations: which surface each wall node lies on, and where (x/c on the blade)
    const auto& wi = gr.wall_identity;
    const auto& curves = c.cascade.active_curves;
    auto xc = [&](T::Surface s, double u) {
        if (!curves) return std::numeric_limits<double>::quiet_NaN();
        return (s == T::Surface::Suction ? curves->suction_unstaggered : curves->pressure_unstaggered).evaluate(u).m;
    };
    auto name = [](T::Surface s) { return s == T::Surface::Suction ? QString("SS") : QString("PS"); };
    QString edges;
    bool ok = !wi.stations.empty();
    const bool touching = gr.touching();
    if (touching && wi.stations.size() >= 3) {
        // the touching stations: both wall nodes on the nose tip / tail end (the point's own curve and u)
        for (int e = 0; e < 2; ++e) {
            const std::size_t b = e == 0 ? 0 : wi.stations.size() - 1;
            const auto s = gr.touching_surface[static_cast<std::size_t>(e)];
            edges += QString("<tr><td>%1 blade station i = %2</td><td colspan=2>both walls on the %3 (on %4 at x/c %5)</td>"
                             "<td><span style='color:#a60'>touching (experiment)</span></td></tr>")
                         .arg(e == 0 ? "first" : "last")
                         .arg(wi.stations[b])
                         .arg(e == 0 ? "nose tip" : "tail end")
                         .arg(name(s))
                         .arg(xc(s, gr.touching_u[static_cast<std::size_t>(e)]), 0, 'f', 5);
        }
    }
    // the stations whose walls must be one suction and one pressure node: the first and the last
    // blade station (inside), the stations next to the touching ones (touching)
    const std::size_t inner = touching && wi.stations.size() >= 3 ? 1 : 0;
    for (const std::size_t b : {inner, wi.stations.empty() ? std::size_t{0} : wi.stations.size() - 1 - inner}) {
        if (wi.stations.empty()) break;
        const bool one_each = wi.lower.surface[b] == T::Surface::Suction && wi.upper.surface[b] == T::Surface::Pressure;
        ok = ok && one_each;
        edges += QString("<tr><td>%1 blade station i = %2</td><td>lower wall on %3 at x/c %4</td><td>upper wall on %5 at x/c %6</td><td>%7</td></tr>")
                     .arg(touching ? (b == inner ? "next to the first" : "next to the last") : (b == 0 ? "first" : "last"))
                     .arg(wi.stations[b])
                     .arg(name(wi.lower.surface[b]))
                     .arg(xc(wi.lower.surface[b], wi.lower.u[b]), 0, 'f', 5)
                     .arg(name(wi.upper.surface[b]))
                     .arg(xc(wi.upper.surface[b], wi.upper.u[b]), 0, 'f', 5)
                     .arg(one_each ? "<span style='color:#080'>one node per surface</span>" : "<span style='color:#b00'>NOT one per surface</span>");
    }
    QString t = QString("<b>%1 x %2 nodes</b> (n_m x n_eta), Δm/c = %3, Δη = %4, blade stations %5 .. %6 (m_LE %7, m_TE %8)")
                    .arg(gr.n_m())
                    .arg(gr.n_eta())
                    .arg(gr.dm() / gr.chord(), 0, 'g', 6)
                    .arg(gr.deta(), 0, 'g', 6)
                    .arg(gr.first_blade_i)
                    .arg(gr.last_blade_i)
                    .arg(gr.m_le(), 0, 'f', 5)
                    .arg(gr.m_te(), 0, 'f', 5);
    t += "<table cellspacing=6>" + edges + "</table>";
    if (c.standard_info) {
        const auto& in = *c.standard_info;
        t += QString("level <b>%1</b>: %2 blade cells (target Δm/c %3), first SS / PS x/c %4 / %5, last %6 / %7, LE gap %8 c, "
                     "TE gap %9 c")
                 .arg(QString::fromStdString(in.level))
                 .arg(in.n_blade_cells)
                 .arg(in.target_dm_over_chord, 0, 'g', 5)
                 .arg(in.first_ss_x_c, 0, 'f', 5)
                 .arg(in.first_ps_x_c, 0, 'f', 5)
                 .arg(in.last_ss_x_c, 0, 'f', 5)
                 .arg(in.last_ps_x_c, 0, 'f', 5)
                 .arg(in.le_gap_over_chord, 0, 'g', 4)
                 .arg(in.te_gap_over_chord, 0, 'g', 4);
        t += QString(", cell aspect %1, actual extensions %2 / %3 chords").arg(in.cell_aspect, 0, 'g', 4).arg(in.upstream_extension_chords, 0, 'g', 5).arg(in.downstream_extension_chords, 0, 'g', 5);
        for (const auto& ch : in.checks) t += "<br>· " + QString::fromStdString(ch).toHtmlEscaped();
    }
    t += "<pre>" + QString::fromStdString(T::format_grid_summary(c)).toHtmlEscaped() + "</pre>";
    mesh_info_->setText(t);
    if (!ok)
        mesh_error_->setText(touching ? "a station next to a touching one does not cut one suction and one pressure surface"
                                      : "the first or last blade station does not cut one suction and one pressure surface");
    status_->setText(QString("mesh %1 x %2").arg(gr.n_m()).arg(gr.n_eta()));
}

// ── Solver ───────────────────────────────────────────────────────────────────

QWidget* TurbolenzaWindow::solverPage()
{
    const T::CaseSettings d;
    auto* inner = new QWidget;
    auto* v = new QVBoxLayout(inner);
    v->addWidget(help("Every solver input of the case (v1 scripts/run_full_solver.py; the defaults are its current values). "
                      "The coupled solve runs at most min(max density, max Kutta) passes.", inner));
    auto* rc = new QGroupBox("Run control", inner);
    auto* rf = new QFormLayout(rc);
    case_name_ = new QLineEdit(QString::fromStdString(d.case_name), rc);
    mode_ = new QComboBox(rc);
    mode_->addItems({"coupled_kutta (density + Kutta)", "fixed_exit (density only, fixed discharge angle)"});
    auto* fx = new QHBoxLayout;
    fixed_is_kappa_ = new QCheckBox("the blade TE angle κ_TE", rc);
    fixed_is_kappa_->setChecked(!d.fixed_exit_angle_deg.has_value());
    fixed_angle_ = dspin(rc, -85.0, 85.0, 4, 0.5, "deg");
    fixed_angle_->setValue(d.fixed_exit_angle_deg.value_or(20.0));
    fx->addWidget(fixed_is_kappa_);
    fx->addWidget(fixed_angle_);
    rf->addRow("case name", case_name_);
    rf->addRow("solver mode", mode_);
    rf->addRow("fixed exit angle", fx);
    v->addWidget(rc);
    auto* dg = new QGroupBox("Density iteration", inner);
    auto* df = new QFormLayout(dg);
    tol_rho_ = sci(dg, d.density_tolerance);
    auto* vw = new QHBoxLayout;
    vel_is_rho_ = new QCheckBox("= density tolerance", dg);
    vel_is_rho_->setChecked(!d.velocity_tolerance.has_value());
    tol_w_ = sci(dg, d.velocity_tolerance.value_or(d.density_tolerance));
    vw->addWidget(vel_is_rho_);
    vw->addWidget(tol_w_);
    max_rho_ = new QSpinBox(dg);
    max_rho_->setRange(1, 10000000);
    max_rho_->setValue(d.max_density_iterations);
    lambda_rho_ = dspin(dg, 0.001, 1.0, 4, 0.05);
    lambda_rho_->setValue(d.density_relaxation);
    stop_sonic_ = new QCheckBox("stop at the first pass with M_rel >= 1 (TransonicFlowError)", dg);
    stop_sonic_->setChecked(d.stop_at_sonic);
    df->addRow("density tolerance", tol_rho_);
    df->addRow("velocity tolerance", vw);
    df->addRow("max density passes", max_rho_);
    df->addRow("relaxation λ_ρ", lambda_rho_);
    df->addRow(stop_sonic_);
    v->addWidget(dg);
    auto* kg = new QGroupBox("Kutta condition", inner);
    auto* kf = new QFormLayout(kg);
    tol_k_ = sci(kg, d.kutta_tolerance);
    max_k_ = new QSpinBox(kg);
    max_k_->setRange(1, 10000000);
    max_k_->setValue(d.max_kutta_iterations);
    lambda_k_ = dspin(kg, 0.001, 1.0, 4, 0.05);
    lambda_k_->setValue(d.kutta_relaxation);
    auto* es = new QHBoxLayout;
    exit_start_ = new QComboBox(kg);
    exit_start_->addItems({"blade_te_angle", "manual"});
    exit_start_->setCurrentIndex(d.initial_exit_angle_mode == T::ExitAngleStart::Manual ? 1 : 0);
    exit_start_deg_ = dspin(kg, -85.0, 85.0, 4, 0.5, "deg");
    exit_start_deg_->setValue(d.initial_exit_angle_deg.value_or(20.0));
    es->addWidget(exit_start_);
    es->addWidget(exit_start_deg_);
    kf->addRow("Kutta tolerance ε_K", tol_k_);
    kf->addRow("max Kutta corrections", max_k_);
    kf->addRow("relaxation λ_K", lambda_k_);
    kf->addRow("initial exit angle", es);
    v->addWidget(kg);
    auto* cg = new QGroupBox("Compute", inner);
    auto* cf = new QFormLayout(cg);
    threads_ = new QSpinBox(cg);
    threads_->setRange(1, std::max(1, T::max_threads()));
    threads_->setValue(std::clamp(d.threads, 1, threads_->maximum()));
    threads_->setToolTip(QString("1 .. %1 (the machine's threads, at most %2)").arg(T::max_threads()).arg(T::kMaxThreads));
    linear_ = new QComboBox(cg);
    linear_->addItems({"Eigen SparseLU", "MKL PARDISO"});
    if (auto* m = qobject_cast<QStandardItemModel*>(linear_->model())) m->item(1)->setEnabled(T::pardiso_available());
    linear_->setItemData(1, T::pardiso_available() ? "parallel sparse direct solver (agrees with SparseLU to round-off)"
                                                    : "not in this build (Intel MKL not found)",
                         Qt::ToolTipRole);
    measure_ = new QCheckBox("measure the run time (wall, CPU, effective cores)", cg);
    measure_->setChecked(d.measure_runtime);
    cf->addRow("cores", threads_);
    cf->addRow("linear solver", linear_);
    cf->addRow(measure_);
    v->addWidget(cg);
    bl_box_ = new QGroupBox("Boundary layer (Aungier Sec. 5.6)", inner);
    auto* bf = new QFormLayout(bl_box_);
    viscous_ = new QCheckBox("compute the blade-surface boundary layer (not in v1)", bl_box_);
    viscous_->setChecked(d.viscous);
    bl_coupling_ = new QComboBox(bl_box_);
    bl_coupling_->addItems({"one_way (Aungier's analysis: the layer from the inviscid flow)",
                            "two_way (transpiration: the equivalent inviscid flow)"});
    bl_coupling_->setCurrentIndex(d.bl_coupling == T::BlCoupling::TwoWay ? 1 : 0);
    bl_cf_ = new QComboBox(bl_box_);
    bl_cf_->addItems({"centrifugal_2000: ½ c_f = μ/(ρWδ') (2 + Λ/6), Eq. 15-20",
                      "axial_2003_as_printed: ½ c_f = μ/(ρWδ') (1 + Λ/6), Eq. 5-123"});
    bl_cf_->setCurrentIndex(d.laminar_cf_form == T::LaminarCfForm::Axial2003AsPrinted ? 1 : 0);
    bl_wake_ = new QComboBox(bl_box_);
    bl_wake_->addItems({"constant (m_wake = m_SS + m_PS at the TE)", "marched (Head / Green with τ_w = 0)"});
    bl_wake_->setCurrentIndex(d.wake_model == T::WakeModel::Marched ? 1 : 0);
    bl_omega_ = dspin(bl_box_, 0.001, 1.0, 4, 0.05);
    bl_omega_->setValue(d.bl_relaxation);
    bl_tol_ = sci(bl_box_, d.bl_tolerance);
    bl_max_ = new QSpinBox(bl_box_);
    bl_max_->setRange(0, 100000);
    bl_max_->setValue(d.max_bl_iterations);
    bl_re_tr_ = dspin(bl_box_, 1.0, 1e6, 2, 10.0);
    bl_re_tr_->setValue(d.transition_re_theta);
    bl_hk_ = dspin(bl_box_, 1.01, 10.0, 3, 0.05);
    bl_hk_->setValue(d.separation_hk);
    bf->addRow(viscous_);
    bf->addRow("coupling", bl_coupling_);
    bf->addRow("laminar c_f form", bl_cf_);
    bf->addRow("wake model", bl_wake_);
    bf->addRow("relaxation ω_BL", bl_omega_);
    bf->addRow("tolerance (max change of m)", bl_tol_);
    bf->addRow("max BL iterations", bl_max_);
    bf->addRow("transition Re_θ (Eq. 5-126)", bl_re_tr_);
    bf->addRow("separation H_k (p. 111)", bl_hk_);
    v->addWidget(bl_box_);
    auto* br = new QHBoxLayout;
    run_ = new QPushButton("Run", inner);
    stop_ = new QPushButton("Stop", inner);
    stop_->setEnabled(false);
    br->addWidget(run_);
    br->addWidget(stop_);
    br->addStretch(1);
    v->addLayout(br);
    v->addStretch(1);
    connect(run_, &QPushButton::clicked, this, [this] { startRun(); });
    connect(stop_, &QPushButton::clicked, this, [this] { stopRun(); });
    connect(mode_, &QComboBox::currentIndexChanged, this, [this] { updateSolverEnabled(); });
    connect(fixed_is_kappa_, &QCheckBox::toggled, this, [this] { updateSolverEnabled(); });
    connect(vel_is_rho_, &QCheckBox::toggled, this, [this] { updateSolverEnabled(); });
    connect(exit_start_, &QComboBox::currentIndexChanged, this, [this] { updateSolverEnabled(); });
    connect(viscous_, &QCheckBox::toggled, this, [this] { updateSolverEnabled(); });
    connect(bl_coupling_, &QComboBox::currentIndexChanged, this, [this] { updateSolverEnabled(); });
    updateSolverEnabled();
    return scrolled(inner);
}

void TurbolenzaWindow::updateSolverEnabled()
{
    const bool fixed = mode_->currentIndex() == 1;
    fixed_is_kappa_->setEnabled(fixed);
    fixed_angle_->setEnabled(fixed && !fixed_is_kappa_->isChecked());
    tol_w_->setEnabled(!vel_is_rho_->isChecked());
    exit_start_deg_->setEnabled(exit_start_->currentIndex() == 1);
    const bool bl = viscous_->isChecked(), two = bl && bl_coupling_->currentIndex() == 1;
    for (QWidget* w : std::initializer_list<QWidget*>{bl_coupling_, bl_cf_, bl_re_tr_, bl_hk_}) w->setEnabled(bl);
    for (QWidget* w : std::initializer_list<QWidget*>{bl_wake_, bl_omega_, bl_tol_, bl_max_}) w->setEnabled(two);
}

T::CaseSettings TurbolenzaWindow::settings() const
{
    T::CaseSettings s;
    s.case_name = case_name_->text().toStdString();
    s.solver_mode = mode_->currentIndex() == 1 ? T::SolverMode::FixedExit : T::SolverMode::CoupledKutta;
    if (!fixed_is_kappa_->isChecked()) s.fixed_exit_angle_deg = fixed_angle_->value();
    s.gamma = gamma_->value();
    s.gas_constant = gas_r_->value();
    s.inlet_total_pressure = pt_->value();
    s.inlet_total_temperature = tt_->value();
    s.omega = omega_->value();
    s.chord = chord_->value();
    s.solidity = solidity_->value();
    s.stagger_deg = stagger_->value();
    s.inlet_mach = mach_given_->isChecked() ? std::optional<double>(mach1_->value()) : std::nullopt;
    s.inlet_flow_angle_deg = beta1_->value();
    s.mass_flow_per_depth_override = mdot_given_->isChecked() ? std::optional<double>(mdot_->value()) : std::nullopt;
    if (level_->currentIndex() < static_cast<int>(T::kGridLevels.size())) {
        s.grid_level = level_->currentText().toStdString();
    } else {
        s.grid_level.reset();
        s.n_passage = n_passage_->value();
        s.n_eta = n_eta_->value();
    }
    s.edge_nodes = edge_nodes_->currentIndex() == 1 ? T::EdgeNodes::Touching : T::EdgeNodes::Inside;
    s.upstream_extension_chords = ext_up_->value();
    s.downstream_extension_chords = ext_down_->value();
    s.density_tolerance = sciValue(tol_rho_);
    if (!vel_is_rho_->isChecked()) s.velocity_tolerance = sciValue(tol_w_);
    s.max_density_iterations = max_rho_->value();
    s.density_relaxation = lambda_rho_->value();
    s.stop_at_sonic = stop_sonic_->isChecked();
    s.kutta_tolerance = sciValue(tol_k_);
    s.max_kutta_iterations = max_k_->value();
    s.kutta_relaxation = lambda_k_->value();
    s.initial_exit_angle_mode = exit_start_->currentIndex() == 1 ? T::ExitAngleStart::Manual : T::ExitAngleStart::BladeTeAngle;
    if (exit_start_->currentIndex() == 1) s.initial_exit_angle_deg = exit_start_deg_->value();
    s.threads = threads_->value();
    s.linear_solver = linear_->currentIndex() == 1 ? T::LinearSolverKind::Pardiso : T::LinearSolverKind::SparseLU;
    s.measure_runtime = measure_->isChecked();
    s.viscous = viscous_->isChecked();
    s.bl_coupling = bl_coupling_->currentIndex() == 1 ? T::BlCoupling::TwoWay : T::BlCoupling::OneWay;
    s.laminar_cf_form = bl_cf_->currentIndex() == 1 ? T::LaminarCfForm::Axial2003AsPrinted : T::LaminarCfForm::Centrifugal2000;
    s.wake_model = bl_wake_->currentIndex() == 1 ? T::WakeModel::Marched : T::WakeModel::Constant;
    s.bl_relaxation = bl_omega_->value();
    s.bl_tolerance = sciValue(bl_tol_);
    s.max_bl_iterations = bl_max_->value();
    s.transition_re_theta = bl_re_tr_->value();
    s.separation_hk = bl_hk_->value();
    return s;
}

// ── Run and post-processing pages ────────────────────────────────────────────

QWidget* TurbolenzaWindow::runPage()
{
    auto* w = new QWidget;
    auto* v = new QVBoxLayout(w);
    auto* top = new QHBoxLayout;
    run_status_ = new QLabel("no run yet", w);
    run2_ = new QPushButton("Run", w);
    stop2_ = new QPushButton("Stop", w);
    stop2_->setEnabled(false);
    top->addWidget(run_status_, 1);
    top->addWidget(run2_);
    top->addWidget(stop2_);
    v->addLayout(top);
    connect(run2_, &QPushButton::clicked, this, [this] { startRun(); });
    connect(stop2_, &QPushButton::clicked, this, [this] { stopRun(); });
    auto* split = new QSplitter(Qt::Vertical, w);
    log_ = new QPlainTextEdit(split);
    log_->setReadOnly(true);
    log_->setLineWrapMode(QPlainTextEdit::NoWrap);
    QFont mono("Consolas");
    mono.setStyleHint(QFont::Monospace);
    log_->setFont(mono);
    split->addWidget(log_);
    auto* plotHost = new QWidget(split);
    run_plot_box_ = new QVBoxLayout(plotHost);
    run_plot_box_->setContentsMargins(0, 0, 0, 0);
    residuals_ = new PlotWidget(plotHost);
    run_plot_note_ = new QLabel("The residual history is shown in the Residuals window.", plotHost);
    run_plot_note_->hide();
    run_plot_box_->addWidget(run_plot_note_);
    run_plot_box_->addWidget(residuals_);
    split->addWidget(plotHost);
    split->setStretchFactor(0, 2);
    split->setStretchFactor(1, 3);
    v->addWidget(split, 1);
    Tl::drawResiduals(residuals_, {}, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                  std::numeric_limits<double>::quiet_NaN(), "no run yet");
    return w;
}

QWidget* TurbolenzaWindow::postPage()
{
    auto* w = new QWidget;
    auto* v = new QVBoxLayout(w);
    auto* row = new QHBoxLayout;
    auto* bres = new QPushButton("Residuals", w);
    auto* bcon = new QPushButton("Contours", w);
    auto* bload = new QPushButton("Blade loading", w);
    auto* bbl = new QPushButton("Boundary layer", w);
    auto* bsave = new QPushButton("Save results...", w);
    for (auto* b : {bres, bcon, bload, bbl, bsave}) row->addWidget(b);
    row->addStretch(1);
    v->addLayout(row);
    connect(bres, &QPushButton::clicked, this, [this] { showResiduals(); });
    connect(bcon, &QPushButton::clicked, this, [this] { showContours(); });
    connect(bload, &QPushButton::clicked, this, [this] { showLoading(); });
    connect(bbl, &QPushButton::clicked, this, [this] { showBoundaryLayer(); });
    connect(bsave, &QPushButton::clicked, this, [this] { saveDialog(); });
    report_ = new QLabel("no run yet", w);
    QFont mono("Consolas");
    mono.setStyleHint(QFont::Monospace);
    report_->setFont(mono);
    report_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    report_->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    v->addWidget(scrolled(report_), 1);
    return w;
}

// ── Mesh independence ────────────────────────────────────────────────────────

QWidget* TurbolenzaWindow::studyPage()
{
    auto* w = new QWidget;
    auto* v = new QVBoxLayout(w);
    v->addWidget(help("Mesh independence: the same case solved on each checked standard grid level, one after the other, "
                      "the coarsest first. Every other input is taken from the Blade, Cascade, Mesh (extensions and the "
                      "first / last blade node layout) and Solver pages; the Mesh page's grid level and explicit "
                      "n_passage x n_eta are not used. A level that fails is recorded and the study goes on.", w));
    auto* lg = new QGroupBox("Grid levels", w);
    auto* lh = new QHBoxLayout(lg);
    for (const char* name : T::kGridLevels) {
        auto* c = new QCheckBox(name, lg);
        c->setChecked(true);
        lh->addWidget(c);
        study_check_.push_back(c);
    }
    lh->addStretch(1);
    v->addWidget(lg);
    auto* br = new QHBoxLayout;
    study_run_ = new QPushButton("Run mesh study", w);
    study_stop_ = new QPushButton("Stop", w);
    study_stop_->setEnabled(false);
    study_save_ = new QPushButton("Save mesh study...", w);
    study_save_->setEnabled(false);
    study_status_ = new QLabel("no mesh study yet", w);
    br->addWidget(study_run_);
    br->addWidget(study_stop_);
    br->addWidget(study_save_);
    br->addWidget(study_status_, 1);
    v->addLayout(br);
    connect(study_run_, &QPushButton::clicked, this, [this] { startStudy(); });
    connect(study_stop_, &QPushButton::clicked, this, [this] { stopStudy(); });
    connect(study_save_, &QPushButton::clicked, this, [this] { saveStudyDialog(); });
    auto* split = new QSplitter(Qt::Vertical, w);
    study_table_ = new QTableWidget(0, 8, split);
    study_table_->setHorizontalHeaderLabels({"level", "n_m x n_eta", "dm/c", "passes", "status", "β_out [deg]", "max Mach", "wall time [s]"});
    study_table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    study_table_->verticalHeader()->hide();
    study_table_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    study_table_->horizontalHeader()->setSectionResizeMode(4, QHeaderView::Stretch);
    split->addWidget(study_table_);
    study_residuals_ = new PlotWidget(split);
    split->addWidget(study_residuals_);
    split->setStretchFactor(0, 2);
    split->setStretchFactor(1, 3);
    v->addWidget(split, 1);
    Tl::drawResiduals(study_residuals_, {}, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                      std::numeric_limits<double>::quiet_NaN(), "no mesh study yet");
    return w;
}

}  // namespace Designer
