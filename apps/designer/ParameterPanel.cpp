#include "ParameterPanel.hpp"

#include <QButtonGroup>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QStandardItemModel>
#include <QTabWidget>
#include <QToolBox>
#include <QVBoxLayout>

#include <cmath>
#include <map>

namespace Designer {

using namespace BladeGen;

Knob::Knob(const QString& label, double lo, double hi, int decimals, QWidget* parent, const QString& tip)
    : QWidget(parent), lo_(lo), hi_(hi)
{
    auto* l = new QHBoxLayout(this);
    l->setContentsMargins(0, 0, 0, 0);
    auto* lab = new QLabel(label, this);
    lab->setMinimumWidth(118);
    lab->setToolTip(tip.isEmpty() ? label : tip);
    slider_ = new QSlider(Qt::Horizontal, this);
    slider_->setRange(0, 1000);
    slider_->setMinimumWidth(90);
    spin_ = new QDoubleSpinBox(this);
    spin_->setRange(lo, hi);
    spin_->setDecimals(decimals);
    spin_->setSingleStep(std::pow(10.0, -decimals) * 10.0 > (hi - lo) / 50.0 ? (hi - lo) / 100.0 : std::pow(10.0, 1 - decimals));
    spin_->setKeyboardTracking(false);
    spin_->setMinimumWidth(78);
    l->addWidget(lab);
    l->addWidget(slider_, 1);
    l->addWidget(spin_);
    connect(slider_, &QSlider::valueChanged, this, [this](int v) {
        const double x = lo_ + (hi_ - lo_) * v / 1000.0;
        spin_->blockSignals(true);
        spin_->setValue(x);
        spin_->blockSignals(false);
        emit changed(spin_->value());
    });
    connect(spin_, &QDoubleSpinBox::valueChanged, this, [this](double v) {
        slider_->blockSignals(true);
        slider_->setValue(static_cast<int>(std::lround(1000.0 * (v - lo_) / (hi_ - lo_))));
        slider_->blockSignals(false);
        emit changed(v);
    });
}

void Knob::setValue(double v)
{
    spin_->blockSignals(true);
    slider_->blockSignals(true);
    if (v < lo_ || v > hi_) {  // values outside the suggested range widen it rather than being clamped
        lo_ = std::min(lo_, v);
        hi_ = std::max(hi_, v);
        spin_->setRange(lo_, hi_);
    }
    spin_->setValue(v);
    slider_->setValue(static_cast<int>(std::lround(1000.0 * (v - lo_) / (hi_ - lo_))));
    spin_->blockSignals(false);
    slider_->blockSignals(false);
}

double Knob::value() const { return spin_->value(); }

// ─────────────────────────────────────────────────────────────────────────────

namespace {

QDoubleSpinBox* coord(QWidget* parent)
{
    auto* s = new QDoubleSpinBox(parent);
    s->setRange(-1e5, 1e5);
    s->setDecimals(6);
    s->setKeyboardTracking(false);
    s->setMinimumWidth(110);
    return s;
}

}  // namespace

ParameterPanel::ParameterPanel(Session& s, QWidget* parent) : QWidget(parent), s_(s)
{
    auto* l = new QVBoxLayout(this);
    l->setContentsMargins(4, 4, 4, 4);
    auto* box = new QToolBox(this);
    box->addItem(meridionalPage(), "Initial Meridional Configuration");
    box->addItem(layerPage(), "Layers and Profiles");
    box->addItem(bladePage(), "Blade and Mesh");
    box->setCurrentIndex(1);
    l->addWidget(box);
    sync();
}

Design::Layer& ParameterPanel::layer()
{
    return s_.design.layers[static_cast<std::size_t>(std::clamp(s_.active_layer, 0, static_cast<int>(s_.design.layers.size()) - 1))];
}

void ParameterPanel::changed()
{
    if (syncing_) return;
    updateEnabled();
    emit designChanged();
}

QWidget* ParameterPanel::meridionalPage()
{
    auto* w = new QWidget(this);
    auto* v = new QVBoxLayout(w);
    auto* help = new QLabel("Axial compressor: hub and casing points at the LE and TE in cylindrical coordinates "
                            "(Z axial, R radius, mm).  Between LE and TE each wall follows Bell & Mehta's fifth-order "
                            "contour (zero slope and curvature at both ends, cylinders upstream / downstream) or a "
                            "straight line.", w);
    help->setWordWrap(true);
    v->addWidget(help);
    auto wall = [&](const QString& name, Wall& W, Blade::WallPoints Blade::MeridionalConfig::*pts, Blade::WallLaw Blade::MeridionalConfig::*law) {
        auto* g = new QGroupBox(name, w);
        auto* f = new QGridLayout(g);
        W.z_le = coord(g);
        W.r_le = coord(g);
        W.z_te = coord(g);
        W.r_te = coord(g);
        f->addWidget(new QLabel("LE   Z:"), 0, 0);
        f->addWidget(W.z_le, 0, 1);
        f->addWidget(new QLabel("R:"), 0, 2);
        f->addWidget(W.r_le, 0, 3);
        f->addWidget(new QLabel("TE   Z:"), 1, 0);
        f->addWidget(W.z_te, 1, 1);
        f->addWidget(new QLabel("R:"), 1, 2);
        f->addWidget(W.r_te, 1, 3);
        W.bell = new QPushButton("Bell && Mehta", g);
        W.linear = new QPushButton("Linear", g);
        W.bell->setCheckable(true);
        W.linear->setCheckable(true);
        W.bell->setToolTip("R(z) = R_LE + (R_TE - R_LE)(6ξ⁵ - 15ξ⁴ + 10ξ³), ξ = (z - z_LE)/(z_TE - z_LE): G2 onto cylinders");
        W.linear->setToolTip("straight line through the LE and TE points, continued straight");
        auto* grp = new QButtonGroup(g);
        grp->setExclusive(true);
        grp->addButton(W.bell);
        grp->addButton(W.linear);
        f->addWidget(W.bell, 2, 0, 1, 2);
        f->addWidget(W.linear, 2, 2, 1, 2);
        v->addWidget(g);
        for (auto* sb : {W.z_le, W.r_le, W.z_te, W.r_te})
            connect(sb, &QDoubleSpinBox::valueChanged, this, [this, &W, pts] {
                if (syncing_) return;
                auto& p = s_.design.meridional.*pts;
                p = {W.z_le->value(), W.r_le->value(), W.z_te->value(), W.r_te->value()};
                changed();
            });
        connect(W.bell, &QPushButton::toggled, this, [this, law](bool on) {
            if (syncing_ || !on) return;
            s_.design.meridional.*law = Blade::WallLaw::BellMehta;
            changed();
        });
        connect(W.linear, &QPushButton::toggled, this, [this, law](bool on) {
            if (syncing_ || !on) return;
            s_.design.meridional.*law = Blade::WallLaw::Linear;
            changed();
        });
    };
    wall("Casing (shroud)", casing_, &Blade::MeridionalConfig::casing, &Blade::MeridionalConfig::casing_law);
    wall("Hub", hub_, &Blade::MeridionalConfig::hub, &Blade::MeridionalConfig::hub_law);
    auto* g = new QGroupBox("Inlet / outlet extensions (× axial LE–TE length)", w);
    auto* f = new QFormLayout(g);
    ext_in_ = new QDoubleSpinBox(g);
    ext_out_ = new QDoubleSpinBox(g);
    for (auto* e : {ext_in_, ext_out_}) {
        e->setRange(0.05, 5.0);
        e->setSingleStep(0.05);
        e->setDecimals(3);
        e->setKeyboardTracking(false);
    }
    f->addRow("inlet", ext_in_);
    f->addRow("outlet", ext_out_);
    connect(ext_in_, &QDoubleSpinBox::valueChanged, this, [this](double x) { s_.design.meridional.inlet_extension = x; changed(); });
    connect(ext_out_, &QDoubleSpinBox::valueChanged, this, [this](double x) { s_.design.meridional.outlet_extension = x; changed(); });
    v->addWidget(g);
    // stacking of the sections (Vos et al. 2025, Sec. 4)
    auto* sg = new QGroupBox("Stacking of the sections", w);
    auto* sl = new QVBoxLayout(sg);
    auto* sh = new QLabel("Prescribed LE / TE lines: each section's LE and TE sit on the meridional LE / TE lines, so "
                          "the chord follows.  Reference point: each section's centre of gravity (or LE, TE) sits on "
                          "the stacking line, and the chord is free (chord scale on the layer page).  Variant 1: every "
                          "reference point at the same Z and θ, a radial line.  Variant 2: each one at the point of its "
                          "stream surface closest to the previous layer's, along the stream-surface normals.", sg);
    sh->setWordWrap(true);
    sl->addWidget(sh);
    stacking_mode_ = new QComboBox(sg);
    stacking_mode_->addItems({"prescribed LE / TE lines", "reference point · variant 1: same Z and θ (radial line)",
                              "reference point · variant 2: closest point to the previous layer"});
    sl->addWidget(stacking_mode_);
    auto* srow = new QHBoxLayout;
    srow->addWidget(new QLabel("reference point", sg));
    stacking_ref_ = new QComboBox(sg);
    stacking_ref_->addItems({"centre of gravity", "leading edge", "trailing edge"});
    srow->addWidget(stacking_ref_, 1);
    sl->addLayout(srow);
    auto* zrow = new QHBoxLayout;
    stacking_z_ = coord(sg);
    stacking_theta_ = coord(sg);
    stacking_theta_->setRange(-180.0, 180.0);
    stacking_z_->setToolTip("axial position of the first (hub) layer's reference point");
    stacking_theta_->setToolTip("angle of every reference point");
    zrow->addWidget(new QLabel("hub Z_ref [mm]", sg));
    zrow->addWidget(stacking_z_);
    zrow->addWidget(new QLabel("θ_ref [deg]", sg));
    zrow->addWidget(stacking_theta_);
    sl->addLayout(zrow);
    stacking_take_ = new QPushButton("Take Z_ref, θ_ref from the current hub section", sg);
    stacking_take_->setToolTip("the hub section's reference point stays where it is now");
    sl->addWidget(stacking_take_);
    v->addWidget(sg);
    connect(stacking_mode_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (syncing_ || i < 0) return;
        auto& st = s_.design.stacking;
        const bool was_edges = st.mode == Design::Stacking::Mode::Edges;
        st.mode = i == 0 ? Design::Stacking::Mode::Edges : Design::Stacking::Mode::ReferencePoint;
        if (i > 0) st.variant = i == 1 ? Design::Stacking::Variant::SameZTheta : Design::Stacking::Variant::ClosestPoint;
        if (was_edges && i > 0) takeStackingPoint();  // the hub section stays where it was
        changed();
    });
    connect(stacking_ref_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (syncing_ || i < 0) return;
        s_.design.stacking.reference = static_cast<Design::Stacking::Reference>(i);
        changed();
    });
    connect(stacking_z_, &QDoubleSpinBox::valueChanged, this, [this](double x) {
        if (syncing_) return;
        s_.design.stacking.z_ref = x;
        changed();
    });
    connect(stacking_theta_, &QDoubleSpinBox::valueChanged, this, [this](double x) {
        if (syncing_) return;
        s_.design.stacking.theta_ref_deg = x;
        changed();
    });
    connect(stacking_take_, &QPushButton::clicked, this, [this] {
        takeStackingPoint();
        changed();
    });
    // re-stagger of the whole blade
    auto* rg = new QGroupBox("Re-stagger the whole blade", w);
    auto* rl = new QVBoxLayout(rg);
    auto* rh = new QLabel("Turns every section by the same angle about a radial axis; every metal angle changes by "
                          "that angle.  + closes the blade (more stall margin), − opens it (more flow).  Rotors: "
                          "about the hub LE.", rg);
    rh->setWordWrap(true);
    rl->addWidget(rh);
    restagger_ = new Knob("re-stagger Δγ [deg]", -15, 15, 3, rg, "added to the stagger of every layer; the sections turn about the axis");
    rl->addWidget(restagger_);
    auto* ar = new QHBoxLayout;
    ar->addWidget(new QLabel("axis through", rg));
    restagger_axis_ = new QComboBox(rg);
    restagger_axis_->addItems({"hub section", "mid-span section", "tip section", "point (Z, θ)"});
    restagger_axis_->setToolTip("the radial line the blade turns about passes through this point");
    ar->addWidget(restagger_axis_, 1);
    rl->addLayout(ar);
    restagger_chord_ = new Knob("at chord x/c", 0.0, 1.0, 3, rg, "pivot on the reference section's chord line: 0 = LE, 1 = TE");
    rl->addWidget(restagger_chord_);
    auto* pr = new QHBoxLayout;
    restagger_z_ = coord(rg);
    restagger_theta_ = coord(rg);
    restagger_theta_->setRange(-180.0, 180.0);
    pr->addWidget(new QLabel("Z [mm]", rg));
    pr->addWidget(restagger_z_);
    pr->addWidget(new QLabel("θ [deg]", rg));
    pr->addWidget(restagger_theta_);
    rl->addLayout(pr);
    v->addWidget(rg);
    connect(restagger_, &Knob::changed, this, [this](double x) { s_.design.restagger.angle_deg = x; changed(); });
    connect(restagger_chord_, &Knob::changed, this, [this](double x) { s_.design.restagger.chord_fraction = x; changed(); });
    connect(restagger_axis_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (syncing_ || i < 0) return;
        s_.design.restagger.axis = static_cast<Design::Restagger::Axis>(i);
        changed();
    });
    connect(restagger_z_, &QDoubleSpinBox::valueChanged, this, [this](double x) {
        if (syncing_) return;
        s_.design.restagger.z = x;
        changed();
    });
    connect(restagger_theta_, &QDoubleSpinBox::valueChanged, this, [this](double x) {
        if (syncing_) return;
        s_.design.restagger.theta_deg = x;
        changed();
    });
    v->addStretch(1);
    return w;
}

QWidget* ParameterPanel::layerPage()
{
    auto* w = new QWidget(this);
    auto* v = new QVBoxLayout(w);
    // layer selector
    auto* row = new QHBoxLayout;
    layers_ = new QButtonGroup(w);
    layers_->setExclusive(true);
    for (int j = 0; j < static_cast<int>(s_.design.layers.size()); ++j) {
        auto* b = new QPushButton(QString::fromStdString(s_.design.layers[static_cast<std::size_t>(j)].name), w);
        b->setCheckable(true);
        layers_->addButton(b, j);
        row->addWidget(b);
    }
    connect(layers_, &QButtonGroup::idClicked, this, [this](int id) {
        s_.active_layer = id;
        sync();
        emit activeLayerChanged();
    });
    v->addLayout(row);
    // method
    auto* mrow = new QHBoxLayout;
    mrow->addWidget(new QLabel("construction:", w));
    method_ = new QComboBox(w);
    method_->addItems({"NACA 65", "Custom", "NACA 65 + Custom", "DCA (planned)", "C4 (planned)"});
    if (auto* model = qobject_cast<QStandardItemModel*>(method_->model()))
        for (int i : {3, 4}) model->item(i)->setEnabled(false);
    mrow->addWidget(method_, 1);
    v->addLayout(mrow);
    connect(method_, &QComboBox::currentIndexChanged, this, [this](int i) {
        if (syncing_ || i > 2) return;
        layer().method = static_cast<Design::Method>(i);
        changed();
    });
    // placement
    auto* pg = new QGroupBox("Placement on the stream surface", w);
    auto* pl = new QVBoxLayout(pg);
    stagger_ = new Knob("stagger γ [deg]", -70, 70, 3, pg, "chord angle from the m' axis toward +θ");
    theta_le_ = new Knob("θ LE [deg] (lean)", -30, 30, 3, pg);
    solidity_ = new Knob("solidity (2-D)", 0.3, 3.0, 3, pg, "2-D cascade solidity for exports; the 3-D solidity follows from the blade count");
    chord_scale_ = new Knob("chord scale (reference-point stacking)", 0.5, 2.0, 3, pg,
                            "c' = scale × the chord of the prescribed LE / TE lines; used when the sections are stacked on a reference point");
    for (auto* k : {stagger_, theta_le_, chord_scale_, solidity_}) pl->addWidget(k);
    connect(chord_scale_, &Knob::changed, this, [this](double x) { layer().chord_scale = x; changed(); });
    connect(stagger_, &Knob::changed, this, [this](double x) { layer().stagger_deg = x; changed(); });
    connect(theta_le_, &Knob::changed, this, [this](double x) { layer().theta_le_deg = x; changed(); });
    connect(solidity_, &Knob::changed, this, [this](double x) { layer().solidity = x; changed(); });
    v->addWidget(pg);
    // NACA 65
    naca_box_ = new QGroupBox("NACA 65-(Cl0 A10) t/c", w);
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
    fit_ = new QPushButton("Fit custom parameters to this NACA 65 section", naca_box_);
    fit_status_ = new QLabel(naca_box_);
    fit_status_->setWordWrap(true);
    nl->addWidget(fit_);
    nl->addWidget(fit_status_);
    connect(cl0_, &Knob::changed, this, [this](double x) { layer().naca.cl0 = x; changed(); });
    connect(tc_, &Knob::changed, this, [this](double x) { layer().naca.t_max_over_c = x; changed(); });
    connect(variant_, &QComboBox::currentTextChanged, this, [this](const QString& t) {
        if (syncing_) return;
        layer().naca.variant = Naca::variant_from_string(t.toStdString());
        changed();
    });
    connect(fit_, &QPushButton::clicked, this, [this] { emit fitRequested(); });
    v->addWidget(naca_box_);
    // custom
    custom_box_ = new QGroupBox("Custom profile parameters", w);
    auto* cl = new QVBoxLayout(custom_box_);
    custom_hint_ = new QLabel(custom_box_);
    custom_hint_->setWordWrap(true);
    cl->addWidget(custom_hint_);
    auto* tabs = new QTabWidget(custom_box_);
    std::map<std::string, QVBoxLayout*> fam;
    const std::map<std::string, QString> short_name = {{"Leading edge", "LE"}, {"Trailing edge", "TE"}, {"Suction side", "SS"},
                                                       {"Thickness distribution", "Thickness"}, {"Pressure side", "PS"}};
    for (const auto& info : Profile::param_table()) {
        if (!fam.count(info.family)) {
            auto* page = new QWidget(tabs);
            auto* pl2 = new QVBoxLayout(page);
            auto* head = new QLabel("<b>" + QString::fromStdString(info.family) + "</b>", page);
            pl2->addWidget(head);
            fam[info.family] = pl2;
            const auto it = short_name.find(info.family);
            tabs->addTab(page, it != short_name.end() ? it->second : QString::fromStdString(info.family));
        }
        auto* k = new Knob(QString::fromStdString(info.label), info.lo, info.hi, 4, custom_box_, QString::fromStdString(info.name));
        fam[info.family]->addWidget(k);
        custom_.emplace_back(info.name, k);
        connect(k, &Knob::changed, this, [this, name = info.name](double x) {
            Profile::set_param(layer().custom, name, x);
            changed();
        });
    }
    if (fam.count("Pressure side")) {
        for (int i = 0; i < 3; ++i) {
            auto* k = new Knob(QString("free PS point %1: m'").arg(i + 1), 0.02, 0.98, 4, custom_box_, "ps.mprime_ctrl");
            fam["Pressure side"]->addWidget(k);
            ps_ctrl_.push_back(k);
            connect(k, &Knob::changed, this, [this, i](double x) {
                auto& c = layer().custom.ps.mprime_ctrl;
                if (static_cast<std::size_t>(i) < c.size()) c[static_cast<std::size_t>(i)] = x;
                changed();
            });
        }
    }
    auto* conn = new QWidget(tabs);
    auto* cn = new QVBoxLayout(conn);
    cn->addWidget(new QLabel("<b>Edge connectors</b>", conn));
    connector_ = new QComboBox(conn);
    connector_->addItems({"Fair", "C2Parametric", "G2Legacy"});
    cn->addWidget(connector_);
    enforce_c1_ = new QCheckBox("thickness: derive shape_te so |dt/dξ| matches at the peak", conn);
    cn->addWidget(enforce_c1_);
    cn->addStretch(1);
    tabs->addTab(conn, "Connector");
    connect(connector_, &QComboBox::currentTextChanged, this, [this](const QString& t) {
        if (syncing_) return;
        layer().custom.connector = Profile::connector_mode_from_string(t.toStdString());
        changed();
    });
    connect(enforce_c1_, &QCheckBox::toggled, this, [this](bool on) {
        if (syncing_) return;
        layer().custom.thickness.enforce_c1 = on;
        changed();
    });
    for (auto& [f, lay] : fam) lay->addStretch(1);
    cl->addWidget(tabs);
    auto* copy = new QPushButton("Copy this layer's profile to all layers", custom_box_);
    connect(copy, &QPushButton::clicked, this, [this] {
        const auto src = layer();
        for (auto& l : s_.design.layers) {
            l.method = src.method;
            l.naca = src.naca;
            l.custom = src.custom;
            l.fitted_to = src.fitted_to;
            l.fit_rms = src.fit_rms;
        }
        changed();
    });
    v->addWidget(custom_box_);
    v->addWidget(copy);
    v->addStretch(1);
    return w;
}

QWidget* ParameterPanel::bladePage()
{
    auto* w = new QWidget(this);
    auto* f = new QFormLayout(w);
    n_blades_ = new QSpinBox(w);
    n_blades_->setRange(1, 500);
    mirror_ = new QCheckBox("suction side toward −θ", w);
    rows_ = new QSpinBox(w);
    rows_->setRange(2, 400);
    n_body_ = new QSpinBox(w);
    n_body_->setRange(8, 1000);
    n_edge_ = new QSpinBox(w);
    n_edge_->setRange(4, 400);
    clustering_ = new QDoubleSpinBox(w);
    clustering_->setRange(0.0, 10.0);
    clustering_->setSingleStep(0.25);
    for (auto* sb : {n_blades_, rows_, n_body_, n_edge_}) sb->setKeyboardTracking(false);
    clustering_->setKeyboardTracking(false);
    f->addRow("number of blades", n_blades_);
    f->addRow("mirror", mirror_);
    f->addRow(new QLabel("<b>3-D mesh</b>", w));
    f->addRow("spanwise rows", rows_);
    f->addRow("points per side body", n_body_);
    f->addRow("points per half edge", n_edge_);
    f->addRow("hub / tip clustering (tanh d)", clustering_);
    connect(n_blades_, &QSpinBox::valueChanged, this, [this](int n) { s_.design.n_blades = n; changed(); });
    connect(mirror_, &QCheckBox::toggled, this, [this](bool on) { s_.design.mirror = on; changed(); });
    connect(rows_, &QSpinBox::valueChanged, this, [this](int n) { s_.design.mesh.n_rows = n; changed(); });
    connect(n_body_, &QSpinBox::valueChanged, this, [this](int n) { s_.design.mesh.n_body = n; changed(); });
    connect(n_edge_, &QSpinBox::valueChanged, this, [this](int n) { s_.design.mesh.n_edge = n; changed(); });
    connect(clustering_, &QDoubleSpinBox::valueChanged, this, [this](double x) { s_.design.mesh.row_clustering = x; changed(); });
    return w;
}

void ParameterPanel::updateEnabled()
{
    const bool point = s_.design.restagger.axis == Design::Restagger::Axis::Point;
    restagger_chord_->setEnabled(!point);
    restagger_z_->setEnabled(point);
    restagger_theta_->setEnabled(point);
    const bool ref = s_.design.stacking.mode == Design::Stacking::Mode::ReferencePoint;
    for (QWidget* x : std::initializer_list<QWidget*>{stacking_ref_, stacking_z_, stacking_theta_, stacking_take_, chord_scale_}) x->setEnabled(ref);
    theta_le_->setEnabled(!ref);
    theta_le_->setToolTip(ref ? "not used: the sections are stacked on the reference point (θ_ref on the meridional page)" : QString());
    const auto& l = layer();
    naca_box_->setEnabled(l.method != Design::Method::Custom);
    fit_->setEnabled(l.method == Design::Method::Naca65Custom);
    custom_box_->setEnabled(l.uses_custom());
    if (l.method == Design::Method::Naca65) custom_hint_->setText("The layer is the NACA 65 section; choose Custom or NACA 65 + Custom to edit these.");
    else if (l.method == Design::Method::Naca65Custom && !l.fitted_to) custom_hint_->setText("Press Fit to start the custom parameters from the NACA 65 section.");
    else if (l.method == Design::Method::Naca65Custom)
        custom_hint_->setText(QString("Fitted to NACA 65-(%1 A10) %2 (%3), RMS %4 % chord; now editable.")
                                  .arg(l.fitted_to->cl0).arg(l.fitted_to->t_max_over_c)
                                  .arg(QString::fromStdString(Naca::to_string(l.fitted_to->variant))).arg(100.0 * l.fit_rms, 0, 'f', 4));
    else custom_hint_->setText("Custom profile: every knob by family.");
}

void ParameterPanel::sync()
{
    syncing_ = true;
    const auto& m = s_.design.meridional;
    auto setWall = [](Wall& W, const Blade::WallPoints& p, Blade::WallLaw law) {
        W.z_le->setValue(p.z_le);
        W.r_le->setValue(p.r_le);
        W.z_te->setValue(p.z_te);
        W.r_te->setValue(p.r_te);
        W.bell->setChecked(law == Blade::WallLaw::BellMehta);
        W.linear->setChecked(law == Blade::WallLaw::Linear);
    };
    setWall(hub_, m.hub, m.hub_law);
    setWall(casing_, m.casing, m.casing_law);
    ext_in_->setValue(m.inlet_extension);
    ext_out_->setValue(m.outlet_extension);
    const auto& rs = s_.design.restagger;
    restagger_->setValue(rs.angle_deg);
    restagger_axis_->setCurrentIndex(static_cast<int>(rs.axis));
    restagger_chord_->setValue(rs.chord_fraction);
    restagger_z_->setValue(rs.z);
    restagger_theta_->setValue(rs.theta_deg);
    const auto& st = s_.design.stacking;
    stacking_mode_->setCurrentIndex(st.mode == Design::Stacking::Mode::Edges ? 0 : st.variant == Design::Stacking::Variant::SameZTheta ? 1 : 2);
    stacking_ref_->setCurrentIndex(static_cast<int>(st.reference));
    stacking_z_->setValue(st.z_ref);
    stacking_theta_->setValue(st.theta_ref_deg);
    if (auto* b = layers_->button(s_.active_layer)) b->setChecked(true);
    const auto& l = layer();
    method_->setCurrentIndex(static_cast<int>(l.method));
    stagger_->setValue(l.stagger_deg);
    theta_le_->setValue(l.theta_le_deg);
    solidity_->setValue(l.solidity);
    chord_scale_->setValue(l.chord_scale);
    cl0_->setValue(l.naca.cl0);
    tc_->setValue(l.naca.t_max_over_c);
    variant_->setCurrentText(QString::fromStdString(Naca::to_string(l.naca.variant)));
    for (auto& [name, k] : custom_) k->setValue(Profile::get_param(l.custom, name));
    for (std::size_t i = 0; i < ps_ctrl_.size(); ++i) {
        ps_ctrl_[i]->setEnabled(i < l.custom.ps.mprime_ctrl.size());
        if (i < l.custom.ps.mprime_ctrl.size()) ps_ctrl_[i]->setValue(l.custom.ps.mprime_ctrl[i]);
    }
    connector_->setCurrentText(QString::fromStdString(Profile::to_string(l.custom.connector)));
    enforce_c1_->setChecked(l.custom.thickness.enforce_c1);
    n_blades_->setValue(s_.design.n_blades);
    mirror_->setChecked(s_.design.mirror);
    rows_->setValue(s_.design.mesh.n_rows);
    n_body_->setValue(s_.design.mesh.n_body);
    n_edge_->setValue(s_.design.mesh.n_edge);
    clustering_->setValue(s_.design.mesh.row_clustering);
    updateEnabled();
    syncing_ = false;
}

void ParameterPanel::takeStackingPoint()
{
    const auto R = s_.result;
    if (!R || R->layers.empty() || !R->layers.front().error.empty()) return;
    const auto& L = R->layers.front();
    auto& st = s_.design.stacking;
    switch (st.reference) {
    case Design::Stacking::Reference::CG:
        st.z_ref = L.cg_xr.x();
        st.theta_ref_deg = Math::rad2deg(L.cg_theta);
        break;
    case Design::Stacking::Reference::LE:
        st.z_ref = L.camber_xr.front().x();
        st.theta_ref_deg = Math::rad2deg(L.theta_le);
        break;
    case Design::Stacking::Reference::TE:
        st.z_ref = L.camber_xr.back().x();
        st.theta_ref_deg = L.theta_te_deg;
        break;
    }
    const bool was = syncing_;
    syncing_ = true;
    stacking_z_->setValue(st.z_ref);
    stacking_theta_->setValue(st.theta_ref_deg);
    syncing_ = was;
}

void ParameterPanel::setFitStatus(const QString& t, bool running)
{
    fit_status_->setText(t);
    fit_->setEnabled(!running && layer().method == Design::Method::Naca65Custom);
}

}  // namespace Designer
