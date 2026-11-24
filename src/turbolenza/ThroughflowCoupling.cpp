/// @file turbolenza/ThroughflowCoupling.cpp
/// @brief Chapter 6 and Chapter 8 as the plug-ins of Chapter 7's through-flow (see the header).
#include "turbolenza/ThroughflowCoupling.hpp"

#include <cmath>
#include <format>
#include <numbers>
#include <string>
#include <utility>

namespace BladeGen::Turbolenza::Throughflow {

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kDeg = 180.0 / kPi;

[[noreturn]] void fail(const std::string& what)
{
    throw ThroughflowError(what);
}

bool positive(double x) noexcept { return std::isfinite(x) && x > 0.0; }

/// The profile family and camberline named by BladeRowSpec::profile_family (see the header).
bool family_of(const std::string& name, ProfileFamily& pf, Camberline& cl)
{
    if (name == "NACA65") {
        pf = ProfileFamily::Naca;
        cl = Camberline::Naca65;
    } else if (name == "C4") {
        pf = ProfileFamily::C4;
        cl = Camberline::CircularArc;
    } else if (name == "DCA") {
        pf = ProfileFamily::DoubleCircularArc;
        cl = Camberline::CircularArc;
    } else if (name == "A4K6") {
        pf = ProfileFamily::Naca;
        cl = Camberline::NacaA4K6;
    } else {
        return false;
    }
    return true;
}

/// R and gamma of a perfect-gas state: P = rho R T, a^2 = gamma R T (ours: one gas for both modules).
Gas gas_from_state(double P, double rho, double T, double a, const char* where)
{
    if (!(positive(P) && positive(rho) && positive(T) && positive(a)))
        fail(std::format("{}: the state (P, rho, T, a) = ({}, {}, {}, {}) is not physical", where, P, rho, T, a));
    Gas g;
    g.gas_constant = P / (rho * T);
    g.gamma = a * a * rho / P;
    return g;
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// CorrelationBladeRow
// ─────────────────────────────────────────────────────────────────────────────

CorrelationBladeRow::CorrelationBladeRow(const BladeRowSpec& row, CorrelationRowOptions options) : opt_(std::move(options))
{
    if (row.n_blades <= 0) fail(row.name + ": CorrelationBladeRow needs the blade count Z > 0");
    for (double r : {row.le_hub_r, row.le_shroud_r, row.te_hub_r, row.te_shroud_r})
        if (!positive(r)) fail(row.name + ": CorrelationBladeRow needs the LE and TE end radii (> 0) for s_bar of Eq. 6-99");
    ProfileFamily pf = ProfileFamily::Naca;
    Camberline cl = Camberline::Naca65;
    const bool known = family_of(row.profile_family, pf, cl);
    if (!known && !(opt_.profile && opt_.camberline))
        fail(std::format("{}: the profile family '{}' has no Ch. 6 model (NACA65, C4, DCA, A4K6); give "
                         "CorrelationRowOptions::profile and camberline",
                         row.name, row.profile_family));
    profile_ = opt_.profile ? *opt_.profile : pf;
    camberline_ = opt_.camberline ? *opt_.camberline : cl;
    // s_bar = 2 pi r_bar / Z, r_bar the mean of the LE and TE end radii (ours; Eqs. 6-99, 6-100)
    const double rbar = 0.25 * (row.le_hub_r + row.le_shroud_r + row.te_hub_r + row.te_shroud_r);
    mean_pitch_ = 2.0 * kPi * rbar / static_cast<double>(row.n_blades);
}

Gas CorrelationBladeRow::gas_of(const BladeRowInlet& in)
{
    return gas_from_state(in.P1, in.rho1, in.T1, in.a1, "CorrelationBladeRow (inlet)");
}

InletState CorrelationBladeRow::inlet_state(const BladeRowInlet& in)
{
    InletState s;
    s.radius = in.r1;
    s.wm = in.Wm1;
    s.wtheta = in.Wtheta1;  // C_theta1 - omega r1 (Eq. 7-44), the row's frame
    s.omega = in.omega;
    s.mach = in.mach1_rel;  // M'_1 = W1 / a1
    s.density = in.rho1;
    s.viscosity = in.mu1;  // Sutherland at the static T1 (the through-flow's gas)
    return s;
}

DischargeEstimate CorrelationBladeRow::discharge_estimate(const BladeRowInlet& in)
{
    return {in.r2, in.Wm2};
}

CascadeGeometry CorrelationBladeRow::cascade_geometry(const BladeRowInlet& in) const
{
    const BladeSectionGeometry& s = in.geometry;
    if (!(std::isfinite(s.kappa1) && std::isfinite(s.kappa2) && std::isfinite(s.camber) && std::isfinite(s.stagger)))
        fail("CorrelationBladeRow: a non-finite section angle");
    if (!(std::abs(s.camber - (s.kappa1 - s.kappa2)) <= 1e-9))
        fail(std::format("CorrelationBladeRow: the section's camber {} rad is not kappa1 - kappa2 = {} rad (Eq. 6-6)",
                         s.camber, s.kappa1 - s.kappa2));
    if (!positive(in.pitch)) fail("CorrelationBladeRow: the pitch must be > 0");
    if (!positive(in.blade_height)) fail("CorrelationBladeRow: the blade height must be > 0 (Eq. 6-99)");
    CascadeGeometry g;
    g.profile = profile_;
    g.camberline = camberline_;
    g.kind = in.is_igv ? BladeRowKind::InletGuideVane : (in.omega != 0.0 ? BladeRowKind::Rotor : BladeRowKind::Stator);
    g.row_sign = in.s_row;
    g.chord = s.chord;
    g.solidity = in.solidity;
    g.camber = s.kappa1 - s.kappa2;
    g.thickness_ratio = s.thickness_ratio;
    // a / c: 0.5 by definition for the circular arc and the NACA 65 equivalent arc (Eq. 4-35, pp. 120, 121), so
    // the section's measured value (e.g. BladeGen's largest camber ordinate) is not used there (ours)
    const bool arc = camberline_ == Camberline::Naca65 || camberline_ == Camberline::CircularArc;
    g.max_camber_position = arc ? 0.5 : s.max_camber_position;
    g.throat_over_pitch = s.throat / in.pitch;  // o at the section, s at r_mean (ours)
    g.pitch_over_height = mean_pitch_ / in.blade_height;
    g.stagger = s.stagger;
    if (opt_.stagger == StaggerSource::MetalAngles) {
        // the cascade frame of Correlations.hpp: s_row, and p. 128's flip of a guide vane of negative camber
        int f = in.s_row;
        if (static_cast<double>(f) * g.camber < 0.0 && in.is_igv) f = -f;
        const double fd = static_cast<double>(f);
        const double theta_c = fd * g.camber * kDeg;
        // a non-IGV row of negative cascade camber is left to evaluate(), which refuses it (p. 128)
        if (theta_c >= 0.0 && theta_c < 90.0) {
            const BladeAngles ba = blade_angles(camberline_, 0.0, theta_c, g.max_camber_position);
            g.stagger = (fd * s.kappa1 * kDeg - ba.chi1_deg) / kDeg / fd;  // gamma_c = kappa1_c - chi1, back to the book
        }
    }
    return g;
}

CascadePerformance CorrelationBladeRow::performance(const BladeRowInlet& in) const
{
    const CascadeGeometry g = cascade_geometry(in);
    CascadePerformance p = Turbolenza::evaluate(g, inlet_state(in), discharge_estimate(in), gas_of(in), opt_.correlations);
    if (opt_.stagger == StaggerSource::MetalAngles && g.stagger != in.geometry.stagger)
        p.notes.push_back(std::format("stagger {} deg from the metal angles (ours); the section's geometric stagger is {} deg",
                                      g.stagger * kDeg, in.geometry.stagger * kDeg));
    return p;
}

BladeRowDischarge CorrelationBladeRow::evaluate(const BladeRowInlet& in) const
{
    const CascadePerformance p = performance(in);
    BladeRowDischarge d;
    d.beta2 = p.beta2;              // book frame (Correlations.hpp)
    d.loss_coefficient = p.loss;    // Eq. 7-49 form
    d.wake_blockage = p.k_w;        // Eq. 7-6
    return d;
}

// ─────────────────────────────────────────────────────────────────────────────
// CorrelationRowLoss
// ─────────────────────────────────────────────────────────────────────────────

Turbolenza::WallState correlation_wall_state(const WallState& w)
{
    Turbolenza::WallState s;
    s.radius = w.r;
    s.density = w.rho;
    s.cm = w.Cm;
    s.ctheta = w.Ctheta;  // absolute, book sign
    s.pressure = w.P;
    s.temperature = w.T;
    return s;
}

CorrelationRowLoss::CorrelationRowLoss(CorrelationRowLossOptions options) : opt_(std::move(options)) {}

namespace {

/// The row's gap (see the header): none, a tip gap at the wall with a clearance, or the seal of a shrouded row.
struct Gap
{
    bool any = false;
    EndWall side = EndWall::Shroud;
    double clearance = 0.0;
};

Gap gap_of(const RowLossInput& in, EndWall seal_wall)
{
    if (!(std::isfinite(in.hub_clearance) && std::isfinite(in.tip_clearance) && in.hub_clearance >= 0.0 &&
          in.tip_clearance >= 0.0))
        fail("CorrelationRowLoss: the clearances must be finite and >= 0");
    Gap g;
    if (in.shrouded) {
        g.any = true;
        g.side = seal_wall;
        g.clearance = seal_wall == EndWall::Hub ? in.hub_clearance : in.tip_clearance;
        const double other = seal_wall == EndWall::Hub ? in.tip_clearance : in.hub_clearance;
        if (other != 0.0)
            fail(std::format("CorrelationRowLoss: a shrouded row has its seal at the {} and no gap at the other wall (got "
                             "{} m)",
                             seal_wall == EndWall::Hub ? "hub" : "shroud", other));
        return g;
    }
    if (in.hub_clearance > 0.0 && in.tip_clearance > 0.0)
        fail("CorrelationRowLoss: clearances at both walls (RowLossResult carries one clearance loss)");
    if (in.tip_clearance > 0.0) g = {true, EndWall::Shroud, in.tip_clearance};
    else if (in.hub_clearance > 0.0) g = {true, EndWall::Hub, in.hub_clearance};
    return g;
}

}  // namespace

RowLeakage CorrelationRowLoss::leakage(const RowLossInput& in) const
{
    const Gap gap = gap_of(in, opt_.seal_wall);
    if (!gap.any) return {};  // no gap, no row-level loss
    const bool hub = gap.side == EndWall::Hub;
    const WallState& w1 = hub ? in.hub1 : in.shroud1;
    const WallState& w2 = hub ? in.hub2 : in.shroud2;
    const BladeSectionGeometry& sec = hub ? in.hub_section : in.shroud_section;
    RowWallData d;
    d.shrouded = in.shrouded;
    d.gap_side = gap.side;
    d.row_sign = in.s_row;
    d.clearance = gap.clearance;
    d.gap_radius = (in.shrouded && opt_.seal.radius) ? *opt_.seal.radius : 0.5 * (w1.r + w2.r);  // ours
    d.blade_count = in.n_blades;
    d.chord = sec.chord;
    d.stagger = sec.stagger;
    d.row_number = in.row_number;
    d.fins = opt_.seal.fins;
    d.fin_thickness = opt_.seal.fin_thickness;
    d.fin_pitch = opt_.seal.fin_pitch;
    const Gas gas = gas_from_state(w1.P, w1.rho, w1.T, w1.a, "CorrelationRowLoss (inlet wall)");
    return endwall_leakage(d, correlation_wall_state(w1), correlation_wall_state(w2), in.mass_flow, gas);
}

RowLossResult CorrelationRowLoss::evaluate(const RowLossInput& in) const
{
    const Gap gap = gap_of(in, opt_.seal_wall);
    RowLossResult r;
    if (!gap.any) return r;
    const RowLeakage L = leakage(in);
    const RowLossResult::Side side = gap.side == EndWall::Hub ? RowLossResult::Side::Hub : RowLossResult::Side::Shroud;
    if (in.shrouded) {
        r.seal_dPt = L.total_pressure_loss;
        r.seal_side = side;
        r.leakage_mass_flow = static_cast<double>(L.seal_flow_direction) * L.leakage_mass_flow;  // + discharge -> inlet
    } else {
        r.clearance_dPt = L.total_pressure_loss;
        r.clearance_side = side;
    }
    return r;
}

BladeRowSpec with_correlations(BladeRowSpec row, const CorrelationRowOptions& row_options,
                               const CorrelationRowLossOptions& loss_options)
{
    row.model = std::make_shared<CorrelationBladeRow>(row, row_options);
    row.row_loss = std::make_shared<CorrelationRowLoss>(loss_options);
    return row;
}

// ─────────────────────────────────────────────────────────────────────────────
// ChapterEightEndwall
// ─────────────────────────────────────────────────────────────────────────────

ChapterEightEndwall::ChapterEightEndwall(EndwallCouplingOptions options) : opt_(std::move(options)) {}

void ChapterEightEndwall::begin_pass(int pass, int n_qn)
{
    if (n_qn < 1) fail("ChapterEightEndwall: a pass needs quasi-normals");
    pass_ = pass;
    rec_.assign(static_cast<std::size_t>(n_qn), Record{});
}

EndwallWallStation ChapterEightEndwall::wall_station(const EndwallEdge& e)
{
    EndwallWallStation w;
    w.z = e.z;
    w.r = e.r;
    w.phi = e.phi;
    w.m = e.m;
    w.wall_omega = e.wall_omega;
    w.Vm = e.Cm;          // C_m = W_m = V_m (Eq. 8-3)
    w.Ctheta = e.Ctheta;  // absolute; Ch. 8 forms V_theta = C_theta - r omega_w itself (Eq. 8-4)
    w.rho = e.rho;
    w.P = e.P;
    w.mu = e.mu;
    return w;
}

Turbolenza::EndwallRowWall ChapterEightEndwall::row_wall(const EndwallRow& row, EndWall wall) const
{
    const EndwallRowWall& w = wall == EndWall::Hub ? row.hub : row.shroud;
    Turbolenza::EndwallRowWall o;
    o.omega = row.omega;
    o.blade_count = row.n_blades;
    o.pitch = w.pitch_mean;      // s_bar (Eq. 8-2)
    o.stagger = w.stagger_mean;  // gamma_bar
    if (row.shrouded && wall == opt_.seal_wall) {
        o.shrouded = true;
        o.clearance = 0.0;  // the band's gap is the seal's (Sec. 8.6; ours, the q3d contract Sec. 2.4)
        o.seal_leakage = row.leakage_mass_flow;
    } else {
        o.clearance = w.clearance;
    }
    return o;
}

EndwallResult ChapterEightEndwall::evaluate(const EndwallStation& st)
{
    if (st.qn < 0 || static_cast<std::size_t>(st.qn) >= rec_.size() || st.pass != pass_)
        fail(std::format("ChapterEightEndwall: quasi-normal {} of pass {} outside the pass begun ({} quasi-normals, pass {})",
                         st.qn, st.pass, rec_.size(), pass_));
    const std::size_t i = static_cast<std::size_t>(st.qn);
    Record r;
    r.hub_wall = wall_station(st.hub);
    r.shroud_wall = wall_station(st.shroud);
    if (i == 0) {  // Sec. 8.8 step 1 at the first station (Eq. 8-72)
        r.hub = endwall_inlet(r.hub_wall, opt_.inlet, st.mass_flow, opt_.settings);
        r.shroud = endwall_inlet(r.shroud_wall, opt_.inlet, st.mass_flow, opt_.settings);
    } else {
        const Record& p = rec_[i - 1];
        if (!p.valid)
            fail(std::format("ChapterEightEndwall: quasi-normal {} before quasi-normal {} in pass {} (Fig. 7-3 marches in "
                             "order)",
                             st.qn, st.qn - 1, st.pass));
        const bool row = st.row_upstream >= 0;
        Turbolenza::EndwallRowWall hw, sw;
        if (row) {
            hw = row_wall(st.upstream, EndWall::Hub);
            sw = row_wall(st.upstream, EndWall::Shroud);
        }
        const EndwallStepResult h =
            endwall_step(p.hub, p.hub_wall, r.hub_wall, row ? &hw : nullptr, st.compressor_omega, opt_.settings);
        const EndwallStepResult s =
            endwall_step(p.shroud, p.shroud_wall, r.shroud_wall, row ? &sw : nullptr, st.compressor_omega, opt_.settings);
        r.hub = h.layer;
        r.shroud = s.layer;
        r.row = h.row;
        r.hub_row = h.row_result;
        r.shroud_row = s.row_result;
    }
    r.blockage = endwall_blockage(r.hub, r.hub_wall, r.shroud, r.shroud_wall, st.mass_flow);  // Eqs. 8-76, 8-77
    if (!(r.blockage.K_B > 0.0 && r.blockage.K_B <= 1.0))
        fail(std::format("ChapterEightEndwall: quasi-normal {}: the end-wall blockage B = {} leaves no flow area (Eq. 8-77)",
                         st.qn, r.blockage.B));
    r.valid = true;
    rec_[i] = r;
    for (std::size_t k = i + 1; k < rec_.size(); ++k) rec_[k].valid = false;  // what followed was marched from the old layer
    EndwallResult out;
    out.blockage_factor = r.blockage.K_B;
    return out;
}

}  // namespace BladeGen::Turbolenza::Throughflow
