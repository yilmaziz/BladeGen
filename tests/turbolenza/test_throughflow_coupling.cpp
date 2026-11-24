// ThroughflowCoupling.hpp: Chapter 6 (Correlations.hpp) as Chapter 7's BladeRowModel and RowLossModel, Chapter 8
// (Endwall.hpp) as its EndwallModel.  The mapping field by field (frames, s_row, the gas, Re_c, o / s, s / h), the
// row-level losses against Eqs. 6-85 .. 6-89 written out here, the end-wall march against Ch. 8's own march and
// Eqs. 8-72, 8-76, 8-77, and the coupled analysis of a synthetic stage (IGV, rotor, shrouded stator) and of
// BladeGen's default rotor: convergence, mass, Euler work = cp dT0 per stream surface, the Ch. 6 values of one
// surface rebuilt by hand from the converged nodes, Eq. 7-49, K_B's effect, the mirrored stage.  Every measured
// number is printed ("MEASURE ..." or "[coupling] ...").
#include "CouplingFixture.hpp"

#include "turbolenza/ThroughflowBladeGen.hpp"
#include "turbolenza/ThroughflowCoupling.hpp"

#include "blade/FlowPath.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>

namespace TF = BladeGen::Turbolenza::Throughflow;
namespace T = BladeGen::Turbolenza;
namespace D = BladeGen::Design;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kDeg = kPi / 180.0;
constexpr double kRh = 0.25, kRs = 0.35;

double measure(const char* name, double v)
{
    std::printf("MEASURE %-60s %.3e\n", name, v);
    return v;
}

/// A section of the plain table: metal angles in degrees (book frame), stagger (kappa1 + kappa2) / 2.
TF::BladeSectionGeometry section(double span, double k1_deg, double k2_deg, double axial_chord, double tb, int Z)
{
    TF::BladeSectionGeometry s;
    s.span = span;
    s.kappa1 = k1_deg * kDeg;
    s.kappa2 = k2_deg * kDeg;
    s.camber = s.kappa1 - s.kappa2;
    s.stagger = 0.5 * (s.kappa1 + s.kappa2);
    s.chord = axial_chord / std::cos(s.stagger);
    s.thickness_ratio = tb;
    s.max_camber_position = 0.5;
    s.max_thickness_position = 0.4;
    const double r = kRh + span * (kRs - kRh);
    s.throat = 0.9 * std::cos(s.stagger) * 2.0 * kPi * r / Z;  // o / s = 0.9 cos(gamma) (ours, a plausible throat)
    return s;
}

TF::BladeRowSpec stage_row(const char* name, double z1, double z2, int Z, int s_row, double omega)
{
    TF::BladeRowSpec r;
    r.name = name;
    r.n_blades = Z;
    r.s_row = s_row;
    r.omega = omega;
    r.profile_family = "NACA65";
    const double gap = 0.008;  // the LE and TE lines 8 mm inside the quasi-normals
    r.le_hub_z = r.le_shroud_z = z1 + gap;
    r.te_hub_z = r.te_shroud_z = z2 - gap;
    r.le_hub_r = r.te_hub_r = kRh;
    r.le_shroud_r = r.te_shroud_r = kRs;
    return r;
}

/// The QN1 state of a mid-span surface of a rotor, by hand (perfect gas gamma, R).
TF::BladeRowInlet rotor_inlet(double gamma = 1.4, double R = 287.05)
{
    TF::BladeRowInlet in;
    in.row = 1;
    in.row_number = 2;
    in.stream_surface = 5;
    in.n_stream_surfaces = 11;
    in.mass_fraction = 0.5;
    in.span = 0.5;
    in.omega = -523.6;
    in.n_blades = 32;
    in.s_row = 1;
    in.profile_family = "NACA65";
    in.blade_height = 0.1;
    in.mass_flow = 22.0;
    in.z1 = 0.12;
    in.r1 = 0.30;
    in.z2 = 0.18;
    in.r2 = 0.302;
    in.pitch = 2.0 * kPi * 0.5 * (in.r1 + in.r2) / in.n_blades;
    in.geometry = section(0.5, 54.0, 34.0, 0.044, 0.07, 32);
    in.solidity = in.geometry.chord / in.pitch;
    in.Ctheta1 = -17.8;
    in.Wm1 = 101.0;
    in.Wtheta1 = in.Ctheta1 - in.r1 * in.omega;
    in.W1 = std::hypot(in.Wm1, in.Wtheta1);
    in.beta1 = std::atan2(in.Wtheta1, in.Wm1);
    in.T1 = 280.0;
    in.P1 = 92000.0;
    in.rho1 = in.P1 / (R * in.T1);
    in.a1 = std::sqrt(gamma * R * in.T1);
    in.mach1_rel = in.W1 / in.a1;
    in.mu1 = 1.716e-5 * std::pow(in.T1 / 273.15, 1.5) * (273.15 + 110.4) / (in.T1 + 110.4);  // Sutherland at T1
    in.reynolds = in.rho1 * in.W1 * in.geometry.chord / in.mu1;
    const double cp = gamma * R / (gamma - 1.0);
    in.h1 = cp * in.T1;
    in.Hrel1 = in.h1 + 0.5 * in.W1 * in.W1;
    in.Tt1_rel = in.Hrel1 / cp;
    in.Pt1_rel = in.P1 * std::pow(in.Tt1_rel / in.T1, gamma / (gamma - 1.0));
    in.Wm2 = 97.0;
    in.rho2 = in.rho1;
    return in;
}

TF::BladeRowSpec rotor_spec()
{
    return stage_row("rotor", 0.12, 0.18, 32, 1, -523.6);
}

/// in mirrored: theta -> -theta (every signed angle, swirl, omega, s_row).
TF::BladeRowInlet mirrored(TF::BladeRowInlet in)
{
    in.omega = -in.omega;
    in.s_row = -in.s_row;
    in.Ctheta1 = -in.Ctheta1;
    in.Wtheta1 = -in.Wtheta1;
    in.beta1 = -in.beta1;
    in.alpha1 = -in.alpha1;
    in.geometry.kappa1 = -in.geometry.kappa1;
    in.geometry.kappa2 = -in.geometry.kappa2;
    in.geometry.camber = -in.geometry.camber;
    in.geometry.stagger = -in.geometry.stagger;
    return in;
}

double max_mass_error(const TF::ThroughflowResult& r)
{
    double e = 0.0;
    for (const auto& s : r.stations) e = std::max(e, std::abs(s.mass_error));
    return e;
}

/// int_0^ys f dy by the trapezoid rule over the nodes.
template <class F>
double trapezoid(const TF::ThroughflowStation& st, F f)
{
    double s = 0.0;
    for (std::size_t j = 0; j + 1 < st.nodes.size(); ++j)
        s += 0.5 * (f(st.nodes[j]) + f(st.nodes[j + 1])) * (st.nodes[j + 1].y - st.nodes[j].y);
    return s;
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// the fixtures of CouplingFixture.hpp
// ─────────────────────────────────────────────────────────────────────────────

namespace coupling {

TF::ThroughflowCase synthetic_stage(const StageOptions& o)
{
    const double f = o.mirror ? -1.0 : 1.0;
    const double omega = f * (-2.0 * kPi * o.rpm / 60.0);  // book sign: rpm > 0 turns toward -theta
    TF::ThroughflowCase c;
    c.model = o.model;
    c.mass_flow = 22.0;
    c.compressor_omega = omega;
    c.options.n_stream_surfaces = o.n_stream_surfaces;
    c.options.max_outer_iterations = 1000;
    const double z[] = {0.0, 0.06, 0.12, 0.18, 0.24, 0.30, 0.38};
    for (int k = 0; k < 7; ++k) {
        TF::QuasiNormal q;
        q.z_hub = q.z_shroud = z[k];
        q.r_hub = kRh;
        q.r_shroud = kRs;
        if (k == 3) q.hub_wall_omega = omega;  // the drum under the rotor (QN2 -> QN3); the stator band is at rest
        c.quasi_normals.push_back(q);
    }
    c.inlet = TF::InletProfile::uniform(101325.0, 288.15);
    // the guide vane: axial in, about -10 deg out (toward the rotation, a free vortex); s_row = +1 (kappa1 > kappa2)
    TF::BladeRowSpec igv = stage_row("igv", 0.06, 0.12, 30, 1, 0.0);
    igv.is_igv = true;
    igv.sections = {section(0.0, 0.0, -15.0, 0.044, 0.08, 30), section(0.5, 0.0, -13.0, 0.044, 0.08, 30),
                    section(1.0, 0.0, -11.5, 0.044, 0.08, 30)};
    // the rotor (a free vortex of about 5.8 kJ/kg by design, ours): beta'_1 about 47 .. 59 deg, beta'_2 about
    // 33 .. 54 deg, kappa_2 6 deg below it; a 0.5 mm tip gap
    TF::BladeRowSpec rot = stage_row("rotor", 0.12, 0.18, 32, 1, omega);
    rot.tip_clearance = 5e-4;
    rot.sections = {section(0.0, 47.3, 27.0, 0.044, 0.07, 32), section(0.5, 54.0, 39.5, 0.044, 0.07, 32),
                    section(1.0, 59.0, 47.5, 0.044, 0.07, 32)};
    // the stator: alpha about -33 .. -25 deg in, axial out: kappa1 < kappa2, so s_row = -1; a hub shroud band, seal.
    // Its loading stays below Lieblein's D_eq = 2: with D_eq near 2 at the hub (44 vanes, alpha_1 -37 deg) the
    // plain substitution of the discharge quasi-normal does not settle (the K_W of Eq. 6-84 against W_m2; measured)
    TF::BladeRowSpec sta = stage_row("stator", 0.18, 0.24, 60, -1, 0.0);
    sta.shrouded = true;
    sta.hub_clearance = 3e-4;
    sta.sections = {section(0.0, -33.0, 5.0, 0.044, 0.08, 60), section(0.5, -28.4, 5.0, 0.044, 0.08, 60),
                    section(1.0, -24.8, 5.0, 0.044, 0.08, 60)};
    igv.qn_inlet = 1;
    rot.qn_inlet = 2;
    sta.qn_inlet = 3;
    TF::CorrelationRowLossOptions lo;
    lo.seal.fins = 3;
    lo.seal.fin_thickness = 5e-4;
    lo.seal.fin_pitch = 4e-3;
    for (TF::BladeRowSpec* r : {&igv, &rot, &sta}) {
        if (o.mirror) {
            r->s_row = -r->s_row;
            for (auto& s : r->sections) {
                s.kappa1 = -s.kappa1;
                s.kappa2 = -s.kappa2;
                s.camber = -s.camber;
                s.stagger = -s.stagger;
            }
        }
        TF::CorrelationRowOptions ro;
        ro.correlations.reynolds_correction = o.reynolds_correction;
        *r = TF::with_correlations(*r, ro, lo);
        if (!o.row_loss) r->row_loss.reset();
    }
    c.rows = {igv, rot, sta};
    if (o.endwall) c.endwall = std::make_shared<TF::ChapterEightEndwall>();
    return c;
}

TF::ThroughflowCase bladegen_rotor_case(TF::ThroughflowModel model, int* qn_before_le)
{
    D::Evaluator ev;
    TF::BladeGenThroughflowOptions opt;
    opt.rpm = 15000.0;
    const TF::BladeGenThroughflowGeometry g = TF::throughflow_geometry(D::BladeDesign::defaults(), ev, opt);
    TF::ThroughflowCase c;
    c.model = model;
    c.quasi_normals = g.quasi_normals;
    c.inlet = TF::InletProfile::uniform(101325.0, 288.15);
    c.mass_flow = 3.0;
    c.compressor_omega = g.compressor_omega;
    c.options.max_outer_iterations = 3000;
    TF::BladeRowSpec row = g.row;
    row.tip_clearance = 3e-4;
    TF::CorrelationRowOptions ro;  // BladeGen's default profile is "custom": a NACA 65-(12)10 fit (ours)
    ro.profile = T::ProfileFamily::Naca;
    ro.camberline = T::Camberline::Naca65;
    c.rows = {TF::with_correlations(row, ro)};
    c.endwall = std::make_shared<TF::ChapterEightEndwall>();
    if (qn_before_le) *qn_before_le = g.qn_before_le;
    return c;
}

TF::ThroughflowCase cone_rotor_case()
{
    const double omega = -2.0 * kPi * 5000.0 / 60.0, t = std::tan(15.0 * kDeg);
    TF::ThroughflowCase c;
    c.model = TF::ThroughflowModel::Approximate;
    c.mass_flow = 22.0;
    c.compressor_omega = omega;
    const double z[] = {0.0, 0.05, 0.10, 0.16, 0.22};
    for (int k = 0; k < 5; ++k) {
        TF::QuasiNormal q;
        q.z_hub = q.z_shroud = z[k];
        q.r_hub = 0.22 + z[k] * t;
        q.r_shroud = 0.32 + z[k] * t;
        if (k == 3) q.hub_wall_omega = omega;  // the drum under the rotor
        c.quasi_normals.push_back(q);
    }
    c.inlet = TF::InletProfile::uniform(101325.0, 288.15, TF::InletProfile::Swirl::Ctheta, 20.0);
    TF::BladeRowSpec rot;
    rot.name = "cone rotor";
    rot.qn_inlet = 2;
    rot.omega = omega;
    rot.n_blades = 32;
    rot.le_hub_z = rot.le_shroud_z = 0.108;
    rot.te_hub_z = rot.te_shroud_z = 0.152;
    rot.le_hub_r = 0.22 + 0.108 * t;
    rot.le_shroud_r = 0.32 + 0.108 * t;
    rot.te_hub_r = 0.22 + 0.152 * t;
    rot.te_shroud_r = 0.32 + 0.152 * t;
    rot.tip_clearance = 5e-4;
    rot.sections = {section(0.0, 50.0, 30.0, 0.044, 0.07, 32), section(1.0, 60.0, 45.0, 0.044, 0.07, 32)};
    rot.model = std::make_shared<TF::SpecifiedBladeRow>(std::vector<double>{0.0, 1.0}, std::vector<double>{0.6, 0.6},
                                                        std::vector<double>{0.03, 0.03});
    c.rows = {rot};
    return c;
}

}  // namespace coupling

// ─────────────────────────────────────────────────────────────────────────────
// Chapter 6 -> Chapter 7, per stream surface
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("coupling: the Ch. 6 row adapter maps the Ch. 7 state field by field", "[turbolenza][coupling]")
{
    const TF::BladeRowSpec spec = rotor_spec();
    const TF::CorrelationBladeRow row(spec);
    const TF::BladeRowInlet in = rotor_inlet();
    // the inlet of Correlations.hpp: Eq. 7-44's W_theta1 in the row frame, M'_1, Sutherland's mu at the static T1
    const T::CascadeInletState is = TF::CorrelationBladeRow::inlet_state(in);
    CHECK(is.radius == 0.30);
    CHECK(is.wm == 101.0);
    CHECK(is.wtheta == -17.8 + 0.30 * 523.6);
    CHECK(is.omega == -523.6);
    CHECK(is.mach == in.W1 / in.a1);
    CHECK(is.density == in.rho1);
    CHECK(is.viscosity == in.mu1);
    const T::DischargeEstimate de = TF::CorrelationBladeRow::discharge_estimate(in);
    CHECK(de.radius == 0.302);
    CHECK(de.wm == 97.0);
    // the gas from the state: R = P / (rho T), gamma = a^2 rho / P, exactly, also for gamma != 1.4
    const T::Gas g14 = TF::CorrelationBladeRow::gas_of(in);
    CHECK(std::abs(g14.gas_constant / 287.05 - 1.0) <= 1e-14);
    CHECK(std::abs(g14.gamma / 1.4 - 1.0) <= 1e-14);
    const T::Gas g13 = TF::CorrelationBladeRow::gas_of(rotor_inlet(1.3, 300.0));
    CHECK(std::abs(g13.gas_constant / 300.0 - 1.0) <= 1e-14);
    CHECK(std::abs(g13.gamma / 1.3 - 1.0) <= 1e-14);
    // the geometry: o / s at the surface's pitch, s_bar / h_bar of the row (s_bar at the mean of the LE / TE radii)
    const T::CascadeGeometry cg = row.cascade_geometry(in);
    CHECK(cg.row_sign == 1);
    CHECK(cg.kind == T::BladeRowKind::Rotor);
    CHECK(cg.profile == T::ProfileFamily::Naca);
    CHECK(cg.camberline == T::Camberline::Naca65);
    CHECK(cg.chord == in.geometry.chord);
    CHECK(cg.solidity == in.solidity);
    CHECK(std::abs(cg.camber - 20.0 * kDeg) <= 1e-15);
    CHECK(cg.thickness_ratio == 0.07);
    CHECK(std::abs(cg.throat_over_pitch - in.geometry.throat / in.pitch) <= 1e-15);
    CHECK(std::abs(row.mean_pitch() - 2.0 * kPi * 0.30 / 32.0) <= 1e-15);
    CHECK(std::abs(cg.pitch_over_height - 2.0 * kPi * 0.30 / 32.0 / 0.1) <= 1e-14);
    // MetalAngles: the NACA 65 equivalent arc (chi1 = chi2 = theta / 2) through the section's kappa: gamma = 44 deg
    CHECK(std::abs(cg.stagger - 44.0 * kDeg) <= 1e-14);
    // the result: Re_c = rho1 W1 c / mu(T1) as the through-flow's, beta'_2 = kappa_2 + delta (6-3), i = beta1 - kappa1
    const T::CascadePerformance p = row.performance(in);
    CHECK(std::abs(p.reynolds / in.reynolds - 1.0) <= 1e-14);
    CHECK(std::abs(p.kappa1_deg - 54.0) <= 1e-12);
    CHECK(std::abs(p.kappa2_deg - 34.0) <= 1e-12);
    CHECK(std::abs(p.incidence_deg - (in.beta1 / kDeg - 54.0)) <= 1e-12);
    CHECK(std::abs(p.beta2 - (in.geometry.kappa2 + p.deviation_deg * kDeg)) <= 1e-14);
    CHECK(p.deviation_deg > 0.0);  // a compressor cascade: the flow leaves short of the metal (beta2 > kappa2)
    // the same as Correlations' evaluate with the inputs built here
    T::CascadeGeometry hand;
    hand.profile = T::ProfileFamily::Naca;
    hand.camberline = T::Camberline::Naca65;
    hand.kind = T::BladeRowKind::Rotor;
    hand.row_sign = 1;
    hand.chord = in.geometry.chord;
    hand.solidity = in.geometry.chord / (2.0 * kPi * 0.301 / 32.0);
    hand.stagger = 44.0 * kDeg;
    hand.camber = 20.0 * kDeg;
    hand.thickness_ratio = 0.07;
    hand.max_camber_position = 0.5;
    hand.throat_over_pitch = 0.9 * std::cos(44.0 * kDeg) * 0.30 / 0.301;
    hand.pitch_over_height = 2.0 * kPi * 0.30 / 32.0 / 0.1;
    T::CascadeInletState hi;
    hi.radius = 0.30;
    hi.wm = 101.0;
    hi.wtheta = -17.8 + 0.30 * 523.6;
    hi.omega = -523.6;
    hi.mach = std::hypot(101.0, hi.wtheta) / std::sqrt(1.4 * 287.05 * 280.0);
    hi.density = 92000.0 / (287.05 * 280.0);
    hi.viscosity = in.mu1;
    const T::CascadePerformance q = T::evaluate(hand, hi, {0.302, 97.0}, {1.4, 287.05});
    CHECK(std::abs(p.beta2 - q.beta2) <= 1e-13);
    CHECK(std::abs(p.loss - q.loss) <= 1e-13);
    CHECK(std::abs(p.k_w - q.k_w) <= 1e-13);
    const TF::BladeRowDischarge d = row.evaluate(in);
    CHECK(d.beta2 == p.beta2);
    CHECK(d.loss_coefficient == p.losses.total);  // the whole of Eq. 6-80 .. 6-82 incl. the shock part
    CHECK(d.wake_blockage == p.k_w);
    std::printf("[coupling] rotor mid surface by hand: i %.3f deg, i* %.3f, delta %.3f, delta* %.3f, beta'2 %.3f deg, "
                "omega %.5f (design %.5f), K_W %.4f, M'1 %.3f, Re_c %.3e\n",
                p.incidence_deg, p.design_incidence_deg, p.deviation_deg, p.design_deviation_deg, p.beta2 / kDeg, p.loss,
                p.losses.design, p.k_w, in.mach1_rel, p.reynolds);
    // the Reynolds correction sees the Sutherland Re_c: K_Re of Eq. 6-101 / 6-102 at in.reynolds
    {
        TF::CorrelationRowOptions ro;
        ro.correlations.reynolds_correction = true;
        const TF::CorrelationBladeRow rr(spec, ro);
        const T::CascadePerformance pr = rr.performance(in);
        const double k2 = 1.0 + hand.pitch_over_height * std::cos(p.design_beta2_deg * kDeg) +
                          0.004 * T::reynolds_factor(in.reynolds) / 0.0073;  // Eqs. 6-99, 6-100
        CHECK(std::abs(pr.losses.k2 - k2) <= 1e-13);
    }
    // the mirrored surface (theta -> -theta, s_row -1): the mirrored discharge angle, the same loss and K_W
    const TF::BladeRowDischarge dm = row.evaluate(mirrored(in));
    CHECK(std::abs(dm.beta2 + d.beta2) <= 1e-15);
    CHECK(std::abs(dm.loss_coefficient - d.loss_coefficient) <= 1e-15);
    CHECK(dm.wake_blockage == d.wake_blockage);
    // Geometric stagger: the section's own chord line; Ch. 4's kappa then differ from the section's
    {
        TF::BladeRowInlet ig = in;
        ig.geometry.stagger = 46.0 * kDeg;  // 2 deg off the metal angles' mean
        TF::CorrelationRowOptions go;
        go.stagger = TF::StaggerSource::Geometric;
        const T::CascadePerformance pg = TF::CorrelationBladeRow(spec, go).performance(ig);
        CHECK(std::abs(pg.kappa2_deg - 36.0) <= 1e-12);
        const T::CascadePerformance pm = row.performance(ig);  // MetalAngles ignores it: the section's kappa
        CHECK(std::abs(pm.kappa2_deg - 34.0) <= 1e-12);
        CHECK(std::abs(pm.beta2 - p.beta2) <= 1e-15);
        bool noted = false;
        for (const auto& n : pm.notes) noted = noted || n.find("from the metal angles") != std::string::npos;
        CHECK(noted);
    }
    // a parabolic-arc camberline (a / c = 0.4): gamma = kappa1 - chi1 keeps both metal angles
    {
        TF::CorrelationRowOptions po;
        po.camberline = T::Camberline::ParabolicArc;
        TF::BladeRowInlet ip = in;
        ip.geometry.max_camber_position = 0.4;
        const T::CascadePerformance pp = TF::CorrelationBladeRow(spec, po).performance(ip);
        CHECK(std::abs(pp.kappa1_deg - 54.0) <= 1e-12);
        CHECK(std::abs(pp.kappa2_deg - 34.0) <= 1e-12);
        CHECK(pp.chi1_deg > pp.chi2_deg);  // front-loaded
    }
}

TEST_CASE("coupling: the Ch. 6 row adapter's families, guide vanes and refusals", "[turbolenza][coupling]")
{
    TF::BladeRowSpec spec = rotor_spec();
    const TF::BladeRowInlet in = rotor_inlet();
    // the families by name
    for (const auto& [name, pf, cl] : {std::tuple{"NACA65", T::ProfileFamily::Naca, T::Camberline::Naca65},
                                       std::tuple{"C4", T::ProfileFamily::C4, T::Camberline::CircularArc},
                                       std::tuple{"DCA", T::ProfileFamily::DoubleCircularArc, T::Camberline::CircularArc},
                                       std::tuple{"A4K6", T::ProfileFamily::Naca, T::Camberline::NacaA4K6}}) {
        spec.profile_family = name;
        const TF::CorrelationBladeRow r(spec);
        CHECK(r.profile() == pf);
        CHECK(r.camberline() == cl);
    }
    spec.profile_family = "custom";  // BladeGen's custom profile: refused unless given
    CHECK_THROWS_AS(TF::CorrelationBladeRow(spec), TF::ThroughflowError);
    TF::CorrelationRowOptions given;
    given.profile = T::ProfileFamily::C4;
    CHECK_THROWS_AS(TF::CorrelationBladeRow(spec, given), TF::ThroughflowError);  // both are needed
    given.camberline = T::Camberline::CircularArc;
    CHECK(TF::CorrelationBladeRow(spec, given).profile() == T::ProfileFamily::C4);
    spec = rotor_spec();
    spec.n_blades = 0;
    CHECK_THROWS_AS(TF::CorrelationBladeRow(spec), TF::ThroughflowError);
    spec = rotor_spec();
    spec.te_shroud_r = 0.0;
    CHECK_THROWS_AS(TF::CorrelationBladeRow(spec), TF::ThroughflowError);
    const TF::CorrelationBladeRow row(rotor_spec());
    // the state and geometry the adapter refuses
    TF::BladeRowInlet bad = in;
    bad.P1 = 0.0;
    CHECK_THROWS_AS(row.evaluate(bad), TF::ThroughflowError);
    bad = in;
    bad.blade_height = 0.0;
    CHECK_THROWS_AS(row.evaluate(bad), TF::ThroughflowError);
    bad = in;
    bad.geometry.camber += 0.01;  // not kappa1 - kappa2
    CHECK_THROWS_AS(row.evaluate(bad), TF::ThroughflowError);
    bad = in;
    bad.geometry.kappa1 = std::nan("");
    CHECK_THROWS_AS(row.evaluate(bad), TF::ThroughflowError);
    // what Ch. 6 refuses: the wrong row sign (negative cascade camber), no throat
    bad = in;
    bad.s_row = -1;
    CHECK_THROWS_AS(row.evaluate(bad), T::CorrelationError);
    bad = in;
    bad.geometry.throat = 0.0;
    CHECK_THROWS_AS(row.evaluate(bad), T::CorrelationError);
    // a guide vane: Howell's deviation, no incidence term (p. 144): beta'_2 does not move with beta'_1
    {
        TF::BladeRowSpec gs = stage_row("igv", 0.06, 0.12, 30, 1, 0.0);
        gs.is_igv = true;
        const TF::CorrelationBladeRow igv(gs);
        TF::BladeRowInlet gi = in;
        gi.omega = 0.0;
        gi.is_igv = true;
        gi.geometry = section(0.5, 0.0, -13.0, 0.044, 0.08, 30);
        gi.pitch = 2.0 * kPi * 0.3 / 30;
        gi.solidity = gi.geometry.chord / gi.pitch;
        gi.Ctheta1 = gi.Wtheta1 = 0.0;
        gi.beta1 = 0.0;
        gi.W1 = gi.Wm1;
        gi.mach1_rel = gi.W1 / gi.a1;
        gi.Wm2 = gi.Wm1;
        const T::CascadeGeometry gg = igv.cascade_geometry(gi);
        CHECK(gg.kind == T::BladeRowKind::InletGuideVane);
        CHECK(std::abs(gg.stagger + 6.5 * kDeg) <= 1e-14);
        const TF::BladeRowDischarge a = igv.evaluate(gi);
        gi.Wtheta1 = gi.Ctheta1 = gi.Wm1 * std::tan(3.0 * kDeg);
        gi.beta1 = 3.0 * kDeg;
        gi.W1 = std::hypot(gi.Wm1, gi.Wtheta1);
        gi.mach1_rel = gi.W1 / gi.a1;
        const TF::BladeRowDischarge b = igv.evaluate(gi);
        CHECK(std::abs(a.beta2 - b.beta2) <= 1e-15);
        CHECK(a.beta2 < 0.0);  // toward the rotation (-theta)
        CHECK(a.beta2 > gi.geometry.kappa2);  // with Howell's deviation, short of the metal
        // the same vane given with s_row = -1 and its camber mirrored (p. 128's flip on top of s_row)
        TF::BladeRowSpec gm = gs;
        gm.s_row = -1;
        TF::BladeRowInlet mi = mirrored(gi);
        const TF::BladeRowDischarge c = TF::CorrelationBladeRow(gm).evaluate(mi);
        CHECK(std::abs(c.beta2 + b.beta2) <= 1e-15);
        // and with s_row = +1 kept: a guide vane of negative camber, p. 128's flip gives the same mirrored answer
        mi.s_row = 1;
        const TF::CorrelationBladeRow igv_plus(gs);
        CHECK(igv_plus.cascade_geometry(mi).stagger == -igv.cascade_geometry(gi).stagger);
        CHECK(std::abs(igv_plus.evaluate(mi).beta2 + b.beta2) <= 1e-15);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Chapter 6 -> Chapter 7, the row-level losses
// ─────────────────────────────────────────────────────────────────────────────

namespace {

TF::WallState wall(double r, double Cm, double Ct, double P, double T)
{
    TF::WallState w;
    w.r = r;
    w.Cm = Cm;
    w.Ctheta = Ct;
    w.P = P;
    w.T = T;
    w.rho = P / (287.05 * T);
    w.a = std::sqrt(1.4 * 287.05 * T);
    w.mu = 1.8e-5;
    return w;
}

TF::RowLossInput rotor_loss_input()
{
    TF::RowLossInput in;
    in.row = 1;
    in.row_number = 2;
    in.omega = -523.6;
    in.n_blades = 32;
    in.s_row = 1;
    in.tip_clearance = 5e-4;
    in.mass_flow = 22.0;
    in.blade_height = 0.1;
    in.hub_section = section(0.0, 47.3, 27.0, 0.044, 0.07, 32);
    in.shroud_section = section(1.0, 59.0, 39.0, 0.050, 0.07, 32);  // an axial chord of its own: c cos(gamma) enters 6-86, 6-88
    in.hub1 = wall(0.25, 100.0, -21.4, 92000.0, 280.0);
    in.shroud1 = wall(0.35, 102.0, -15.3, 93500.0, 281.0);
    in.hub2 = wall(0.25, 96.0, -67.8, 99000.0, 287.0);
    in.shroud2 = wall(0.354, 98.0, -85.8, 101000.0, 292.0);  // a rising casing: the gap radius is the mean (ours)
    return in;
}

}  // namespace

TEST_CASE("coupling: the Ch. 6 row-level loss adapter (Eqs. 6-85 .. 6-91) and its gaps", "[turbolenza][coupling]")
{
    const TF::CorrelationRowLoss model;
    TF::RowLossInput in = rotor_loss_input();
    // the rotor's tip gap at the shroud, Eqs. 6-85 .. 6-89 written out (cascade signs: s_row = +1 here; |tau|)
    {
        const TF::WallState &a = in.shroud1, &b = in.shroud2;
        const double dc = 5e-4, c = in.shroud_section.chord, cg = std::cos(in.shroud_section.stagger), Z = 32.0;
        const double tau = kPi * dc * (a.r * a.rho * a.Cm + b.r * b.rho * b.Cm) * std::abs(b.r * b.Ctheta - a.r * a.Ctheta);
        const double rtip = 0.5 * (a.r + b.r);  // the mean of the wall radii (ours)
        const double dP = tau / (Z * rtip * dc * c * cg);
        const double rho = 0.5 * (a.rho + b.rho);
        const double Uc = 0.816 * std::sqrt(2.0 * dP / rho) / std::pow(2.0, 0.2);  // N_row = 2
        const double mc = rho * Uc * Z * dc * c * cg;
        const double dPt = dP * mc / 22.0;
        const TF::RowLossResult r = model.evaluate(in);
        std::printf("[coupling] rotor tip gap by hand: tau %.4f N m, dP %.1f Pa, U_c %.2f m/s, mdot_c %.4f kg/s, dPt %.2f Pa\n",
                    tau, dP, Uc, mc, dPt);
        CHECK(std::abs(r.clearance_dPt / dPt - 1.0) <= 1e-13);
        CHECK(r.clearance_side == TF::RowLossResult::Side::Shroud);
        CHECK(r.seal_dPt == 0.0);
        CHECK(r.leakage_mass_flow == 0.0);  // a tip gap sends nothing to Ch. 8's seal model
        CHECK(std::abs(model.leakage(in).leakage_mass_flow / mc - 1.0) <= 1e-13);
        // the mirrored row (s_row = -1, every swirl sign-changed) loses the same
        TF::RowLossInput mi = in;
        mi.s_row = -1;
        mi.omega = -mi.omega;
        for (TF::WallState* w : {&mi.hub1, &mi.hub2, &mi.shroud1, &mi.shroud2}) w->Ctheta = -w->Ctheta;
        for (auto* s : {&mi.hub_section, &mi.shroud_section}) s->stagger = -s->stagger;
        CHECK(std::abs(model.evaluate(mi).clearance_dPt / dPt - 1.0) <= 1e-13);
    }
    // a cantilevered stator's hub gap: the hub wall and section
    {
        TF::RowLossInput h = in;
        h.tip_clearance = 0.0;
        h.hub_clearance = 4e-4;
        const TF::RowLossResult r = model.evaluate(h);
        CHECK(r.clearance_side == TF::RowLossResult::Side::Hub);
        CHECK(r.clearance_dPt > 0.0);
        T::RowWallData d;
        d.gap_side = T::EndWall::Hub;
        d.clearance = 4e-4;
        d.gap_radius = 0.25;
        d.blade_count = 32;
        d.chord = h.hub_section.chord;
        d.stagger = h.hub_section.stagger;
        d.row_number = 2;
        const T::RowLeakage L = T::endwall_leakage(d, TF::correlation_wall_state(h.hub1), TF::correlation_wall_state(h.hub2),
                                                   22.0, {1.4, 287.05});
        CHECK(std::abs(r.clearance_dPt / L.total_pressure_loss - 1.0) <= 1e-13);
    }
    // a shrouded stator, the seal at the hub: Egli's leakage with its sign (+ discharge -> inlet when P2 > P1)
    {
        TF::CorrelationRowLossOptions so;
        so.seal.fins = 3;
        so.seal.fin_thickness = 5e-4;
        so.seal.fin_pitch = 4e-3;
        const TF::CorrelationRowLoss seal(so);
        TF::RowLossInput s = in;
        s.omega = 0.0;
        s.s_row = -1;
        s.shrouded = true;
        s.tip_clearance = 0.0;
        s.hub_clearance = 3e-4;
        const TF::RowLossResult r = seal.evaluate(s);
        CHECK(r.seal_side == TF::RowLossResult::Side::Hub);
        CHECK(r.clearance_dPt == 0.0);
        CHECK(r.seal_dPt > 0.0);
        CHECK(r.leakage_mass_flow > 0.0);  // P2 = 99000 > P1 = 92000 at the hub
        T::RowWallData d;
        d.shrouded = true;
        d.gap_side = T::EndWall::Hub;
        d.row_sign = -1;
        d.clearance = 3e-4;
        d.gap_radius = 0.25;
        d.blade_count = 32;
        d.chord = s.hub_section.chord;
        d.stagger = s.hub_section.stagger;
        d.row_number = 2;
        d.fins = 3;
        d.fin_thickness = 5e-4;
        d.fin_pitch = 4e-3;
        const T::RowLeakage L = T::endwall_leakage(d, TF::correlation_wall_state(s.hub1), TF::correlation_wall_state(s.hub2),
                                                   22.0, {1.4, 287.05});
        CHECK(std::abs(r.seal_dPt / L.total_pressure_loss - 1.0) <= 1e-13);
        CHECK(std::abs(r.leakage_mass_flow / L.leakage_mass_flow - 1.0) <= 1e-13);
        std::swap(s.hub1, s.hub2);  // P2 < P1: the seal flows the other way
        CHECK(seal.evaluate(s).leakage_mass_flow < 0.0);
        so.seal.radius = 0.24;  // a given seal radius scales Eq. 6-91
        std::swap(s.hub1, s.hub2);
        CHECK(std::abs(TF::CorrelationRowLoss(so).evaluate(s).leakage_mass_flow / r.leakage_mass_flow - 0.24 / 0.25) <= 1e-13);
        // refusals: a shrouded row with a gap at its other wall; seal data missing (Ch. 6 refuses)
        TF::RowLossInput t = s;
        t.tip_clearance = 1e-4;
        CHECK_THROWS_AS(seal.evaluate(t), TF::ThroughflowError);
        CHECK_THROWS_AS(model.evaluate(s), T::CorrelationError);  // no fins
    }
    // no gap: no loss; gaps at both walls and negative gaps are refused
    in.tip_clearance = 0.0;
    const TF::RowLossResult z = model.evaluate(in);
    CHECK(z.clearance_dPt == 0.0);
    CHECK(z.seal_dPt == 0.0);
    in.tip_clearance = 5e-4;
    in.hub_clearance = 5e-4;
    CHECK_THROWS_AS(model.evaluate(in), TF::ThroughflowError);
    in.hub_clearance = -1e-4;
    CHECK_THROWS_AS(model.evaluate(in), TF::ThroughflowError);
}

// ─────────────────────────────────────────────────────────────────────────────
// Chapter 8 -> Chapter 7
// ─────────────────────────────────────────────────────────────────────────────

namespace {

/// The wall geometry of a march input, per station: arc length along each wall and the contour angle.
struct WallGeometry
{
    std::vector<double> m_hub, m_shroud, phi_hub, phi_shroud;
};

/// Straight walls between the stations (a cylinder, a cone): m the sum of the segment lengths, phi = atan2(dr, dz)
/// of the segment (the same for every segment of a cone; the first segment's at station 0).
WallGeometry straight_walls(const TF::ThroughflowResult& r)
{
    WallGeometry g;
    for (int side = 0; side < 2; ++side) {
        std::vector<double>& m = side == 0 ? g.m_hub : g.m_shroud;
        std::vector<double>& phi = side == 0 ? g.phi_hub : g.phi_shroud;
        for (std::size_t k = 0; k < r.stations.size(); ++k) {
            const auto& a = side == 0 ? r.stations[k == 0 ? 0 : k - 1].nodes.front() : r.stations[k == 0 ? 0 : k - 1].nodes.back();
            const auto& b = side == 0 ? r.stations[k == 0 ? 1 : k].nodes.front() : r.stations[k == 0 ? 1 : k].nodes.back();
            m.push_back(k == 0 ? 0.0 : m.back() + std::hypot(b.z - a.z, b.r - a.r));
            phi.push_back(std::atan2(b.r - a.r, b.z - a.z));
        }
    }
    return g;
}

/// The Ch. 8 march input of a converged through-flow: its wall edges, the wall geometry given, its rows (from the
/// case, s_bar at the wall's radii on the two quasi-normals), its mdot; a shrouded row's seal leakage recovered from
/// the adapter's record (Eq. 8-59: mdot_leak = Delta+ 2 pi r rho_e V_me).
T::EndwallInput march_input(const TF::ThroughflowCase& c, const TF::ThroughflowResult& r, const TF::ChapterEightEndwall& ew,
                            const WallGeometry& wg)
{
    T::EndwallInput in;
    in.mass_flow = c.mass_flow;
    in.compressor_omega = c.compressor_omega;
    for (std::size_t k = 0; k < r.stations.size(); ++k) {
        const auto& st = r.stations[k];
        T::EndwallStation s;
        for (int side = 0; side < 2; ++side) {
            const TF::ThroughflowNode& nd = side == 0 ? st.nodes.front() : st.nodes.back();
            T::EndwallWallStation& w = side == 0 ? s.hub : s.shroud;
            w.z = nd.z;
            w.r = nd.r;
            w.phi = side == 0 ? wg.phi_hub[k] : wg.phi_shroud[k];
            w.m = side == 0 ? wg.m_hub[k] : wg.m_shroud[k];
            w.Vm = nd.Wm;
            w.Ctheta = nd.Ctheta;
            w.rho = nd.rho;
            w.P = nd.P;
            w.mu = nd.mu;
            w.wall_omega = side == 0 ? c.quasi_normals[k].hub_wall_omega : c.quasi_normals[k].shroud_wall_omega;
        }
        s.quasi_normal_length = st.y_s;
        in.stations.push_back(s);
    }
    for (const auto& row : c.rows) {
        T::EndwallBladeRow b;
        b.station_in = row.qn_inlet;
        b.station_out = row.qn_inlet + 1;
        for (int side = 0; side < 2; ++side) {
            T::EndwallRowWall& w = side == 0 ? b.hub : b.shroud;
            const auto& sec = side == 0 ? row.sections.front() : row.sections.back();
            const auto& n1 = r.stations[static_cast<std::size_t>(b.station_in)].nodes;
            const auto& n2 = r.stations[static_cast<std::size_t>(b.station_out)].nodes;
            const double r1 = side == 0 ? n1.front().r : n1.back().r, r2 = side == 0 ? n2.front().r : n2.back().r;
            w.omega = row.omega;
            w.blade_count = row.n_blades;
            w.pitch = 2.0 * kPi * 0.5 * (r1 + r2) / row.n_blades;
            w.stagger = sec.stagger;
            w.clearance = side == 0 ? row.hub_clearance : row.tip_clearance;
            if (row.shrouded && side == 0) {
                w.shrouded = true;
                w.clearance = 0.0;
                const auto& rec = ew.records()[static_cast<std::size_t>(b.station_out)];
                const auto& wi = ew.records()[static_cast<std::size_t>(b.station_in)].hub_wall;
                w.seal_leakage = rec.hub_row.leakage.delta_plus * 2.0 * kPi * wi.r * wi.rho * wi.Vm;
            }
        }
        in.rows.push_back(b);
    }
    return in;
}

}  // namespace

TEST_CASE("coupling: the Ch. 8 adapter in an annulus without rows (Eqs. 8-72, 8-76, 8-77; the march, purity, order)",
          "[turbolenza][coupling]")
{
    coupling::StageOptions so;
    TF::ThroughflowCase c = coupling::synthetic_stage(so);
    c.rows.clear();
    for (auto& q : c.quasi_normals) q.hub_wall_omega = 0.0;
    auto ew = std::make_shared<TF::ChapterEightEndwall>();
    c.endwall = ew;
    const TF::ThroughflowResult r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    CHECK(max_mass_error(r) <= 1e-11);
    REQUIRE(ew->records().size() == r.stations.size());
    // QN 0: B = 0.02 shared by the walls (Eq. 8-72), K_B = 0.98 exactly (Eq. 8-77)
    CHECK(std::abs(r.stations[0].K_B_endwall - 0.98) <= 1e-15);
    double kb_min = 1.0;
    for (std::size_t k = 0; k < r.stations.size(); ++k) {
        const auto& st = r.stations[k];
        const auto& rec = ew->records()[k];
        REQUIRE(rec.valid);
        // the station's K_B is Eq. 8-76 / 8-77 of ITS OWN walls: B = sum 2 pi r rho_e V_me delta1* / mdot
        const auto &h = st.nodes.front(), &s = st.nodes.back();
        const double B = 2.0 * kPi * (h.r * h.rho * h.Wm * rec.hub.delta1 + s.r * s.rho * s.Wm * rec.shroud.delta1) / c.mass_flow;
        CHECK(std::abs(st.K_B_endwall - (1.0 - B)) <= 1e-13);
        // the edge handed over is the end-wall stream surface's
        CHECK(rec.hub_wall.Vm == h.Wm);
        CHECK(rec.shroud_wall.Ctheta == s.Ctheta);
        CHECK(rec.hub_wall.mu == h.mu);
        // Eq. 7-5 with K_B: the node data (trapezoid) carry the case's mass flow
        const double mdot = 2.0 * kPi * st.K_B_endwall * trapezoid(st, [](const TF::ThroughflowNode& n) {
                                return n.K_W * n.r * n.rho * n.Wm * std::cos(n.epsilon);
                            });
        CHECK(std::abs(mdot / c.mass_flow - 1.0) <= 1e-4);  // the trapezoid against the solver's three-point rule
        kb_min = std::min(kb_min, st.K_B_endwall);
    }
    CHECK(kb_min < 0.98);  // the layers grow along the walls
    // the adapter's march is Ch. 8's own march on the converged edges, to round-off
    const T::EndwallResult m = T::analyze_endwalls(march_input(c, r, *ew, straight_walls(r)));
    double dk = 0.0;
    for (std::size_t k = 0; k < r.stations.size(); ++k) dk = std::max(dk, std::abs(m.stations[k].blockage.K_B - r.stations[k].K_B_endwall));
    CHECK(measure("Ch. 8 adapter vs analyze_endwalls, annulus: max |dK_B|", dk) <= 1e-12);
    std::printf("[coupling] annulus K_B:");
    for (const auto& st : r.stations) std::printf(" %.5f", st.K_B_endwall);
    std::printf("\n");
    // purity and order: the same station twice gives the same K_B; a station out of order is refused
    TF::ChapterEightEndwall e2;
    TF::EndwallStation s0;
    s0.qn = 0;
    s0.n_qn = 3;
    s0.mass_flow = 22.0;
    for (TF::EndwallEdge* w : {&s0.hub, &s0.shroud}) {
        w->r = w == &s0.hub ? kRh : kRs;
        w->Cm = 100.0;
        w->rho = 1.1;
        w->P = 95000.0;
        w->mu = 1.8e-5;
        w->m = 0.0;
    }
    CHECK_THROWS_AS(e2.evaluate(s0), TF::ThroughflowError);  // no pass begun
    e2.begin_pass(0, 3);
    TF::EndwallStation s2 = s0;
    s2.qn = 2;
    s2.hub.m = s2.shroud.m = 0.1;
    CHECK_THROWS_AS(e2.evaluate(s2), TF::ThroughflowError);  // QN 2 before QN 1
    const double k0 = e2.evaluate(s0).blockage_factor;
    CHECK(e2.evaluate(s0).blockage_factor == k0);
    TF::EndwallStation s1 = s0;
    s1.qn = 1;
    s1.hub.m = s1.shroud.m = 0.05;
    const double k1 = e2.evaluate(s1).blockage_factor;
    CHECK(e2.evaluate(s1).blockage_factor == k1);
    CHECK(k1 < k0);
    s1.pass = 1;
    CHECK_THROWS_AS(e2.evaluate(s1), TF::ThroughflowError);  // not the pass begun
}

TEST_CASE("coupling: the Ch. 8 adapter's row walls (s_bar, gamma_bar, clearance, a shrouded band)", "[turbolenza][coupling]")
{
    TF::EndwallRow row;
    row.row = 2;
    row.omega = 0.0;
    row.n_blades = 44;
    row.shrouded = true;
    row.leakage_mass_flow = -0.031;  // P2 < P1: the seal flows from the inlet to the discharge (Ch. 8 takes the sign)
    row.hub.pitch_mean = 0.0357;
    row.hub.stagger_mean = -0.25;
    row.hub.clearance = 3e-4;
    row.shroud.pitch_mean = 0.05;
    row.shroud.stagger_mean = -0.3;
    row.shroud.clearance = 0.0;
    const TF::ChapterEightEndwall hub_band;
    const T::EndwallRowWall h = hub_band.row_wall(row, T::EndWall::Hub);
    CHECK(h.shrouded);
    CHECK(h.clearance == 0.0);  // the band's gap is the seal's (Sec. 8.6)
    CHECK(h.seal_leakage == -0.031);
    CHECK(h.pitch == 0.0357);
    CHECK(h.stagger == -0.25);
    CHECK(h.blade_count == 44);
    const T::EndwallRowWall s = hub_band.row_wall(row, T::EndWall::Shroud);
    CHECK_FALSE(s.shrouded);
    CHECK(s.seal_leakage == 0.0);
    CHECK(s.pitch == 0.05);
    TF::EndwallCouplingOptions so;
    so.seal_wall = T::EndWall::Shroud;
    const T::EndwallRowWall s2 = TF::ChapterEightEndwall(so).row_wall(row, T::EndWall::Shroud);
    CHECK(s2.shrouded);
    CHECK(s2.seal_leakage == -0.031);
    const T::EndwallRowWall h2 = TF::ChapterEightEndwall(so).row_wall(row, T::EndWall::Hub);
    CHECK_FALSE(h2.shrouded);
    CHECK(h2.clearance == 3e-4);
    row.shrouded = false;  // an unshrouded rotor: its tip gap
    row.omega = -523.6;
    row.shroud.clearance = 5e-4;
    const T::EndwallRowWall t = hub_band.row_wall(row, T::EndWall::Shroud);
    CHECK(t.clearance == 5e-4);
    CHECK(t.omega == -523.6);
    CHECK_FALSE(hub_band.row_wall(row, T::EndWall::Hub).shrouded);
}

// ─────────────────────────────────────────────────────────────────────────────
// the coupled analysis: a synthetic stage
// ─────────────────────────────────────────────────────────────────────────────

namespace {

/// The Ch. 6 inputs of surface j of row k rebuilt by hand from the converged nodes and evaluated with Correlations.
T::CascadePerformance hand_ch6(const TF::ThroughflowCase& c, const TF::ThroughflowResult& r, int k, std::size_t j,
                            const T::CorrelationOptions& opt = {})
{
    const TF::BladeRowSpec& row = c.rows[static_cast<std::size_t>(k)];
    const auto& n1 = r.stations[static_cast<std::size_t>(row.qn_inlet)].nodes[j];
    const auto& n2 = r.stations[static_cast<std::size_t>(row.qn_inlet + 1)].nodes[j];
    // the span: the mean of the LE and TE crossings of the segment n1 -> n2 (radial edge lines: the radius there)
    const double tl = (row.le_hub_z - n1.z) / (n2.z - n1.z), tt = (row.te_hub_z - n1.z) / (n2.z - n1.z);
    const double rl = n1.r + tl * (n2.r - n1.r), rt = n1.r + tt * (n2.r - n1.r);
    const double span = 0.5 * ((rl - kRh) + (rt - kRh)) / (kRs - kRh);
    // the section by linear interpolation in span (three sections at 0, 0.5, 1)
    const auto& S = row.sections;
    const std::size_t a = span < 0.5 ? 0 : 1;
    const double t = (span - S[a].span) / (S[a + 1].span - S[a].span);
    auto lin = [&](double TF::BladeSectionGeometry::*f) { return S[a].*f + t * (S[a + 1].*f - S[a].*f); };
    const double k1 = lin(&TF::BladeSectionGeometry::kappa1), k2 = lin(&TF::BladeSectionGeometry::kappa2);
    const double chord = lin(&TF::BladeSectionGeometry::chord), o = lin(&TF::BladeSectionGeometry::throat);
    const double pitch = 2.0 * kPi * 0.5 * (n1.r + n2.r) / row.n_blades;
    T::CascadeGeometry g;
    g.profile = T::ProfileFamily::Naca;
    g.camberline = T::Camberline::Naca65;
    g.kind = row.is_igv ? T::BladeRowKind::InletGuideVane : (row.omega != 0.0 ? T::BladeRowKind::Rotor : T::BladeRowKind::Stator);
    g.row_sign = row.s_row;
    g.chord = chord;
    g.solidity = chord / pitch;
    g.stagger = 0.5 * (k1 + k2);  // the equivalent circular arc through the metal angles
    g.camber = k1 - k2;
    g.thickness_ratio = lin(&TF::BladeSectionGeometry::thickness_ratio);
    g.max_camber_position = 0.5;
    g.throat_over_pitch = o / pitch;
    g.pitch_over_height = 2.0 * kPi * 0.5 * (kRh + kRs) / row.n_blades / (kRs - kRh);
    T::CascadeInletState in;
    in.radius = n1.r;
    in.wm = n1.Wm;
    in.wtheta = n1.Ctheta - row.omega * n1.r;  // Eq. 7-44 (QN1 is solved in the frame of the row before)
    in.omega = row.omega;
    in.mach = std::hypot(in.wm, in.wtheta) / n1.a;
    in.density = n1.rho;
    in.viscosity = 1.716e-5 * std::pow(n1.T / 273.15, 1.5) * (273.15 + 110.4) / (n1.T + 110.4);
    return T::evaluate(g, in, {n2.r, n2.Wm}, {1.4, 287.05}, opt);
}

}  // namespace

TEST_CASE("coupling: a synthetic stage (IGV, rotor, shrouded stator) with Ch. 6 and Ch. 8", "[turbolenza][coupling]")
{
    // approximate, full, and approximate with Ch. 6's Reynolds correction (K_Re of the node's Re_c = rho1 W1 c / mu(T1))
    for (int v = 0; v < 3; ++v) {
        const auto model = v == 1 ? TF::ThroughflowModel::FullNormalEquilibrium : TF::ThroughflowModel::Approximate;
        coupling::StageOptions so;
        so.model = model;
        so.reynolds_correction = v == 2;
        T::CorrelationOptions copt;
        copt.reynolds_correction = so.reynolds_correction;
        const TF::ThroughflowCase c = coupling::synthetic_stage(so);
        const TF::ThroughflowResult r = TF::solve_throughflow(c);
        INFO(r.message);
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        const char* mn = v == 0 ? "approximate" : v == 1 ? "full" : "approximate, K_Re";
        std::printf("[coupling] stage, %s: %zu passes, max mass error %.2e\n", mn, r.history.size(), max_mass_error(r));
        CHECK(max_mass_error(r) <= 1e-10);  // measured 1.1e-11 (approximate and full)
        const double cp = c.gas.cp();
        const auto& q1 = r.stations[2];  // rotor inlet
        const auto& q2 = r.stations[3];  // rotor discharge
        const auto& q3 = r.stations[4];  // stator discharge
        double euler = 0.0, rothalpy = 0.0;
        for (std::size_t j = 0; j < q1.nodes.size(); ++j) {
            const auto &a = q1.nodes[j], &b = q2.nodes[j], &s = q3.nodes[j];
            // Euler work = cp dT0 on every stream surface (Eqs. 7-47, 7-51): omega (r2 C_theta2 - r1 C_theta1)
            const double work = c.rows[1].omega * (b.rCtheta - a.rCtheta);
            euler = std::max(euler, std::abs(cp * (b.Tt - a.Tt) / work - 1.0));
            rothalpy = std::max(rothalpy, std::abs((b.H - c.rows[1].omega * b.rCtheta) - (a.H - c.rows[1].omega * a.rCtheta)));
            CHECK(work > 0.0);                      // a compressor rotor (omega < 0, r C_theta falls)
            CHECK(b.s > a.s);                       // Ch. 6's loss
            CHECK(std::abs(s.H - b.H) <= 1e-9 * b.H);  // the stator does no work
            CHECK(std::abs(s.alpha) < 5.0 * kDeg);  // the stator turns the flow back to about axial
            CHECK(b.Pt > a.Pt);
        }
        CHECK(measure("stage: max |cp dT0 / (omega d(r C_theta)) - 1| over the surfaces", euler) <= 1e-12);
        CHECK(measure("stage: max |I2 - I1| [J/kg] (Eq. 7-47)", rothalpy) <= 1e-8);
        // the Ch. 6 values of the mid surface rebuilt by hand from the converged nodes, every row
        const std::size_t jm = q1.nodes.size() / 2;
        for (int k = 0; k < 3; ++k) {
            const auto& nd2 = r.stations[static_cast<std::size_t>(c.rows[static_cast<std::size_t>(k)].qn_inlet + 1)].nodes[jm];
            const T::CascadePerformance p = hand_ch6(c, r, k, jm, copt);
            std::printf("[coupling] stage %s, %s mid surface: beta'2 %.4f deg (hand %.4f), omega %.5f (hand %.5f), "
                        "row loss %.5f, K_W %.4f (hand %.4f), i %.2f, delta %.2f deg, D_eq %.3f, M'1 %.3f\n",
                        mn, c.rows[static_cast<std::size_t>(k)].name.c_str(), nd2.beta / kDeg, p.beta2 / kDeg,
                        nd2.loss_coefficient - nd2.row_loss_coefficient, p.loss, nd2.row_loss_coefficient, nd2.K_W, p.k_w,
                        p.incidence_deg, p.deviation_deg, p.d_eq, std::hypot(nd2.Wm, nd2.Wtheta) / nd2.a);
            CHECK(std::abs(nd2.beta - p.beta2) <= 1e-7);  // W_m2 converged to the flow tolerance
            CHECK(std::abs((nd2.loss_coefficient - nd2.row_loss_coefficient) - p.loss) <= 1e-8);
            CHECK(std::abs(nd2.K_W - p.k_w) <= 1e-8);
        }
        // Eq. 7-49 on every surface of the rotor: P't2 = P't2,id - omegabar (P't1 - P1)
        double e749 = 0.0;
        for (std::size_t j = 0; j < q1.nodes.size(); ++j) {
            const auto &a = q1.nodes[j], &b = q2.nodes[j];
            const double w = c.rows[1].omega;
            const double Wt1 = a.Ctheta - w * a.r, Hrel1 = cp * a.T + 0.5 * (a.Wm * a.Wm + Wt1 * Wt1);
            const double Pt1 = a.P * std::pow(Hrel1 / (cp * a.T), 3.5);
            const double Hrel2 = (a.H - w * a.rCtheta) + 0.5 * (w * b.r) * (w * b.r);  // Eq. 7-48
            const double Pt2id = Pt1 * std::pow(Hrel2 / Hrel1, 3.5);                    // at s1
            const double Pt2 = b.P * std::pow(Hrel2 / (cp * b.T), 3.5);
            e749 = std::max(e749, std::abs((Pt2id - Pt2) / (Pt1 - a.P) - b.loss_coefficient));
        }
        CHECK(measure("stage: Eq. 7-49 from the nodes, max |omegabar - loss|", e749) <= 1e-9);
        // the row-level losses: the rotor's tip gap is zero at the hub (p. 147), the stator's hub seal zero at the shroud
        CHECK(q2.nodes.front().row_loss_coefficient == 0.0);
        CHECK(q2.nodes.back().row_loss_coefficient > 0.0);
        CHECK(q3.nodes.back().row_loss_coefficient == 0.0);
        CHECK(q3.nodes.front().row_loss_coefficient > 0.0);
        // the end walls: K_B of every station from its own walls, and Ch. 8's march on the converged edges
        const auto ew = std::dynamic_pointer_cast<TF::ChapterEightEndwall>(c.endwall);
        REQUIRE(ew);
        const T::EndwallResult m = T::analyze_endwalls(march_input(c, r, *ew, straight_walls(r)));
        double dk = 0.0;
        std::printf("[coupling] stage %s K_B:", mn);
        for (std::size_t k = 0; k < r.stations.size(); ++k) {
            std::printf(" %.5f", r.stations[k].K_B_endwall);
            dk = std::max(dk, std::abs(m.stations[k].blockage.K_B - r.stations[k].K_B_endwall));
        }
        std::printf("\n");
        CHECK(measure("stage: Ch. 8 adapter vs analyze_endwalls on the converged edges, max |dK_B|", dk) <= 1e-9);
        CHECK(r.stations[0].K_B_endwall == 1.0 - 0.02);
        CHECK(r.stations.back().K_B_endwall < r.stations[0].K_B_endwall);
    }
}

TEST_CASE("coupling: K_B's effect, and the mirrored stage (theta -> -theta, s_row -> -s_row)", "[turbolenza][coupling]")
{
    coupling::StageOptions so;
    const TF::ThroughflowCase c = coupling::synthetic_stage(so);
    const TF::ThroughflowResult r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    so.endwall = false;
    const TF::ThroughflowCase c0 = coupling::synthetic_stage(so);
    const TF::ThroughflowResult r0 = TF::solve_throughflow(c0);
    REQUIRE(r0.status == TF::ThroughflowResult::Status::Converged);
    // Eq. 7-5: the blockage raises the mean mass flux by 1 / K_B (K_W differs a little between the two)
    std::printf("[coupling] K_B effect, station: K_B, int(rho W_m) ratio, 1 / K_B\n");
    double worst = 0.0;
    for (std::size_t k = 0; k < r.stations.size(); ++k) {
        auto I = [](const TF::ThroughflowStation& st) {
            return trapezoid(st, [](const TF::ThroughflowNode& n) { return n.K_W * n.r * n.rho * n.Wm * std::cos(n.epsilon); });
        };
        const double ratio = I(r.stations[k]) / I(r0.stations[k]);
        std::printf("[coupling]   %zu: %.5f %.6f %.6f\n", k, r.stations[k].K_B_endwall, ratio, 1.0 / r.stations[k].K_B_endwall);
        worst = std::max(worst, std::abs(ratio * r.stations[k].K_B_endwall - 1.0));
        CHECK(r0.stations[k].K_B_endwall == 1.0);
    }
    CHECK(measure("K_B effect: max |ratio K_B - 1|", worst) <= 2e-4);
    // the mirrored stage: every angle and swirl sign-changed, the same W_m, losses, K_W and K_B
    so.endwall = true;
    so.mirror = true;
    const TF::ThroughflowResult rm = TF::solve_throughflow(coupling::synthetic_stage(so));
    REQUIRE(rm.status == TF::ThroughflowResult::Status::Converged);
    double dw = 0.0, db = 0.0, dl = 0.0, dk = 0.0;
    for (std::size_t k = 0; k < r.stations.size(); ++k) {
        dk = std::max(dk, std::abs(rm.stations[k].K_B_endwall - r.stations[k].K_B_endwall));
        for (std::size_t j = 0; j < r.stations[k].nodes.size(); ++j) {
            const auto &a = r.stations[k].nodes[j], &b = rm.stations[k].nodes[j];
            dw = std::max(dw, std::abs(b.Wm / a.Wm - 1.0));
            db = std::max(db, std::abs(b.beta + a.beta) + std::abs(b.alpha + a.alpha));
            dl = std::max(dl, std::abs(b.loss_coefficient - a.loss_coefficient) + std::abs(b.K_W - a.K_W));
        }
    }
    CHECK(measure("mirrored stage: max |W_m' / W_m - 1|", dw) <= 1e-12);
    CHECK(measure("mirrored stage: max |beta' + beta| + |alpha' + alpha| [rad]", db) <= 1e-12);
    CHECK(measure("mirrored stage: max |d loss| + |d K_W|", dl) <= 1e-12);
    CHECK(measure("mirrored stage: max |d K_B|", dk) <= 1e-12);
    CHECK(rm.stations[3].nodes[5].alpha > 0.0);  // the mirrored rotor swirls toward +theta
}

TEST_CASE("coupling: BladeGen's default rotor with Ch. 6 and Ch. 8 (simple, approximate, full)", "[turbolenza][coupling]")
{
    for (const auto model : {TF::ThroughflowModel::SimpleNonIsentropic, TF::ThroughflowModel::Approximate,
                             TF::ThroughflowModel::FullNormalEquilibrium}) {
        int qb = 0;
        const TF::ThroughflowCase c = coupling::bladegen_rotor_case(model, &qb);
        const TF::ThroughflowResult r = TF::solve_throughflow(c);
        INFO(r.message);
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        CHECK(max_mass_error(r) <= 1e-11);
        const auto& a = r.stations[static_cast<std::size_t>(qb)];
        const auto& b = r.stations[static_cast<std::size_t>(qb + 1)];
        const double cp = c.gas.cp(), w = c.rows[0].omega;
        double euler = 0.0;
        for (std::size_t j = 0; j < a.nodes.size(); ++j) {
            const double work = w * (b.nodes[j].rCtheta - a.nodes[j].rCtheta);
            euler = std::max(euler, std::abs(cp * (b.nodes[j].Tt - a.nodes[j].Tt) / work - 1.0));
            CHECK(work > 0.0);
            CHECK(b.nodes[j].s > a.nodes[j].s);
        }
        CHECK(euler <= 1e-12);
        const auto& m1 = a.nodes[5];
        const auto& m2 = b.nodes[5];
        std::printf("[coupling] BladeGen rotor, %s: %zu passes; mid surface beta'1 %.2f -> beta'2 %.2f deg (kappa2 %.2f), "
                    "omega %.4f (tip gap part %.4f at the tip), K_W %.4f, Pt2/Pt1 %.4f; K_B before / after the row %.4f / %.4f; "
                    "Euler vs cp dT0 %.1e\n",
                    TF::to_string(model).c_str(), r.history.size(), std::atan2(m1.Ctheta - w * m1.r, m1.Wm) / kDeg,
                    m2.beta / kDeg, c.rows[0].sections[2].kappa2 / kDeg, m2.loss_coefficient,
                    b.nodes.back().row_loss_coefficient, m2.K_W, m2.Pt / m1.Pt, a.K_B_endwall, b.K_B_endwall, euler);
        CHECK(m2.beta > c.rows[0].sections[2].kappa2);  // Ch. 6: the flow leaves short of the metal (delta > 0)
        // across this rotor B falls (0.0204 -> 0.017 measured: the contracting, accelerating channel); not bounded here
        CHECK(b.K_B_endwall < 1.0);
        CHECK(b.K_B_endwall > 0.9);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// the wake blockage K_W (Eqs. 6-83, 6-84) and its use in Eq. 7-5
// ─────────────────────────────────────────────────────────────────────────────

namespace {

/// A test wrapper: the adapter evaluated at a FROZEN discharge estimate W_m2 = W_m1 (rho2 = rho1), so that K_W no
/// longer follows the iterate.  With D_eq > 2 the plain substitution of K_W against W_m2 has a loop gain near -1
/// (coupling.md, open items) and does not settle; frozen, the case converges and K_W's way into Eq. 7-5 can be
/// checked.  Records the last output and D_eq per surface.
class FrozenDischarge final : public TF::BladeRowModel
{
public:
    FrozenDischarge(std::shared_ptr<const TF::CorrelationBladeRow> m, int n)
        : m_(std::move(m)), kw(static_cast<std::size_t>(n)), deq(static_cast<std::size_t>(n))
    {
    }
    TF::BladeRowDischarge evaluate(const TF::BladeRowInlet& in) const override
    {
        TF::BladeRowInlet f = in;
        f.Wm2 = in.Wm1;
        f.rho2 = in.rho1;
        const T::CascadePerformance p = m_->performance(f);
        const TF::BladeRowDischarge d = m_->evaluate(f);
        kw[static_cast<std::size_t>(in.stream_surface)] = d.wake_blockage;
        deq[static_cast<std::size_t>(in.stream_surface)] = p.d_eq;
        return d;
    }
    std::shared_ptr<const TF::CorrelationBladeRow> m_;
    mutable std::vector<double> kw, deq;
};

}  // namespace

TEST_CASE("coupling: the wake blockage K_W = (2 / D_eq)^0.9 (Eq. 6-84) through the adapter and into Eq. 7-5",
          "[turbolenza][coupling]")
{
    // by hand: the rotor surface with a strong deceleration W_m 101 -> 55 m/s, D_eq > 2
    {
        const TF::CorrelationBladeRow row(rotor_spec());
        TF::BladeRowInlet in = rotor_inlet();
        in.Wm2 = 55.0;
        const T::CascadePerformance p = row.performance(in);
        REQUIRE(p.d_eq > 2.0);
        const TF::BladeRowDischarge d = row.evaluate(in);
        std::printf("[coupling] K_W by hand: D_eq %.4f, K_W %.5f\n", p.d_eq, d.wake_blockage);
        CHECK(std::abs(d.wake_blockage - std::pow(2.0 / p.d_eq, 0.9)) <= 1e-14);  // 1 - B_wake, Eqs. 6-84, 7-6
        CHECK(d.wake_blockage < 0.95);
    }
    // a single stator with alpha_1 = -40 deg and 22 vanes (D_eq above 2 over most of the span), the discharge frozen
    TF::ThroughflowCase c;
    c.model = TF::ThroughflowModel::Approximate;
    c.mass_flow = 22.0;
    for (double z : {0.0, 0.06, 0.12, 0.18}) {
        TF::QuasiNormal q;
        q.z_hub = q.z_shroud = z;
        q.r_hub = kRh;
        q.r_shroud = kRs;
        c.quasi_normals.push_back(q);
    }
    c.inlet = TF::InletProfile::uniform(101325.0, 288.15, TF::InletProfile::Swirl::Alpha, -40.0 * kDeg);
    TF::BladeRowSpec sta = stage_row("stator", 0.06, 0.12, 22, -1, 0.0);
    sta.qn_inlet = 1;
    sta.sections = {section(0.0, -40.0, 0.0, 0.044, 0.08, 22), section(1.0, -40.0, 0.0, 0.044, 0.08, 22)};
    auto frozen = std::make_shared<FrozenDischarge>(std::make_shared<TF::CorrelationBladeRow>(sta), c.options.n_stream_surfaces);
    sta.model = frozen;
    c.rows = {sta};
    const TF::ThroughflowResult r = TF::solve_throughflow(c);
    INFO(r.message);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    const auto& q2 = r.stations[2];
    double kmin = 1.0;
    for (std::size_t j = 0; j < q2.nodes.size(); ++j) {
        CHECK(q2.nodes[j].K_W == frozen->kw[j]);  // the adapter's K_W is the node's, on the discharge quasi-normal only
        if (frozen->deq[j] > 2.0) CHECK(std::abs(frozen->kw[j] - std::pow(2.0 / frozen->deq[j], 0.9)) <= 1e-14);
        else CHECK(frozen->kw[j] == 1.0);  // Eq. 6-83
        CHECK(r.stations[1].nodes[j].K_W == 1.0);
        CHECK(r.stations[3].nodes[j].K_W == 1.0);
        kmin = std::min(kmin, q2.nodes[j].K_W);
    }
    CHECK(measure("K_W case: the smallest K_W on the stator's discharge quasi-normal", kmin) < 0.95);
    // Eq. 7-5 on QN2: mdot = 2 pi K_B int K_W r rho W_m cos(eps) dy; without K_W the node data carry more
    const double with = 2.0 * kPi * trapezoid(q2, [](const TF::ThroughflowNode& n) { return n.K_W * n.r * n.rho * n.Wm * std::cos(n.epsilon); });
    const double without = 2.0 * kPi * trapezoid(q2, [](const TF::ThroughflowNode& n) { return n.r * n.rho * n.Wm * std::cos(n.epsilon); });
    CHECK(measure("K_W case: Eq. 7-5 with K_W from the nodes, |mdot / 22 - 1|", std::abs(with / 22.0 - 1.0)) <= 1e-4);
    CHECK(measure("K_W case: the same without K_W", std::abs(without / 22.0 - 1.0)) > 1e-2);
}

// ─────────────────────────────────────────────────────────────────────────────
// Chapter 8 on sloped and curved walls
// ─────────────────────────────────────────────────────────────────────────────

TEST_CASE("coupling: the Ch. 8 adapter on a 15-degree cone (sloped walls, swirl) against Ch. 8's march",
          "[turbolenza][coupling]")
{
    // the cone without its rotor (Ch. 8's row step does not converge on this rotor: 'the blade force defect did
    // not converge in 200 iterations', measured), the inlet swirl kept: sin(phi) acts in Eq. 8-6
    TF::ThroughflowCase c = coupling::cone_rotor_case();
    c.rows.clear();
    for (auto& q : c.quasi_normals) q.hub_wall_omega = 0.0;
    c.endwall = std::make_shared<TF::ChapterEightEndwall>();
    const TF::ThroughflowResult r = TF::solve_throughflow(c);
    INFO(r.message);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    const auto ew = std::dynamic_pointer_cast<TF::ChapterEightEndwall>(c.endwall);
    REQUIRE(ew);
    // the march input with m = the straight segment lengths and phi = atan2(dr, dz) of the contours, by hand
    const WallGeometry wg = straight_walls(r);
    for (std::size_t k = 0; k < r.stations.size(); ++k) {
        CHECK(std::abs(wg.phi_hub[k] - 15.0 * kDeg) <= 1e-12);
        CHECK(std::abs(ew->records()[k].hub_wall.phi - 15.0 * kDeg) <= 1e-12);  // what the adapter handed over
        CHECK(std::abs(ew->records()[k].shroud_wall.m - wg.m_shroud[k]) <= 1e-12);
    }
    const T::EndwallResult m = T::analyze_endwalls(march_input(c, r, *ew, wg));
    double dk = 0.0;
    std::printf("[coupling] cone K_B:");
    for (std::size_t k = 0; k < r.stations.size(); ++k) {
        std::printf(" %.6f", r.stations[k].K_B_endwall);
        dk = std::max(dk, std::abs(m.stations[k].blockage.K_B - r.stations[k].K_B_endwall));
    }
    std::printf("\n");
    CHECK(measure("cone: Ch. 8 adapter vs analyze_endwalls (m, phi of the contours), max |dK_B|", dk) <= 1e-12);
}

namespace {

/// u of a BladeGen wall spline at axial position z [mm] (z(u) is monotone on every piece: bisection).
double u_at(const BladeGen::Math::BSpline2& w, double z)
{
    double a = w.t0(), b = w.t1();
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (a + b);
        if (w.eval(mid).x() < z) a = mid;
        else b = mid;
    }
    return 0.5 * (a + b);
}

/// The arc length of a wall spline from u0 to u1 (composite 5-point Gauss-Legendre on 400 pieces) [mm].
double arc(const BladeGen::Math::BSpline2& w, double u0, double u1)
{
    static const double x[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
    static const double g[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
                                0.2369268850561891};
    const int n = 400;
    double s = 0.0;
    for (int p = 0; p < n; ++p) {
        const double a = u0 + (u1 - u0) * p / n, b = u0 + (u1 - u0) * (p + 1) / n;
        for (int q = 0; q < 5; ++q) {
            const auto d = w.derivs(0.5 * (a + b) + 0.5 * (b - a) * x[q], 1);
            s += 0.5 * (b - a) * g[q] * std::hypot(d[1].x(), d[1].y());
        }
    }
    return s;
}

}  // namespace

TEST_CASE("coupling: the Ch. 8 adapter on BladeGen's curved walls against Ch. 8's march from the contours",
          "[turbolenza][coupling]")
{
    const TF::ThroughflowCase c = coupling::bladegen_rotor_case(TF::ThroughflowModel::Approximate);
    const TF::ThroughflowResult r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    const auto ew = std::dynamic_pointer_cast<TF::ChapterEightEndwall>(c.endwall);
    REQUIRE(ew);
    // the wall geometry integrated independently along BladeGen's own hub and casing splines (mm -> m)
    const BladeGen::Blade::FlowPath fp(D::BladeDesign::defaults().meridional);
    WallGeometry wg;
    for (int side = 0; side < 2; ++side) {
        const auto& w = side == 0 ? fp.hub() : fp.casing();
        const double u0 = u_at(w, 1e3 * r.stations.front().nodes.front().z);
        for (const auto& st : r.stations) {
            const double u = u_at(w, 1e3 * st.nodes.front().z);  // radial quasi-normals
            const auto d = w.derivs(u, 1);
            (side == 0 ? wg.m_hub : wg.m_shroud).push_back(1e-3 * arc(w, u0, u));
            (side == 0 ? wg.phi_hub : wg.phi_shroud).push_back(std::atan2(d[1].y(), d[1].x()));
        }
    }
    double dm = 0.0, dphi = 0.0, mlen = 0.0;
    for (std::size_t k = 0; k < r.stations.size(); ++k) {
        const auto& rec = ew->records()[k];
        dm = std::max({dm, std::abs(rec.hub_wall.m - wg.m_hub[k]), std::abs(rec.shroud_wall.m - wg.m_shroud[k])});
        dphi = std::max({dphi, std::abs(rec.hub_wall.phi - wg.phi_hub[k]), std::abs(rec.shroud_wall.phi - wg.phi_shroud[k])});
        mlen = std::max(mlen, wg.m_hub[k]);
    }
    std::printf("[coupling] BladeGen walls: hub length %.5f m over %zu stations\n", mlen, r.stations.size());
    // the through-flow's wall m is Eq. 7-58 on the quasi-normal points only (two across the row): short of the
    // contour's arc where the hub rises inside the row (measured 5.3e-4 m of 0.042 m); the end-wall march takes it
    CHECK(measure("BladeGen walls: |m handed to Ch. 8 - the contour's arc length| [m]", dm) <= 5e-3);
    CHECK(measure("BladeGen walls: |phi handed to Ch. 8 - atan2(dr, dz) of the contour| [rad]", dphi) <= 1e-9);
    const T::EndwallResult m = T::analyze_endwalls(march_input(c, r, *ew, wg));
    double dk = 0.0;
    for (std::size_t k = 0; k < r.stations.size(); ++k)
        dk = std::max(dk, std::abs(m.stations[k].blockage.K_B - r.stations[k].K_B_endwall));
    CHECK(measure("BladeGen walls: Ch. 8 adapter vs analyze_endwalls from the contours, max |dK_B|", dk) <= 4e-4);  // measured 3.8e-5: the m difference above
    CHECK(dk > 0.0);
}
