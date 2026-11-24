// ThroughflowSheet.hpp: a Chapter 7 stream surface as the blade-to-blade solver's stream sheet (Aungier 2003
// Ch. 12 p. 287, the simple one-way approach): r(m) and b(m) against the through-flow at its quasi-normals and
// linear between them, the LE / TE crossings, the solver's StreamSheet, the inlet state on the surface, the case
// settings (motion, frames, a mirrored stator), and one blade-to-blade solve on BladeGen's default rotor whose
// exit angle is compared (measured, not bounded) with Chapter 6's discharge angle.  Every measured number is
// printed ("MEASURE ..." or "[sheet] ...").
#include "CascadeFixture.hpp"
#include "CouplingFixture.hpp"

#include "turbolenza/BladeGenSource.hpp"
#include "turbolenza/Case.hpp"
#include "turbolenza/ThroughflowSheet.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>

namespace TF = BladeGen::Turbolenza::Throughflow;
namespace T = BladeGen::Turbolenza;
namespace D = BladeGen::Design;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kDeg = kPi / 180.0;

double measure(const char* name, double v)
{
    std::printf("MEASURE %-60s %.3e\n", name, v);
    return v;
}

const TF::ThroughflowCase& stage()
{
    static const TF::ThroughflowCase c = coupling::synthetic_stage();
    return c;
}

const TF::ThroughflowResult& stage_result()
{
    static const TF::ThroughflowResult r = TF::solve_throughflow(stage());
    return r;
}

T::CaseSettings b2b_settings()
{
    T::CaseSettings s;
    s.grid_level = std::string("very_coarse");
    s.density_tolerance = 1e-6;
    s.kutta_tolerance = 1e-6;
    s.density_relaxation = 0.5;
    s.kutta_relaxation = 0.5;
    s.max_density_iterations = 3000;
    s.max_kutta_iterations = 3000;
    return s;
}

/// The synthetic blade by hand (no sheet): its LE angle near the rotor's beta'_1.
T::StreamSurfaceBlade hand_blade(double stagger_deg, int Z, double chord)
{
    T::StreamSurfaceBlade b;
    b.section = fixture::synthetic_section();
    b.stagger_deg = stagger_deg;
    b.chord = chord;
    b.n_blades = Z;
    return b;
}

}  // namespace

TEST_CASE("sheet: a through-flow surface as r(m), b(m): the quasi-normals, linear between them (p. 287)", "[turbolenza][sheet][coupling]")
{
    const TF::ThroughflowCase& c = stage();
    const TF::ThroughflowResult& r = stage_result();
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    const int j = 5;
    const TF::ThroughflowSurfaceSheet s = TF::throughflow_surface_sheet(c, r, coupling::kRotor, j);
    const TF::ThroughflowStreamSurface& ss = r.stream_surfaces[static_cast<std::size_t>(j)];
    CHECK(s.qn_inlet == 2);
    CHECK(s.qn_outlet == 3);
    CHECK(s.qn_free == 2);  // the guide vane ends at QN 2
    CHECK(s.s_row == 1);
    CHECK(s.omega == c.rows[1].omega);
    // at every quasi-normal the through-flow's own r, b (the Eq. 7-32 area term), m, and the flow
    for (std::size_t k = 0; k < ss.m.size(); ++k) {
        CHECK(s.m[k] == ss.m[k]);
        CHECK(s.radius(ss.m[k]) == ss.r[k]);
        CHECK(s.thickness_at(ss.m[k]) == ss.b[k]);
        CHECK(s.Pt[k] == r.stations[k].nodes[5].Pt);
        CHECK(s.rCtheta[k] == r.stations[k].nodes[5].rCtheta);
    }
    // linear between QN1 and QN2 of the row (p. 287), and between the other quasi-normals (ours); constant beyond
    const double mh = 0.5 * (ss.m[2] + ss.m[3]);
    CHECK(std::abs(s.radius(mh) - 0.5 * (ss.r[2] + ss.r[3])) <= 1e-15);
    CHECK(std::abs(s.thickness_at(mh) - 0.5 * (ss.b[2] + ss.b[3])) <= 1e-15);
    CHECK(std::abs(s.radius(0.25 * ss.m[0] + 0.75 * ss.m[1]) - (0.25 * ss.r[0] + 0.75 * ss.r[1])) <= 1e-15);
    CHECK(s.radius(ss.m.front() - 1.0) == ss.r.front());
    CHECK(s.radius(ss.m.back() + 1.0) == ss.r.back());
    CHECK(s.thickness_at(ss.m.back() + 1.0) == ss.b.back());
    std::printf("[sheet] stage rotor surface 5: r at the quasi-normals");
    for (double x : ss.r) std::printf(" %.6f", x);
    std::printf(" m; b");
    for (double x : ss.b) std::printf(" %.6f", x);
    std::printf(" m\n");
    CHECK(std::abs(ss.r[3] - ss.r[2]) > 1e-4);  // the surface moves radially across the rotor: the rule is exercised
    // the LE / TE crossings of the radial edge lines 8 mm inside the quasi-normals: t = 8 mm / 60 mm along the segment
    const double tl = 0.008 / 0.06, tt = 0.052 / 0.06;
    CHECK(std::abs(s.z_le - 0.128) <= 1e-15);
    CHECK(std::abs(s.z_te - 0.172) <= 1e-15);
    CHECK(std::abs(s.m_le - (ss.m[2] + tl * (ss.m[3] - ss.m[2]))) <= 1e-14);
    CHECK(std::abs(s.m_te - (ss.m[2] + tt * (ss.m[3] - ss.m[2]))) <= 1e-14);
    CHECK(std::abs(s.r_le - (ss.r[2] + tl * (ss.r[3] - ss.r[2]))) <= 1e-15);
    CHECK(s.beta2 == r.stations[3].nodes[5].beta);
    // the row span: the radial edge lines of a cylindrical annulus, (r - r_h) / (r_s - r_h) at the crossings
    const double r_te = ss.r[2] + tt * (ss.r[3] - ss.r[2]);
    CHECK(std::abs(s.span - 0.5 * ((s.r_le - 0.25) + (r_te - 0.25)) / 0.1) <= 1e-13);
    CHECK(TF::throughflow_surface_sheet(c, r, coupling::kRotor, 0).span == 0.0);
    CHECK(TF::throughflow_surface_sheet(c, r, coupling::kRotor, 10).span == 1.0);
    // Eq. 12-39's thickness on request
    const TF::ThroughflowSurfaceSheet sm = TF::throughflow_surface_sheet(c, r, coupling::kRotor, j, TF::SheetThickness::Mass);
    for (std::size_t k = 0; k < ss.m.size(); ++k) CHECK(sm.b[k] == ss.b_mass[k]);
    // Eq. 12-39 is the thickness that carries the surface tube's mass: 2 pi r rho W_m K_B K_W b_mass = mdot / (n_ss - 1)
    double e1239 = 0.0;
    for (std::size_t k = 0; k < ss.m.size(); ++k) {
        const auto& nd = r.stations[k].nodes[5];
        e1239 = std::max(e1239, std::abs(2.0 * kPi * nd.r * nd.rho * nd.Wm * r.stations[k].K_B_endwall * nd.K_W * sm.b[k] /
                                             (c.mass_flow / 10.0) -
                                         1.0));
    }
    CHECK(measure("sheet: Eq. 12-39 b_mass carries the tube's mass, max rel. error", e1239) <= 1e-12);
    // the solver's sheet: M = (m - m_LE) / L, r / L, b; linear inside the row away from its two kinks
    const double L = 0.06, rref = s.r_le / L;
    const auto sh = TF::solver_stream_sheet(s, L, rref, -0.5, 2.0, 0.0025);
    CHECK(sh->m_front() == -0.5);
    CHECK(sh->m_back() == 2.0);
    CHECK(sh->r_ref() == rref);
    CHECK(std::abs(sh->radius(0.0) - s.r_le / L) <= 1e-15 * rref);  // the LE crossing is M = 0
    const double Mi = 0.5 * (s.m_te - s.m_le) / L;                  // mid-chord, 0.3 chords from the kinks at QN1 / QN2
    const double slope = (ss.r[3] - ss.r[2]) / (ss.m[3] - ss.m[2]);
    CHECK(measure("solver sheet: r inside the row against the linear rule (rel.)",
                  std::abs(sh->radius(Mi) * L / s.radius(s.m_le + Mi * L) - 1.0)) <= 1e-12);
    CHECK(measure("solver sheet: sin(phi) = dr/dm inside the row against the segment's slope",
                  std::abs(sh->sin_phi(Mi) - slope)) <= 1e-9);
    CHECK(std::abs(sh->thickness(Mi) / s.thickness_at(s.m_le + Mi * L) - 1.0) <= 1e-12);
    // refusals
    CHECK_THROWS_AS(TF::throughflow_surface_sheet(c, r, 3, j), TF::ThroughflowError);
    CHECK_THROWS_AS(TF::throughflow_surface_sheet(c, r, -1, j), TF::ThroughflowError);
    CHECK_THROWS_AS(TF::throughflow_surface_sheet(c, r, 1, 11), TF::ThroughflowError);
    CHECK_THROWS_AS(TF::throughflow_surface_sheet(c, r, 1, -1), TF::ThroughflowError);
    TF::ThroughflowResult cut = r;
    cut.stations.pop_back();  // a choked pass that ended early
    CHECK_THROWS_AS(TF::throughflow_surface_sheet(c, cut, 1, j), TF::ThroughflowError);
    TF::ThroughflowCase off = c;
    off.rows[1].le_hub_z = off.rows[1].le_shroud_z = 0.11;  // the LE line before the row's QN1
    CHECK_THROWS_AS(TF::throughflow_surface_sheet(off, r, 1, j), TF::ThroughflowError);
    CHECK_THROWS_AS(TF::solver_stream_sheet(s, 0.0, rref, -0.5, 2.0), std::invalid_argument);
    CHECK_THROWS_AS(TF::solver_stream_sheet(s, L, 0.0, -0.5, 2.0), std::invalid_argument);
    CHECK_THROWS_AS(TF::solver_stream_sheet(s, L, rref, 2.0, -0.5), std::invalid_argument);
    CHECK_THROWS_AS(TF::solver_stream_sheet(s, L, rref, -0.5, 2.0, 0.0), std::invalid_argument);
}

TEST_CASE("sheet: the inlet state on the surface (the blade-free stretch, Eqs. 7-40 .. 7-42)", "[turbolenza][sheet][coupling]")
{
    const TF::ThroughflowCase& c = stage();
    const TF::ThroughflowResult& r = stage_result();
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    // the guide vane's surface: the stretch QN 0 .. QN 1
    const TF::ThroughflowSurfaceSheet s = TF::throughflow_surface_sheet(c, r, coupling::kIgv, 3);
    CHECK(s.qn_free == 0);
    for (int k = 0; k <= 1; ++k) {
        const auto& nd = r.stations[static_cast<std::size_t>(k)].nodes[3];
        const TF::SurfaceInletState in = TF::surface_inlet_state(s, c.gas, s.m[static_cast<std::size_t>(k)]);
        CHECK(std::abs(in.mach - nd.mach) <= 1e-13);
        CHECK(std::abs(in.alpha - nd.alpha) <= 1e-15);
        CHECK(std::abs(in.T - nd.T) <= 1e-11);
        CHECK(in.Pt == nd.Pt);
    }
    // halfway: P_t, T_t, W_m and r C_theta linear in m, C_theta = r C_theta / r
    const double mh = 0.5 * (s.m[0] + s.m[1]);
    const TF::SurfaceInletState h = TF::surface_inlet_state(s, c.gas, mh);
    CHECK(std::abs(h.Wm - 0.5 * (s.Wm[0] + s.Wm[1])) <= 1e-13);
    CHECK(std::abs(h.Ctheta - 0.5 * (s.rCtheta[0] + s.rCtheta[1]) / s.radius(mh)) <= 1e-13);
    // before QN 0: QN 0's state; between QN1 and the LE: QN1's (with C_theta = (r C_theta)_1 / r there)
    CHECK(TF::surface_inlet_state(s, c.gas, s.m[0] - 0.5).Wm == s.Wm[0]);
    const TF::SurfaceInletState b = TF::surface_inlet_state(s, c.gas, 0.5 * (s.m[1] + s.m_le));
    CHECK(b.Wm == s.Wm[1]);
    CHECK(b.Pt == s.Pt[1]);
    // the rotor's surface: the stretch is QN 2 alone; a far inlet inside the guide vane is refused
    const TF::ThroughflowSurfaceSheet sr = TF::throughflow_surface_sheet(c, r, coupling::kRotor, 3);
    const TF::SurfaceInletState q2 = TF::surface_inlet_state(sr, c.gas, sr.m[2]);
    CHECK(std::abs(q2.alpha - r.stations[2].nodes[3].alpha) <= 1e-15);
    CHECK(q2.alpha < -5.0 * kDeg);  // the guide vane's swirl toward -theta
    CHECK_THROWS_AS(TF::surface_inlet_state(sr, c.gas, sr.m[2] - 0.01), TF::ThroughflowError);
    TF::ThroughflowSurfaceSheet cold = s;
    for (double& t : cold.Tt) t = 1.0;  // no static state
    CHECK_THROWS_AS(TF::surface_inlet_state(cold, c.gas, s.m[0]), TF::ThroughflowError);
}

TEST_CASE("sheet: the blade-to-blade case on the stage's rotor and (mirrored) stator surfaces", "[turbolenza][sheet][coupling]")
{
    const TF::ThroughflowCase& c = stage();
    const TF::ThroughflowResult& r = stage_result();
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    T::CaseSettings base = b2b_settings();
    // the rows are packed: the far inlet of a 1-chord extension lies in the guide vane, refused
    const double L = 0.06;
    CHECK_THROWS_AS(TF::blade_to_blade_on_surface(c, r, coupling::kRotor, 5, hand_blade(37.0, 32, L), base), TF::ThroughflowError);
    base.upstream_extension_chords = 0.12;
    base.downstream_extension_chords = 0.12;
    const TF::BladeToBladeOnSurface o = TF::blade_to_blade_on_surface(c, r, coupling::kRotor, 5, hand_blade(37.0, 32, L), base);
    CHECK(o.frame_sign == 1);
    CHECK(o.settings.stream_surface == T::StreamSurfaceKind::Layer);
    CHECK(o.settings.row == T::RowMotion::Rotating);
    CHECK(std::abs(o.settings.rpm / 5000.0 - 1.0) <= 1e-14);  // omega = -2 pi rpm / 60 and back
    CHECK(o.settings.inlet_total_pressure == o.inlet.Pt);
    CHECK(o.settings.inlet_absolute_mach == o.inlet.mach);
    CHECK(o.settings.inlet_absolute_flow_angle_deg == o.inlet.alpha / kDeg);
    CHECK(o.settings.gamma == c.gas.gamma);
    // a blade without a sheet: its pitch at the through-flow's r_LE
    CHECK(std::abs(o.blade.sheet->r_ref() - o.sheet.r_le / L) <= 1e-15);
    CHECK(std::abs(o.blade.solidity - L * 32.0 / (2.0 * kPi * o.sheet.r_le)) <= 1e-14);
    CHECK(o.r_le_ratio == 1.0);
    CHECK(o.blade.sheet->m_front() <= -base.upstream_extension_chords - 1.0 + 1e-12);
    const T::Case cs = T::build_case(o.settings, o.blade);
    REQUIRE(cs.triangle);
    // the case's relative inlet: beta'_1 at the far-inlet radius from the through-flow's QN1 state
    const double rin = cs.triangle->radius, U = -cs.triangle->omega * rin;
    const double beta = std::atan2(o.inlet.Ctheta + U, o.inlet.Wm);  // the triangle of the absolute state at r_in
    std::printf("[sheet] stage rotor surface 5: case beta1 %.4f deg (through-flow QN1 %.4f deg), r_in %.5f m (QN1 %.5f m)\n",
                cs.inlet_angle / kDeg, std::atan2(r.stations[2].nodes[5].Ctheta - c.rows[1].omega * r.stations[2].nodes[5].r,
                                                  r.stations[2].nodes[5].Wm) / kDeg,
                rin, r.stations[2].nodes[5].r);
    CHECK(std::abs(cs.inlet_angle - beta) <= 1e-12);
    // the stator: s_row = -1, stationary, mirrored (theta -> -theta): the case sees a positive inlet angle
    const TF::BladeToBladeOnSurface st = TF::blade_to_blade_on_surface(c, r, coupling::kStator, 5, hand_blade(11.0, 60, 0.046), base);
    CHECK(st.frame_sign == -1);
    CHECK(st.settings.row == T::RowMotion::Stationary);
    CHECK(st.settings.rpm == 0.0);
    CHECK(st.inlet.alpha < 0.0);
    CHECK(st.settings.inlet_flow_angle_deg == -st.inlet.alpha / kDeg);
    CHECK(st.settings.inlet_mach == st.inlet.mach);
    // the far inlet lies beyond the stator's QN1: QN1's state, C_theta = (r C_theta)_1 / r there
    CHECK(st.inlet.m > st.sheet.m[3]);
    CHECK(st.inlet.Wm == r.stations[3].nodes[5].Wm);
    CHECK(st.inlet.Pt == r.stations[3].nodes[5].Pt);
    CHECK(std::abs(st.inlet.Ctheta * st.inlet.r - r.stations[3].nodes[5].rCtheta) <= 1e-14);
    const T::Case cst = T::build_case(st.settings, st.blade);
    CHECK(cst.inlet_angle > 0.0);
    // the mirrored stage (theta -> -theta): its rotor has s_row = -1 and runs mirrored, the very case of the stage
    {
        coupling::StageOptions mo;
        mo.mirror = true;
        const TF::ThroughflowCase cm = coupling::synthetic_stage(mo);
        const TF::ThroughflowResult rm = TF::solve_throughflow(cm);
        REQUIRE(rm.status == TF::ThroughflowResult::Status::Converged);
        const TF::BladeToBladeOnSurface om = TF::blade_to_blade_on_surface(cm, rm, coupling::kRotor, 5, hand_blade(37.0, 32, L), base);
        CHECK(om.frame_sign == -1);
        CHECK(om.inlet.alpha == -o.inlet.alpha);  // the book frame's
        CHECK(om.settings.rpm == o.settings.rpm);
        CHECK(om.settings.inlet_absolute_flow_angle_deg == o.settings.inlet_absolute_flow_angle_deg);
        const T::Case csm = T::build_case(om.settings, om.blade);
        CHECK(csm.inlet_angle == cs.inlet_angle);
    }
    // refusals
    CHECK_THROWS_AS(TF::blade_to_blade_on_surface(c, r, coupling::kRotor, 5, hand_blade(37.0, 0, L), base), std::invalid_argument);
    CHECK_THROWS_AS(TF::blade_to_blade_on_surface(c, r, coupling::kRotor, 5, hand_blade(37.0, 32, 0.0), base), std::invalid_argument);
}

TEST_CASE("sheet: blade-to-blade on BladeGen's default rotor (mid surface) against Chapter 6 (measured)", "[turbolenza][sheet][coupling][bladegen]")
{
    int qb = 0;
    const TF::ThroughflowCase c = coupling::bladegen_rotor_case(TF::ThroughflowModel::Approximate, &qb);
    const TF::ThroughflowResult r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    D::BladeDesign design = D::BladeDesign::defaults();
    D::Evaluator ev;
    const T::StreamSurfaceBlade blade = T::stream_surface_blade(design, 2, design.layers[2], ev);
    // the surface whose row span is nearest the Mid layer's (0.5): not psi = 0.5, which lies at a larger radius
    int j = 0;
    for (int k = 1; k + 1 < c.options.n_stream_surfaces; ++k)
        if (std::abs(TF::throughflow_surface_sheet(c, r, 0, k).span - design.layers[2].span) <
            std::abs(TF::throughflow_surface_sheet(c, r, 0, j).span - design.layers[2].span))
            j = k;
    const TF::BladeToBladeOnSurface o = TF::blade_to_blade_on_surface(c, r, 0, j, blade, b2b_settings());
    std::printf("[sheet] BladeGen rotor surface %d (psi %.2f, span %.4f): r_LE through-flow / layer = %.6f, b(TE)/b(LE) %.5f (layer's own %.5f), "
                "far inlet alpha %.3f deg M %.4f\n",
                j, o.sheet.mass_fraction, o.sheet.span, o.r_le_ratio, o.blade.sheet->thickness((o.sheet.m_te - o.sheet.m_le) / blade.chord) /
                                                            o.blade.sheet->thickness(0.0),
                blade.sheet->thickness(std::cos(blade.stagger_deg * kDeg)) / blade.sheet->thickness(0.0), o.inlet.alpha / kDeg,
                o.inlet.mach);
    CHECK(std::abs(o.r_le_ratio - 1.0) <= 0.05);  // the layer and the through-flow surface nearest its span (placement mismatch)
    CHECK(std::abs(o.settings.rpm / 15000.0 - 1.0) <= 1e-14);
    const T::Case cs = T::build_case(o.settings, o.blade);
    const T::RunResult rr = T::run_case(cs);
    INFO(rr.error);
    REQUIRE(rr.converged);
    REQUIRE(rr.solution);
    const T::SolutionData& d = *rr.solution;
    // the solve's stream sheet is the through-flow's: r at every station from the sheet
    double dr = 0.0;
    for (int i = 0; i < d.grid.n_m(); ++i) {
        const double M = d.grid.m[static_cast<std::size_t>(i)];
        dr = std::max(dr, std::abs(d.station_radius[static_cast<std::size_t>(i)] / (o.sheet.radius(o.sheet.m_le + M * blade.chord)) - 1.0));
    }
    CHECK(measure("B2B station radius against the through-flow's piecewise-linear r(m)", dr) <= 1e-6);
    // the measured comparison: the blade-to-blade exit angle (Kutta, far exit) against Ch. 6's beta'_2 on QN2
    const double b2b = static_cast<double>(o.frame_sign) * d.exit_angle;
    const auto& n1 = r.stations[static_cast<std::size_t>(qb)].nodes[static_cast<std::size_t>(j)];
    const auto& n2 = r.stations[static_cast<std::size_t>(qb + 1)].nodes[static_cast<std::size_t>(j)];
    const double w = c.rows[0].omega;
    std::printf("[sheet] BladeGen rotor mid surface: beta'1 through-flow %.3f deg, B2B %.3f deg; beta'2 Ch. 6 %.3f deg, B2B exit "
                "%.3f deg (difference %.3f deg); kappa_TE %.3f deg; %d passes\n",
                std::atan2(n1.Ctheta - w * n1.r, n1.Wm) / kDeg, cs.inlet_angle / kDeg, o.sheet.beta2 / kDeg, b2b / kDeg,
                (b2b - o.sheet.beta2) / kDeg, cs.cascade.kappa_te / kDeg, rr.iterations);
    if (rr.summary && rr.summary->rotor) {
        const T::RotorSummary& rs = *rr.summary->rotor;
        std::printf("[sheet] BladeGen rotor mid surface: Euler work B2B %.1f J/kg (cp dT0 %.1f), through-flow %.1f J/kg; Pt ratio "
                    "B2B %.4f (inviscid), through-flow %.4f\n",
                    rs.euler_work, rs.total_enthalpy_rise, w * (n2.rCtheta - n1.rCtheta), rs.total_pressure_ratio, n2.Pt / n1.Pt);
        CHECK(std::abs(rs.euler_work / rs.total_enthalpy_rise - 1.0) <= 1e-6);
    }
    CHECK(std::isfinite(b2b));
}

TEST_CASE("sheet: the BladeGen rotor's blade-to-blade exit angle by grid level (report)", "[.b2bgrid][turbolenza][sheet]")
{
    const TF::ThroughflowCase c = coupling::bladegen_rotor_case(TF::ThroughflowModel::Approximate);
    const TF::ThroughflowResult r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    D::BladeDesign design = D::BladeDesign::defaults();
    D::Evaluator ev;
    const T::StreamSurfaceBlade blade = T::stream_surface_blade(design, 2, design.layers[2], ev);
    for (const char* level : {"very_coarse", "coarse", "medium"}) {
        T::CaseSettings s = b2b_settings();
        s.grid_level = std::string(level);
        s.max_density_iterations = 8000;
        s.max_kutta_iterations = 8000;
        const TF::BladeToBladeOnSurface o = TF::blade_to_blade_on_surface(c, r, 0, 4, blade, s);
        const T::RunResult rr = T::run_case(T::build_case(o.settings, o.blade));
        std::printf("[sheet] grid %s: converged %d, B2B exit %.3f deg against Ch. 6 %.3f deg\n", level, rr.converged,
                    rr.solution ? rr.solution->exit_angle / kDeg : 0.0, o.sheet.beta2 / kDeg);
    }
}
