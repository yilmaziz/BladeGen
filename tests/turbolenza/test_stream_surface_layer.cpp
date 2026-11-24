// A BladeGen layer on its stream surface (BladeGenSource.hpp stream_surface_blade, not in v1): the
// placed section in the solver's frame, BladeGen's blade angles, the pitch, r(m) and b(m); a cylindrical
// channel against the planar cascade of the same pitch and chord; and (hidden, [.rotor]) one rotating
// case on the default design's Mid layer.
#include "SolverRef.hpp"

#include "turbolenza/BladeGenSource.hpp"
#include "turbolenza/Case.hpp"
#include "turbolenza/Detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>

namespace T = BladeGen::Turbolenza;
namespace D = BladeGen::Design;

namespace {

constexpr double kPi = std::numbers::pi;

T::CaseSettings layer_settings()
{
    T::CaseSettings s;
    s.stream_surface = T::StreamSurfaceKind::Layer;
    s.grid_level = std::string("very_coarse");
    s.density_tolerance = 1e-6;
    s.kutta_tolerance = 1e-6;
    s.density_relaxation = 0.5;
    s.kutta_relaxation = 0.5;
    s.max_density_iterations = 3000;
    s.max_kutta_iterations = 3000;
    s.inlet_total_pressure = 101325.0;
    s.inlet_total_temperature = 288.15;
    return s;
}

}  // namespace

TEST_CASE("a BladeGen layer on its stream surface: the placed section, angles, pitch, r(m), b(m)", "[turbolenza][bladegen][sheet]")
{
    D::BladeDesign design = D::BladeDesign::defaults();
    D::Evaluator ev;
    const int k = 2;  // Mid
    T::StreamSurfaceReport rep;
    const T::StreamSurfaceBlade b = T::stream_surface_blade(design, k, design.layers[2], ev, {}, &rep);
    const D::DesignResult res = ev.evaluate(design);
    const D::LayerResult& L = res.layers[2];
    REQUIRE(L.error.empty());
    // the section: LE (0, 0), TE (1, 0), fitted within the fit's tolerance
    CHECK(std::hypot(b.section.suction.evaluate(0.0).m, b.section.suction.evaluate(0.0).y) <= 1e-15);
    CHECK(std::abs(b.section.suction.evaluate(1.0).m - 1.0) <= 1e-15);
    CHECK(b.max_fit_residual <= 1e-7);
    // our exact map against BladeGen's StreamSurface table (Hermite, 128 intervals per piece)
    CHECK(sref::measure("layer map m(m') vs BladeGen's StreamSurface [chords]", rep.map_max_difference) <= 1e-8);
    // the true blade angles: kappa_le, kappa_te are BladeGen's metal angles beta_LE, beta_TE (custom layers)
    const T::Cascade c = T::make_cascade(b);
    CHECK(std::abs(T::detail::degrees(c.kappa_le) - L.beta_le) <= 1e-9);
    CHECK(std::abs(T::detail::degrees(c.kappa_te) - L.beta_te) <= 1e-9);
    // the pitch r_LE 2 pi / N, the radius at the LE r_LE (sigma = 1 there)
    CHECK(std::abs(b.r_le - 1e-3 * L.r_le) <= 1e-12 * b.r_le);
    // the sheet's cubic spline through its stations (step 0.0025 chords) against the exact streamline
    CHECK(sref::measure("Mid layer sheet r(LE) vs the streamline (relative)",
                        std::abs(b.sheet->radius(0.0) * b.chord / b.r_le - 1.0)) <= 1e-9);
    CHECK(std::abs(b.sheet->scale(0.0) - 1.0) <= 1e-9);
    CHECK(std::abs(1.0 / b.solidity - b.sheet->r_ref() * 2.0 * kPi / design.n_blades) <= 1e-12);
    CHECK(b.n_blades == design.n_blades);
    // r at the TE point (the sheet at the TE's m = cos(stagger) chords) is BladeGen's r_TE
    const double m_te = std::cos(T::detail::radians(b.stagger_deg));
    CHECK(sref::measure("Mid layer sheet r(TE) vs BladeGen r_TE (relative)", std::abs(b.sheet->radius(m_te) * b.chord / (1e-3 * L.r_te) - 1.0)) <= 1e-9);
    // b(m) of the channel family: positive, smooth, b_te / b_le reported
    CHECK(rep.b_te > 0.0);
    sref::measure("Mid layer b_TE / b_LE", rep.b_te);
    sref::measure("Mid layer r_TE / r_LE", b.r_te / b.r_le);
    CHECK_THROWS_AS(T::stream_surface_blade(design, 7, design.layers[2], ev), std::invalid_argument);
    D::BladeDesign mirrored = design;
    mirrored.mirror = true;
    CHECK_THROWS_AS(T::stream_surface_blade(mirrored, k, design.layers[2], ev), std::invalid_argument);
}

TEST_CASE("a cylindrical channel: the layer on its stream surface = the planar cascade of the same pitch and chord", "[turbolenza][bladegen][sheet]")
{
    // straight cylindrical walls: r = r0 on the Mid streamline, b = casing - hub radius, constant
    D::BladeDesign design = D::BladeDesign::naca65_defaults();
    design.meridional.hub = {-15.0, 60.0, 15.0, 60.0};
    design.meridional.casing = {-15.0, 100.0, 15.0, 100.0};
    design.meridional.hub_law = design.meridional.casing_law = BladeGen::Blade::WallLaw::Linear;
    design.n_blades = 30;
    D::Evaluator ev;
    const D::Layer& layer = design.layers[2];
    const T::StreamSurfaceBlade b = T::stream_surface_blade(design, 2, layer, ev);
    const double r0 = 80.0;  // mid span
    for (double m : {-3.0, 0.0, 0.7, 4.0}) {
        CHECK(std::abs(b.sheet->radius(m) * b.chord / 1e-3 - r0) <= 1e-11 * r0);
        CHECK(std::abs(b.sheet->sin_phi(m)) <= 1e-12);
        CHECK(std::abs(b.sheet->thickness(m) - 40.0) <= 1e-10);
    }
    const D::DesignResult res = ev.evaluate(design);
    const D::LayerResult& L = res.layers[2];
    // the planar cascade: chord r0 c', pitch 2 pi r0 / N, BladeGen's stagger
    T::CaseSettings planar = layer_settings();
    planar.stream_surface = T::StreamSurfaceKind::Planar;
    planar.chord = 1e-3 * r0 * L.chord_mprime;
    planar.solidity = L.chord_mprime * design.n_blades / (2.0 * kPi);
    planar.stagger_deg = T::detail::degrees(L.gamma);
    planar.inlet_mach = 0.5;
    planar.inlet_flow_angle_deg = planar.stagger_deg + 10.0;
    T::CaseSettings onsheet = planar;
    onsheet.stream_surface = T::StreamSurfaceKind::Layer;
    const T::Case cp = T::build_case(planar, T::section_from_layer(layer, ev));
    const T::Case cs = T::build_case(onsheet, b);
    CHECK(sref::measure("cylinder: chord vs r0 c' (relative)", std::abs(cs.settings.chord / planar.chord - 1.0)) <= 1e-12);
    CHECK(sref::measure("cylinder: stagger vs BladeGen's [deg]", std::abs(cs.settings.stagger_deg - planar.stagger_deg)) <= 1e-10);
    CHECK(std::abs(cs.settings.solidity / planar.solidity - 1.0) <= 1e-12);
    CHECK(std::abs(cs.cascade.kappa_le - cp.cascade.kappa_le) <= 1e-12);
    CHECK(std::abs(cs.cascade.kappa_te - cp.cascade.kappa_te) <= 1e-12);
    CHECK(std::abs(cs.mass_flow_per_depth / cp.mass_flow_per_depth - 1.0) <= 1e-12);
    const T::RunResult a = T::run_case(cp), r = T::run_case(cs);
    REQUIRE(a.converged);
    REQUIRE(r.converged);
    const double d_beta = std::abs(a.solution->exit_angle - r.solution->exit_angle);
    CHECK(sref::measure("cylinder layer vs planar cascade: beta_out [rad]", d_beta) <= 1e-8);
    double dm = 0.0;
    for (std::size_t q = 0; q < a.solution->suction.mach.size() && q < r.solution->suction.mach.size(); ++q)
        dm = std::max(dm, std::abs(a.solution->suction.mach[q] - r.solution->suction.mach[q]));
    CHECK(a.solution->suction.mach.size() == r.solution->suction.mach.size());
    CHECK(sref::measure("cylinder layer vs planar cascade: SS Mach", dm) <= 1e-7);
}

TEST_CASE("a rotating row on the default design's Mid layer (report)", "[.rotor][turbolenza][bladegen][sheet]")
{
    D::BladeDesign design = D::BladeDesign::defaults();
    D::Evaluator ev;
    const T::StreamSurfaceBlade b = T::stream_surface_blade(design, 2, design.layers[2], ev);
    const T::Cascade c = T::make_cascade(b);
    T::CaseSettings s = layer_settings();
    s.grid_level = std::string("medium");
    s.density_tolerance = 1e-5;
    s.kutta_tolerance = 1e-5;
    s.density_relaxation = 0.3;
    s.kutta_relaxation = 0.3;
    s.max_density_iterations = 6000;
    s.max_kutta_iterations = 6000;
    s.row = T::RowMotion::Rotating;
    s.inlet_absolute_mach = 0.3;
    s.inlet_absolute_flow_angle_deg = 0.0;
    // rpm for an incidence of 2 deg at the nominal far inlet: tan(beta1) = Omega r_in / C_m
    const double r_in = b.chord * b.sheet->radius(-s.upstream_extension_chords);
    const double T1 = s.inlet_total_temperature / (1.0 + 0.2 * 0.09);
    const double C = 0.3 * std::sqrt(1.4 * 287.05 * T1);
    const double beta = c.kappa_le + T::detail::radians(2.0);
    s.rpm = std::round(C * std::tan(beta) / r_in * 60.0 / (2.0 * kPi));
    const T::Case cs = T::build_case(s, b);
    std::printf("%s\n%s\n", T::format_case_summary(cs).c_str(), T::format_grid_summary(cs).c_str());
    const T::RunResult r = T::run_case(cs);
    std::printf("%s\n", T::format_run_report(r).c_str());
    CHECK(r.converged);
    REQUIRE(r.summary);
    REQUIRE(r.summary->rotor);
    const T::RotorSummary& o = *r.summary->rotor;
    // rothalpy: cp (T0_2 - T0_1) = U2 C_theta2 - U1 C_theta1 at every node, so in the mass averages too
    CHECK(std::abs(o.total_enthalpy_rise - o.euler_work) <= 1e-6 * std::abs(o.euler_work));
    CHECK(o.total_pressure_ratio > 1.0);  // a compressor rotor
}
