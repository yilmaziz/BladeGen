// Postprocess (port of v1 postprocess.py): the blade surfaces, pressure coefficients and summary
// of the fixed-angle 45 x 21 solution against v1's; the chordwise coordinate and the relation of
// the two pressure coefficients in closed form; the refusals.
#include "SolverRef.hpp"

#include "turbolenza/Postprocess.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen::Turbolenza;

TEST_CASE("postprocess: surfaces and summary of the fixed-angle solution against v1", "[turbolenza][postprocess][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "fixed_exit_explicit_45x21");
    const Solution s = solve_at_exit_angle(sref::grid(d.at("first_grid")), sref::config(d),
                                           ref::num(d.at("settings").at("chord")));
    const BladeSurfaces b = blade_surfaces(s);
    const ref::json& v = d.at("surfaces");
    CHECK(b.ss_j == v.at("ss_j").get<int>());
    CHECK(b.ps_j == v.at("ps_j").get<int>());
    // x/c and m/c are geometry: bit for bit (measured 0); the flow values carry the solve's
    // round-off (measured below 1e-14)
    CHECK(sref::measure("postprocess x_c_ss", ref::rel_diff(b.x_c_ss, ref::vec(v.at("x_c_ss")))) <= 1e-15);
    CHECK(sref::measure("postprocess x_c_ps", ref::rel_diff(b.x_c_ps, ref::vec(v.at("x_c_ps")))) <= 1e-15);
    CHECK(sref::measure("postprocess m_c", ref::rel_diff(b.m_c, ref::vec(v.at("m_c")))) <= 1e-15);
    CHECK(sref::measure("postprocess mach_ss", ref::rel_diff(b.mach_ss, ref::vec(v.at("mach_ss")))) <= 1e-12);
    CHECK(sref::measure("postprocess mach_ps", ref::rel_diff(b.mach_ps, ref::vec(v.at("mach_ps")))) <= 5e-13);  // measured 4.7e-14
    CHECK(sref::measure("postprocess pressure_ss", ref::rel_diff(b.pressure_ss, ref::vec(v.at("pressure_ss")))) <= 1e-12);
    const SurfaceCoefficient S = surface_pressure_coefficient(s);
    const ref::json& vs = d.at("surface_pressure_coefficient");
    CHECK(S.label == vs.at("label").get<std::string>());
    CHECK(sref::measure("postprocess Dunavant S ss", ref::rel_diff(S.values_ss, ref::vec(vs.at("values_ss")))) <= 2e-13);  // measured 1.2e-14
    CHECK(sref::measure("postprocess Dunavant S ps", ref::rel_diff(S.values_ps, ref::vec(vs.at("values_ps")))) <= 1e-12);  // measured 8.4e-14
    const SolutionSummary sum = summarize(s);
    const ref::json& vm = d.at("summary");
    CHECK(sum.limitation == vm.at("limitation").get<std::string>());
    CHECK(sum.mass_flow_source == vm.at("mass_flow_source").get<std::string>());
    CHECK(sum.mach_max_node[0] == vm.at("mach_max_node")[0].get<int>());
    CHECK(sum.mach_max_node[1] == vm.at("mach_max_node")[1].get<int>());
    CHECK(sref::measure("postprocess exit velocity angle",
                        sref::rel(sum.exit_velocity_angle_deg, ref::num(vm.at("exit_velocity_angle_deg")))) <= 1e-12);
    CHECK(sum.discharge_angle_deg == ref::num(vm.at("discharge_angle_deg")));
    CHECK(sref::measure("postprocess ss peak", sref::rel(sum.ss_peak_mach, ref::num(vm.at("ss_peak_mach")))) <= 1e-12);
    CHECK(sum.ss_peak_x_c == ref::num(vm.at("ss_peak_x_c")));
    CHECK(sref::measure("postprocess pressure_min", sref::rel(sum.pressure_min, ref::num(vm.at("pressure_min")))) <= 1e-12);
    CHECK(!format_summary(sum).empty());
    // the two conventions: S = (Pt1 - p1) / q1 - Cp
    const SurfaceCoefficient Cp = surface_pressure_coefficient(s, PressureConvention::InletStatic);
    const double offset = (s.inlet->total_pressure - s.inlet->static_pressure) / s.inlet->dynamic_pressure;
    for (std::size_t k = 0; k < Cp.values_ss.size(); ++k)
        CHECK(std::abs(S.values_ss[k] - (offset - Cp.values_ss[k])) <= 1e-12);
}

TEST_CASE("postprocess: chordwise coordinate and refusals", "[turbolenza][postprocess]")
{
    const Point le{0.2, 0.1}, te{1.0, 0.7};
    const std::vector<double> x = chordwise_coordinate({le, te, {0.6, 0.4}, {0.6 - 0.3, 0.4 + 0.4}}, le, te);
    CHECK(x[0] == 0.0);
    CHECK(x[1] == 1.0);
    CHECK(std::abs(x[2] - 0.5) <= 1e-15);
    CHECK(std::abs(x[3] - 0.5) <= 1e-15);  // a point off the chord line, normal to it
    CHECK_THROWS_AS(pressure_coefficient_values({1.0}, std::nullopt, PressureConvention::Dunavant), PostprocessError);
    InletState zero;
    CHECK_THROWS_AS(pressure_coefficient_values({1.0}, zero, PressureConvention::Dunavant), PostprocessError);
    // a wall whose surface changes along the blade cannot be labeled
    Solution s;
    Grid g = sref::channel(31, 9, 8, 20, 0.2);
    g.wall_identity.lower.surface[3] = Surface::Pressure;
    s.grid = std::make_shared<const Grid>(g);
    CHECK_THROWS_AS(blade_surfaces(s), PostprocessError);
}
