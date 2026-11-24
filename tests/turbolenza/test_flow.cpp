// Flow (port of v1 flow.py): the velocities of v1 pass 0 / 1 from v1's psi, density and
// coefficients; the inlet state in closed form and against v1; uniform flow in a channel; the
// station mass flow; the refusals.
#include "SolverRef.hpp"

#include "turbolenza/Flow.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen::Turbolenza;

TEST_CASE("flow: v1 velocities from v1's psi", "[turbolenza][flow][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    for (int k = 0; k < 2; ++k) {
        const ref::json& p = d.at("passes").at(k);
        const auto g = sref::grid(p.at("grid"));
        const GridMetrics m = sref::metrics(p.at("metrics"), g);
        const BoundaryConditions b = build_boundary_conditions(g);
        const CoefficientFields c = sref::coefficients(p.at("coefficients"));
        const FlowField f = velocity_components(ref::field(p.at("psi")), ref::field(p.at("density_resident")), c, m, b);
        const FlowField v1 = sref::flow(p.at("flow"));
        const std::string tag = " pass " + std::to_string(k);
        // measured 0 for every field (the same operations as NumPy); bound 1e-15
        CHECK(sref::measure("flow psi_m" + tag, ref::rel_diff(f.psi_m, v1.psi_m)) <= 1e-15);
        CHECK(sref::measure("flow psi_eta" + tag, ref::rel_diff(f.psi_eta, v1.psi_eta)) <= 1e-15);
        CHECK(sref::measure("flow W_m" + tag, ref::rel_diff(f.W_m, v1.W_m)) <= 1e-15);
        CHECK(sref::measure("flow W_theta" + tag, ref::rel_diff(f.W_theta, v1.W_theta)) <= 1e-15);
        CHECK(sref::measure("flow W" + tag, ref::rel_diff(f.W, v1.W)) <= 1e-15);
        const std::vector<double> mass = station_mass_flow(f, ref::field(p.at("density_resident")), c, m);
        double worst = 0.0, median_ref = ref::num(p.at("record").at("mass_flow_error_max"));
        for (double q : mass) worst = std::max(worst, std::abs(q / c.mass_flow_per_depth - 1.0));
        CHECK(sref::measure("flow mass-flow error max vs v1 record" + tag, sref::rel(worst, median_ref)) <= 1e-12);
    }
}

TEST_CASE("flow: the inlet state", "[turbolenza][flow]")
{
    FlowConditions fc;
    fc.inlet_mach = 0.65;
    fc.inlet_angle_deg = 45.0;
    const PerfectGas gas(1.4, 287.05);
    const InletState s = inlet_state(fc, gas, 0.1);
    // closed form of the perfect gas: T1 = Tt / (1 + 0.2 M^2), W1 = M1 sqrt(gamma R T1)
    const double T1 = 288.15 / (1.0 + 0.2 * 0.65 * 0.65);
    CHECK(std::abs(s.static_temperature - T1) <= 1e-12 * T1);
    CHECK(std::abs(s.velocity - 0.65 * std::sqrt(1.4 * 287.05 * T1)) <= 1e-12 * s.velocity);
    CHECK(std::abs(s.velocity / s.sound_speed - 0.65) <= 1e-15);
    CHECK(s.mass_flow_per_depth == s.static_density * s.meridional_velocity * 0.1);
    CHECK(s.dynamic_pressure == 0.5 * s.static_density * s.velocity * s.velocity);
    fc.inlet_mach = 0.0;
    CHECK(inlet_state(fc, gas, 0.1).mass_flow_per_depth == 0.0);
    fc.inlet_mach.reset();
    CHECK_THROWS_AS(inlet_state(fc, gas, 0.1), InletStateError);
    fc.inlet_mach = 0.5;
    fc.inlet_angle_deg = 90.0;
    CHECK_THROWS_AS(inlet_state(fc, gas, 0.1), InletStateError);
    fc.inlet_angle_deg = 30.0;
    CHECK_THROWS_AS(inlet_state(fc, gas, 0.0), InletStateError);
    fc.angular_speed = 100.0;
    CHECK_THROWS_AS(inlet_state(fc, gas, 0.1), InletStateError);
}

TEST_CASE("flow: the inlet state of the v1 cases", "[turbolenza][flow][v1ref]")
{
    for (const char* name : {"m065_very_coarse_12", "m030_very_coarse_converged", "fixed_exit_explicit_45x21"}) {
        auto d = ref::load(name);
        if (!d) SKIP("reference data missing");
        const sref::T::SolverConfig cfg = sref::config(*d);
        const ref::json& v1 = d->at("solution").at("inlet");
        const InletState s = inlet_state(cfg.flow, cfg.gas, ref::num(v1.at("pitch")));
        // the bisection lands on the same double (measured 0); bound 1 ulp-ish
        CHECK(sref::measure(std::string("flow inlet W1 vs v1 ") + name, sref::rel(s.velocity, ref::num(v1.at("velocity")))) <=
              1e-15);
        CHECK(sref::measure(std::string("flow inlet mdot vs v1 ") + name,
                            sref::rel(s.mass_flow_per_depth, ref::num(v1.at("mass_flow_per_depth")))) <= 1e-15);
        CHECK(sref::rel(s.dynamic_pressure, ref::num(v1.at("dynamic_pressure"))) <= 1e-15);
        CHECK(sref::rel(s.static_pressure, ref::num(v1.at("static_pressure"))) <= 1e-15);
    }
}

TEST_CASE("flow: uniform flow in a channel; the refusals", "[turbolenza][flow]")
{
    const double angle = 0.3, pitch = 0.8, c_ref = 0.1, mdot = 6.0, rho = 1.15;
    const auto g = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, angle, pitch));
    const GridMetrics m = build_metrics(g);
    const BoundaryConditions b = build_boundary_conditions(g);
    const CoefficientFields c = build_coefficients(m, Field(31, 9, rho), mdot, c_ref);
    Field psi(31, 9);
    for (int i = 0; i < 31; ++i)
        for (int j = 0; j < 9; ++j) psi(i, j) = g->eta[static_cast<std::size_t>(j)];
    const FlowField f = velocity_components(psi, Field(31, 9, rho), c, m, b);
    const double Wm = mdot / (rho * pitch * c_ref);
    for (std::size_t q = 0; q < f.W.size(); ++q) {
        CHECK(std::abs(f.W_m.v[q] - Wm) <= 1e-12 * Wm);
        CHECK(std::abs(f.W_theta.v[q] / f.W_m.v[q] - std::tan(angle)) <= 1e-11);
        CHECK(std::abs(f.W.v[q] - Wm / std::cos(angle)) <= 1e-12 * Wm);
    }
    for (double q : station_mass_flow(f, Field(31, 9, rho), c, m)) CHECK(std::abs(q / mdot - 1.0) <= 1e-13);
    // a density other than the coefficients' own
    CHECK_THROWS_AS(velocity_components(psi, Field(31, 9, 1.2), c, m, b), FlowRecoveryError);
    Field bad = psi;
    bad(3, 3) = NAN;
    CHECK_THROWS_AS(velocity_components(bad, Field(31, 9, rho), c, m, b), FlowRecoveryError);
    const BoundaryConditions other = build_boundary_conditions(std::make_shared<const Grid>(*g));
    CHECK_THROWS_AS(velocity_components(psi, Field(31, 9, rho), c, m, other), FlowRecoveryError);
}
