// Eos (port of v1 eos.py): the perfect-gas relations in closed form, the refusals, and the v1
// static state / relaxed density of pass 0 of the M1 0.65 case from v1's own speeds.
#include "SolverRef.hpp"

#include "turbolenza/Eos.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen::Turbolenza;

TEST_CASE("eos: the perfect-gas state at W satisfies its defining relations", "[turbolenza][eos]")
{
    const PerfectGas gas(1.4, 287.05);
    const TotalState total(101325.0, 288.15);
    const BladeRowReference ref = blade_row_reference(total, gas);
    CHECK(gas.cp() == 1.4 * 287.05 / (1.4 - 1.0));
    CHECK(ref.relative_total_enthalpy() == gas.cp() * 288.15);
    // W = 0: the total state itself, rho = Pt / (R Tt)
    const ScalarState s0 = static_state_at(0.0, gas, ref);
    CHECK(std::abs(s0.temperature - 288.15) <= 1e-12 * 288.15);
    CHECK(std::abs(s0.pressure - 101325.0) <= 1e-12 * 101325.0);
    CHECK(std::abs(initial_density(gas, ref) - 101325.0 / (287.05 * 288.15)) <= 1e-15);
    // any W: T = Tt - W^2 / (2 cp), the isentropic Tt / T = (1 + (g-1)/2 M^2), p = rho R T
    for (double W : {50.0, 150.0, 250.0, 400.0}) {
        const ScalarState s = static_state_at(W, gas, ref);
        const double M = W / s.sound_speed;
        CHECK(std::abs(s.temperature - (288.15 - W * W / (2.0 * gas.cp()))) <= 1e-12 * 288.15);
        CHECK(std::abs(288.15 / s.temperature - (1.0 + 0.2 * M * M)) <= 1e-13);
        CHECK(std::abs(101325.0 / s.pressure - std::pow(1.0 + 0.2 * M * M, 3.5)) <= 1e-12);
        CHECK(std::abs(s.pressure - s.density * 287.05 * s.temperature) <= 1e-10 * s.pressure);
    }
    // the field version gives the same bits as the scalar one
    Field W(2, 3);
    for (std::size_t k = 0; k < W.size(); ++k) W.v[k] = 37.5 * static_cast<double>(k);
    const StaticState f = static_state_from_velocity(W, gas, ref);
    for (std::size_t k = 0; k < W.size(); ++k) {
        const ScalarState s = static_state_at(W.v[k], gas, ref);
        CHECK(f.density.v[k] == s.density);
        CHECK(f.pressure.v[k] == s.pressure);
        CHECK(f.sound_speed.v[k] == s.sound_speed);
    }
}

TEST_CASE("eos: refusals (nothing is clipped)", "[turbolenza][eos]")
{
    CHECK_THROWS_AS(PerfectGas(1.0, 287.0), ThermodynamicStateError);
    CHECK_THROWS_AS(PerfectGas(1.4, 0.0), ThermodynamicStateError);
    CHECK_THROWS_AS(TotalState(0.0, 288.0), ThermodynamicStateError);
    CHECK_THROWS_AS(TotalState(1e5, NAN), ThermodynamicStateError);
    const PerfectGas gas;
    const BladeRowReference ref = blade_row_reference(TotalState(101325.0, 288.15), gas);
    CHECK_THROWS_AS(blade_row_reference(TotalState(101325.0, 288.15), gas, 10.0), ThermodynamicStateError);
    const double W_max = std::sqrt(2.0 * ref.relative_total_enthalpy());
    CHECK_THROWS_AS(static_state_at(W_max, gas, ref), ThermodynamicStateError);
    CHECK_THROWS_AS(static_state_at(-1.0, gas, ref), ThermodynamicStateError);
    Field W(1, 3, 100.0);
    W.v[1] = NAN;
    CHECK_THROWS_AS(static_state_from_velocity(W, gas, ref), ThermodynamicStateError);
    const Field a(1, 3, 1.0), b(1, 3, 2.0);
    CHECK_THROWS_AS(relax_density(a, b, 0.0), ThermodynamicStateError);
    CHECK_THROWS_AS(relax_density(a, b, 1.5), ThermodynamicStateError);
    CHECK_THROWS_AS(relax_density(a, Field(1, 2, 1.0), 0.5), ThermodynamicStateError);
    CHECK_THROWS_AS(relax_density(a, Field(1, 3, -1.0), 0.5), ThermodynamicStateError);
    // lambda = 1 is Aungier's direct replacement, exactly the target
    CHECK(relax_density(a, b, 1.0).v == b.v);
    CHECK(relax_density(a, b, 0.5).v == Field(1, 3, 1.5).v);
}

TEST_CASE("eos: static state and relaxed density of v1 pass 0", "[turbolenza][eos][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    const ref::json& p = d.at("passes").at(0);
    const sref::T::SolverConfig cfg = sref::config(d);
    const BladeRowReference ref = blade_row_reference(
        TotalState(cfg.flow.inlet_absolute_total_pressure, cfg.flow.inlet_absolute_total_temperature), cfg.gas);
    const Field W = ref::field(p.at("flow").at("W"));
    const StaticState s = static_state_from_velocity(W, cfg.gas, ref);
    const ref::json& v1 = p.at("static_state");
    // measured 0 (bit for bit: the same UCRT pow as NumPy on Windows); bound a few ulp
    CHECK(sref::measure("eos temperature vs v1", ref::rel_diff(s.temperature, ref::field(v1.at("temperature")))) <= 1e-15);
    CHECK(sref::measure("eos pressure vs v1", ref::rel_diff(s.pressure, ref::field(v1.at("pressure")))) <= 1e-15);
    CHECK(sref::measure("eos density vs v1", ref::rel_diff(s.density, ref::field(v1.at("density")))) <= 1e-15);
    CHECK(sref::measure("eos sound speed vs v1", ref::rel_diff(s.sound_speed, ref::field(v1.at("sound_speed")))) <= 1e-15);
    const Field next = relax_density(ref::field(p.at("density_resident")), s.density, cfg.iteration.density_relaxation);
    CHECK(sref::measure("eos density_next vs v1", ref::rel_diff(next, ref::field(p.at("density_next")))) <= 1e-15);
    CHECK(sref::measure("eos initial density vs v1",
                        sref::rel(initial_density(cfg.gas, ref), ref::field(p.at("density_resident")).v[0])) <= 1e-15);
}
