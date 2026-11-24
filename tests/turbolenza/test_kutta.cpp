// Kutta (port of v1 kutta.py, production rule): v1's residual and correction from v1's flow;
// the sign of the correction, the zero step, eps_K at zero speed, the refusals.
#include "SolverRef.hpp"

#include "turbolenza/Kutta.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen::Turbolenza;

TEST_CASE("kutta: v1 residual and correction from v1's flow", "[turbolenza][kutta][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    for (int k = 0; k < 2; ++k) {
        const ref::json& p = d.at("passes").at(static_cast<std::size_t>(k));
        const auto g = sref::grid(p.at("grid"));
        const GridMetrics m = sref::metrics(p.at("metrics"), g);
        const FlowField f = sref::flow(p.at("flow"));
        const KuttaResidual r = kutta_residual(f, *g);
        const ExitAngleCorrection c = exit_angle_correction(r, f, m);
        const ref::json& vr = p.at("kutta_residual");
        const ref::json& vc = p.at("exit_angle_correction");
        const std::string tag = " pass " + std::to_string(k);
        CHECK(r.station.i == vr.at("station").at("i").get<int>());
        CHECK(r.station.ss_j == vr.at("station").at("ss_j").get<int>());
        // measured 0 (residual, eps_K) and <= 1 ulp (the trapezoidal integral); bounds 1e-15
        CHECK(sref::measure("kutta R_K" + tag, sref::rel(r.residual, ref::num(vr.at("residual")))) <= 1e-15);
        CHECK(sref::measure("kutta eps_K" + tag, sref::rel(r.normalized, ref::num(vr.at("normalized")))) <= 1e-15);
        CHECK(sref::measure("kutta exit integral" + tag,
                            sref::rel(c.exit_meridional_integral, ref::num(vc.at("exit_meridional_integral")))) <= 1e-15);
        CHECK(sref::measure("kutta delta_tan" + tag, sref::rel(c.delta_tan, ref::num(vc.at("delta_tan")))) <= 1e-15);
        const double next = next_exit_angle(g->exit_angle(), c, ref::num(d.at("config").at("iteration").at("kutta_relaxation")));
        CHECK(sref::measure("kutta next angle" + tag, sref::rel(next, ref::num(p.at("record").at("exit_angle_next")))) <=
              1e-15);
    }
}

TEST_CASE("kutta: sign, zero step, zero speed, refusals", "[turbolenza][kutta]")
{
    const auto g = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, 0.3, 0.8));
    const GridMetrics m = build_metrics(g);
    FlowField f{Field(31, 9, 0.0), Field(31, 9, 1.0), Field(31, 9, 100.0), Field(31, 9, 30.0), Field(31, 9), Field(31, 9)};
    // equal TE speeds: R_K = 0 and the angle stays bit for bit
    KuttaResidual r = kutta_residual(f, *g);
    CHECK(r.residual == 0.0);
    CHECK(r.normalized == 0.0);
    ExitAngleCorrection c = exit_angle_correction(r, f, m);
    CHECK(c.delta_tan == 0.0);
    CHECK(next_exit_angle(0.3, c, 1.0) == 0.3);
    // the eta = 0 wall (SS here) faster: Delta(tan beta_out) > 0, beta_out up (kutta.py section 5)
    f.W_m(20, 0) = 110.0;
    r = kutta_residual(f, *g);
    CHECK(r.residual > 0.0);
    CHECK(r.wall_difference > 0.0);
    c = exit_angle_correction(r, f, m);
    CHECK(c.delta_tan > 0.0);
    CHECK(next_exit_angle(0.3, c, 0.5) > 0.3);
    CHECK(next_exit_angle(0.3, c, 0.5) < next_exit_angle(0.3, c, 1.0));
    // the closed form of Eqs. 13-57, 13-59 on this channel
    const double cosb = std::cos(0.3);
    const double expected = 2.0 * r.wall_difference * g->dm() / (0.8 * 2.0 * cosb) / 100.0;
    CHECK(std::abs(c.delta_tan - expected) <= 1e-13 * expected);
    // both speeds zero: equal, eps_K = 0
    FlowField z{Field(31, 9), Field(31, 9), Field(31, 9, 0.0), Field(31, 9, 0.0), Field(31, 9), Field(31, 9)};
    CHECK(kutta_residual(z, *g).normalized == 0.0);
    CHECK_THROWS_AS(exit_angle_correction(kutta_residual(z, *g), z, m), KuttaError);  // no flow leaves the exit
    CHECK_THROWS_AS(next_exit_angle(0.3, c, 0.0), KuttaError);
    CHECK_THROWS_AS(next_exit_angle(0.3, c, 1.5), KuttaError);
    CHECK_THROWS_AS(next_exit_angle(1.6, c, 1.0), KuttaError);
    // a TE station with two suction nodes is inadmissible
    Grid bad = *g;
    bad.wall_identity.upper.surface.back() = Surface::Suction;
    CHECK_THROWS_AS(kutta_residual(f, bad), KuttaError);
}
