// Coefficients (port of v1 coefficients.py): A..E of v1 pass 0 and pass 1 from v1's metrics and
// resident density; the defining equations (5-22 .. 5-26) on a channel; the scaling with the
// chord; the refusals.
#include "SolverRef.hpp"

#include "turbolenza/Coefficients.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen::Turbolenza;

TEST_CASE("coefficients: v1 A..E from v1 metrics and density", "[turbolenza][coefficients][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    for (int k = 0; k < 2; ++k) {
        const ref::json& p = d.at("passes").at(k);
        const auto g = sref::grid(p.at("grid"));
        const GridMetrics metrics = sref::metrics(p.at("metrics"), g);
        const CoefficientFields v1 = sref::coefficients(p.at("coefficients"));
        const CoefficientFields c = build_coefficients(metrics, ref::field(p.at("density_resident")),
                                                       v1.mass_flow_per_depth, v1.c_ref);
        const std::string tag = " pass " + std::to_string(k);
        // measured 0 for all five (the same operations in the same order as NumPy); bound 1e-15
        CHECK(sref::measure("coefficients A" + tag, ref::rel_diff(c.A, v1.A)) <= 1e-15);
        CHECK(sref::measure("coefficients B" + tag, ref::rel_diff(c.B, v1.B)) <= 1e-15);
        CHECK(sref::measure("coefficients C" + tag, ref::rel_diff(c.C, v1.C)) <= 1e-15);
        CHECK(sref::measure("coefficients D" + tag, ref::rel_diff(c.D, v1.D)) <= 1e-15);
        CHECK(sref::measure("coefficients E" + tag, ref::rel_diff(c.E, v1.E)) <= 1e-15);
    }
}

TEST_CASE("coefficients: Eqs. 5-22 to 5-26 on a channel, and the chord scaling", "[turbolenza][coefficients]")
{
    const auto g = std::make_shared<const Grid>(sref::channel(31, 7, 8, 20, 0.0, 0.7));
    const GridMetrics m = build_metrics(g);
    const Field rho(31, 7, 1.2);
    const double mdot = 10.0;
    const CoefficientFields c = build_coefficients(m, rho, mdot, 0.1);
    // beta = 0, S = 0.07 m: A = mdot / (S rho), B = 0, C = mdot S / rho, D = E = 0 exactly
    for (std::size_t q = 0; q < c.A.size(); ++q) {
        CHECK(std::abs(c.A.v[q] - mdot / (0.07 * 1.2)) <= 1e-12 * c.A.v[q]);
        // B = 0 up to the round-off of the one-sided slopes of the constant side lines (v1 identical)
        CHECK(std::abs(c.B.v[q]) <= 1e-13 * mdot / 1.2);
        CHECK(std::abs(c.C.v[q] - mdot * 0.07 / 1.2) <= 1e-12 * c.C.v[q]);
        CHECK(std::abs(c.D.v[q]) <= 1e-9 * c.A.v[q]);
        CHECK(std::abs(c.E.v[q]) <= 1e-9 * c.C.v[q]);
    }
    // fixed mdot', chord x lambda: A ~ 1/lambda, B ~ 1, C ~ lambda (coefficients.py "Scaling")
    const auto sheared = std::make_shared<const Grid>(sref::channel(31, 7, 8, 20, 0.3, 0.7));
    const GridMetrics ms = build_metrics(sheared);
    Field r(31, 7);
    for (int i = 0; i < 31; ++i)
        for (int j = 0; j < 7; ++j) r(i, j) = 1.0 + 0.01 * i + (i >= 8 && i <= 20 ? 0.02 * j : 0.0);
    const CoefficientFields a = build_coefficients(ms, r, mdot, 0.1), b = build_coefficients(ms, r, mdot, 0.3);
    for (std::size_t q = 0; q < a.A.size(); ++q) {
        CHECK(std::abs(b.A.v[q] * 3.0 - a.A.v[q]) <= 1e-13 * a.A.v[q]);
        CHECK(b.B.v[q] == a.B.v[q]);
        CHECK(std::abs(b.C.v[q] - 3.0 * a.C.v[q]) <= 1e-13 * std::abs(b.C.v[q]));
        CHECK(std::abs(b.D.v[q] * 3.0 - a.D.v[q]) <= 1e-12 * (std::abs(a.D.v[q]) + a.A.v[q]));
    }
}

TEST_CASE("coefficients: refusals", "[turbolenza][coefficients]")
{
    const auto g = std::make_shared<const Grid>(sref::channel(31, 7, 8, 20, 0.2));
    const GridMetrics m = build_metrics(g);
    CHECK_THROWS_AS(build_coefficients(m, Field(31, 7, 1.0), 0.0), CoefficientError);
    CHECK_THROWS_AS(build_coefficients(m, Field(31, 7, 1.0), 1.0, -0.1), CoefficientError);
    CHECK_THROWS_AS(build_coefficients(m, Field(31, 6, 1.0), 1.0), CoefficientError);
    Field neg(31, 7, 1.0);
    neg(10, 3) = -1.0;
    CHECK_THROWS_AS(build_coefficients(m, neg, 1.0), CoefficientError);
    Field jump(31, 7, 1.0);
    jump(2, 6) = 1.001;  // not periodic on the inlet station 2
    CHECK_THROWS_AS(build_coefficients(m, jump, 1.0), CoefficientError);
    Field blade(31, 7, 1.0);
    blade(10, 6) = 1.001;  // a blade station may differ across the passage
    CHECK_NOTHROW(build_coefficients(m, blade, 1.0));
    Grid scaled = sref::channel(31, 7, 8, 20, 0.2);
    scaled.boundaries.chord = 0.1;  // a metre grid: the chord would be applied twice
    const GridMetrics ms = build_metrics(std::make_shared<const Grid>(scaled));
    CHECK_THROWS_AS(build_coefficients(ms, Field(31, 7, 1.0), 1.0), CoefficientError);
}
