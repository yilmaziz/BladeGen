// Tests of math/BSpline.hpp and math/Interpolation.hpp against closed forms.
#include "math/BSpline.hpp"
#include "math/Interpolation.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <random>

using namespace BladeGen::Math;
using Catch::Approx;

namespace {
double bernstein(int n, int i, double t)
{
    double c = 1.0;
    for (int k = 1; k <= i; ++k) c = c * (n - k + 1) / k;
    return c * std::pow(t, i) * std::pow(1.0 - t, n - i);
}
}  // namespace

TEST_CASE("clamped knot vector layout", "[bspline]")
{
    const auto U = clamped_uniform_knots(3, 5);
    REQUIRE(U == std::vector<double>{0, 0, 0, 0, 0.5, 1, 1, 1, 1});
    REQUIRE_THROWS_AS(clamped_uniform_knots(3, 3), std::invalid_argument);
    REQUIRE_THROWS_AS(clamped_uniform_knots(0, 3), std::invalid_argument);
}

TEST_CASE("a single-span clamped B-spline is the Bezier curve", "[bspline]")
{
    std::vector<Vec2> P = {{0, 0}, {1, 2}, {3, 3}, {4, 0}};
    const auto c = BSpline2::clamped(3, P);
    for (double t : {0.0, 0.1, 0.37, 0.5, 0.81, 1.0}) {
        Vec2 ref = Vec2::Zero();
        for (int i = 0; i < 4; ++i) ref += bernstein(3, i, t) * P[static_cast<std::size_t>(i)];
        REQUIRE((c.eval(t) - ref).norm() < 1e-14);
    }
    // end derivatives of a Bezier curve: n (P1 - P0), n (n-1) (P2 - 2 P1 + P0)
    const auto d0 = c.derivs(0.0, 2);
    REQUIRE((d0[1] - 3.0 * (P[1] - P[0])).norm() < 1e-13);
    REQUIRE((d0[2] - 6.0 * (P[2] - 2.0 * P[1] + P[0])).norm() < 1e-12);
}

TEST_CASE("basis functions form a partition of unity and derivatives match finite differences", "[bspline]")
{
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<Vec3> P;
    for (int i = 0; i < 9; ++i) P.emplace_back(u(rng), u(rng), u(rng));
    std::vector<double> U = {0, 0, 0, 0, 0.1, 0.3, 0.3, 0.6, 0.9, 1, 1, 1, 1};
    const BSpline3 c(3, U, P);
    for (int k = 0; k <= 40; ++k) {
        const double t = k / 40.0;
        const int span = find_span(8, 3, t, U);
        const auto N = ders_basis_funs(span, t, 3, 0, U);
        double sum = 0.0;
        for (double v : N[0]) sum += v;
        REQUIRE(sum == Approx(1.0).epsilon(1e-14));
        if (t > 0.01 && t < 0.99 && std::abs(t - 0.3) > 0.01) {
            const double h = 1e-6;
            const auto d = c.derivs(t, 2);
            const Vec3 fd1 = (c.eval(t + h) - c.eval(t - h)) / (2 * h);
            const Vec3 fd2 = (c.derivs(t + h, 1)[1] - c.derivs(t - h, 1)[1]) / (2 * h);
            REQUIRE((d[1] - fd1).norm() < 1e-6 * (1 + d[1].norm()));
            REQUIRE((d[2] - fd2).norm() < 1e-4 * (1 + d[2].norm()));
        }
    }
}

TEST_CASE("reversal keeps the point set", "[bspline]")
{
    const auto c = BSpline2::clamped(3, {{0, 0}, {1, 1}, {2, -1}, {3, 0.5}, {4, 0}});
    const auto r = c.reversed();
    for (double t : {0.0, 0.2, 0.5, 0.77, 1.0}) REQUIRE((c.eval(t) - r.eval(1.0 - t)).norm() < 1e-14);
}

TEST_CASE("invalid B-spline input throws", "[bspline]")
{
    REQUIRE_THROWS_AS(BSpline2(3, {0, 0, 0, 1, 1, 1}, {{0, 0}, {1, 0}, {2, 0}}), std::invalid_argument);
    REQUIRE_THROWS_AS(BSpline2(2, {0, 0, 0, 1, 0.5, 1}, {{0, 0}, {1, 0}, {2, 0}}), std::invalid_argument);
    REQUIRE_THROWS_AS(BSpline2(2, {0, 0, 0, 0, 0, 0}, {{0, 0}, {1, 0}, {2, 0}}), std::invalid_argument);
    REQUIRE_THROWS_AS(BSpline2::clamped(3, {{0, 0}, {1, NAN}, {2, 0}, {3, 0}}), std::invalid_argument);
}

TEST_CASE("not-a-knot cubic interpolation reproduces a cubic polynomial exactly", "[interp]")
{
    auto f = [](double x) { return Vec2(x, 1.0 - 2.0 * x + 0.5 * x * x * x); };
    std::vector<double> t = {0.0, 0.13, 0.4, 0.55, 0.9, 1.3, 2.0};
    std::vector<Vec2> pts;
    for (double x : t) pts.push_back(f(x));
    const auto c = interpolate_cubic<2>(t, pts);
    for (double x = 0.0; x <= 2.0; x += 0.05) REQUIRE((c.eval(x) - f(x)).norm() < 1e-12);
}

TEST_CASE("clamped and natural end conditions hold", "[interp]")
{
    std::vector<double> t = {0, 1, 2, 3, 4};
    std::vector<Vec2> pts = {{0, 0}, {1, 1}, {2, 0}, {3, 1}, {4, 0}};
    const auto nat = interpolate_cubic<2>(t, pts, EndCondition<2>::natural(), EndCondition<2>::natural());
    REQUIRE(nat.derivs(0.0, 2)[2].norm() < 1e-12);
    REQUIRE(nat.derivs(4.0, 2)[2].norm() < 1e-12);
    const auto cl = interpolate_cubic<2>(t, pts, EndCondition<2>::derivative({1, 2}), EndCondition<2>::derivative({1, -3}));
    REQUIRE((cl.derivs(0.0, 1)[1] - Vec2(1, 2)).norm() < 1e-12);
    REQUIRE((cl.derivs(4.0, 1)[1] - Vec2(1, -3)).norm() < 1e-12);
    for (std::size_t i = 0; i < t.size(); ++i) REQUIRE((cl.eval(t[i]) - pts[i]).norm() < 1e-12);
    REQUIRE_THROWS_AS(interpolate_cubic<2>({0, 1, 1}, {{0, 0}, {1, 1}, {2, 2}}), std::invalid_argument);
}

TEST_CASE("solve_bspline meets position and derivative constraints (two-span quintic)", "[interp]")
{
    const double um = 0.37;
    std::vector<double> U = {0, 0, 0, 0, 0, 0, um, 1, 1, 1, 1, 1, 1};
    std::vector<Constraint<2>> cons = {{0.0, 0, {0, -0.01}}, {0.0, 1, {-0.02, 0.03}}, {0.0, 2, {0.4, -0.1}},
                                       {um, 0, {-0.005, 0.0}}, {1.0, 0, {0.01, 0.012}}, {1.0, 1, {0.03, 0.01}},
                                       {1.0, 2, {-0.2, 0.3}}};
    const auto c = solve_bspline<2>(5, U, cons);
    for (const auto& k : cons) REQUIRE((c.derivs(k.t, k.order)[static_cast<std::size_t>(k.order)] - k.value).norm() < 1e-12);
    // a quintic with a simple knot is C4 there
    const double h = 1e-7;
    const auto l = c.derivs(um - h, 4), r = c.derivs(um + h, 4);
    for (int k = 0; k <= 4; ++k) REQUIRE((l[static_cast<std::size_t>(k)] - r[static_cast<std::size_t>(k)]).norm() < 1e-3 * (1 + l[static_cast<std::size_t>(k)].norm()));
}
