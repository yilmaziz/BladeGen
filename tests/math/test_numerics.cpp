// Tests of math/Numerics.hpp and math/LeastSquares.hpp.
#include "math/LeastSquares.hpp"
#include "math/Numerics.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen::Math;
using Catch::Approx;

TEST_CASE("10-point Gauss-Legendre is exact to degree 19", "[quad]")
{
    for (int deg = 0; deg <= 19; ++deg) {
        const double I = gauss_legendre10([deg](double x) { return std::pow(x, deg); }, -0.3, 1.7);
        const double ref = (std::pow(1.7, deg + 1) - std::pow(-0.3, deg + 1)) / (deg + 1);
        REQUIRE(I == Approx(ref).epsilon(1e-13));
    }
}

TEST_CASE("adaptive integration of non-polynomial integrands", "[quad]")
{
    REQUIRE(integrate([](double x) { return std::sin(x); }, 0.0, kPi) == Approx(2.0).epsilon(1e-13));
    REQUIRE(integrate([](double x) { return std::sqrt(x); }, 0.0, 1.0, 1e-12) == Approx(2.0 / 3.0).epsilon(1e-10));
}

TEST_CASE("Brent finds bracketed roots and rejects unbracketed ones", "[roots]")
{
    REQUIRE(brent_root([](double x) { return x * x - 2.0; }, 0.0, 2.0) == Approx(std::sqrt(2.0)).epsilon(1e-15));
    REQUIRE(brent_root([](double x) { return std::cos(x) - x; }, 0.0, 1.0) == Approx(0.7390851332151607).epsilon(1e-15));
    REQUIRE_THROWS_AS(brent_root([](double x) { return x * x + 1.0; }, -1.0, 1.0), std::invalid_argument);
    REQUIRE(std::isnan(first_root([](double x) { return x * x + 1.0; }, -1.0, 1.0)));
    REQUIRE(first_root([](double x) { return std::sin(x); }, 0.5, 10.0) == Approx(kPi).epsilon(1e-14));
}

TEST_CASE("arc length of a B-spline line and of a parabola", "[quad]")
{
    // a straight line with a non-uniform parametrisation
    const auto line = BSpline2::clamped(3, {{0, 0}, {0.1, 0.2}, {1.5, 3.0}, {2.0, 4.0}});
    REQUIRE(arc_length(line, 0.0, 1.0) == Approx(std::sqrt(20.0)).epsilon(1e-13));
    // y = x^2 on [0, 1] as a quadratic Bezier: P0 (0,0), P1 (0.5, 0), P2 (1, 1)
    const auto par = BSpline2::clamped(2, {{0, 0}, {0.5, 0}, {1, 1}});
    const double ref = 0.5 * std::sqrt(5.0) + 0.25 * std::asinh(2.0);
    REQUIRE(arc_length(par, 0.0, 1.0) == Approx(ref).epsilon(1e-13));
}

TEST_CASE("Levenberg-Marquardt solves Rosenbrock and respects bounds", "[lsq]")
{
    auto rosen = [](const Eigen::VectorXd& x) {
        Eigen::VectorXd r(2);
        r << 10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0];
        return r;
    };
    Eigen::VectorXd x0(2), lb(2), ub(2);
    x0 << -1.2, 1.0;
    lb << -5, -5;
    ub << 5, 5;
    auto res = levenberg_marquardt(rosen, x0, lb, ub);
    REQUIRE(res.x[0] == Approx(1.0).margin(1e-8));
    REQUIRE(res.x[1] == Approx(1.0).margin(1e-8));
    // with x0 <= 0.5 the constrained minimum is on the bound x0 = 0.5, x1 = 0.25
    ub << 0.5, 5;
    res = levenberg_marquardt(rosen, x0, lb, ub);
    REQUIRE(res.x[0] == Approx(0.5).margin(1e-9));
    REQUIRE(res.x[1] == Approx(0.25).margin(1e-6));
    REQUIRE_THROWS_AS(levenberg_marquardt(rosen, x0, ub, lb), std::invalid_argument);
}

TEST_CASE("Levenberg-Marquardt fits an exponential decay", "[lsq]")
{
    std::vector<double> t, y;
    for (int i = 0; i < 30; ++i) {
        t.push_back(0.1 * i);
        y.push_back(2.5 * std::exp(-1.3 * t.back()) + 0.4);
    }
    auto r = [&](const Eigen::VectorXd& p) {
        Eigen::VectorXd out(static_cast<Eigen::Index>(t.size()));
        for (std::size_t i = 0; i < t.size(); ++i)
            out[static_cast<Eigen::Index>(i)] = p[0] * std::exp(-p[1] * t[i]) + p[2] - y[i];
        return out;
    };
    Eigen::VectorXd x0(3), lb(3), ub(3);
    x0 << 1, 0.5, 0;
    lb << 0, 0, -1;
    ub << 10, 5, 1;
    const auto res = levenberg_marquardt(r, x0, lb, ub);
    REQUIRE(res.x[0] == Approx(2.5).epsilon(1e-8));
    REQUIRE(res.x[1] == Approx(1.3).epsilon(1e-8));
    REQUIRE(res.x[2] == Approx(0.4).epsilon(1e-8));
}
