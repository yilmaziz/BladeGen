// ThroughflowNumerics.hpp: Aungier (2003) Eqs. 7-67 to 7-73 (three-point differences and integrals on
// non-uniform points) and Eqs. 7-34 to 7-36 (the linear ODE and its solution), against closed forms
// and a direct numerical integration.
#include "turbolenza/ThroughflowNumerics.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace TF = BladeGen::Turbolenza::Throughflow;

namespace {

/// Non-uniform points on [a, b]: x_k = a + (b - a) (k / (n - 1))^1.6 (strongly graded).
std::vector<double> graded(double a, double b, int n)
{
    std::vector<double> x(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) x[static_cast<std::size_t>(k)] = a + (b - a) * std::pow(k / (n - 1.0), 1.6);
    return x;
}

std::vector<double> map(const std::vector<double>& x, double (*f)(double))
{
    std::vector<double> y(x.size());
    for (std::size_t k = 0; k < x.size(); ++k) y[k] = f(x[k]);
    return y;
}

double quad(double x) { return 3.0 - 2.0 * x + 5.0 * x * x; }
double dquad(double x) { return -2.0 + 10.0 * x; }
double quad_integral(double a, double b)
{
    auto F = [](double x) { return 3.0 * x - x * x + 5.0 * x * x * x / 3.0; };
    return F(b) - F(a);
}

}  // namespace

TEST_CASE("throughflow numerics: Eqs. 7-67 to 7-71 are exact for quadratics on non-uniform points",
          "[turbolenza][throughflow][numerics]")
{
    const auto x = graded(0.2, 1.7, 9);
    const auto f = map(x, quad);
    const auto d3 = TF::derivative(x, f, TF::EndRule::ThreePoint);
    const auto d2 = TF::second_derivative(x, f);
    double e1 = 0.0, e2 = 0.0;
    for (std::size_t k = 0; k < x.size(); ++k) {
        e1 = std::max(e1, std::abs(d3[k] - dquad(x[k])));
        if (k > 0 && k + 1 < x.size()) e2 = std::max(e2, std::abs(d2[k] - 10.0));
    }
    std::printf("[numerics] quadratic: max |f' - exact| = %.3e (Eqs. 7-67..7-69), max |f'' - 10| = %.3e (Eq. 7-71)\n", e1, e2);
    CHECK(e1 <= 1e-12);  // measured 3.6e-15
    CHECK(e2 <= 1e-11);  // measured 3.6e-14
    CHECK(std::isnan(d2.front()));
    CHECK(std::isnan(d2.back()));
    // Eq. 7-70 (two-point ends) is exact for lines only: on the quadratic its end error is the
    // first-order Taylor remainder 5 (x2 - x1) (the interior points keep Eq. 7-67)
    const auto d2p = TF::derivative(x, f, TF::EndRule::TwoPoint);
    CHECK(std::abs(d2p.front() - (dquad(x[0]) + 5.0 * (x[1] - x[0]))) <= 1e-12);
    CHECK(std::abs(d2p.back() - (dquad(x[8]) - 5.0 * (x[8] - x[7]))) <= 1e-12);
    CHECK(d2p[4] == d3[4]);
    const auto lin = map(x, [](double v) { return 4.0 - 7.0 * v; });
    for (double v : TF::derivative(x, lin, TF::EndRule::TwoPoint)) CHECK(std::abs(v + 7.0) <= 1e-12);
}

TEST_CASE("throughflow numerics: Eqs. 7-72 / 7-73 integrate quadratics exactly, segment by segment",
          "[turbolenza][throughflow][numerics]")
{
    for (int n : {3, 4, 9}) {
        const auto x = graded(-0.5, 2.0, n);
        const auto f = map(x, quad);
        const auto I = TF::cumulative_integral(x, f);
        double e = 0.0;
        for (std::size_t k = 0; k < x.size(); ++k) e = std::max(e, std::abs(I[k] - quad_integral(x[0], x[k])));
        CHECK(e <= 1e-12);  // measured <= 1.8e-15
        for (int k = 0; k + 1 < n; ++k) {
            const auto ku = static_cast<std::size_t>(k);
            CHECK(std::abs(TF::segment_integral(x, f, k) - quad_integral(x[ku], x[ku + 1])) <= 1e-13);
            // the partial integral is the same series: exact for quadratics anywhere in the segment
            const double xi = x[ku] + 0.37 * (x[ku + 1] - x[ku]);
            CHECK(std::abs(TF::partial_segment_integral(x, f, k, xi) - quad_integral(x[ku], xi)) <= 1e-13);
            const double full = TF::segment_integral(x, f, k);
            CHECK(std::abs(TF::partial_segment_integral(x, f, k, x[ku + 1]) - full) <= 1e-14 * std::abs(full));  // round-off
        }
    }
}

TEST_CASE("throughflow numerics: the three-point integral of a smooth function converges at third order or better",
          "[turbolenza][throughflow][numerics]")
{
    // int_0^1 exp(2 y) dy = (e^2 - 1) / 2
    const double exact = 0.5 * (std::exp(2.0) - 1.0);
    double prev = 0.0;
    for (int n : {6, 11, 21, 41}) {
        const auto x = graded(0.0, 1.0, n);
        const double e = std::abs(TF::integral(x, map(x, [](double v) { return std::exp(2.0 * v); })) - exact);
        std::printf("[numerics] int exp(2y), n = %2d: error %.3e%s\n", n, e, prev > 0.0 ? "" : " (first)");
        if (prev > 0.0) CHECK(prev / e >= 3.0);  // measured ratios 3.6, 3.6, 5.9 on these graded points
        prev = e;
    }
    prev = 0.0;
    for (int n : {6, 11, 21, 41}) {  // equal spacing
        std::vector<double> x(static_cast<std::size_t>(n));
        for (int k = 0; k < n; ++k) x[static_cast<std::size_t>(k)] = k / (n - 1.0);
        const double e = std::abs(TF::integral(x, map(x, [](double v) { return std::exp(2.0 * v); })) - exact);
        std::printf("[numerics] int exp(2y), equal spacing, n = %2d: error %.3e\n", n, e);
        if (prev > 0.0) CHECK(prev / e >= 4.5);  // measured ratios 4.7, 6.4, 7.2 (toward 8: third order)
        prev = e;
    }
}

TEST_CASE("throughflow numerics: Eqs. 7-35 / 7-36 against a direct numerical integration of Eq. 7-34",
          "[turbolenza][throughflow][numerics]")
{
    // dW/dy = f1 W + f4 with f1 = 3 sin(2 y) - 1, f4 = 40 y^2 - 20 y + 7, W(0) = 50 on [0, 0.8]
    auto f1 = [](double y) { return 3.0 * std::sin(2.0 * y) - 1.0; };
    auto f4 = [](double y) { return 40.0 * y * y - 20.0 * y + 7.0; };
    const double w0 = 50.0, ys = 0.8;
    // reference: classical RK4 with 200000 steps
    const int nrk = 200000;
    const double h = ys / nrk;
    std::vector<double> yr(nrk + 1), wr(nrk + 1);
    double w = w0;
    for (int k = 0; k <= nrk; ++k) {
        const double y = k * h;
        yr[static_cast<std::size_t>(k)] = y;
        wr[static_cast<std::size_t>(k)] = w;
        if (k == nrk) break;
        auto rhs = [&](double yy, double ww) { return f1(yy) * ww + f4(yy); };
        const double k1 = rhs(y, w), k2 = rhs(y + h / 2, w + h / 2 * k1), k3 = rhs(y + h / 2, w + h / 2 * k2),
                     k4 = rhs(y + h, w + h * k3);
        w += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    double prev = 0.0;
    for (int n : {6, 11, 21, 41, 81}) {
        const auto y = graded(0.0, ys, n);
        std::vector<double> a(y.size()), b(y.size());
        for (std::size_t k = 0; k < y.size(); ++k) {
            a[k] = f1(y[k]);
            b[k] = f4(y[k]);
        }
        const auto s = TF::solve_linear_ode(y, a, b);
        const auto W = s.evaluate(w0);
        CHECK(W.front() == w0);
        double e = 0.0;
        for (std::size_t k = 0; k < y.size(); ++k) {
            const auto idx = static_cast<std::size_t>(std::lround(y[k] / h));
            if (std::abs(static_cast<double>(idx) * h - y[k]) < 1e-12) e = std::max(e, std::abs(W[k] - wr[idx]) / wr[idx]);
            else {
                // interpolate the RK4 reference linearly (h = 4e-6: error ~1e-11 relative)
                const auto i0 = static_cast<std::size_t>(y[k] / h);
                const double t = (y[k] - yr[i0]) / h;
                const double ref = wr[i0] + t * (wr[i0 + 1] - wr[i0]);
                e = std::max(e, std::abs(W[k] - ref) / ref);
            }
        }
        std::printf("[numerics] ODE 7-34 by 7-35/7-36, n = %2d: max relative error %.3e\n", n, e);
        if (prev > 0.0) CHECK(prev / e >= 6.0);  // measured ratios 9.0 .. 15.6
        prev = e;
    }
    CHECK(prev <= 1e-6);  // n = 81: measured 2.1e-7
}

TEST_CASE("throughflow numerics: refusals", "[turbolenza][throughflow][numerics]")
{
    const std::vector<double> x = {0.0, 1.0, 2.0}, f = {1.0, 2.0, 3.0};
    CHECK_THROWS_AS(TF::derivative({0.0, 1.0}, {1.0, 2.0}), std::invalid_argument);
    CHECK_THROWS_AS(TF::derivative({0.0, 1.0, 1.0}, f), std::invalid_argument);
    CHECK_THROWS_AS(TF::derivative({0.0, 2.0, 1.0}, f), std::invalid_argument);
    CHECK_THROWS_AS(TF::derivative(x, {1.0, 2.0}), std::invalid_argument);
    CHECK_THROWS_AS(TF::integral(x, {1.0, std::numeric_limits<double>::quiet_NaN(), 3.0}), std::invalid_argument);
    CHECK_THROWS_AS(TF::segment_integral(x, f, 2), std::invalid_argument);
    CHECK_THROWS_AS(TF::segment_integral(x, f, -1), std::invalid_argument);
    CHECK_THROWS_AS(TF::partial_segment_integral(x, f, 0, 1.5), std::invalid_argument);
    CHECK_THROWS_AS(TF::solve_linear_ode(x, f, {1.0, 1.0}), std::invalid_argument);
}
