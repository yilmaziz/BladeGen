/// @file math/Numerics.hpp
/// @brief Quadrature (Gauss–Legendre), scalar root finding (Brent) and the
///        arc length of B-spline curves.
#pragma once

#include "math/BSpline.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace BladeGen::Math {

/// 10-point Gauss–Legendre rule on [a, b] (exact for polynomials of degree <= 19).
template <class F>
[[nodiscard]] double gauss_legendre10(F&& f, double a, double b)
{
    static constexpr std::array<double, 5> x = {0.1488743389816312, 0.4333953941292472,
                                                0.6794095682990244, 0.8650633666889845,
                                                0.9739065285171717};
    static constexpr std::array<double, 5> w = {0.2955242247147529, 0.2692667193099963,
                                                0.2190863625159820, 0.1494513491505806,
                                                0.0666713443086881};
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    double s = 0.0;
    for (std::size_t i = 0; i < 5; ++i) s += w[i] * (f(c - h * x[i]) + f(c + h * x[i]));
    return s * h;
}

/// Adaptive Gauss–Legendre: bisects until two halves agree with the whole to `tol`.
template <class F>
[[nodiscard]] double integrate(F&& f, double a, double b, double tol = 1e-13, int depth = 40)
{
    const double whole = gauss_legendre10(f, a, b);
    const double m = 0.5 * (a + b);
    const double left = gauss_legendre10(f, a, m), right = gauss_legendre10(f, m, b);
    if (depth <= 0 || std::abs(left + right - whole) <= tol * std::max(1.0, std::abs(whole)))
        return left + right;
    return integrate(f, a, m, 0.5 * tol, depth - 1) + integrate(f, m, b, 0.5 * tol, depth - 1);
}

/// Brent's method (R. P. Brent, "Algorithms for Minimization without Derivatives",
/// 1973, ch. 4) for a root of f in [a, b]; f(a), f(b) must have opposite signs.
template <class F>
[[nodiscard]] double brent_root(F&& f, double a, double b, double xtol = 1e-14, int max_iter = 200)
{
    double fa = f(a), fb = f(b);
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if ((fa > 0.0) == (fb > 0.0)) throw std::invalid_argument("brent_root: root not bracketed");
    double c = a, fc = fa, d = b - a, e = d;
    for (int it = 0; it < max_iter; ++it) {
        if ((fb > 0.0) == (fc > 0.0)) { c = a; fc = fa; d = e = b - a; }
        if (std::abs(fc) < std::abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
        const double tol1 = 2.0 * std::numeric_limits<double>::epsilon() * std::abs(b) + 0.5 * xtol;
        const double xm = 0.5 * (c - b);
        if (std::abs(xm) <= tol1 || fb == 0.0) return b;
        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
            double p, q, r;
            const double s = fb / fa;
            if (a == c) { p = 2.0 * xm * s; q = 1.0 - s; }
            else {
                q = fa / fc; r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = std::abs(p);
            if (2.0 * p < std::min(3.0 * xm * q - std::abs(tol1 * q), std::abs(e * q))) { e = d; d = p / q; }
            else { d = xm; e = d; }
        } else { d = xm; e = d; }
        a = b; fa = fb;
        b += (std::abs(d) > tol1) ? d : (xm > 0.0 ? tol1 : -tol1);
        fb = f(b);
    }
    return b;
}

/// First root of f on [a, b]: scans n sub-intervals for a sign change, then Brent.
/// Returns NaN when no sign change is found.
template <class F>
[[nodiscard]] double first_root(F&& f, double a, double b, int n = 64, double xtol = 1e-14)
{
    double x0 = a, f0 = f(a);
    if (f0 == 0.0) return a;
    for (int i = 1; i <= n; ++i) {
        const double x1 = a + (b - a) * i / n, f1 = f(x1);
        if (f1 == 0.0) return x1;
        if ((f0 > 0.0) != (f1 > 0.0)) return brent_root(f, x0, x1, xtol);
        x0 = x1; f0 = f1;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

/// Arc length of a B-spline between parameters ta < tb (per knot span, adaptive GL).
template <int D>
[[nodiscard]] double arc_length(const BSpline<D>& c, double ta, double tb, double tol = 1e-13)
{
    if (tb < ta) return -arc_length(c, tb, ta, tol);
    auto speed = [&](double t) { return c.derivs(t, 1)[1].norm(); };
    double s = 0.0, lo = ta;
    for (double k : c.knots()) {
        if (k > ta && k < tb) { s += integrate(speed, lo, k, tol); lo = k; }
    }
    return s + integrate(speed, lo, tb, tol);
}

/// Cumulative arc-length table s(t) of a B-spline at the parameters t (ascending).
template <int D>
[[nodiscard]] std::vector<double> arc_length_table(const BSpline<D>& c, const std::vector<double>& t)
{
    std::vector<double> s(t.size(), 0.0);
    for (std::size_t i = 1; i < t.size(); ++i) s[i] = s[i - 1] + arc_length(c, t[i - 1], t[i]);
    return s;
}

/// n parameters of a clustered distribution on [0, 1]: x = (1 - cos(pi*u))/2 for both
/// ends clustered, or uniform.
[[nodiscard]] inline std::vector<double> cosine_spacing(int n, bool clustered = true)
{
    if (n < 2) throw std::invalid_argument("cosine_spacing: need n >= 2");
    std::vector<double> u(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double s = static_cast<double>(i) / (n - 1);
        u[static_cast<std::size_t>(i)] = clustered ? 0.5 * (1.0 - std::cos(kPi * s)) : s;
    }
    u.front() = 0.0;
    u.back() = 1.0;
    return u;
}

}  // namespace BladeGen::Math
