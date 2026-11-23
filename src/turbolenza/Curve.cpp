// Port of v1 geometry.py ParametricBladeCurve / SurfaceCurves (see Curve.hpp).
#include "turbolenza/Curve.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace BladeGen::Turbolenza {

namespace {

std::size_t sz(int k) { return static_cast<std::size_t>(k); }

}  // namespace

ParametricBladeCurve::ParametricBladeCurve(std::vector<double> knots, std::vector<Point> control, double reference_length,
                                           int degree)
    : t_(std::move(knots)), c_(std::move(control)), L_(reference_length), p_(degree)
{
    const int k = p_;
    const int n = static_cast<int>(c_.size());
    if (k < 1 || n < k + 1 || static_cast<int>(t_.size()) != n + k + 1)
        throw std::invalid_argument("ParametricBladeCurve: needs (n, 2) control points and n + degree + 1 knots");
    for (int a = 0; a <= k; ++a)
        if (t_[sz(a)] != 0.0 || t_[t_.size() - 1 - sz(a)] != 1.0)
            throw std::invalid_argument("ParametricBladeCurve: the knot vector must be clamped on [0, 1]");
    for (std::size_t a = sz(k); a + 1 < t_.size() - sz(k); ++a)
        if (!(t_[a + 1] > t_[a]))
            throw std::invalid_argument("ParametricBladeCurve: the interior knots must be simple and increasing");
    if (!(std::isfinite(L_) && L_ > 0.0)) throw std::invalid_argument("ParametricBladeCurve: reference_length must be > 0");
}

std::vector<double> ParametricBladeCurve::interior_knots() const
{
    return {t_.begin() + p_ + 1, t_.end() - p_ - 1};
}

int ParametricBladeCurve::span(double u) const
{
    // t_l <= u < t_{l+1}, l in [p, n - 1]; u = 1 takes the last span (scipy BSpline convention)
    const int n = static_cast<int>(c_.size());
    if (u >= t_[sz(n)]) return n - 1;
    const auto it = std::upper_bound(t_.begin() + p_, t_.begin() + n + 1, u);
    return static_cast<int>(it - t_.begin()) - 1;
}

std::vector<Point> ParametricBladeCurve::derivatives(double u, int order) const
{
    if (order < 0 || order > p_)
        throw std::invalid_argument("ParametricBladeCurve: derivative order must lie in [0, " + std::to_string(p_) + "]");
    if (!std::isfinite(u) || u < -kCurveParameterSlack || u > 1.0 + kCurveParameterSlack)
        throw std::invalid_argument("ParametricBladeCurve: the curve parameter u must lie in [0, 1]");
    u = std::clamp(u, 0.0, 1.0);
    const int p = p_;
    const int l = span(u);
    // The NURBS Book A2.3: basis functions and their derivatives on span l
    std::vector<std::vector<double>> ndu(sz(p + 1), std::vector<double>(sz(p + 1), 0.0));
    std::vector<double> left(sz(p + 1)), right(sz(p + 1));
    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[sz(j)] = u - t_[sz(l + 1 - j)];
        right[sz(j)] = t_[sz(l + j)] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[sz(j)][sz(r)] = right[sz(r + 1)] + left[sz(j - r)];
            const double temp = ndu[sz(r)][sz(j - 1)] / ndu[sz(j)][sz(r)];
            ndu[sz(r)][sz(j)] = saved + right[sz(r + 1)] * temp;
            saved = left[sz(j - r)] * temp;
        }
        ndu[sz(j)][sz(j)] = saved;
    }
    std::vector<std::vector<double>> ders(sz(order + 1), std::vector<double>(sz(p + 1), 0.0));
    for (int j = 0; j <= p; ++j) ders[0][sz(j)] = ndu[sz(j)][sz(p)];
    std::vector<std::vector<double>> a(2, std::vector<double>(sz(p + 1), 0.0));
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= order; ++k) {
            double d = 0.0;
            const int rk = r - k, pk = p - k;
            if (r >= k) {
                a[sz(s2)][0] = a[sz(s1)][0] / ndu[sz(pk + 1)][sz(rk)];
                d = a[sz(s2)][0] * ndu[sz(rk)][sz(pk)];
            }
            const int j1 = rk >= -1 ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (int j = j1; j <= j2; ++j) {
                a[sz(s2)][sz(j)] = (a[sz(s1)][sz(j)] - a[sz(s1)][sz(j - 1)]) / ndu[sz(pk + 1)][sz(rk + j)];
                d += a[sz(s2)][sz(j)] * ndu[sz(rk + j)][sz(pk)];
            }
            if (r <= pk) {
                a[sz(s2)][sz(k)] = -a[sz(s1)][sz(k - 1)] / ndu[sz(pk + 1)][sz(r)];
                d += a[sz(s2)][sz(k)] * ndu[sz(r)][sz(pk)];
            }
            ders[sz(k)][sz(r)] = d;
            std::swap(s1, s2);
        }
    }
    double factor = static_cast<double>(p);
    for (int k = 1; k <= order; ++k) {
        for (int j = 0; j <= p; ++j) ders[sz(k)][sz(j)] *= factor;
        factor *= static_cast<double>(p - k);
    }
    std::vector<Point> out(sz(order + 1));
    for (int k = 0; k <= order; ++k) {
        double m = 0.0, y = 0.0;
        for (int j = 0; j <= p; ++j) {
            const Point& P = c_[sz(l - p + j)];
            m += ders[sz(k)][sz(j)] * P.m;
            y += ders[sz(k)][sz(j)] * P.y;
        }
        out[sz(k)] = {m, y};
    }
    // a clamped end is the end control point exactly (v1: de Boor is exact there only to round-off)
    if (u == 0.0) out[0] = c_.front();
    if (u == 1.0) out[0] = c_.back();
    return out;
}

Point ParametricBladeCurve::derivative(double u, int order) const { return derivatives(u, order)[sz(order)]; }

Point ParametricBladeCurve::unit_tangent(double u) const
{
    const Point d = derivative(u, 1);
    const double n = std::hypot(d.m, d.y);
    return {d.m / n, d.y / n};
}

double ParametricBladeCurve::curvature(double u) const
{
    const auto d = derivatives(u, 2);
    return signed_curvature(d[1], d[2]);
}

std::array<double, 2> ParametricBladeCurve::project(const Point& q, double u_guess) const
{
    double u = u_guess;
    if (u_guess < 0.0) {
        double best = std::numeric_limits<double>::infinity();
        for (int k = 0; k < 4097; ++k) {
            const double s = static_cast<double>(k) / 4096.0;
            const Point x = evaluate(s);
            const double d = std::hypot(q.m - x.m, q.y - x.y);
            if (d < best) {   // first minimum, as numpy argmin
                best = d;
                u = s;
            }
        }
    }
    for (int it = 0; it < kProjectionIterations; ++it) {
        const auto d = derivatives(u, 2);
        const double gm = d[0].m - q.m, gy = d[0].y - q.y;
        const double slope = gm * d[1].m + gy * d[1].y;
        const double curv = d[1].m * d[1].m + d[1].y * d[1].y + gm * d[2].m + gy * d[2].y;
        u = std::clamp(u - slope / curv, 0.0, 1.0);
    }
    const Point tan = unit_tangent(u);
    const Point x = evaluate(u);
    const double gm = q.m - x.m, gy = q.y - x.y;
    return {u, tan.m * gy - tan.y * gm};
}

ParametricBladeCurve ParametricBladeCurve::scaled(double factor) const
{
    std::vector<Point> c = c_;
    for (auto& P : c) P = {P.m * factor, P.y * factor};
    return {t_, std::move(c), L_ * factor, p_};
}

ParametricBladeCurve ParametricBladeCurve::rotated(double angle) const
{
    std::vector<Point> c = c_;
    for (auto& P : c) P = rotate_by_stagger(P, angle);
    return {t_, std::move(c), L_, p_};
}

SurfaceCurves SurfaceCurves::scaled(double factor) const
{
    return {suction.scaled(factor),      pressure.scaled(factor), suction_unstaggered.scaled(factor),
            pressure_unstaggered.scaled(factor), suction_parameters, pressure_parameters};
}

double signed_curvature(const Point& d1, const Point& d2)
{
    const double h = std::hypot(d1.m, d1.y);
    return (d1.m * d2.y - d1.y * d2.m) / (h * h * h);
}

double slope_min_cos() noexcept { return std::sqrt(std::numeric_limits<double>::epsilon()); }

std::array<double, 2> wall_slopes(const Point& d1, const Point& d2, double min_cos)
{
    if (!(min_cos > 0.0 && min_cos < 1.0)) throw std::invalid_argument("wall_slopes: min_cos must lie in (0, 1)");
    const double m1 = d1.m, y1 = d1.y, m2 = d2.m, y2 = d2.y;
    if (!(std::abs(m1) >= min_cos * std::hypot(m1, y1))) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan, nan};
    }
    return {y1 / m1, (y2 * m1 - y1 * m2) / (m1 * m1 * m1)};
}

Point rotate_by_stagger(const Point& p, double stagger)
{
    // v1: points @ [[c, s], [-s, c]]  ->  X = x c - y s, Y = x s + y c
    const double c = std::cos(stagger), s = std::sin(stagger);
    return {p.m * c - p.y * s, p.m * s + p.y * c};
}

}  // namespace BladeGen::Turbolenza
