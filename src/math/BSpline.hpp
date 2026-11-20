/// @file math/BSpline.hpp
/// @brief Non-rational B-spline curves in D dimensions with exact derivatives.
///
/// Algorithms: L. Piegl & W. Tiller, "The NURBS Book", 2nd ed., Springer 1997
///   A2.1 FindSpan, A2.2 BasisFuns, A2.3 DersBasisFuns, A3.2 CurveDerivsAlg1.
///
/// A curve of degree p with n+1 control points P_0..P_n needs m+1 = n+p+2
/// knots U = {u_0 <= ... <= u_m}.  The parameter domain is [u_p, u_{n+1}].
#pragma once

#include "math/Vec.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Math {

/// Clamped knot vector with uniformly spaced interior knots on [a, b]:
/// (p+1) copies of a, n_ctrl-p-1 interior knots, (p+1) copies of b.
[[nodiscard]] inline std::vector<double> clamped_uniform_knots(int degree, int n_ctrl,
                                                               double a = 0.0, double b = 1.0)
{
    if (degree < 1) throw std::invalid_argument("clamped_uniform_knots: degree must be >= 1");
    if (n_ctrl <= degree)
        throw std::invalid_argument("clamped_uniform_knots: need n_ctrl > degree");
    const int n_int = n_ctrl - degree - 1;
    std::vector<double> U;
    U.reserve(static_cast<std::size_t>(n_ctrl + degree + 1));
    for (int i = 0; i <= degree; ++i) U.push_back(a);
    for (int i = 1; i <= n_int; ++i) U.push_back(a + (b - a) * i / (n_int + 1));
    for (int i = 0; i <= degree; ++i) U.push_back(b);
    return U;
}

/// Knot span index i with U[i] <= t < U[i+1] (A2.1); the last non-empty span at t = u_{n+1}.
[[nodiscard]] inline int find_span(int n, int p, double t, const std::vector<double>& U)
{
    const auto un1 = static_cast<std::size_t>(n + 1);
    if (t >= U[un1]) {
        int k = n;  // last non-degenerate span
        while (k > p && U[static_cast<std::size_t>(k)] >= U[un1]) --k;
        return k;
    }
    if (t <= U[static_cast<std::size_t>(p)]) {
        int k = p;
        while (k < n && U[static_cast<std::size_t>(k + 1)] <= U[static_cast<std::size_t>(p)]) ++k;
        return k;
    }
    int low = p, high = n + 1;
    int mid = (low + high) / 2;
    while (t < U[static_cast<std::size_t>(mid)] || t >= U[static_cast<std::size_t>(mid + 1)]) {
        if (t < U[static_cast<std::size_t>(mid)]) high = mid;
        else low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

/// Non-zero basis functions and their derivatives up to order nd at t (A2.3).
/// Returns ders[k][j] = d^k N_{span-p+j,p}(t) / dt^k, k = 0..nd, j = 0..p.
[[nodiscard]] inline std::vector<std::vector<double>>
ders_basis_funs(int span, double t, int p, int nd, const std::vector<double>& U)
{
    const auto P = static_cast<std::size_t>(p);
    std::vector<std::vector<double>> ndu(P + 1, std::vector<double>(P + 1, 0.0));
    std::vector<double> left(P + 1, 0.0), right(P + 1, 0.0);
    ndu[0][0] = 1.0;
    for (std::size_t j = 1; j <= P; ++j) {
        left[j] = t - U[static_cast<std::size_t>(span) + 1 - j];
        right[j] = U[static_cast<std::size_t>(span) + j] - t;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    const auto ND = static_cast<std::size_t>(nd);
    std::vector<std::vector<double>> ders(ND + 1, std::vector<double>(P + 1, 0.0));
    for (std::size_t j = 0; j <= P; ++j) ders[0][j] = ndu[j][P];
    std::vector<std::vector<double>> a(2, std::vector<double>(P + 1, 0.0));
    for (int r = 0; r <= p; ++r) {
        std::size_t s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= nd; ++k) {
            double d = 0.0;
            const int rk = r - k, pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[static_cast<std::size_t>(pk + 1)][static_cast<std::size_t>(rk)];
                d = a[s2][0] * ndu[static_cast<std::size_t>(rk)][static_cast<std::size_t>(pk)];
            }
            const int j1 = (rk >= -1) ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (int j = j1; j <= j2; ++j) {
                const auto uj = static_cast<std::size_t>(j);
                a[s2][uj] = (a[s1][uj] - a[s1][uj - 1]) /
                            ndu[static_cast<std::size_t>(pk + 1)][static_cast<std::size_t>(rk + j)];
                d += a[s2][uj] * ndu[static_cast<std::size_t>(rk + j)][static_cast<std::size_t>(pk)];
            }
            if (r <= pk) {
                a[s2][static_cast<std::size_t>(k)] =
                    -a[s1][static_cast<std::size_t>(k - 1)] /
                    ndu[static_cast<std::size_t>(pk + 1)][static_cast<std::size_t>(r)];
                d += a[s2][static_cast<std::size_t>(k)] *
                     ndu[static_cast<std::size_t>(r)][static_cast<std::size_t>(pk)];
            }
            ders[static_cast<std::size_t>(k)][static_cast<std::size_t>(r)] = d;
            std::swap(s1, s2);
        }
    }
    double fac = p;
    for (int k = 1; k <= nd; ++k) {
        for (std::size_t j = 0; j <= P; ++j) ders[static_cast<std::size_t>(k)][j] *= fac;
        fac *= (p - k);
    }
    return ders;
}

/// B-spline curve C(t) = sum_i N_{i,p}(t) P_i in D dimensions.
template <int D>
class BSpline
{
public:
    using Point = Vec<D>;

    BSpline() = default;

    BSpline(int degree, std::vector<double> knots, std::vector<Point> ctrl)
        : p_(degree), U_(std::move(knots)), P_(std::move(ctrl))
    {
        validate();
    }

    /// Clamped curve on [a, b] with uniformly spaced interior knots.
    [[nodiscard]] static BSpline clamped(int degree, std::vector<Point> ctrl, double a = 0.0, double b = 1.0)
    {
        auto U = clamped_uniform_knots(degree, static_cast<int>(ctrl.size()), a, b);
        return BSpline(degree, std::move(U), std::move(ctrl));
    }

    [[nodiscard]] int degree() const noexcept { return p_; }
    [[nodiscard]] const std::vector<double>& knots() const noexcept { return U_; }
    [[nodiscard]] const std::vector<Point>& ctrl() const noexcept { return P_; }
    [[nodiscard]] int n_ctrl() const noexcept { return static_cast<int>(P_.size()); }
    [[nodiscard]] bool empty() const noexcept { return P_.empty(); }
    [[nodiscard]] double t0() const { return U_[static_cast<std::size_t>(p_)]; }
    [[nodiscard]] double t1() const { return U_[P_.size()]; }

    /// Point at parameter t (clamped into the domain).
    [[nodiscard]] Point eval(double t) const { return derivs(t, 0)[0]; }

    /// C(t), C'(t), ..., C^(nd)(t) (A3.2).  Derivatives above the degree are zero.
    [[nodiscard]] std::vector<Point> derivs(double t, int nd) const
    {
        if (P_.empty()) throw std::logic_error("BSpline: empty curve");
        if (nd < 0) throw std::invalid_argument("BSpline::derivs: nd must be >= 0");
        t = std::clamp(t, t0(), t1());
        const int n = n_ctrl() - 1;
        const int du = std::min(nd, p_);
        const int span = find_span(n, p_, t, U_);
        const auto nders = ders_basis_funs(span, t, p_, du, U_);
        std::vector<Point> out(static_cast<std::size_t>(nd + 1), Point::Zero());
        for (int k = 0; k <= du; ++k)
            for (int j = 0; j <= p_; ++j)
                out[static_cast<std::size_t>(k)] +=
                    nders[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)] *
                    P_[static_cast<std::size_t>(span - p_ + j)];
        return out;
    }

    /// Curve with the parameter direction reversed (same point set, domain mirrored).
    [[nodiscard]] BSpline reversed() const
    {
        std::vector<double> U(U_.size());
        const double a = U_.front(), b = U_.back();
        for (std::size_t i = 0; i < U_.size(); ++i) U[i] = a + b - U_[U_.size() - 1 - i];
        std::vector<Point> P(P_.rbegin(), P_.rend());
        return BSpline(p_, std::move(U), std::move(P));
    }

private:
    void validate() const
    {
        if (p_ < 1) throw std::invalid_argument("BSpline: degree must be >= 1");
        if (static_cast<int>(P_.size()) <= p_)
            throw std::invalid_argument("BSpline: need more control points than the degree");
        if (U_.size() != P_.size() + static_cast<std::size_t>(p_) + 1)
            throw std::invalid_argument("BSpline: knot count must be n_ctrl + degree + 1, got " +
                                        std::to_string(U_.size()));
        for (std::size_t i = 1; i < U_.size(); ++i)
            if (!(U_[i] >= U_[i - 1])) throw std::invalid_argument("BSpline: knots must be non-decreasing");
        if (!(t1() > t0())) throw std::invalid_argument("BSpline: empty parameter domain");
        for (const auto& q : P_)
            if (!q.allFinite()) throw std::invalid_argument("BSpline: non-finite control point");
    }

    int p_ = 0;
    std::vector<double> U_;
    std::vector<Point> P_;
};

using BSpline2 = BSpline<2>;
using BSpline3 = BSpline<3>;

}  // namespace BladeGen::Math
