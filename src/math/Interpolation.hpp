/// @file math/Interpolation.hpp
/// @brief B-spline interpolation from point / derivative constraints.
///
/// Everything here reduces to one square linear system: for a fixed degree and
/// knot vector, every constraint "d^k C / dt^k (t_i) = v_i" is one row of basis
/// function derivatives (The NURBS Book §9.2).  The system is banded and solved
/// with Eigen's sparse LU.
#pragma once

#include "math/BSpline.hpp"

#include <Eigen/Sparse>

#include <stdexcept>
#include <vector>

namespace BladeGen::Math {

/// One interpolation condition: the order-th derivative at parameter t.
template <int D>
struct Constraint
{
    double t;
    int order;  ///< 0 = position, 1 = first derivative, ...
    Vec<D> value;
};

/// Control points of the degree-p B-spline on `knots` that meets every constraint.
/// The number of constraints must equal the number of control points.
template <int D>
[[nodiscard]] BSpline<D> solve_bspline(int degree, std::vector<double> knots,
                                       const std::vector<Constraint<D>>& cons)
{
    const auto n_ctrl = static_cast<int>(knots.size()) - degree - 1;
    if (n_ctrl <= degree) throw std::invalid_argument("solve_bspline: too few knots");
    if (static_cast<int>(cons.size()) != n_ctrl)
        throw std::invalid_argument("solve_bspline: need exactly one constraint per control point");
    const int n = n_ctrl - 1;
    std::vector<Eigen::Triplet<double>> trip;
    Eigen::Matrix<double, Eigen::Dynamic, D> rhs(n_ctrl, D);
    for (int r = 0; r < n_ctrl; ++r) {
        const auto& c = cons[static_cast<std::size_t>(r)];
        if (c.order < 0 || c.order > degree)
            throw std::invalid_argument("solve_bspline: derivative order out of range");
        const double tt = std::clamp(c.t, knots[static_cast<std::size_t>(degree)],
                                     knots[static_cast<std::size_t>(n + 1)]);
        const int span = find_span(n, degree, tt, knots);
        const auto N = ders_basis_funs(span, tt, degree, c.order, knots);
        for (int j = 0; j <= degree; ++j) {
            const double v = N[static_cast<std::size_t>(c.order)][static_cast<std::size_t>(j)];
            if (v != 0.0) trip.emplace_back(r, span - degree + j, v);
        }
        rhs.row(r) = c.value.transpose();
    }
    Eigen::SparseMatrix<double> A(n_ctrl, n_ctrl);
    A.setFromTriplets(trip.begin(), trip.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(A);
    if (lu.info() != Eigen::Success)
        throw std::runtime_error("solve_bspline: singular interpolation system");
    const Eigen::Matrix<double, Eigen::Dynamic, D> X = lu.solve(rhs);
    if (lu.info() != Eigen::Success || !X.allFinite())
        throw std::runtime_error("solve_bspline: interpolation solve failed");
    std::vector<Vec<D>> P(static_cast<std::size_t>(n_ctrl));
    for (int i = 0; i < n_ctrl; ++i) P[static_cast<std::size_t>(i)] = X.row(i).transpose();
    return BSpline<D>(degree, std::move(knots), std::move(P));
}

/// End condition of a cubic interpolating spline.
template <int D>
struct EndCondition
{
    enum class Kind { NotAKnot, Natural, FirstDerivative } kind = Kind::NotAKnot;
    Vec<D> value = Vec<D>::Zero();  ///< dC/dt at the end for FirstDerivative

    [[nodiscard]] static EndCondition not_a_knot() { return {}; }
    [[nodiscard]] static EndCondition natural() { return {Kind::Natural, Vec<D>::Zero()}; }
    [[nodiscard]] static EndCondition derivative(const Vec<D>& d) { return {Kind::FirstDerivative, d}; }
};

/// C2 cubic spline through points[i] at strictly increasing params[i].
/// Not-a-knot ends drop the knots at params[1] / params[n-2] (needs >= 4 points;
/// with 2 or 3 points a lower-degree / natural fallback is used).
template <int D>
[[nodiscard]] BSpline<D> interpolate_cubic(const std::vector<double>& t, const std::vector<Vec<D>>& pts,
                                           EndCondition<D> left = EndCondition<D>::not_a_knot(),
                                           EndCondition<D> right = EndCondition<D>::not_a_knot())
{
    const std::size_t n = t.size();
    if (n != pts.size()) throw std::invalid_argument("interpolate_cubic: size mismatch");
    if (n < 2) throw std::invalid_argument("interpolate_cubic: need at least 2 points");
    for (std::size_t i = 1; i < n; ++i)
        if (!(t[i] > t[i - 1])) throw std::invalid_argument("interpolate_cubic: parameters must increase");
    using K = typename EndCondition<D>::Kind;
    if (n < 4) {  // too few points for not-a-knot: treat as natural
        if (left.kind == K::NotAKnot) left = EndCondition<D>::natural();
        if (right.kind == K::NotAKnot) right = EndCondition<D>::natural();
    }
    std::vector<double> U;
    for (int i = 0; i < 4; ++i) U.push_back(t.front());
    const std::size_t i0 = (left.kind == K::NotAKnot) ? 2 : 1;
    const std::size_t i1 = (right.kind == K::NotAKnot) ? n - 2 : n - 1;
    for (std::size_t i = i0; i < i1; ++i) U.push_back(t[i]);
    for (int i = 0; i < 4; ++i) U.push_back(t.back());

    std::vector<Constraint<D>> cons;
    auto end_row = [&](const EndCondition<D>& e, double tt) {
        if (e.kind == K::Natural) cons.push_back({tt, 2, Vec<D>::Zero()});
        else if (e.kind == K::FirstDerivative) cons.push_back({tt, 1, e.value});
    };
    cons.push_back({t.front(), 0, pts.front()});
    end_row(left, t.front());
    for (std::size_t i = 1; i + 1 < n; ++i) cons.push_back({t[i], 0, pts[i]});
    end_row(right, t.back());
    cons.push_back({t.back(), 0, pts.back()});
    return solve_bspline<D>(3, std::move(U), cons);
}

/// Chord-length parameters of a polyline, starting at 0.
template <int D>
[[nodiscard]] std::vector<double> chord_length_params(const std::vector<Vec<D>>& pts)
{
    std::vector<double> t(pts.size(), 0.0);
    for (std::size_t i = 1; i < pts.size(); ++i) t[i] = t[i - 1] + (pts[i] - pts[i - 1]).norm();
    return t;
}

}  // namespace BladeGen::Math
