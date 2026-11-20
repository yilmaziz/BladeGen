/// @file math/FairCurve.hpp
/// @brief Planar B-spline of minimum third-derivative energy under linear constraints.
///
/// Minimises E(P) = \int_{t0}^{t1} |C'''(t)|^2 dt over the control points P of a
/// planar B-spline with fixed degree and knots, subject to linear equality
/// constraints a_x . P_x + a_y . P_y = b.  E is a quadratic form in P, so the
/// minimiser solves one symmetric KKT system (Lagrange multipliers):
///
///     [ 2 Q   A^T ] [P]   [0]
///     [  A     0  ] [l] = [b]
///
/// The third-derivative ("jerk") energy is a standard fairing functional:
/// minimal-energy splines of this kind give slowly varying curvature (see e.g.
/// H. P. Moreton & C. H. Séquin, "Functional optimization for fair surface
/// design", SIGGRAPH 1992, and G. Farin, "Curves and Surfaces for CAGD", §23).
#pragma once

#include "math/BSpline.hpp"

#include <Eigen/Dense>

#include <stdexcept>
#include <vector>

namespace BladeGen::Math {

/// One scalar linear condition on a planar B-spline: w . d^order C/dt^order (t) = value.
struct LinearCondition
{
    double t;
    int order;
    Vec2 w;      ///< direction the derivative is projected on
    double value;
};

/// Gram matrix Q_ij = \int N_i^(r) N_j^(r) dt (10-point Gauss on every knot span).
[[nodiscard]] inline Eigen::MatrixXd derivative_gram(int degree, const std::vector<double>& U, int r)
{
    static constexpr double x[5] = {0.1488743389816312, 0.4333953941292472, 0.6794095682990244,
                                    0.8650633666889845, 0.9739065285171717};
    static constexpr double w[5] = {0.2955242247147529, 0.2692667193099963, 0.2190863625159820,
                                    0.1494513491505806, 0.0666713443086881};
    const int n_ctrl = static_cast<int>(U.size()) - degree - 1;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(n_ctrl, n_ctrl);
    const int n = n_ctrl - 1;
    for (int s = degree; s <= n; ++s) {
        const double a = U[static_cast<std::size_t>(s)], b = U[static_cast<std::size_t>(s + 1)];
        if (!(b > a)) continue;
        const double c = 0.5 * (a + b), h = 0.5 * (b - a);
        for (int q = 0; q < 10; ++q) {
            const double xi = (q < 5) ? -x[q] : x[q - 5];
            const double wi = w[q < 5 ? q : q - 5] * h;
            const auto N = ders_basis_funs(s, c + h * xi, degree, r, U);
            for (int i = 0; i <= degree; ++i)
                for (int j = 0; j <= degree; ++j)
                    Q(s - degree + i, s - degree + j) +=
                        wi * N[static_cast<std::size_t>(r)][static_cast<std::size_t>(i)] *
                        N[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)];
        }
    }
    return Q;
}

/// Minimum-energy planar B-spline (energy order `r`, default 3) meeting all conditions.
[[nodiscard]] inline BSpline2 fair_bspline(int degree, std::vector<double> U,
                                           const std::vector<LinearCondition>& cons, int r = 3)
{
    const int n_ctrl = static_cast<int>(U.size()) - degree - 1;
    if (n_ctrl <= degree) throw std::invalid_argument("fair_bspline: too few knots");
    const int m = static_cast<int>(cons.size());
    const int nv = 2 * n_ctrl;
    if (m > nv) throw std::invalid_argument("fair_bspline: more conditions than unknowns");
    const Eigen::MatrixXd Q = derivative_gram(degree, U, r);
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nv + m, nv + m);
    K.block(0, 0, n_ctrl, n_ctrl) = 2.0 * Q;
    K.block(n_ctrl, n_ctrl, n_ctrl, n_ctrl) = 2.0 * Q;
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nv + m);
    const int n = n_ctrl - 1;
    for (int k = 0; k < m; ++k) {
        const auto& c = cons[static_cast<std::size_t>(k)];
        if (c.order < 0 || c.order > degree) throw std::invalid_argument("fair_bspline: bad derivative order");
        const double tt = std::clamp(c.t, U[static_cast<std::size_t>(degree)], U[static_cast<std::size_t>(n + 1)]);
        const int span = find_span(n, degree, tt, U);
        const auto N = ders_basis_funs(span, tt, degree, c.order, U);
        for (int j = 0; j <= degree; ++j) {
            const double v = N[static_cast<std::size_t>(c.order)][static_cast<std::size_t>(j)];
            const int col = span - degree + j;
            K(nv + k, col) = c.w.x() * v;
            K(nv + k, n_ctrl + col) = c.w.y() * v;
            K(col, nv + k) = c.w.x() * v;
            K(n_ctrl + col, nv + k) = c.w.y() * v;
        }
        rhs[nv + k] = c.value;
    }
    const Eigen::VectorXd X = K.fullPivLu().solve(rhs);
    if (!X.allFinite() || (K * X - rhs).norm() > 1e-8 * (1.0 + rhs.norm()))
        throw std::runtime_error("fair_bspline: inconsistent or singular conditions");
    std::vector<Vec2> P(static_cast<std::size_t>(n_ctrl));
    for (int i = 0; i < n_ctrl; ++i) P[static_cast<std::size_t>(i)] = Vec2(X[i], X[n_ctrl + i]);
    return BSpline2(degree, std::move(U), std::move(P));
}

}  // namespace BladeGen::Math
