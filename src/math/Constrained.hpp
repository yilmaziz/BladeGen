/// @file math/Constrained.hpp
/// @brief Constrained linear least squares and Gauss–Legendre rules of any order.
///
/// C. L. Lawson & R. J. Hanson, "Solving Least Squares Problems", Prentice-Hall 1974:
///   nnls  Ch. 23, Algorithm NNLS:  min |A x - b|  subject to x >= 0
///   ldp   Ch. 23:                  min |z|        subject to G z >= h   (through NNLS)
///   lsei  Chs. 20-23:              min |A x - b|  subject to E x = f, x >= 0
/// (the equalities are eliminated through a null-space basis of E, the rest is an LDP).
/// Used by the NACA 65 leading-edge fair (port of TURBOlenza's _lsei / _ldp).
#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace BladeGen::Math {

/// n-point Gauss–Legendre nodes and weights on [-1, 1] (Newton on P_n).
inline void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w)
{
    if (n < 1) throw std::invalid_argument("gauss_legendre: n >= 1");
    x.assign(static_cast<std::size_t>(n), 0.0);
    w.assign(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z = std::cos(3.14159265358979323846 * (i + 0.75) / (n + 0.5)), dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = z;
            for (int k = 2; k <= n; ++k) {
                const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
                p0 = p1;
                p1 = p2;
            }
            dp = n * (z * p1 - p0) / (z * z - 1.0);
            const double dz = p1 / dp;
            z -= dz;
            if (std::abs(dz) < 1e-16) break;
        }
        const auto a = static_cast<std::size_t>(i), b = static_cast<std::size_t>(n - 1 - i);
        x[a] = -z;
        x[b] = z;
        w[a] = w[b] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
}

/// Lawson–Hanson NNLS: min |A x - b| subject to x >= 0.  An index whose least-squares value
/// on entering the passive set is not positive is blocked until the next step that makes
/// progress (the L&H safeguard against cycling on degenerate problems).
[[nodiscard]] inline Eigen::VectorXd nnls(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, int max_iter = 0)
{
    const Eigen::Index n = A.cols();
    if (max_iter <= 0) max_iter = static_cast<int>(3 * n + 50);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    std::vector<bool> passive(static_cast<std::size_t>(n), false), blocked(static_cast<std::size_t>(n), false);
    const double tol = 10.0 * std::numeric_limits<double>::epsilon() * A.norm() * std::max<double>(1.0, static_cast<double>(n));
    auto solve_passive = [&](std::vector<Eigen::Index>& P) {
        P.clear();
        for (Eigen::Index j = 0; j < n; ++j) if (passive[static_cast<std::size_t>(j)]) P.push_back(j);
        if (P.empty()) return Eigen::VectorXd();
        Eigen::MatrixXd AP(A.rows(), static_cast<Eigen::Index>(P.size()));
        for (std::size_t k = 0; k < P.size(); ++k) AP.col(static_cast<Eigen::Index>(k)) = A.col(P[k]);
        return Eigen::VectorXd(AP.colPivHouseholderQr().solve(b));
    };
    std::vector<Eigen::Index> P;
    for (int outer = 0; outer < max_iter; ++outer) {
        const Eigen::VectorXd w = A.transpose() * (b - A * x);
        Eigen::Index t = -1;
        double wmax = tol;
        for (Eigen::Index j = 0; j < n; ++j) {
            const auto uj = static_cast<std::size_t>(j);
            if (!passive[uj] && !blocked[uj] && w[j] > wmax) { wmax = w[j]; t = j; }
        }
        if (t < 0) break;
        passive[static_cast<std::size_t>(t)] = true;
        Eigen::VectorXd zP = solve_passive(P);
        const auto pos_t = static_cast<Eigen::Index>(std::find(P.begin(), P.end(), t) - P.begin());
        if (!(zP[pos_t] > 0.0)) {  // t cannot enter: block it, try the next candidate
            passive[static_cast<std::size_t>(t)] = false;
            blocked[static_cast<std::size_t>(t)] = true;
            continue;
        }
        std::fill(blocked.begin(), blocked.end(), false);
        for (int inner = 0; inner < max_iter; ++inner) {
            bool all_pos = true;
            for (std::size_t k = 0; k < P.size(); ++k) if (zP[static_cast<Eigen::Index>(k)] <= 0.0) all_pos = false;
            if (all_pos) {
                x.setZero();
                for (std::size_t k = 0; k < P.size(); ++k) x[P[k]] = zP[static_cast<Eigen::Index>(k)];
                break;
            }
            double alpha = 1.0;
            for (std::size_t k = 0; k < P.size(); ++k) {
                const double zk = zP[static_cast<Eigen::Index>(k)];
                if (zk <= 0.0) alpha = std::min(alpha, x[P[k]] / (x[P[k]] - zk));
            }
            for (std::size_t k = 0; k < P.size(); ++k) x[P[k]] += alpha * (zP[static_cast<Eigen::Index>(k)] - x[P[k]]);
            for (std::size_t k = 0; k < P.size(); ++k)
                if (x[P[k]] <= tol) { x[P[k]] = 0.0; passive[static_cast<std::size_t>(P[k])] = false; }
            zP = solve_passive(P);
            if (P.empty()) break;
        }
    }
    return x;
}

/// Least distance: min |z| subject to G z >= h; nullopt when infeasible (L&H Ch. 23).
[[nodiscard]] inline std::optional<Eigen::VectorXd> ldp(const Eigen::MatrixXd& G, const Eigen::VectorXd& h, double tol = 1e-12)
{
    const Eigen::Index n = G.cols();
    Eigen::MatrixXd system(n + 1, G.rows());
    system.topRows(n) = G.transpose();
    system.row(n) = h.transpose();
    Eigen::VectorXd target = Eigen::VectorXd::Zero(n + 1);
    target[n] = 1.0;
    const Eigen::VectorXd u = nnls(system, target, static_cast<int>(50 * (system.cols() + n + 10)));
    const Eigen::VectorXd r = system * u - target;
    if (r.norm() <= tol || !(std::abs(r[n]) > tol)) return std::nullopt;
    return Eigen::VectorXd(-r.head(n) / r[n]);
}

/// min |A x - b| subject to E x = f and x >= 0; nullopt when the constraints admit no x.
[[nodiscard]] inline std::optional<Eigen::VectorXd> lsei(const Eigen::MatrixXd& A_in, const Eigen::VectorXd& b_in,
                                                         const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    const Eigen::Index n = A_in.cols();
    // a tall A is replaced by its triangular factor: |A x - b|^2 = |R x - (Q^T b)_n|^2 + const
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    if (A_in.rows() > n) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(A_in);
        A = qr.matrixQR().topRows(n).triangularView<Eigen::Upper>();
        b = (qr.householderQ().transpose() * b_in).head(n);
    } else {
        A = A_in;
        b = b_in;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> se(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const auto& sv = se.singularValues();
    const double cut = (sv.size() ? sv[0] : 0.0) * std::max<double>(static_cast<double>(E.rows()), static_cast<double>(n)) *
                       std::numeric_limits<double>::epsilon();
    Eigen::Index rank = 0;
    for (Eigen::Index i = 0; i < sv.size(); ++i) if (sv[i] > cut) ++rank;
    // minimum-norm particular solution and the null space of E
    Eigen::VectorXd particular = Eigen::VectorXd::Zero(n);
    for (Eigen::Index i = 0; i < rank; ++i)
        particular += se.matrixV().col(i) * (se.matrixU().col(i).dot(f) / sv[i]);
    const Eigen::MatrixXd basis = se.matrixV().rightCols(n - rank);
    if (basis.cols() == 0) {
        if ((particular.array() < -1e-14).any()) return std::nullopt;
        return particular;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> sa(A * basis, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& s = sa.singularValues();
    Eigen::Index keep = 0;
    for (Eigen::Index i = 0; i < s.size(); ++i) if (s[i] > s[0] * 1e-13) ++keep;
    if (keep == 0) return std::nullopt;
    const Eigen::MatrixXd U = sa.matrixU().leftCols(keep), V = sa.matrixV().leftCols(keep);
    const Eigen::VectorXd offset = U.transpose() * (b - A * particular);
    const Eigen::MatrixXd back = basis * V * s.head(keep).cwiseInverse().asDiagonal();  // x = p + back (z + offset)
    const auto z = ldp(back, -particular - back * offset);
    if (!z) return std::nullopt;
    return Eigen::VectorXd(particular + back * (*z + offset));
}

}  // namespace BladeGen::Math
