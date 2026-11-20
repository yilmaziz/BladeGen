/// @file math/LeastSquares.hpp
/// @brief Bound-constrained nonlinear least squares (Levenberg–Marquardt).
///
/// Minimises 0.5 * ||r(x)||^2 subject to lb <= x <= ub.  The method is the
/// classical Levenberg–Marquardt iteration with Marquardt's diagonal scaling and
/// Nielsen's damping update (H. B. Nielsen, "Damping parameter in Marquardt's
/// method", IMM-REP-1999-05, DTU), applied in the scaled variables
/// z = (x - lb) / (ub - lb).  Bounds are handled by projecting each trial step
/// onto the box and freezing variables that sit on an active bound with the
/// gradient pointing outward (a projected / active-set LM).
///
/// The Jacobian is formed by central differences unless the caller provides one.
#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>

namespace BladeGen::Math {

struct LMOptions
{
    int max_iterations = 200;
    double ftol = 1e-12;   ///< stop when the relative cost decrease is below this
    double xtol = 1e-10;   ///< stop when the scaled step is below this
    double gtol = 1e-12;   ///< stop when the scaled projected gradient is below this
    double fd_step = 1e-6; ///< central-difference step in scaled variables
    double lambda0 = 1e-3; ///< initial damping relative to max(diag(J^T J))
};

struct LMResult
{
    Eigen::VectorXd x;
    Eigen::VectorXd residuals;
    double cost = 0.0;      ///< 0.5 * ||r||^2
    int iterations = 0;
    int evaluations = 0;
    std::string status;
};

using ResidualFn = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

[[nodiscard]] inline LMResult levenberg_marquardt(const ResidualFn& fun, Eigen::VectorXd x0,
                                                  const Eigen::VectorXd& lb, const Eigen::VectorXd& ub,
                                                  const LMOptions& opt = {})
{
    const Eigen::Index n = x0.size();
    if (lb.size() != n || ub.size() != n) throw std::invalid_argument("levenberg_marquardt: bound size");
    for (Eigen::Index i = 0; i < n; ++i)
        if (!(ub[i] > lb[i])) throw std::invalid_argument("levenberg_marquardt: need lb < ub");
    const Eigen::VectorXd span = ub - lb;
    auto to_x = [&](const Eigen::VectorXd& z) { return Eigen::VectorXd(lb + span.cwiseProduct(z)); };
    Eigen::VectorXd z = ((x0 - lb).cwiseQuotient(span)).cwiseMax(0.0).cwiseMin(1.0);

    LMResult res;
    auto eval = [&](const Eigen::VectorXd& zz) {
        ++res.evaluations;
        return fun(to_x(zz));
    };
    Eigen::VectorXd r = eval(z);
    const Eigen::Index m = r.size();
    double cost = 0.5 * r.squaredNorm();

    auto jacobian = [&](const Eigen::VectorXd& zz) {
        Eigen::MatrixXd J(m, n);
        for (Eigen::Index j = 0; j < n; ++j) {
            const double h = opt.fd_step;
            Eigen::VectorXd zp = zz, zm = zz;
            zp[j] = std::min(1.0, zz[j] + h);
            zm[j] = std::max(0.0, zz[j] - h);
            const double dz = zp[j] - zm[j];
            J.col(j) = (eval(zp) - eval(zm)) / dz;
        }
        return J;
    };

    double lambda = -1.0, nu = 2.0;
    res.status = "max iterations";
    for (int it = 0; it < opt.max_iterations; ++it) {
        res.iterations = it + 1;
        const Eigen::MatrixXd J = jacobian(z);
        Eigen::VectorXd g = J.transpose() * r;
        // active set: on a bound with the descent direction pointing outside
        Eigen::Array<bool, Eigen::Dynamic, 1> free(n);
        for (Eigen::Index i = 0; i < n; ++i)
            free[i] = !((z[i] <= 0.0 && g[i] > 0.0) || (z[i] >= 1.0 && g[i] < 0.0));
        double gmax = 0.0;
        for (Eigen::Index i = 0; i < n; ++i) if (free[i]) gmax = std::max(gmax, std::abs(g[i]));
        if (gmax <= opt.gtol) { res.status = "gradient tolerance"; break; }

        Eigen::MatrixXd A = J.transpose() * J;
        if (lambda < 0.0) lambda = opt.lambda0 * std::max(A.diagonal().maxCoeff(), 1e-300);
        bool accepted = false;
        for (int inner = 0; inner < 50 && !accepted; ++inner) {
            Eigen::MatrixXd Ad = A;
            for (Eigen::Index i = 0; i < n; ++i) {
                Ad(i, i) += lambda * std::max(A(i, i), 1e-12);
                if (!free[i]) {  // freeze the variable
                    Ad.row(i).setZero(); Ad.col(i).setZero(); Ad(i, i) = 1.0;
                }
            }
            Eigen::VectorXd gf = g;
            for (Eigen::Index i = 0; i < n; ++i) if (!free[i]) gf[i] = 0.0;
            const Eigen::VectorXd step = -Ad.ldlt().solve(gf);
            if (!step.allFinite()) { lambda *= nu; nu *= 2.0; continue; }
            const Eigen::VectorXd zt = (z + step).cwiseMax(0.0).cwiseMin(1.0);
            const Eigen::VectorXd dz = zt - z;
            if (dz.norm() <= opt.xtol * (opt.xtol + z.norm())) {
                res.status = "step tolerance";
                accepted = true;
                it = opt.max_iterations;  // leave the outer loop too
                break;
            }
            const Eigen::VectorXd rt = eval(zt);
            const double cost_t = 0.5 * rt.squaredNorm();
            // gain ratio against the linear model of the projected step
            const double pred = -(g.dot(dz) + 0.5 * dz.dot(A * dz));
            const double rho = (pred > 0.0) ? (cost - cost_t) / pred : -1.0;
            if (std::isfinite(cost_t) && cost_t < cost) {
                const double rel = (cost - cost_t) / std::max(cost, 1e-300);
                z = zt; r = rt; cost = cost_t;
                lambda *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * std::clamp(rho, 0.0, 1.0) - 1.0, 3));
                nu = 2.0;
                accepted = true;
                if (rel <= opt.ftol) { res.status = "cost tolerance"; it = opt.max_iterations; }
            } else {
                lambda *= nu;
                nu *= 2.0;
            }
        }
        if (!accepted) { res.status = "damping exhausted"; break; }
    }
    res.x = to_x(z);
    res.residuals = r;
    res.cost = cost;
    return res;
}

}  // namespace BladeGen::Math
