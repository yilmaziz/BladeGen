// Port of v1 convergence.py (see Convergence.hpp).  The maxima are serial (exact in any order).
#include "turbolenza/Convergence.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BladeGen::Turbolenza {

double change_norm(const Field& a, const Field& b)
{
    if (!a.same_shape(b))
        throw std::invalid_argument("the two fields have shapes (" + std::to_string(a.n_m) + ", " +
                                    std::to_string(a.n_eta) + ") and (" + std::to_string(b.n_m) + ", " +
                                    std::to_string(b.n_eta) + ")");
    double scale = 0.0, change = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        if (!std::isfinite(a.v[k]) || !std::isfinite(b.v[k]))
            throw std::invalid_argument("a field compared between passes is not finite");
        scale = std::max(scale, std::abs(b.v[k]));
        change = std::max(change, std::abs(a.v[k] - b.v[k]));
    }
    if (!(scale > 0.0)) throw std::invalid_argument("the previous field is zero everywhere: its relative change is undefined");
    return change / scale;
}

double velocity_change(const Field& W_m_new, const Field& W_theta_new, const Field& W_m_old, const Field& W_theta_old)
{
    if (!W_m_new.same_shape(W_m_old) || !W_theta_new.same_shape(W_theta_old) || !W_m_new.same_shape(W_theta_new))
        throw std::invalid_argument("the two velocity fields have different shapes");
    double scale = 0.0, change = 0.0;
    for (std::size_t k = 0; k < W_m_new.size(); ++k) {
        const double a = W_m_new.v[k], b = W_theta_new.v[k], c = W_m_old.v[k], d = W_theta_old.v[k];
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d))
            throw std::invalid_argument("a velocity field compared between passes is not finite");
        scale = std::max(scale, std::hypot(c, d));
        change = std::max(change, std::hypot(a - c, b - d));
    }
    if (!(scale > 0.0))
        throw std::invalid_argument("the previous velocity is zero everywhere: its relative change is undefined");
    return change / scale;
}

double velocity_tolerance(const IterationSettings& s) noexcept
{
    return s.velocity_tolerance ? *s.velocity_tolerance : s.density_tolerance;
}

bool density_converged(const History& h, const IterationSettings& s) noexcept
{
    if (h.density_change.empty() || h.velocity_change.empty()) return false;
    const double d = h.density_change.back(), w = h.velocity_change.back();
    if (std::isnan(d) || std::isnan(w)) return false;
    return d <= s.density_tolerance && w <= velocity_tolerance(s);
}

bool kutta_converged(const History& h, const IterationSettings& s) noexcept
{
    if (h.kutta_error.empty()) return false;
    const double e = h.kutta_error.back();
    return !std::isnan(e) && e <= s.kutta_tolerance;
}

bool coupled_converged(const History& h, const IterationSettings& s) noexcept
{
    return density_converged(h, s) && kutta_converged(h, s);
}

}  // namespace BladeGen::Turbolenza
