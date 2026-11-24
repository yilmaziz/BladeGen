/// @file turbolenza/Convergence.hpp
/// @brief Bookkeeping of the nonlinear iteration: change norms, history, stop tests.
///
/// Port of v1 convergence.py (TURBOlenza v1 commit 67599a6):
///     change_norm(new, old) = max |new - old| / max |old|,
///     velocity_change = max |(W_m, W_theta)_new - (W_m, W_theta)_old| / max |(W_m, W_theta)_old|;
/// density_converged: the last density change <= density_tolerance AND its velocity change <=
/// the velocity tolerance (density_tolerance when unset); a NaN (pass 0) never passes.
/// kutta_converged: the last eps_K <= kutta_tolerance; coupled_converged: both on the same pass.
#pragma once

#include "turbolenza/Config.hpp"
#include "turbolenza/Field.hpp"

#include <vector>

namespace BladeGen::Turbolenza {

/// v1 History: one entry per pass (exit_angle and kutta_error in the coupled loop only).
struct History
{
    std::vector<int> kutta_pass;
    std::vector<double> density_change, psi_change, velocity_change, residual;
    std::vector<double> exit_angle, kutta_error;
};

/// v1 change_norm (std::invalid_argument for shapes that differ, non-finite values, old == 0).
[[nodiscard]] double change_norm(const Field& new_field, const Field& old_field);
/// v1 velocity_change.
[[nodiscard]] double velocity_change(const Field& W_m_new, const Field& W_theta_new, const Field& W_m_old,
                                     const Field& W_theta_old);
[[nodiscard]] double velocity_tolerance(const IterationSettings& settings) noexcept;
[[nodiscard]] bool density_converged(const History& history, const IterationSettings& settings) noexcept;
[[nodiscard]] bool kutta_converged(const History& history, const IterationSettings& settings) noexcept;
[[nodiscard]] bool coupled_converged(const History& history, const IterationSettings& settings) noexcept;

}  // namespace BladeGen::Turbolenza
