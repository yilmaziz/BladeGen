/// @file turbolenza/Flow.hpp
/// @brief From the stream function to the flow at the nodes; the far-inlet state.
///
/// Port of v1 flow.py (TURBOlenza v1 commit 67599a6; Aungier 2003 Eqs. 5-18, 5-19):
///     W_m = mu psi_eta / S,   W_theta = mu [(tan(beta) / S) psi_eta - psi_m],   mu = mdot' / rho,
/// S = L_ref S*, psi_m in the physical m (dm = L_ref grid.dm).  The derivatives:
///  * psi_eta: central inside; one-sided into the passage at the blade walls; at both periodic
///    images central across the line with the continuation psi[i, -1] = psi[i, J-1] - 1,
///    psi[i, J+1] = psi[i, 1] + 1, i.e. (psi[i, 1] - psi[i, J-1] + 1) / (2 deta);
///  * psi_m: by m region (metrics.py's runs and stencils: central inside a run, one-sided at its
///    ends, so backward at f - 1, forward at f, backward at l, forward at l + 1), and exactly 0
///    on the blade walls (Dirichlet).
/// Mass flow (station_mass_flow, a diagnostic): the trapezoidal integral over eta of rho W_m S.
/// Inlet closure (inlet_state, stationary only): W1 the root of W / a(W) = M1 by bisection to
/// adjacent doubles on eos's sound speed; mdot' = rho1 W1 cos(beta1) S1.
/// The per-node maps run on the solve's threads.
#pragma once

#include "turbolenza/BoundaryConditions.hpp"
#include "turbolenza/Coefficients.hpp"
#include "turbolenza/Config.hpp"
#include "turbolenza/Eos.hpp"
#include "turbolenza/Field.hpp"
#include "turbolenza/Metrics.hpp"

#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 FlowRecoveryError.
class FlowRecoveryError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 InletStateError.
class InletStateError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 FlowField.
struct FlowField
{
    Field psi_m, psi_eta, W_m, W_theta, W2, W;
};

/// v1 InletState (SI with SI totals; mdot' in kg/(m s)).
struct InletState
{
    double mach = 0.0, flow_angle = 0.0, total_pressure = 0.0, total_temperature = 0.0;
    double static_temperature = 0.0, static_pressure = 0.0, static_density = 0.0, sound_speed = 0.0;
    double velocity = 0.0, meridional_velocity = 0.0, tangential_velocity = 0.0, dynamic_pressure = 0.0;
    double pitch = 0.0, mass_flow_per_depth = 0.0;
};

/// v1 inlet_state (InletStateError for what it cannot close).
[[nodiscard]] InletState inlet_state(const FlowConditions& conditions, const PerfectGas& gas, double inlet_pitch);
/// v1 passage_mass_flow: rho1 W_m1 S1.
[[nodiscard]] double passage_mass_flow(double density, double meridional_velocity, double inlet_pitch) noexcept;

/// v1 velocity_components (FlowRecoveryError for broken or mismatched inputs).
[[nodiscard]] FlowField velocity_components(const Field& psi, const Field& density,
                                            const CoefficientFields& coefficients, const GridMetrics& metrics,
                                            const BoundaryConditions& boundary_conditions);

/// v1 station_mass_flow: (n_m) trapezoidal integrals of rho W_m S over eta.
[[nodiscard]] std::vector<double> station_mass_flow(const FlowField& flow, const Field& density,
                                                    const CoefficientFields& coefficients, const GridMetrics& metrics);

}  // namespace BladeGen::Turbolenza
