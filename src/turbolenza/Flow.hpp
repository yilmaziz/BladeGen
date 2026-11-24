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
///    on the blade walls (Dirichlet).  Not in v1: with boundary conditions built with a
///    Transpiration (the equivalent inviscid flow of the boundary-layer coupling) the wall psi
///    varies along m; psi_m on the walls is then d(psi_wall)/dm, the same run differences along
///    the wall, and the periodic continuation carries each station's jump.  Where that jump varies
///    along the exit (the marched wake: a source sheet on the wake cut), psi_m of the two periodic
///    images differs by d(jump)/dm; both take their mean, so W and the density stay periodic (ours).
/// Mass flow (station_mass_flow, a diagnostic): the trapezoidal integral over eta of rho W_m S.
/// Inlet closure (inlet_state; v1 stationary only): W1 the root of W / a(W) = M1 by bisection to
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
    // not in v1: a rotating row (mach, flow_angle, velocity and dynamic_pressure are relative)
    bool rotating = false;
    double blade_speed = 0.0;                   ///< U = omega r_in
    double absolute_tangential_velocity = 0.0;  ///< C_theta
    double relative_total_pressure = 0.0;       ///< P0' = p (T0' / T)^(gamma / (gamma - 1)), T0' = T + W^2 / (2 cp)
    /// The total pressure of the pressure coefficients: P0' of a rotating row, else total_pressure.
    [[nodiscard]] double reference_total_pressure() const noexcept { return rotating ? relative_total_pressure : total_pressure; }
};

/// v1 inlet_state (InletStateError for what it cannot close).  Not in v1: a rotating row
/// (angular_speed != 0) closes with its rothalpy I = cp Tt - U C_theta at stream_surface_radius
/// (Eos.hpp): the inlet Mach number and angle are the RELATIVE ones (velocity_triangle derives
/// them), W1 the root of W / a(W) = M1_rel on the relative total enthalpy H' = I + U^2 / 2.
[[nodiscard]] InletState inlet_state(const FlowConditions& conditions, const PerfectGas& gas, double inlet_pitch);
/// v1 passage_mass_flow: rho1 W_m1 S1.
[[nodiscard]] double passage_mass_flow(double density, double meridional_velocity, double inlet_pitch) noexcept;

/// Not in v1: the inlet velocity triangle of a rotating row from the absolute inlet state.
/// The sign convention (ours; TURBOlenza's frame, every angle and tangential component positive
/// toward +theta = +y, the side of the reference blade's suction surface, eta = 0 at the suction wall):
///   omega = -2 pi rpm / 60 (Aungier's omega of Eqs. 5-11 .. 5-21 and 5-51, U = omega r):
///   positive rpm moves the blades toward -theta, the way their pressure surfaces face (a compressor
///   rotor, whose blades push the gas with their pressure side);
///   T1 = T0 / (1 + (gamma - 1) M1^2 / 2), a1 = sqrt(gamma R T1), C = M1 a1 (isentropic, absolute);
///   C_m = C cos(alpha1), C_theta = C sin(alpha1), U = omega r_in, W_theta = C_theta - U,
///   beta1 = atan(W_theta / C_m), W = sqrt(C_m^2 + W_theta^2), M1_rel = W / a1, I = cp T0 - U C_theta.
/// So rpm > 0 with alpha1 = 0 gives W_theta = 2 pi rpm r_in / 60 > 0, beta1 > 0; pre-swirl in the
/// direction of rotation is a negative alpha1.  InletStateError for inputs that make no state.
struct VelocityTriangle
{
    double rpm = 0.0, omega = 0.0, radius = 0.0;          ///< omega [rad/s] signed as above; radius [m]
    double alpha = 0.0, mach_absolute = 0.0;              ///< absolute angle [rad] and Mach number
    double total_pressure = 0.0, total_temperature = 0.0; ///< absolute totals
    double static_temperature = 0.0, static_pressure = 0.0, static_density = 0.0, sound_speed = 0.0;
    double C = 0.0, C_m = 0.0, C_theta = 0.0, U = 0.0, W_theta = 0.0, W = 0.0;
    double beta = 0.0, mach_relative = 0.0;               ///< relative angle [rad] and Mach number
    double rothalpy = 0.0;                                ///< I = cp T0 - U C_theta [J/kg]
};
[[nodiscard]] VelocityTriangle velocity_triangle(double rpm, double alpha1_deg, double mach_absolute, double total_pressure,
                                                 double total_temperature, const PerfectGas& gas, double radius);

/// v1 velocity_components (FlowRecoveryError for broken or mismatched inputs).  Not in v1: on a
/// stream sheet mu = mdot' / (rho b) [Eqs. 5-18, 5-19], b = metrics.thickness.
[[nodiscard]] FlowField velocity_components(const Field& psi, const Field& density,
                                            const CoefficientFields& coefficients, const GridMetrics& metrics,
                                            const BoundaryConditions& boundary_conditions);

/// v1 station_mass_flow: (n_m) trapezoidal integrals of rho W_m S over eta (on a stream sheet
/// rho W_m S b, Eq. 5-17: each equals mdot' of the passage).
[[nodiscard]] std::vector<double> station_mass_flow(const FlowField& flow, const Field& density,
                                                    const CoefficientFields& coefficients, const GridMetrics& metrics);

}  // namespace BladeGen::Turbolenza
