/// @file turbolenza/GeneralSolver.hpp
/// @brief The density-coupled (and Kutta-coupled) potential-flow solve on a general multi-block mesh:
///        the entry point solve_general, separate from the Aungier solver (Solver.hpp).  Not a port of v1.
///
/// One pass n (the order of the Aungier solver, Solver.hpp, so the two can be compared pass by pass):
///  1 the coefficient mdot' / (rho(n) b) at every node, the rotation source omega (r^2)_x J (fixed);
///  2 the assembled system of GeneralOperator, solved by LinearSolver.hpp's solve (SparseLU / PARDISO);
///  3 psi, its gradient, W_m = mdot' psi_y / (rho b r), W_theta = -mdot' psi_x / (rho b r), W;
///  4 the static state of W by Eos.hpp (rothalpy, isentropic to the inlet totals; per node the radius
///    r for a rotating row: h = I + (omega r)^2 / 2 - W^2 / 2), rho_target, the Mach number; M >= 1
///    anywhere stops the solve (unless stop_at_sonic is off); rho(n+1) = relax_density (lambda_rho);
///  5 change norms and the stop tests of Convergence.hpp (density and velocity change; with the Kutta
///    condition also eps_K, all on one pass); at convergence one more evaluation at rho(n+1) (and
///    beta_out(n+1)) is reported, as in the Aungier solver.
/// The Kutta condition (Aungier 2003 Eq. 5-34): equal relative speeds at the suction- and pressure-side
/// trailing-edge nodes, found by their tags (not by index); the discharge angle enters only as the outlet
/// T_f (no mesh rebuild).  The correction is Aungier 2000 Eqs. 13-57 and 13-59 as in Kutta.cpp,
/// Delta tan(beta_out) = 2 (W_+ - W_-) dm / (S (cos beta_+ + cos beta_-) mean W_m(outlet)), written for a
/// general mesh (ours, documented): W_+ / W_- at the TE node on the +theta- / -theta-facing side of the
/// blade (Aungier's W0 - W1: the passage's eta = 0 wall is the reference blade's +theta-facing side);
/// dm the m distance from each TE node to its upstream wall neighbour (the grid dm on Aungier's grid),
/// averaged; cos beta of the wall tangent from the nodal difference along the wall; S = r (theta_- +
/// pitch - theta_+) the passage width at the TE (both theta on the same blade); mean W_m the
/// trapezoidal integral of W_m over theta along the outlet divided by the outlet's theta extent
/// (Aungier's int_0^1 W_m d eta at the far exit).  On Aungier's grid this is Kutta.cpp's correction
/// except cos beta (there from the C2 curve's tangent).
///
/// Rotation sign (the cross-module contract of the quasi-3D work, "frames and signs" 0.1 and 0.3): the
/// operator's omega is the book's, omega_book, with W_theta = C_theta - omega r (Aungier Eqs. 3-1, 5-11:
/// omega > 0 turns toward +theta).  TURBOlenza's frame (beta from the meridional direction toward +theta,
/// eta = 0 the suction wall; a positive rpm with axial inflow gives a positive relative inlet angle) has
/// theta_book = theta_T and omega_book = -Omega, Omega = 2 pi rpm / 60.  FlowConditions::angular_speed is
/// Omega (positive for positive rpm) and RotationSense::AgainstTheta (the default) applies the contract:
/// omega = -angular_speed in the source f = omega d(r^2)/dx and in the rothalpy I = cp Tt - omega r1
/// C_theta1 (C_theta1 = inlet_absolute_tangential_velocity in the T = book frame), so the relative total
/// enthalpy I + (omega r)^2 / 2 at every node's r.  RotationSense::AlongTheta takes omega = +angular_speed
/// (angular_speed is then omega_book itself).  r and b are plain per-node inputs of the mesh (the
/// quasi-3D loop will give r(m), b(m) of each stream sheet).
#pragma once

#include "turbolenza/Config.hpp"
#include "turbolenza/Convergence.hpp"
#include "turbolenza/GeneralOperator.hpp"
#include "turbolenza/LinearSolver.hpp"

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

enum class RotationSense { AgainstTheta, AlongTheta };

/// omega in the +theta sense of the operator for the user's angular speed (see the file comment).
[[nodiscard]] double theta_angular_speed(double angular_speed, RotationSense sense) noexcept;

struct GeneralConfig
{
    FlowConditions flow;           ///< inlet Mach, angle, totals; exit_angle_guess_deg = beta_out(0) (or the fixed angle)
    PerfectGas gas;
    IterationSettings iteration;   ///< tolerances, pass limits, relaxations, stop_at_sonic
    RotationSense rotation = RotationSense::AgainstTheta;
    bool kutta = true;             ///< false: the outlet angle stays exit_angle_guess_deg
    std::optional<double> mass_flow_per_depth;  ///< kg/(m s); needed for a rotating row, else from the inlet Mach
    LinearSolverKind linear_solver = LinearSolverKind::SparseLU;
    int threads = 1;
};

/// The numbers of one pass.
struct GeneralPass
{
    int index = 0;
    bool final = false;
    double density_change = 0.0, velocity_change = 0.0, psi_change = 0.0;
    double mach_max = 0.0, linear_residual = 0.0;
    double exit_angle = 0.0, exit_angle_next = 0.0;  ///< rad
    double kutta_residual = 0.0, kutta_normalized = 0.0, tan_correction = 0.0;  ///< NaN without the Kutta loop
    double inlet_mass_flow_error = 0.0, outlet_mass_flow_error = 0.0;  ///< |int rho W_m b r dtheta / mdot' - 1|
};

/// The two trailing-edge nodes of the Kutta condition and their geometry.
struct GeneralTrailingEdge
{
    int suction = -1, pressure = -1;   ///< nodes
    int plus = -1, minus = -1;         ///< the +theta- and -theta-facing sides (one of suction / pressure each)
    double passage_width = 0.0;        ///< S = r (theta_- + pitch - theta_+) [m]
    double dm = 0.0;                   ///< mean m distance to the upstream wall neighbours [m]
    double cos_plus = 0.0, cos_minus = 0.0;
};

struct GeneralKutta
{
    double W_ss = 0.0, W_ps = 0.0, residual = 0.0, normalized = 0.0;  ///< R_K = W_ss - W_ps, eps_K
    double wall_difference = 0.0;   ///< W_+ - W_-
    double outlet_mean_W_m = 0.0;
    double delta_tan = 0.0;
};

struct GeneralSolution
{
    std::shared_ptr<const GeneralOperator> op;
    std::vector<double> psi, density, density_target, W_m, W_theta, W, mach, pressure, temperature;
    double mass_flow_per_depth = 0.0;
    double inlet_angle = 0.0, exit_angle = 0.0;  ///< rad
    std::optional<GeneralKutta> kutta;
    History history;
    std::vector<GeneralPass> passes;
    bool converged = false;
    int iterations = 0;
};

/// The solve did not converge, went sonic or a module failed in a pass; solution: the last pass (or null).
class GeneralSolverError : public std::runtime_error
{
public:
    GeneralSolverError(const std::string& what, std::shared_ptr<const GeneralSolution> solution_)
        : std::runtime_error(what), solution(std::move(solution_))
    {
    }
    std::shared_ptr<const GeneralSolution> solution;
};

/// The trailing-edge nodes from the mesh tags (GeneralOperatorError unless exactly one suction and one
/// pressure trailing-edge tag on different nodes, each with an upstream wall neighbour).
[[nodiscard]] GeneralTrailingEdge general_trailing_edge(const GeneralOperator& op);

/// The Kutta residual and the correction of tan(beta_out) at the given flow.
[[nodiscard]] GeneralKutta general_kutta(const GeneralOperator& op, const GeneralTrailingEdge& te,
                                         const std::vector<double>& W_m, const std::vector<double>& W);

/// The solve (see the file comment).  The mesh needs inlet and outlet segments; config.kutta needs the
/// trailing-edge tags.
[[nodiscard]] GeneralSolution solve_general(std::shared_ptr<const GeneralOperator> op, const GeneralConfig& config);
[[nodiscard]] GeneralSolution solve_general(const GeneralMesh& mesh, const GeneralConfig& config);

}  // namespace BladeGen::Turbolenza
