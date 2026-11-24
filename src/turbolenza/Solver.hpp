/// @file turbolenza/Solver.hpp
/// @brief Runs the method by calling the other modules in order: the density pass, the coupled
///        density + Kutta loop and the fixed-angle loop.
///
/// Port of v1 solver.py, production paths (TURBOlenza v1 commit 67599a6).  The density pass n:
///  3 coefficients from the resident rho(n)   4 F Psi = Q   5 Psi, psi(n)   6 W_m, W_theta, W
///  7 static state of W(n), rho_target = rho_EOS(W(n)); Mach >= 1 anywhere: TransonicFlowError
///    (unless stop_at_sonic is off); rho(n+1) = (1 - lambda_rho) rho(n) + lambda_rho rho_target
///  8 change norms and the density stop test (pass 0 has no previous W: never converged).
/// The coupled pass n (solve / solve_coupled): grid(n) = rebuild_exit(grid(0), beta_out(n)) (only
/// the exit stations' nodes may move, the grid must stay admissible, else ExitRebuildError),
/// metrics and boundary relations rebuilt, the density pass, the Kutta residual and Aungier's
/// correction, beta_out(n+1); converged when the density, velocity and Kutta tests hold on one
/// pass.  At convergence one more evaluation is made at (rho(n+1), beta_out(n+1)) and reported
/// (final = true).  No convergence within the pass limit: DensityNotConvergedError (fixed angle)
/// or CoupledNotConvergedError (coupled), carrying the last pass.
/// Errors of the called modules keep their type; the last completed pass and the failed pass's
/// index go to SolveOptions::failure (v1 attaches them to the exception as .solution and
/// .failed_pass).  Observers (RunObserver) get every PassRecord once, the final evaluation once
/// (marked final); an observer that throws is switched off, nothing it does reaches the solve.
/// SolveOptions::cancel is polled between passes: a cancelled solve ends like an unconverged one
/// (CancelledError, the last completed pass).
/// v2 additions: the thread count of the per-node loops (and PARDISO), the linear solver kind,
/// the grid operations as a seam (tests inject the v1 grids), and the solve starting from a
/// grid that is already built (Case.cpp passes build_case's grid instead of building it again).
#pragma once

#include "turbolenza/BoundaryConditions.hpp"
#include "turbolenza/Case.hpp"
#include "turbolenza/Coefficients.hpp"
#include "turbolenza/Config.hpp"
#include "turbolenza/Convergence.hpp"
#include "turbolenza/Eos.hpp"
#include "turbolenza/Flow.hpp"
#include "turbolenza/Kutta.hpp"
#include "turbolenza/LinearSolver.hpp"
#include "turbolenza/Metrics.hpp"

#include <atomic>
#include <exception>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 LIMITATION, COUPLED_LIMITATION, MASS_FLOW_* and MASS_FLOW_CONSISTENCY.
extern const char* const kLimitation;
extern const char* const kCoupledLimitation;
extern const char* const kMassFlowFromInlet;
extern const char* const kMassFlowOverride;
extern const char* const kMassFlowChecked;
inline constexpr double kMassFlowConsistency = 1e-12;

/// v1 Solution: everything of the reported evaluation (one pass).
struct Solution
{
    std::shared_ptr<const Grid> grid;
    GridMetrics metrics;
    Field psi;
    FlowField flow;
    History history;
    Field density, density_target;
    StaticState static_state;
    Field mach;
    CoefficientFields coefficients;
    LinearSolveResult linear;
    std::vector<PassRecord> passes;
    bool converged = false;
    int iterations = 0;
    double mass_flow_per_depth = 0.0, c_ref = 0.0, density_relaxation = 0.0;
    double exit_angle = 0.0;
    std::optional<KuttaResidual> kutta;
    double kutta_residual = 0.0, normalized_kutta_residual = 0.0;
    std::optional<double> kutta_relaxation;
    std::string limitation;
    std::optional<InletState> inlet;
    std::string mass_flow_source;
};

/// v1 SolverInputError.
class SolverInputError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 DensityIterationError: the iteration stopped without a converged solution; solution is the
/// last pass (converged false), or null.
class DensityIterationError : public std::runtime_error
{
public:
    DensityIterationError(const std::string& what, std::shared_ptr<const Solution> solution_)
        : std::runtime_error(what), solution(std::move(solution_))
    {
    }
    std::shared_ptr<const Solution> solution;
};
class DensityNotConvergedError : public DensityIterationError
{
public:
    using DensityIterationError::DensityIterationError;
};
class CoupledNotConvergedError : public DensityIterationError
{
public:
    using DensityIterationError::DensityIterationError;
};
class TransonicFlowError : public DensityIterationError
{
public:
    using DensityIterationError::DensityIterationError;
};
/// v2: the run was cancelled between passes (RunResult::error "Cancelled: ...").
class CancelledError : public DensityIterationError
{
public:
    using DensityIterationError::DensityIterationError;
};
/// v1 ExitRebuildError.
class ExitRebuildError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

/// The v1 class name of an exception of the port ("TransonicFlowError", "LinearSolveError", ...;
/// "Cancelled" for CancelledError), else "RuntimeError"-like std names.
[[nodiscard]] std::string error_type_name(const std::exception& error);

/// The fields of one pass for an inspector (v1 PassFields): views valid during the call only.
struct PassFields
{
    const PassRecord& record;
    std::shared_ptr<const Grid> grid;
    const GridMetrics& metrics;
    const Field& density_resident;
    const Field& density_target;
    const Field& density_next;
    const CoefficientFields& coefficients;
    const LinearSolveResult& linear;
    const FlowField& flow;
    const StaticState& static_state;
    const Field& mach;
    double mass_flow_per_depth = 0.0, c_ref = 0.0;
};

/// The grid operations of the coupled loop; empty functions mean Grid.hpp's rebuild_exit and
/// validate_uniform_grid_topology (a seam for tests that inject v1's grids).
struct GridOps
{
    std::function<Grid(const Grid&, double)> rebuild_exit;
    std::function<void(const Grid&)> validate_topology;
};

/// Filled when a called module fails inside a run (v1's .solution and .failed_pass).
struct SolveFailure
{
    std::shared_ptr<const Solution> solution;
    int failed_pass = -1;
};

struct SolveOptions
{
    RunObserver* observer = nullptr;
    const std::atomic<bool>* cancel = nullptr;
    std::optional<Field> initial_density;   ///< the resident density of pass 0 (default eos.initial_density)
    int threads = 1;                        ///< 1 .. min(kMaxThreads, hardware threads)
    LinearSolverKind linear_solver = LinearSolverKind::SparseLU;
    GridOps grid_ops;
    std::function<void(const PassFields&)> inspect;  ///< every pass's fields (diagnostics, tests)
    SolveFailure* failure = nullptr;
};

/// v1 solve: make_grid (admissible) at the inlet angle and exit_angle_guess_deg, then the coupled loop.
[[nodiscard]] Solution solve(const Cascade& cascade, const SolverConfig& config, double c_ref,
                             std::optional<double> mass_flow_per_depth = std::nullopt, const SolveOptions& options = {});
/// v1 solve_density_coupled: make_grid at the fixed angle, then solve_at_exit_angle.
[[nodiscard]] Solution solve_density_coupled(const Cascade& cascade, const SolverConfig& config, double c_ref,
                                             std::optional<double> mass_flow_per_depth = std::nullopt,
                                             const SolveOptions& options = {});
/// v1 solve_at_exit_angle: the density loop at the angle grid was built for (no Kutta condition).
[[nodiscard]] Solution solve_at_exit_angle(std::shared_ptr<const Grid> grid, const SolverConfig& config, double c_ref,
                                           std::optional<double> mass_flow_per_depth = std::nullopt,
                                           const SolveOptions& options = {});
/// v1 _solve_coupled: the coupled loop from a built grid, whose exit angle is beta_out(0).
[[nodiscard]] Solution solve_coupled(std::shared_ptr<const Grid> grid, const SolverConfig& config, double c_ref,
                                     std::optional<double> mass_flow_per_depth = std::nullopt,
                                     const SolveOptions& options = {});

}  // namespace BladeGen::Turbolenza
