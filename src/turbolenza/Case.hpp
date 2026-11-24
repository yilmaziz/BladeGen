/// @file turbolenza/Case.hpp
/// @brief The run layer TURBOlenza's window calls: flat settings -> case -> one timed solve ->
///        report, plot data and files.  Port of v1 controller.py (no numerics here: it calls
///        make_cascade, make_grid / build_grid, the inlet state, the solver and postprocess).
///
/// The settings and their defaults are those of v1 scripts/run_full_solver.py (the user's
/// current values, 2026-11-24) minus the blade inputs, which come from BladeGen as a
/// BladeSection; plus the number of threads and the sparse direct solver.
#pragma once

#include "turbolenza/BoundaryLayer.hpp"
#include "turbolenza/Cascade.hpp"
#include "turbolenza/Field.hpp"
#include "turbolenza/Grid.hpp"
#include "turbolenza/GridData.hpp"

#include <atomic>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

enum class SolverMode { CoupledKutta, FixedExit };          ///< v1 "coupled_kutta", "fixed_exit"
enum class ExitAngleStart { BladeTeAngle, Manual };         ///< v1 "blade_te_angle", "manual"
enum class LinearSolverKind { SparseLU, Pardiso };          ///< Eigen SparseLU (serial) or MKL PARDISO (parallel)
/// Not in v1: the boundary layer's coupling to the potential flow (Viscous.hpp).
enum class BlCoupling { OneWay, TwoWay };                   ///< "one_way" (Aungier's analysis), "two_way"
enum class WakeModel { Constant, Marched };                 ///< "constant", "marched"
[[nodiscard]] const char* bl_coupling_name(BlCoupling c) noexcept;
[[nodiscard]] const char* wake_model_name(WakeModel w) noexcept;

/// The most cores the solver may use (the user's safety limit).
inline constexpr int kMaxThreads = 10;

/// v1 controller.CaseSettings (blade inputs replaced by the BladeSection of build_case).
struct CaseSettings
{
    // 0. run control
    std::string case_name = "turbolenza";
    SolverMode solver_mode = SolverMode::CoupledKutta;
    std::optional<double> fixed_exit_angle_deg;      ///< fixed_exit only; empty: kappa_te
    // 1. gas / thermodynamics
    double gamma = 1.4;
    double gas_constant = 287.05;                    ///< J/(kg K)
    double inlet_total_pressure = 101325.0;          ///< Pa, absolute frame
    double inlet_total_temperature = 288.15;         ///< K, absolute frame
    double omega = 0.0;                              ///< rad/s; only 0 is solved (a rotating row is refused)
    // 2. cascade (the blade is the BladeSection)
    double chord = 0.5;                              ///< m: the true chord, the solver's c_ref
    double solidity = 1.0;                           ///< chord / pitch
    double stagger_deg = 33.0;                       ///< from the axial direction
    // 3. inlet flow
    std::optional<double> inlet_mach = 0.65;         ///< far-inlet relative Mach number
    double inlet_flow_angle_deg = 45.0;              ///< from the axial direction, positive toward +y
    std::optional<double> mass_flow_per_depth_override;  ///< kg/(m s); empty: derived from M1
    // 4. grid
    std::optional<std::string> grid_level = std::string("medium");  ///< empty: explicit n_passage x n_eta
    std::optional<int> n_passage;
    std::optional<int> n_eta;
    double upstream_extension_chords = 1.0;
    double downstream_extension_chords = 1.0;
    /// not in v1: the first / last blade nodes inside the passage (Aungier, default) or on the LE
    /// and TE (touching, Fig. 5-5 "outside": an EXPERIMENT; GridData.hpp, EdgeNodes)
    EdgeNodes edge_nodes = EdgeNodes::Inside;
    // 5. density
    double density_tolerance = 1e-5;
    std::optional<double> velocity_tolerance;        ///< empty: density_tolerance
    int max_density_iterations = 5000;
    double density_relaxation = 0.3;                 ///< lambda_rho (0.5 v1 baseline; 1.0 Aungier)
    bool stop_at_sonic = true;
    // 6. Kutta
    double kutta_tolerance = 1e-5;
    int max_kutta_iterations = 6000;
    double kutta_relaxation = 0.3;                   ///< lambda_K (1.0 v1 baseline)
    ExitAngleStart initial_exit_angle_mode = ExitAngleStart::BladeTeAngle;
    std::optional<double> initial_exit_angle_deg;    ///< manual only
    // 7. compute
    int threads = 1;                                 ///< 1 .. min(kMaxThreads, hardware threads)
    LinearSolverKind linear_solver = LinearSolverKind::SparseLU;
    bool measure_runtime = true;
    // 8. boundary layer (not in v1: Aungier 2003 Sec. 5.6, coupled by transpiration; Viscous.hpp)
    bool viscous = false;                            ///< off: the inviscid solve of v1, bit for bit
    BlCoupling bl_coupling = BlCoupling::TwoWay;
    LaminarCfForm laminar_cf_form = LaminarCfForm::Centrifugal2000;
    WakeModel wake_model = WakeModel::Constant;
    double bl_relaxation = 0.5;                      ///< omega_BL, 0 < omega <= 1
    double bl_tolerance = 1e-4;                      ///< max relative change of m on the blade surfaces
    int max_bl_iterations = 30;                      ///< coupled iterations after the inviscid one (0: one way)
    double transition_re_theta = 250.0;              ///< Aungier Eq. 5-126
    double separation_hk = 2.4;                      ///< Aungier 2003 p. 111
    SutherlandLaw sutherland;                        ///< mu(T) of the gas (air)
};

/// v1 CaseError: the settings do not make a case.
class CaseError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 solver.PassRecord: the numbers of one pass (no fields).
struct PassRecord
{
    int index = 0;
    bool final = false;                 ///< the evaluation at the converged state, not a loop pass
    double density_change = 0.0;
    double velocity_change = 0.0;       ///< NaN on pass 0
    double psi_change = 0.0;            ///< NaN on pass 0
    double target_mismatch = 0.0;
    double density_min = 0.0, density_max = 0.0;
    double speed_min = 0.0, speed_max = 0.0;
    double mach_max = 0.0;
    double linear_residual = 0.0;
    double mass_flow_error_max = 0.0;
    double mass_flow_error_blade_median = 0.0;
    double exit_angle = 0.0;            ///< rad
    double kutta_residual = 0.0;        ///< NaN if not formed
    double kutta_normalized = 0.0;      ///< eps_K; NaN if not formed
    double exit_angle_next = 0.0;       ///< rad
    double tan_correction = 0.0;        ///< NaN at a fixed angle
    int bl_iteration = 0;               ///< not in v1: the boundary-layer iteration k (0: the inviscid solve)
};

/// v1 monitor.SolveInfo.
struct SolveInfo
{
    std::string mode;                   ///< "coupled" or "fixed angle"
    double density_tolerance = 0.0, velocity_tolerance = 0.0;
    std::optional<double> kutta_tolerance;
    int max_passes = 0;
    double exit_angle_start = 0.0, inlet_angle = 0.0;
    double mass_flow_per_depth = 0.0;
    std::string mass_flow_source;
    int n_m = 0, n_eta = 0;
    int threads = 1;
    std::string linear_solver;
    std::optional<int> bl_iteration;    ///< not in v1: the boundary-layer iteration of this solve (viscous runs)
};

/// Observer of a solve (v1 monitor.py protocol): observation only, the result is the same with
/// or without it.  update() is called from the solver thread after every pass.
class RunObserver
{
public:
    virtual ~RunObserver() = default;
    virtual void start(const SolveInfo&) {}
    virtual void update(const PassRecord&) {}
    virtual void finish(bool /*converged*/, const std::string& /*error*/) {}
};

/// v1 console line of one pass (monitor.format_record).
[[nodiscard]] std::string format_record(const PassRecord& r);

/// One blade surface's values at the wall nodes of the blade stations (v1 postprocess.BladeSurfaces).
struct SurfaceData
{
    std::vector<double> x_c, m_c, mach, pressure;
    std::vector<double> dunavant_S;     ///< (Pt1 - p) / q1; empty without an inlet state
    std::vector<double> cp_static;      ///< (p - p1) / q1
};

/// v1 flow.InletState (the numbers).
struct InletData
{
    double mach = 0.0, angle = 0.0, speed = 0.0, density = 0.0, pressure = 0.0, temperature = 0.0;
    double sound_speed = 0.0, mass_flow_per_depth = 0.0, total_pressure = 0.0, dynamic_pressure = 0.0;
};

/// The solved (or last) state, in plot-ready form.
struct SolutionData
{
    Grid grid;                          ///< the reported evaluation's grid (its exit angle), chord-normalized
    double c_ref = 1.0;                 ///< true chord [m]: physical lengths = c_ref * grid lengths
    Field psi, density, pressure, temperature, mach, W, W_m, W_theta;
    SurfaceData suction, pressure_side;
    std::optional<InletData> inlet;
    double exit_angle = 0.0;            ///< rad
    double kutta_residual = 0.0, normalized_kutta_residual = 0.0;
    double mass_flow_per_depth = 0.0;
    std::string limitation;
};

struct RuntimeStats
{
    double wall_time_s = 0.0, process_cpu_time_s = 0.0, effective_cpu_cores = 0.0;
    int logical_cpus_available = 0, threads = 1;
    std::string linear_solver;
};

/// Not in v1: one boundary-layer iteration of a viscous run.
struct BlIteration
{
    int k = 0;                          ///< 0: the inviscid solve
    int inviscid_passes = 0;            ///< the passes of its inviscid solve (the final evaluation included)
    double beta_out = 0.0;              ///< rad, of its inviscid solve
    double max_change_m = 0.0;          ///< max |m_BL - m_used| / max |m_BL| on the blade surfaces
    double transition_x_c_ss = 0.0, transition_x_c_ps = 0.0;  ///< NaN: none
    double separation_x_c_ss = 0.0, separation_x_c_ps = 0.0;  ///< turbulent H_k >= the limit; NaN: none
    double laminar_separation_x_c_ss = 0.0, laminar_separation_x_c_ps = 0.0;  ///< Lambda = -12; NaN: none
    double loss = 0.0;                  ///< omega of Eq. 5-139
};

/// One surface's layer at its stations (plot and file data).
struct BoundaryLayerSurfaceData
{
    std::vector<double> x_c, x, W_e, rho_e, T_e, mu, theta, delta_star, delta, H, H_k, c_f, Lambda, Re_theta;
    std::vector<double> mass_deficit;
    std::vector<std::string> regime;
    double transition_x_c = 0.0, separation_x_c = 0.0, laminar_separation_x_c = 0.0;  ///< NaN: none
};

/// The viscous part of a run (RunResult::viscous).
struct ViscousData
{
    std::string coupling, laminar_cf_form, wake_model;
    bool converged = false;
    std::vector<BlIteration> history;
    /// The reported layers: those of BL iteration layers_iteration (-1: none were computed).  When a later
    /// inviscid solve failed they are the last layers computed, of the previous iteration's solution
    /// (layers_of_last_solve false), not of the solution the run reports.
    int layers_iteration = -1;
    bool layers_of_last_solve = true;
    BoundaryLayerSurfaceData suction, pressure_side;
    std::vector<double> wake_m_c, wake_mass_deficit;   ///< the exit stations' m_wake asked by the reported layers
    std::vector<double> wake_theta, wake_delta_star, wake_H;  ///< marched model only
    double loss_omega = 0.0, loss_Theta = 0.0, loss_Delta_star = 0.0;
    double te_jump_mismatch = 0.0;      ///< |jump(first exit) - (psi_up - psi_low)(last blade)| of the solve's transpiration
    double inviscid_exit_angle = 0.0;   ///< rad, k = 0
    SurfaceData inviscid_suction, inviscid_pressure_side;  ///< the loading of k = 0
};

/// A built case (v1 controller.Case).
struct Case
{
    CaseSettings settings;
    BladeSection section;
    Cascade cascade;
    Grid grid;                          ///< the first grid (coupled) or the fixed-angle grid
    std::optional<StandardGridInfo> standard_info;  ///< a standard level: what build_grid found
    double inlet_angle = 0.0, exit_angle_start = 0.0;  ///< rad
    std::optional<InletData> inlet;
    double mass_flow_per_depth = 0.0;
    std::string mass_flow_source;
};

/// The numbers of v1 postprocess.SolutionSummary for the report and the run summary (v2
/// addition: RunResult carries them, so the report needs no Solution).
struct RunSummary
{
    std::string limitation;
    std::optional<double> inlet_mach;
    double inlet_angle_deg = 0.0, discharge_angle_deg = 0.0, exit_velocity_angle_deg = 0.0;
    double exit_angle_mismatch_deg = 0.0, turning_deg = 0.0;
    double surface_mach_max = 0.0, surface_mach_min = 0.0;
    double ss_peak_mach = 0.0, ss_peak_x_c = 0.0, ps_peak_mach = 0.0, ps_peak_x_c = 0.0;
    double te_ss_mach = 0.0, te_ps_mach = 0.0;
    double kutta_residual = 0.0, normalized_kutta_residual = 0.0;
    double mach_max = 0.0;
    int mach_max_i = 0, mach_max_j = 0;
    double density_min = 0.0, density_max = 0.0, pressure_min = 0.0, pressure_max = 0.0;
    double temperature_min = 0.0, temperature_max = 0.0;
    double mass_flow_per_depth = 0.0;
    std::string mass_flow_source;
};

/// v1 controller.RunResult.
struct RunResult
{
    std::string case_name;
    bool converged = false;
    int iterations = 0;
    std::string error;                  ///< empty: converged; else "<ErrorType>: message"
    std::vector<PassRecord> records;    ///< every pass, the final evaluation last (final = true)
    std::optional<SolutionData> solution;  ///< the reported evaluation, or the last completed pass
    RuntimeStats runtime;
    std::string started;                ///< ISO time of the solve start
    // v2 additions
    std::string status;                 ///< v1 RunResult.status: "converged" or "failed: <ErrorType>"
    std::string mode;                   ///< v1 solver_mode string: "coupled_kutta" or "fixed_exit"
    bool stop_at_sonic = true;          ///< the case's sonic stop (for the report)
    std::optional<int> failed_pass;     ///< a called module failed in this pass (v1 error.failed_pass)
    std::optional<RunSummary> summary;  ///< the reported evaluation's numbers (empty without a solution
                                        ///< or when its walls cannot be labeled)
    std::string postprocess_error;      ///< why summary / surfaces are empty, if they are
    std::optional<double> chord_reynolds;  ///< rho1 W1 c / mu(T1) of the inlet state (Sutherland), if any
    std::optional<ViscousData> viscous;    ///< a viscous run: its layers and the coupling history
};

/// v1 build_case: cascade, grid (admissibility checked before any solve), inlet state and the
/// derived values.  Throws CaseError (with the production message) when the settings do not make
/// a case.
[[nodiscard]] Case build_case(const CaseSettings& settings, const BladeSection& section);

/// v1 run_case: the timed solve; failures are returned in RunResult::error, never swallowed.
/// `cancel` (may be null) is polled between passes: a cancelled solve ends like an unconverged one
/// ("Cancelled: ..."), with the last completed pass reported.
[[nodiscard]] RunResult run_case(const Case& c, RunObserver* observer = nullptr, const std::atomic<bool>* cancel = nullptr);

/// v1 format_case_summary / format_inlet_state / format_grid_summary / format_run_report.
[[nodiscard]] std::string format_case_summary(const Case& c);
[[nodiscard]] std::string format_grid_summary(const Case& c);
[[nodiscard]] std::string format_run_report(const RunResult& r);

/// v1 export_run without the figures: surface_mach.csv, surface_pressure.csv,
/// convergence_history.csv, fields.csv (every node: m, y, psi, rho, p, T, M, W, W_m, W_theta)
/// and run_summary.json into `directory` (created; an existing file is overwritten only when
/// `overwrite`).  Returns the names of the files written, relative to `directory`.
std::vector<std::string> export_run(const RunResult& r, const Case& c, const std::string& directory, bool overwrite = false);

/// The chord Reynolds number rho1 W1 c / mu(T1) of an inlet state at Mach number `mach` with the
/// totals P0, T0 (stationary perfect gas, isentropic: T1 = T0 / (1 + (gamma - 1)/2 M^2),
/// p1 = P0 (T1/T0)^(gamma/(gamma - 1))), mu by Sutherland at the inlet static temperature.
[[nodiscard]] double chord_reynolds_number(double mach, double total_pressure, double total_temperature, double chord,
                                           double gamma, double gas_constant, const SutherlandLaw& law = {});
/// The total pressure P0 that gives the chord Reynolds number `reynolds` at Mach `mach` with the
/// total temperature T0: the same relations solved for P0 (closed form; the tunnel holds P0 and
/// T0).  CaseError for inputs that do not make a state.
[[nodiscard]] double total_pressure_for_reynolds(double reynolds, double mach, double total_temperature, double chord,
                                                 double gamma, double gas_constant, const SutherlandLaw& law = {});

/// Threads the machine offers, capped at kMaxThreads.
[[nodiscard]] int max_threads();
/// True when the build has MKL PARDISO.
[[nodiscard]] bool pardiso_available() noexcept;

}  // namespace BladeGen::Turbolenza
