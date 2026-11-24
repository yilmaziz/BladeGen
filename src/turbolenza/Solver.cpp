// Port of v1 solver.py, production paths (see Solver.hpp): solve, solve_density_coupled,
// solve_at_exit_angle, _solve_coupled, _Observers, _observed, _first_grid, _limit, _start_density,
// _context, _mass_flow, _at_angle, _density_pass, _kutta_step, _checked_pass, _record, _attach,
// _solution.
#include "turbolenza/Solver.hpp"

#include "turbolenza/Detail.hpp"
#include "turbolenza/Grid.hpp"
#include "turbolenza/MatrixAssembly.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <typeinfo>
#include <utility>

namespace BladeGen::Turbolenza {

using detail::nan;
using detail::sz;

const char* const kLimitation = "DENSITY-COUPLED / FIXED-DISCHARGE-ANGLE SOLUTION; NOT FINAL AUNGIER KUTTA-CONVERGED SOLUTION";
const char* const kCoupledLimitation =
    "DENSITY + KUTTA COUPLED SOLUTION (Aungier 2003 Eq. 5-34 at the last blade station, exit angle by Aungier 2000 "
    "Eqs. 13-57 and 13-59); VERIFIED ON A CONSTRUCTED CASE ONLY, NOT YET VALIDATED AGAINST EXPERIMENT OR PUBLISHED "
    "RESULTS";
const char* const kMassFlowFromInlet = "derived from the inlet Mach number";
const char* const kMassFlowOverride = "override (no inlet Mach number)";
const char* const kMassFlowChecked = "override, checked against the inlet Mach number";

std::string error_type_name(const std::exception& e)
{
#define TL_NAME(T, name)                                                                                              \
    if (dynamic_cast<const T*>(&e)) return name
    TL_NAME(CancelledError, "Cancelled");
    TL_NAME(TransonicFlowError, "TransonicFlowError");
    TL_NAME(CoupledNotConvergedError, "CoupledNotConvergedError");
    TL_NAME(DensityNotConvergedError, "DensityNotConvergedError");
    TL_NAME(DensityIterationError, "DensityIterationError");
    TL_NAME(ExitRebuildError, "ExitRebuildError");
    TL_NAME(SolverInputError, "SolverInputError");
    TL_NAME(ThermodynamicStateError, "ThermodynamicStateError");
    TL_NAME(LinearSolveError, "LinearSolveError");
    TL_NAME(MatrixAssemblyError, "MatrixAssemblyError");
    TL_NAME(FlowRecoveryError, "FlowRecoveryError");
    TL_NAME(InletStateError, "InletStateError");
    TL_NAME(KuttaError, "KuttaError");
    TL_NAME(DiscretizationError, "DiscretizationError");
    TL_NAME(BoundaryConditionError, "BoundaryConditionError");
    TL_NAME(CoefficientError, "CoefficientError");
    TL_NAME(MetricsError, "MetricsError");
    TL_NAME(CaseError, "CaseError");
    TL_NAME(BoundaryLayerError, "BoundaryLayerError");
    TL_NAME(GridTopologyError, "GridTopologyError");
    TL_NAME(GridGenerationError, "GridGenerationError");
    TL_NAME(GridError, "GridError");
    TL_NAME(GeometryError, "GeometryError");
    TL_NAME(std::invalid_argument, "ValueError");
    TL_NAME(std::logic_error, "logic_error");
    TL_NAME(std::bad_alloc, "MemoryError");
#undef TL_NAME
    return "RuntimeError";
}

namespace {

// v1 _Observers: an observer that throws is switched off; nothing reaches the solve.
class Observers
{
public:
    explicit Observers(RunObserver* o) : o_(o) {}
    void start(const SolveInfo& info) { call([&] { o_->start(info); }); }
    void update(const PassRecord& r) { call([&] { o_->update(r); }); }
    void finish(bool converged, const std::string& error) { call([&] { o_->finish(converged, error); }); }

private:
    template <class F>
    void call(F&& f)
    {
        if (!o_) return;
        try {
            f();
        } catch (...) {
            o_ = nullptr;
        }
    }
    RunObserver* o_;
};

struct Context
{
    std::shared_ptr<const Grid> grid;
    std::shared_ptr<const GridMetrics> metrics;
    std::shared_ptr<const BoundaryConditions> bcs;
    PerfectGas gas;
    BladeRowReference reference;
    double mass_flow_per_depth = 0.0, c_ref = 0.0, relaxation = 0.0;
    int f = 0, l = 0;
    std::optional<double> kutta_relaxation;  // empty: fixed angle, no correction
    std::optional<InletState> inlet;
    std::string mass_flow_source;
    bool stop_at_sonic = true;
    LinearSolverKind linear_solver = LinearSolverKind::SparseLU;
    int threads = 1;
    std::shared_ptr<const Transpiration> transpiration;  // the equivalent inviscid flow (null: v1)
    std::optional<int> bl_iteration;
    double omega = 0.0;  // not in v1: Aungier's omega of a rotating row (the source of Eq. 5-21)
};

struct Pass
{
    Field density_resident;
    CoefficientFields coefficients;
    LinearSolveResult linear;
    FlowField flow;
    StaticState state;
    Field mach;
    Field density_target, density_next;
    std::vector<double> mass_flow_error;
    std::optional<KuttaResidual> kutta;
    std::optional<ExitAngleCorrection> correction;
    double exit_angle = nan(), exit_angle_next = nan();
};

struct Loop
{
    const SolverConfig& config;
    const SolveOptions& options;
    Observers& observers;
    History history;
    std::vector<PassRecord> records;
};

int limit_of(int value, const char* name)
{
    if (value < 1)
        throw SolverInputError(std::format("{} = {}: at least one pass is needed", name, value));
    return value;
}

// v1 _mass_flow: (inlet state, mdot', where mdot' came from).
void mass_flow(Context& c, const Grid& grid, const SolverConfig& config, std::optional<double> given)
{
    if (config.flow.inlet_mach) {
        try {
            c.inlet = inlet_state(config.flow, config.gas, c.c_ref * inlet_passage_width(grid));
        } catch (const InletStateError& error) {
            throw SolverInputError(std::string("the inlet state cannot be derived: ") + error.what());
        }
    }
    if (!given) {
        if (!c.inlet)
            throw SolverInputError("no mass flow: give config.flow.inlet_mach (normal use) or mass_flow_per_depth "
                                   "(verification override)");
        if (!(c.inlet->mass_flow_per_depth > 0.0))
            throw SolverInputError(std::format("the inlet Mach number {} gives mdot' = {}: the solver refuses a zero flow "
                                               "(coefficients.py needs a positive mass flow)",
                                               detail::py_repr(c.inlet->mach),
                                               detail::py_repr(c.inlet->mass_flow_per_depth)));
        c.mass_flow_per_depth = c.inlet->mass_flow_per_depth;
        c.mass_flow_source = kMassFlowFromInlet;
        return;
    }
    const double mdot = *given;
    if (!c.inlet) {
        c.mass_flow_per_depth = mdot;
        c.mass_flow_source = kMassFlowOverride;
        return;
    }
    const double derived = c.inlet->mass_flow_per_depth;
    if (!(std::abs(mdot - derived) <= kMassFlowConsistency * std::abs(derived)))
        throw SolverInputError(std::format("mass_flow_per_depth = {} and the inlet Mach number {} (mdot' = {}) disagree "
                                           "by more than {:g} relative: give one of them, or inlet_mach None for a pure "
                                           "override",
                                           detail::py_repr(mdot), detail::py_repr(c.inlet->mach),
                                           detail::py_repr(derived), kMassFlowConsistency));
    c.mass_flow_per_depth = mdot;
    c.mass_flow_source = kMassFlowChecked;
}

// v1 _context.
Context make_context(std::shared_ptr<const Grid> grid, const SolverConfig& config, std::optional<double> mdot,
                     double c_ref, const SolveOptions& options)
{
    if (!grid) throw SolverInputError("the solver takes a Grid, got None");
    if (!(std::isfinite(c_ref) && c_ref > 0.0))
        throw SolverInputError(std::format("c_ref = {}: the coupled solve needs the true chord as a positive length "
                                           "(metres in SI); the normalized mode of the lower modules is no physical "
                                           "solution",
                                           detail::py_repr(c_ref)));
    const double omega = config.flow.angular_speed;
    if (omega != 0.0) {
        // not in v1: a rotating row lies on a stream sheet with its radius, and its rothalpy is taken
        // at the far-inlet radius (Eqs. 5-51, 13-56)
        if (!grid->boundaries.sheet)
            throw SolverInputError("a rotating row (angular_speed != 0) needs the grid of a stream surface: its radius "
                                   "r(m) carries the rotation (Aungier 2003 Eqs. 5-21, 5-51); a planar cascade has none");
        if (!config.flow.stream_surface_radius || !config.flow.inlet_absolute_tangential_velocity)
            throw SolverInputError("a rotating row needs stream_surface_radius (the far-inlet radius) and "
                                   "inlet_absolute_tangential_velocity for its rothalpy");
        const double r_in = c_ref * grid->boundaries.sheet->radius(grid->m.front());
        if (!(std::abs(*config.flow.stream_surface_radius - r_in) <= 1e-12 * r_in))
            throw SolverInputError(std::format("stream_surface_radius = {} m is not the far-inlet radius of the grid's "
                                               "stream sheet, {} m: the rothalpy must be taken where the inlet state is",
                                               detail::py_repr(*config.flow.stream_surface_radius), detail::py_repr(r_in)));
    }
    if (options.threads < 1 || options.threads > kMaxThreads || options.threads > max_threads())
        throw SolverInputError(std::format("threads = {}: 1 .. {} on this machine (at most {})", options.threads,
                                           max_threads(), kMaxThreads));
    if (options.linear_solver == LinearSolverKind::Pardiso && !pardiso_available())
        throw SolverInputError("the PARDISO linear solver was asked for, but this build has no Intel MKL");
    Context c;
    c.grid = grid;
    c.c_ref = c_ref;
    c.gas = config.gas;
    const TotalState total(config.flow.inlet_absolute_total_pressure, config.flow.inlet_absolute_total_temperature);
    c.reference = blade_row_reference(total, config.gas, config.flow.angular_speed, config.flow.stream_surface_radius,
                                      config.flow.inlet_absolute_tangential_velocity);
    mass_flow(c, *grid, config, mdot);
    c.relaxation = config.iteration.density_relaxation;
    c.f = grid->first_blade_i;
    c.l = grid->last_blade_i;
    c.stop_at_sonic = config.iteration.stop_at_sonic;
    c.linear_solver = options.linear_solver;
    c.threads = options.threads;
    c.transpiration = options.transpiration;
    c.bl_iteration = options.bl_iteration;
    c.metrics = std::make_shared<const GridMetrics>(build_metrics(grid));
    c.bcs = std::make_shared<const BoundaryConditions>(build_boundary_conditions(grid, c.transpiration.get()));
    c.omega = omega;
    if (omega != 0.0) {  // not in v1: U_i = omega r_i at every station (Eq. 5-51)
        c.reference.station_blade_speed.resize(sz(grid->n_m()));
        for (int i = 0; i < grid->n_m(); ++i)
            c.reference.station_blade_speed[sz(i)] = omega * (c_ref * c.metrics->radius_at(i));
    }
    return c;
}

Field start_density(const Grid& grid, const Context& c, const SolveOptions& options)
{
    if (!options.initial_density) return initial_density(c.gas, c.reference, grid.n_m(), grid.n_eta());
    const Field& rho = *options.initial_density;
    if (rho.n_m != grid.n_m() || rho.n_eta != grid.n_eta())
        throw SolverInputError(std::format("the start density has shape ({}, {}), the grid ({}, {})", rho.n_m, rho.n_eta,
                                           grid.n_m(), grid.n_eta()));
    return rho;
}

Grid rebuild(const SolveOptions& o, const Grid& g, double angle)
{
    return o.grid_ops.rebuild_exit ? o.grid_ops.rebuild_exit(g, angle) : rebuild_exit(g, angle);
}
void validate(const SolveOptions& o, const Grid& g)
{
    if (o.grid_ops.validate_topology) o.grid_ops.validate_topology(g);
    else validate_uniform_grid_topology(g);
}

// v1 _at_angle: the grid at angle (start itself at its own angle) with its own metrics and relations.
Context at_angle(const std::shared_ptr<const Grid>& start, double angle, const Context& base, const SolveOptions& o)
{
    std::shared_ptr<const Grid> grid = start;
    if (angle != start->exit_angle()) {
        grid = std::make_shared<const Grid>(rebuild(o, *start, angle));
        const int l = start->last_blade_i;
        bool same = grid->m == start->m && grid->eta == start->eta && grid->first_blade_i == start->first_blade_i &&
                    grid->last_blade_i == start->last_blade_i && grid->node_m.same_shape(start->node_m);
        for (int i = 0; same && i <= l; ++i)
            for (int j = 0; same && j < start->n_eta(); ++j)
                same = grid->node_m(i, j) == start->node_m(i, j) && grid->node_y(i, j) == start->node_y(i, j);
        if (!same)
            throw ExitRebuildError(std::format("the exit rebuild at {:.6f} deg changed the stations, eta, the blade run or "
                                               "an inlet or blade node: the Kutta loop may move the exit nodes only",
                                               detail::degrees(angle)));
        try {
            validate(o, *grid);
        } catch (const std::invalid_argument& error) {
            throw ExitRebuildError(std::format("the exit rebuild at {:.6f} deg left an inadmissible grid: {}",
                                               detail::degrees(angle), error.what()));
        }
    }
    Context c = base;
    c.grid = grid;
    c.metrics = std::make_shared<const GridMetrics>(build_metrics(grid));
    c.bcs = std::make_shared<const BoundaryConditions>(build_boundary_conditions(grid, c.transpiration.get()));
    return c;
}

// v1 _density_pass: steps 3 to 7 with the resident density rho.
Pass density_pass(const Field& rho, const Context& c)
{
    Pass p;
    p.density_resident = rho;
    p.coefficients = build_coefficients(*c.metrics, rho, c.mass_flow_per_depth, c.c_ref, c.omega);
    const AssemblyResult assembly = assemble_from_coefficients(p.coefficients, *c.metrics, *c.bcs);
    p.linear = solve_linear_system(assembly, c.linear_solver, c.threads);
    p.flow = velocity_components(p.linear.psi, rho, p.coefficients, *c.metrics, *c.bcs);
    p.state = static_state_from_velocity(p.flow.W, c.gas, c.reference);
    p.mach = Field(rho.n_m, rho.n_eta);
    const int n = static_cast<int>(rho.size());
    TURBOLENZA_PARALLEL_FOR
    for (int k = 0; k < n; ++k) {
        const std::size_t q = static_cast<std::size_t>(k);
        p.mach.v[q] = p.flow.W.v[q] / p.state.sound_speed.v[q];
    }
    p.density_target = p.state.density;
    p.density_next = relax_density(rho, p.density_target, c.relaxation);
    const std::vector<double> mass = station_mass_flow(p.flow, rho, p.coefficients, *c.metrics);
    p.mass_flow_error.resize(mass.size());
    // against the station's own flow: mdot' (v1), or mdot' (psi_upper - psi_lower) of the equivalent
    // inviscid flow with transpiration (station_flow is exactly 1 without it)
    for (std::size_t k = 0; k < mass.size(); ++k)
        p.mass_flow_error[k] =
            std::abs(mass[k] / (c.mass_flow_per_depth * c.bcs->station_flow(static_cast<int>(k))) - 1.0);
    return p;
}

// v1 _kutta_step (production rule): the residual and, coupled, the corrected next angle.
void kutta_step(Pass& p, const Context& c)
{
    const Grid& grid = *c.grid;
    const double angle = grid.exit_angle();
    p.exit_angle = angle;
    if (!c.kutta_relaxation) {
        try {
            p.kutta = kutta_residual(p.flow, grid);
        } catch (const KuttaError&) {
            p.kutta.reset();
        }
        p.exit_angle_next = angle;
        return;
    }
    p.kutta = kutta_residual(p.flow, grid);
    p.correction = exit_angle_correction(*p.kutta, p.flow, *c.metrics);
    p.exit_angle_next = next_exit_angle(angle, *p.correction, *c.kutta_relaxation);
}

double max_of(const Field& f)
{
    double m = f.v[0];
    for (double v : f.v) m = std::max(m, v);
    return m;
}
double min_of(const Field& f)
{
    double m = f.v[0];
    for (double v : f.v) m = std::min(m, v);
    return m;
}

// v1 _solution.
std::shared_ptr<Solution> make_solution(const Context& c, const Pass& p, const Loop& loop, bool converged, int iterations)
{
    auto s = std::make_shared<Solution>();
    s->grid = c.grid;
    s->metrics = *c.metrics;
    s->psi = p.linear.psi;
    s->flow = p.flow;
    s->history = loop.history;
    s->density = p.density_resident;
    s->density_target = p.density_target;
    s->static_state = p.state;
    s->mach = p.mach;
    s->coefficients = p.coefficients;
    s->linear = p.linear;
    s->passes = loop.records;
    s->converged = converged;
    s->iterations = iterations;
    s->mass_flow_per_depth = c.mass_flow_per_depth;
    s->c_ref = c.c_ref;
    s->density_relaxation = c.relaxation;
    s->exit_angle = c.grid->exit_angle();
    s->kutta = p.kutta;
    s->kutta_residual = p.kutta ? p.kutta->residual : nan();
    s->normalized_kutta_residual = p.kutta ? p.kutta->normalized : nan();
    s->kutta_relaxation = c.kutta_relaxation;
    s->limitation = c.kutta_relaxation ? kCoupledLimitation : kLimitation;
    if (!c.stop_at_sonic) {
        double peak = max_of(p.mach);
        if (!loop.records.empty()) {
            peak = loop.records.front().mach_max;
            for (const PassRecord& r : loop.records) peak = std::max(peak, r.mach_max);
        }
        if (peak >= 1.0)
            s->limitation += std::format("; SONIC STOP OFF: M_rel >= 1 REACHED (peak {:.4f}), OUTSIDE THE SUBSONIC METHOD",
                                         peak);
    }
    s->inlet = c.inlet;
    s->mass_flow_source = c.mass_flow_source;
    return s;
}

// v1 _record: the PassRecord, History's entries from it, then the observers.
void record(const Pass& p, const Pass* previous, int n, bool final, const Context& c, Loop& loop)
{
    PassRecord r;
    r.index = n;
    r.final = final;
    r.bl_iteration = c.bl_iteration.value_or(0);
    r.density_change = change_norm(p.density_next, p.density_resident);
    if (!previous) {
        r.velocity_change = r.psi_change = nan();
    } else {
        r.velocity_change = velocity_change(p.flow.W_m, p.flow.W_theta, previous->flow.W_m, previous->flow.W_theta);
        r.psi_change = change_norm(p.linear.psi, previous->linear.psi);
    }
    r.target_mismatch = change_norm(p.density_target, p.density_resident);
    r.density_min = min_of(p.density_resident);
    r.density_max = max_of(p.density_resident);
    r.speed_min = min_of(p.flow.W);
    r.speed_max = max_of(p.flow.W);
    r.mach_max = max_of(p.mach);
    r.linear_residual = p.linear.residual_inf;
    r.mass_flow_error_max = *std::max_element(p.mass_flow_error.begin(), p.mass_flow_error.end());
    r.mass_flow_error_blade_median = detail::median(std::vector<double>(p.mass_flow_error.begin() + c.f,
                                                                        p.mass_flow_error.begin() + c.l + 1));
    r.exit_angle = p.exit_angle;
    r.kutta_residual = p.kutta ? p.kutta->residual : nan();
    r.kutta_normalized = p.kutta ? p.kutta->normalized : nan();
    r.exit_angle_next = p.exit_angle_next;
    r.tan_correction = p.correction ? p.correction->delta_tan : nan();
    loop.records.push_back(r);
    const bool coupled = c.kutta_relaxation.has_value();
    loop.history.kutta_pass.push_back(coupled ? r.index : 0);
    loop.history.density_change.push_back(r.density_change);
    loop.history.velocity_change.push_back(r.velocity_change);
    loop.history.psi_change.push_back(r.psi_change);
    loop.history.residual.push_back(r.linear_residual);
    if (coupled) {
        loop.history.exit_angle.push_back(r.exit_angle);
        loop.history.kutta_error.push_back(r.kutta_normalized);
    }
    loop.observers.update(r);
    if (loop.options.inspect) {
        try {
            loop.options.inspect(PassFields{loop.records.back(), c.grid, *c.metrics, p.density_resident,
                                            p.density_target, p.density_next, p.coefficients, p.linear, p.flow, p.state,
                                            p.mach, c.mass_flow_per_depth, c.c_ref});
        } catch (...) {  // an inspector never stops or changes the solve
        }
    }
}

// v1 _checked_pass: one pass, recorded; Mach >= 1 stops the loop (sonic stop on).
Pass checked_pass(const Field& rho, const Context& c, Loop& loop, const Pass* previous, int n, bool final)
{
    Pass p = density_pass(rho, c);
    kutta_step(p, c);
    record(p, previous, n, final, c, loop);
    const double peak = max_of(p.mach);
    if (!(peak < 1.0)) {
        if (!c.stop_at_sonic) return p;
        std::size_t at = 0;
        for (std::size_t k = 0; k < p.mach.size(); ++k)
            if (p.mach.v[k] > p.mach.v[at]) at = k;
        const int i = static_cast<int>(at) / p.mach.n_eta, j = static_cast<int>(at) % p.mach.n_eta;
        throw TransonicFlowError(std::format("pass {}: the relative Mach number reaches {:.4f} at node ({}, {}), m = "
                                             "{:.4g}: the subsonic density iteration does not apply and no transonic "
                                             "extension is implemented",
                                             n, peak, i, j, c.grid->m[sz(i)]),
                                 make_solution(c, p, loop, false, n));
    }
    return p;
}

void attach(const SolveOptions& o, const Context* last, const Pass* previous, const Loop& loop, int failed_pass)
{
    if (!o.failure) return;
    try {
        o.failure->solution = (previous && last) ? make_solution(*last, *previous, loop, false, failed_pass) : nullptr;
        o.failure->failed_pass = failed_pass;
    } catch (...) {
    }
}

bool cancelled(const SolveOptions& o) { return o.cancel && o.cancel->load(); }

SolveInfo info_of(const char* mode, const IterationSettings& s, int limit, const Grid& grid, const Context& c,
                  std::optional<double> kutta_tolerance)
{
    SolveInfo info;
    info.mode = mode;
    info.density_tolerance = s.density_tolerance;
    info.velocity_tolerance = velocity_tolerance(s);
    info.kutta_tolerance = kutta_tolerance;
    info.max_passes = limit;
    info.exit_angle_start = grid.exit_angle();
    info.inlet_angle = grid.inlet_angle();
    info.mass_flow_per_depth = c.mass_flow_per_depth;
    info.mass_flow_source = c.mass_flow_source;
    info.n_m = grid.n_m();
    info.n_eta = grid.n_eta();
    info.threads = c.threads;
    info.linear_solver = linear_solver_name(c.linear_solver);
    info.bl_iteration = c.bl_iteration;
    return info;
}

// v1 _observed: start, the loop, finish.
template <class F>
Solution observed(Observers& obs, const SolveInfo& info, F&& loop)
{
    obs.start(info);
    try {
        Solution s = loop();
        obs.finish(s.converged, "");
        return s;
    } catch (const std::exception& error) {
        obs.finish(false, error_type_name(error) + ": " + error.what());
        throw;
    }
}

}  // namespace

Solution solve_at_exit_angle(std::shared_ptr<const Grid> grid, const SolverConfig& config, double c_ref,
                             std::optional<double> mass_flow_per_depth, const SolveOptions& options)
{
    const detail::ThreadScope threads(options.threads);
    const Context context = make_context(grid, config, mass_flow_per_depth, c_ref, options);
    const Field rho0 = start_density(*grid, context, options);
    const IterationSettings& settings = config.iteration;
    const int limit = limit_of(settings.max_density_iterations, "max_density_iterations");
    Observers observers(options.observer);
    Loop loop{config, options, observers, {}, {}};
    return observed(observers, info_of("fixed angle", settings, limit, *grid, context, std::nullopt), [&]() -> Solution {
        std::optional<Pass> previous;
        Field rho = rho0;
        int n = 0;
        try {
            for (n = 0; n < limit; ++n) {
                if (cancelled(options))
                    throw CancelledError(std::format("the run was cancelled after {} pass(es)", n),
                                         previous ? make_solution(context, *previous, loop, false, n) : nullptr);
                Pass current = checked_pass(rho, context, loop, previous ? &*previous : nullptr, n, false);
                previous = std::move(current);
                if (density_converged(loop.history, settings)) {
                    ++n;
                    const Pass final = checked_pass(previous->density_next, context, loop, &*previous, n, true);
                    return *make_solution(context, final, loop, true, n);
                }
                rho = previous->density_next;
            }
        } catch (const DensityIterationError&) {
            throw;
        } catch (const std::exception&) {
            attach(options, &context, previous ? &*previous : nullptr, loop, n);
            throw;
        }
        throw DensityNotConvergedError(
            std::format("the density loop did not converge in {} passes (last density change {:.3e}, velocity change "
                        "{:.3e}, tolerance {:g})",
                        limit, loop.history.density_change.back(), loop.history.velocity_change.back(),
                        settings.density_tolerance),
            make_solution(context, *previous, loop, false, limit));
    });
}

Solution solve_coupled(std::shared_ptr<const Grid> grid, const SolverConfig& config, double c_ref,
                       std::optional<double> mass_flow_per_depth, const SolveOptions& options)
{
    const detail::ThreadScope threads(options.threads);
    Context base = make_context(grid, config, mass_flow_per_depth, c_ref, options);
    const IterationSettings& settings = config.iteration;
    base.kutta_relaxation = settings.kutta_relaxation;
    validate(options, *grid);
    const Field rho0 = start_density(*grid, base, options);
    const int limit = std::min(limit_of(settings.max_density_iterations, "max_density_iterations"),
                               limit_of(settings.max_kutta_iterations, "max_kutta_iterations"));
    Observers observers(options.observer);
    Loop loop{config, options, observers, {}, {}};
    return observed(observers, info_of("coupled", settings, limit, *grid, base, settings.kutta_tolerance),
                    [&]() -> Solution {
        std::optional<Pass> previous;
        std::optional<Context> last;
        Field rho = rho0;
        double beta = grid->exit_angle();
        int n = 0;
        try {
            for (n = 0; n < limit; ++n) {
                if (cancelled(options))
                    throw CancelledError(std::format("the run was cancelled after {} pass(es)", n),
                                         previous ? make_solution(*last, *previous, loop, false, n) : nullptr);
                const Context context = at_angle(grid, beta, base, options);
                Pass current = checked_pass(rho, context, loop, previous ? &*previous : nullptr, n, false);
                previous = std::move(current);
                last = context;
                if (coupled_converged(loop.history, settings)) {
                    ++n;
                    const Context final_context = at_angle(grid, previous->exit_angle_next, base, options);
                    const Pass final = checked_pass(previous->density_next, final_context, loop, &*previous, n, true);
                    return *make_solution(final_context, final, loop, true, n);
                }
                rho = previous->density_next;
                beta = previous->exit_angle_next;
            }
        } catch (const DensityIterationError&) {
            throw;
        } catch (const std::exception&) {
            attach(options, last ? &*last : nullptr, previous ? &*previous : nullptr, loop, n);
            throw;
        }
        throw CoupledNotConvergedError(
            std::format("the coupled density + Kutta loop did not converge in {} passes (last density change {:.3e} "
                        "against {:g}, velocity change {:.3e} against {:g}, eps_K {:.3e} against {:g}; last exit angle "
                        "{:.6f} deg)",
                        limit, loop.history.density_change.back(), settings.density_tolerance,
                        loop.history.velocity_change.back(), velocity_tolerance(settings),
                        loop.history.kutta_error.back(), settings.kutta_tolerance,
                        detail::degrees(previous->exit_angle)),
            make_solution(*last, *previous, loop, false, limit));
    });
}

Solution solve(const Cascade& cascade, const SolverConfig& config, double c_ref, std::optional<double> mass_flow_per_depth,
               const SolveOptions& options)
{
    auto grid = std::make_shared<const Grid>(make_grid(cascade, config.grid, detail::radians(config.flow.inlet_angle_deg),
                                                       detail::radians(config.flow.exit_angle_guess_deg), true));
    return solve_coupled(grid, config, c_ref, mass_flow_per_depth, options);
}

Solution solve_density_coupled(const Cascade& cascade, const SolverConfig& config, double c_ref,
                               std::optional<double> mass_flow_per_depth, const SolveOptions& options)
{
    auto grid = std::make_shared<const Grid>(make_grid(cascade, config.grid, detail::radians(config.flow.inlet_angle_deg),
                                                       detail::radians(config.flow.exit_angle_guess_deg), true));
    return solve_at_exit_angle(grid, config, c_ref, mass_flow_per_depth, options);
}

}  // namespace BladeGen::Turbolenza
