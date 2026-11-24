// Port of v1 controller.py (see Case.hpp): build_case, run_case (timed: wall and process CPU
// time), the console text (format_record of monitor.py, format_case_summary, format_grid_summary,
// format_run_report) and export_run (the v1 CSV / JSON layouts plus fields.csv).  No numerics.
#include "turbolenza/Case.hpp"

#include "turbolenza/Detail.hpp"
#include "turbolenza/Postprocess.hpp"
#include "turbolenza/Solver.hpp"
#include "turbolenza/Viscous.hpp"

#include <nlohmann/json.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <memory>
#include <thread>

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

namespace BladeGen::Turbolenza {

using detail::degrees;
using detail::py_repr;
using detail::radians;
using detail::sz;

int max_threads()
{
    const unsigned hw = std::thread::hardware_concurrency();
    const int n = hw == 0 ? 1 : static_cast<int>(hw);
    return std::max(1, std::min(kMaxThreads, n));
}

const char* bl_coupling_name(BlCoupling c) noexcept { return c == BlCoupling::OneWay ? "one_way" : "two_way"; }
const char* wake_model_name(WakeModel w) noexcept { return w == WakeModel::Constant ? "constant" : "marched"; }
const char* stream_surface_name(StreamSurfaceKind k) noexcept { return k == StreamSurfaceKind::Planar ? "planar" : "layer"; }
const char* row_motion_name(RowMotion r) noexcept { return r == RowMotion::Stationary ? "stationary" : "rotating"; }

double chord_reynolds_number(double mach, double P0, double T0, double chord, double gamma, double R, const SutherlandLaw& law)
{
    if (!(std::isfinite(mach) && mach > 0.0 && std::isfinite(P0) && P0 > 0.0 && std::isfinite(T0) && T0 > 0.0 &&
          std::isfinite(chord) && chord > 0.0 && std::isfinite(gamma) && gamma > 1.0 && std::isfinite(R) && R > 0.0))
        throw CaseError("the chord Reynolds number needs M > 0, P0, T0, chord, R > 0 and gamma > 1, all finite");
    const double ratio = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;  // T0 / T1
    const double T1 = T0 / ratio;
    const double p1 = P0 * std::pow(ratio, -gamma / (gamma - 1.0));
    const double rho1 = p1 / (R * T1);
    const double W1 = mach * std::sqrt(gamma * R * T1);
    double mu = 0.0;
    try {
        mu = sutherland_viscosity(T1, law);
    } catch (const BoundaryLayerError& e) {
        throw CaseError(std::string("the chord Reynolds number: ") + e.what());
    }
    return rho1 * W1 * chord / mu;
}

double total_pressure_for_reynolds(double reynolds, double mach, double T0, double chord, double gamma, double R,
                                   const SutherlandLaw& law)
{
    if (!(std::isfinite(reynolds) && reynolds > 0.0))
        throw CaseError(std::format("the reference Reynolds number {} must be positive", py_repr(reynolds)));
    // Re is proportional to P0 at a fixed Mach number and T0: P0 = Re / Re(P0 = 1 Pa)
    return reynolds / chord_reynolds_number(mach, 1.0, T0, chord, gamma, R, law);
}

bool pardiso_available() noexcept
{
#ifdef TURBOLENZA_HAVE_PARDISO
    return true;
#else
    return false;
#endif
}

namespace {

const char* mode_name(SolverMode m) { return m == SolverMode::CoupledKutta ? "coupled_kutta" : "fixed_exit"; }

std::string opt_repr(const std::optional<double>& v) { return v ? py_repr(*v) : "None"; }
std::string opt_repr(const std::optional<int>& v) { return v ? std::to_string(*v) : "None"; }
std::string opt_repr(const std::optional<std::string>& v) { return v ? "'" + *v + "'" : "None"; }
std::string bool_repr(bool b) { return b ? "True" : "False"; }

// The first (coupled) or fixed (fixed_exit) discharge angle and its source (v1 build_case).
std::pair<double, std::string> exit_angle_of(const CaseSettings& s, double kappa_te_deg)
{
    if (s.solver_mode == SolverMode::FixedExit)
        return s.fixed_exit_angle_deg
                   ? std::pair<double, std::string>{*s.fixed_exit_angle_deg, "fixed: user input"}
                   : std::pair<double, std::string>{kappa_te_deg,
                                                    "fixed: the blade TE angle kappa_te (fixed_exit_angle_deg None)"};
    if (s.initial_exit_angle_mode == ExitAngleStart::BladeTeAngle)
        return {kappa_te_deg, "first Kutta angle: the blade TE angle kappa_te"};
    if (!s.initial_exit_angle_deg) throw CaseError("initial_exit_angle_mode 'manual' needs initial_exit_angle_deg");
    return {*s.initial_exit_angle_deg, "first Kutta angle: user input (manual)"};
}

int pass_limit(const CaseSettings& s)
{
    return s.solver_mode == SolverMode::CoupledKutta ? std::min(s.max_density_iterations, s.max_kutta_iterations)
                                                     : s.max_density_iterations;
}

GridSettings grid_settings(const CaseSettings& s)
{
    GridSettings g;
    g.n_passage = s.n_passage;
    g.n_eta = s.n_eta;
    g.upstream_extension_chords = s.upstream_extension_chords;
    g.downstream_extension_chords = s.downstream_extension_chords;
    g.level = s.grid_level;
    g.edge_nodes = s.edge_nodes;
    return g;
}

FlowConditions flow_conditions(const CaseSettings& s, double exit_angle_deg,
                               const std::optional<VelocityTriangle>& triangle = std::nullopt)
{
    FlowConditions f;
    f.inlet_mach = s.inlet_mach;
    f.inlet_angle_deg = s.inlet_flow_angle_deg;
    f.exit_angle_guess_deg = exit_angle_deg;
    f.inlet_absolute_total_pressure = s.inlet_total_pressure;
    f.inlet_absolute_total_temperature = s.inlet_total_temperature;
    f.angular_speed = s.omega;
    if (triangle) {  // not in v1: the rotating row's rothalpy at the far-inlet radius
        f.stream_surface_radius = triangle->radius;
        f.inlet_absolute_tangential_velocity = triangle->C_theta;
    }
    return f;
}

IterationSettings iteration_settings(const CaseSettings& s)
{
    IterationSettings it;
    it.density_tolerance = s.density_tolerance;
    it.max_density_iterations = s.max_density_iterations;
    it.kutta_tolerance = s.kutta_tolerance;
    it.max_kutta_iterations = s.max_kutta_iterations;
    it.density_relaxation = s.density_relaxation;
    it.velocity_tolerance = s.velocity_tolerance;
    it.kutta_relaxation = s.kutta_relaxation;
    it.stop_at_sonic = s.stop_at_sonic;
    return it;
}

InletData inlet_data(const InletState& i)
{
    InletData d;
    d.mach = i.mach;
    d.angle = i.flow_angle;
    d.speed = i.velocity;
    d.density = i.static_density;
    d.pressure = i.static_pressure;
    d.temperature = i.static_temperature;
    d.sound_speed = i.sound_speed;
    d.mass_flow_per_depth = i.mass_flow_per_depth;
    d.total_pressure = i.total_pressure;
    d.dynamic_pressure = i.dynamic_pressure;
    return d;
}

double process_cpu_seconds()
{
#ifdef _WIN32
    FILETIME creation, exit, kernel, user;
    if (!GetProcessTimes(GetCurrentProcess(), &creation, &exit, &kernel, &user)) return 0.0;
    auto seconds = [](const FILETIME& t) {
        const unsigned long long ticks =
            (static_cast<unsigned long long>(t.dwHighDateTime) << 32) | static_cast<unsigned long long>(t.dwLowDateTime);
        return static_cast<double>(ticks) * 1e-7;
    };
    return seconds(kernel) + seconds(user);
#else
    timespec ts{};
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return static_cast<double>(ts.tv_sec) + 1e-9 * static_cast<double>(ts.tv_nsec);
#endif
}

std::string iso_now()
{
    const std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return buf;
}

// The run's observer: keeps every record and the time spent in the user's observer; an
// observer that throws is switched off (it never reaches the solve).
class Recorder final : public RunObserver
{
public:
    explicit Recorder(RunObserver* user) : user_(user) {}
    void start(const SolveInfo& info) override { call([&] { user_->start(info); }); }
    void update(const PassRecord& r) override
    {
        records.push_back(r);
        call([&] { user_->update(r); });
    }
    void finish(bool converged, const std::string& error) override { call([&] { user_->finish(converged, error); }); }
    std::vector<PassRecord> records;

private:
    template <class F>
    void call(F&& f)
    {
        if (!user_) return;
        try {
            f();
        } catch (...) {
            user_ = nullptr;
        }
    }
    RunObserver* user_;
};

SurfaceData surface_data(const std::vector<double>& x_c, const std::vector<double>& m_c, const std::vector<double>& mach,
                         const std::vector<double>& pressure, const std::optional<InletState>& inlet)
{
    SurfaceData d{x_c, m_c, mach, pressure, {}, {}};
    if (inlet) {
        d.dunavant_S = pressure_coefficient_values(pressure, inlet, PressureConvention::Dunavant);
        d.cp_static = pressure_coefficient_values(pressure, inlet, PressureConvention::InletStatic);
    }
    return d;
}

std::string rule(const std::string& title = "")
{
    if (title.empty()) return std::string(64, '-');
    std::string s = "--- " + title + " ";
    if (s.size() < 64) s += std::string(64 - s.size(), '-');
    return s;
}

ViscousSettings viscous_settings(const CaseSettings& s)
{
    ViscousSettings v;
    v.coupling = s.bl_coupling;
    v.wake = s.wake_model;
    v.layer.laminar_cf_form = s.laminar_cf_form;
    v.layer.transition_re_theta = s.transition_re_theta;
    v.layer.separation_hk = s.separation_hk;
    v.viscosity = s.sutherland;
    v.relaxation = s.bl_relaxation;
    v.tolerance = s.bl_tolerance;
    v.max_iterations = s.max_bl_iterations;
    return v;
}

BoundaryLayerSurfaceData layer_data(const WallLayer& w)
{
    const BoundaryLayerProfile& p = w.profile;
    BoundaryLayerSurfaceData d;
    d.x_c = w.x_c;
    d.x = p.x;
    d.W_e = p.W;
    d.rho_e = p.rho;
    d.T_e = p.T;
    d.mu = p.mu;
    d.theta = p.theta;
    d.delta_star = p.delta_star;
    d.delta = p.delta;
    d.H = p.H;
    d.H_k = p.Hk;
    d.c_f = p.cf;
    d.Lambda = p.Lambda;
    d.Re_theta = p.Re_theta;
    d.mass_deficit = p.mass_deficit;
    for (BlRegime g : p.regime) d.regime.emplace_back(bl_regime_name(g));
    auto at = [&](int k) { return k < 0 ? detail::nan() : w.x_c[sz(k)]; };
    d.transition_x_c = at(p.transition_index);
    d.separation_x_c = at(p.separation_index);
    d.laminar_separation_x_c = at(p.laminar_separation_index);
    return d;
}

ViscousData viscous_data(const ViscousSolve& v, const CaseSettings& s)
{
    ViscousData d;
    d.coupling = bl_coupling_name(s.bl_coupling);
    d.laminar_cf_form = laminar_cf_form_name(s.laminar_cf_form);
    d.wake_model = wake_model_name(s.wake_model);
    d.converged = v.converged;
    d.history = v.history;
    d.loss_omega = d.loss_Theta = d.loss_Delta_star = d.te_jump_mismatch = d.inviscid_exit_angle = detail::nan();
    const std::shared_ptr<const Solution> owner = v.analysis_solution ? v.analysis_solution : v.solution;
    if (v.analysis && owner) {
        const ViscousAnalysis& a = *v.analysis;
        const Grid& g = *owner->grid;
        d.layers_iteration = v.analysis_iteration;
        d.layers_of_last_solve = owner == v.solution;
        d.suction = layer_data(a.suction());
        d.pressure_side = layer_data(a.pressure());
        const int l = g.last_blade_i;
        for (std::size_t q = 1; q < a.wake.stations.size(); ++q) {
            const int i = a.wake.stations[q];
            d.wake_m_c.push_back(a.wake.m_c[q]);
            d.wake_mass_deficit.push_back(a.deficit.wake[sz(i)]);
            if (a.wake.marched) {
                d.wake_theta.push_back(a.wake.profile.theta[q]);
                d.wake_delta_star.push_back(a.wake.profile.delta_star[q]);
                d.wake_H.push_back(a.wake.profile.H[q]);
            }
        }
        d.loss_omega = a.loss.omega;
        d.loss_Theta = a.loss.Theta;
        d.loss_Delta_star = a.loss.Delta_star;
        if (!v.used.wake.empty() && l + 1 < g.n_m())
            d.te_jump_mismatch = std::abs(v.used.wake[sz(l + 1)] - (v.used.lower[sz(l)] + v.used.upper[sz(l)])) /
                                 owner->mass_flow_per_depth;
    }
    if (v.inviscid) {
        d.inviscid_exit_angle = v.inviscid->exit_angle;
        try {
            const BladeSurfaces b = blade_surfaces(*v.inviscid);
            d.inviscid_suction = surface_data(b.x_c_ss, b.m_c, b.mach_ss, b.pressure_ss, v.inviscid->inlet);
            d.inviscid_pressure_side = surface_data(b.x_c_ps, b.m_c, b.mach_ps, b.pressure_ps, v.inviscid->inlet);
        } catch (const std::exception&) {
        }
    }
    return d;
}

// Not in v1: the absolute frame of a rotating row, C = W + U with U = omega r_i at station i.
AbsoluteFrameData absolute_frame(const SolutionData& d, const PerfectGas& gas)
{
    const int n_m = d.grid.n_m(), n_eta = d.grid.n_eta();
    AbsoluteFrameData a{Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta),
                        Field(n_m, n_eta)};
    const double cp = gas.cp(), exponent = gas.gamma / (gas.gamma - 1.0);
    for (int i = 0; i < n_m; ++i) {
        const double U = d.omega * d.station_radius[sz(i)];
        for (int j = 0; j < n_eta; ++j) {
            const double ct = d.W_theta(i, j) + U, cm = d.W_m(i, j), T = d.temperature(i, j);
            const double C = std::hypot(cm, ct);
            const double T0 = T + 0.5 * C * C / cp;
            a.C_theta(i, j) = ct;
            a.C(i, j) = C;
            a.alpha(i, j) = std::atan2(ct, cm);
            a.mach(i, j) = C / std::sqrt(gas.gamma * gas.gas_constant * T);
            a.total_temperature(i, j) = T0;
            a.total_pressure(i, j) = d.pressure(i, j) * std::pow(T0 / T, exponent);
        }
    }
    return a;
}

// Not in v1: the far exit of a rotating row, mass averages over eta (weights rho W_m; ours).
RotorSummary rotor_summary(const SolutionData& d, const VelocityTriangle& in, const PerfectGas& gas)
{
    const AbsoluteFrameData& a = *d.absolute;
    const int i = d.grid.n_m() - 1, n_eta = d.grid.n_eta();
    std::vector<double> w(sz(n_eta)), x(sz(n_eta));
    for (int j = 0; j < n_eta; ++j) w[sz(j)] = d.density(i, j) * d.W_m(i, j);
    const double total = detail::trapezoid(w.data(), d.grid.eta);
    auto mean = [&](const Field& f) {
        for (int j = 0; j < n_eta; ++j) x[sz(j)] = w[sz(j)] * f(i, j);
        return detail::trapezoid(x.data(), d.grid.eta) / total;
    };
    RotorSummary r;
    r.rpm = in.rpm;
    r.omega = in.omega;
    r.U_in = in.U;
    r.U_out = d.omega * d.station_radius[sz(i)];
    r.C_m_out = mean(d.W_m);
    r.C_theta_out = mean(a.C_theta);
    r.alpha_out_deg = degrees(std::atan2(r.C_theta_out, r.C_m_out));
    r.mach_out = mean(a.mach);
    r.total_pressure_out = mean(a.total_pressure);
    r.total_temperature_out = mean(a.total_temperature);
    r.total_pressure_ratio = r.total_pressure_out / in.total_pressure;
    r.total_temperature_ratio = r.total_temperature_out / in.total_temperature;
    r.euler_work = r.U_out * r.C_theta_out - r.U_in * in.C_theta;
    r.total_enthalpy_rise = gas.cp() * (r.total_temperature_out - in.total_temperature);
    return r;
}

}  // namespace

// --- case construction -----------------------------------------------------------------------

namespace {

// v1 build_case, for a planar cascade (stream == nullptr) or a blade on its stream surface (not in v1).
Case build_case_on(const CaseSettings& settings, const BladeSection& section, const StreamSurfaceBlade* stream)
{
    const CaseSettings& s = settings;
    if (s.threads < 1 || s.threads > max_threads())
        throw CaseError(std::format("threads {}: 1 .. {} (the machine offers {} logical CPUs, the limit is {})", s.threads,
                                    max_threads(), std::thread::hardware_concurrency(), kMaxThreads));
    if (s.linear_solver == LinearSolverKind::Pardiso && !pardiso_available())
        throw CaseError("linear solver PARDISO: this build has no Intel MKL (pardiso_available() is false)");
    // not in v1: the stream surface and the row
    if ((s.stream_surface == StreamSurfaceKind::Layer) != (stream != nullptr))
        throw CaseError(s.stream_surface == StreamSurfaceKind::Layer
                            ? "stream surface 'layer' needs the blade on its stream surface (a StreamSurfaceBlade)"
                            : "a StreamSurfaceBlade was given with stream surface 'planar'");
    const bool rotating = s.row == RowMotion::Rotating;
    if (rotating && !stream)
        throw CaseError("a rotating row needs the BladeGen layer's stream surface (stream surface 'layer'): its radius "
                        "r(m) carries the rotation (Aungier 2003 Eqs. 5-21, 5-51); a planar cascade has none");
    if (!rotating && s.omega != 0.0)
        throw CaseError(std::format("omega = {} rad/s with a stationary row: set the row to 'rotating' (rpm and the "
                                    "absolute inlet state)",
                                    py_repr(s.omega)));
    if (rotating) {
        if (!std::isfinite(s.rpm)) throw CaseError("rpm must be finite");
        if (!(std::isfinite(s.inlet_absolute_mach) && s.inlet_absolute_mach > 0.0))
            throw CaseError(std::format("the absolute inlet Mach number {} must be positive", py_repr(s.inlet_absolute_mach)));
        if (!(std::isfinite(s.inlet_absolute_flow_angle_deg) && std::abs(s.inlet_absolute_flow_angle_deg) < 90.0))
            throw CaseError(std::format("the absolute inlet angle {} deg must lie in (-90, 90)",
                                        py_repr(s.inlet_absolute_flow_angle_deg)));
    }
    const double chord = stream ? stream->chord : s.chord;
    if (!(std::isfinite(chord) && chord > 0.0))
        throw CaseError(std::format("chord {}: a positive length in metres", py_repr(chord)));
    if (s.viscous) {
        if (!(s.bl_relaxation > 0.0 && s.bl_relaxation <= 1.0))
            throw CaseError(std::format("bl_relaxation {}: 0 < omega_BL <= 1", py_repr(s.bl_relaxation)));
        if (!(std::isfinite(s.bl_tolerance) && s.bl_tolerance > 0.0))
            throw CaseError(std::format("bl_tolerance {}: a positive number", py_repr(s.bl_tolerance)));
        if (s.max_bl_iterations < 0)
            throw CaseError(std::format("max_bl_iterations {}: 0 or more", s.max_bl_iterations));
        if (!(std::isfinite(s.transition_re_theta) && s.transition_re_theta > 0.0))
            throw CaseError(std::format("transition_re_theta {}: a positive number", py_repr(s.transition_re_theta)));
        if (!(std::isfinite(s.separation_hk) && s.separation_hk > 1.0))
            throw CaseError(std::format("separation_hk {}: a number > 1", py_repr(s.separation_hk)));
        try {
            (void)sutherland_viscosity(s.inlet_total_temperature, s.sutherland);
        } catch (const BoundaryLayerError& e) {
            throw CaseError(std::string("the viscosity law: ") + e.what());
        }
        if (s.edge_nodes == EdgeNodes::Touching)
            throw CaseError("the boundary layer needs the first and last blade nodes inside the passage (edge nodes "
                            "'inside'): the touching experiment puts them on the nose tip and the tail end");
        if (s.omega != 0.0 || rotating)
            throw CaseError("the boundary layer is set up for a stationary cascade (omega = 0) only");
        // not in v1 (ours): the march of BoundaryLayer.cpp has b = 1 and no r; Aungier's Eq. 5-106,
        // (1/b) d(b rho_e W_e^2 theta)/dx + ..., carries b(x).  On a cylinder (r and b constant) the
        // sheet is the planar cascade and the layer applies unchanged; otherwise it is refused
        if (stream && !stream->sheet->cylinder())
            throw CaseError("the boundary layer (Aungier 2003 Sec. 5.6) is implemented with b = 1 and no radius: its "
                            "momentum integral Eq. 5-106 carries b(x), which this march does not have. On a stream "
                            "sheet with r(m) or b(m) it is refused; it runs on the planar cascade and on a cylinder");
    }
    Case c;
    c.settings = s;
    if (stream) {  // the chord, solidity and stagger of the blade in the (m, r_LE theta) plane
        c.settings.chord = stream->chord;
        c.settings.solidity = stream->solidity;
        c.settings.stagger_deg = stream->stagger_deg;
        c.stream = *stream;
    }
    c.section = section;
    try {
        c.cascade = stream ? make_cascade(*stream) : make_cascade(section, s.solidity, s.stagger_deg);
    } catch (const std::invalid_argument& error) {
        throw CaseError(std::string("the cascade cannot be built: ") + error.what());
    }
    PerfectGas gas;
    try {
        gas = PerfectGas(s.gamma, s.gas_constant);
    } catch (const ThermodynamicStateError& error) {
        throw CaseError(std::string("the gas cannot be formed: ") + error.what());
    }
    // a rotating row: the relative inlet from the absolute one at the far-inlet radius (velocity_triangle).
    // The far inlet lies upstream_extension_chords ahead of the LE (a standard level rounds it to whole
    // cells): the triangle is taken there first and again at the grid's own m_in (below)
    auto triangle_at = [&](double m_in) {
        try {
            const double r_in = chord * c.cascade.sheet->radius(m_in);
            VelocityTriangle t = velocity_triangle(s.rpm, s.inlet_absolute_flow_angle_deg, s.inlet_absolute_mach,
                                                   s.inlet_total_pressure, s.inlet_total_temperature, gas, r_in);
            if (!(std::abs(t.beta) < radians(89.0)))
                throw CaseError(std::format("the relative inlet angle beta1 = {:.4f} deg: the side lines cannot lie along "
                                            "it",
                                            degrees(t.beta)));
            c.triangle = t;
            c.settings.omega = t.omega;
            c.settings.inlet_mach = t.mach_relative;
            c.settings.inlet_flow_angle_deg = degrees(t.beta);
        } catch (const std::invalid_argument& error) {
            throw CaseError(std::string("the inlet velocity triangle: ") + error.what());
        }
    };
    if (rotating) triangle_at(c.cascade.active.leading_edge_point.m - s.upstream_extension_chords * c.cascade.active.chord);
    const double kappa_te = degrees(c.cascade.kappa_te);
    const auto [angle, source] = exit_angle_of(s, kappa_te);
    (void)source;
    if (s.grid_level && (s.n_passage || s.n_eta))
        throw CaseError(std::format("grid_level '{}' sets the resolution: give n_passage = n_eta = None with it (got {}, "
                                    "{}), or grid_level None for the explicit n_passage x n_eta grid",
                                    *s.grid_level, opt_repr(s.n_passage), opt_repr(s.n_eta)));
    const GridSettings gs = grid_settings(s);
    c.exit_angle_start = radians(angle);
    auto build_grid_now = [&]() {
        c.inlet_angle = radians(c.settings.inlet_flow_angle_deg);
        try {
            if (!s.grid_level) {
                c.grid = make_grid(c.cascade, gs, c.inlet_angle, c.exit_angle_start, true);
            } else {  // the same grid make_grid(level) gives the solver, with its metadata
                StandardGrid standard = build_grid(c.cascade, *s.grid_level, c.inlet_angle, c.exit_angle_start,
                                                   {s.upstream_extension_chords, s.downstream_extension_chords}, s.edge_nodes);
                // v1 _standard_checks' last check, "metrics with positive S", needs the metrics, which come
                // after the grid layer in v2: run here, with v1's message
                try {
                    const GridMetrics metrics = build_metrics(std::make_shared<const Grid>(standard.grid));
                    for (double S : metrics.S)
                        if (!(S > 0.0)) throw GridError("metrics with positive S: S is not positive");
                } catch (const MetricsError& error) {
                    throw GridError(std::string("metrics: MetricsError: ") + error.what());
                }
                standard.info.checks.emplace_back("metrics with positive S");
                c.grid = std::move(standard.grid);
                c.standard_info = std::move(standard.info);
            }
        } catch (const std::invalid_argument& error) {
            const std::string which = s.grid_level ? "grid_level " + *s.grid_level
                                                   : std::format("n_passage {}, n_eta {}", opt_repr(s.n_passage),
                                                                 opt_repr(s.n_eta));
            throw CaseError(std::format("the grid cannot be used ({}): {}", which, error.what()));
        }
    };
    build_grid_now();
    if (rotating && c.grid.m.front() != c.cascade.active.leading_edge_point.m - s.upstream_extension_chords) {
        // the grid's own far inlet (a standard level's whole cells): the triangle there, the grid again at its
        // beta1 (only the side lines' angles change, not the stations)
        const double m_in = c.grid.m.front();
        triangle_at(m_in);
        build_grid_now();
        if (c.grid.m.front() != m_in)
            throw CaseError("the grid's far inlet moved with the relative inlet angle: the triangle cannot be closed");
    }
    const FlowConditions conditions = flow_conditions(c.settings, angle, c.triangle);
    const double pitch = chord * inlet_passage_width(c.grid);
    std::optional<InletState> inlet;
    if (c.settings.inlet_mach) {
        try {
            inlet = inlet_state(conditions, gas, pitch);
        } catch (const std::invalid_argument& error) {
            throw CaseError(std::string("the inlet state cannot be derived: ") + error.what());
        }
        c.inlet = inlet_data(*inlet);
    }
    const auto& override_ = s.mass_flow_per_depth_override;
    if (!override_) {
        if (!inlet) throw CaseError("no mass flow: give inlet_mach (normal use) or mass_flow_per_depth_override");
        c.mass_flow_per_depth = inlet->mass_flow_per_depth;
        c.mass_flow_source = kMassFlowFromInlet;
    } else if (!inlet) {
        c.mass_flow_per_depth = *override_;
        c.mass_flow_source = kMassFlowOverride;
    } else {
        if (!(std::abs(*override_ - inlet->mass_flow_per_depth) <= kMassFlowConsistency * std::abs(inlet->mass_flow_per_depth)))
            throw CaseError(std::format("mass_flow_per_depth_override {} contradicts the inlet Mach number {} (mdot' {}); "
                                        "set inlet_mach None for a pure override (solver.MASS_FLOW_CONSISTENCY = {:g})",
                                        py_repr(*override_), opt_repr(c.settings.inlet_mach),
                                        py_repr(inlet->mass_flow_per_depth), kMassFlowConsistency));
        c.mass_flow_per_depth = *override_;
        c.mass_flow_source = kMassFlowChecked;
    }
    return c;
}

}  // namespace

Case build_case(const CaseSettings& settings, const BladeSection& section)
{
    return build_case_on(settings, section, nullptr);
}

Case build_case(const CaseSettings& settings, const StreamSurfaceBlade& blade)
{
    return build_case_on(settings, blade.section, &blade);
}

// --- the timed solve -------------------------------------------------------------------------

RunResult run_case(const Case& c, RunObserver* observer, const std::atomic<bool>* cancel)
{
    const CaseSettings& s = c.settings;
    SolverConfig config;
    config.grid = grid_settings(s);
    config.flow = flow_conditions(s, degrees(c.exit_angle_start), c.triangle);
    config.iteration = iteration_settings(s);
    RunResult r;
    r.case_name = s.case_name;
    r.mode = mode_name(s.solver_mode);
    r.stop_at_sonic = s.stop_at_sonic;
    Recorder recorder(observer);
    SolveFailure failure;
    SolveOptions options;
    options.observer = &recorder;
    options.cancel = cancel;
    options.threads = s.threads;
    options.linear_solver = s.linear_solver;
    options.failure = &failure;
    std::shared_ptr<const Solution> solution;
    std::optional<ViscousSolve> viscous;
    r.started = iso_now();
    const auto wall0 = std::chrono::steady_clock::now();
    const double cpu0 = process_cpu_seconds();
    try {
        config.gas = PerfectGas(s.gamma, s.gas_constant);
        auto grid = std::make_shared<const Grid>(c.grid);
        const std::optional<double> override_ = s.mass_flow_per_depth_override;
        if (!s.viscous) {
            solution = std::make_shared<const Solution>(s.solver_mode == SolverMode::CoupledKutta
                                                            ? solve_coupled(grid, config, s.chord, override_, options)
                                                            : solve_at_exit_angle(grid, config, s.chord, override_, options));
        } else {
            options.failure = nullptr;  // solve_viscous keeps its own per solve
            viscous = solve_viscous(grid, config, s.chord, override_, options, viscous_settings(s),
                                    s.solver_mode == SolverMode::CoupledKutta);
            solution = viscous->solution;
            r.error = viscous->error;
            if (viscous->failed_pass) r.failed_pass = viscous->failed_pass;
        }
    } catch (const DensityIterationError& error) {
        solution = error.solution;
        r.error = error_type_name(error) + ": " + error.what();
    } catch (const std::exception& error) {
        solution = failure.solution;
        if (failure.failed_pass >= 0) r.failed_pass = failure.failed_pass;
        r.error = error_type_name(error) + ": " + error.what();
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    const double cpu = process_cpu_seconds() - cpu0;
    if (s.measure_runtime) {
        r.runtime.wall_time_s = wall;
        r.runtime.process_cpu_time_s = cpu;
        r.runtime.effective_cpu_cores = wall > 0.0 ? cpu / wall : detail::nan();
    }
    r.runtime.logical_cpus_available = static_cast<int>(std::thread::hardware_concurrency());
    r.runtime.threads = s.threads;
    r.runtime.linear_solver = linear_solver_name(s.linear_solver);
    r.records = std::move(recorder.records);
    r.converged = r.error.empty() && solution && solution->converged;
    r.status = r.converged ? "converged" : "failed: " + r.error.substr(0, r.error.find(':'));
    if (solution) {
        const Solution& sol = *solution;
        r.iterations = sol.iterations;
        SolutionData d;
        d.grid = *sol.grid;
        d.c_ref = sol.c_ref;
        d.psi = sol.psi;
        d.density = sol.density;
        d.pressure = sol.static_state.pressure;
        d.temperature = sol.static_state.temperature;
        d.mach = sol.mach;
        d.W = sol.flow.W;
        d.W_m = sol.flow.W_m;
        d.W_theta = sol.flow.W_theta;
        if (sol.inlet) d.inlet = inlet_data(*sol.inlet);
        d.exit_angle = sol.exit_angle;
        d.kutta_residual = sol.kutta_residual;
        d.normalized_kutta_residual = sol.normalized_kutta_residual;
        d.mass_flow_per_depth = sol.mass_flow_per_depth;
        d.limitation = sol.limitation;
        if (sol.metrics.on_sheet)  // not in v1: the stream sheet at the stations
            for (int i = 0; i < sol.grid->n_m(); ++i) {
                d.station_radius.push_back(sol.c_ref * sol.metrics.radius_at(i));
                d.station_thickness.push_back(sol.metrics.thickness_at(i));
            }
        d.omega = sol.coefficients.omega;
        if (d.omega != 0.0 && !d.station_radius.empty()) d.absolute = absolute_frame(d, config.gas);
        try {
            const BladeSurfaces b = blade_surfaces(sol);
            d.suction = surface_data(b.x_c_ss, b.m_c, b.mach_ss, b.pressure_ss, sol.inlet);
            d.pressure_side = surface_data(b.x_c_ps, b.m_c, b.mach_ps, b.pressure_ps, sol.inlet);
            const SolutionSummary p = summarize(sol);
            RunSummary sum;
            sum.limitation = p.limitation;
            sum.inlet_mach = p.inlet_mach;
            sum.inlet_angle_deg = p.inlet_angle_deg;
            sum.discharge_angle_deg = p.discharge_angle_deg;
            sum.exit_velocity_angle_deg = p.exit_velocity_angle_deg;
            sum.exit_angle_mismatch_deg = p.exit_angle_mismatch_deg;
            sum.turning_deg = p.turning_deg;
            sum.surface_mach_max = p.surface_mach_max;
            sum.surface_mach_min = p.surface_mach_min;
            sum.ss_peak_mach = p.ss_peak_mach;
            sum.ss_peak_x_c = p.ss_peak_x_c;
            sum.ps_peak_mach = p.ps_peak_mach;
            sum.ps_peak_x_c = p.ps_peak_x_c;
            sum.te_ss_mach = b.mach_ss.back();
            sum.te_ps_mach = b.mach_ps.back();
            sum.kutta_residual = p.kutta_residual;
            sum.normalized_kutta_residual = p.normalized_kutta_residual;
            sum.mach_max = p.mach_max;
            sum.mach_max_i = p.mach_max_node[0];
            sum.mach_max_j = p.mach_max_node[1];
            sum.density_min = p.density_min;
            sum.density_max = p.density_max;
            sum.pressure_min = p.pressure_min;
            sum.pressure_max = p.pressure_max;
            sum.temperature_min = p.temperature_min;
            sum.temperature_max = p.temperature_max;
            sum.mass_flow_per_depth = p.mass_flow_per_depth;
            sum.mass_flow_source = p.mass_flow_source;
            if (d.absolute && c.triangle) sum.rotor = rotor_summary(d, *c.triangle, config.gas);
            r.summary = sum;
        } catch (const std::exception& error) {
            r.postprocess_error = error_type_name(error) + ": " + error.what();
        }
        r.solution = std::move(d);
    }
    if (c.inlet) {
        try {
            r.chord_reynolds = c.inlet->density * c.inlet->speed * s.chord / sutherland_viscosity(c.inlet->temperature, s.sutherland);
        } catch (const BoundaryLayerError&) {
        }
    }
    if (viscous) r.viscous = viscous_data(*viscous, s);
    return r;
}

// --- console text -----------------------------------------------------------------------------

std::string format_record(const PassRecord& r)
{
    const std::string label = r.final ? " final" : std::format("{:6d}", r.index);
    return std::format("iter{} | rho {:.1e} | W {:.1e} | Kutta {:.1e} | beta2 {:.4f} deg | Mmax {:.3f} | lin {:.1e}", label,
                       r.density_change, r.velocity_change, r.kutta_normalized, degrees(r.exit_angle), r.mach_max,
                       r.linear_residual);
}

std::string format_case_summary(const Case& c)
{
    const CaseSettings& s = c.settings;
    std::vector<std::string> lines{rule("CASE INPUTS")};
    auto row = [&](const std::string& label, const std::string& value, const std::string& tag = "[user]",
                   const std::string& note = "") {
        lines.push_back(std::format("  {:34s}: {:>24s}  {}{}", label, value, tag, note.empty() ? "" : "  " + note));
    };
    auto tag = [](bool is_default) { return is_default ? "[user = project default]" : "[user]"; };
    auto head = [&](const std::string& text) { lines.push_back(" " + text); };
    head("0. run control");
    row("case name", s.case_name);
    row("solver mode", mode_name(s.solver_mode));
    if (s.solver_mode == SolverMode::FixedExit) row("fixed exit angle [deg]", opt_repr(s.fixed_exit_angle_deg));
    head("1. gas / thermodynamics");
    row("gamma", py_repr(s.gamma));
    row("gas constant R [J/(kg K)]", py_repr(s.gas_constant));
    row("Pt,in absolute [Pa]", py_repr(s.inlet_total_pressure));
    row("Tt,in absolute [K]", py_repr(s.inlet_total_temperature));
    const bool rotating = s.row == RowMotion::Rotating;
    row("omega [rad/s]", py_repr(s.omega), rotating ? "[derived]" : "[user]",
        rotating ? "-2 pi rpm / 60 (positive toward +theta)" : "0: stationary row");
    head("2. blade / cascade geometry");
    row("blade", c.section.designation, "[BladeGen]");
    row("chi_le, chi_te [deg]", std::format("{:.4f}, {:.4f}", degrees(c.section.chi_le), degrees(c.section.chi_te)),
        "[BladeGen]", c.stream ? "from the (m, r_LE theta) chord line" : "");
    const char* geo = c.stream ? "[BladeGen]" : "[user]";
    row("chord c [m]", py_repr(s.chord), geo, c.stream ? "LE-TE in the (m, r_LE theta) plane" : "");
    row("solidity sigma", py_repr(s.solidity), geo, c.stream ? "chord / (r_LE 2 pi / N)" : "");
    row("stagger [deg]", py_repr(s.stagger_deg), geo, c.stream ? "chord line in the (m, r_LE theta) plane" : "");
    row("pitch [m]", std::format("{:.9g}", s.chord * c.grid.pitch()), "[derived]",
        c.stream ? "r_LE 2 pi / N (the y-pitch; S = r Theta varies with r)" : "chord / solidity");
    row("angle of attack [deg]", std::format("{:.6g}", s.inlet_flow_angle_deg - s.stagger_deg), "[derived]",
        "beta1 - stagger");
    row("kappa_le, kappa_te [deg]", std::format("{:.4f}, {:.4f}", degrees(c.cascade.kappa_le), degrees(c.cascade.kappa_te)),
        "[derived]");
    head("3. inlet flow");
    row("inlet Mach number M1", opt_repr(s.inlet_mach));
    row("inlet flow angle beta1 [deg]", py_repr(s.inlet_flow_angle_deg), "[user]", "from the axial direction");
    row("mass-flow override [kg/(m s)]", opt_repr(s.mass_flow_per_depth_override));
    row("mass flow per depth [kg/(m s)]", std::format("{:.9g}", c.mass_flow_per_depth), "[derived]", c.mass_flow_source);
    head("4. grid");
    row("grid level", opt_repr(s.grid_level), "[user]",
        s.grid_level ? "geometry-driven standard grid (grid.py section 9)" : "None: explicit n_passage x n_eta");
    row("n_passage", opt_repr(s.n_passage));
    row("n_eta", opt_repr(s.n_eta));
    row("upstream extension [chords]", py_repr(s.upstream_extension_chords), tag(s.upstream_extension_chords == 1.0));
    row("downstream extension [chords]", py_repr(s.downstream_extension_chords), tag(s.downstream_extension_chords == 1.0));
    head("5. density / nonlinear");
    row("density tolerance", py_repr(s.density_tolerance));
    row("velocity tolerance", opt_repr(s.velocity_tolerance), tag(!s.velocity_tolerance), "None: the density tolerance");
    row("max density iterations", std::to_string(s.max_density_iterations));
    row("density relaxation lambda_rho", py_repr(s.density_relaxation), tag(s.density_relaxation == kDensityRelaxation),
        "0.5 ours; 1.0 = Aungier direct replacement");
    row("stop at M_rel >= 1", bool_repr(s.stop_at_sonic), tag(s.stop_at_sonic),
        s.stop_at_sonic ? "" : "OFF: the run goes on past sonic, outside the subsonic method");
    head("6. Kutta");
    row("Kutta tolerance (eps_K)", py_repr(s.kutta_tolerance));
    row("max Kutta iterations", std::to_string(s.max_kutta_iterations));
    row("Kutta relaxation lambda_K", py_repr(s.kutta_relaxation), tag(s.kutta_relaxation == kKuttaRelaxation),
        "1.0 ours: the undamped Aungier correction");
    row("initial exit-angle mode", s.initial_exit_angle_mode == ExitAngleStart::BladeTeAngle ? "blade_te_angle" : "manual");
    row("initial exit angle [deg]", opt_repr(s.initial_exit_angle_deg));
    const auto [angle, source] = exit_angle_of(s, degrees(c.cascade.kappa_te));
    row("exit angle used [deg]", std::format("{:.6f}", angle), "[derived]", source);
    row("pass limit", std::to_string(pass_limit(s)), "[derived]",
        s.solver_mode == SolverMode::CoupledKutta ? "min(max density, max Kutta)" : "max density iterations");
    head("7. compute / runtime");
    row("threads", std::to_string(s.threads), "[user]", std::format("1 .. {} on this machine", max_threads()));
    row("linear solver", linear_solver_name(s.linear_solver));
    row("measure runtime", bool_repr(s.measure_runtime));
    row("logical CPUs available", std::to_string(std::thread::hardware_concurrency()), "[machine]");
    head("8. boundary layer (Aungier 2003 Sec. 5.6; not in v1)");
    row("viscous", bool_repr(s.viscous), tag(!s.viscous), s.viscous ? "" : "off: the inviscid solve");
    if (s.viscous) {
        row("coupling", bl_coupling_name(s.bl_coupling), "[user]",
            s.bl_coupling == BlCoupling::OneWay ? "Aungier's analysis: the layer from the inviscid edge state"
                                                : "transpiration (equivalent inviscid flow, Lighthill 1958)");
        row("laminar c_f form", laminar_cf_form_name(s.laminar_cf_form), "[user]",
            s.laminar_cf_form == LaminarCfForm::Centrifugal2000 ? "2 + Lambda/6 (Aungier 2000 Eq. 15-20)"
                                                               : "1 + Lambda/6 (Aungier 2003 Eq. 5-123 as printed)");
        row("wake model", wake_model_name(s.wake_model));
        row("BL relaxation omega_BL", py_repr(s.bl_relaxation));
        row("BL tolerance (max dm/m)", py_repr(s.bl_tolerance));
        row("max BL iterations", std::to_string(s.max_bl_iterations));
        row("transition Re_theta", py_repr(s.transition_re_theta), tag(s.transition_re_theta == 250.0), "Eq. 5-126");
        row("separation H_k", py_repr(s.separation_hk), tag(s.separation_hk == 2.4), "2003 p. 111");
        row("Sutherland mu_ref [Pa s]", py_repr(s.sutherland.mu_ref));
        row("Sutherland T_ref, S [K]", std::format("{}, {}", py_repr(s.sutherland.t_ref), py_repr(s.sutherland.s)));
    }
    head("9. stream surface and row (not in v1)");
    row("stream surface", stream_surface_name(s.stream_surface), "[user]",
        c.stream ? "r(m), b(m) of the BladeGen layer (Aungier Eqs. 5-9 .. 5-26)" : "planar cascade: r = b = 1");
    if (c.stream) {
        const StreamSurfaceBlade& b = *c.stream;
        const StreamSheet& sh = *b.sheet;
        row("source", b.source, "[BladeGen]");
        row("number of blades N", std::to_string(b.n_blades), "[BladeGen]");
        row("BladeGen stagger [deg]", std::format("{:.6f}", b.layer_stagger_deg), "[BladeGen]", "in (m', theta)");
        row("conformal chord c' [rad]", std::format("{:.9g}", b.conformal_chord), "[BladeGen]");
        row("theta_LE [deg]", std::format("{:.6f}", b.theta_le_deg), "[BladeGen]");
        row("r_LE, r_TE [m]", std::format("{:.6g}, {:.6g}", b.r_le, b.r_te), "[BladeGen]");
        const double m0 = c.grid.m.front(), m1 = c.grid.m.back();
        row("r far inlet, far exit [m]", std::format("{:.6g}, {:.6g}", s.chord * sh.radius(m0), s.chord * sh.radius(m1)),
            "[derived]");
        row("b far exit / b far inlet", std::format("{:.6f}", sh.thickness(m1) / sh.thickness(m0)), "[derived]",
            "our model: the channel's streamline spacing");
        row("sin(phi) far inlet, far exit", std::format("{:.6f}, {:.6f}", sh.sin_phi(m0), sh.sin_phi(m1)), "[derived]");
        row("fit residual [chords]", std::format("{:.3g}", b.max_fit_residual), "[derived]");
    }
    row("row", row_motion_name(s.row), "[user]", rotating ? "" : "stationary: beta1 and M1 relative (= absolute)");
    if (rotating && c.triangle) {
        const VelocityTriangle& t = *c.triangle;
        row("rpm", py_repr(s.rpm), "[user]", "> 0: the blades move toward -theta");
        row("alpha1 absolute [deg]", py_repr(s.inlet_absolute_flow_angle_deg), "[user]", "0: axial, no swirl");
        row("M1 absolute", py_repr(s.inlet_absolute_mach), "[user]");
        row("r_in [m]", std::format("{:.9g}", t.radius), "[derived]", "the far inlet");
        row("U = omega r_in [m/s]", std::format("{:.6f}", t.U), "[derived]");
        row("C, C_m, C_theta [m/s]", std::format("{:.6f}, {:.6f}, {:.6f}", t.C, t.C_m, t.C_theta), "[derived]");
        row("W_theta = C_theta - U [m/s]", std::format("{:.6f}", t.W_theta), "[derived]");
        row("beta1 relative [deg]", std::format("{:.6f}", degrees(t.beta)), "[derived]");
        row("M1 relative", std::format("{:.6f}", t.mach_relative), "[derived]");
        row("rothalpy I [J/kg]", std::format("{:.6f}", t.rothalpy), "[derived]", "cp T0 - U C_theta");
    }
    if (c.inlet) {
        try {
            const double re =
                c.inlet->density * c.inlet->speed * s.chord / sutherland_viscosity(c.inlet->temperature, s.sutherland);
            row("chord Reynolds number", std::format("{:.6g}", re), "[derived]", "rho1 W1 c / mu(T1), Sutherland");
        } catch (const BoundaryLayerError&) {
        }
    }
    lines.push_back(rule());
    std::string out;
    for (std::size_t k = 0; k < lines.size(); ++k) out += (k ? "\n" : "") + lines[k];
    return out;
}

std::string format_grid_summary(const Case& c)
{
    const Grid& g = c.grid;
    const CaseSettings& s = c.settings;
    std::vector<std::pair<std::string, std::string>> rows;
    if (c.standard_info) {
        const StandardGridInfo& info = *c.standard_info;
        rows.emplace_back("level", std::format("{} (target dm/c {:.5f}, actual {:.5f}; {} cells between the LE and TE "
                                               "faces)",
                                               info.level, info.target_dm_over_chord, info.dm_over_chord,
                                               info.n_blade_cells));
        rows.emplace_back("(m_f - m_LE)/c", std::format("{:.6f}", info.le_gap_over_chord));
        rows.emplace_back("(m_TE - m_l)/c", std::format("{:.6f}", info.te_gap_over_chord));
        rows.emplace_back("extensions [chords]", std::format("{:.4f} up, {:.4f} down", info.upstream_extension_chords,
                                                             info.downstream_extension_chords));
    }
    rows.emplace_back("n_m", std::to_string(g.n_m()));
    rows.emplace_back("n_eta", std::to_string(g.n_eta()));
    rows.emplace_back("total nodes", std::to_string(g.n_m() * g.n_eta()));
    rows.emplace_back("first blade station", std::to_string(g.first_blade_i));
    rows.emplace_back("last blade station", std::to_string(g.last_blade_i));
    rows.emplace_back("blade stations", std::to_string(g.last_blade_i - g.first_blade_i + 1));
    rows.emplace_back("dm", std::format("{:.6g} chords = {:.6g} m", g.dm(), g.dm() * s.chord));
    rows.emplace_back("deta", std::format("{:.6g}", g.deta()));
    if (g.touching()) {
        rows.emplace_back("edge nodes", "touching (EXPERIMENT, Aungier 2003 Fig. 5-5 'outside', 2000 Fig. 13-3): the first "
                                        "and last blade stations on the nose tip and the tail end, both wall nodes on "
                                        "that point; the side lines end there");
        rows.emplace_back("touching stations", std::format("i = {} at m = {:.6f}, i = {} at m = {:.6f} (chords)", g.first_blade_i,
                                                           g.m_first_blade(), g.last_blade_i, g.m_last_blade()));
        rows.emplace_back("admissibility", "admissible: the stations next to the touching ones each cut one suction and "
                                           "one pressure surface");
    } else {
        rows.emplace_back("edge nodes", "inside (Aungier, v1): the first and last blade stations inside the passage, one "
                                        "node on each surface");
        rows.emplace_back("admissibility", "admissible: first and last blade station each cut one suction and one pressure "
                                           "surface");
    }
    std::string out = rule("GRID");
    for (const auto& [label, value] : rows) out += std::format("\n  {:20s}: {}", label, value);
    return out + "\n" + rule();
}

std::string format_run_report(const RunResult& r)
{
    std::vector<std::string> out{std::string(64, '='), r.converged ? "SOLVER RESULT" : "SOLVER FAILURE",
                                 std::string(64, '='), ""};
    auto block = [&](const std::string& title, const std::vector<std::pair<std::string, std::string>>& rows) {
        out.push_back(title);
        for (const auto& [label, value] : rows) out.push_back(std::format("  {:26s}: {}", label, value));
        out.emplace_back("");
    };
    {
        std::vector<std::pair<std::string, std::string>> rows{{"name", r.case_name}, {"mode", r.mode}, {"status", r.status}};
        if (r.chord_reynolds) rows.emplace_back("chord Reynolds number", std::format("{:.6g}", *r.chord_reynolds));
        block("Case", rows);
    }
    const PassRecord* last = r.records.empty() ? nullptr : &r.records.back();
    if (!r.converged) {
        std::vector<std::pair<std::string, std::string>> rows{{"termination reason", r.error}};
        if (r.failed_pass)
            rows.emplace_back("failed in pass", std::format("{} (the report shows the last completed pass)", *r.failed_pass));
        if (!r.stop_at_sonic) rows.emplace_back("sonic stop", "OFF (stop_at_sonic False)");
        if (last) {
            double peak = last->mach_max;
            for (const PassRecord& p : r.records) peak = std::max(peak, p.mach_max);
            rows.emplace_back("last pass", std::format("{}{}", last->index, last->final ? " (final evaluation)" : ""));
            rows.emplace_back("last beta_out", std::format("{:.6f} deg", degrees(last->exit_angle)));
            rows.emplace_back("last density change", std::format("{:.3e}", last->density_change));
            rows.emplace_back("last velocity change", std::format("{:.3e}", last->velocity_change));
            rows.emplace_back("last Kutta error eps_K", std::format("{:.3e}", last->kutta_normalized));
            rows.emplace_back("peak Mach (last pass)", std::format("{:.4f}", last->mach_max));
            rows.emplace_back("peak Mach (all passes)", std::format("{:.4f}", peak));
        } else {
            rows.emplace_back("last pass", "none recorded");
        }
        block("Failure", rows);
    }
    if (r.solution && r.summary && last) {
        const RunSummary& s = *r.summary;
        const SolutionData& sol = *r.solution;
        if (r.converged)
            block("Convergence", {{"converged", "True"},
                                  {r.mode == "coupled_kutta" ? "coupled passes" : "density passes",
                                   std::format("{} (+ final evaluation)", r.iterations)},
                                  {"density error", std::format("{:.3e}", last->density_change)},
                                  {"velocity error", std::format("{:.3e}", last->velocity_change)},
                                  {"Kutta error eps_K", std::format("{:.3e}", sol.normalized_kutta_residual)},
                                  {"Kutta residual", std::format("{:+.3e} m/s", sol.kutta_residual)}});
        block(std::string("Flow") + (r.converged ? "" : " (last pass, NOT CONVERGED)"),
              {{"beta_in", std::format("{:.4f} deg", s.inlet_angle_deg)},
               {"beta_out", std::format("{:.6f} deg", s.discharge_angle_deg)},
               {"turning", std::format("{:.4f} deg", s.turning_deg)},
               {"max Mach", std::format("{:.4f} at node ({}, {})", s.mach_max, s.mach_max_i, s.mach_max_j)},
               {"min/max density", std::format("{:.6g} / {:.6g} kg/m^3", s.density_min, s.density_max)},
               {"min/max pressure", std::format("{:.6g} / {:.6g} Pa", s.pressure_min, s.pressure_max)},
               {"min/max temperature", std::format("{:.6g} / {:.6g} K", s.temperature_min, s.temperature_max)}});
        block("Surface", {{"SS peak Mach", std::format("{:.4f} at x/c = {:.4f}", s.ss_peak_mach, s.ss_peak_x_c)},
                          {"PS peak Mach", std::format("{:.4f} at x/c = {:.4f}", s.ps_peak_mach, s.ps_peak_x_c)},
                          {"TE SS Mach", std::format("{:.4f}", s.te_ss_mach)},
                          {"TE PS Mach", std::format("{:.4f}", s.te_ps_mach)}});
        block("Linear algebra", {{"||F Psi - Q||_inf", std::format("{:.3e}", last->linear_residual)},
                                 {"solver", r.runtime.linear_solver}});
        block("Mass flow", {{"mdot/depth", std::format("{:.6f} kg/(m s)", s.mass_flow_per_depth)},
                            {"source", s.mass_flow_source},
                            {"max station error", std::format("{:.3e}", last->mass_flow_error_max)},
                            {"blade median error", std::format("{:.3e}", last->mass_flow_error_blade_median)}});
        if (s.rotor) {
            const RotorSummary& o = *s.rotor;
            block("Rotor, absolute frame (far exit, mass averages)",
                  {{"rpm, omega", std::format("{:g}, {:.6f} rad/s", o.rpm, o.omega)},
                   {"U in, out", std::format("{:.4f}, {:.4f} m/s", o.U_in, o.U_out)},
                   {"C_m, C_theta out", std::format("{:.4f}, {:.4f} m/s", o.C_m_out, o.C_theta_out)},
                   {"alpha out", std::format("{:.4f} deg", o.alpha_out_deg)},
                   {"M out (absolute)", std::format("{:.4f}", o.mach_out)},
                   {"P0 out, P0 out / P0 in", std::format("{:.6g} Pa, {:.6f}", o.total_pressure_out, o.total_pressure_ratio)},
                   {"T0 out, T0 out / T0 in", std::format("{:.6g} K, {:.6f}", o.total_temperature_out, o.total_temperature_ratio)},
                   {"Euler work U2 Ct2 - U1 Ct1", std::format("{:.4f} J/kg", o.euler_work)},
                   {"cp (T0 out - T0 in)", std::format("{:.4f} J/kg", o.total_enthalpy_rise)}});
        }
    } else if (r.solution && !r.postprocess_error.empty()) {
        block("Postprocess", {{"not available", r.postprocess_error}});
    }
    if (r.viscous) {
        const ViscousData& v = *r.viscous;
        auto xc = [](double x) { return std::isfinite(x) ? std::format("{:.4f}", x) : std::string("none"); };
        std::vector<std::pair<std::string, std::string>> rows{
            {"coupling", v.coupling},
            {"laminar c_f form", v.laminar_cf_form},
            {"wake model", v.wake_model},
            {"BL iterations", std::format("{} ({})", v.history.size(), v.converged ? "converged" : "NOT converged")},
            {"layers reported", v.layers_iteration < 0      ? std::string("none (no layer was computed)")
                                : v.layers_of_last_solve ? std::format("BL iteration {}", v.layers_iteration)
                                                         : std::format("BL iteration {}, the last one computed (a later "
                                                                       "inviscid solve failed)",
                                                                       v.layers_iteration)}};
        if (!v.history.empty()) {
            const BlIteration& h = v.history.back();
            rows.emplace_back("last max change of m", std::format("{:.3e}", h.max_change_m));
            rows.emplace_back("beta_out inviscid", std::format("{:.6f} deg", degrees(v.inviscid_exit_angle)));
            rows.emplace_back("beta_out reported", std::format("{:.6f} deg", degrees(h.beta_out)));
            rows.emplace_back("SS transition x/c", xc(h.transition_x_c_ss));
            rows.emplace_back("PS transition x/c", xc(h.transition_x_c_ps));
            rows.emplace_back("SS separation x/c", xc(h.separation_x_c_ss));
            rows.emplace_back("PS separation x/c", xc(h.separation_x_c_ps));
            if (std::isfinite(h.laminar_separation_x_c_ss) || std::isfinite(h.laminar_separation_x_c_ps))
                rows.emplace_back("laminar separation x/c",
                                  std::format("SS {}, PS {}", xc(h.laminar_separation_x_c_ss), xc(h.laminar_separation_x_c_ps)));
        }
        rows.emplace_back("loss omega (Eq. 5-139)", std::format("{:.5f}", v.loss_omega));
        rows.emplace_back("Theta, Delta*", std::format("{:.5f}, {:.5f}", v.loss_Theta, v.loss_Delta_star));
        block("Boundary layer (Aungier 2003 Sec. 5.6)", rows);
    }
    const RuntimeStats& rt = r.runtime;
    block("Runtime", {{"solver wall time", std::format("{:.3f} s", rt.wall_time_s)},
                      {"process CPU time", std::format("{:.3f} s", rt.process_cpu_time_s)},
                      {"effective CPU cores", std::format("{:.2f}", rt.effective_cpu_cores)},
                      {"threads", std::to_string(rt.threads)},
                      {"logical CPUs available", std::to_string(rt.logical_cpus_available)}});
    out.push_back(std::string(64, '='));
    std::string text;
    for (std::size_t k = 0; k < out.size(); ++k) text += (k ? "\n" : "") + out[k];
    return text;
}

// --- files ------------------------------------------------------------------------------------

namespace {

using nlohmann::json;

json num(double v) { return std::isfinite(v) ? json(v) : json(nullptr); }
json num(const std::optional<double>& v) { return v ? num(*v) : json(nullptr); }

json record_json(const PassRecord& p)
{
    return {{"index", p.index},
            {"final", p.final},
            {"density_change", num(p.density_change)},
            {"velocity_change", num(p.velocity_change)},
            {"psi_change", num(p.psi_change)},
            {"target_mismatch", num(p.target_mismatch)},
            {"density_min", num(p.density_min)},
            {"density_max", num(p.density_max)},
            {"speed_min", num(p.speed_min)},
            {"speed_max", num(p.speed_max)},
            {"mach_max", num(p.mach_max)},
            {"linear_residual", num(p.linear_residual)},
            {"mass_flow_error_max", num(p.mass_flow_error_max)},
            {"mass_flow_error_blade_median", num(p.mass_flow_error_blade_median)},
            {"exit_angle", num(p.exit_angle)},
            {"kutta_residual", num(p.kutta_residual)},
            {"kutta_normalized", num(p.kutta_normalized)},
            {"exit_angle_next", num(p.exit_angle_next)},
            {"tan_correction", num(p.tan_correction)}};
}

json settings_json(const CaseSettings& s)
{
    auto opt_s = [](const std::optional<std::string>& v) { return v ? json(*v) : json(nullptr); };
    auto opt_i = [](const std::optional<int>& v) { return v ? json(*v) : json(nullptr); };
    return {{"case_name", s.case_name},
            {"solver_mode", mode_name(s.solver_mode)},
            {"fixed_exit_angle_deg", num(s.fixed_exit_angle_deg)},
            {"gamma", num(s.gamma)},
            {"gas_constant", num(s.gas_constant)},
            {"inlet_total_pressure", num(s.inlet_total_pressure)},
            {"inlet_total_temperature", num(s.inlet_total_temperature)},
            {"omega", num(s.omega)},
            {"chord", num(s.chord)},
            {"solidity", num(s.solidity)},
            {"stagger_deg", num(s.stagger_deg)},
            {"inlet_mach", num(s.inlet_mach)},
            {"inlet_flow_angle_deg", num(s.inlet_flow_angle_deg)},
            {"mass_flow_per_depth_override", num(s.mass_flow_per_depth_override)},
            {"grid_level", opt_s(s.grid_level)},
            {"n_passage", opt_i(s.n_passage)},
            {"n_eta", opt_i(s.n_eta)},
            {"upstream_extension_chords", num(s.upstream_extension_chords)},
            {"downstream_extension_chords", num(s.downstream_extension_chords)},
            {"density_tolerance", num(s.density_tolerance)},
            {"velocity_tolerance", num(s.velocity_tolerance)},
            {"max_density_iterations", s.max_density_iterations},
            {"density_relaxation", num(s.density_relaxation)},
            {"stop_at_sonic", s.stop_at_sonic},
            {"kutta_tolerance", num(s.kutta_tolerance)},
            {"max_kutta_iterations", s.max_kutta_iterations},
            {"kutta_relaxation", num(s.kutta_relaxation)},
            {"initial_exit_angle_mode",
             s.initial_exit_angle_mode == ExitAngleStart::BladeTeAngle ? "blade_te_angle" : "manual"},
            {"initial_exit_angle_deg", num(s.initial_exit_angle_deg)},
            {"threads", s.threads},
            {"linear_solver", s.linear_solver == LinearSolverKind::Pardiso ? "pardiso" : "sparselu"},
            {"measure_runtime", s.measure_runtime},
            {"stream_surface", stream_surface_name(s.stream_surface)},
            {"row", row_motion_name(s.row)},
            {"rpm", num(s.rpm)},
            {"inlet_absolute_mach", num(s.inlet_absolute_mach)},
            {"inlet_absolute_flow_angle_deg", num(s.inlet_absolute_flow_angle_deg)}};
}

std::string csv_line(const std::vector<std::string>& cells)
{
    std::string s;
    for (std::size_t k = 0; k < cells.size(); ++k) s += (k ? "," : "") + cells[k];
    return s + "\r\n";  // Python's csv.writer line terminator
}

}  // namespace

std::vector<std::string> export_run(const RunResult& r, const Case& c, const std::string& directory, bool overwrite)
{
    namespace fs = std::filesystem;
    const fs::path dir(directory);
    std::vector<std::string> files;
    const bool sol = r.solution.has_value();
    const bool surfaces = sol && !r.solution->suction.x_c.empty();
    if (surfaces) files.emplace_back("surface_mach.csv");
    if (surfaces && r.solution->inlet) files.emplace_back("surface_pressure.csv");
    if (sol) {
        files.emplace_back("convergence_history.csv");
        files.emplace_back("fields.csv");
    }
    const bool layers = r.viscous && !r.viscous->suction.x_c.empty();
    if (layers) {
        files.emplace_back("bl_suction.csv");
        files.emplace_back("bl_pressure.csv");
    }
    if (r.viscous) files.emplace_back("bl_history.csv");
    files.emplace_back("run_summary.json");
    if (!overwrite)
        for (const std::string& f : files)
            if (fs::exists(dir / f))
                throw std::runtime_error(std::format("{} exists: export_run overwrites only when asked", (dir / f).string()));
    fs::create_directories(dir);
    auto open = [&](const std::string& name) {
        std::ofstream out(dir / name, std::ios::binary);
        if (!out) throw std::runtime_error("cannot write " + (dir / name).string());
        return out;
    };
    if (surfaces) {
        const SolutionData& d = *r.solution;
        std::ofstream out = open("surface_mach.csv");
        out << csv_line({"surface", "x_c", "m_c", "mach", "pressure", "dunavant_S"});
        for (const auto& [name, s] : {std::pair<const char*, const SurfaceData*>{"suction", &d.suction},
                                      std::pair<const char*, const SurfaceData*>{"pressure", &d.pressure_side}})
            for (std::size_t k = 0; k < s->x_c.size(); ++k)
                out << csv_line({name, py_repr(s->x_c[k]), py_repr(s->m_c[k]), py_repr(s->mach[k]),
                                 py_repr(s->pressure[k]), s->dunavant_S.empty() ? "" : py_repr(s->dunavant_S[k])});
        if (d.inlet) {
            std::ofstream p = open("surface_pressure.csv");
            p << csv_line({"surface", "x_c", "dunavant_S", "inlet_static_cp"});
            for (const auto& [name, s] : {std::pair<const char*, const SurfaceData*>{"suction", &d.suction},
                                          std::pair<const char*, const SurfaceData*>{"pressure", &d.pressure_side}})
                for (std::size_t k = 0; k < s->x_c.size(); ++k)
                    p << csv_line({name, py_repr(s->x_c[k]), py_repr(s->dunavant_S[k]), py_repr(s->cp_static[k])});
        }
    }
    if (sol) {
        std::ofstream out = open("convergence_history.csv");
        std::vector<std::string> head{"iteration", "final", "density_change", "velocity_change", "psi_change", "kutta_error",
                                      "kutta_residual", "exit_angle_deg", "exit_angle_next_deg", "tan_correction",
                                      "linear_residual", "peak_mach", "target_mismatch", "density_min", "density_max",
                                      "speed_min", "speed_max", "mass_flow_error_max", "mass_flow_error_blade_median"};
        if (r.viscous) head.emplace_back("bl_iteration");  // viscous runs only: the inviscid layout is v1's
        out << csv_line(head);
        for (const PassRecord& p : r.records) {
            std::vector<std::string> cells{std::to_string(p.index), p.final ? "1" : "0"};
            for (double v : {p.density_change, p.velocity_change, p.psi_change, p.kutta_normalized, p.kutta_residual,
                             degrees(p.exit_angle), degrees(p.exit_angle_next), p.tan_correction, p.linear_residual,
                             p.mach_max, p.target_mismatch, p.density_min, p.density_max, p.speed_min, p.speed_max,
                             p.mass_flow_error_max, p.mass_flow_error_blade_median})
                cells.push_back(py_repr(v));
            if (r.viscous) cells.push_back(std::to_string(p.bl_iteration));
            out << csv_line(cells);
        }
        // every node of the reported evaluation; m and y in metres (c_ref times the chord-normalized grid)
        const SolutionData& d = *r.solution;
        std::ofstream f = open("fields.csv");
        std::vector<std::string> head_f{"i", "j", "m", "y", "psi", "rho", "p", "T", "M", "W", "W_m", "W_theta"};
        const bool sheet = !d.station_radius.empty();  // not in v1: y = r_LE theta, and r, theta, b of the sheet
        if (sheet) head_f.insert(head_f.end(), {"r", "theta", "b"});
        if (d.absolute) head_f.insert(head_f.end(), {"C_theta", "C", "alpha_deg", "M_abs", "P0_abs", "T0_abs"});
        f << csv_line(head_f);
        const double r_ref = sheet && c.stream ? c.stream->sheet->r_ref() : 1.0;
        for (int i = 0; i < d.grid.n_m(); ++i)
            for (int j = 0; j < d.grid.n_eta(); ++j) {
                std::vector<std::string> cells{std::to_string(i), std::to_string(j), py_repr(d.c_ref * d.grid.node_m(i, j)),
                                               py_repr(d.c_ref * d.grid.node_y(i, j)), py_repr(d.psi(i, j)),
                                               py_repr(d.density(i, j)), py_repr(d.pressure(i, j)),
                                               py_repr(d.temperature(i, j)), py_repr(d.mach(i, j)), py_repr(d.W(i, j)),
                                               py_repr(d.W_m(i, j)), py_repr(d.W_theta(i, j))};
                if (sheet) {
                    cells.push_back(py_repr(d.station_radius[sz(i)]));
                    cells.push_back(py_repr(d.grid.node_y(i, j) / r_ref));
                    cells.push_back(py_repr(d.station_thickness[sz(i)]));
                }
                if (d.absolute) {
                    const AbsoluteFrameData& a = *d.absolute;
                    for (double v : {a.C_theta(i, j), a.C(i, j), degrees(a.alpha(i, j)), a.mach(i, j),
                                     a.total_pressure(i, j), a.total_temperature(i, j)})
                        cells.push_back(py_repr(v));
                }
                f << csv_line(cells);
            }
    }
    if (layers) {
        const ViscousData& v = *r.viscous;
        for (const auto& [name, d] : {std::pair<const char*, const BoundaryLayerSurfaceData*>{"bl_suction.csv", &v.suction},
                                      std::pair<const char*, const BoundaryLayerSurfaceData*>{"bl_pressure.csv", &v.pressure_side}}) {
            std::ofstream out = open(name);
            out << csv_line({"x_c", "x", "W_e", "rho_e", "T_e", "mu", "theta", "delta_star", "delta", "H", "H_k", "c_f",
                             "Lambda", "Re_theta", "regime", "mass_deficit"});
            for (std::size_t k = 0; k < d->x_c.size(); ++k) {
                std::vector<std::string> cells;
                for (const std::vector<double>* col : {&d->x_c, &d->x, &d->W_e, &d->rho_e, &d->T_e, &d->mu, &d->theta,
                                                       &d->delta_star, &d->delta, &d->H, &d->H_k, &d->c_f, &d->Lambda,
                                                       &d->Re_theta})
                    cells.push_back(py_repr((*col)[k]));
                cells.push_back(d->regime[k]);
                cells.push_back(py_repr(d->mass_deficit[k]));
                out << csv_line(cells);
            }
        }
    }
    if (r.viscous) {
        std::ofstream out = open("bl_history.csv");
        out << csv_line({"k", "inviscid_passes", "beta_out_deg", "max_change_m", "transition_x_c_ss", "transition_x_c_ps",
                         "separation_x_c_ss", "separation_x_c_ps", "laminar_separation_x_c_ss", "laminar_separation_x_c_ps",
                         "loss_omega"});
        for (const BlIteration& h : r.viscous->history) {
            std::vector<std::string> cells{std::to_string(h.k), std::to_string(h.inviscid_passes)};
            for (double x : {degrees(h.beta_out), h.max_change_m, h.transition_x_c_ss, h.transition_x_c_ps,
                             h.separation_x_c_ss, h.separation_x_c_ps, h.laminar_separation_x_c_ss,
                             h.laminar_separation_x_c_ps, h.loss})
                cells.push_back(py_repr(x));
            out << csv_line(cells);
        }
    }
    // run_summary.json (v1 controller.run_summary)
    const CaseSettings& s = c.settings;
    const Grid& g = c.grid;
    const auto [angle, source] = exit_angle_of(s, degrees(c.cascade.kappa_te));
    json inlet = nullptr;
    if (c.inlet)
        inlet = {{"mach", num(c.inlet->mach)},
                 {"flow_angle", num(c.inlet->angle)},
                 {"flow_angle_deg", num(degrees(c.inlet->angle))},
                 {"total_pressure", num(c.inlet->total_pressure)},
                 {"static_temperature", num(c.inlet->temperature)},
                 {"static_pressure", num(c.inlet->pressure)},
                 {"static_density", num(c.inlet->density)},
                 {"sound_speed", num(c.inlet->sound_speed)},
                 {"velocity", num(c.inlet->speed)},
                 {"dynamic_pressure", num(c.inlet->dynamic_pressure)},
                 {"pitch", num(s.chord * g.pitch())},
                 {"mass_flow_per_depth", num(c.inlet->mass_flow_per_depth)}};
    json standard = nullptr;
    if (c.standard_info) {
        const StandardGridInfo& i = *c.standard_info;
        standard = {{"level", i.level},
                    {"target_dm_over_chord", num(i.target_dm_over_chord)},
                    {"n_blade_cells", i.n_blade_cells},
                    {"dm", num(i.dm)},
                    {"dm_over_chord", num(i.dm_over_chord)},
                    {"deta", num(i.deta)},
                    {"n_m", i.n_m},
                    {"n_eta", i.n_eta},
                    {"total_nodes", i.total_nodes},
                    {"first_blade_i", i.first_blade_i},
                    {"last_blade_i", i.last_blade_i},
                    {"n_blade_stations", i.n_blade_stations},
                    {"le_gap_over_chord", num(i.le_gap_over_chord)},
                    {"te_gap_over_chord", num(i.te_gap_over_chord)},
                    {"upstream_extension_chords", num(i.upstream_extension_chords)},
                    {"downstream_extension_chords", num(i.downstream_extension_chords)},
                    {"cell_aspect", num(i.cell_aspect)},
                    {"checks", i.checks}};
    }
    json results = nullptr;
    if (r.solution && r.summary && !r.records.empty()) {
        const RunSummary& m = *r.summary;
        double peak = r.records.front().mach_max;
        for (const PassRecord& p : r.records) peak = std::max(peak, p.mach_max);
        const PassRecord& last = r.records.back();
        results = {{"iterations", r.iterations},
                   {"passes_recorded", r.records.size()},
                   {"limitation", m.limitation},
                   {"beta_in_deg", num(m.inlet_angle_deg)},
                   {"beta_out_deg", num(m.discharge_angle_deg)},
                   {"exit_velocity_angle_deg", num(m.exit_velocity_angle_deg)},
                   {"turning_deg", num(m.turning_deg)},
                   {"kutta_residual_m_s", num(m.kutta_residual)},
                   {"normalized_kutta_residual", num(m.normalized_kutta_residual)},
                   {"surface",
                    {{"ss_peak_mach", num(m.ss_peak_mach)},
                     {"ss_peak_x_c", num(m.ss_peak_x_c)},
                     {"ps_peak_mach", num(m.ps_peak_mach)},
                     {"ps_peak_x_c", num(m.ps_peak_x_c)},
                     {"te_ss_mach", num(m.te_ss_mach)},
                     {"te_ps_mach", num(m.te_ps_mach)},
                     {"min_mach", num(m.surface_mach_min)},
                     {"max_mach", num(m.surface_mach_max)}}},
                   {"field",
                    {{"max_mach", num(m.mach_max)},
                     {"max_mach_node", {m.mach_max_i, m.mach_max_j}},
                     {"density_min", num(m.density_min)},
                     {"density_max", num(m.density_max)},
                     {"pressure_min", num(m.pressure_min)},
                     {"pressure_max", num(m.pressure_max)},
                     {"temperature_min", num(m.temperature_min)},
                     {"temperature_max", num(m.temperature_max)}}},
                   {"mass_flow",
                    {{"per_depth_kg_m_s", num(m.mass_flow_per_depth)},
                     {"source", m.mass_flow_source},
                     {"max_station_error", num(last.mass_flow_error_max)},
                     {"blade_median_error", num(last.mass_flow_error_blade_median)}}},
                   {"final_pass", record_json(last)},
                   {"peak_mach_all_passes", num(peak)}};
    }
    const std::string type = r.error.empty() ? "" : r.error.substr(0, r.error.find(':'));
    const std::string message = r.error.empty() ? "" : r.error.substr(std::min(r.error.size(), type.size() + 2));
    json data = {
        {"format", "turbolenza.b2b run summary, version 1"},
        {"producer", "TURBOlenza v2 (BladeGen, C++ port of turbolenza.b2b v1)"},
        {"started", r.started},
        {"status", r.status},
        {"termination",
         {{"converged", r.converged},
          {"error_type", r.error.empty() ? json(nullptr) : json(type)},
          {"error_message", r.error.empty() ? json(nullptr) : json(message)}}},
        {"inputs", settings_json(s)},
        {"derived",
         {{"blade_name", c.section.designation},
          {"pitch_m", num(s.chord * g.pitch())},
          {"angle_of_attack_deg", num(s.inlet_flow_angle_deg - s.stagger_deg)},
          {"kappa_le_deg", num(degrees(c.cascade.kappa_le))},
          {"kappa_te_deg", num(degrees(c.cascade.kappa_te))},
          {"exit_angle_deg", num(angle)},
          {"exit_angle_source", source},
          {"mass_flow_per_depth_kg_m_s", num(c.mass_flow_per_depth)},
          {"mass_flow_source", c.mass_flow_source},
          {"pass_limit", pass_limit(s)},
          {"inlet_state", inlet},
          {"chord_reynolds_number", num(r.chord_reynolds)}}},
        {"grid",
         {{"grid_level", s.grid_level ? json(*s.grid_level) : json(nullptr)},
          {"standard_grid", standard},
          {"n_passage", s.n_passage ? json(*s.n_passage) : json(nullptr)},
          {"n_eta_setting", s.n_eta ? json(*s.n_eta) : json(nullptr)},
          {"n_m", g.n_m()},
          {"n_eta", g.n_eta()},
          {"total_nodes", g.n_m() * g.n_eta()},
          {"first_blade_station", g.first_blade_i},
          {"last_blade_station", g.last_blade_i},
          {"blade_stations", g.last_blade_i - g.first_blade_i + 1},
          {"dm_chords", num(g.dm())},
          {"dm_m", num(g.dm() * s.chord)},
          {"deta", num(g.deta())},
          {"admissible", true}}},
        {"results", results},
        {"runtime",
         {{"wall_time_s", num(r.runtime.wall_time_s)},
          {"process_cpu_time_s", num(r.runtime.process_cpu_time_s)},
          {"effective_cpu_cores", num(r.runtime.effective_cpu_cores)},
          {"logical_cpus_available", r.runtime.logical_cpus_available},
          {"threads", r.runtime.threads},
          {"linear_solver", r.runtime.linear_solver}}},
        {"environment",
         {{"eigen", std::format("{}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION)},
          {"openmp",
#ifdef TURBOLENZA_HAVE_OPENMP
           true
#else
           false
#endif
          },
          {"pardiso", pardiso_available()},
          {"compiler",
#if defined(_MSC_VER)
           std::format("MSVC {}", _MSC_VER)
#elif defined(__clang__)
           std::format("clang {}.{}", __clang_major__, __clang_minor__)
#elif defined(__GNUC__)
           std::format("gcc {}.{}", __GNUC__, __GNUC_MINOR__)
#else
           "unknown"
#endif
          }}},
    };
    if (s.viscous) {
        data["inputs"]["boundary_layer"] = {{"viscous", true},
                                            {"bl_coupling", bl_coupling_name(s.bl_coupling)},
                                            {"laminar_cf_form", laminar_cf_form_name(s.laminar_cf_form)},
                                            {"wake_model", wake_model_name(s.wake_model)},
                                            {"bl_relaxation", num(s.bl_relaxation)},
                                            {"bl_tolerance", num(s.bl_tolerance)},
                                            {"max_bl_iterations", s.max_bl_iterations},
                                            {"transition_re_theta", num(s.transition_re_theta)},
                                            {"separation_hk", num(s.separation_hk)},
                                            {"sutherland_mu_ref", num(s.sutherland.mu_ref)},
                                            {"sutherland_t_ref", num(s.sutherland.t_ref)},
                                            {"sutherland_s", num(s.sutherland.s)}};
    }
    if (r.viscous) {
        const ViscousData& v = *r.viscous;
        json bl = {{"converged", v.converged},
                   {"iterations", v.history.size()},
                   {"loss_omega", num(v.loss_omega)},
                   {"loss_Theta", num(v.loss_Theta)},
                   {"loss_Delta_star", num(v.loss_Delta_star)},
                   {"beta_out_inviscid_deg", num(degrees(v.inviscid_exit_angle))},
                   {"te_jump_mismatch", num(v.te_jump_mismatch)}};
        if (!v.history.empty()) {
            const BlIteration& h = v.history.back();
            bl["beta_out_deg"] = num(degrees(h.beta_out));
            bl["max_change_m"] = num(h.max_change_m);
            bl["transition_x_c"] = {{"suction", num(h.transition_x_c_ss)}, {"pressure", num(h.transition_x_c_ps)}};
            bl["separation_x_c"] = {{"suction", num(h.separation_x_c_ss)}, {"pressure", num(h.separation_x_c_ps)}};
            bl["laminar_separation_x_c"] = {{"suction", num(h.laminar_separation_x_c_ss)},
                                            {"pressure", num(h.laminar_separation_x_c_ps)}};
        }
        data["boundary_layer"] = bl;
    }
    // not in v1: the stream surface and the rotating row
    if (c.stream) {
        const StreamSurfaceBlade& b = *c.stream;
        const StreamSheet& sh = *b.sheet;
        data["derived"]["stream_surface"] = {
            {"source", b.source},
            {"n_blades", b.n_blades},
            {"chord_m", num(b.chord)},
            {"solidity", num(b.solidity)},
            {"stagger_deg", num(b.stagger_deg)},
            {"bladegen_stagger_deg", num(b.layer_stagger_deg)},
            {"conformal_chord_rad", num(b.conformal_chord)},
            {"theta_le_deg", num(b.theta_le_deg)},
            {"r_le_m", num(b.r_le)},
            {"r_te_m", num(b.r_te)},
            {"r_ref_chords", num(sh.r_ref())},
            {"r_far_inlet_m", num(s.chord * sh.radius(g.m.front()))},
            {"r_far_exit_m", num(s.chord * sh.radius(g.m.back()))},
            {"b_far_exit_over_far_inlet", num(sh.thickness(g.m.back()) / sh.thickness(g.m.front()))},
            {"cylinder", sh.cylinder()},
            {"max_fit_residual_chords", num(b.max_fit_residual)}};
    }
    if (c.triangle) {
        const VelocityTriangle& t = *c.triangle;
        data["derived"]["inlet_triangle"] = {{"rpm", num(t.rpm)},
                                             {"omega_rad_s", num(t.omega)},
                                             {"r_in_m", num(t.radius)},
                                             {"U_m_s", num(t.U)},
                                             {"alpha1_deg", num(degrees(t.alpha))},
                                             {"M1_absolute", num(t.mach_absolute)},
                                             {"C_m_s", num(t.C)},
                                             {"C_m_m_s", num(t.C_m)},
                                             {"C_theta_m_s", num(t.C_theta)},
                                             {"W_theta_m_s", num(t.W_theta)},
                                             {"W_m_s", num(t.W)},
                                             {"beta1_deg", num(degrees(t.beta))},
                                             {"M1_relative", num(t.mach_relative)},
                                             {"rothalpy_J_kg", num(t.rothalpy)}};
    }
    if (r.summary && r.summary->rotor && results.is_object()) {
        const RotorSummary& o = *r.summary->rotor;
        data["results"]["rotor_absolute"] = {{"U_in_m_s", num(o.U_in)},
                                             {"U_out_m_s", num(o.U_out)},
                                             {"C_m_out_m_s", num(o.C_m_out)},
                                             {"C_theta_out_m_s", num(o.C_theta_out)},
                                             {"alpha_out_deg", num(o.alpha_out_deg)},
                                             {"mach_out", num(o.mach_out)},
                                             {"total_pressure_out_pa", num(o.total_pressure_out)},
                                             {"total_temperature_out_k", num(o.total_temperature_out)},
                                             {"total_pressure_ratio", num(o.total_pressure_ratio)},
                                             {"total_temperature_ratio", num(o.total_temperature_ratio)},
                                             {"euler_work_j_kg", num(o.euler_work)},
                                             {"total_enthalpy_rise_j_kg", num(o.total_enthalpy_rise)}};
    }
    json listed = json::array();
    for (const std::string& f : files) listed.push_back(f);
    data["files"] = listed;
    std::ofstream out = open("run_summary.json");
    out << data.dump(2);
    return files;
}

}  // namespace BladeGen::Turbolenza
