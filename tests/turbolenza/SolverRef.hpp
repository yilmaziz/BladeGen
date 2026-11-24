// tests/turbolenza/SolverRef.hpp — helpers of the solver-port tests: the v1 objects of a dumped
// pass as port objects (so that each module is checked from v1's own inputs), the measured
// differences printed as "MEASURE <name> <value>" lines (measure before tolerance, PORT.md rule
// 3), and a manufactured straight channel grid for the closed-form tests that need no v1 data.
#pragma once

#include "RefData.hpp"

#include "turbolenza/Coefficients.hpp"
#include "turbolenza/Flow.hpp"
#include "turbolenza/Metrics.hpp"
#include "turbolenza/Solver.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace sref {

namespace T = BladeGen::Turbolenza;
using ref::json;

/// Prints the measured value (grep "MEASURE" in the test output) and returns it.
inline double measure(const std::string& name, double value)
{
    std::printf("MEASURE %-48s %.3e\n", name.c_str(), value);
    return value;
}

inline std::shared_ptr<const T::Grid> grid(const json& g) { return std::make_shared<const T::Grid>(ref::grid(g)); }

/// v1 GridMetrics as the port's struct (the v1 fields themselves).
inline T::GridMetrics metrics(const json& m, std::shared_ptr<const T::Grid> g)
{
    T::GridMetrics out;
    out.grid = std::move(g);
    out.S = ref::vec(m.at("S"));
    out.dS_dm = ref::vec(m.at("dS_dm"));
    out.tan_beta0 = ref::vec(m.at("tan_beta0"));
    out.tan_beta1 = ref::vec(m.at("tan_beta1"));
    out.tan_beta = ref::field(m.at("tan_beta"));
    out.beta = ref::field(m.at("beta"));
    out.sin_beta = ref::field(m.at("sin_beta"));
    out.cos_beta = ref::field(m.at("cos_beta"));
    out.dtan_beta0_dm = ref::vec(m.at("dtan_beta0_dm"));
    out.dtan_beta1_dm = ref::vec(m.at("dtan_beta1_dm"));
    out.dtan_beta_dm = ref::field(m.at("dtan_beta_dm"));
    return out;
}

/// v1 CoefficientFields.
inline T::CoefficientFields coefficients(const json& c)
{
    T::CoefficientFields out{ref::field(c.at("A")), ref::field(c.at("B")), ref::field(c.at("C")), ref::field(c.at("D")),
                             ref::field(c.at("E")), ref::num(c.at("mass_flow_per_depth")), std::nullopt};
    if (!c.at("c_ref").is_null()) out.c_ref = ref::num(c.at("c_ref"));
    return out;
}

/// v1 FlowField.
inline T::FlowField flow(const json& f)
{
    return {ref::field(f.at("psi_m")), ref::field(f.at("psi_eta")), ref::field(f.at("W_m")),
            ref::field(f.at("W_theta")), ref::field(f.at("W2")),     ref::field(f.at("W"))};
}

/// max |a - b| over the entries.
inline double abs_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    REQUIRE(a.size() == b.size());
    double m = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) m = std::max(m, std::abs(a[k] - b[k]));
    return m;
}

/// |a - b| / max(|b|, tiny) for two numbers (NaN == NaN counts as equal).
inline double rel(double a, double b)
{
    if (std::isnan(a) && std::isnan(b)) return 0.0;
    const double d = std::abs(a - b);
    return std::abs(b) > 0.0 ? d / std::abs(b) : d;
}

/// The worst difference of the numbers of two PassRecords (NaN at the same place equal): relative
/// for the levels (densities, speeds, Mach, angles); absolute for the dimensionless small
/// quantities that go to zero at convergence (the change norms, mass-flow errors, eps_K, the
/// tan correction), whose relative difference is the solve's round-off over a vanishing value;
/// R_K relative to the largest speed.
inline double record_diff(const T::PassRecord& a, const json& b, std::string* worst = nullptr)
{
    const std::pair<const char*, double> relative[] = {{"density_min", a.density_min}, {"density_max", a.density_max},
                                                       {"speed_min", a.speed_min},     {"speed_max", a.speed_max},
                                                       {"mach_max", a.mach_max},       {"exit_angle", a.exit_angle},
                                                       {"exit_angle_next", a.exit_angle_next}};
    const std::pair<const char*, double> absolute[] = {{"density_change", a.density_change},
                                                       {"velocity_change", a.velocity_change},
                                                       {"psi_change", a.psi_change},
                                                       {"target_mismatch", a.target_mismatch},
                                                       {"mass_flow_error_max", a.mass_flow_error_max},
                                                       {"mass_flow_error_blade_median", a.mass_flow_error_blade_median},
                                                       {"kutta_normalized", a.kutta_normalized},
                                                       {"tan_correction", a.tan_correction}};
    double m = 0.0;
    auto take = [&](const char* name, double d) {
        if (std::isnan(d)) d = INFINITY;
        if (d > m) {
            m = d;
            if (worst) *worst = name;
        }
    };
    for (const auto& [name, v] : relative) take(name, rel(v, ref::num(b.at(name))));
    for (const auto& [name, v] : absolute) {
        const double r = ref::num(b.at(name));
        take(name, std::isnan(v) && std::isnan(r) ? 0.0 : std::abs(v - r));
    }
    const double rk = ref::num(b.at("kutta_residual"));
    take("kutta_residual", std::isnan(a.kutta_residual) && std::isnan(rk) ? 0.0
                                                                          : std::abs(a.kutta_residual - rk) / a.speed_max);
    return m;
}

/// The v1 settings of a dump as a SolverConfig.
inline T::SolverConfig config(const json& d)
{
    const json& c = d.at("config");
    T::SolverConfig out;
    const json& f = c.at("flow");
    if (!f.at("inlet_mach").is_null()) out.flow.inlet_mach = ref::num(f.at("inlet_mach"));
    out.flow.inlet_angle_deg = ref::num(f.at("inlet_angle_deg"));
    out.flow.exit_angle_guess_deg = ref::num(f.at("exit_angle_guess_deg"));
    out.flow.inlet_absolute_total_pressure = ref::num(f.at("inlet_absolute_total_pressure"));
    out.flow.inlet_absolute_total_temperature = ref::num(f.at("inlet_absolute_total_temperature"));
    out.flow.angular_speed = ref::num(f.at("angular_speed"));
    out.gas = T::PerfectGas(ref::num(c.at("gas").at("gamma")), ref::num(c.at("gas").at("gas_constant")));
    const json& it = c.at("iteration");
    out.iteration.density_tolerance = ref::num(it.at("density_tolerance"));
    out.iteration.max_density_iterations = it.at("max_density_iterations").get<int>();
    out.iteration.kutta_tolerance = ref::num(it.at("kutta_tolerance"));
    out.iteration.max_kutta_iterations = it.at("max_kutta_iterations").get<int>();
    out.iteration.density_relaxation = ref::num(it.at("density_relaxation"));
    if (!it.at("velocity_tolerance").is_null()) out.iteration.velocity_tolerance = ref::num(it.at("velocity_tolerance"));
    out.iteration.kutta_relaxation = ref::num(it.at("kutta_relaxation"));
    out.iteration.stop_at_sonic = it.at("stop_at_sonic").get<bool>();
    const json& g = c.at("grid");
    out.grid.upstream_extension_chords = ref::num(g.at("upstream_extension_chords"));
    out.grid.downstream_extension_chords = ref::num(g.at("downstream_extension_chords"));
    if (!g.at("level").is_null()) out.grid.level = g.at("level").get<std::string>();
    if (!g.at("n_passage").is_null()) out.grid.n_passage = g.at("n_passage").get<int>();
    if (!g.at("n_eta").is_null()) out.grid.n_eta = g.at("n_eta").get<int>();
    return out;
}

/// A straight channel: side lines and walls on y0 = tan(angle) m, y1 = y0 + pitch (a cascade of
/// infinitely thin flat plates along the flow), stations 0 .. n_m - 1 with blade stations f .. l.
/// The C2 wall curves are the straight walls X(u) = (m_le + u L, tan (m_le + u L)).
inline T::Grid channel(int n_m, int n_eta, int f, int l, double angle, double pitch = 1.0)
{
    T::Grid g;
    const double t = std::tan(angle);
    const double dm = 3.0 / static_cast<double>(n_m - 1);
    g.m.resize(static_cast<std::size_t>(n_m));
    for (int i = 0; i < n_m; ++i) g.m[static_cast<std::size_t>(i)] = -1.0 + dm * i;
    g.m.back() = 2.0;
    g.eta.resize(static_cast<std::size_t>(n_eta));
    for (int j = 0; j < n_eta; ++j) g.eta[static_cast<std::size_t>(j)] = static_cast<double>(j) / (n_eta - 1);
    g.node_m = T::Field(n_m, n_eta);
    g.node_y = T::Field(n_m, n_eta);
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j) {
            const double m = g.m[static_cast<std::size_t>(i)];
            g.node_m(i, j) = m;
            g.node_y(i, j) = t * m + g.eta[static_cast<std::size_t>(j)] * pitch;
        }
    g.first_blade_i = f;
    g.last_blade_i = l;
    auto& B = g.boundaries;
    B.chord = 1.0;
    B.pitch = pitch;
    B.periodic_translation = {0.0, pitch};
    B.inlet_angle = angle;
    B.exit_angle = angle;
    B.m_in = g.m.front();
    B.m_out = g.m.back();
    B.m_le = g.m[static_cast<std::size_t>(f)] - 0.5 * dm;
    B.m_te = g.m[static_cast<std::size_t>(l)] + 0.5 * dm;
    B.kappa_le = B.kappa_te = angle;
    B.lower_blade_boundary = {{B.m_le, t * B.m_le}, {B.m_te, t * B.m_te}};
    B.upper_blade_boundary = {{B.m_le, t * B.m_le + pitch}, {B.m_te, t * B.m_te + pitch}};
    const double L = B.m_te - B.m_le;
    for (int i = f; i <= l; ++i) {
        const double m = g.m[static_cast<std::size_t>(i)];
        const double u = (m - B.m_le) / L;
        g.wall_identity.stations.push_back(i);
        for (T::WallSide* w : {&g.wall_identity.lower, &g.wall_identity.upper}) {
            w->surface.push_back(w == &g.wall_identity.lower ? T::Surface::Suction : T::Surface::Pressure);
            w->parameter.push_back(u);
            w->u.push_back(u);
            w->contour.push_back(u * L);
            w->hidden.push_back(0);
        }
        g.wall_d[0].lower.push_back({m, g.node_y(i, 0)});
        g.wall_d[0].upper.push_back({m, g.node_y(i, n_eta - 1)});
        for (auto* side : {&g.wall_d[1].lower, &g.wall_d[1].upper}) side->push_back({L, t * L});
        for (auto* side : {&g.wall_d[2].lower, &g.wall_d[2].upper}) side->push_back({0.0, 0.0});
    }
    return g;
}

}  // namespace sref
