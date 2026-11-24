// Port of v1 flow.py (see Flow.hpp): inlet_state, passage_mass_flow, _speed_at_mach,
// velocity_components, station_mass_flow, _derivatives, _velocities, _check_density_of.
#include "turbolenza/Flow.hpp"

#include "turbolenza/Detail.hpp"

#include <cmath>
#include <format>
#include <limits>
#include <numbers>
#include <string>

namespace BladeGen::Turbolenza {

using detail::py_repr;
using detail::sz;

double passage_mass_flow(double density, double meridional_velocity, double inlet_pitch) noexcept
{
    return density * meridional_velocity * inlet_pitch;
}

namespace {

// v1 _speed_at_mach: W with W / a(W) = mach by bisection on [0, sqrt(2 H')).
double speed_at_mach(double mach, const PerfectGas& gas, const BladeRowReference& reference)
{
    if (mach == 0.0) return 0.0;
    auto excess = [&](double W) {
        try {
            return W - mach * static_state_at(W, gas, reference).sound_speed;
        } catch (const ThermodynamicStateError&) {
            return std::numeric_limits<double>::infinity();
        }
    };
    double low = 0.0, high = std::sqrt(2.0 * reference.relative_total_enthalpy());
    while (true) {
        const double middle = 0.5 * (low + high);
        if (!(low < middle && middle < high)) break;
        if (excess(middle) < 0.0) low = middle;
        else high = middle;
    }
    return std::abs(excess(low)) <= std::abs(excess(high)) ? low : high;
}

}  // namespace

InletState inlet_state(const FlowConditions& conditions, const PerfectGas& gas, double inlet_pitch)
{
    if (!conditions.inlet_mach || !(std::isfinite(*conditions.inlet_mach) && *conditions.inlet_mach >= 0.0))
        throw InletStateError(std::format("the inlet Mach number is {}: a finite number >= 0 is needed to derive the "
                                          "inlet state",
                                          conditions.inlet_mach ? py_repr(*conditions.inlet_mach) : "None"));
    const double mach = *conditions.inlet_mach;
    const double beta = detail::radians(conditions.inlet_angle_deg);
    if (!(std::isfinite(beta) && std::abs(beta) < 0.5 * std::numbers::pi))
        throw InletStateError(
            std::format("the inlet angle {} deg is not within (-90, 90)", py_repr(conditions.inlet_angle_deg)));
    if (!(std::isfinite(inlet_pitch) && inlet_pitch > 0.0))
        throw InletStateError(std::format("the inlet pitch {} must be a positive length", py_repr(inlet_pitch)));
    const TotalState total(conditions.inlet_absolute_total_pressure, conditions.inlet_absolute_total_temperature);
    const bool rotating = conditions.angular_speed != 0.0;
    BladeRowReference reference;
    if (!rotating) {
        reference = blade_row_reference(total, gas);
    } else {  // not in v1: the rothalpy of the rotating row at the inlet radius (Eq. 5-51)
        try {
            reference = blade_row_reference(total, gas, conditions.angular_speed, conditions.stream_surface_radius,
                                            conditions.inlet_absolute_tangential_velocity);
        } catch (const ThermodynamicStateError& error) {
            throw InletStateError(std::string("the rotating row's rothalpy: ") + error.what());
        }
    }
    const double W = speed_at_mach(mach, gas, reference);
    const ScalarState state = static_state_at(W, gas, reference);
    InletState s;
    if (rotating) {
        s.rotating = true;
        s.blade_speed = reference.blade_speed();
        s.absolute_tangential_velocity = *conditions.inlet_absolute_tangential_velocity;
        // the relative total state at the inlet: T0' = T + W^2 / (2 cp), isentropic
        const double T0r = state.temperature + 0.5 * (W * W) / gas.cp();
        s.relative_total_pressure = state.pressure * std::pow(T0r / state.temperature, gas.gamma / (gas.gamma - 1.0));
    }
    s.mach = mach;
    s.flow_angle = beta;
    s.total_pressure = total.pressure;
    s.total_temperature = total.temperature;
    s.static_temperature = state.temperature;
    s.static_pressure = state.pressure;
    s.static_density = state.density;
    s.sound_speed = state.sound_speed;
    s.velocity = W;
    s.meridional_velocity = W * std::cos(beta);
    s.tangential_velocity = W * std::sin(beta);
    s.dynamic_pressure = 0.5 * state.density * W * W;
    s.pitch = inlet_pitch;
    s.mass_flow_per_depth = passage_mass_flow(state.density, s.meridional_velocity, inlet_pitch);
    return s;
}

VelocityTriangle velocity_triangle(double rpm, double alpha1_deg, double mach_absolute, double P0, double T0,
                                   const PerfectGas& gas, double radius)
{
    // not in v1 (see Flow.hpp for the sign convention, ours)
    if (!std::isfinite(rpm)) throw InletStateError("the rotation speed rpm must be finite");
    if (!(std::isfinite(alpha1_deg) && std::abs(alpha1_deg) < 90.0))
        throw InletStateError(std::format("the absolute inlet angle {} deg is not within (-90, 90)", py_repr(alpha1_deg)));
    if (!(std::isfinite(mach_absolute) && mach_absolute > 0.0))
        throw InletStateError(std::format("the absolute inlet Mach number {} must be positive", py_repr(mach_absolute)));
    if (!(std::isfinite(P0) && P0 > 0.0 && std::isfinite(T0) && T0 > 0.0))
        throw InletStateError("the absolute inlet totals P0, T0 must be positive");
    if (!(std::isfinite(radius) && radius > 0.0))
        throw InletStateError(std::format("the inlet radius {} must be a positive length", py_repr(radius)));
    VelocityTriangle t;
    t.rpm = rpm;
    t.omega = -2.0 * std::numbers::pi * rpm / 60.0;
    t.radius = radius;
    t.alpha = detail::radians(alpha1_deg);
    t.mach_absolute = mach_absolute;
    t.total_pressure = P0;
    t.total_temperature = T0;
    const double g = gas.gamma, R = gas.gas_constant;
    t.static_temperature = T0 / (1.0 + 0.5 * (g - 1.0) * mach_absolute * mach_absolute);
    t.static_pressure = P0 * std::pow(t.static_temperature / T0, g / (g - 1.0));
    t.static_density = t.static_pressure / (R * t.static_temperature);
    t.sound_speed = std::sqrt(g * R * t.static_temperature);
    t.C = mach_absolute * t.sound_speed;
    t.C_m = t.C * std::cos(t.alpha);
    t.C_theta = t.C * std::sin(t.alpha);
    t.U = t.omega * radius;
    t.W_theta = t.C_theta - t.U;
    t.beta = std::atan(t.W_theta / t.C_m);
    t.W = std::hypot(t.C_m, t.W_theta);
    t.mach_relative = t.W / t.sound_speed;
    t.rothalpy = gas.cp() * T0 - t.U * t.C_theta;
    return t;
}

namespace {

void check_field(const Field& f, int n_m, int n_eta, const char* name)
{
    if (f.n_m != n_m || f.n_eta != n_eta)
        throw FlowRecoveryError(
            std::format("{} has shape ({}, {}), the grid needs ({}, {})", name, f.n_m, f.n_eta, n_m, n_eta));
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (!std::isfinite(f(i, j))) throw FlowRecoveryError(std::format("{} is not finite at node ({}, {})", name, i, j));
}

// v1 _check_density_of: C = mdot' S / rho and B = mdot' tan(beta) / rho to round-off (rho b on a
// stream sheet, not in v1).
void check_density_of(const CoefficientFields& co, const Field& rho, const std::vector<double>& S, const Field& tan_beta,
                      double mdot, const GridMetrics& metrics)
{
    for (int pass = 0; pass < 2; ++pass)
        for (int i = 0; i < rho.n_m; ++i)
            for (int j = 0; j < rho.n_eta; ++j) {
                const double rb = metrics.on_sheet ? rho(i, j) * metrics.thickness_at(i) : rho(i, j);
                const double off =
                    pass == 0 ? std::abs(co.C(i, j) * rb - mdot * S[sz(i)]) / (mdot * S[sz(i)])
                              : std::abs(co.B(i, j) * rb - mdot * tan_beta(i, j)) /
                                    (mdot * std::max(1.0, std::abs(tan_beta(i, j))));
                if (!(off <= kRoundOff))
                    throw FlowRecoveryError(std::format(
                        "{} fails at node ({}, {}) by {:.3g} (relative): the density, the metrics and the coefficients "
                        "do not belong together. Recover the velocities with the density and the coefficients psi was "
                        "solved with",
                        pass == 0 ? "C = mdot' S / rho" : "B = mdot' tan(beta) / rho", i, j, off));
            }
}

}  // namespace

FlowField velocity_components(const Field& psi, const Field& density, const CoefficientFields& co,
                              const GridMetrics& metrics, const BoundaryConditions& bcs)
{
    if (!metrics.grid) throw FlowRecoveryError("velocity_components takes GridMetrics with their grid");
    if (bcs.grid != metrics.grid)
        throw FlowRecoveryError("the boundary conditions were built from another grid than the metrics");
    const Grid& grid = *metrics.grid;
    const int n_m = grid.n_m(), n_eta = grid.n_eta();
    check_field(psi, n_m, n_eta, "psi");
    check_field(density, n_m, n_eta, "the density");
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (!(density(i, j) > 0.0))
                throw FlowRecoveryError(
                    std::format("the density is {} at node ({}, {}): it must be positive", py_repr(density(i, j)), i, j));
    if (co.A.n_m != n_m || co.A.n_eta != n_eta)
        throw FlowRecoveryError(std::format("the coefficient fields are ({}, {}), the grid has ({}, {}) nodes", co.A.n_m,
                                            co.A.n_eta, n_m, n_eta));
    const double length = co.length_scale();
    const double mdot = co.mass_flow_per_depth;
    std::vector<double> S(sz(n_m));
    for (int i = 0; i < n_m; ++i) {
        S[sz(i)] = length * metrics.S[sz(i)];
        if (!(std::isfinite(S[sz(i)]) && S[sz(i)] > 0.0))
            throw FlowRecoveryError("S must be finite and positive at every station");
    }
    for (double t : metrics.tan_beta.v)
        if (!std::isfinite(t))
            throw FlowRecoveryError(std::format("tan(beta) must be a finite ({}, {}) field", n_m, n_eta));
    const double dm = length * grid.dm(), deta = grid.deta();
    if (!(std::isfinite(dm) && dm > 0.0))
        throw FlowRecoveryError(std::format("dm = {}: a grid spacing must be finite and positive", py_repr(dm)));
    if (!(std::isfinite(deta) && deta > 0.0))
        throw FlowRecoveryError(std::format("deta = {}: a grid spacing must be finite and positive", py_repr(deta)));
    check_density_of(co, density, S, metrics.tan_beta, mdot, metrics);

    // v1 _derivatives
    if (bcs.n_m != n_m || bcs.n_eta != n_eta)
        throw FlowRecoveryError(std::format("psi is ({}, {}), the boundary conditions are for ({}, {})", n_m, n_eta,
                                            bcs.n_m, bcs.n_eta));
    if (n_m < 3 || n_eta < 3)
        throw FlowRecoveryError(
            std::format("{} x {} nodes: the second-order differences need three a side", n_m, n_eta));
    const int J = n_eta - 1;
    std::vector<int> blade;
    for (int i = 0; i < n_m; ++i) {
        if (bcs.lower_wall[bcs.at(i, 0)] != bcs.upper_wall[bcs.at(i, J)])
            throw FlowRecoveryError("the lower and upper walls are not on the same stations");
        if (bcs.lower_wall[bcs.at(i, 0)]) blade.push_back(i);
    }
    if (blade.empty() || blade.back() - blade.front() + 1 != static_cast<int>(blade.size()))
        throw FlowRecoveryError("the blade stations of the boundary conditions are not one run");
    FlowField out{Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta),
                  Field(n_m, n_eta)};
    Field& psi_eta = out.psi_eta;
    for (int i = 0; i < n_m; ++i)
        for (int j = 1; j < J; ++j) psi_eta(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / (2.0 * deta);
    for (int i : blade) {
        psi_eta(i, 0) = (-3.0 * psi(i, 0) + 4.0 * psi(i, 1) + -1.0 * psi(i, 2)) / (2.0 * deta);
        psi_eta(i, J) = (3.0 * psi(i, J) + -4.0 * psi(i, J - 1) + 1.0 * psi(i, J - 2)) / (2.0 * deta);
    }
    for (int i = 0; i < n_m; ++i) {
        if (!bcs.periodic_lower[bcs.at(i, 0)]) continue;
        const double jump = bcs.station_jump(i);  // 1 (v1), or the wake's jump of the equivalent inviscid flow
        const double below = bcs.periodic_ghost(i, -1, jump).value(psi);     // psi[i, -1] = psi[i, J - 1] - 1
        const double above = bcs.periodic_ghost(i, n_eta, jump).value(psi);  // psi[i, J + 1] = psi[i, 1] + 1
        psi_eta(i, 0) = (psi(i, 1) - below) / (2.0 * deta);
        psi_eta(i, J) = (above - psi(i, J - 1)) / (2.0 * deta);
    }
    std::vector<MRun> regions;
    try {
        regions = m_runs(n_m, blade.front(), blade.back());
    } catch (const MetricsError& error) {
        throw FlowRecoveryError(std::string("psi_m by m region: ") + error.what());
    }
    out.psi_m.v = d_dm_runs(psi.v, n_eta, dm, regions);
    // v1: psi_m = 0 on the walls (Dirichlet psi = 0, 1).  With transpiration the wall psi varies along
    // m and psi_m there is d(psi_wall)/dm, the differences of the blade run along the wall itself
    if (!bcs.transpiration)
        for (int i = 0; i < n_m; ++i)
            for (int j = 0; j < n_eta; ++j)
                if (bcs.blade_wall(i, j)) out.psi_m(i, j) = 0.0;  // the walls take precedence
    if (bcs.transpiration) {
        // a jump that varies along the exit (the marched wake) is a source sheet on the wake cut: psi_m
        // of the two images differs by d(jump)/dm.  Ours: both images take the mean (the sheet's normal
        // velocity split equally between its sides), so W and the density stay periodic, as the
        // coefficients require.  A constant jump (v1, the constant wake) is left as it is.
        const int first_exit = blade.back() + 1;
        bool varying = false;
        for (int i = first_exit + 1; i < n_m; ++i) varying = varying || bcs.station_jump(i) != bcs.station_jump(first_exit);
        if (varying)
            for (int i = first_exit; i < n_m; ++i) {
                const double mean = 0.5 * (out.psi_m(i, 0) + out.psi_m(i, J));
                out.psi_m(i, 0) = mean;
                out.psi_m(i, J) = mean;
            }
    }

    // v1 _velocities: Eqs. 5-18, 5-19 with mu = mdot' / rho
    const int n = n_m * n_eta;
    TURBOLENZA_PARALLEL_FOR
    for (int k = 0; k < n; ++k) {
        const std::size_t q = static_cast<std::size_t>(k);
        const double s = S[sz(k / n_eta)];
        // mu = mdot' / (rho b) on a stream sheet (Eqs. 5-18, 5-19), v1's mdot' / rho otherwise
        const double mu = metrics.on_sheet ? mdot / (density.v[q] * metrics.thickness_at(k / n_eta)) : mdot / density.v[q];
        const double wm = mu * out.psi_eta.v[q] / s;
        const double wt = mu * (metrics.tan_beta.v[q] / s * out.psi_eta.v[q] - out.psi_m.v[q]);
        const double w2 = wm * wm + wt * wt;
        out.W_m.v[q] = wm;
        out.W_theta.v[q] = wt;
        out.W2.v[q] = w2;
        out.W.v[q] = std::sqrt(w2);
    }
    const std::pair<const char*, const Field*> all[] = {{"psi_m", &out.psi_m}, {"psi_eta", &out.psi_eta},
                                                        {"W_m", &out.W_m},     {"W_theta", &out.W_theta},
                                                        {"W2", &out.W2},       {"W", &out.W}};
    for (const auto& [name, f] : all)
        for (double v : f->v)
            if (!std::isfinite(v)) throw FlowRecoveryError(std::format("the recovered {} is not finite", name));
    return out;
}

std::vector<double> station_mass_flow(const FlowField& flow, const Field& density, const CoefficientFields& co,
                                      const GridMetrics& metrics)
{
    const Grid& grid = *metrics.grid;
    const int n_m = grid.n_m(), n_eta = grid.n_eta();
    check_field(density, n_m, n_eta, "the density");
    if (flow.W_m.n_m != n_m || flow.W_m.n_eta != n_eta)
        throw FlowRecoveryError(std::format("the flow field is ({}, {}), the grid has ({}, {}) nodes", flow.W_m.n_m,
                                            flow.W_m.n_eta, n_m, n_eta));
    const double deta = grid.deta();
    std::vector<double> out(sz(n_m));
    std::vector<double> row(sz(n_eta));
    for (int i = 0; i < n_m; ++i) {
        const double S = co.length_scale() * metrics.S[sz(i)];
        for (int j = 0; j < n_eta; ++j)  // rho W_m S b on a stream sheet (Eq. 5-17), v1's rho W_m S otherwise
            row[sz(j)] = metrics.on_sheet ? density(i, j) * flow.W_m(i, j) * S * metrics.thickness_at(i)
                                          : density(i, j) * flow.W_m(i, j) * S;
        out[sz(i)] = deta * (detail::pairwise_sum(row) - 0.5 * (row.front() + row.back()));
    }
    return out;
}

}  // namespace BladeGen::Turbolenza
