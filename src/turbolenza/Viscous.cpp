// The boundary layer on a solved flow and its coupling to the potential flow (see Viscous.hpp).
#include "turbolenza/Viscous.hpp"

#include "turbolenza/Curve.hpp"
#include "turbolenza/Detail.hpp"
#include "turbolenza/Postprocess.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <utility>

namespace BladeGen::Turbolenza {

using detail::nan;
using detail::sz;

namespace {

// 8-point Gauss-Legendre on [-1, 1] (Abramowitz and Stegun, Table 25.4)
constexpr double kGaussX[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
constexpr double kGaussW[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
constexpr int kGaussParts = 8;  // parts of every knot span between two nodes

double speed_of(const ParametricBladeCurve& c, double u)
{
    const Point d = c.derivative(u, 1);
    return std::hypot(d.m, d.y);
}

double gauss(const ParametricBladeCurve& c, double a, double b)
{
    const double h = 0.5 * (b - a), mid = 0.5 * (a + b);
    double s = 0.0;
    for (int k = 0; k < 4; ++k) s += kGaussW[k] * (speed_of(c, mid - h * kGaussX[k]) + speed_of(c, mid + h * kGaussX[k]));
    return h * s;
}

// int_a^b |X'(u)| du, split at the curve's interior knots (|X'| is smooth inside a span)
double arc(const ParametricBladeCurve& c, double a, double b)
{
    std::vector<double> cuts{a};
    for (double t : c.interior_knots())
        if (a < t && t < b) cuts.push_back(t);
    cuts.push_back(b);
    double s = 0.0;
    for (std::size_t k = 0; k + 1 < cuts.size(); ++k)
        for (int q = 0; q < kGaussParts; ++q) {
            const double lo = cuts[k] + (cuts[k + 1] - cuts[k]) * q / kGaussParts;
            const double hi = cuts[k] + (cuts[k + 1] - cuts[k]) * (q + 1) / kGaussParts;
            s += gauss(c, lo, hi);
        }
    return s;
}

const char* side_name(int j) { return j == 0 ? "lower (eta = 0)" : "upper (eta = 1)"; }

double x_c_at(const std::vector<double>& x_c, int index) { return index < 0 ? nan() : x_c[sz(index)]; }

}  // namespace

std::vector<double> blade_station_fraction(const Grid& g)
{
    const int first = g.first_blade_i, last = g.last_blade_i;
    if (!(0 <= first && first < last && last < static_cast<int>(g.m.size())))
        throw BoundaryLayerError(std::format("the grid's blade stations {} .. {} are not two or more stations of its {}", first,
                                             last, g.m.size()));
    const double m0 = g.m[sz(first)], m1 = g.m[sz(last)];
    std::vector<double> f;
    for (int i = first; i <= last; ++i) {
        if (i > first && !(g.m[sz(i)] > g.m[sz(i - 1)]))
            throw BoundaryLayerError(std::format("the grid's m does not increase at blade station {}", i));
        f.push_back((g.m[sz(i)] - m0) / (m1 - m0));
    }
    f.back() = 1.0;  // exactly, whatever the rounding of the division
    return f;
}

std::vector<double> wall_arc_length(const Grid& g, int j)
{
    if (!g.boundaries.blade_curves) throw BoundaryLayerError("the grid has no blade surface curves");
    if (g.touching())
        throw BoundaryLayerError("the boundary layer needs the first and last blade nodes inside the passage (Aungier's "
                                 "layout): the touching edge-node experiment puts them on the nose tip and the tail end");
    const bool lower = j == 0;
    if (!lower && j != g.n_eta() - 1) throw BoundaryLayerError(std::format("j = {} is not a wall", j));
    const WallSide& side = lower ? g.wall_identity.lower : g.wall_identity.upper;
    const std::vector<Point>& nodes = lower ? g.wall_d[0].lower : g.wall_d[0].upper;
    const SurfaceCurves& curves = *g.boundaries.blade_curves;
    const Point shift = lower ? Point{0.0, 0.0} : g.boundaries.periodic_translation;
    std::vector<double> x(side.u.size());
    double u_prev = 0.0, s = 0.0;
    for (std::size_t b = 0; b < side.u.size(); ++b) {
        if (side.surface[b] != side.surface[0])
            throw BoundaryLayerError(std::format("the {} wall changes its surface along the blade", side_name(j)));
        const ParametricBladeCurve& c = side.surface[b] == Surface::Suction ? curves.suction : curves.pressure;
        const double u = side.u[b];
        // the wall node must be its curve's point at u (the topology sources agree)
        const Point X = c.evaluate(u);
        const double off = std::hypot(X.m + shift.m - nodes[b].m, X.y + shift.y - nodes[b].y);
        if (!(off <= 1e-9 * g.chord()))
            throw BoundaryLayerError(std::format("the {} wall node of station {} is {:.3g} chords off its surface curve at "
                                                 "u = {}",
                                                 side_name(j), g.first_blade_i + static_cast<int>(b), off / g.chord(), u));
        if (b > 0 && !(u > u_prev))
            throw BoundaryLayerError(std::format("the {} wall's curve parameter does not increase at station {}",
                                                 side_name(j), g.first_blade_i + static_cast<int>(b)));
        s += arc(c, u_prev, u);
        x[b] = s;
        u_prev = u;
    }
    return x;
}

namespace {

WallLayer wall_layer(const Solution& sol, int j, const std::vector<double>& x_c, const PerfectGas& gas, double Tt,
                     const ViscousSettings& settings)
{
    const Grid& g = *sol.grid;
    WallLayer w;
    w.j = j;
    w.surface = (j == 0 ? g.wall_identity.lower : g.wall_identity.upper).surface.front();
    w.x_c = x_c;
    w.f = blade_station_fraction(g);
    const std::vector<double> s = wall_arc_length(g, j);
    EdgeConditions& e = w.edge;
    e.cp = gas.cp();
    e.total_temperature = Tt;
    for (int i = g.first_blade_i; i <= g.last_blade_i; ++i) {
        const std::size_t b = sz(i - g.first_blade_i);
        w.stations.push_back(i);
        e.x.push_back(sol.c_ref * s[b]);
        e.W.push_back(sol.flow.W(i, j));
        e.rho.push_back(sol.static_state.density(i, j));
        e.T.push_back(sol.static_state.temperature(i, j));
        e.mu.push_back(sutherland_viscosity(e.T.back(), settings.viscosity));
    }
    w.profile = march_surface(e, settings.layer);
    return w;
}

// the marched wake along j = 0 from the last blade station (see Viscous.hpp)
WakeLayer wake_layer(const Solution& sol, const WallLayer& lo, const WallLayer& up, const PerfectGas& gas, double Tt,
                     const ViscousSettings& settings)
{
    const Grid& g = *sol.grid;
    const int l = g.last_blade_i, J = g.n_eta() - 1, n_m = g.n_m();
    WakeLayer w;
    w.marched = true;
    EdgeConditions& e = w.edge;
    e.cp = gas.cp();
    e.total_temperature = Tt;
    double x = 0.0;
    for (int i = l; i < n_m; ++i) {
        w.stations.push_back(i);
        w.m_c.push_back((g.m[sz(i)] - g.m_te()) / g.chord());
        if (i > l) x += sol.c_ref * std::hypot(g.node_m(i, 0) - g.node_m(i - 1, 0), g.node_y(i, 0) - g.node_y(i - 1, 0));
        e.x.push_back(x);
        const double W = i == l ? 0.5 * (sol.flow.W(i, 0) + sol.flow.W(i, J)) : sol.flow.W(i, 0);
        const double rho = i == l ? 0.5 * (sol.static_state.density(i, 0) + sol.static_state.density(i, J))
                                  : sol.static_state.density(i, 0);
        const double T = i == l ? 0.5 * (sol.static_state.temperature(i, 0) + sol.static_state.temperature(i, J))
                                : sol.static_state.temperature(i, 0);
        e.W.push_back(W);
        e.rho.push_back(rho);
        e.T.push_back(T);
        e.mu.push_back(sutherland_viscosity(T, settings.viscosity));
    }
    const BoundaryLayerProfile& a = lo.profile;
    const BoundaryLayerProfile& b = up.profile;
    const double theta = a.theta.back() + b.theta.back();
    const double q = (a.delta.back() - a.delta_star.back()) + (b.delta.back() - b.delta_star.back());
    const double ds = a.delta_star.back() + b.delta_star.back();
    w.profile = march_wake(e, theta, q, ds, settings.layer);
    return w;
}

}  // namespace

ViscousAnalysis analyze_boundary_layers(const Solution& sol, const PerfectGas& gas, double Tt, const ViscousSettings& settings)
{
    if (!sol.grid) throw BoundaryLayerError("the solution has no grid");
    const Grid& g = *sol.grid;
    const BladeSurfaces surfaces = blade_surfaces(sol);
    const int J = g.n_eta() - 1;
    const auto& xc_lower = surfaces.ss_j == 0 ? surfaces.x_c_ss : surfaces.x_c_ps;
    const auto& xc_upper = surfaces.ss_j == 0 ? surfaces.x_c_ps : surfaces.x_c_ss;
    ViscousAnalysis a;
    a.lower = wall_layer(sol, 0, xc_lower, gas, Tt, settings);
    a.upper = wall_layer(sol, J, xc_upper, gas, Tt, settings);
    const int n_m = g.n_m(), f = g.first_blade_i, l = g.last_blade_i;
    MassDeficit& d = a.deficit;
    d.lower.assign(sz(n_m), 0.0);
    d.upper.assign(sz(n_m), 0.0);
    d.wake.assign(sz(n_m), 0.0);
    for (int i = f; i <= l; ++i) {
        d.lower[sz(i)] = a.lower.profile.mass_deficit[sz(i - f)];
        d.upper[sz(i)] = a.upper.profile.mass_deficit[sz(i - f)];
    }
    const double te = d.lower[sz(l)] + d.upper[sz(l)];
    if (settings.wake == WakeModel::Constant) {
        for (int i = l + 1; i < n_m; ++i) d.wake[sz(i)] = te;
        a.wake.marched = false;
        for (int i = l; i < n_m; ++i) {
            a.wake.stations.push_back(i);
            a.wake.m_c.push_back((g.m[sz(i)] - g.m_te()) / g.chord());
        }
    } else {
        a.wake = wake_layer(sol, a.lower, a.upper, gas, Tt, settings);
        for (int i = l + 1; i < n_m - 1; ++i) d.wake[sz(i)] = a.wake.profile.mass_deficit[sz(i - l)];
        // the far exit: what its psi_m = 0 relation implies (its psi is eliminated by that relation)
        d.wake[sz(n_m - 1)] = 4.0 / 3.0 * d.wake[sz(n_m - 2)] - 1.0 / 3.0 * d.wake[sz(n_m - 3)];
    }
    a.loss = profile_loss(a.lower.profile.theta.back() + a.upper.profile.theta.back(),
                          a.lower.profile.delta_star.back() + a.upper.profile.delta_star.back(), sol.c_ref * g.pitch(),
                          g.inlet_angle(), g.exit_angle());
    return a;
}

Transpiration transpiration_of(const Grid& g, const MassDeficit& d, double mdot)
{
    if (!(std::isfinite(mdot) && mdot > 0.0)) throw BoundaryLayerError("the transpiration needs mdot' > 0");
    const int n_m = g.n_m(), f = g.first_blade_i, l = g.last_blade_i;
    if (d.lower.size() != sz(n_m) || d.upper.size() != sz(n_m) || d.wake.size() != sz(n_m))
        throw BoundaryLayerError("the mass deficit does not have one value per station");
    Transpiration t;
    t.lower_wall.assign(sz(n_m), kPsiLowerWall);
    t.upper_wall.assign(sz(n_m), kPsiUpperWall);
    t.jump.assign(sz(n_m), kPsiPeriodicJump);
    for (int i = f; i <= l; ++i) {
        t.lower_wall[sz(i)] = -d.lower[sz(i)] / mdot;
        t.upper_wall[sz(i)] = 1.0 + d.upper[sz(i)] / mdot;
    }
    for (int i = l + 1; i < n_m; ++i) t.jump[sz(i)] = 1.0 + d.wake[sz(i)] / mdot;
    return t;
}

double deficit_change(const Grid& g, const MassDeficit& a, const MassDeficit& b)
{
    double num = 0.0, den = 0.0;
    for (int i = g.first_blade_i; i <= g.last_blade_i; ++i) {
        num = std::max({num, std::abs(a.lower[sz(i)] - b.lower[sz(i)]), std::abs(a.upper[sz(i)] - b.upper[sz(i)])});
        den = std::max({den, std::abs(a.lower[sz(i)]), std::abs(a.upper[sz(i)])});
    }
    return den > 0.0 ? num / den : (num > 0.0 ? detail::nan() : 0.0);
}

namespace {

BlIteration summary_of(int k, const Solution& sol, const ViscousAnalysis& a, double change, int passes)
{
    BlIteration h;
    h.k = k;
    h.inviscid_passes = passes;
    h.beta_out = sol.exit_angle;
    h.max_change_m = change;
    const WallLayer& ss = a.suction();
    const WallLayer& ps = a.pressure();
    h.transition_x_c_ss = x_c_at(ss.x_c, ss.profile.transition_index);
    h.transition_x_c_ps = x_c_at(ps.x_c, ps.profile.transition_index);
    h.separation_x_c_ss = x_c_at(ss.x_c, ss.profile.separation_index);
    h.separation_x_c_ps = x_c_at(ps.x_c, ps.profile.separation_index);
    h.laminar_separation_x_c_ss = x_c_at(ss.x_c, ss.profile.laminar_separation_index);
    h.laminar_separation_x_c_ps = x_c_at(ps.x_c, ps.profile.laminar_separation_index);
    h.loss = a.loss.omega;
    return h;
}

MassDeficit zero_deficit(int n_m)
{
    return {std::vector<double>(sz(n_m), 0.0), std::vector<double>(sz(n_m), 0.0), std::vector<double>(sz(n_m), 0.0)};
}

}  // namespace

ViscousSolve solve_viscous(std::shared_ptr<const Grid> grid, const SolverConfig& config, double c_ref,
                           std::optional<double> mdot_override, const SolveOptions& base, const ViscousSettings& settings,
                           bool coupled)
{
    ViscousSolve out;
    if (!grid) {
        out.error = "SolverInputError: the solver takes a Grid, got None";
        return out;
    }
    if (!(settings.relaxation > 0.0 && settings.relaxation <= 1.0) || !(settings.tolerance > 0.0) ||
        settings.max_iterations < 0) {
        out.error = std::format("ValueError: the coupling needs 0 < omega_BL <= 1 (got {}), a tolerance > 0 (got {}) and "
                                "max_bl_iterations >= 0 (got {})",
                                detail::py_repr(settings.relaxation), detail::py_repr(settings.tolerance),
                                settings.max_iterations);
        return out;
    }
    const BladeRowReference reference = blade_row_reference(
        TotalState(config.flow.inlet_absolute_total_pressure, config.flow.inlet_absolute_total_temperature), config.gas,
        config.flow.angular_speed, config.flow.stream_surface_radius, config.flow.inlet_absolute_tangential_velocity);
    const double Tt = reference.relative_total_enthalpy() / config.gas.cp();  // T't
    const int limit = settings.coupling == BlCoupling::OneWay ? 0 : settings.max_iterations;
    MassDeficit used = zero_deficit(grid->n_m());
    std::shared_ptr<const Transpiration> transpiration;
    std::shared_ptr<const Solution> previous;
    for (int k = 0;; ++k) {
        SolveOptions options = base;
        SolveFailure failure;
        options.failure = &failure;
        options.bl_iteration = k;
        options.transpiration = transpiration;
        if (previous) options.initial_density = previous->density;
        const std::shared_ptr<const Grid> start = previous ? previous->grid : grid;
        std::shared_ptr<const Solution> sol;
        // a failed inviscid solve: its last completed pass is reported (else the previous iteration's
        // solution); in an iteration k >= 1 the coupling failed.  The last layers computed (iteration
        // k - 1) are kept with the solution they belong to (analysis_solution), so they can still be shown.
        auto fail = [&](const std::exception& e, std::shared_ptr<const Solution> partial, int failed_pass) {
            const std::string type = error_type_name(e);
            out.error = (k == 0 || type == "Cancelled")
                            ? type + ": " + e.what()
                            : std::format("BoundaryLayerCouplingError: the coupled iteration failed in BL iteration {} "
                                          "({}: {})",
                                          k, type, e.what());
            if (partial) out.solution = std::move(partial);
            if (failed_pass >= 0) out.failed_pass = failed_pass;
        };
        try {
            sol = std::make_shared<const Solution>(coupled ? solve_coupled(start, config, c_ref, mdot_override, options)
                                                           : solve_at_exit_angle(start, config, c_ref, mdot_override, options));
        } catch (const DensityIterationError& e) {
            fail(e, e.solution, -1);
            return out;
        } catch (const std::exception& e) {
            fail(e, failure.solution, failure.failed_pass);
            return out;
        }
        if (k == 0) out.inviscid = sol;
        ViscousAnalysis a;
        try {
            a = analyze_boundary_layers(*sol, config.gas, Tt, settings);
        } catch (const std::exception& e) {
            out.solution = sol;   // the layers of this solution failed; analysis keeps iteration k - 1's, if any
            out.used = used;
            const std::string what = error_type_name(e) + ": " + e.what();
            out.error = k == 0 ? what
                               : std::format("BoundaryLayerCouplingError: the coupled iteration failed in BL iteration {} "
                                             "({})",
                                             k, what);
            return out;
        }
        const double change = deficit_change(*grid, a.deficit, used);
        out.history.push_back(summary_of(k, *sol, a, change, static_cast<int>(sol->passes.size())));
        out.solution = sol;
        out.used = used;
        out.analysis = a;
        out.analysis_solution = sol;
        out.analysis_iteration = k;
        if (limit == 0 || change <= settings.tolerance) {
            out.converged = true;
            return out;
        }
        if (k >= limit) {
            out.error = std::format("BoundaryLayerNotConvergedError: the coupled boundary-layer iteration did not converge in "
                                    "{} iteration(s) (max relative change of m {:.3e} against {:g})",
                                    limit, change, settings.tolerance);
            return out;
        }
        if (base.cancel && base.cancel->load()) {
            out.error = std::format("Cancelled: the run was cancelled after BL iteration {}", k);
            return out;
        }
        const double w = settings.relaxation;
        for (auto [u, n] : {std::pair{&used.lower, &a.deficit.lower}, std::pair{&used.upper, &a.deficit.upper},
                            std::pair{&used.wake, &a.deficit.wake}})
            for (std::size_t q = 0; q < u->size(); ++q) (*u)[q] += w * ((*n)[q] - (*u)[q]);
        transpiration = std::make_shared<const Transpiration>(transpiration_of(*grid, used, sol->mass_flow_per_depth));
        previous = sol;
    }
}

}  // namespace BladeGen::Turbolenza
