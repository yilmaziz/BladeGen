// The density- and Kutta-coupled solve on a general mesh (see GeneralSolver.hpp).  Not a port of v1;
// the pass order, relaxation, stop tests and the Kutta correction are the Aungier solver's (Solver.cpp,
// Kutta.cpp), reused through Eos, Convergence, Kutta and LinearSolver.
#include "turbolenza/GeneralSolver.hpp"

#include "turbolenza/Detail.hpp"
#include "turbolenza/Eos.hpp"
#include "turbolenza/Flow.hpp"
#include "turbolenza/Kutta.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>
#include <utility>

namespace BladeGen::Turbolenza {

using detail::nan;

namespace {

std::size_t z(int k) noexcept { return static_cast<std::size_t>(k); }

Field as_field(const std::vector<double>& v)
{
    Field f(static_cast<int>(v.size()), 1);
    f.v = v;
    return f;
}

}  // namespace

double theta_angular_speed(double angular_speed, RotationSense sense) noexcept
{
    return sense == RotationSense::AlongTheta ? angular_speed : -angular_speed;
}

GeneralTrailingEdge general_trailing_edge(const GeneralOperator& op)
{
    const GeneralMesh& mesh = op.mesh();
    const EdgeTag* ss = nullptr;
    const EdgeTag* ps = nullptr;
    int n_ss = 0, n_ps = 0;
    for (const EdgeTag& t : mesh.tags) {
        if (t.edge != EdgeKind::Trailing) continue;
        if (t.surface == Surface::Suction) {
            ss = &t;
            ++n_ss;
        } else {
            ps = &t;
            ++n_ps;
        }
    }
    if (n_ss != 1 || n_ps != 1)
        throw GeneralOperatorError(std::format("the Kutta condition needs exactly one suction-side and one pressure-side "
                                               "trailing-edge tag; the mesh has {} and {}",
                                               n_ss, n_ps));
    GeneralTrailingEdge te;
    te.suction = op.node(ss->block, ss->i, ss->j);
    te.pressure = op.node(ps->block, ps->i, ps->j);
    if (te.suction == te.pressure)
        throw GeneralOperatorError("the suction and pressure trailing-edge tags name one node: Eq. 5-34 compares two nodes");
    const double pitch = mesh.pitch;
    const double th_ss = op.y(te.suction) - static_cast<double>(ss->blade) * pitch;
    const double th_ps = op.y(te.pressure) - static_cast<double>(ps->blade) * pitch;
    te.plus = th_ss > th_ps ? te.suction : te.pressure;
    te.minus = th_ss > th_ps ? te.pressure : te.suction;
    const double th_plus = std::max(th_ss, th_ps), th_minus = std::min(th_ss, th_ps);
    const double r_te = 0.5 * (op.r(te.suction) + op.r(te.pressure));
    te.passage_width = r_te * (th_minus + pitch - th_plus);
    if (!(te.passage_width > 0.0))
        throw GeneralOperatorError(std::format("the passage width at the trailing edge is {:.3e} <= 0", te.passage_width));
    auto along = [&](int g, double& dm, double& cosb) {
        const auto l = op.location(g);
        const MeshBlock& B = mesh.blocks[z(l[0])];
        const int axis = (l[1] == 0 || l[1] == B.ni - 1) && !(l[2] == 0 || l[2] == B.nj - 1) ? 1 : 0;
        double best = nan();
        for (int s : {-1, 1}) {
            const auto q = op.step(g, axis, s);
            if (!q || !op.wall(q->node)) continue;
            const double d = op.m(g) - op.m(q->node);
            if (d > 0.0 && !(best <= d)) best = d;
        }
        if (!(best > 0.0))
            throw GeneralOperatorError(std::format("the trailing-edge node {} has no upstream neighbour on its wall", g));
        dm = best;
        const NodalMetric& mt = op.metric(g);
        const double xt = axis == 0 ? mt.x_z : mt.x_e, yt = axis == 0 ? mt.y_z : mt.y_e;
        cosb = std::abs(xt) / std::hypot(xt, yt);
    };
    double dm_p = 0.0, dm_m = 0.0;
    along(te.plus, dm_p, te.cos_plus);
    along(te.minus, dm_m, te.cos_minus);
    te.dm = 0.5 * (dm_p + dm_m);
    return te;
}

GeneralKutta general_kutta(const GeneralOperator& op, const GeneralTrailingEdge& te, const std::vector<double>& W_m,
                           const std::vector<double>& W)
{
    GeneralKutta k;
    k.W_ss = W[z(te.suction)];
    k.W_ps = W[z(te.pressure)];
    if (!std::isfinite(k.W_ss) || !std::isfinite(k.W_ps)) throw KuttaError("a trailing-edge speed is not finite");
    k.residual = k.W_ss - k.W_ps;
    const double mean = 0.5 * (std::abs(k.W_ss) + std::abs(k.W_ps));
    k.normalized = mean > 0.0 ? std::abs(k.residual) / mean : 0.0;
    k.wall_difference = W[z(te.plus)] - W[z(te.minus)];
    double integral = 0.0, extent = 0.0;
    for (std::size_t s = 0; s < op.mesh().segments.size(); ++s) {
        if (op.mesh().segments[s].kind != SegmentKind::Outlet) continue;
        // each segment in the direction of increasing theta (its node order may run either way)
        const std::vector<int> nodes = op.segment_nodes(static_cast<int>(s));
        const double sense = op.y(nodes.back()) >= op.y(nodes.front()) ? 1.0 : -1.0;
        for (std::size_t q = 0; q + 1 < nodes.size(); ++q) {
            const double dth = sense * (op.y(nodes[q + 1]) - op.y(nodes[q]));
            integral += 0.5 * (W_m[z(nodes[q])] + W_m[z(nodes[q + 1])]) * dth;
            extent += dth;
        }
    }
    if (!(extent > 0.0)) throw KuttaError("the mesh has no outlet extent in theta to carry the exit-angle correction");
    k.outlet_mean_W_m = integral / extent;
    if (!(k.outlet_mean_W_m > 0.0))
        throw KuttaError(std::format("the mean W_m at the outlet is {}: no flow leaves, the exit angle cannot carry an "
                                     "angular-momentum correction (Eq. 13-59)",
                                     k.outlet_mean_W_m));
    const double delta_c = 2.0 * k.wall_difference * te.dm / (te.passage_width * (te.cos_plus + te.cos_minus));  // Eq. 13-57
    k.delta_tan = delta_c / k.outlet_mean_W_m;                                                                   // Eq. 13-59
    return k;
}

double general_boundary_flux(const GeneralOperator& op, SegmentKind kind, const std::vector<double>& rho,
                             const std::vector<double>& W_m, const std::vector<double>& W_theta)
{
    const std::size_t n = static_cast<std::size_t>(op.n_nodes());
    if (rho.size() != n || W_m.size() != n || W_theta.size() != n)
        throw GeneralOperatorError("general_boundary_flux: the fields must have one value per node");
    double flux = 0.0;
    for (std::size_t s = 0; s < op.mesh().segments.size(); ++s) {
        if (op.mesh().segments[s].kind != kind) continue;
        const std::vector<int> nodes = op.segment_nodes(static_cast<int>(s));
        const double sense = op.y(nodes.back()) >= op.y(nodes.front()) ? 1.0 : -1.0;  // toward +theta
        auto f = [&](int g, const std::vector<double>& w) { return rho[z(g)] * op.b(g) * op.r(g) * w[z(g)]; };
        for (std::size_t q = 0; q + 1 < nodes.size(); ++q) {
            const int a = nodes[q], b = nodes[q + 1];
            flux += sense * (0.5 * (f(a, W_m) + f(b, W_m)) * (op.y(b) - op.y(a)) -
                             0.5 * (f(a, W_theta) + f(b, W_theta)) * (op.x(b) - op.x(a)));
        }
    }
    return flux;
}

namespace {

struct Setup
{
    std::shared_ptr<const GeneralOperator> op;
    PerfectGas gas;
    BladeRowReference reference;
    double omega = 0.0, mdot = 0.0, relaxation = 0.0, inlet_tan = 0.0;
    std::vector<double> source;
    std::vector<int> inlet, outlet;
    std::optional<GeneralTrailingEdge> te;
    bool stop_at_sonic = true;
    LinearSolverKind kind = LinearSolverKind::SparseLU;
    int threads = 1;
};

struct PassState
{
    std::vector<double> rho, psi, W_m, W_theta, W, mach, target, next, pressure, temperature;
    double linear_residual = 0.0;
    double exit_angle = 0.0, exit_angle_next = 0.0;
    std::optional<GeneralKutta> kutta;
    double inlet_error = 0.0, outlet_error = 0.0;
};

ScalarState state_at(const Setup& c, int g, double W)
{
    if (c.omega == 0.0) return static_state_at(W, c.gas, c.reference);
    BladeRowReference ref = c.reference;
    ref.radius = c.op->r(g);
    return static_state_at(W, c.gas, ref);
}

// |flux / mdot' - 1| through the segments of one kind.
double boundary_flow_error(const Setup& c, SegmentKind kind, const PassState& p)
{
    return std::abs(general_boundary_flux(*c.op, kind, p.rho, p.W_m, p.W_theta) / c.mdot - 1.0);
}

PassState density_pass(const Setup& c, const std::vector<double>& rho, double exit_angle)
{
    const GeneralOperator& op = *c.op;
    const int n = op.n_nodes();
    PassState p;
    p.rho = rho;
    p.exit_angle = exit_angle;
    OperatorFields fields;
    fields.coefficient.resize(z(n));
    for (int g = 0; g < n; ++g) fields.coefficient[z(g)] = c.mdot / (rho[z(g)] * op.b(g));
    fields.source = c.source;
    fields.flow_angle_tan.assign(z(n), 0.0);
    for (int g : c.inlet) fields.flow_angle_tan[z(g)] = c.inlet_tan;
    for (int g : c.outlet) fields.flow_angle_tan[z(g)] = std::tan(exit_angle);
    const GeneralSystem sys = op.assemble(fields);
    const std::vector<double> u = solve(sys.matrix, sys.rhs, c.kind, c.threads);
    p.linear_residual = residual_norm(sys.matrix, u, sys.rhs);
    p.psi = op.expand(u, fields);
    const StreamGradient G = op.gradient(p.psi);
    p.W_m.resize(z(n));
    p.W_theta.resize(z(n));
    p.W.resize(z(n));
    p.mach.resize(z(n));
    p.target.resize(z(n));
    p.pressure.resize(z(n));
    p.temperature.resize(z(n));
    for (int g = 0; g < n; ++g) {
        const double f = c.mdot / (rho[z(g)] * op.b(g) * op.r(g));
        p.W_m[z(g)] = f * G.psi_y[z(g)];
        p.W_theta[z(g)] = -f * G.psi_x[z(g)];
        p.W[z(g)] = std::hypot(p.W_m[z(g)], p.W_theta[z(g)]);
        const ScalarState s = state_at(c, g, p.W[z(g)]);
        p.target[z(g)] = s.density;
        p.pressure[z(g)] = s.pressure;
        p.temperature[z(g)] = s.temperature;
        p.mach[z(g)] = p.W[z(g)] / s.sound_speed;
    }
    p.next = relax_density(as_field(rho), as_field(p.target), c.relaxation).v;
    p.inlet_error = boundary_flow_error(c, SegmentKind::Inlet, p);
    p.outlet_error = boundary_flow_error(c, SegmentKind::Outlet, p);
    p.exit_angle_next = exit_angle;
    return p;
}

std::shared_ptr<GeneralSolution> make_solution(const Setup& c, const PassState& p, const History& h,
                                               const std::vector<GeneralPass>& passes, bool converged, int iterations,
                                               double inlet_angle)
{
    auto s = std::make_shared<GeneralSolution>();
    s->op = c.op;
    s->psi = p.psi;
    s->density = p.rho;
    s->density_target = p.target;
    s->W_m = p.W_m;
    s->W_theta = p.W_theta;
    s->W = p.W;
    s->mach = p.mach;
    s->pressure = p.pressure;
    s->temperature = p.temperature;
    s->mass_flow_per_depth = c.mdot;
    s->inlet_angle = inlet_angle;
    s->exit_angle = p.exit_angle;
    s->kutta = p.kutta;
    s->history = h;
    s->passes = passes;
    s->converged = converged;
    s->iterations = iterations;
    return s;
}

}  // namespace

GeneralSolution solve_general(const GeneralMesh& mesh, const GeneralConfig& config)
{
    return solve_general(std::make_shared<const GeneralOperator>(mesh), config);
}

GeneralSolution solve_general(std::shared_ptr<const GeneralOperator> op, const GeneralConfig& config)
{
    if (!op) throw GeneralOperatorError("solve_general needs an operator");
    const detail::ThreadScope threads(config.threads);
    Setup c;
    c.op = op;
    c.gas = config.gas;
    c.kind = config.linear_solver;
    c.threads = config.threads;
    c.stop_at_sonic = config.iteration.stop_at_sonic;
    c.relaxation = config.iteration.density_relaxation;
    c.inlet = op->segment_nodes(SegmentKind::Inlet);
    c.outlet = op->segment_nodes(SegmentKind::Outlet);
    if (c.inlet.empty() || c.outlet.empty())
        throw GeneralOperatorError("solve_general needs an inlet and an outlet segment (the flow angles)");
    const double inlet_angle = detail::radians(config.flow.inlet_angle_deg);
    const double exit0 = detail::radians(config.flow.exit_angle_guess_deg);
    for (double a : {inlet_angle, exit0})
        if (!(std::isfinite(a) && std::abs(a) < 0.5 * std::numbers::pi))
            throw GeneralOperatorError(std::format("the flow angle {} rad is not within (-pi/2, pi/2)", a));
    c.inlet_tan = std::tan(inlet_angle);
    c.omega = theta_angular_speed(config.flow.angular_speed, config.rotation);
    const TotalState total(config.flow.inlet_absolute_total_pressure, config.flow.inlet_absolute_total_temperature);
    c.reference = blade_row_reference(total, config.gas, c.omega, config.flow.stream_surface_radius,
                                      config.flow.inlet_absolute_tangential_velocity);
    if (config.mass_flow_per_depth) {
        c.mdot = *config.mass_flow_per_depth;
    } else {
        if (c.omega != 0.0)
            throw GeneralOperatorError("a rotating row needs mass_flow_per_depth (the inlet closure is stationary only)");
        const int g0 = c.inlet.front();
        for (int g : c.inlet)
            if (std::abs(op->r(g) - op->r(g0)) > 1e-12 * op->r(g0) || std::abs(op->b(g) - op->b(g0)) > 1e-12 * op->b(g0))
                throw GeneralOperatorError("the inlet must lie at one radius and thickness to close the inlet state");
        const InletState in = inlet_state(config.flow, config.gas, op->r(g0) * op->mesh().pitch * op->b(g0));
        c.mdot = in.mass_flow_per_depth;
    }
    if (!(std::isfinite(c.mdot) && c.mdot > 0.0))
        throw GeneralOperatorError(std::format("the mass flow per depth is {}: it must be positive", c.mdot));
    c.source = op->rotation_source(c.omega);
    if (config.kutta) c.te = general_trailing_edge(*op);
    const IterationSettings& it = config.iteration;
    if (it.max_density_iterations < 1) throw GeneralOperatorError("max_density_iterations must be at least 1");
    if (config.kutta && it.max_kutta_iterations < 1) throw GeneralOperatorError("max_kutta_iterations must be at least 1");
    const int limit = config.kutta ? std::min(it.max_density_iterations, it.max_kutta_iterations) : it.max_density_iterations;

    std::vector<double> rho(z(op->n_nodes()));
    for (int g = 0; g < op->n_nodes(); ++g) rho[z(g)] = state_at(c, g, 0.0).density;  // Aungier's first guess

    History history;
    std::vector<GeneralPass> passes;
    std::optional<PassState> previous;
    double beta = exit0;
    auto evaluate = [&](const std::vector<double>& density, double angle, int index, bool final) {
        PassState p = density_pass(c, density, angle);
        if (c.te) {
            p.kutta = general_kutta(*op, *c.te, p.W_m, p.W);
            if (config.kutta) {
                ExitAngleCorrection corr;
                corr.delta_tan = p.kutta->delta_tan;
                p.exit_angle_next = next_exit_angle(angle, corr, it.kutta_relaxation);
            }
        }
        GeneralPass r;
        r.index = index;
        r.final = final;
        r.density_change = change_norm(as_field(p.next), as_field(p.rho));
        r.velocity_change = previous ? velocity_change(as_field(p.W_m), as_field(p.W_theta), as_field(previous->W_m),
                                                       as_field(previous->W_theta))
                                     : nan();
        r.psi_change = previous ? change_norm(as_field(p.psi), as_field(previous->psi)) : nan();
        r.mach_max = *std::max_element(p.mach.begin(), p.mach.end());
        r.linear_residual = p.linear_residual;
        r.exit_angle = angle;
        r.exit_angle_next = p.exit_angle_next;
        r.kutta_residual = p.kutta ? p.kutta->residual : nan();
        r.kutta_normalized = p.kutta ? p.kutta->normalized : nan();
        r.tan_correction = (p.kutta && config.kutta) ? p.kutta->delta_tan : nan();
        r.inlet_mass_flow_error = p.inlet_error;
        r.outlet_mass_flow_error = p.outlet_error;
        passes.push_back(r);
        history.kutta_pass.push_back(config.kutta ? index : 0);
        history.density_change.push_back(r.density_change);
        history.velocity_change.push_back(r.velocity_change);
        history.psi_change.push_back(r.psi_change);
        history.residual.push_back(r.linear_residual);
        if (config.kutta) {
            history.exit_angle.push_back(angle);
            history.kutta_error.push_back(r.kutta_normalized);
        }
        if (!(r.mach_max < 1.0) && c.stop_at_sonic)
            throw GeneralSolverError(std::format("pass {}: the relative Mach number reaches {:.4f}: the subsonic density "
                                                 "iteration does not apply",
                                                 index, r.mach_max),
                                     make_solution(c, p, history, passes, false, index, inlet_angle));
        return p;
    };
    for (int n = 0; n < limit; ++n) {
        PassState p = evaluate(rho, beta, n, false);
        previous = std::move(p);
        const bool done = config.kutta ? coupled_converged(history, it) : density_converged(history, it);
        if (done) {
            const PassState final = evaluate(previous->next, previous->exit_angle_next, n + 1, true);
            return *make_solution(c, final, history, passes, true, n + 1, inlet_angle);
        }
        rho = previous->next;
        beta = previous->exit_angle_next;
    }
    throw GeneralSolverError(std::format("the general solve did not converge in {} passes (density change {:.3e}, velocity "
                                         "change {:.3e}{})",
                                         limit, history.density_change.back(), history.velocity_change.back(),
                                         config.kutta ? std::format(", eps_K {:.3e}", history.kutta_error.back()) : ""),
                             make_solution(c, *previous, history, passes, false, limit, inlet_angle));
}

}  // namespace BladeGen::Turbolenza
