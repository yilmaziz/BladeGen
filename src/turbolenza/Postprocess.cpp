// Port of v1 postprocess.py (see Postprocess.hpp): blade_surfaces, chordwise_coordinate,
// surface_pressure_coefficient, _coefficient, discharge_angle_deg, exit_velocity_angle_deg,
// summarize, format_summary, _one_surface.
#include "turbolenza/Postprocess.hpp"

#include "turbolenza/Detail.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <set>

namespace BladeGen::Turbolenza {

using detail::sz;

const char* pressure_label(PressureConvention c) noexcept
{
    return c == PressureConvention::Dunavant ? "S = (Pt1 - p) / q1" : "Cp = (p - p1) / q1";
}

namespace {

const char* surface_name(Surface s) { return s == Surface::Suction ? "suction" : "pressure"; }

Surface one_surface(const std::vector<Surface>& names, const char* wall)
{
    std::set<int> found;
    for (Surface s : names) found.insert(static_cast<int>(s));
    if (found.size() != 1) {
        std::string list;
        for (int s : found) list += std::string(list.empty() ? "" : " and ") + surface_name(static_cast<Surface>(s));
        throw PostprocessError(std::format("the {} wall lies on {} along the blade: its surface changes, so it cannot be "
                                           "labeled one side (an inadmissible grid; grid.py sections 3 and 8)",
                                           wall, list));
    }
    return static_cast<Surface>(*found.begin());
}

std::size_t argmax(const std::vector<double>& v)
{
    std::size_t at = 0;
    for (std::size_t k = 0; k < v.size(); ++k)
        if (v[k] > v[at]) at = k;
    return at;
}

}  // namespace

std::vector<double> chordwise_coordinate(const std::vector<Point>& points, const Point& le, const Point& te)
{
    const double dx = te.m - le.m, dy = te.y - le.y;
    const double dd = dx * dx + dy * dy;
    std::vector<double> out(points.size());
    for (std::size_t k = 0; k < points.size(); ++k)
        out[k] = ((points[k].m - le.m) * dx + (points[k].y - le.y) * dy) / dd;
    return out;
}

BladeSurfaces blade_surfaces(const Solution& solution)
{
    const Grid& grid = *solution.grid;
    const Surface lower = one_surface(grid.wall_identity.lower.surface, "lower (eta = 0)");
    const Surface upper = one_surface(grid.wall_identity.upper.surface, "upper (eta = 1)");
    const int last = grid.n_eta() - 1;
    BladeSurfaces s;
    if (lower == Surface::Suction && upper == Surface::Pressure) {
        s.ss_j = 0;
        s.ps_j = last;
    } else if (lower == Surface::Pressure && upper == Surface::Suction) {
        s.ss_j = last;
        s.ps_j = 0;
    } else {
        throw PostprocessError(std::format("the walls lie on the {} and the {} surface: not one suction and one pressure "
                                           "surface",
                                           surface_name(lower), surface_name(upper)));
    }
    const PassageBoundaries& b = grid.boundaries;
    if (b.lower_blade_boundary.empty()) throw PostprocessError("the grid has no blade boundary (LE and TE points)");
    const Point le = b.lower_blade_boundary.front(), te = b.lower_blade_boundary.back();
    for (int i = grid.first_blade_i; i <= grid.last_blade_i; ++i) {
        s.stations.push_back(i);
        s.m.push_back(grid.m[sz(i)]);
        s.m_c.push_back((grid.m[sz(i)] - grid.m_le()) / grid.chord());
    }
    auto x_c = [&](int j) {
        std::vector<Point> pts;
        for (int i : s.stations) {
            Point p{grid.node_m(i, j), grid.node_y(i, j)};
            if (j == last) {  // the upper wall lies on the +pitch neighbour
                p.m = p.m - b.periodic_translation.m;
                p.y = p.y - b.periodic_translation.y;
            }
            pts.push_back(p);
        }
        return chordwise_coordinate(pts, le, te);
    };
    s.x_c_ss = x_c(s.ss_j);
    s.x_c_ps = x_c(s.ps_j);
    for (int i : s.stations) {
        s.mach_ss.push_back(solution.mach(i, s.ss_j));
        s.mach_ps.push_back(solution.mach(i, s.ps_j));
        s.pressure_ss.push_back(solution.static_state.pressure(i, s.ss_j));
        s.pressure_ps.push_back(solution.static_state.pressure(i, s.ps_j));
    }
    return s;
}

std::vector<double> pressure_coefficient_values(const std::vector<double>& p, const std::optional<InletState>& inlet,
                                                PressureConvention convention)
{
    if (!inlet)
        throw PostprocessError("the solution has no inlet state (mdot' given without an inlet Mach number): no reference "
                               "for a pressure coefficient");
    const double q1 = inlet->dynamic_pressure;
    if (!(std::isfinite(q1) && q1 > 0.0))
        throw PostprocessError(std::format("the inlet dynamic pressure is {}: a pressure coefficient needs q1 > 0",
                                           detail::py_repr(q1)));
    std::vector<double> out(p.size());
    for (std::size_t k = 0; k < p.size(); ++k)
        out[k] = convention == PressureConvention::Dunavant ? (inlet->total_pressure - p[k]) / q1
                                                            : (p[k] - inlet->static_pressure) / q1;
    return out;
}

SurfaceCoefficient surface_pressure_coefficient(const Solution& solution, PressureConvention convention)
{
    const BladeSurfaces s = blade_surfaces(solution);
    return SurfaceCoefficient{convention,
                              pressure_label(convention),
                              s.x_c_ss,
                              pressure_coefficient_values(s.pressure_ss, solution.inlet, convention),
                              s.x_c_ps,
                              pressure_coefficient_values(s.pressure_ps, solution.inlet, convention)};
}

double discharge_angle_deg(const Solution& solution) { return detail::degrees(solution.grid->exit_angle()); }

double exit_velocity_angle_deg(const Solution& solution)
{
    const Grid& g = *solution.grid;
    const std::size_t row = solution.flow.W_m.index(g.n_m() - 1, 0);
    const double W_m = detail::trapezoid(solution.flow.W_m.v.data() + row, g.eta);
    const double W_theta = detail::trapezoid(solution.flow.W_theta.v.data() + row, g.eta);
    return detail::degrees(std::atan2(W_theta, W_m));
}

SolutionSummary summarize(const Solution& solution)
{
    const BladeSurfaces s = blade_surfaces(solution);
    const Grid& grid = *solution.grid;
    if (solution.passes.empty()) throw PostprocessError("the solution has no pass record");
    SolutionSummary r;
    r.limitation = solution.limitation;
    r.converged = solution.converged;
    r.iterations = solution.iterations;
    r.mass_flow_per_depth = solution.mass_flow_per_depth;
    r.mass_flow_source = solution.mass_flow_source;
    if (solution.inlet) r.inlet_mach = solution.inlet->mach;
    r.inlet_angle_deg = detail::degrees(grid.inlet_angle());
    r.discharge_angle_deg = discharge_angle_deg(solution);
    r.exit_velocity_angle_deg = exit_velocity_angle_deg(solution);
    r.exit_angle_mismatch_deg = r.exit_velocity_angle_deg - r.discharge_angle_deg;
    r.turning_deg = r.inlet_angle_deg - r.discharge_angle_deg;
    std::vector<double> surface = s.mach_ss;
    surface.insert(surface.end(), s.mach_ps.begin(), s.mach_ps.end());
    r.surface_mach_max = *std::max_element(surface.begin(), surface.end());
    r.surface_mach_min = *std::min_element(surface.begin(), surface.end());
    r.ss_peak_mach = s.mach_ss[argmax(s.mach_ss)];
    r.ss_peak_x_c = s.x_c_ss[argmax(s.mach_ss)];
    r.ps_peak_mach = s.mach_ps[argmax(s.mach_ps)];
    r.ps_peak_x_c = s.x_c_ps[argmax(s.mach_ps)];
    r.kutta_residual = solution.kutta_residual;
    r.normalized_kutta_residual = solution.normalized_kutta_residual;
    const std::size_t at = argmax(solution.mach.v);
    r.mach_max = solution.mach.v[at];
    r.mach_max_node = {static_cast<int>(at) / solution.mach.n_eta, static_cast<int>(at) % solution.mach.n_eta};
    auto lo = [](const Field& f) { return *std::min_element(f.v.begin(), f.v.end()); };
    auto hi = [](const Field& f) { return *std::max_element(f.v.begin(), f.v.end()); };
    r.density_min = lo(solution.density);
    r.density_max = hi(solution.density);
    r.pressure_min = lo(solution.static_state.pressure);
    r.pressure_max = hi(solution.static_state.pressure);
    r.temperature_min = lo(solution.static_state.temperature);
    r.temperature_max = hi(solution.static_state.temperature);
    const PassRecord& last = solution.passes.back();
    r.mass_flow_error_max = last.mass_flow_error_max;
    r.mass_flow_error_blade_median = last.mass_flow_error_blade_median;
    r.linear_residual = last.linear_residual;
    r.density_change = last.density_change;
    r.velocity_change = last.velocity_change;
    return r;
}

std::string format_summary(const SolutionSummary& s)
{
    const std::string state = s.converged ? "converged" : "NOT CONVERGED";
    const std::string inlet = s.inlet_mach ? std::format("{:.4f}", *s.inlet_mach) : "-";
    std::string out;
    out += std::format("{} after {} passes; {}\n", state, s.iterations, s.limitation);
    out += std::format("inlet: M1 {}, beta1 {:.4f} deg, mdot' {:.6g} ({})\n", inlet, s.inlet_angle_deg,
                       s.mass_flow_per_depth, s.mass_flow_source);
    out += std::format("beta_out {:.6f} deg (far-exit velocity {:.6f} deg, mismatch {:+.2e} deg); turning {:.4f} deg\n",
                       s.discharge_angle_deg, s.exit_velocity_angle_deg, s.exit_angle_mismatch_deg, s.turning_deg);
    out += std::format("surface Mach {:.4f} .. {:.4f}; SS peak {:.4f} at x/c {:.4f}; PS peak {:.4f} at x/c {:.4f}\n",
                       s.surface_mach_min, s.surface_mach_max, s.ss_peak_mach, s.ss_peak_x_c, s.ps_peak_mach,
                       s.ps_peak_x_c);
    out += std::format("Kutta: R_K {:+.3e}, eps_K {:.3e}\n", s.kutta_residual, s.normalized_kutta_residual);
    out += std::format("field: max M_rel {:.4f} at node ({}, {}); rho {:.6g} .. {:.6g}; p {:.6g} .. {:.6g}; T {:.6g} .. "
                       "{:.6g}\n",
                       s.mach_max, s.mach_max_node[0], s.mach_max_node[1], s.density_min, s.density_max,
                       s.pressure_min, s.pressure_max, s.temperature_min, s.temperature_max);
    out += std::format("mass flow: station error max {:.2e}, blade median {:.2e}; linear residual {:.2e}; last changes "
                       "rho {:.2e}, W {:.2e}",
                       s.mass_flow_error_max, s.mass_flow_error_blade_median, s.linear_residual, s.density_change,
                       s.velocity_change);
    return out;
}

}  // namespace BladeGen::Turbolenza
