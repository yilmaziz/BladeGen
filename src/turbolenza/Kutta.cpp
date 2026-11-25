// Port of v1 kutta.py, production rule (see Kutta.hpp): trailing_edge_station, kutta_residual,
// exit_angle_correction, next_exit_angle.
#include "turbolenza/Kutta.hpp"

#include "turbolenza/Detail.hpp"

#include <cmath>
#include <format>
#include <numbers>

namespace BladeGen::Turbolenza {

using detail::sz;

TrailingEdgeStation trailing_edge_station(const Grid& grid)
{
    const int i = grid.last_blade_i;
    const int k = i - grid.first_blade_i;
    const WallIdentity& id = grid.wall_identity;
    if (k < 0 || sz(k) >= id.lower.surface.size() || sz(k) >= id.upper.surface.size())
        throw KuttaError(std::format("the wall identity has no entry for the TE station {}", i));
    const Surface lower = id.lower.surface[sz(k)], upper = id.upper.surface[sz(k)];
    const int last = grid.n_eta() - 1;
    TrailingEdgeStation s;
    s.i = i;
    s.m = grid.m[sz(i)];
    s.m_offset = grid.m[sz(i)] - grid.m_te();
    if (lower == Surface::Suction && upper == Surface::Pressure) {
        s.ss_j = 0;
        s.ps_j = last;
    } else if (lower == Surface::Pressure && upper == Surface::Suction) {
        s.ss_j = last;
        s.ps_j = 0;
    } else {
        auto name = [](Surface x) { return x == Surface::Suction ? "suction" : "pressure"; };
        throw KuttaError(std::format("the TE station {} (m = {:.6g}) has its wall nodes on the {} and the {} surface, not "
                                     "one suction and one pressure node: an inadmissible TE station (grid.py section 8)",
                                     i, grid.m[sz(i)], name(lower), name(upper)));
    }
    return s;
}

KuttaResidual kutta_residual(const FlowField& flow, const Grid& grid)
{
    const TrailingEdgeStation station = trailing_edge_station(grid);
    const int n_m = grid.n_m(), n_eta = grid.n_eta();
    if (flow.W_m.n_m != n_m || flow.W_m.n_eta != n_eta || !flow.W_m.same_shape(flow.W_theta))
        throw KuttaError(
            std::format("the flow has shape ({}, {}), the grid ({}, {})", flow.W_m.n_m, flow.W_m.n_eta, n_m, n_eta));
    const int i = station.i, last = n_eta - 1;
    const double s0 = std::hypot(flow.W_m(i, 0), flow.W_theta(i, 0));
    const double s1 = std::hypot(flow.W_m(i, last), flow.W_theta(i, last));
    if (!std::isfinite(s0) || !std::isfinite(s1))
        throw KuttaError(std::format("a wall speed at the TE station {} is not finite", i));
    KuttaResidual r;
    r.station = station;
    r.W_ss = station.ss_j == 0 ? s0 : s1;
    r.W_ps = station.ps_j == 0 ? s0 : s1;
    r.residual = r.W_ss - r.W_ps;
    const double mean = 0.5 * (std::abs(r.W_ss) + std::abs(r.W_ps));
    r.normalized = mean > 0.0 ? std::abs(r.residual) / mean : 0.0;  // both speeds zero: Eq. 5-34 holds
    r.wall_difference = s0 - s1;
    return r;
}

ExitAngleCorrection exit_angle_correction(const KuttaResidual& residual, const FlowField& flow, const GridMetrics& metrics)
{
    const Grid& grid = *metrics.grid;
    const int i = residual.station.i;
    if (i != grid.last_blade_i)
        throw KuttaError(std::format("the residual is at station {}, the TE station of these metrics is {}", i,
                                     grid.last_blade_i));
    if (flow.W_m.n_m != grid.n_m() || flow.W_m.n_eta != grid.n_eta())
        throw KuttaError(std::format("the flow has shape ({}, {}), the grid ({}, {})", flow.W_m.n_m, flow.W_m.n_eta,
                                     grid.n_m(), grid.n_eta()));
    const double tan0 = metrics.tan_beta0[sz(i)], tan1 = metrics.tan_beta1[sz(i)];
    const double width = metrics.S[sz(i)], dm = grid.dm();
    const double integral = detail::trapezoid(flow.W_m.v.data() + flow.W_m.index(grid.n_m() - 1, 0), grid.eta);
    for (double v : {tan0, tan1, width, dm, integral, residual.wall_difference})
        if (!std::isfinite(v))
            throw KuttaError("a TE slope, width, spacing, speed or the far-exit meridional integral is not finite");
    if (!(width > 0.0 && dm > 0.0))
        throw KuttaError(std::format("the TE passage width {} and the station spacing {} must be positive",
                                     detail::py_repr(width), detail::py_repr(dm)));
    if (!(integral > 0.0))
        throw KuttaError(std::format("int W_m d(eta) at the far exit is {}: no flow leaves through it, the exit angle "
                                     "cannot carry an angular-momentum correction (Eq. 13-59)",
                                     detail::py_repr(integral)));
    ExitAngleCorrection c;
    c.cos_beta0 = 1.0 / std::sqrt(1.0 + tan0 * tan0);
    c.cos_beta1 = 1.0 / std::sqrt(1.0 + tan1 * tan1);
    c.wall_difference = residual.wall_difference;
    c.delta_c_theta = 2.0 * residual.wall_difference * dm / (width * (c.cos_beta0 + c.cos_beta1));  // Eq. 13-57 over r
    if (metrics.on_sheet) {
        // not in v1: Eq. 13-57 carries r at the TE station, Delta(r C_U) = 2 (W0 - W1) r_te dm / (S (cos beta0 +
        // cos beta1)), and Eq. 13-59 divides by r at the far exit: delta_c_theta = Delta(r C_U) / r_out
        c.radius_ratio = metrics.radius_at(i) / metrics.radius_at(grid.n_m() - 1);
        c.delta_c_theta = c.delta_c_theta * c.radius_ratio;
    }
    c.exit_meridional_integral = integral;
    c.passage_width = width;
    c.dm = dm;
    c.delta_tan = c.delta_c_theta / integral;  // Eq. 13-59
    return c;
}

double next_exit_angle(double exit_angle, const ExitAngleCorrection& correction, double relaxation)
{
    if (!(0.0 < relaxation && relaxation <= 1.0))
        throw KuttaError(std::format("the Kutta relaxation {} must lie in (0, 1]", detail::py_repr(relaxation)));
    if (!(std::isfinite(exit_angle) && std::abs(exit_angle) < 0.5 * std::numbers::pi))
        throw KuttaError(std::format("the exit angle {} rad is not within (-pi/2, pi/2)", detail::py_repr(exit_angle)));
    const double step = relaxation * correction.delta_tan;
    if (!std::isfinite(step))
        throw KuttaError(std::format("the correction of tan(beta_out) is {}: not finite",
                                     detail::py_repr(correction.delta_tan)));
    if (step == 0.0) return exit_angle;  // equal TE speeds: the angle stays exactly
    return std::atan(std::tan(exit_angle) + step);
}

}  // namespace BladeGen::Turbolenza
