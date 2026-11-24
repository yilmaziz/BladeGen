/// @file turbolenza/ThroughflowSheet.cpp
/// @brief A through-flow stream surface as the blade-to-blade stream sheet (see the header).
#include "turbolenza/ThroughflowSheet.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>

namespace BladeGen::Turbolenza::Throughflow {

namespace {

constexpr double kPi = std::numbers::pi;

[[noreturn]] void fail(const std::string& what)
{
    throw ThroughflowError(what);
}

/// f(m) linear between the quasi-normals k .. k + 1 that bracket m, the end values beyond (the header's rule).
double piecewise_linear(const std::vector<double>& m, const std::vector<double>& f, double x)
{
    if (x <= m.front()) return f.front();
    if (x >= m.back()) return f.back();
    const std::size_t k = static_cast<std::size_t>(std::upper_bound(m.begin(), m.end(), x) - m.begin()) - 1;
    const double t = (x - m[k]) / (m[k + 1] - m[k]);
    return f[k] + t * (f[k + 1] - f[k]);
}

}  // namespace

double ThroughflowSurfaceSheet::radius(double x) const
{
    return piecewise_linear(m, r, x);
}

double ThroughflowSurfaceSheet::thickness_at(double x) const
{
    return piecewise_linear(m, b, x);
}

ThroughflowSurfaceSheet throughflow_surface_sheet(const ThroughflowCase& c, const ThroughflowResult& result, int row,
                                                  int surface, SheetThickness thickness)
{
    if (row < 0 || static_cast<std::size_t>(row) >= c.rows.size())
        fail(std::format("throughflow_surface_sheet: row {} of {} rows", row, c.rows.size()));
    const BladeRowSpec& br = c.rows[static_cast<std::size_t>(row)];
    const int nqn = static_cast<int>(c.quasi_normals.size());
    const int qi = br.qn_inlet, qo = br.qn_outlet < 0 ? br.qn_inlet + 1 : br.qn_outlet;
    if (qi < 0 || qo <= qi || qo >= nqn) fail("throughflow_surface_sheet: the row's quasi-normals are out of range");
    if (static_cast<int>(result.stations.size()) != nqn || result.stream_surfaces.empty())
        fail(std::format("throughflow_surface_sheet: the result has {} of {} quasi-normals (a choked or failed pass)",
                         result.stations.size(), nqn));
    if (surface < 0 || static_cast<std::size_t>(surface) >= result.stream_surfaces.size())
        fail(std::format("throughflow_surface_sheet: stream surface {} of {}", surface, result.stream_surfaces.size()));
    const std::size_t j = static_cast<std::size_t>(surface);
    const ThroughflowStreamSurface& ss = result.stream_surfaces[j];
    ThroughflowSurfaceSheet s;
    s.row = row;
    s.surface = surface;
    s.qn_inlet = qi;
    s.qn_outlet = qo;
    for (const auto& other : c.rows) {  // the blade-free stretch upstream of the row ends at qi
        const int out = other.qn_outlet < 0 ? other.qn_inlet + 1 : other.qn_outlet;
        if (out <= qi) s.qn_free = std::max(s.qn_free, out);
    }
    s.s_row = br.s_row;
    s.omega = br.omega;
    s.mass_fraction = ss.mass_fraction;
    s.thickness = thickness;
    s.m = ss.m;
    s.z = ss.z;
    s.r = ss.r;
    s.b = thickness == SheetThickness::Geometric ? ss.b : ss.b_mass;
    for (std::size_t k = 1; k < s.m.size(); ++k)
        if (!(s.m[k] > s.m[k - 1])) fail("throughflow_surface_sheet: m does not increase along the surface");
    for (const auto& st : result.stations) {
        const ThroughflowNode& nd = st.nodes[j];
        s.Pt.push_back(nd.Pt);
        s.Tt.push_back(nd.Tt);
        s.rCtheta.push_back(nd.rCtheta);
        s.Wm.push_back(nd.Wm);
    }
    // the LE and TE where the surface's straight segment QN1 -> QN2 crosses the row's edge lines (p. 287)
    const std::size_t a = static_cast<std::size_t>(qi), e = static_cast<std::size_t>(qo);
    auto cross = [&](double hz, double hr, double tz, double tr, const char* what) {
        const double t = line_crossing_fraction(s.z[a], s.r[a], s.z[e], s.r[e], hz, hr, tz, tr);
        if (!(t > 0.0 && t < 1.0))
            fail(std::format("throughflow_surface_sheet: surface {} meets the {} line outside its row segment (t = {})", surface,
                             what, t));
        return t;
    };
    const double tl = cross(br.le_hub_z, br.le_hub_r, br.le_shroud_z, br.le_shroud_r, "LE");
    const double tt = cross(br.te_hub_z, br.te_hub_r, br.te_shroud_z, br.te_shroud_r, "TE");
    if (!(tt > tl)) fail("throughflow_surface_sheet: the TE crossing lies upstream of the LE crossing");
    s.m_le = s.m[a] + tl * (s.m[e] - s.m[a]);
    s.m_te = s.m[a] + tt * (s.m[e] - s.m[a]);
    s.z_le = s.z[a] + tl * (s.z[e] - s.z[a]);
    s.r_le = s.r[a] + tl * (s.r[e] - s.r[a]);
    s.z_te = s.z[a] + tt * (s.z[e] - s.z[a]);
    s.r_te = s.r[a] + tt * (s.r[e] - s.r[a]);
    s.beta2 = result.stations[e].nodes[j].beta;  // in the frame of QN2, the row's
    if (j == 0) s.span = 0.0;  // the wall surfaces are the walls (as the through-flow)
    else if (j + 1 == result.stream_surfaces.size()) s.span = 1.0;
    else
        s.span = 0.5 * (line_crossing_fraction(br.le_hub_z, br.le_hub_r, br.le_shroud_z, br.le_shroud_r, s.z[a], s.r[a], s.z[e], s.r[e]) +
                        line_crossing_fraction(br.te_hub_z, br.te_hub_r, br.te_shroud_z, br.te_shroud_r, s.z[a], s.r[a], s.z[e], s.r[e]));
    return s;
}

std::shared_ptr<const StreamSheet> solver_stream_sheet(const ThroughflowSurfaceSheet& s, double L, double r_ref, double M_lo,
                                                       double M_hi, double step)
{
    if (!(std::isfinite(L) && L > 0.0)) throw std::invalid_argument("solver_stream_sheet: the length L must be > 0");
    if (!(std::isfinite(r_ref) && r_ref > 0.0)) throw std::invalid_argument("solver_stream_sheet: r_ref must be > 0");
    if (!(std::isfinite(M_lo) && std::isfinite(M_hi) && M_lo < M_hi))
        throw std::invalid_argument("solver_stream_sheet: the range must have M_lo < M_hi");
    if (!(std::isfinite(step) && step > 0.0)) throw std::invalid_argument("solver_stream_sheet: the step must be > 0");
    const int n = std::max(1, static_cast<int>(std::ceil((M_hi - M_lo) / step - 1e-9)));
    std::vector<double> M(static_cast<std::size_t>(n + 1)), r(M.size()), b(M.size());
    for (int q = 0; q <= n; ++q) {
        const std::size_t k = static_cast<std::size_t>(q);
        M[k] = q == n ? M_hi : M_lo + (M_hi - M_lo) * static_cast<double>(q) / static_cast<double>(n);
        const double m = s.m_le + M[k] * L;
        r[k] = s.radius(m) / L;
        b[k] = s.thickness_at(m);
    }
    return std::make_shared<const StreamSheet>(std::move(M), std::move(r), std::move(b), r_ref);
}

SurfaceInletState surface_inlet_state(const ThroughflowSurfaceSheet& s, const ThroughflowGas& gas, double m)
{
    // linear between the quasi-normals of the blade-free stretch qn_free .. qn_inlet (rC_theta, H, s conserved,
    // Eqs. 7-40 .. 7-42); QN 0's state before QN 0, QN1's beyond QN1 (ours); inside an upstream row: refused
    const std::size_t k0 = static_cast<std::size_t>(s.qn_free), qi = static_cast<std::size_t>(s.qn_inlet);
    if (k0 > 0 && m < s.m[k0])
        fail(std::format("surface_inlet_state: m = {} m lies upstream of quasi-normal {} (m = {} m), inside the row before: "
                         "shorten the upstream extension",
                         m, k0, s.m[k0]));
    const auto b = static_cast<std::ptrdiff_t>(k0), e = static_cast<std::ptrdiff_t>(qi + 1);
    const std::vector<double> mu(s.m.begin() + b, s.m.begin() + e);
    auto at = [&](const std::vector<double>& f) {
        const std::vector<double> fu(f.begin() + b, f.begin() + e);
        return fu.size() == 1 ? fu.front() : piecewise_linear(mu, fu, m);
    };
    SurfaceInletState o;
    o.m = m;
    o.r = s.radius(m);
    o.Pt = at(s.Pt);
    o.Tt = at(s.Tt);
    o.Wm = at(s.Wm);
    o.Ctheta = at(s.rCtheta) / o.r;
    const double cp = gas.cp(), C2 = o.Wm * o.Wm + o.Ctheta * o.Ctheta;
    const double h = cp * o.Tt - 0.5 * C2;
    if (!(h > 0.0)) fail(std::format("surface_inlet_state: no static state at m = {} (h = {} J/kg)", m, h));
    o.T = h / cp;
    o.a = std::sqrt(gas.gamma * gas.gas_constant * o.T);
    o.mach = std::sqrt(C2) / o.a;
    o.alpha = std::atan2(o.Ctheta, o.Wm);
    return o;
}

BladeToBladeOnSurface blade_to_blade_on_surface(const ThroughflowCase& c, const ThroughflowResult& result, int row,
                                                int surface, StreamSurfaceBlade blade, CaseSettings base,
                                                const BladeToBladeOptions& opt)
{
    if (blade.n_blades <= 0) throw std::invalid_argument("blade_to_blade_on_surface: the blade needs n_blades > 0");
    if (!(std::isfinite(blade.chord) && blade.chord > 0.0))
        throw std::invalid_argument("blade_to_blade_on_surface: the blade needs a chord > 0 (its normalizing length)");
    if (!(opt.margin_chords >= 0.0)) throw std::invalid_argument("blade_to_blade_on_surface: margin_chords must be >= 0");
    BladeToBladeOnSurface o;
    o.sheet = throughflow_surface_sheet(c, result, row, surface, opt.thickness);
    const ThroughflowSurfaceSheet& s = o.sheet;
    const double L = blade.chord;
    double r_ref = 0.0;
    if (blade.sheet) {
        r_ref = blade.sheet->r_ref();
    } else {  // a blade by hand: the pitch at the through-flow's r_LE
        r_ref = s.r_le / L;
        blade.solidity = L * static_cast<double>(blade.n_blades) / (2.0 * kPi * s.r_le);
    }
    o.r_le_ratio = s.r_le / (r_ref * L);
    // the sheet covers the grid's nominal far inlet and exit, plus the margin
    const double M_in = -base.upstream_extension_chords;
    const double M_te = std::max(std::cos(blade.stagger_deg * kPi / 180.0), (s.m_te - s.m_le) / L);
    const double M_out = M_te + base.downstream_extension_chords;
    blade.sheet = solver_stream_sheet(s, L, r_ref, M_in - opt.margin_chords, M_out + opt.margin_chords, opt.table_step_chords);
    if (blade.source.empty()) blade.source = "through-flow stream surface";
    blade.source += std::format(" on through-flow surface {} (psi {:.4g})", surface, s.mass_fraction);
    o.blade = std::move(blade);
    // the inlet at the nominal far inlet, and the row's motion; mirrored for s_row = -1 (the header)
    o.inlet = surface_inlet_state(s, c.gas, s.m_le + M_in * L);
    o.frame_sign = s.s_row;
    const double fs = static_cast<double>(s.s_row);
    CaseSettings& st = base;
    st.stream_surface = StreamSurfaceKind::Layer;
    st.gamma = c.gas.gamma;
    st.gas_constant = c.gas.gas_constant;
    st.inlet_total_pressure = o.inlet.Pt;
    st.inlet_total_temperature = o.inlet.Tt;
    const double alpha_deg = fs * o.inlet.alpha * 180.0 / kPi;
    if (s.omega != 0.0) {
        st.row = RowMotion::Rotating;
        o.rpm = fs * (-60.0 * s.omega / (2.0 * kPi));
        st.rpm = o.rpm;
        st.omega = 0.0;  // derived by build_case from rpm
        st.inlet_absolute_mach = o.inlet.mach;
        st.inlet_absolute_flow_angle_deg = alpha_deg;
    } else {
        st.row = RowMotion::Stationary;
        st.omega = 0.0;
        st.rpm = 0.0;
        st.inlet_mach = o.inlet.mach;
        st.inlet_flow_angle_deg = alpha_deg;
    }
    o.settings = std::move(st);
    return o;
}

}  // namespace BladeGen::Turbolenza::Throughflow
