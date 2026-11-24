/// @file turbolenza/ThroughflowBladeGen.cpp
/// @brief A BladeGen design -> the through-flow's quasi-normals and blade row (see the header).
#include "turbolenza/ThroughflowBladeGen.hpp"

#include "turbolenza/Cascade.hpp"

#include "blade/FlowPath.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace BladeGen::Turbolenza::Throughflow {

namespace {

using Math::Vec2;

/// u on [0, 3] where the wall has axial position z (z is linear in u on every piece, so z(u) is
/// monotone: bisection).
double u_at_z(const Math::BSpline2& w, double z)
{
    double a = w.t0(), b = w.t1();
    if (z < w.eval(a).x() - 1e-9 || z > w.eval(b).x() + 1e-9)
        throw std::invalid_argument("throughflow_geometry: a radial quasi-normal misses a wall");
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (a + b);
        if (w.eval(mid).x() < z) a = mid;
        else b = mid;
    }
    return 0.5 * (a + b);
}

/// The throat o (Eq. 6-53) of a placed section, normalized by its conformal chord c': the geometric
/// throat of the cascade module (closest_points of Cascade.hpp, as make_cascade: the suction side
/// against the neighbour's pressure side) on BladeGen's own side samples, placed as the design places
/// them (stagger gamma, the mirror sign) at the conformal pitch (2 pi / Z) / c'; the neighbour on
/// either side, the smaller distance (ours).
double throat_over_chord(const Design::LayerResult& lr, int n_blades, bool mirror)
{
    const auto& q = lr.section->quality;
    const double g = lr.gamma, sgn = mirror ? -1.0 : 1.0;
    auto placed = [&](const std::vector<Math::Vec2>& side, double shift) {
        std::vector<Point> out;
        out.reserve(side.size());
        for (const auto& p : side)
            out.push_back({p.x() * std::cos(g) - sgn * p.y() * std::sin(g), p.x() * std::sin(g) + sgn * p.y() * std::cos(g) + shift});
        return out;
    };
    const double pitch = 2.0 * std::numbers::pi / static_cast<double>(n_blades) / lr.chord_mprime;
    const auto ss = placed(q.ss.xy, 0.0);
    double best = std::numeric_limits<double>::infinity();
    for (double side : {1.0, -1.0}) best = std::min(best, closest_points(ss, placed(q.ps.xy, side * pitch)).distance);
    return best;
}

/// The chordwise position of the largest camber-line ordinate (ours: a / c).
double max_camber_position(const Profile::CamberLine& c)
{
    double best = -1.0, x = 0.5;
    for (const auto& p : c.point)
        if (std::abs(p.y()) > best) {
            best = std::abs(p.y());
            x = p.x();
        }
    return x;
}

}  // namespace

void wall_slope_curvature(const Math::BSpline2& wall, double u, double& phi, double& kappa)
{
    const auto d = wall.derivs(u, 2);
    const double z1 = d[1].x(), r1 = d[1].y(), z2 = d[2].x(), r2 = d[2].y();
    const double s = std::hypot(z1, r1);
    if (!(s > 0.0)) throw std::invalid_argument("wall_slope_curvature: degenerate wall tangent");
    phi = std::atan2(r1, z1);                       // Eq. 7-2
    kappa = -(z1 * r2 - r1 * z2) / (s * s * s);     // Eq. 7-9: -dphi/dm
}

BladeGenThroughflowGeometry throughflow_geometry(const Design::BladeDesign& d, Design::Evaluator& ev,
                                                 const BladeGenThroughflowOptions& opt)
{
    if (opt.n_inlet < 0 || opt.n_exit < 0) throw std::invalid_argument("throughflow_geometry: negative station count");
    if (!(opt.edge_gap > 0.0 && opt.edge_gap < 0.5)) throw std::invalid_argument("throughflow_geometry: edge_gap must be in (0, 0.5)");
    if (!(opt.length_scale > 0.0)) throw std::invalid_argument("throughflow_geometry: length_scale must be > 0");
    if (!std::isfinite(opt.rpm)) throw std::invalid_argument("throughflow_geometry: non-finite rpm");
    d.validate();
    const Blade::FlowPath fp(d.meridional);
    const auto& hub = fp.hub();
    const auto& cas = fp.casing();
    const double L = opt.length_scale;
    const auto& cfg = d.meridional;

    BladeGenThroughflowGeometry out;
    // the stations: per station either (z) radial or (u) equal-parameter
    std::vector<double> at;
    if (opt.layout == BladeGenThroughflowOptions::Layout::Radial) {
        const double z_in = std::max(hub.eval(hub.t0()).x(), cas.eval(cas.t0()).x());
        const double z_out = std::min(hub.eval(hub.t1()).x(), cas.eval(cas.t1()).x());
        const double gap = opt.edge_gap * (cfg.hub.z_te - cfg.hub.z_le);
        const double z_b = std::min(cfg.hub.z_le, cfg.casing.z_le) - gap;
        const double z_a = std::max(cfg.hub.z_te, cfg.casing.z_te) + gap;
        if (!(z_in < z_b) || !(z_a < z_out))
            throw std::invalid_argument("throughflow_geometry: the extensions are too short for radial quasi-normals");
        at.push_back(z_in);
        for (int k = 1; k <= opt.n_inlet; ++k) at.push_back(z_in + (z_b - z_in) * k / (opt.n_inlet + 1.0));
        at.push_back(z_b);
        out.qn_before_le = static_cast<int>(at.size()) - 1;
        at.push_back(z_a);
        out.qn_after_te = out.qn_before_le + 1;
        for (int k = 1; k <= opt.n_exit; ++k) at.push_back(z_a + (z_out - z_a) * k / (opt.n_exit + 1.0));
        at.push_back(z_out);
    } else {
        const double u_b = 1.0 - opt.edge_gap, u_a = 2.0 + opt.edge_gap;
        at.push_back(0.0);
        for (int k = 1; k <= opt.n_inlet; ++k) at.push_back(u_b * k / (opt.n_inlet + 1.0));
        at.push_back(u_b);
        out.qn_before_le = static_cast<int>(at.size()) - 1;
        at.push_back(u_a);
        out.qn_after_te = out.qn_before_le + 1;
        for (int k = 1; k <= opt.n_exit; ++k) at.push_back(u_a + (3.0 - u_a) * k / (opt.n_exit + 1.0));
        at.push_back(3.0);
    }
    const double omega = -2.0 * std::numbers::pi * opt.rpm / 60.0;  // book sign (header)
    for (std::size_t k = 0; k < at.size(); ++k) {
        double uh = at[k], uc = at[k];
        if (opt.layout == BladeGenThroughflowOptions::Layout::Radial) {
            uh = u_at_z(hub, at[k]);
            uc = u_at_z(cas, at[k]);
        }
        const Vec2 ph = hub.eval(uh), pc = cas.eval(uc);
        QuasiNormal q;
        q.z_hub = L * ph.x();
        q.r_hub = L * ph.y();
        q.z_shroud = opt.layout == BladeGenThroughflowOptions::Layout::Radial ? q.z_hub : L * pc.x();  // radial: the same z
        q.r_shroud = L * pc.y();
        double fh = 0.0, kh = 0.0, fs = 0.0, ks = 0.0;
        wall_slope_curvature(hub, uh, fh, kh);
        wall_slope_curvature(cas, uc, fs, ks);
        q.phi_hub = fh;
        q.phi_shroud = fs;
        q.kappa_hub = kh / L;  // curvature in 1/m
        q.kappa_shroud = ks / L;
        const int ik = static_cast<int>(k);
        if (ik == out.qn_after_te) q.hub_wall_omega = omega;  // the hub segment across the rotor
        out.quasi_normals.push_back(q);
    }

    // the blade row from the placed layers
    const Design::DesignResult res = ev.evaluate(d);
    if (!res.errors.empty()) throw std::invalid_argument("throughflow_geometry: " + res.errors.front());
    if (res.layers.size() != d.layers.size() || res.layers.empty())
        throw std::invalid_argument("throughflow_geometry: the design has no layers");
    for (const auto& lr : res.layers)
        if (!lr.error.empty()) throw std::invalid_argument("throughflow_geometry: layer: " + lr.error);
    const std::size_t mid = res.layers.size() / 2;
    BladeRowSpec& row = out.row;
    row.name = d.name;
    row.qn_inlet = out.qn_before_le;
    row.qn_outlet = out.qn_after_te;
    row.omega = omega;
    out.compressor_omega = omega;
    row.n_blades = d.n_blades;
    row.s_row = res.layers[mid].beta_le > res.layers[mid].beta_te ? 1 : -1;
    row.profile_family = d.layers[mid].uses_custom() ? "custom" : "NACA65";
    row.le_hub_z = L * cfg.hub.z_le;
    row.le_hub_r = L * cfg.hub.r_le;
    row.le_shroud_z = L * cfg.casing.z_le;
    row.le_shroud_r = L * cfg.casing.r_le;
    row.te_hub_z = L * cfg.hub.z_te;
    row.te_hub_r = L * cfg.hub.r_te;
    row.te_shroud_z = L * cfg.casing.z_te;
    row.te_shroud_r = L * cfg.casing.r_te;
    const double deg = std::numbers::pi / 180.0;
    for (std::size_t k = 0; k < res.layers.size(); ++k) {
        const auto& lr = res.layers[k];
        const auto& cam = lr.section->camber;
        BladeSectionGeometry g;
        g.span = d.layers[k].span;
        g.chord = L * lr.chord;
        g.stagger = lr.gamma;
        g.kappa1 = lr.beta_le * deg;
        g.kappa2 = lr.beta_te * deg;
        g.camber = g.kappa1 - g.kappa2;
        g.thickness_ratio = cam.max_thickness;
        g.max_thickness_position = cam.max_thickness_x;
        g.max_camber_position = max_camber_position(cam);
        g.throat = L * lr.r_le * lr.chord_mprime * throat_over_chord(lr, d.n_blades, d.mirror);  // at the LE radius
        row.sections.push_back(g);
    }
    return out;
}

}  // namespace BladeGen::Turbolenza::Throughflow
