#include "design/DesignSolid.hpp"

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRep_Builder.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax1.hxx>

#include <cmath>
#include <fstream>
#include <stdexcept>

namespace BladeGen::Design {

namespace {

Blade::MeridionalCurve scaled(const Math::BSpline2& c, double s)
{
    std::vector<Vec2> P;
    for (const auto& q : c.ctrl()) P.push_back(s * q);
    return Blade::MeridionalCurve(Math::BSpline2(c.degree(), c.knots(), std::move(P)));
}

/// The (x, r) B-spline placed in the meridional plane theta: (x, r) -> (x, r sin theta, r cos theta)
/// is linear, so mapping the control points gives the same curve exactly.
Handle(Geom_BSplineCurve) in_plane(const Math::BSpline2& c, double scale, double theta)
{
    const auto& U = c.knots();
    std::vector<double> k;
    std::vector<int> m;
    for (double u : U) {
        if (!k.empty() && u == k.back()) ++m.back();
        else {
            k.push_back(u);
            m.push_back(1);
        }
    }
    TColStd_Array1OfReal knots(1, static_cast<int>(k.size()));
    TColStd_Array1OfInteger mults(1, static_cast<int>(m.size()));
    for (std::size_t i = 0; i < k.size(); ++i) {
        knots(static_cast<int>(i) + 1) = k[i];
        mults(static_cast<int>(i) + 1) = m[i];
    }
    TColgp_Array1OfPnt poles(1, c.n_ctrl());
    for (int i = 0; i < c.n_ctrl(); ++i) {
        const auto& q = c.ctrl()[static_cast<std::size_t>(i)];
        const Vec3 p = Blade::cartesian(scale * q.x(), scale * q.y(), theta);
        poles(i + 1) = gp_Pnt(p.x(), p.y(), p.z());
    }
    return new Geom_BSplineCurve(poles, knots, mults, c.degree());
}

/// Face of revolution of a wall contour from theta0 over `sweep` toward +theta.
TopoDS_Shape wall(const Math::BSpline2& c, double scale, double theta0, double sweep)
{
    const TopoDS_Edge e = BRepBuilderAPI_MakeEdge(in_plane(c, scale, theta0)).Edge();
    // a rotation about -x increases theta (same axis as occt/OcctBlade's sector)
    BRepPrimAPI_MakeRevol rev(e, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(-1, 0, 0)), sweep);
    rev.Build();
    if (!rev.IsDone()) throw std::runtime_error("wall of revolution failed");
    return rev.Shape();
}

}  // namespace

gp_Trsf theta_rotation(double angle)
{
    gp_Trsf t;
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(-1, 0, 0)), angle);
    return t;
}

Scene build_scene(const BladeDesign& d, Evaluator& ev, const SceneOptions& opt)
{
    if (d.n_blades < 1) throw std::invalid_argument("build_scene: n_blades must be >= 1");
    constexpr double scale = 1e-3;  // design mm -> m (the unit build_design_solid is tested in)
    Scene s;
    s.shapes = build_design_solid(d, ev, scale, opt.rows, 0.0, true);
    s.n_blades = d.n_blades;
    s.pitch = 2.0 * Math::kPi / d.n_blades;
    // theta of the blade: the solid's inside point (mid row, mid chord), the sector is centred on it
    const auto r = ev.evaluate(d);
    const int j = r.blade.nj / 2, n = r.blade.ni - 1;
    const Vec3 c = 0.5 * (r.blade.at(n / 4, j) + r.blade.at(3 * n / 4, j));
    s.theta_blade = std::atan2(c.y(), c.z());
    s.copies = opt.annulus ? d.n_blades : 1;
    s.sweep = opt.annulus ? 2.0 * Math::kPi : s.pitch;
    s.theta_start = opt.annulus ? s.theta_blade - Math::kPi : s.theta_blade - 0.5 * s.pitch;

    BRep_Builder b;
    TopoDS_Compound all;
    b.MakeCompound(all);
    for (int k = 0; k < s.copies; ++k)
        b.Add(all, k == 0 ? s.shapes.blade : s.shapes.blade.Moved(TopLoc_Location(theta_rotation(k * s.pitch))));
    s.blades = all;
    const Blade::FlowPath fp(d.meridional);
    s.hub = wall(fp.hub(), scale, s.theta_start, s.sweep);
    s.casing = wall(fp.casing(), scale, s.theta_start, s.sweep);
    return s;
}

Occt::BladeShapes build_design_solid(const BladeDesign& d, Evaluator& ev, double scale, int rows, double fillet_radius, bool blade_only)
{
    if (!(scale > 0.0)) throw std::invalid_argument("build_design_solid: scale must be > 0");
    if (rows < 3) throw std::invalid_argument("build_design_solid: rows >= 3");
    auto e = d;
    e.mesh.n_rows = rows;
    e.mesh.row_clustering = 0.0;
    if (!blade_only) {  // the export resolution; the viewer uses the design mesh
        e.mesh.n_edge = std::max(d.mesh.n_edge, 30);
        e.mesh.n_body = std::max(d.mesh.n_body, 100);
    }
    const auto r = ev.evaluate(e);
    if (!r.ok()) throw std::invalid_argument("build_design_solid: the design has errors");
    auto loops = section_loops(r);
    for (auto& l : loops)
        for (auto& p : l) p *= scale;
    loops = Occt::extend_loops(loops, 1.0, 1.0);  // cross hub and casing; trimmed there
    const Blade::FlowPath fp(e.meridional);
    Occt::SolidOptions opt;
    opt.hub = scaled(fp.hub(), scale);
    opt.casing = scaled(fp.casing(), scale);
    opt.fillet_radius = fillet_radius;
    opt.blade_only = blade_only;
    // inside point: between the suction and pressure side at mid chord of the middle row
    const int j = r.blade.nj / 2, n = r.blade.ni - 1;
    const Vec3 inside = 0.5 * scale * (r.blade.at(n / 4, j) + r.blade.at(3 * n / 4, j));
    return Occt::build_shapes(loops, opt, inside);
}

void write_obj(const DesignResult& r, const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f.precision(10);
    f << "# BladeGen designer mesh: blade, hub, casing\n";
    int base = 1;
    auto grid = [&](const Grid& g, const char* name) {
        if (g.p.empty()) return;
        f << "o " << name << "\n";
        for (const auto& p : g.p) f << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
        for (int j = 0; j + 1 < g.nj; ++j)
            for (int i = 0; i + 1 < g.ni; ++i) {
                const int a = base + i + g.ni * j;
                f << "f " << a << " " << a + 1 << " " << a + 1 + g.ni << " " << a + g.ni << "\n";
            }
        base += static_cast<int>(g.p.size());
    };
    grid(r.blade, "blade");
    grid(r.hub_wall, "hub");
    grid(r.casing_wall, "casing");
}

void write_layer_csv(const BladeDesign& d, const DesignResult& r, const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f.precision(12);
    f << "layer,span,index,piece,x_c,y_c,mprime,theta,x,y,z\n";
    for (std::size_t k = 0; k < r.layers.size() && k < d.layers.size(); ++k) {
        const auto& L = r.layers[k];
        if (!L.error.empty()) continue;
        const auto& loop = L.section->loop;
        for (std::size_t i = 0; i < loop.xy.size(); ++i)
            f << d.layers[k].name << "," << d.layers[k].span << "," << i << ","
              << Profile::to_string(static_cast<Profile::Piece>(loop.piece[i])) << "," << loop.xy[i].x() << ","
              << loop.xy[i].y() << "," << L.conformal[i].x() << "," << L.conformal[i].y() << "," << L.xyz[i].x() << ","
              << L.xyz[i].y() << "," << L.xyz[i].z() << "\n";
    }
}

}  // namespace BladeGen::Design
