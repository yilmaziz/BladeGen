#include "design/DesignSolid.hpp"

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

}  // namespace

Occt::BladeShapes build_design_solid(const BladeDesign& d, Evaluator& ev, double scale, int rows, double fillet_radius)
{
    if (!(scale > 0.0)) throw std::invalid_argument("build_design_solid: scale must be > 0");
    if (rows < 3) throw std::invalid_argument("build_design_solid: rows >= 3");
    auto e = d;
    e.mesh.n_rows = rows;
    e.mesh.row_clustering = 0.0;
    e.mesh.n_edge = std::max(d.mesh.n_edge, 30);
    e.mesh.n_body = std::max(d.mesh.n_body, 100);
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
