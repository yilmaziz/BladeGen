// Tests of design/DesignSolid build_scene: the OCCT 3-D view's geometry (one blade or the whole
// annulus, hub and casing surfaces of revolution over one pitch or 360 deg).
#include "design/DesignSolid.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <BRepAdaptor_Surface.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace BladeGen;
using namespace BladeGen::Design;
using Catch::Approx;

namespace {

std::vector<TopoDS_Shape> solids(const TopoDS_Shape& s)
{
    std::vector<TopoDS_Shape> out;
    for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next()) out.push_back(ex.Current());
    return out;
}

gp_Pnt centroid(const TopoDS_Shape& s)
{
    GProp_GProps g;
    BRepGProp::VolumeProperties(s, g, 1e-7);  // adaptive, as Occt::volume
    return g.CentreOfMass();
}

double theta_of(const gp_Pnt& p) { return std::atan2(p.Y(), p.Z()); }

/// Largest distance of face points from the wall contour (x, r) polyline.
double wall_error(const TopoDS_Shape& face, const std::vector<Math::Vec2>& xr)
{
    double worst = 0.0;
    for (TopExp_Explorer ex(face, TopAbs_FACE); ex.More(); ex.Next()) {
        BRepAdaptor_Surface S(TopoDS::Face(ex.Current()));
        for (int i = 0; i <= 8; ++i)
            for (int j = 0; j <= 4; ++j) {
                const double u = S.FirstUParameter() + (S.LastUParameter() - S.FirstUParameter()) * i / 8.0;
                const double v = S.FirstVParameter() + (S.LastVParameter() - S.FirstVParameter()) * j / 4.0;
                const gp_Pnt p = S.Value(u, v);
                const Math::Vec2 q(p.X(), std::hypot(p.Y(), p.Z()));
                double d = 1e300;
                for (std::size_t k = 1; k < xr.size(); ++k) {
                    const Math::Vec2 a = xr[k - 1], ab = xr[k] - a;
                    const double t = std::clamp((q - a).dot(ab) / ab.squaredNorm(), 0.0, 1.0);
                    d = std::min(d, (q - (a + t * ab)).norm());
                }
                worst = std::max(worst, d);
            }
    }
    return worst;
}

}  // namespace

TEST_CASE("Scene: one blade on a pitch of hub and casing, the whole annulus", "[design][scene]")
{
    auto d = BladeDesign::defaults();
    d.n_blades = 12;  // fewer copies: a fast annulus
    Evaluator ev;
    const auto one = build_scene(d, ev, {false, 9});
    const auto all = build_scene(d, ev, {true, 9});
    REQUIRE_FALSE(one.shapes.blade.IsNull());
    CHECK(Occt::is_valid(one.shapes.blade));
    CHECK(one.copies == 1);
    CHECK(all.copies == 12);
    CHECK(one.pitch == Approx(2.0 * Math::kPi / 12).epsilon(1e-15));
    CHECK(one.sweep == Approx(one.pitch).epsilon(1e-15));
    CHECK(all.sweep == Approx(2.0 * Math::kPi).epsilon(1e-15));
    // the single-blade sector is centred on the blade
    CHECK(one.theta_start + 0.5 * one.sweep == Approx(one.theta_blade).margin(1e-15));

    const auto b1 = solids(one.blades), bn = solids(all.blades);
    REQUIRE(b1.size() == 1);
    REQUIRE(bn.size() == 12);
    const double v0 = Occt::volume(one.shapes.blade);
    CHECK(v0 > 0.0);
    // copy k is the blade located by the rotation of k pitches toward +theta (shared geometry)
    for (std::size_t k = 0; k < bn.size(); ++k) {
        CHECK(bn[k].TShape() == all.shapes.blade.TShape());
        const gp_Trsf t = bn[k].Location().Transformation(), e = theta_rotation(static_cast<double>(k) * all.pitch);
        for (int r = 1; r <= 3; ++r)
            for (int c = 1; c <= 4; ++c) CHECK(t.Value(r, c) == Approx(e.Value(r, c)).margin(1e-14));
    }
    // ... with the same volume, its centroid one pitch further (integrated for three copies)
    const gp_Pnt c0 = centroid(bn[0]);
    for (const std::size_t k : {std::size_t{1}, std::size_t{11}}) {
        CHECK(Occt::volume(bn[k]) == Approx(v0).epsilon(1e-9));
        const gp_Pnt ck = centroid(bn[k]);
        const double dt = std::remainder(theta_of(ck) - theta_of(c0) - static_cast<double>(k) * all.pitch, 2.0 * Math::kPi);
        CHECK(std::abs(dt) < 1e-9);
        CHECK(std::hypot(ck.Y(), ck.Z()) == Approx(std::hypot(c0.Y(), c0.Z())).epsilon(1e-9));
        CHECK(ck.X() == Approx(c0.X()).epsilon(1e-9));
    }
    // theta_rotation turns toward +theta: (x, r sin t, r cos t) -> theta + a
    const gp_Pnt turned = gp_Pnt(0.1, 0.0, 1.0).Transformed(theta_rotation(0.3));
    CHECK(theta_of(turned) == Approx(0.3).margin(1e-15));
    CHECK(turned.X() == Approx(0.1).margin(1e-15));
    // neighbouring blades do not touch
    CHECK(Occt::distance(bn[0], bn[1]) > 1e-4);

    // the walls are surfaces of revolution of the hub and casing contours (metres)
    const Blade::FlowPath fp(d.meridional);
    auto scaled = [](std::vector<Math::Vec2> v) {
        for (auto& q : v) q *= 1e-3;
        return v;
    };
    const auto hub = scaled(fp.surface(0.0).polyline(400)), cas = scaled(fp.surface(1.0).polyline(400));
    CHECK(wall_error(one.hub, hub) < 1e-7);
    CHECK(wall_error(one.casing, cas) < 1e-7);
    CHECK(wall_error(all.hub, hub) < 1e-7);
    // area of a surface of revolution is sweep * integral of r ds: 360 deg = n_blades pitches
    CHECK(Occt::area(all.hub) == Approx(12.0 * Occt::area(one.hub)).epsilon(1e-9));
    CHECK(Occt::area(all.casing) == Approx(12.0 * Occt::area(one.casing)).epsilon(1e-9));
    CHECK(Occt::area(one.casing) > Occt::area(one.hub));  // larger radius
}

TEST_CASE("Scene: invalid input", "[design][scene]")
{
    Evaluator ev;
    auto d = BladeDesign::defaults();
    d.n_blades = 0;
    CHECK_THROWS_AS(build_scene(d, ev), std::invalid_argument);
    d = BladeDesign::defaults();
    CHECK_THROWS_AS(build_scene(d, ev, {false, 2}), std::invalid_argument);  // rows >= 3
    d.layers[2].method = Method::Custom;
    d.layers[2].custom.le.radius = 0.0;  // a layer error
    CHECK_THROWS_AS(build_scene(d, ev), std::invalid_argument);
}

namespace {

/// Signed volume enclosed by a closed, outward-oriented triangle mesh (divergence theorem).
double enclosed_volume(const TriangleMesh& m)
{
    double v = 0.0;
    for (const auto& [a, b, c] : m.tri) {
        const auto& p = m.p[static_cast<std::size_t>(a)];
        const auto& q = m.p[static_cast<std::size_t>(b)];
        const auto& s = m.p[static_cast<std::size_t>(c)];
        v += p.dot(q.cross(s)) / 6.0;
    }
    return v;
}

}  // namespace

TEST_CASE("Scene mesh: the design rows in triangles, closed, outward, close to the exact skin", "[design][scene][mesh]")
{
    const auto d = BladeDesign::defaults();
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    const auto m = build_scene_mesh(d, r, false);
    const int ni = r.blade.ni, nj = r.blade.nj;
    // every design mesh point is a vertex (metres)
    REQUIRE(m.blade.p.size() >= r.blade.p.size());
    for (std::size_t k = 0; k < r.blade.p.size(); k += 7) CHECK((m.blade.p[k] - 1e-3 * r.blade.p[k]).norm() < 1e-15);
    REQUIRE(m.blade.n.size() == m.blade.p.size());
    for (const auto& n : m.blade.n) CHECK(n.norm() == Approx(1.0).margin(1e-12));
    for (const auto& t : m.blade.tri)
        for (const int v : t) CHECK((v >= 0 && static_cast<std::size_t>(v) < m.blade.p.size()));

    // closed and outward: the enclosed volume is positive and equals the exact OCCT solid's
    const double v_mesh = enclosed_volume(m.blade);
    const auto exact = build_design_solid(d, ev, 1e-3, 11, 0.0, true);
    const double v_exact = Occt::volume(exact.blade);
    CHECK(v_mesh > 0.0);
    CHECK(v_mesh == Approx(v_exact).epsilon(1e-2));  // measured +0.47 %: the caps are flat, the walls curved
    // the normals at the mid row point away from the blade's inside
    const int N = ni - 1, jm = nj / 2;
    const Math::Vec3 inside = 0.5e-3 * (r.blade.at(N / 4, jm) + r.blade.at(3 * N / 4, jm));
    for (const int i : {N / 4, 3 * N / 4}) {
        const auto k = static_cast<std::size_t>(i + ni * jm);
        CHECK(m.blade.n[k].dot(m.blade.p[k] - inside) > 0.0);
    }

    // flat triangles against the C2 skin through the same rows: largest gap at the cell centres
    const auto loops = [&] {
        auto l = section_loops(r);
        for (auto& row : l)
            for (auto& q : row) q *= 1e-3;
        return l;
    }();
    const auto skin = Occt::skin_loops(loops);
    double gap = 0.0;
    for (int j = 0; j + 1 < nj; ++j)
        for (int i = 0; i + 1 < ni; ++i) {
            const auto P = [&](int a, int b) { return m.blade.p[static_cast<std::size_t>(a + ni * b)]; };
            const Math::Vec3 mid = 0.25 * (P(i, j) + P(i + 1, j) + P(i, j + 1) + P(i + 1, j + 1));
            const double u = 0.5 * (skin.u[static_cast<std::size_t>(i)] + skin.u[static_cast<std::size_t>(i + 1)]);
            const double v = 0.5 * (skin.v[static_cast<std::size_t>(j)] + skin.v[static_cast<std::size_t>(j + 1)]);
            const gp_Pnt s = skin.surface->Value(u, v);
            gap = std::max(gap, (mid - Math::Vec3(s.X(), s.Y(), s.Z())).norm());
        }
    INFO("largest gap between the display mesh and the exact skin: " << gap * 1e3 << " mm");
    CHECK(gap < 0.1e-3);  // measured 0.054 mm (1/800 of the chord): below a pixel at whole-blade zoom

    // walls: on their contours, over one pitch centred on the blade, normals into the passage
    const Blade::FlowPath fp(d.meridional);
    auto scaled = [](std::vector<Math::Vec2> v) {
        for (auto& q : v) q *= 1e-3;
        return v;
    };
    const auto hub = scaled(fp.surface(0.0).polyline(400)), cas = scaled(fp.surface(1.0).polyline(400));
    auto on = [](const std::vector<Math::Vec2>& poly, const Math::Vec3& p) {
        const Math::Vec2 q(p.x(), std::hypot(p.y(), p.z()));
        double best = 1e300;
        for (std::size_t k = 1; k < poly.size(); ++k) {
            const Math::Vec2 a = poly[k - 1], ab = poly[k] - a;
            const double t = std::clamp((q - a).dot(ab) / ab.squaredNorm(), 0.0, 1.0);
            best = std::min(best, (q - (a + t * ab)).norm());
        }
        return best;
    };
    // 1e-7 m: the sag of the 400-segment reference polyline (the vertices are exact surface points)
    for (std::size_t k = 0; k < m.hub.p.size(); k += 11) CHECK(on(hub, m.hub.p[k]) < 1e-7);
    for (std::size_t k = 0; k < m.casing.p.size(); k += 11) CHECK(on(cas, m.casing.p[k]) < 1e-7);
    double tlo = 1e9, thi = -1e9;
    for (const auto& p : m.hub.p) {
        const double t = m.theta_blade + std::remainder(std::atan2(p.y(), p.z()) - m.theta_blade, 2.0 * Math::kPi);
        tlo = std::min(tlo, t);
        thi = std::max(thi, t);
    }
    CHECK(tlo == Approx(m.theta_blade - 0.5 * m.pitch).margin(1e-12));
    CHECK(thi == Approx(m.theta_blade + 0.5 * m.pitch).margin(1e-12));
    for (std::size_t k = 0; k < m.hub.p.size(); k += 13) {
        const Math::Vec3 rad(0.0, m.hub.p[k].y(), m.hub.p[k].z());
        CHECK(m.hub.n[k].dot(rad) > 0.0);
        CHECK(m.casing.n[k].dot(Math::Vec3(0.0, m.casing.p[k].y(), m.casing.p[k].z())) < 0.0);
    }
    CHECK_FALSE(m.blade_lines.empty());
    CHECK(m.blade_edges.size() == 4);
    CHECK(m.copies == 1);
    // the annulus: 360 deg walls, n_blades copies
    const auto a = build_scene_mesh(d, r, true);
    CHECK(a.copies == d.n_blades);
    CHECK(a.sweep == Approx(2.0 * Math::kPi).epsilon(1e-15));
    CHECK(a.blade.p.size() == m.blade.p.size());  // one blade; the viewer turns the copies
}

TEST_CASE("Scene mesh: invalid input", "[design][scene][mesh]")
{
    auto d = BladeDesign::defaults();
    CHECK_THROWS_AS(build_scene_mesh(d, DesignResult{}, false), std::invalid_argument);
    Evaluator ev;
    const auto r = ev.evaluate(d);
    d.n_blades = 0;
    CHECK_THROWS_AS(build_scene_mesh(d, r, false), std::invalid_argument);
}
