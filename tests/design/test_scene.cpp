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
    const gp_Pnt q = gp_Pnt(0.1, 0.0, 1.0).Transformed(theta_rotation(0.3));
    CHECK(theta_of(q) == Approx(0.3).margin(1e-15));
    CHECK(q.X() == Approx(0.1).margin(1e-15));
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
