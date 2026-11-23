// Tests of the centres of gravity (design/BladeDesign.hpp): polygon and closed-surface centroids
// against closed forms, the section centre of gravity on cylindrical layers against the chord-frame
// centroid mapped by the placement, and the blade's centre of gravity against the exact OCCT solid.
#include "design/BladeDesign.hpp"
#include "design/DesignSolid.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Design;
using Catch::Approx;
using Math::Vec2;
using Math::Vec3;

TEST_CASE("Polygon centroid: rectangle and triangle in closed form, either orientation", "[design][cg]")
{
    const std::vector<Vec2> rect = {{1, 2}, {3, 2}, {3, 6}, {1, 6}};
    auto r = polygon_centroid(rect);
    CHECK(r.area == Approx(8.0).epsilon(1e-15));
    CHECK((r.centroid - Vec2(2, 4)).norm() < 1e-14);
    const std::vector<Vec2> tri = {{0, 0}, {0, 3}, {4, 0}, {0, 0}};  // clockwise, closing point repeated
    r = polygon_centroid(tri);
    CHECK(r.area == Approx(6.0).epsilon(1e-15));
    CHECK((r.centroid - Vec2(4.0 / 3.0, 1.0)).norm() < 1e-14);
    // far from the origin: no cancellation
    std::vector<Vec2> far = rect;
    for (auto& p : far) p += Vec2(1e6, -2e6);
    r = polygon_centroid(far);
    CHECK(r.area == Approx(8.0).epsilon(1e-9));
    CHECK((r.centroid - Vec2(2 + 1e6, 4 - 2e6)).norm() < 1e-8);
    CHECK_THROWS_AS(polygon_centroid({{0, 0}, {1, 1}}), std::invalid_argument);
    CHECK_THROWS_AS(polygon_centroid({{0, 0}, {1, 1}, {2, 2}}), std::invalid_argument);  // zero area
}

TEST_CASE("Closed-surface centroid: polygonal prism and a sheared prism in closed form", "[design][cg]")
{
    const int N = 40, nj = 4;
    const double a = 3.0, x0 = 5.0, y0 = -2.0, z0 = 1.0, h = 6.0;
    for (const double shear : {0.0, 0.4}) {
        Grid g;
        g.ni = N + 1;
        g.nj = nj;
        for (int j = 0; j < nj; ++j) {
            const double z = z0 + h * j / (nj - 1);
            for (int i = 0; i <= N; ++i) {
                const double f = 2.0 * Math::kPi * (i % N) / N;
                g.p.emplace_back(x0 + shear * (z - z0) + a * std::cos(f), y0 + a * std::sin(f), z);
            }
        }
        const auto m = mass_properties(g);
        const double area = 0.5 * N * a * a * std::sin(2.0 * Math::kPi / N);  // the inscribed N-gon
        CHECK(m.volume == Approx(area * h).epsilon(1e-12));  // Cavalieri: the shear keeps the volume
        CHECK((m.centroid - Vec3(x0 + shear * 0.5 * h, y0, z0 + 0.5 * h)).norm() < 1e-12);
    }
    Grid bad;
    bad.ni = 5;
    bad.nj = 2;
    bad.p.assign(10, Vec3::Zero());
    CHECK_THROWS_AS(mass_properties(bad), std::invalid_argument);
}

TEST_CASE("Section centre of gravity on cylindrical layers: the chord-frame centroid, fixed by a re-stagger about it", "[design][cg]")
{
    auto d = BladeDesign::defaults();
    d.meridional.hub = {-15.0, 50.0, 15.0, 50.0};
    d.meridional.casing = {-15.0, 80.0, 15.0, 80.0};
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    for (std::size_t k = 0; k < r.layers.size(); ++k) {
        const auto& L = r.layers[k];
        REQUIRE(L.surface.cylindrical());
        // on a cylinder the (m, r theta) plane is the chord frame scaled by c' r: centroid and area map exactly
        const auto c = polygon_centroid(L.section->loop.xy);
        const auto conf = placement(d.layers[k].stagger_deg, d.layers[k].theta_le_deg, L.surface.mprime_te(), d.mirror);
        const Vec2 cg = conf(c.centroid);
        CHECK(L.cg_mprime == Approx(cg.x()).margin(1e-12));
        CHECK(L.cg_theta == Approx(cg.y()).margin(1e-12));
        const double s = L.chord_mprime * L.surface.radius();
        CHECK(L.section_area == Approx(c.area * s * s).epsilon(1e-12));
        CHECK(std::hypot(L.cg_xyz.y(), L.cg_xyz.z()) == Approx(L.surface.radius()).epsilon(1e-12));
        CHECK(L.cg_xr.x() == Approx(L.cg_xyz.x()).margin(1e-12));
    }
    // a re-stagger about the radial axis through the mid layer's centre of gravity leaves it in place
    const auto& M = r.layers[2];
    auto e = d;
    e.restagger = {12.0, Restagger::Axis::Point, 0.0, M.cg_xr.x(), Math::rad2deg(M.cg_theta)};
    const auto q = ev.evaluate(e);
    REQUIRE(q.ok());
    CHECK((q.layers[2].cg_xyz - M.cg_xyz).norm() < 1e-9);
    CHECK((q.layers[0].cg_xyz - r.layers[0].cg_xyz).norm() > 1e-3);  // the other layers' do move
}

TEST_CASE("Blade centre of gravity: the closed mesh against the exact OCCT solid", "[design][cg]")
{
    const auto d = BladeDesign::defaults();
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    REQUIRE(r.blade_volume > 0.0);
    // the exact solid (trimmed by the hub and casing surfaces), metres -> mm
    const auto s = build_design_solid(d, ev, 1e-3, 25, 0.0, true);
    GProp_GProps g;
    BRepGProp::VolumeProperties(s.blade, g, 1e-7);
    const double v_exact = g.Mass() * 1e9;
    const gp_Pnt c = g.CentreOfMass();
    const Vec3 cg_exact(1e3 * c.X(), 1e3 * c.Y(), 1e3 * c.Z());
    INFO("volume mesh " << r.blade_volume << " exact " << v_exact << " mm^3; CG gap " << (r.blade_cg - cg_exact).norm() << " mm");
    CHECK(r.blade_volume == Approx(v_exact).epsilon(1e-2));  // flat caps on curved walls
    CHECK((r.blade_cg - cg_exact).norm() < 0.05);            // mm (measured 0.016 mm; chord ~45 mm, span ~60 mm)
    // the centre of gravity lies between hub and casing and between the LE and TE lines
    const double rad = std::hypot(r.blade_cg.y(), r.blade_cg.z());
    CHECK(rad > r.layers.front().r_le);
    CHECK(rad < r.layers.back().r_le);
    CHECK(r.blade_cg.x() > d.meridional.hub.z_le);
    CHECK(r.blade_cg.x() < d.meridional.hub.z_te);
}
