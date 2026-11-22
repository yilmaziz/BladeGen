// Tests of the re-stagger of the whole blade (design/BladeDesign.hpp, "Re-stagger"): every
// section turns rigidly in its conformal plane about the point where the radial axis crosses
// its stream surface, so the pivot stays put, every distance to it is kept, and every metal
// angle changes by exactly the re-stagger angle.
#include "design/BladeDesign.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Design;
using Catch::Approx;
using Math::Vec2;

namespace {

/// The pivot of layer k: where the axis (x_a, theta_a) crosses its stream surface.
Vec2 pivot_of(const LayerResult& L, const Vec2& axis) { return {mprime_at_x(L.surface, axis.x()), axis.y()}; }

/// Nominal vs re-staggered layer: rigid rotation by d about P (conformal plane), angles + d.
void check_rigid(const LayerResult& a, const LayerResult& b, const Vec2& P, double d_deg)
{
    const double d = Math::deg2rad(d_deg), c = std::cos(d), s = std::sin(d);
    REQUIRE(a.conformal.size() == b.conformal.size());
    for (std::size_t i = 0; i < a.conformal.size(); ++i) {
        const Vec2 v = a.conformal[i] - P;
        const Vec2 rotated = P + Vec2(c * v.x() - s * v.y(), s * v.x() + c * v.y());
        CHECK((b.conformal[i] - rotated).norm() < 1e-12);
        CHECK((b.conformal[i] - P).norm() == Approx((a.conformal[i] - P).norm()).margin(1e-12));
    }
    REQUIRE(a.beta_deg.size() == b.beta_deg.size());
    for (std::size_t i = 0; i < a.beta_deg.size(); ++i) CHECK(b.beta_deg[i] == Approx(a.beta_deg[i] + d_deg).margin(1e-11));
    CHECK(b.chord_mprime == Approx(a.chord_mprime).epsilon(1e-14));
    CHECK(Math::rad2deg(b.gamma) == Approx(Math::rad2deg(a.gamma) + d_deg).margin(1e-12));
}

}  // namespace

TEST_CASE("Placement: the similarity map, re-stagger about a conformal point", "[design][restagger]")
{
    const auto n = placement(30.0, 2.0, 0.8, false);
    CHECK(n.m0 == 0.0);
    CHECK(n.t0 == Approx(Math::deg2rad(2.0)).margin(1e-15));
    CHECK(n.c == Approx(0.8 / std::cos(Math::deg2rad(30.0))).epsilon(1e-15));
    CHECK(n({1.0, 0.0}).x() == Approx(0.8).margin(1e-15));  // TE on the TE line: c' cos(gamma) = m'_TE
    const Vec2 P(0.3, 0.1);
    for (bool mirror : {false, true}) {
        const auto a = placement(30.0, 2.0, 0.8, mirror), b = placement(30.0, 2.0, 0.8, mirror, -7.0, P);
        const double d = Math::deg2rad(-7.0);
        for (const Vec2& q : {Vec2(0.0, 0.0), Vec2(1.0, 0.0), Vec2(0.4, 0.06), Vec2(0.7, -0.03)}) {
            const Vec2 v = a(q) - P;
            const Vec2 r = P + Vec2(std::cos(d) * v.x() - std::sin(d) * v.y(), std::sin(d) * v.x() + std::cos(d) * v.y());
            CHECK((b(q) - r).norm() < 1e-14);
        }
    }
    // a pivot on the chord line stays where it is
    const Vec2 Q = n({0.25, 0.0});
    CHECK((placement(30.0, 2.0, 0.8, false, 12.0, Q)({0.25, 0.0}) - Q).norm() < 1e-14);
}

TEST_CASE("Re-stagger about the hub LE (rotor): pivot fixed, angles + d on every layer, mesh consistent", "[design][restagger]")
{
    auto d = BladeDesign::defaults();
    Evaluator ev;
    const auto a = ev.evaluate(d);
    REQUIRE(a.ok());
    d.restagger.angle_deg = 4.0;  // closed by 4 deg
    const auto b = ev.evaluate(d);
    REQUIRE(b.ok());
    const Blade::FlowPath fp(d.meridional);
    const Vec2 axis = restagger_axis(d, fp);
    // the axis passes through the hub LE point
    CHECK(axis.x() == Approx(d.meridional.hub.z_le).margin(1e-12));
    CHECK(axis.y() == Approx(Math::deg2rad(d.layers[0].theta_le_deg)).margin(1e-15));
    CHECK(b.layers[0].mprime_le == Approx(0.0).margin(1e-12));
    CHECK(b.layers[0].theta_le == Approx(a.layers[0].theta_le).margin(1e-15));
    for (std::size_t k = 0; k < d.layers.size(); ++k) check_rigid(a.layers[k], b.layers[k], pivot_of(a.layers[k], axis), 4.0);
    // closing turns the sections about the axis: the TE leaves the TE line (here by more than
    // 0.5 mm axially), the chord fractions still run 0 -> 1
    CHECK(std::abs(b.layers[4].camber_xr.back().x() - a.layers[4].camber_xr.back().x()) > 0.5);
    CHECK(b.layers[4].mfrac.front() == Approx(0.0).margin(1e-12));
    CHECK(b.layers[4].mfrac.back() == Approx(1.0).margin(1e-9));
    // the mesh rows at the layer spans are the placed layers (rows interpolate exactly at the layers)
    const auto& L0 = b.layers[0];
    for (int i = 0; i < b.blade.ni; i += 5) CHECK((b.blade.at(i, 0) - L0.xyz[static_cast<std::size_t>(i)]).norm() < 1e-9);
    const auto& L4 = b.layers[4];
    for (int i = 0; i < b.blade.ni; i += 5) CHECK((b.blade.at(i, b.blade.nj - 1) - L4.xyz[static_cast<std::size_t>(i)]).norm() < 1e-9);
    // the metal-angle field: at the placed LE of the hub row, beta_LE + d (inside the passage the
    // field is sampled at meridional fractions, and turning a section moves its chord points in m')
    REQUIRE(a.beta.value.size() == b.beta.value.size());
    CHECK(b.beta.value.front() == Approx(a.layers[0].beta_le + 4.0).margin(1e-9));
}

TEST_CASE("Re-stagger about the mid span, the tip, a chord point and an arbitrary point", "[design][restagger]")
{
    auto base = BladeDesign::defaults();
    Evaluator ev;
    const auto a = ev.evaluate(base);
    REQUIRE(a.ok());
    const Blade::FlowPath fp(base.meridional);
    struct Case
    {
        Restagger::Axis axis;
        double xp, z, th;
    };
    for (const auto& c : {Case{Restagger::Axis::Mid, 0.5, 0, 0}, Case{Restagger::Axis::Tip, 0.0, 0, 0}, Case{Restagger::Axis::Hub, 1.0, 0, 0},
                          Case{Restagger::Axis::Point, 0.0, 2.5, -3.0}}) {
        auto d = base;
        d.restagger = {-6.0, c.axis, c.xp, c.z, c.th};  // opened by 6 deg
        const auto b = ev.evaluate(d);
        REQUIRE(b.ok());
        const Vec2 axis = restagger_axis(d, fp);
        for (std::size_t k = 0; k < d.layers.size(); ++k) check_rigid(a.layers[k], b.layers[k], pivot_of(a.layers[k], axis), -6.0);
        if (c.axis == Restagger::Axis::Point) {
            CHECK(axis.x() == 2.5);
            CHECK(axis.y() == Approx(Math::deg2rad(-3.0)).margin(1e-15));
        }
        if (c.axis == Restagger::Axis::Tip) {  // the tip LE stays
            CHECK(b.layers[4].mprime_le == Approx(0.0).margin(1e-12));
            CHECK(b.layers[4].theta_le == Approx(a.layers[4].theta_le).margin(1e-14));
        }
        if (c.axis == Restagger::Axis::Hub && c.xp == 1.0) {  // the hub TE stays
            CHECK(b.layers[0].camber_xr.back().x() == Approx(a.layers[0].camber_xr.back().x()).margin(1e-9));
            CHECK(b.layers[0].theta_deg.back() == Approx(a.layers[0].theta_deg.back()).margin(1e-9));
        }
    }
}

TEST_CASE("Re-stagger: invalid input, JSON round trip", "[design][restagger]")
{
    auto d = BladeDesign::defaults();
    d.restagger.angle_deg = 31.0;
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    d.restagger.angle_deg = std::nan("");
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    d.restagger = {};
    d.restagger.chord_fraction = 1.5;
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    d.restagger = {};
    d.layers[4].stagger_deg = 78.0;
    d.restagger.angle_deg = 5.0;  // 83 deg in total
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    // an axis outside the channel is a design error (reported, not thrown)
    d = BladeDesign::defaults();
    d.restagger = {3.0, Restagger::Axis::Point, 0.0, 1e4, 0.0};
    Evaluator ev;
    const auto r = ev.evaluate(d);
    CHECK_FALSE(r.ok());
    REQUIRE_FALSE(r.errors.empty());
    CHECK(r.errors.front().find("re-stagger") != std::string::npos);
    // zero angle: the axis is not needed
    d.restagger.angle_deg = 0.0;
    CHECK(ev.evaluate(d).ok());
    // JSON
    d.restagger = {-2.5, Restagger::Axis::Point, 0.3, 4.0, 1.5};
    const auto e = design_from_json(to_json(d));
    CHECK(e.restagger.angle_deg == -2.5);
    CHECK(e.restagger.axis == Restagger::Axis::Point);
    CHECK(e.restagger.chord_fraction == 0.3);
    CHECK(e.restagger.z == 4.0);
    CHECK(e.restagger.theta_deg == 1.5);
    auto j = to_json(d);
    j.erase("restagger");  // designs saved before the option: no re-stagger
    CHECK(design_from_json(j).restagger.angle_deg == 0.0);
    j = to_json(d);
    j["restagger"]["axis"] = "shaft";
    CHECK_THROWS_AS(design_from_json(j), std::invalid_argument);
    for (auto a : {Restagger::Axis::Hub, Restagger::Axis::Mid, Restagger::Axis::Tip, Restagger::Axis::Point})
        CHECK(restagger_axis_from_string(to_string(a)) == a);
}
