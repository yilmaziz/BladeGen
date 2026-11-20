// Tests of blade/Meridional and blade/Blade against closed forms.
#include "blade/Blade.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Blade;
using Catch::Approx;
using Math::Vec2;
using Math::Vec3;

TEST_CASE("m' of a conical streamline has the closed form L ln(r/r0) / (r1 - r0)", "[meridional]")
{
    const Vec2 a(0.0, 0.2), b(0.3, 0.35);
    const auto sl = MeridionalCurve::line(a, b);
    const double L = (b - a).norm();
    REQUIRE(sl.length() == Approx(L).epsilon(1e-13));
    auto mp = [&](double m) { return L / (b.y() - a.y()) * std::log((a.y() + (b.y() - a.y()) * m / L) / a.y()); };
    REQUIRE(sl.mprime_length() == Approx(mp(L)).epsilon(1e-12));
    for (double m : {0.0, 0.01, 0.1, 0.2, L}) {
        REQUIRE(sl.mprime_of_m(m) == Approx(mp(m)).margin(1e-13));
        REQUIRE(sl.m_of_mprime(sl.mprime_of_m(m)) == Approx(m).margin(1e-12));
    }
    // straight extrapolation past both ends keeps the same closed form
    for (double m : {-0.05, L + 0.04}) {
        REQUIRE(sl.mprime_of_m(m) == Approx(mp(m)).margin(1e-12));
        REQUIRE(sl.m_of_mprime(sl.mprime_of_m(m)) == Approx(m).margin(1e-12));
        const Vec2 p = sl.at_m(m);
        REQUIRE((p - (a + (b - a) * m / L)).norm() < 1e-12);
    }
}

TEST_CASE("m' on a cylinder is m / r", "[meridional]")
{
    const auto sl = MeridionalCurve::line({0.0, 0.5}, {0.2, 0.5});
    REQUIRE(sl.mprime_length() == Approx(0.4).epsilon(1e-13));
    REQUIRE(sl.mprime_of_m(-0.1) == Approx(-0.2).epsilon(1e-12));
}

TEST_CASE("a curved streamline: m is the arc length and m' integrates dm / r", "[meridional]")
{
    std::vector<Vec2> pts;
    for (int i = 0; i <= 20; ++i) {
        const double x = 0.01 * i;
        pts.emplace_back(x, 0.3 + 0.2 * x * x);
    }
    const auto sl = MeridionalCurve::through_points(pts);
    // compare with a fine polyline of the same curve
    double m = 0.0, mp = 0.0;
    Vec2 prev = sl.curve().eval(sl.curve().t0());
    for (int i = 1; i <= 20000; ++i) {
        const Vec2 q = sl.curve().eval(sl.curve().t0() + (sl.curve().t1() - sl.curve().t0()) * i / 20000.0);
        const double ds = (q - prev).norm();
        m += ds;
        mp += ds / (0.5 * (q.y() + prev.y()));
        prev = q;
    }
    REQUIRE(sl.length() == Approx(m).epsilon(1e-8));
    REQUIRE(sl.mprime_length() == Approx(mp).epsilon(1e-8));
    REQUIRE_THROWS_AS(MeridionalCurve::line({0, 0.1}, {0, 0.1}), std::invalid_argument);
    REQUIRE_THROWS_AS(MeridionalCurve::line({0, -0.1}, {1, 0.1}), std::invalid_argument);
}

TEST_CASE("channel streamlines interpolate hub and shroud", "[meridional]")
{
    const auto ch = MeridionalChannel::straight(0.0, 0.1, 0.2, 0.22, 0.3, 0.29);
    const auto s0 = ch.streamline(0.0), s1 = ch.streamline(1.0), sh = ch.streamline(0.5);
    REQUIRE((s0.at_m(0.0) - Vec2(0.0, 0.2)).norm() < 1e-12);
    REQUIRE((s1.at_m(s1.length()) - Vec2(0.1, 0.29)).norm() < 1e-12);
    REQUIRE((sh.at_m(0.0) - Vec2(0.0, 0.25)).norm() < 1e-12);
    const auto below = ch.streamline(-0.1);
    REQUIRE(below.at_m(0.0).y() == Approx(0.19).epsilon(1e-12));
}

TEST_CASE("section placement: LE / TE points, chord and conformality", "[blade]")
{
    const auto sl = MeridionalCurve::line({0.0, 0.3}, {0.1, 0.3});  // cylinder r = 0.3
    SectionPlacement pl{0.2, 0.8, 35.0, 5.0, false};
    const auto map = make_section_map(sl, pl);
    const double g = Math::deg2rad(35.0);
    REQUIRE(map.mprime_le == Approx(0.02 / 0.3).epsilon(1e-12));
    REQUIRE(map.chord * std::cos(g) == Approx(0.06 / 0.3).epsilon(1e-12));
    const Vec2 te = map.conformal(Vec2(1.0, 0.0));
    REQUIRE(te.x() == Approx(0.08 / 0.3).epsilon(1e-12));
    REQUIRE(te.y() == Approx(Math::deg2rad(5.0) + map.chord * std::sin(g)).epsilon(1e-12));
    // on a cylinder the chord is a helix: its 3-D angle to the axis equals the stagger
    const Vec3 a = map.cartesian(map.conformal(Vec2(0.3, 0.0))), b = map.cartesian(map.conformal(Vec2(0.3001, 0.0)));
    const Vec3 d = (b - a).normalized();
    REQUIRE(std::acos(d.x()) == Approx(g).epsilon(1e-6));
    // mirrored placement puts the suction side on the other side of the chord
    SectionPlacement mi = pl;
    mi.mirror = true;
    const auto mm = make_section_map(sl, mi);
    REQUIRE((mm.conformal(Vec2(0.5, 0.05)) - map.conformal(Vec2(0.5, -0.05))).norm() < 1e-14);
    REQUIRE_THROWS_AS(make_section_map(sl, SectionPlacement{0.8, 0.2, 0.0, 0.0, false}), std::invalid_argument);
}

TEST_CASE("demo rotor: stations reproduced, sections closed and on their stream surfaces", "[blade]")
{
    const BladeGen::Blade::Blade blade(BladeDefinition::demo_rotor());
    const auto& st = blade.definition().stations;
    for (const auto& s : st) {
        const auto p = blade.profile_at(s.span);
        REQUIRE(p.le.beta_deg == Approx(s.profile.le.beta_deg).epsilon(1e-12));
        REQUIRE(blade.placement_at(s.span).stagger_deg == Approx(s.placement.stagger_deg).epsilon(1e-12));
    }
    const auto secs = blade.sections(20, 60);
    REQUIRE(secs.size() == static_cast<std::size_t>(blade.definition().n_sections + 2));
    for (const auto& s : secs) {
        REQUIRE((s.xyz.front() - s.xyz.back()).norm() < 1e-12);
        // every point lies on the stream surface: its (x, r) is on the streamline
        for (std::size_t i = 0; i < s.xyz.size(); i += 17) {
            const Vec3& q = s.xyz[i];
            const double r = std::hypot(q.y(), q.z());
            const Vec2 xr = s.streamline.at_m(s.streamline.m_of_mprime(s.conformal[i].x()));
            REQUIRE(std::abs(xr.x() - q.x()) < 1e-12);
            REQUIRE(std::abs(xr.y() - r) < 1e-12);
        }
    }
    REQUIRE(secs.front().span < 0.0);
    REQUIRE(secs.back().span > 1.0);
}

TEST_CASE("invalid blade definitions throw instead of misbehaving", "[blade]")
{
    auto d = BladeDefinition::demo_rotor();
    d.channel = MeridionalChannel{};  // no channel, no stream lines
    REQUIRE_THROWS_AS(d.validate(), std::invalid_argument);
    d = BladeDefinition::demo_rotor();
    d.stations[1].streamline = d.channel.streamline(0.5);  // only one station with its own stream line
    REQUIRE_THROWS_AS(d.validate(), std::invalid_argument);
    d = BladeDefinition::demo_rotor();
    d.stations[2].placement.mirror = true;
    REQUIRE_THROWS_AS(d.validate(), std::invalid_argument);
}
