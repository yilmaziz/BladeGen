// Tests of profile/Camber: inscribed-circle camber line and thickness against closed forms.
#include "naca/Naca65.hpp"
#include "profile/Camber.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Profile;
using Catch::Approx;
using Math::Vec2;

TEST_CASE("Camber of two concentric arcs: the mid arc, thickness 2h, angle = arc tangent", "[camber]")
{
    // centre (0.5, -Y), mid radius R; sides at R + h (suction, above) and R - h
    const double R = 1.2, h = 0.03, Y = std::sqrt(R * R - 0.25);
    const double a0 = std::atan2(Y, -0.5), a1 = std::atan2(Y, 0.5);  // angles of the chord ends
    auto arc = [&](double rad) {
        std::vector<Vec2> p;
        for (int i = 0; i <= 600; ++i) {
            const double a = a0 + (a1 - a0) * i / 600.0;
            p.emplace_back(0.5 + rad * std::cos(a), -Y + rad * std::sin(a));
        }
        return p;
    };
    // round edges: quarter circles of radius h around the mid-arc end points close each side
    // to the LE / TE apex on the mid-arc tangent (a profile of constant thickness 2h)
    auto side = [&](double rad, double sgn) {
        std::vector<Vec2> p;
        const Vec2 e0(0.0, 0.0), e1(1.0, 0.0);
        const Vec2 r0 = (e0 - Vec2(0.5, -Y)).normalized(), r1 = (e1 - Vec2(0.5, -Y)).normalized();
        const Vec2 t0(r0.y(), -r0.x()), t1(r1.y(), -r1.x());  // clockwise tangents (LE -> TE)
        for (int i = 0; i <= 40; ++i) {  // LE cap: apex e0 - h t0 -> arc start e0 + sgn h r0
            const double a = 0.5 * Math::kPi * i / 40.0;
            p.push_back(e0 + h * (-std::cos(a) * t0 + sgn * std::sin(a) * r0));
        }
        auto body = arc(rad);
        p.insert(p.end(), body.begin() + 1, body.end() - 1);
        for (int i = 0; i <= 40; ++i) {  // TE cap
            const double a = 0.5 * Math::kPi * i / 40.0;
            p.push_back(e1 + h * (std::sin(a) * t1 + sgn * std::cos(a) * r1));
        }
        return p;
    };
    const auto c = camber_line(side(R + h, 1.0), side(R - h, -1.0), 61);
    REQUIRE(c.point.size() == 61);
    // away from the sharp edges, where the inscribed circles are those of the two arcs
    int checked = 0;
    for (std::size_t i = 1; i + 1 < c.point.size(); ++i) {
        if (c.point[i].x() < 0.15 || c.point[i].x() > 0.85) continue;
        const Vec2 q = c.point[i] - Vec2(0.5, -Y);
        CHECK(q.norm() == Approx(R).margin(5e-5));
        CHECK(c.thickness[i] == Approx(2 * h).margin(4e-6));
        // the tangent of the arc at q is perpendicular to q
        CHECK(c.angle_deg[i] == Approx(Math::rad2deg(std::atan2(-q.x(), q.y()))).margin(0.05));
        ++checked;
    }
    CHECK(checked > 25);
    CHECK(c.max_thickness == Approx(2 * h).margin(4e-6));
}

TEST_CASE("Camber of a symmetric profile is the chord line", "[camber]")
{
    std::vector<Vec2> ss, ps;
    for (int i = 0; i <= 400; ++i) {
        const double x = (i / 400.0) * (i / 400.0);  // x = t^2 resolves the sqrt nose
        const double y = 0.6 * (0.2969 * std::sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1036 * x * x * x * x);
        ss.emplace_back(x, y);
        ps.emplace_back(x, -y);
    }
    const auto c = camber_line(ss, ps, 81);
    for (std::size_t i = 0; i < c.point.size(); ++i) {
        CHECK(c.point[i].y() == Approx(0.0).margin(1e-6));  // polyline discretisation of the sides
        CHECK(c.angle_deg[i] == Approx(0.0).margin(5e-3));
    }
    // NACA 4-digit with t = 0.12 (0.6 * 0.2): the inscribed circle is at most the vertical thickness
    CHECK(c.max_thickness == Approx(0.12).epsilon(0.02));
    CHECK(c.max_thickness_x == Approx(0.30).margin(0.03));
}

TEST_CASE("Camber of a NACA 65 section: close to its A10 mean line, t_max = t/c", "[camber]")
{
    const Naca::Naca65Section s({1.2, 0.10, Naca::Variant::KovachSandercock});
    const auto q = Naca::check_quality(s);
    const auto c = camber_line(q.ss.xy, q.ps.xy, 121);
    CHECK(c.max_thickness == Approx(0.10).epsilon(0.02));
    double dev = 0.0;
    for (std::size_t i = 1; i + 1 < c.point.size(); ++i) {
        const double x = c.point[i].x();
        if (x < 0.02 || x > 0.97) continue;
        dev = std::max(dev, std::abs(c.point[i].y() - s.camber_point(x).y()));
    }
    CHECK(dev < 2e-3);
    const auto eq = Naca::equivalent_arc_angles(1.2);
    CHECK(c.angle_deg.front() > 0.5 * Math::rad2deg(eq[1]));   // camber rises at the LE
    CHECK(c.angle_deg.back() < -0.5 * Math::rad2deg(eq[2]));   // and falls at the TE
}

TEST_CASE("Camber line: degenerate input is rejected; helpers", "[camber]")
{
    const std::vector<Vec2> few = {{0, 0}, {1, 0}};
    CHECK_THROWS_AS(camber_line(few, few), std::invalid_argument);
    std::vector<Vec2> a = {{0, 0}, {0.3, 0.1}, {0.6, 0.1}, {1, 0}};
    CHECK_THROWS_AS(camber_line(a, a, 3), std::invalid_argument);
    auto bad = a;
    bad[1].x() = std::nan("");
    CHECK_THROWS_AS(camber_line(bad, a), std::invalid_argument);
    CHECK_THROWS_AS(camber_line(a, a), std::invalid_argument);  // coincident sides: no circle
    CHECK(polygon_area({{0, 0}, {1, 0}, {1, 1}, {0, 1}}) == Approx(1.0));
    CHECK(polygon_area({{0, 0}, {0, 1}, {1, 1}, {1, 0}}) == Approx(-1.0));
    CHECK(polyline_distance({0.5, 1.0}, {{0, 0}, {1, 0}}) == Approx(1.0));
    CHECK(polyline_distance({2.0, 0.0}, {{0, 0}, {1, 0}}) == Approx(1.0));
    CHECK_THROWS_AS(polyline_distance({0, 0}, {}), std::invalid_argument);
}
