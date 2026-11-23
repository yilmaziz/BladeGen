// Cascade.hpp: cascade generation from a blade section (port of v1 geometry.py's cascade part),
// the geometric throat and the not-a-knot cubic spline behind it.
#include "CascadeFixture.hpp"

#include "turbolenza/Cascade.hpp"
#include "turbolenza/CubicSpline.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <iostream>

namespace T = BladeGen::Turbolenza;

namespace {

double max_diff(const std::vector<T::Point>& a, const std::vector<T::Point>& b)
{
    REQUIRE(a.size() == b.size());
    double d = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) d = std::max({d, std::abs(a[k].m - b[k].m), std::abs(a[k].y - b[k].y)});
    return d;
}

T::Point pt(const ref::json& a)
{
    const auto v = ref::vec(a);
    return {v[0], v[1]};
}

}  // namespace

TEST_CASE("the not-a-knot cubic spline reproduces a cubic and the scipy end cases", "[turbolenza][cascade]")
{
    // not-a-knot reproduces every cubic exactly (Boor, A Practical Guide to Splines, Ch. IV)
    auto f = [](double x) { return T::Point{1.0 - 2.0 * x + 0.5 * x * x * x, 0.3 * x * x - x * x * x}; };
    auto f1 = [](double x) { return T::Point{-2.0 + 1.5 * x * x, 0.6 * x - 3.0 * x * x}; };
    std::vector<double> x{0.0, 0.1, 0.35, 0.4, 0.8, 1.3, 1.31, 2.0};
    std::vector<T::Point> y;
    for (double v : x) y.push_back(f(v));
    const T::detail::CubicSpline s(x, y);
    for (double v : {0.0, 0.05, 0.2, 0.37, 1.0, 1.305, 1.7, 2.0}) {
        CHECK(s(v).m == Catch::Approx(f(v).m).margin(1e-13));
        CHECK(s(v).y == Catch::Approx(f(v).y).margin(1e-13));
        CHECK(s(v, 1).m == Catch::Approx(f1(v).m).margin(1e-12));
        CHECK(s(v, 1).y == Catch::Approx(f1(v).y).margin(1e-12));
        CHECK(s(v, 3).m == Catch::Approx(3.0).margin(1e-9));
    }
    // n = 3: the parabola through the points; n = 2: the line
    const T::detail::CubicSpline p({0.0, 1.0, 3.0}, {{0.0, 1.0}, {1.0, 0.0}, {9.0, 4.0}});
    CHECK(p(2.0).m == Catch::Approx(4.0).margin(1e-13));
    CHECK(p(2.0).y == Catch::Approx(1.0).margin(1e-13));  // y = (x - 1)^2, m = x^2
    const T::detail::CubicSpline l({0.0, 2.0}, {{0.0, 0.0}, {2.0, 1.0}});
    CHECK(l(0.5).m == Catch::Approx(0.5));
    CHECK(l(0.5).y == Catch::Approx(0.25));
    CHECK_THROWS_AS(T::detail::CubicSpline({0.0, 0.0, 1.0}, {{0, 0}, {1, 1}, {2, 2}}), std::invalid_argument);
    CHECK_THROWS_AS(T::detail::CubicSpline({0.0}, {{0, 0}}), std::invalid_argument);
}

TEST_CASE("closest_points finds the gap between two circles", "[turbolenza][cascade]")
{
    // circles of radius 1 about (0, 0) and (3, 0.5): gap |c2 - c1| - 2 along the centre line
    std::vector<T::Point> a, b;
    for (int k = 0; k <= 60; ++k) {
        const double t = -1.0 + 2.0 * k / 60.0;  // arcs facing each other
        a.push_back({std::cos(t), std::sin(t)});
        b.push_back({3.0 - std::cos(t), 0.5 + std::sin(t)});
    }
    const auto c = T::closest_points(a, b);
    const double d = std::hypot(3.0, 0.5);
    // the splines through 61 points of a circle are off the circle by ~1e-7 (h^4)
    CHECK(c.distance == Catch::Approx(d - 2.0).margin(1e-6));
    CHECK(std::atan2(c.a.y, c.a.m) == Catch::Approx(std::atan2(0.5, 3.0)).margin(1e-5));
    CHECK(std::hypot(c.a.m - c.b.m, c.a.y - c.b.y) == Catch::Approx(c.distance).margin(1e-15));
}

TEST_CASE("make_cascade refuses invalid input and overlapping blades", "[turbolenza][cascade]")
{
    const auto s = fixture::synthetic_section();
    CHECK_THROWS_AS(T::make_cascade(s, 0.0, 30.0), std::invalid_argument);
    CHECK_THROWS_AS(T::make_cascade(s, -1.0, 30.0), std::invalid_argument);
    CHECK_THROWS_AS(T::make_cascade(s, std::nan(""), 30.0), std::invalid_argument);
    CHECK_THROWS_AS(T::make_cascade(s, 1.0, 90.0), std::invalid_argument);
    CHECK_THROWS_AS(T::make_cascade(s, 1.0, -90.0), std::invalid_argument);
    // pitch 0.05 < the blade's 0.09 thickness at zero stagger: the neighbours overlap
    CHECK_THROWS_AS(T::make_cascade(s, 20.0, 0.0), T::GeometryError);
    // sample parameters that do not start at u = 0
    auto bad = s;
    bad.suction_parameters.front() = 0.001;
    CHECK_THROWS_AS(T::make_cascade(bad, 1.0, 30.0), T::GeometryError);
    // a blade whose suction side is on -y (surfaces swapped)
    auto swapped = s;
    std::swap(swapped.suction, swapped.pressure);
    CHECK_THROWS_AS(T::make_cascade(swapped, 1.0, 30.0), T::GeometryError);
}

TEST_CASE("make_cascade of the synthetic blade: frame, curves and throat", "[turbolenza][cascade]")
{
    const auto s = fixture::synthetic_section();
    const double stagger_deg = 30.0;
    const auto c = T::make_cascade(s, 1.25, stagger_deg);
    const double st = stagger_deg * fixture::kPi / 180.0;
    CHECK(c.stagger == st);
    CHECK(c.active.pitch == 0.8);
    CHECK(c.active.chord == 1.0);
    CHECK(c.kappa_le == st + 0.3);
    CHECK(c.kappa_te == st - 0.25);
    CHECK(c.camber_angle == Catch::Approx(0.55));
    CHECK(c.active.leading_edge_point.m == 0.0);
    CHECK(c.active.trailing_edge_point.m == Catch::Approx(std::cos(st)).margin(1e-15));
    CHECK(c.active.trailing_edge_point.y == Catch::Approx(std::sin(st)).margin(1e-15));
    // the cascade-frame curves are the blade-frame ones turned: X(u) = R X_blade(u)
    for (double u : {0.0, 0.13, 0.5, 0.77, 1.0}) {
        const auto a = c.active_curves->suction.evaluate(u);
        const auto b = T::rotate_by_stagger(s.suction.evaluate(u), st);
        CHECK(std::abs(a.m - b.m) <= 1e-15);
        CHECK(std::abs(a.y - b.y) <= 1e-15);
    }
    // the throat: a suction point and a pressure point of the neighbour, no closer pair of samples
    const auto& t = c.active;
    CHECK(t.throat_opening > 0.0);
    CHECK(t.throat_opening < t.pitch);
    double nearest = 1e9;
    for (const auto& p : t.suction_surface)
        for (const auto& q : t.pressure_surface) nearest = std::min(nearest, std::hypot(p.m - q.m, p.y - (q.y + t.pitch)));
    CHECK(t.throat_opening <= nearest + 1e-12);
    CHECK(t.throat_opening >= nearest - 2e-3);  // the samples are 1/80 apart in u
    CHECK(c.throat_over_pitch == t.throat_opening / t.pitch);
}

TEST_CASE("make_cascade matches v1 for the NACA 65-(12A10)10 at 33 deg", "[turbolenza][cascade][v1ref]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    const T::Cascade& c = R->cascade;
    const ref::json& v = R->family.at("cascade");
    const ref::json& a = v.at("active");
    CHECK(c.stagger == ref::num(v.at("stagger")));
    CHECK(c.kappa_le == ref::num(v.at("kappa_le")));
    CHECK(c.kappa_te == ref::num(v.at("kappa_te")));
    CHECK(c.camber_angle == Catch::Approx(ref::num(v.at("camber_angle"))).margin(1e-15));
    CHECK(c.active.chord == ref::num(a.at("chord")));
    CHECK(c.active.pitch == ref::num(a.at("pitch")));
    const double d_su = max_diff(c.active.suction_surface_unstaggered, ref::points(a.at("suction_surface_unstaggered")));
    const double d_pu = max_diff(c.active.pressure_surface_unstaggered, ref::points(a.at("pressure_surface_unstaggered")));
    const double d_s = max_diff(c.active.suction_surface, ref::points(a.at("suction_surface")));
    const double d_p = max_diff(c.active.pressure_surface, ref::points(a.at("pressure_surface")));
    const double d_le = max_diff({c.active.leading_edge_point, c.active.trailing_edge_point},
                                 {pt(a.at("leading_edge_point")), pt(a.at("trailing_edge_point"))});
    const double d_throat = std::abs(c.active.throat_opening - ref::num(a.at("throat_opening")));
    const double d_tp = max_diff({c.active.throat_point_reference, c.active.throat_point_neighbor},
                                 {pt(a.at("throat_point_reference")), pt(a.at("throat_point_neighbor"))});
    const auto cv = ref::surface_curves(v.at("active_curves"));
    const double d_cs = max_diff(c.active_curves->suction.control_points(), cv.suction.control_points());
    const double d_cp = max_diff(c.active_curves->pressure.control_points(), cv.pressure.control_points());
    INFO("cascade vs v1: surfaces unstaggered " << d_su << " " << d_pu << ", staggered " << d_s << " " << d_p
                                                << ", LE/TE " << d_le << ", throat " << d_throat << ", throat points " << d_tp
                                                << ", curves " << d_cs << " " << d_cp);
    // measured (MSVC 19.51): every surface sample, landmark and curve bit for bit (0); the throat
    // opening and its end points 1.1e-16 (the not-a-knot splines and the Newton refinement against
    // scipy's CubicSpline and LAPACK).  Bounds: 1e-15 where 0 was measured (an ulp of slack for
    // another compiler's libm), about 10x the measurement otherwise.
    CHECK(d_su <= 1e-15);
    CHECK(d_pu <= 1e-15);
    CHECK(d_s <= 1e-15);
    CHECK(d_p <= 1e-15);
    CHECK(d_le <= 1e-15);
    CHECK(d_throat <= 2e-15);
    CHECK(d_tp <= 2e-15);
    CHECK(d_cs <= 1e-15);
    CHECK(d_cp <= 1e-15);
    CHECK(std::abs(c.throat_over_pitch - ref::num(v.at("throat_over_pitch"))) <= 2e-15);
}
