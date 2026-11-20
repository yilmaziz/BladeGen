// Tests of blade/FlowPath: the Bell & Mehta wall law, the exact wall B-splines, layer
// streamlines, the conformal map (cylindrical wrap and surface of revolution).
#include "blade/FlowPath.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Blade;
using Catch::Approx;

TEST_CASE("Bell & Mehta fifth-order contour: ends, symmetry, zero slope and curvature at both ends", "[flowpath]")
{
    CHECK(bell_mehta(0.0) == 0.0);
    CHECK(bell_mehta(1.0) == Approx(1.0).margin(1e-15));
    CHECK(bell_mehta(0.5) == Approx(0.5).margin(1e-15));
    for (double x : {0.1, 0.3, 0.7}) CHECK(bell_mehta(x) + bell_mehta(1.0 - x) == Approx(1.0).margin(1e-14));
    CHECK(bell_mehta_d1(0.0) == 0.0);
    CHECK(bell_mehta_d1(1.0) == 0.0);
    CHECK(bell_mehta_d2(0.0) == 0.0);
    CHECK(bell_mehta_d2(1.0) == 0.0);
    const double h = 1e-6;
    for (double x : {0.2, 0.5, 0.8}) {
        CHECK(bell_mehta_d1(x) == Approx((bell_mehta(x + h) - bell_mehta(x - h)) / (2 * h)).epsilon(1e-8));
        CHECK(bell_mehta_d2(x) == Approx((bell_mehta_d1(x + h) - bell_mehta_d1(x - h)) / (2 * h)).margin(1e-6));
    }
}

TEST_CASE("Bell & Mehta wall: the exact quintic between LE and TE, cylinders outside, G2 joins", "[flowpath]")
{
    const WallPoints w{-15.0, 40.0, 15.0, 60.0};
    const auto c = wall_spline(w, WallLaw::BellMehta, 0.5, 0.75);
    CHECK(c.t0() == 0.0);
    CHECK(c.t1() == 3.0);
    for (int i = 0; i <= 50; ++i) {
        const double u = 1.0 + i / 50.0;
        const auto p = c.eval(u);
        const double xi = (p.x() - w.z_le) / (w.z_te - w.z_le);
        CHECK(xi == Approx(u - 1.0).margin(1e-13));  // z linear in u
        CHECK(p.y() == Approx(w.r_le + (w.r_te - w.r_le) * bell_mehta(xi)).margin(1e-12));
    }
    CHECK(c.eval(0.0).x() == Approx(-15.0 - 0.5 * 30.0));
    CHECK(c.eval(3.0).x() == Approx(15.0 + 0.75 * 30.0));
    for (double u : {0.0, 0.4, 0.9}) CHECK(c.eval(u).y() == Approx(40.0).margin(1e-12));
    for (double u : {2.1, 2.6, 3.0}) CHECK(c.eval(u).y() == Approx(60.0).margin(1e-12));
    for (double u : {1.0, 2.0}) {
        // slope and curvature continuous across the joins (both are zero there)
        const auto l = c.derivs(u - 1e-12, 2), r = c.derivs(u + 1e-12, 2);
        CHECK(l[1].y() == Approx(0.0).margin(1e-9));
        CHECK(r[1].y() == Approx(0.0).margin(1e-9));
        CHECK(Math::signed_curvature(l[1], l[2]) == Approx(0.0).margin(1e-9));
        CHECK(Math::signed_curvature(r[1], r[2]) == Approx(0.0).margin(1e-9));
    }
}

TEST_CASE("Linear wall: the straight line through the LE and TE points, continued straight", "[flowpath]")
{
    const WallPoints w{-10.0, 50.0, 20.0, 62.0};
    const auto c = wall_spline(w, WallLaw::Linear, 0.5, 0.5);
    const double k = (w.r_te - w.r_le) / (w.z_te - w.z_le);
    for (int i = 0; i <= 60; ++i) {
        const auto p = c.eval(3.0 * i / 60.0);
        CHECK(p.y() == Approx(w.r_le + k * (p.x() - w.z_le)).margin(1e-12));
    }
    CHECK((c.eval(1.0) - Math::Vec2(w.z_le, w.r_le)).norm() < 1e-13);
    CHECK((c.eval(2.0) - Math::Vec2(w.z_te, w.r_te)).norm() < 1e-13);
}

TEST_CASE("Meridional configuration: invalid input is rejected", "[flowpath]")
{
    MeridionalConfig ok;
    REQUIRE_NOTHROW(ok.validate());
    auto bad = ok;
    bad.hub.z_te = bad.hub.z_le;
    CHECK_THROWS_AS(FlowPath(bad), std::invalid_argument);
    bad = ok;
    bad.casing.r_le = 30.0;  // below the hub
    CHECK_THROWS_AS(FlowPath(bad), std::invalid_argument);
    bad = ok;
    bad.hub.r_le = -1.0;
    CHECK_THROWS_AS(bad.validate(), std::invalid_argument);
    bad = ok;
    bad.inlet_extension = 0.0;
    CHECK_THROWS_AS(bad.validate(), std::invalid_argument);
    bad = ok;
    bad.hub_law = WallLaw::Linear;
    bad.hub = {-15.0, 5.0, 15.0, 60.0};  // steep line: r < 0 upstream
    CHECK_THROWS_AS(bad.validate(), std::invalid_argument);
    bad = ok;
    bad.hub.r_te = std::nan("");
    CHECK_THROWS_AS(bad.validate(), std::invalid_argument);
    CHECK_THROWS_AS(wall_law_from_string("spline"), std::invalid_argument);
    CHECK(wall_law_from_string(to_string(WallLaw::BellMehta)) == WallLaw::BellMehta);
    CHECK(wall_law_from_string(to_string(WallLaw::Linear)) == WallLaw::Linear);
}

TEST_CASE("Layer streamlines: exact blend, LE / TE on the straight LE / TE lines", "[flowpath]")
{
    MeridionalConfig cfg;
    cfg.hub = {-12.0, 40.0, 15.0, 60.0};
    cfg.casing = {-15.0, 80.0, 18.0, 78.0};
    cfg.casing_law = WallLaw::Linear;
    const FlowPath fp(cfg);
    for (double s : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto sl = fp.streamline(s);
        CHECK((sl.eval(1.0) - fp.le_point(s)).norm() < 1e-12);
        CHECK((sl.eval(2.0) - fp.te_point(s)).norm() < 1e-12);
        for (double u : {0.3, 1.4, 2.7}) {
            const Math::Vec2 p = sl.eval(u), q = (1 - s) * fp.hub().eval(u) + s * fp.casing().eval(u);
            CHECK((p - q).norm() < 1e-12);
        }
        const auto S = fp.surface(s);
        CHECK((S.le_point() - fp.le_point(s)).norm() < 1e-12);
        CHECK(S.xr(0.0).isApprox(fp.le_point(s), 1e-12));
        CHECK((S.xr(S.mprime_te()) - fp.te_point(s)).norm() < 1e-9);
    }
    CHECK_THROWS_AS(fp.streamline(std::nan("")), std::invalid_argument);
}

TEST_CASE("Cylindrical wrapping: constant radius gives x = x_LE + r m', r = const", "[flowpath]")
{
    MeridionalConfig cfg;
    cfg.hub = {-15.0, 50.0, 15.0, 50.0};
    cfg.casing = {-15.0, 80.0, 15.0, 80.0};
    const FlowPath fp(cfg);
    for (double s : {0.0, 0.5, 1.0}) {
        const auto S = fp.surface(s);
        REQUIRE(S.cylindrical());
        const double r = 50.0 + 30.0 * s;
        CHECK(S.radius() == Approx(r));
        CHECK(S.mprime_te() == Approx(30.0 / r).epsilon(1e-13));
        for (double mp : {-0.2, 0.0, 0.1, 0.3}) {
            CHECK(S.xr(mp).x() == Approx(-15.0 + r * mp).margin(1e-12));
            CHECK(S.xr(mp).y() == Approx(r));
            CHECK(S.m_of_mprime(mp) == Approx(r * mp).margin(1e-12));
        }
        // the unrolled cylinder: (r m', r theta) is the true surface distance
        const auto a = S.cartesian(0.1, 0.0), b = S.cartesian(0.1, 0.02);
        CHECK((a - b).norm() == Approx(2.0 * r * std::sin(0.01)).epsilon(1e-12));
    }
}

TEST_CASE("Surface of revolution: m' of a cone in closed form, Hermite map vs the adaptive integral", "[flowpath]")
{
    MeridionalConfig cfg;
    cfg.hub = {-15.0, 40.0, 15.0, 60.0};
    cfg.hub_law = WallLaw::Linear;  // a cone between LE and TE
    cfg.casing = {-15.0, 80.0, 15.0, 80.0};
    const FlowPath fp(cfg);
    const auto S = fp.surface(0.0);
    REQUIRE_FALSE(S.cylindrical());
    // cone: r = r_LE + m sin(phi), m' = ln(r / r_LE) / sin(phi)
    const double L = std::hypot(30.0, 20.0), sphi = 20.0 / L;
    CHECK(S.mprime_te() == Approx(std::log(60.0 / 40.0) / sphi).epsilon(1e-11));
    CHECK(S.m_te() == Approx(L).epsilon(1e-12));
    // a Bell & Mehta layer against the arc-length tables of MeridionalCurve
    const FlowPath fb(MeridionalConfig{});
    for (double s : {0.0, 0.6}) {
        const auto T = fb.surface(s);
        const MeridionalCurve mc(fb.streamline(s));
        // m of the LE point along mc: integrate up to u = 1
        double m_le = 0.0;
        {
            double lo = 0.0, hi = mc.length();
            for (int it = 0; it < 80; ++it) {
                const double mid = 0.5 * (lo + hi);
                (mc.at_m(mid).x() < T.le_point().x() ? lo : hi) = mid;
            }
            m_le = 0.5 * (lo + hi);
        }
        const double mp_le = mc.mprime_of_m(m_le);
        for (double f : {0.0, 0.2, 0.5, 0.9, 1.0}) {
            const double mp = f * T.mprime_te();
            const auto exact = mc.at_m(mc.m_of_mprime(mp_le + mp));
            CHECK((T.xr(mp) - exact).norm() < 1e-6);  // ~1e-8 relative to r
            CHECK(T.m_of_mprime(mp) == Approx(mc.m_of_mprime(mp_le + mp) - m_le).margin(1e-6));
        }
    }
}

TEST_CASE("Conformal map: an angle drawn in (m', theta) is the true angle on the surface", "[flowpath]")
{
    const FlowPath fp(MeridionalConfig{});
    const auto S = fp.surface(0.3);
    for (double beta_deg : {-30.0, 0.0, 25.0, 55.0})
        for (double mp : {0.05, 0.2, 0.35}) {
            const double b = Math::deg2rad(beta_deg), th = 0.1, h = 1e-6;
            const Math::Vec3 t = (S.cartesian(mp + h * std::cos(b), th + h * std::sin(b)) -
                                  S.cartesian(mp - h * std::cos(b), th - h * std::sin(b)));
            const Math::Vec3 em = (S.cartesian(mp + h, th) - S.cartesian(mp - h, th)).normalized();
            const Math::Vec3 et = (S.cartesian(mp, th + h) - S.cartesian(mp, th - h)).normalized();
            CHECK(Math::rad2deg(std::atan2(t.dot(et), t.dot(em))) == Approx(beta_deg).margin(1e-5));
        }
}

TEST_CASE("Streamline curvature: zero on the cylinders, Bell & Mehta extremes inside the blade", "[flowpath]")
{
    const FlowPath fp(MeridionalConfig{});
    const auto S = fp.surface(0.0);
    CHECK(S.curvature(S.mprime_inlet() * 0.5) == Approx(0.0).margin(1e-12));
    CHECK(S.curvature(0.5 * (S.mprime_te() + S.mprime_outlet())) == Approx(0.0).margin(1e-12));
    CHECK(std::abs(S.curvature(0.2 * S.mprime_te())) > 1e-3);
    const auto poly = S.polyline(10);
    CHECK(poly.size() == 31u);
    CHECK_THROWS_AS(S.polyline(0), std::invalid_argument);
}
