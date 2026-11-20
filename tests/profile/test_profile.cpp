// Tests of profile/Profile2D: every construction step against its defining equation.
#include "math/Numerics.hpp"
#include "profile/Profile2D.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Profile;
using Catch::Approx;
using Math::Vec2;

TEST_CASE("edge characteristic points follow the ellipse construction", "[profile]")
{
    EdgeParams e;
    e.radius = 0.02;
    e.stretch_ab = 1.5;
    e.asym = 0.004;
    e.beta_deg = 20.0;
    const auto g = build_edge_geometry(Vec2(0, 0), e, +1.0);
    const double b = Math::deg2rad(20.0);
    const Vec2 V(std::cos(b), std::sin(b)), N(-std::sin(b), std::cos(b));
    REQUIRE((g.v_axis - V).norm() < 1e-15);
    REQUIRE((g.p_asym - 0.004 * N).norm() < 1e-15);
    REQUIRE((g.p_mid - (0.004 * N + 0.03 * V)).norm() < 1e-15);
    REQUIRE((g.p_ss - g.p_mid).norm() == Approx(0.02));
    REQUIRE((g.p_ps - g.p_mid).norm() == Approx(0.02));
    REQUIRE(std::abs((g.p_ss - g.p_ps).dot(V)) < 1e-15);  // segment SS-PS normal to V
    REQUIRE(g.p_ss.y() > g.p_ps.y());
    // TE: axis points to the left
    const auto t = build_edge_geometry(Vec2(1, 0), e, -1.0);
    REQUIRE(t.v_axis.x() < 0.0);
}

TEST_CASE("thickness law end values, peak and C1 option", "[profile]")
{
    ThicknessParams t{0.1, 0.4, 0.03, 0.005, 0.7, 1.4, false};
    const double xm = 0.37;
    REQUIRE(thickness_law(0.0, t, xm) == Approx(0.03));
    REQUIRE(thickness_law(xm, t, xm) == Approx(0.1));
    REQUIRE(thickness_law(1.0, t, xm) == Approx(0.005));
    // blade2d's enforce_c1 equalises the slope MAGNITUDES at the peak (the law is a
    // maximum there, so the slopes have opposite signs): a symmetric peak, not C1.
    t.shape_le = 1.3;
    t.enforce_c1 = true;
    const double h = 1e-7;
    const double dl = (thickness_law(xm, t, xm) - thickness_law(xm - h, t, xm)) / h;
    const double dr = (thickness_law(xm + h, t, xm) - thickness_law(xm, t, xm)) / h;
    REQUIRE(dl > 0.0);
    REQUIRE(dr < 0.0);
    REQUIRE(dl == Approx(-dr).epsilon(1e-4));
}

TEST_CASE("NACA 65-(12)10 default: sides, coupling and connectors", "[profile]")
{
    const Profile2D p(ProfileParams::naca65_1210());
    const auto& prm = p.params();
    // suction side ends and tangent directions (handles along V_LE_SS / V_TE_SS)
    REQUIRE((p.suction_side().eval(0.0) - p.p_le_ss()).norm() < 1e-15);
    REQUIRE((p.suction_side().eval(1.0) - p.p_te_ss()).norm() < 1e-15);
    const Vec2 t0 = Math::unit<2>(p.suction_side().derivs(0.0, 1)[1]);
    REQUIRE((t0 - unit_from_beta(prm.le.beta_deg + prm.le.d_beta_ss_deg, +1)).norm() < 1e-12);
    const Vec2 t1 = Math::unit<2>(p.suction_side().derivs(1.0, 1)[1]);
    REQUIRE((-t1 - unit_from_beta(prm.te.beta_deg + prm.te.d_beta_ss_deg, -1)).norm() < 1e-12);
    // pressure side ends
    REQUIRE((p.pressure_side().eval(0.0) - p.p_le_ps()).norm() < 1e-15);
    REQUIRE((p.pressure_side().eval(1.0) - p.p_te_ps()).norm() < 1e-15);
    // free PS control points sit at the law thickness along the SS normal
    for (std::size_t i = 0; i < p.xi_ctrl().size(); ++i) {
        const double u = p.ss_u_at_mprime(prm.ps.mprime_ctrl[i]);
        const Vec2 q = p.suction_side().eval(u);
        REQUIRE(q.x() == Approx(prm.ps.mprime_ctrl[i]).margin(1e-12));
        const double d = (p.ps_free_points()[i] - q).norm();
        REQUIRE(d == Approx(thickness_law(p.xi_ctrl()[i], prm.thickness, p.xi_max())).epsilon(1e-12));
        const double xi = Math::arc_length(p.suction_side(), 0.0, u) / p.ss_length();
        REQUIRE(xi == Approx(p.xi_ctrl()[i]).margin(1e-12));
    }
    // connectors through the edge points
    REQUIRE((p.le_connector().eval(0.0) - p.p_le_ps()).norm() < 1e-13);
    REQUIRE(p.le_connector().eval(p.le_u_mid()).norm() < 1e-13);
    REQUIRE((p.le_connector().eval(1.0) - p.p_le_ss()).norm() < 1e-13);
    REQUIRE((p.te_connector().eval(p.te_u_mid()) - Vec2(1, 0)).norm() < 1e-13);
    REQUIRE((p.point(p.t_le_point())).norm() < 1e-13);
    REQUIRE((p.point(p.t_te_point()) - Vec2(1, 0)).norm() < 1e-13);
}

TEST_CASE("C2Parametric mode: the closed loop is parametrically C2 at all four junctions", "[profile]")
{
    for (auto prm : {ProfileParams::naca65_1210(), ProfileParams::cda_start()}) {
        prm.connector = ConnectorMode::C2Parametric;
        prm.le.tension = 0.28;
        prm.te.tension = 0.2;
        const Profile2D p(prm);
        for (const auto& j : p.junctions()) {
            INFO(j.name);
            REQUIRE(j.gap < 1e-13);
            REQUIRE(j.d1_jump < 1e-10);
            REQUIRE(j.d2_jump < 1e-7);
            REQUIRE(j.tangent_angle < 1e-9);
        }
        // and also periodic: the end of the PS meets the start of the LE connector
        const auto a = p.derivs(p.period(), true), b = p.derivs(0.0, false);
        REQUIRE((a[0] - b[0]).norm() < 1e-13);
        REQUIRE((a[1] - b[1]).norm() < 1e-10 * a[1].norm());
        REQUIRE((a[2] - b[2]).norm() < 1e-7 * a[1].squaredNorm());
    }
}

TEST_CASE("Fair (default) and G2Legacy modes: position, tangent and curvature continuous (G2)", "[profile]")
{
    for (const auto& prm : {ProfileParams::naca65_1210(), ProfileParams::naca65_1210_blade2d(),
                            ProfileParams::cda_start()}) {
        const Profile2D p(prm);
        for (const auto& j : p.junctions()) {
            INFO(to_string(prm.connector) << " " << j.name);
            REQUIRE(j.gap < 1e-13);
            REQUIRE(j.tangent_angle < 1e-9);
            REQUIRE(j.curvature_jump < 1e-8 * std::max(1.0, std::abs(p.curvature(j.t))));
        }
    }
}

TEST_CASE("the Fair connector is the minimum-jerk curve among those meeting its conditions", "[profile]")
{
    const Profile2D p(ProfileParams::naca65_1210());
    const auto& c = p.le_connector();
    REQUIRE(c.degree() == 5);
    REQUIRE(c.n_ctrl() == 10);
    auto energy = [](const Math::BSpline2& b) {
        double e = 0.0;
        for (int i = 0; i < 4000; ++i) e += b.derivs((i + 0.5) / 4000.0, 3)[3].squaredNorm() / 4000.0;
        return e;
    };
    const double e0 = energy(c);
    // admissible perturbations: P0..P2 and P7..P9 carry the end conditions; moving
    // two of P3..P6 with weights that cancel at u_mid keeps C(u_mid) = P_LE
    const int span = Math::find_span(9, 5, p.le_u_mid(), c.knots());
    const auto N = Math::ders_basis_funs(span, p.le_u_mid(), 5, 0, c.knots());
    auto weight = [&](int i) {
        const int j = i - (span - 5);
        return (j >= 0 && j <= 5) ? N[0][static_cast<std::size_t>(j)] : 0.0;
    };
    for (int i = 3; i <= 5; ++i) {
        auto P = c.ctrl();
        const Vec2 d(1e-4, -2e-4);
        P[static_cast<std::size_t>(i)] += weight(i + 1) * d;
        P[static_cast<std::size_t>(i + 1)] -= weight(i) * d;
        const Math::BSpline2 q(5, c.knots(), P);
        REQUIRE((q.eval(p.le_u_mid()) - c.eval(p.le_u_mid())).norm() < 1e-14);
        REQUIRE(energy(q) >= e0 * (1.0 - 1e-9));
    }
}

TEST_CASE("side sampling runs from the LE point to the TE point", "[profile]")
{
    const Profile2D p(ProfileParams::naca65_1210());
    for (auto side : {Side::Suction, Side::Pressure}) {
        const auto c = p.sample_side(side, 30, 100);
        REQUIRE(c.xy.size() == 30u + 99u + 29u);
        REQUIRE(c.xy.front().norm() < 1e-13);
        REQUIRE((c.xy.back() - Vec2(1, 0)).norm() < 1e-13);
        // suction side above the chord, pressure side mostly below the SS
        const Vec2 mid = c.xy[c.xy.size() / 2];
        if (side == Side::Suction) REQUIRE(mid.y() > 0.0);
        for (std::size_t i = 1; i < c.s.size(); ++i) REQUIRE(c.s[i] >= c.s[i - 1]);
    }
    const auto loop = p.sample_loop(20, 60);
    REQUIRE(loop.xy.size() == 4u * 19u + 2u * 59u + 1u);
    REQUIRE((loop.xy.front() - loop.xy.back()).norm() < 1e-13);
}

TEST_CASE("invalid profile parameters throw", "[profile]")
{
    auto bad = ProfileParams::naca65_1210();
    bad.le.radius = 0.0;
    REQUIRE_THROWS_AS(Profile2D(bad), std::invalid_argument);
    bad = ProfileParams::naca65_1210();
    bad.ps.mprime_ctrl = {0.5, 0.3};
    REQUIRE_THROWS_AS(Profile2D(bad), std::invalid_argument);
    bad = ProfileParams::naca65_1210();
    bad.thickness.t_max = std::nan("");
    REQUIRE_THROWS_AS(Profile2D(bad), std::invalid_argument);
    bad = ProfileParams::naca65_1210();
    bad.ss.handle_le_frac = -0.1;
    REQUIRE_THROWS_AS(Profile2D(bad), std::invalid_argument);
    REQUIRE_THROWS_AS(set_param(bad, "no.such", 1.0), std::invalid_argument);
}

TEST_CASE("generic parameter access round-trips", "[profile]")
{
    auto p = ProfileParams::naca65_1210();
    const auto v = to_vector(p);
    REQUIRE(v.size() == param_table().size());
    set_param(p, "le.radius", 0.02);
    REQUIRE(p.le.radius == 0.02);
    REQUIRE(get_param(p, "thickness.t_max") == Approx(0.11027));
    auto q = p;
    from_vector(q, v);
    REQUIRE(q.le.radius == Approx(0.012073));
    const auto m = blend(ProfileParams::naca65_1210(), ProfileParams::cda_start(), 0.5);
    REQUIRE(m.le.radius == Approx(0.5 * (0.012073 + 0.010)));
}
