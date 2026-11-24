// The empirical cascade models of Aungier (2003) Ch. 6 (Correlations.hpp): each correlation against
// hand-computed values of the equations as read (independent arithmetic written out here), the book's
// own figure checks (the "Equation (6-xx)" circles of Figs. 6-2 .. 6-18, read from the page images to
// about +-0.02 of the plotted unit), the limits the book states, the monotonic trends, the two
// corrected misprints, and the refusals.  The Dunavant cascade is in test_correlations_dunavant.cpp.
#include "turbolenza/Correlations.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

using namespace BladeGen::Turbolenza;
using Catch::Approx;

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;
double tand(double a) { return std::tan(a * kDeg); }
double cosd(double a) { return std::cos(a * kDeg); }
double sind(double a) { return std::sin(a * kDeg); }

CascadeGeometry naca_cascade()
{
    CascadeGeometry g;
    g.chord = 0.05;
    g.solidity = 1.2;
    g.stagger = 30.0 * kDeg;
    g.camber = 25.0 * kDeg;
    g.thickness_ratio = 0.08;
    g.throat_over_pitch = 0.7;
    g.pitch_over_height = 0.4;
    return g;
}

/// The QN1 state of a stream surface: W_theta1 = W_m1 tan(beta'_1).
InletState inlet(double beta1_rad, double mach, double wm = 100.0, double r = 0.3, double omega = 0.0)
{
    InletState s;
    s.radius = r;
    s.wm = wm;
    s.wtheta = wm * std::tan(beta1_rad);
    s.omega = omega;
    s.mach = mach;
    return s;
}

DischargeEstimate discharge(const InletState& in, double wm_ratio = 1.0, double radius_ratio = 1.0)
{
    return {in.radius * radius_ratio, in.wm * wm_ratio};
}

}  // namespace

TEST_CASE("Ch. 4 blade angles and Eq. 6-7", "[turbolenza][correlations]")
{
    // circular arc: chi1 = chi2 = theta/2 (4-35), kappa1 = gamma + chi1, kappa2 = gamma - chi2 (4-32, 4-33)
    const BladeAngles c = blade_angles(Camberline::CircularArc, 30.0, 40.0, 0.5);
    CHECK(c.chi1_deg == 20.0);
    CHECK(c.chi2_deg == 20.0);
    CHECK(c.kappa1_deg == 50.0);
    CHECK(c.kappa2_deg == 10.0);
    // the parabolic arc: chi1 + chi2 = theta (Eq. 4-1), the circular values at a/c = 0.5 (p. 76),
    // front-loaded (a/c < 0.5) has chi1 > chi2
    for (const double ac : {0.3, 0.4, 0.5, 0.6, 0.7}) {
        const BladeAngles p = blade_angles(Camberline::ParabolicArc, 10.0, 35.0, ac);
        CHECK(p.chi1_deg + p.chi2_deg == Approx(35.0).margin(1e-12));
        CHECK(p.kappa1_deg - p.kappa2_deg == Approx(35.0).margin(1e-12));
        if (ac < 0.5) CHECK(p.chi1_deg > p.chi2_deg);
        if (ac > 0.5) CHECK(p.chi1_deg < p.chi2_deg);
    }
    const BladeAngles h = blade_angles(Camberline::ParabolicArc, 10.0, 35.0, 0.5);
    CHECK(h.chi1_deg == Approx(17.5).margin(1e-12));
    // b/c of Eq. 4-19 by hand at a/c = 0.4, theta = 35 deg: then tan chi1 = 4b/(4a - c) (4-17)
    {
        const double t4 = 4.0 * tand(35.0);
        const double bc = (std::sqrt(1.0 + t4 * t4 * (0.4 - 0.16 - 0.1875)) - 1.0) / t4;
        CHECK(blade_angles(Camberline::NacaA4K6, 0.0, 35.0, 0.4).chi1_deg ==
              Approx(std::atan(4.0 * bc / 0.6) / kDeg).margin(1e-12));
    }
    CHECK(blade_angles(Camberline::ParabolicArc, 0.0, 0.0, 0.4).chi1_deg == 0.0);
    CHECK_THROWS_AS(blade_angles(Camberline::ParabolicArc, 0.0, 20.0, 0.25), CorrelationError);  // 4-20
    CHECK_THROWS_AS(blade_angles(Camberline::ParabolicArc, 0.0, 20.0, 0.75), CorrelationError);
    CHECK_THROWS_AS(blade_angles(Camberline::Naca65, 0.0, 20.0, 0.4), CorrelationError);  // a/c = 0.5 only
    CHECK_THROWS_AS(blade_angles(Camberline::CircularArc, 0.0, -1.0, 0.5), CorrelationError);
    CHECK_THROWS_AS(blade_angles(Camberline::CircularArc, 90.0, 10.0, 0.5), CorrelationError);

    // Eq. 6-7: tan(theta/4) = 0.05515 C_l0/(a/c); C_l0 = 1.2, a/c = 0.5: theta = 4 atan(0.13236)
    CHECK(naca_camber_from_cl0_deg(1.2, 0.5) == Approx(30.15937).margin(1e-5));
    CHECK(naca_camber_from_cl0_deg(0.0, 0.5) == 0.0);
}

TEST_CASE("Sec. 6.2 design incidence: Eqs. 6-8 .. 6-18", "[turbolenza][correlations]")
{
    CHECK(profile_shape_factor(ProfileFamily::Naca) == 1.0);
    CHECK(profile_shape_factor(ProfileFamily::C4) == 1.1);
    CHECK(profile_shape_factor(ProfileFamily::DoubleCircularArc) == 0.7);

    // K_t,i (6-10, 6-11): 1 at t_b/c = 0.1, 0 at 0; Fig. 6-3: 0.68 at 0.05, 1.19 at 0.15
    CHECK(incidence_thickness_factor(0.1) == Approx(1.0).margin(1e-15));
    CHECK(incidence_thickness_factor(0.0) == 0.0);
    CHECK(incidence_thickness_factor(0.05) == Approx(0.68).margin(0.02));
    CHECK(incidence_thickness_factor(0.15) == Approx(1.19).margin(0.02));
    CHECK(incidence_thickness_factor(0.06) > incidence_thickness_factor(0.04));
    CHECK_THROWS_AS(incidence_thickness_factor(-0.01), CorrelationError);

    // alpha* (6-8, 6-9): hand value for 65-(12)10, sigma 1: 3.6 + 0.3532 * 30.15937 * 0.5^0.25
    {
        const double th = 30.15937;
        CHECK(design_attack_angle_deg(1.0, 1.0, th, 0.5, 1.0) == Approx(3.6 + 0.3532 * th * std::pow(0.5, 0.25)).margin(1e-12));
        CHECK(design_attack_angle_deg(1.0, 1.0, th, 0.5, 1.0) == Approx(12.557).margin(1e-3));
        // the Fig. 6-2 circles: 65-(24)10 at sigma 0.4 ~ 13.0, 65-(27)10 at sigma 0.6 ~ 17.9 (read)
        CHECK(design_attack_angle_deg(1.0, 1.0, naca_camber_from_cl0_deg(2.4, 0.5), 0.5, 0.4) == Approx(13.0).margin(0.3));
        CHECK(design_attack_angle_deg(1.0, 1.0, naca_camber_from_cl0_deg(2.7, 0.5), 0.5, 0.6) == Approx(17.9).margin(0.3));
        // zero camber: alpha* = 3.6 K_sh K_t,i sigma^0.65
        CHECK(design_attack_angle_deg(0.7, 0.9, 0.0, 0.5, 1.5) == Approx(3.6 * 0.63 * std::pow(1.5, 0.65)).margin(1e-12));
        // rises with solidity (Fig. 6-2)
        CHECK(design_attack_angle_deg(1.0, 1.0, th, 0.5, 1.4) > design_attack_angle_deg(1.0, 1.0, th, 0.5, 1.0));
    }

    // Lieblein's zero-camber incidence (6-13, 6-14): the Fig. 6-4 circles at beta1 = 70 deg
    CHECK(lieblein_zero_camber_incidence_deg(70.0, 0.4) == Approx(2.08).margin(0.05));
    CHECK(lieblein_zero_camber_incidence_deg(70.0, 1.2) == Approx(6.26).margin(0.05));
    CHECK(lieblein_zero_camber_incidence_deg(70.0, 2.0) == Approx(10.20).margin(0.05));
    {  // hand: sigma 1, beta1 40: p = 0.92025, 40^p / (5 + 46 e^-2.3) - 0.1 e^-7.5
        const double p = 0.914 + 1.0 / 160.0;
        CHECK(lieblein_zero_camber_incidence_deg(40.0, 1.0) ==
              Approx(std::pow(40.0, p) / (5.0 + 46.0 * std::exp(-2.3)) - 0.1 * std::exp(-7.5)).margin(1e-12));
    }
    CHECK(std::abs(lieblein_zero_camber_incidence_deg(0.0, 1.0)) < 1e-8);  // passes through 0 (Fig. 6-4)
    CHECK_THROWS_AS(lieblein_zero_camber_incidence_deg(-1.0, 1.0), CorrelationError);
    // slope factor n (6-15): Fig. 6-5 at sigma 0.4: -0.46 at 70 deg, -0.05 at 0 deg
    CHECK(lieblein_incidence_slope(70.0, 0.4) == Approx(-0.462).margin(0.01));
    CHECK(lieblein_incidence_slope(0.0, 0.4) == Approx(-0.05).margin(1e-12));

    // the Lieblein design incidence satisfies its defining equation (6-12) at beta1 = kappa1_bar + i*_bar,
    // and for a circular arc with gamma >= 0 (kappa1_bar = kappa1, a/c = 0.5) i* is that i*_bar
    {
        const double kti = incidence_thickness_factor(0.08);
        const BladeAngles ba = blade_angles(Camberline::CircularArc, 30.0, 25.0, 0.5);
        const LiebleinIncidence li = lieblein_design_incidence(1.1, kti, 30.0, 25.0, 0.5, ba.kappa1_deg, 1.2);
        const double b1 = li.pseudo_kappa1_deg + li.pseudo_incidence_deg;
        CHECK(li.pseudo_incidence_deg ==
              Approx(1.1 * kti * lieblein_zero_camber_incidence_deg(b1, 1.2) + lieblein_incidence_slope(b1, 1.2) * 25.0)
                  .margin(1e-11));
        CHECK(li.pseudo_kappa1_deg == Approx(ba.kappa1_deg).margin(1e-12));
        CHECK(li.incidence_deg == Approx(li.pseudo_incidence_deg).margin(1e-12));
        CHECK(li.iterations > 1);
        // a/c != 0.5: Eq. 6-17 shifts alpha* by (a/c - 0.5) theta
        const LiebleinIncidence l4 = lieblein_design_incidence(1.1, kti, 30.0, 25.0, 0.4, ba.kappa1_deg, 1.2);
        CHECK(l4.attack_deg - li.attack_deg == Approx(-0.1 * 25.0).margin(1e-12));
    }
}

TEST_CASE("Sec. 6.3 design deviation: Eqs. 6-19 .. 6-26", "[turbolenza][correlations]")
{
    // (delta*_0)_10 (6-20): Fig. 6-6 circles at 70 deg
    CHECK(lieblein_zero_camber_deviation_deg(70.0, 0.4) == Approx(1.06).margin(0.03));
    CHECK(lieblein_zero_camber_deviation_deg(70.0, 1.0) == Approx(2.57).margin(0.03));
    CHECK(lieblein_zero_camber_deviation_deg(70.0, 2.0) == Approx(4.73).margin(0.05));
    CHECK(lieblein_zero_camber_deviation_deg(0.0, 1.0) == 0.0);
    // m_1.0 (Fig. 6-7) at 70 deg and sigma 1: 0.310 NACA 65, 0.3445 circular arc; 0.17 / 0.249 at 0
    CHECK(lieblein_deviation_slope(Camberline::Naca65, 70.0, 1.0) == Approx(0.17 - 0.0333 * 0.7 + 0.333 * 0.49).margin(1e-14));
    CHECK(lieblein_deviation_slope(Camberline::Naca65, 70.0, 1.0) == Approx(0.310).margin(0.002));
    CHECK(lieblein_deviation_slope(Camberline::CircularArc, 70.0, 1.0) == Approx(0.3445).margin(0.002));
    CHECK(lieblein_deviation_slope(Camberline::Naca65, 0.0, 1.0) == Approx(0.17).margin(1e-15));
    // the exponent b (6-24, Fig. 6-8: 0.55 at 70 deg): m = m_1.0 / sigma^b
    CHECK(lieblein_deviation_slope(Camberline::Naca65, 70.0, 2.0) ==
          Approx(0.310 / std::pow(2.0, 0.9625 - 0.17 * 0.7 - 0.85 * 0.343)).margin(1e-3));
    CHECK_THROWS_AS(lieblein_deviation_slope(Camberline::ParabolicArc, 30.0, 1.0), CorrelationError);
    // K_t,delta (6-25): 1 at 0.1, 1.29 at 0.12 (Fig. 6-9)
    CHECK(deviation_thickness_factor(0.1) == Approx(1.0).margin(1e-15));
    CHECK(deviation_thickness_factor(0.12) == Approx(1.29).margin(1e-12));
    // 6-19: zero camber leaves the thickness term only
    CHECK(lieblein_design_deviation_deg(Camberline::Naca65, 0.7, 0.8, 40.0, 0.0, 1.0) ==
          Approx(0.56 * lieblein_zero_camber_deviation_deg(40.0, 1.0)).margin(1e-14));
    // Howell (6-26): hand value, a/c 0.5, kappa2 10, theta 30, sigma 1.44, NACA 10 % (no thickness term)
    {
        const double rs = 1.2;
        const double hand = (0.92 * 0.25 + 0.002 * 10.0) / (1.0 - 0.002 * 30.0 / rs) * 30.0 / rs;
        CHECK(howell_design_deviation_deg(1.0, 1.0, 0.5, 10.0, 30.0, -5.0, 1.44) == Approx(hand).margin(1e-13));
        // a C.4 profile adds (K_sh K_t,delta - 1)(delta*_0)_10
        CHECK(howell_design_deviation_deg(1.1, 1.0, 0.5, 10.0, 30.0, 40.0, 1.44) ==
              Approx(hand + 0.1 * lieblein_zero_camber_deviation_deg(40.0, 1.44)).margin(1e-13));
        // its defining form: Howell's delta = [0.23 (2a/c)^2 + beta2/500] theta/sqrt(sigma) with beta2 = kappa2 + delta
        const double d = howell_design_deviation_deg(1.0, 1.0, 0.4, 12.0, 30.0, 0.0, 1.44);
        CHECK(d == Approx((0.23 * 0.64 + (12.0 + d) / 500.0) * 30.0 / rs).margin(1e-12));
    }
}

TEST_CASE("Sec. 6.4 diffusion factors and design loss: Eqs. 6-27 .. 6-46, 6-99 .. 6-102", "[turbolenza][correlations]")
{
    // D (6-31) from velocity components: W1 = W_m (1, tan b1), W2 = (W_m2/W_m1) W_m (1, tan b2)
    {
        const double wm = 100.0, r = 0.9, b1 = 50.0, b2 = 20.0, sg = 1.3;
        const double w1 = wm / cosd(b1), w2 = r * wm / cosd(b2);
        const double hand = 1.0 - w2 / w1 + (wm * tand(b1) - r * wm * tand(b2)) / (2.0 * sg * w1);
        CHECK(diffusion_factor(b1, b2, sg, r) == Approx(hand).margin(1e-14));
        CHECK(velocity_ratio_w2_w1(b1, b2, r) == Approx(w2 / w1).margin(1e-14));
        CHECK(diffusion_factor(b1, b1, sg, 1.0) == Approx(0.0).margin(1e-15));  // no turning, no diffusion
        CHECK(diffusion_factor_from_velocities(1.3, 1.0, 0.8) == Approx(0.5).margin(1e-15));    // 6-29
        CHECK(equivalent_diffusion_from_velocities(1.3, 0.65) == Approx(2.0).margin(1e-15));    // 6-30
    }
    // W_max/W1: Eq. 6-39 is Eq. 6-35 for a planar stationary cascade at i = i*; rotation cancels for r2 = r1
    {
        const double p = lieblein_wmax_ratio_design(45.0, 20.0, 1.1);
        CHECK(p == Approx(1.12 + 0.61 * cosd(45.0) * cosd(45.0) / 1.1 * (1.0 - tand(20.0))).margin(1e-14));
        CHECK(wmax_ratio(45.0, 20.0, 1.1, 1.0, 1.0, 0.0, 0.0, 0.0117) == Approx(p).margin(1e-14));
        CHECK(wmax_ratio(45.0, 20.0, 1.1, 1.0, 1.0, -1.7, 0.0, 0.0117) == Approx(p).margin(1e-14));
        // hand value of Klapproth's term with r2/r1 = 1.1, W_m2/W_m1 = 0.9, U1/W_m1 = 1.5:
        // (r1 C_t1 - r2 C_t2)/(r1 C_m1) with C_t = W_m tan b + omega r
        const double rr = 1.1, vr = 0.9, u = 1.5;
        const double ct1 = tand(45.0) + u, ct2 = vr * tand(20.0) + u * rr;
        const double hand = 1.12 + 0.61 * cosd(45.0) * cosd(45.0) / 1.1 * (ct1 - rr * ct2) + 0.0117 * std::pow(3.0, 1.43);
        CHECK(wmax_ratio(45.0, 20.0, 1.1, vr, rr, u, 3.0, 0.0117) == Approx(hand).margin(1e-14));
        CHECK_THROWS_AS(wmax_ratio(45.0, 20.0, 1.1, 1.0, 1.0, 0.0, -1.0, 0.0117), CorrelationError);  // i < i*
    }
    CHECK(lieblein_alpha(ProfileFamily::Naca) == 0.0117);
    CHECK(lieblein_alpha(ProfileFamily::C4) == 0.007);
    CHECK_THROWS_AS(lieblein_alpha(ProfileFamily::DoubleCircularArc), CorrelationError);
    // D_eq (6-36) = (W_max/W1) cos b2 W_m1 / (cos b1 W_m2)
    CHECK(equivalent_diffusion(1.4, 50.0, 20.0, 0.8) == Approx(1.4 * cosd(20.0) / (cosd(50.0) * 0.8)).margin(1e-14));

    // the loss correlations: Fig. 6-11 at D* = 0.6: omega* cos b2 / (2 sigma) = 0.02763; Fig. 6-12 at
    // D_eq* = 2: omega* cos b2 (W1/W2)^2/(2 sigma) = 0.0180; 0.004 at D_eq* = 1
    {
        const double sg = 1.2, b2 = 15.0;
        CHECK(design_loss_d_factor(0.6, b2, sg, 0.0035, 1.0) * cosd(b2) / (2.0 * sg) == Approx(0.02763).margin(2e-5));
        CHECK(design_loss_d_factor(0.0, b2, sg, 0.0035, 1.0) * cosd(b2) / (2.0 * sg) == Approx(0.0035).margin(1e-15));
        CHECK_THROWS_AS(design_loss_d_factor(0.61, b2, sg, 0.0035, 1.0), CorrelationError);  // the diffusion limit
        const double w12 = 1.3;
        CHECK(design_loss_equivalent_diffusion(2.0, b2, sg, w12, 0.004, 1.0) * cosd(b2) * w12 * w12 / (2.0 * sg) ==
              Approx(0.004 * 4.5).margin(1e-15));
        CHECK(design_loss_equivalent_diffusion(1.0, b2, sg, w12, 0.004, 1.0) * cosd(b2) * w12 * w12 / (2.0 * sg) ==
              Approx(0.004).margin(1e-16));
        // increasing in D_eq* above 1 (Fig. 6-12)
        CHECK(design_loss_equivalent_diffusion(1.8, b2, sg, w12, 0.004, 1.0) >
              design_loss_equivalent_diffusion(1.6, b2, sg, w12, 0.004, 1.0));
    }
    // the wake form (6-27) and Howell's (6-40 .. 6-44)
    CHECK(wake_loss_coefficient(0.01, 1.2, 20.0, 0.8) == Approx(2.0 * 0.01 * 1.2 / cosd(20.0) * 0.64).margin(1e-15));
    {
        const double bm = std::atan(0.5 * (tand(50.0) + tand(20.0))) / kDeg;
        CHECK(howell_mean_angle_deg(50.0, 20.0) == Approx(bm).margin(1e-12));
        const double cl = 2.0 * cosd(bm) * (tand(50.0) - tand(20.0)) / 1.2;
        CHECK(howell_lift_coefficient(50.0, 20.0, 1.2) == Approx(cl).margin(1e-14));
        CHECK(howell_secondary_drag(cl) == Approx(0.18 * cl * cl).margin(1e-15));
        CHECK(howell_endwall_drag(0.5) == Approx(0.01).margin(1e-16));
        CHECK(howell_loss_from_drag(0.02, 1.2, 50.0, 20.0) ==
              Approx(0.02 * 1.2 * cosd(50.0) * cosd(50.0) / std::pow(cosd(bm), 3)).margin(1e-15));
    }
    // K2 (6-99, 6-100) and K_Re (6-101, 6-102): K_Re = 0 at 2.5e5 from both sides, 1 at 6.25e4
    CHECK(aspect_ratio_k2(0.5, 20.0, 0.0, 0.0073) == Approx(1.0 + 0.5 * cosd(20.0)).margin(1e-15));
    CHECK(aspect_ratio_k2(0.5, 20.0, 0.3, 0.0073) == Approx(1.0 + 0.5 * cosd(20.0) + 0.004 * 0.3 / 0.0073).margin(1e-14));
    CHECK(reynolds_factor(2.5e5) == Approx(0.0).margin(1e-15));
    CHECK(reynolds_factor(2.5e5 * (1.0 - 1e-12)) == Approx(0.0).margin(1e-11));
    CHECK(reynolds_factor(6.25e4) == Approx(1.0).margin(1e-14));
    CHECK(reynolds_factor(2.5e6) == Approx(std::pow(std::log10(2.5e5) / std::log10(2.5e6), 2.58) - 1.0).margin(1e-14));
    CHECK(reynolds_factor(2.5e6) < 0.0);  // less skin friction above the basis Reynolds number
    CHECK_THROWS_AS(reynolds_factor(0.0), CorrelationError);
}

TEST_CASE("Sec. 6.5 stall incidence: Eqs. 6-47 .. 6-50", "[turbolenza][correlations]")
{
    // Fig. 6-14: beta_1c = 30 deg gives -9 for any camber; (60, 60) -4.93; positive (30, 45) 15.77, (70, 35) 3.61
    CHECK(negative_stall_attack_offset_deg(30.0, 0.0) == -9.0);
    CHECK(negative_stall_attack_offset_deg(30.0, 55.0) == Approx(-9.0).margin(1e-14));
    CHECK(negative_stall_attack_offset_deg(60.0, 60.0) == Approx(-4.93).margin(0.01));
    CHECK(positive_stall_attack_offset_deg(30.0, 45.0) == Approx(15.77).margin(0.01));
    CHECK(positive_stall_attack_offset_deg(70.0, 35.0) == Approx(3.61).margin(0.01));
    CHECK(positive_stall_attack_offset_deg(45.0, 0.0) == 10.3);
    CHECK_THROWS_AS(negative_stall_attack_offset_deg(19.9, 20.0), CorrelationError);  // beta_1c >= 20 (p. 135)

    // the solved ranges satisfy beta_1 = gamma + alpha (p. 135) with Eqs. 6-47 and 6-48
    const StallRanges r = stall_ranges(10.0, 30.0, 30.0);
    CHECK(r.beta1s_deg == Approx(30.0 + 10.0 + positive_stall_attack_offset_deg(r.beta1s_deg, 30.0)).margin(1e-11));
    CHECK(r.beta1c_deg == Approx(30.0 + 10.0 + negative_stall_attack_offset_deg(r.beta1c_deg, 30.0)).margin(1e-11));
    CHECK(r.positive_deg == Approx(r.beta1s_deg - 40.0).margin(1e-11));
    CHECK(r.negative_deg == Approx(40.0 - r.beta1c_deg).margin(1e-11));
    CHECK_FALSE(r.beta1c_limited);
    // a low-stagger cascade: beta_1c < 20 deg, Eq. 6-47 at 20 deg (Aungier's limit)
    const StallRanges lo = stall_ranges(5.0, 5.0, 20.0);
    CHECK(lo.beta1c_limited);
    CHECK(lo.negative_deg == Approx(-negative_stall_attack_offset_deg(20.0, 20.0)).margin(1e-12));
}

TEST_CASE("Sec. 6.6 Mach number effects: Eqs. 6-51 .. 6-60", "[turbolenza][correlations]")
{
    // 6-51, 6-52: ranges shrink by 1 + 0.5 M^3; K_sh limited to 1 (a C.4, 1.1, acts as 1); DCA (0.7) shrinks i_s less
    CHECK(mach_negative_stall_incidence_deg(-2.0, 8.0, 0.8) == Approx(-2.0 - 8.0 / (1.0 + 0.5 * 0.512)).margin(1e-14));
    CHECK(mach_positive_stall_incidence_deg(-2.0, 8.0, 0.8, 1.1) == Approx(-2.0 + 8.0 / (1.0 + 0.5 * 0.512)).margin(1e-14));
    CHECK(mach_positive_stall_incidence_deg(-2.0, 8.0, 0.8, 0.7) ==
          Approx(-2.0 + 8.0 / (1.0 + 0.5 * std::pow(0.56, 3))).margin(1e-14));
    CHECK(mach_positive_stall_incidence_deg(0.0, 8.0, 0.9, 1.0) < mach_positive_stall_incidence_deg(0.0, 8.0, 0.5, 1.0));
    // 6-55: i_m = i* when both ranges shrink alike (p. 138), i_m > i* when i_s shrinks less
    {
        const double ic = mach_negative_stall_incidence_deg(-2.0, 8.0, 0.8), is = mach_positive_stall_incidence_deg(-2.0, 6.0, 0.8, 1.0);
        CHECK(minimum_loss_incidence_deg(ic, is, 8.0, 6.0) == Approx(-2.0).margin(1e-13));
        const double is7 = mach_positive_stall_incidence_deg(-2.0, 6.0, 0.8, 0.7);
        CHECK(minimum_loss_incidence_deg(ic, is7, 8.0, 6.0) > -2.0);
    }
    // the perfect-gas mass flux: 1 at M = 1, the A*/A of the isentropic table at M = 0.65 (A/A* = 1.1356, gamma 1.4)
    CHECK(mass_flux_ratio_to_sonic(1.0, 1.4) == Approx(1.0).margin(1e-15));
    CHECK(1.0 / mass_flux_ratio_to_sonic(0.65, 1.4) == Approx(1.1356).margin(1e-4));
    CHECK(sonic_velocity_ratio(1.0, 1.4) == Approx(1.0).margin(1e-15));
    CHECK(sonic_velocity_ratio(0.5, 1.4) == Approx(std::sqrt(2.1 / 0.6)).margin(1e-14));
    // 6-53: cos(beta_choke) = (o/s) / (rho1 W1 / rho* W*); no choke when that is >= 1
    CHECK(choke_inlet_angle_deg(0.7, 0.65, 1.4) == Approx(std::acos(0.7 * 1.1356) / kDeg).margin(5e-3));
    CHECK(std::isnan(choke_inlet_angle_deg(0.7, 0.3, 1.4)));
    CHECK(choke_inlet_angle_deg(0.6, 1.0, 1.4) == Approx(std::acos(0.6) / kDeg).margin(1e-12));
    // 6-59 and 6-58 / 6-60
    CHECK(critical_mach_number(0.7, 1.2, 1.4) == Approx(0.6).margin(1e-15));
    CHECK(minimum_loss(0.02, 1.0, -1.0, 8.0, 1.0, 0.5, 0.7, 1.3) == Approx(0.02 * (1.0 + 4.0 / 64.0)).margin(1e-16));
    CHECK(minimum_loss(0.02, -1.0, -1.0, 8.0, 0.7, 0.8, 0.7, 1.1) == Approx(0.02 + 0.7 * std::pow((0.8 / 0.7 - 1.0) * 1.1, 2)).margin(1e-15));
    // M' limited to 1 (p. 138): M'_1 = 1.3 and 1.0 give the same Mach term
    CHECK(minimum_loss(0.02, -1.0, -1.0, 8.0, 1.0, 1.3, 0.7, 0.9) == Approx(minimum_loss(0.02, -1.0, -1.0, 8.0, 1.0, 1.0, 0.7, 0.9)).margin(1e-16));
}

TEST_CASE("Sec. 6.7 shock loss: Eqs. 6-61 .. 6-75 (perfect gas)", "[turbolenza][correlations]")
{
    // the corrected Eq. 6-61 (c/(2 R_u) = sin(theta_u/2), Eq. 4-7): R_u >= c/2 for every camber, the arc
    // spans the chord; the printed 2R_u/c = sin(theta_u/2) would give R_u = c sin(theta_u/2)/2 < c/2
    for (const double th : {2.0, 10.0, 30.0, 60.0}) {
        const ShockGeometry g = shock_expansion_geometry(th, 0.05, 40.0, 1.5);
        const double ru_c = g.radius_over_pitch / 1.5;
        CHECK(ru_c >= 0.5);
        CHECK(ru_c * 2.0 * std::sin(g.suction_camber_deg * kDeg / 2.0) == Approx(1.0).margin(1e-14));
        CHECK(0.5 * std::sin(g.suction_camber_deg * kDeg / 2.0) < 0.5);  // the printed form's radius
        // 6-65: tan(theta_u/4) = tan(theta/4) + t_b/c; 6-68; 6-69
        CHECK(std::tan(g.suction_camber_deg * kDeg / 4.0) == Approx(std::tan(th * kDeg / 4.0) + 0.05).margin(1e-14));
        CHECK(g.psi_deg == Approx(90.0 - g.suction_camber_deg / 2.0 - 40.0).margin(1e-12));
        CHECK(std::tan(g.expansion_deg * kDeg) ==
              Approx(cosd(g.psi_deg) / (sind(g.psi_deg) + g.radius_over_pitch)).margin(1e-13));
    }
    // a flat, thin plate: no expansion
    CHECK(shock_expansion_geometry(0.0, 0.0, 50.0, 1.0).expansion_deg == 0.0);
    // the expansion grows with camber and thickness
    CHECK(shock_expansion_geometry(10.0, 0.06, 50.0, 1.0).expansion_deg > shock_expansion_geometry(10.0, 0.03, 50.0, 1.0).expansion_deg);

    // Prandtl-Meyer: nu(1) = 0, nu(2) = 26.3798 deg (gamma 1.4, the standard table); Eq. 6-70 integrated
    // numerically (dW/W = dM / (M (1 + (g-1)/2 M^2)) on the isentrope) gives the closed form
    CHECK(prandtl_meyer_angle(1.0, 1.4) == 0.0);
    CHECK(prandtl_meyer_angle(2.0, 1.4) / kDeg == Approx(26.3798).margin(1e-4));
    {
        const int n = 20000;
        double phi = 0.0;
        for (int k = 0; k < n; ++k) {  // midpoint rule in M from 1.2 to 1.6
            const double m = 1.2 + 0.4 * (k + 0.5) / n;
            phi += std::sqrt(m * m - 1.0) / (m * (1.0 + 0.2 * m * m)) * (0.4 / n);
        }
        CHECK(phi == Approx(prandtl_meyer_angle(1.6, 1.4) - prandtl_meyer_angle(1.2, 1.4)).margin(1e-9));
    }
    CHECK(prandtl_meyer_mach(prandtl_meyer_angle(1.7, 1.4), 1.4) == Approx(1.7).margin(1e-12));
    CHECK_THROWS_AS(prandtl_meyer_mach(2.3, 1.4), CorrelationError);  // beyond nu_max = 2.27 rad
    // normal shock P_t2/P_t1: 1 at M = 1, 0.72087 at M = 2 (standard table)
    CHECK(normal_shock_total_pressure_ratio(1.0, 1.4) == Approx(1.0).margin(1e-15));
    CHECK(normal_shock_total_pressure_ratio(2.0, 1.4) == Approx(0.72087).margin(1e-5));
    // the shock loss: 0 subsonic, positive and growing above M = 1; omega_s = (1 - P_t2/P_t1)/(1 - P1/P_t1)
    CHECK(shock_loss(0.9, 10.0, 0.05, 50.0, 1.2, 1.4).loss == 0.0);
    const ShockLoss s = shock_loss(1.3, 10.0, 0.05, 50.0, 1.2, 1.4);
    CHECK(s.surface_mach > 1.3);
    CHECK(s.shock_mach == Approx(std::sqrt(1.3 * s.surface_mach)).margin(1e-15));  // 6-71: the geometric mean
    CHECK(s.loss == Approx((1.0 - normal_shock_total_pressure_ratio(s.shock_mach, 1.4)) /
                           (1.0 - std::pow(1.0 + 0.2 * 1.69, -3.5)))
                        .margin(1e-15));
    CHECK(shock_loss(1.4, 10.0, 0.05, 50.0, 1.2, 1.4).loss > s.loss);
    CHECK(s.loss > 0.0);
}

TEST_CASE("Sec. 6.8 off-design: Eqs. 6-76 .. 6-84", "[turbolenza][correlations]")
{
    // 6-76: 1 at sigma 0 (Fig. 6-18), exp(-3.1) at sigma 1, beta1 0; decreasing in solidity
    CHECK(deviation_incidence_slope(40.0, 1e-12) == Approx(1.0).margin(1e-9));
    CHECK(deviation_incidence_slope(0.0, 1.0) == Approx(std::exp(-3.1)).margin(1e-15));
    CHECK(deviation_incidence_slope(60.0, 1.0) > deviation_incidence_slope(40.0, 1.0));
    CHECK(deviation_incidence_slope(60.0, 1.2) < deviation_incidence_slope(60.0, 0.8));
    // 6-77: delta* at i*, the AVDR term, no incidence term for a guide vane
    CHECK(off_design_deviation_deg(6.0, 0.1, 2.0, 2.0, 1.0, false) == 6.0);
    CHECK(off_design_deviation_deg(6.0, 0.1, 5.0, 2.0, 0.9, false) == Approx(6.0 + 0.3 + 1.0).margin(1e-14));
    CHECK(off_design_deviation_deg(6.0, 0.1, 5.0, 2.0, 1.0, true) == 6.0);
    // 6-78 .. 6-82: xi = -1 at i_c and +1 at i_s, where the loss is 2 omega_m (Fig. 6-13); C1 joins at -2, 1
    CHECK(normalized_incidence(-8.0, -1.0, -8.0, 5.0) == -1.0);
    CHECK(normalized_incidence(5.0, -1.0, -8.0, 5.0) == 1.0);
    CHECK(off_design_loss_factor(0.0) == 1.0);
    CHECK(off_design_loss_factor(1.0) == 2.0);
    CHECK(off_design_loss_factor(-1.0) == 2.0);
    CHECK(off_design_loss_factor(-2.0) == 5.0);
    CHECK(off_design_loss_factor(-2.0 - 1e-9) == Approx(5.0).margin(1e-8));
    CHECK(off_design_loss_factor(3.0) == 6.0);
    CHECK(off_design_loss_factor(-3.0) == 9.0);
    CHECK((off_design_loss_factor(1.0 + 1e-6) - off_design_loss_factor(1.0)) / 1e-6 == Approx(2.0).margin(1e-5));
    CHECK_THROWS_AS(normalized_incidence(0.0, 5.0, 5.0, 6.0), CorrelationError);  // no range
    // the wake blockage: 0 up to D_eq = 2, continuous, 1 - (2/D_eq)^0.9 above
    CHECK(wake_blockage(2.0) == 0.0);
    CHECK(wake_blockage(2.0 + 1e-9) == Approx(0.0).margin(1e-9));
    CHECK(wake_blockage(2.5) == Approx(1.0 - std::pow(0.8, 0.9)).margin(1e-15));
}

TEST_CASE("Secs. 6.9, 6.10 tip clearance and shroud seal: Eqs. 6-85 .. 6-96", "[turbolenza][correlations]")
{
    TipClearanceInput t;
    t.clearance = 0.0005;
    t.tip_radius = 0.25;
    t.blade_count = 30;
    t.chord = 0.05;
    t.stagger_rad = 40.0 * kDeg;
    t.row_number = 2;
    t.r1 = 0.249;
    t.r2 = 0.248;
    t.rho1 = 1.2;
    t.rho2 = 1.3;
    t.cm1 = 150.0;
    t.cm2 = 145.0;
    t.ctheta1 = 20.0;
    t.ctheta2 = 180.0;
    t.mean_density = 1.25;
    t.mass_flow = 30.0;
    const TipClearanceLoss a = tip_clearance_loss(t);
    const double tau = std::numbers::pi * 0.0005 * (0.249 * 1.2 * 150.0 + 0.248 * 1.3 * 145.0) * (0.248 * 180.0 - 0.249 * 20.0);
    const double ax = 0.05 * cosd(40.0);
    const double dp = tau / (30.0 * 0.25 * 0.0005 * ax);
    const double uc = 0.816 * std::sqrt(2.0 * dp / 1.25) / std::pow(2.0, 0.2);
    const double mc = 1.25 * uc * 30.0 * 0.0005 * ax;
    CHECK(a.torque == Approx(tau).margin(1e-12));
    CHECK(a.pressure_difference == Approx(dp).margin(1e-9));
    CHECK(a.leakage_velocity == Approx(uc).margin(1e-11));
    CHECK(a.leakage_mass_flow == Approx(mc).margin(1e-14));
    CHECK(a.total_pressure_loss == Approx(dp * mc / 30.0).margin(1e-10));
    // a stator (swirl removed): the same gap flow with tau < 0 gives the same loss (|tau|, ours)
    TipClearanceInput st = t;
    std::swap(st.ctheta1, st.ctheta2);
    st.r1 = t.r2;
    st.r2 = t.r1;
    std::swap(st.rho1, st.rho2);
    std::swap(st.cm1, st.cm2);
    const TipClearanceLoss b = tip_clearance_loss(st);
    CHECK(b.torque == Approx(-a.torque).margin(1e-12));
    CHECK(b.total_pressure_loss == Approx(a.total_pressure_loss).margin(1e-10));
    // the later row leaks less (N_row^-0.2)
    st.row_number = 6;
    CHECK(tip_clearance_loss(st).total_pressure_loss < b.total_pressure_loss);
    TipClearanceInput bad = t;
    bad.blade_count = 0;
    CHECK_THROWS_AS(tip_clearance_loss(bad), CorrelationError);

    // C_r (6-92): ~1 for a tiny clearance, -> 2/3 for a large one (Aungier 2000 Fig. 4-7)
    CHECK(seal_contraction_ratio(0.0) == Approx(1.0 - 1.0 / (3.0 + std::pow(54.3, 3.45))).margin(1e-15));
    CHECK(seal_contraction_ratio(100.0) == Approx(2.0 / 3.0).margin(1e-4));
    // C_t: Aungier (2000) Eq. 4-35; Fig. 4-8 at P_R = 0 reads ~0.94 (N 1), ~0.71 (N 2), ~0.52 (N 4)
    CHECK(seal_throttling_coefficient(1, 0.0) == Approx(0.944).margin(0.005));
    CHECK(seal_throttling_coefficient(2, 0.0) == Approx(0.711).margin(0.005));
    CHECK(seal_throttling_coefficient(4, 0.0) == Approx(0.517).margin(0.005));
    for (int n = 1; n < 30; ++n) CHECK(seal_throttling_coefficient(n + 1, 0.3) < seal_throttling_coefficient(n, 0.3));
    CHECK(seal_throttling_coefficient(5, 0.5) == Approx(2.143 * (std::log(5.0) - 1.464) / (5.0 - 4.322) * std::pow(0.5, 0.1875)).margin(1e-14));
    // the 2003 print, 2.143 ln(N) - 1.464 over N - 4.322: negative for N = 2 .. 4, a pole at 4.322
    {
        const auto printed = [](double n) { return (2.143 * std::log(n) - 1.464) / (n - 4.322); };
        CHECK(printed(3.0) < 0.0);
        CHECK(std::abs(printed(4.3225)) > 1e3);  // the pole at N = 4.322
        // the bracketed form: its numerator vanishes at N = e^1.464 = 4.3230, 0.001 from the pole, so it
        // stays finite there (a near-cancellation; N is an integer, so it never matters)
        CHECK(std::abs(std::exp(1.464) - 4.322) < 2e-3);
        CHECK(std::abs(2.143 * (std::log(4.3225) - 1.464) / (4.3225 - 4.322)) < 1.0);
    }
    // C_c (6-94 .. 6-96): 1 at zero clearance; refused beyond X2 - 1 (p. 149); N > 12 by Aungier 2000 Eqs. 4-36, 4-37
    CHECK(seal_carryover_coefficient(4, 0.0) == 1.0);
    {
        const double x1 = 15.1 - 0.05255 * std::exp(0.507 * 8.0), x2 = 1.058 + 0.0218 * 4.0;
        CHECK(seal_carryover_coefficient(4, 0.05) == Approx(1.0 + x1 * (0.05 - x2 * std::log(1.05)) / (1.0 - x2)).margin(1e-14));
        CHECK_THROWS_AS(seal_carryover_coefficient(4, x2 - 1.0 + 1e-9), CorrelationError);
    }
    CHECK(seal_carryover_coefficient(15, 0.05) ==
          Approx(1.0 + (13.15 + 0.1625 * 15.0) * (0.05 - 1.32 * std::log(1.05)) / (1.0 - 1.32)).margin(1e-14));
    CHECK(seal_carryover_coefficient(4, 0.05) > 1.0);
    CHECK_THROWS_AS(seal_throttling_coefficient(0, 0.5), CorrelationError);
    CHECK_THROWS_AS(seal_throttling_coefficient(3, 1.0), CorrelationError);

    ShroudSealInput s;
    s.seal_radius = 0.15;
    s.clearance = 0.0004;
    s.fins = 3;
    s.fin_thickness = 0.0005;
    s.fin_pitch = 0.004;
    s.pressure_ratio = 0.9;
    s.upstream_density = 1.4;
    s.upstream_temperature = 330.0;
    s.gas_constant = 287.0;
    s.pressure_difference = 8000.0;
    s.mass_flow = 30.0;
    const ShroudSealLoss l = shroud_seal_loss(s);
    const double ml = 2.0 * std::numbers::pi * 0.15 * 0.0004 * seal_throttling_coefficient(3, 0.9) *
                      seal_carryover_coefficient(3, 0.1) * seal_contraction_ratio(0.8) * 1.4 * std::sqrt(287.0 * 330.0);
    CHECK(l.leakage_mass_flow == Approx(ml).margin(1e-14));
    CHECK(l.total_pressure_loss == Approx(8000.0 * ml / 30.0).margin(1e-10));
}

TEST_CASE("evaluate: design point, frame sign, guide vanes, trends, refusals", "[turbolenza][correlations]")
{
    CorrelationOptions opt;  // Aungier's: alpha* (6-8), Eq. 6-46 with K1 = 0.0073, K2 of 6-99
    const Gas gas;
    CascadeGeometry g = naca_cascade();
    InletState in = inlet(0.0, 0.4);

    // at beta1 = beta1*, W_m2 = W_m1, low Mach: i = i*, delta = delta*, xi = 0, omega = omega*
    const CascadePerformance d0 = evaluate(g, in, discharge(in), gas, opt);
    in = inlet(d0.design_beta1_deg * kDeg, 0.4);
    const CascadePerformance d = evaluate(g, in, discharge(in), gas, opt);
    CHECK(d.incidence_deg == Approx(d.design_incidence_deg).margin(1e-12));
    CHECK(d.deviation_deg == Approx(d.design_deviation_deg).margin(1e-12));
    CHECK(d.minimum_loss_incidence_deg == Approx(d.design_incidence_deg).margin(1e-12));  // K_sh = 1, no choke
    CHECK(d.xi == Approx(0.0).margin(1e-12));
    CHECK(d.loss == Approx(d.losses.design).margin(1e-15));
    CHECK(d.d_eq == Approx(d.d_eq_design).margin(1e-12));
    CHECK(d.losses.k1 == 0.0073);
    CHECK(d.losses.k2 == Approx(1.0 + 0.4 * cosd(d.design_beta2_deg)).margin(1e-14));
    CHECK(d.beta2 == Approx(d.beta2_deg * kDeg).margin(1e-15));
    // the design values by hand: kappa1 = 42.5, alpha* (6-8), i* = alpha* + gamma - kappa1 (6-18)
    const double kti = std::pow(0.8, 0.28 / (0.1 + std::pow(0.08, 0.3)));
    const double astar = (3.6 * kti + 0.3532 * 25.0 * std::pow(0.5, 0.25)) * std::pow(1.2, 0.65 - 0.05);
    CHECK(d.design_attack_deg == Approx(astar).margin(1e-12));
    CHECK(d.design_incidence_deg == Approx(astar + 30.0 - 42.5).margin(1e-12));
    CHECK(d.frame_sign == 1);
    CHECK(d.reynolds == 0.0);  // no density / viscosity given

    // the loss rises away from i_m on both sides and the deviation with incidence (Fig. 6-13, Eq. 6-77)
    double prev_dev = -1e9;
    for (const double di : {-6.0, -3.0, 0.0, 3.0, 6.0}) {
        const InletState ii = inlet((d.design_beta1_deg + di) * kDeg, 0.4);
        const CascadePerformance p = evaluate(g, ii, discharge(ii), gas, opt);
        CHECK(p.loss >= d.loss - 1e-15);
        CHECK(p.deviation_deg > prev_dev);
        prev_dev = p.deviation_deg;
        CHECK(p.positive_stall_margin_deg == Approx(p.positive_stall_incidence_deg - p.incidence_deg).margin(1e-12));
        if (di > 0.0) CHECK(p.wmax_ratio > wmax_ratio(p.beta1_deg, p.beta2_deg, g.solidity, 1.0, 1.0, 0.0, 0.0, 0.0));
    }
    // the stall points carry twice the minimum loss (Fig. 6-13)
    {
        const InletState ii = inlet((d.kappa1_deg + d.positive_stall_incidence_deg) * kDeg, 0.4);
        CHECK(evaluate(g, ii, discharge(ii), gas, opt).loss == Approx(2.0 * d.losses.minimum).margin(1e-12));
    }

    // s_row = -1: the mirrored row (all signed inputs negated: stagger, camber, W_theta1, omega) gives the
    // mirrored answer; here a rotor with r2/r1 = 1.05, W_m2/W_m1 = 0.95, U1/W_m1 = 1.2
    CascadeGeometry gr = g;
    gr.kind = BladeRowKind::Rotor;
    gr.row_sign = -1;
    gr.stagger = -g.stagger;
    gr.camber = -g.camber;
    const InletState ip = inlet((d.design_beta1_deg + 2.0) * kDeg, 0.4, 100.0, 0.3, 1.2 * 100.0 / 0.3);
    InletState im = ip;
    im.wtheta = -ip.wtheta;
    im.omega = -ip.omega;
    const CascadePerformance ps = evaluate(g, ip, discharge(ip, 0.95, 1.05), gas, opt);
    const CascadePerformance pr = evaluate(gr, im, discharge(im, 0.95, 1.05), gas, opt);
    CHECK(pr.frame_sign == -1);
    CHECK(pr.beta2 == Approx(-ps.beta2).margin(1e-15));
    CHECK(pr.loss == Approx(ps.loss).margin(1e-15));
    CHECK(pr.d_eq == Approx(ps.d_eq).margin(1e-15));
    // Klapproth's term sees the rotation: U1/W_m1 = 1.2 with r2 = 1.05 r1 lowers W_max/W1 by
    // 0.61 cos^2(b1)/sigma * 1.2 (1 - 1.05^2) relative to omega = 0
    {
        InletState i0 = ip;
        i0.omega = 0.0;
        const CascadePerformance p0 = evaluate(g, i0, discharge(i0, 0.95, 1.05), gas, opt);
        CHECK(ps.wmax_ratio - p0.wmax_ratio ==
              Approx(0.61 * cosd(ps.beta1_deg) * cosd(ps.beta1_deg) / g.solidity * 1.2 * (1.0 - 1.05 * 1.05)).margin(1e-13));
        CHECK(ps.beta2 == p0.beta2);  // the deviation does not depend on rotation
    }
    // the AVDR term of Eq. 6-77: W_m2/W_m1 = 0.95 adds 0.5 deg of deviation
    CHECK(ps.deviation_deg - evaluate(g, ip, discharge(ip, 1.0, 1.05), gas, opt).deviation_deg == Approx(0.5).margin(1e-12));
    // a non-IGV row of negative camber after s_row is refused (contract Sec. 0.2)
    CascadeGeometry gwrong = g;
    gwrong.row_sign = -1;
    CHECK_THROWS_AS(evaluate(gwrong, ip, discharge(ip), gas, opt), CorrelationError);
    gwrong.row_sign = 0;
    CHECK_THROWS_AS(evaluate(gwrong, ip, discharge(ip), gas, opt), CorrelationError);

    // an inlet guide vane of negative camber: Aungier's flip (p. 128), Howell's deviation, no incidence term
    CascadeGeometry gv;
    gv.kind = BladeRowKind::InletGuideVane;
    gv.camberline = Camberline::NacaA4K6;
    gv.max_camber_position = 0.4;
    gv.chord = 0.06;
    gv.solidity = 1.0;
    gv.stagger = 12.0 * kDeg;
    gv.camber = -25.0 * kDeg;
    gv.thickness_ratio = 0.10;
    gv.throat_over_pitch = 0.85;
    gv.pitch_over_height = 0.3;
    const InletState v0i = inlet(0.0, 0.3), v1i = inlet(8.0 * kDeg, 0.3);
    const CascadePerformance v0 = evaluate(gv, v0i, discharge(v0i), gas, opt);
    const CascadePerformance v1 = evaluate(gv, v1i, discharge(v1i), gas, opt);
    CHECK(v0.frame_sign == -1);
    CHECK(v0.deviation_slope == 0.0);
    CHECK(v0.deviation_deg == Approx(v1.deviation_deg).margin(1e-13));  // incidence does not change it (p. 144)
    CHECK(v0.beta2 > 0.0);                                              // the vane turns the flow to positive swirl
    const BladeAngles va = blade_angles(Camberline::NacaA4K6, -12.0, 25.0, 0.4);
    CHECK(v0.kappa2_deg == Approx(va.kappa2_deg).margin(1e-13));
    CHECK(v0.design_deviation_deg ==
          Approx(howell_design_deviation_deg(1.0, 1.0, 0.4, va.kappa2_deg, 25.0, v0.design_beta1_deg, 1.0)).margin(1e-13));

    // zero camber, zero stagger, axial inflow
    CascadeGeometry gz = g;
    gz.camber = 0.0;
    gz.stagger = 0.0;
    const InletState iz = inlet(0.0, 0.3);
    const CascadePerformance z = evaluate(gz, iz, discharge(iz), gas, opt);
    CHECK(z.frame_sign == 1);
    CHECK(z.kappa1_deg == 0.0);

    // Mach number: i_c rises toward i*, the loss at the same incidence grows; supersonic: shock loss
    const InletState i3 = inlet(d.design_beta1_deg * kDeg, 0.3), i85 = inlet(d.design_beta1_deg * kDeg, 0.85);
    const double l3 = evaluate(g, i3, discharge(i3), gas, opt).loss;
    const CascadePerformance m85 = evaluate(g, i85, discharge(i85), gas, opt);
    CHECK(m85.negative_stall_incidence_deg > d.negative_stall_incidence_deg);
    CHECK(m85.loss > l3);
    CHECK(m85.losses.shock == 0.0);
    CascadeGeometry gs = g;
    gs.throat_over_pitch = 0.95;
    const InletState i115 = inlet(d.design_beta1_deg * kDeg, 1.15);
    const CascadePerformance m115 = evaluate(gs, i115, discharge(i115), gas, opt);
    CHECK(m115.losses.shock > 0.0);
    CHECK(m115.loss == Approx(m115.losses.shock + m115.losses.minimum + m115.losses.incidence).margin(1e-15));
    // a narrow throat at high Mach: the choke limit 6-54 raises i_c
    CascadeGeometry gc = g;
    gc.throat_over_pitch = 0.69;
    const CascadePerformance mc = evaluate(gc, i85, discharge(i85), gas, opt);
    CHECK(mc.choke_limit_active);
    CHECK(mc.negative_stall_incidence_deg == Approx(mc.choke_beta1_deg - mc.kappa1_deg + 1.0).margin(1e-12));
    CHECK(mc.minimum_loss_incidence_deg > mc.design_incidence_deg);  // Fig. 6-15: i_m > i* when 6-54 acts
    CHECK(mc.losses.minimum > mc.losses.design);                     // Eq. 6-58
    // narrower still: i_c of Eq. 6-54 passes i_s, no low-loss range is left: refused
    gc.throat_over_pitch = 0.62;
    CHECK_THROWS_AS(evaluate(gc, i85, discharge(i85), gas, opt), CorrelationError);

    // high loading: D_eq > 2 gives a wake blockage and K_W < 1 (Eqs. 6-84, 7-6)
    CascadeGeometry gh = naca_cascade();
    gh.solidity = 0.6;
    gh.camber = 45.0 * kDeg;
    const InletState ih = inlet(60.0 * kDeg, 0.3);
    const CascadePerformance h = evaluate(gh, ih, discharge(ih), gas, opt);
    CHECK(h.d_eq > 2.0);
    CHECK(h.k_w == Approx(std::pow(2.0 / h.d_eq, 0.9)).margin(1e-15));
    CHECK(h.k_w < 1.0);

    // options: the 2D calibration, the D-factor model, Reynolds, DCA without alpha
    CorrelationOptions o2;
    o2.calibration = LossCalibration::TwoDimensionalCascade;
    g = naca_cascade();
    const InletState i4 = inlet(d.design_beta1_deg * kDeg, 0.4);
    const CascadePerformance c2 = evaluate(g, i4, discharge(i4), gas, o2);
    CHECK(c2.losses.design * cosd(c2.design_beta2_deg) / (2.0 * g.solidity) *
              std::pow(cosd(c2.design_beta2_deg) / cosd(c2.design_beta1_deg), 2) ==
          Approx(0.004 * (1.0 + 3.1 * std::pow(c2.d_eq_design - 1.0, 2) + 0.4 * std::pow(c2.d_eq_design - 1.0, 8))).margin(1e-15));
    o2.loss = ProfileLossModel::DFactor;
    CHECK(evaluate(g, i4, discharge(i4), gas, o2).losses.k1 == 0.0035);
    CorrelationOptions bad = opt;
    bad.loss = ProfileLossModel::DFactor;
    CHECK_THROWS_AS(evaluate(g, i4, discharge(i4), gas, bad), CorrelationError);  // K1 = 0.0073 belongs to Eq. 6-46
    CorrelationOptions re = opt;
    re.reynolds_correction = true;
    CHECK_THROWS_AS(evaluate(g, i4, discharge(i4), gas, re), CorrelationError);  // no Re_c without rho1, mu1
    InletState ir = i4;
    ir.density = 1.2;
    ir.viscosity = 1.8e-5;
    const double re_c = 1.2 * std::hypot(ir.wm, ir.wtheta) * g.chord / 1.8e-5;  // Re_c = rho1 W1 c / mu1
    const CascadePerformance pre = evaluate(g, ir, discharge(ir), gas, re);
    CHECK(pre.reynolds == Approx(re_c).margin(1e-9));
    CHECK(pre.losses.k2 ==
          Approx(1.0 + 0.4 * cosd(d.design_beta2_deg) + 0.004 * reynolds_factor(re_c) / 0.0073).margin(1e-13));
    CascadeGeometry gd = g;
    gd.profile = ProfileFamily::DoubleCircularArc;
    gd.camberline = Camberline::CircularArc;
    const InletState i5 = inlet((d.design_beta1_deg + 3.0) * kDeg, 0.4);
    CHECK_THROWS_AS(evaluate(gd, i5, discharge(i5), gas, opt), CorrelationError);  // alpha of Eq. 6-38 unknown for DCA
    CorrelationOptions oa = opt;
    oa.lieblein_alpha = 0.007;
    CHECK_NOTHROW(evaluate(gd, i5, discharge(i5), gas, oa));

    // refusals and notes
    CascadeGeometry gb = naca_cascade();
    gb.solidity = 0.0;
    CHECK_THROWS_AS(evaluate(gb, i4, discharge(i4), gas, opt), CorrelationError);
    gb = naca_cascade();
    gb.throat_over_pitch = 0.0;
    CHECK_THROWS_AS(evaluate(gb, i4, discharge(i4), gas, opt), CorrelationError);
    gb = naca_cascade();
    gb.pitch_over_height = 0.0;
    CHECK_THROWS_AS(evaluate(gb, i4, discharge(i4), gas, opt), CorrelationError);  // Eq. 6-99 needs s/h
    InletState ib = i4;
    ib.mach = 0.0;
    CHECK_THROWS_AS(evaluate(g, ib, discharge(i4), gas, opt), CorrelationError);
    ib = i4;
    ib.wm = 0.0;
    CHECK_THROWS_AS(evaluate(g, ib, discharge(i4), gas, opt), CorrelationError);
    DischargeEstimate ob = discharge(i4);
    ob.wm = -1.0;
    CHECK_THROWS_AS(evaluate(g, i4, ob, gas, opt), CorrelationError);
    gb = naca_cascade();
    gb.solidity = 2.5;
    const CascadePerformance x = evaluate(gb, i4, discharge(i4), gas, opt);
    CHECK(x.extrapolated);
    CHECK_FALSE(x.notes.empty());
    CHECK_FALSE(d.extrapolated);
}

TEST_CASE("the row-level leakage call and the p. 147 distribution", "[turbolenza][correlations]")
{
    RowWallData row;
    row.clearance = 0.0005;
    row.gap_radius = 0.25;
    row.blade_count = 30;
    row.chord = 0.05;
    row.stagger = 40.0 * kDeg;
    row.row_number = 2;
    WallState a{0.249, 1.2, 150.0, 20.0, 9.5e4, 300.0}, b{0.248, 1.3, 145.0, 180.0, 1.05e5, 315.0};
    const RowLeakage t = endwall_leakage(row, a, b, 30.0);
    TipClearanceInput ti;
    ti.clearance = 0.0005;
    ti.tip_radius = 0.25;
    ti.blade_count = 30;
    ti.chord = 0.05;
    ti.stagger_rad = 40.0 * kDeg;
    ti.row_number = 2;
    ti.r1 = 0.249;
    ti.r2 = 0.248;
    ti.rho1 = 1.2;
    ti.rho2 = 1.3;
    ti.cm1 = 150.0;
    ti.cm2 = 145.0;
    ti.ctheta1 = 20.0;
    ti.ctheta2 = 180.0;
    ti.mean_density = 1.25;  // (rho1 + rho2)/2
    ti.mass_flow = 30.0;
    CHECK(t.total_pressure_loss == Approx(tip_clearance_loss(ti).total_pressure_loss).margin(1e-12));
    CHECK(t.seal_flow_direction == 0);
    // s_row = -1 with mirrored C_theta: the same loss
    RowWallData rm = row;
    rm.row_sign = -1;
    WallState am = a, bm = b;
    am.ctheta = -a.ctheta;
    bm.ctheta = -b.ctheta;
    CHECK(endwall_leakage(rm, am, bm, 30.0).total_pressure_loss == Approx(t.total_pressure_loss).margin(1e-12));

    // the seal: P_R = 9.5e4/1.05e5, Delta P = 1e4 Pa, the higher-pressure side (QN2) state
    row.shrouded = true;
    row.gap_side = EndWall::Hub;
    row.gap_radius = 0.15;
    row.clearance = 0.0004;
    row.fins = 3;
    row.fin_thickness = 0.0005;
    row.fin_pitch = 0.004;
    const Gas gas{1.4, 287.0};
    const RowLeakage s = endwall_leakage(row, a, b, 30.0, gas);
    ShroudSealInput si;
    si.seal_radius = 0.15;
    si.clearance = 0.0004;
    si.fins = 3;
    si.fin_thickness = 0.0005;
    si.fin_pitch = 0.004;
    si.pressure_ratio = 9.5e4 / 1.05e5;
    si.upstream_density = 1.3;
    si.upstream_temperature = 315.0;
    si.gas_constant = 287.0;
    si.pressure_difference = 1.0e4;
    si.mass_flow = 30.0;
    CHECK(s.leakage_mass_flow == Approx(shroud_seal_loss(si).leakage_mass_flow).margin(1e-15));
    CHECK(s.total_pressure_loss == Approx(shroud_seal_loss(si).total_pressure_loss).margin(1e-11));
    CHECK(s.seal_flow_direction == 1);  // from the discharge back to the inlet

    // the distribution: 0 at the opposite wall, 2 Delta P_t at the gap, mass-weighted mean Delta P_t
    CHECK(leakage_loss_at(100.0, 0.0, EndWall::Shroud) == 0.0);
    CHECK(leakage_loss_at(100.0, 1.0, EndWall::Shroud) == 200.0);
    CHECK(leakage_loss_at(100.0, 1.0, EndWall::Hub) == 0.0);
    CHECK(leakage_loss_at(100.0, 0.5, EndWall::Hub) == 100.0);  // half the gap value at mid-passage (p. 147)
    {
        double mean = 0.0;
        const int n = 1000;
        for (int k = 0; k < n; ++k) mean += leakage_loss_at(100.0, (k + 0.5) / n, EndWall::Hub) / n;
        CHECK(mean == Approx(100.0).margin(1e-10));
    }
    CHECK_THROWS_AS(leakage_loss_at(100.0, 1.1, EndWall::Hub), CorrelationError);
    CHECK(loss_coefficient_from_pressure(500.0, 1.2e5, 1.0e5) == Approx(0.025).margin(1e-16));
    CHECK_THROWS_AS(loss_coefficient_from_pressure(500.0, 1.0e5, 1.0e5), CorrelationError);
}
