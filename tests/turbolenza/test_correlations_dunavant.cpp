// The Aungier (2003) Ch. 6 models (Correlations.hpp) on the NACA 65-(12A10)10 cascade of Dunavant et
// al. (1955) as Aungier uses it in Ch. 5: beta_in = 45 deg, alpha = 12 deg (stagger 33 deg), solidity 1,
// M_in = 0.65 (Fig. 5-6, p. 91) and a low-speed M_in = 0.3.  The throat o/s is measured on BladeGen's
// NACA 65-(12A10)10 section (section_from_layer + make_cascade, the geometric throat of Sec. 4.9).
//
// What the book states for this cascade (read from the page images):
//  * p. 91: "For this particular cascade test, those correlations [Ch. 6] indicate that the angle of
//    attack is within 0.5 deg of the optimum value" (alpha = 12 deg);
//  * Fig. 5-6, p. 91: measured discharge angle beta_out = 25.0 deg at M_in = 0.65 (the potential flow
//    gives 24.4 deg); Fig. 5-7, p. 92: 25.1 deg at M_in = 0.757.
// Ch. 6 itself gives no number for this cascade except Fig. 6-2 (65-(12)10 at sigma = 1, where the
// Eq. 6-8 circle lies on Herrig's line; not digitized).
#include "turbolenza/BladeGenSource.hpp"
#include "turbolenza/Cascade.hpp"
#include "turbolenza/Correlations.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>

using namespace BladeGen::Turbolenza;
namespace D = BladeGen::Design;

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;

double dunavant_throat_over_pitch()
{
    D::Layer l;
    l.method = D::Method::Naca65;
    l.naca.cl0 = 1.2;
    l.naca.t_max_over_c = 0.10;
    D::Evaluator ev;
    const BladeSection b = section_from_layer(l, ev);
    const Cascade c = make_cascade(b, 1.0, 33.0);
    return c.throat_over_pitch;
}

CascadeGeometry dunavant_geometry(double o_s)
{
    CascadeGeometry g;
    g.profile = ProfileFamily::Naca;
    g.camberline = Camberline::Naca65;
    g.kind = BladeRowKind::Stator;
    g.chord = 0.127;  // 5 in; the loss does not depend on it
    g.solidity = 1.0;
    g.stagger = 33.0 * kDeg;
    g.camber = naca_camber_from_cl0_deg(1.2, 0.5) * kDeg;  // Eq. 6-7: 30.159 deg
    g.thickness_ratio = 0.10;
    g.max_camber_position = 0.5;
    g.throat_over_pitch = o_s;
    g.pitch_over_height = 0.0;  // a two-dimensional cascade: no end-wall term
    return g;
}

void print(const char* tag, const CascadePerformance& p)
{
    std::printf("[correlations] %s: alpha* %.4f i* %.4f delta* %.4f beta1* %.4f beta2* %.4f | i %.4f delta %.4f beta2 %.4f deg\n"
                "               D* %.4f D_eq* %.4f D %.4f D_eq %.4f | omega* %.5f omega_m %.5f omega %.5f | i_c %.4f i_m %.4f "
                "i_s %.4f choke beta1 %.4f M'_c %.4f xi %.4f K_W %.4f\n",
                tag, p.design_attack_deg, p.design_incidence_deg, p.design_deviation_deg, p.design_beta1_deg,
                p.design_beta2_deg, p.incidence_deg, p.deviation_deg, p.beta2_deg, p.d_factor_design, p.d_eq_design,
                p.d_factor, p.d_eq, p.losses.design, p.losses.minimum, p.loss, p.negative_stall_incidence_deg,
                p.minimum_loss_incidence_deg, p.positive_stall_incidence_deg, p.choke_beta1_deg, p.critical_mach, p.xi, p.k_w);
    for (const auto& n : p.notes) std::printf("               note: %s\n", n.c_str());
}

}  // namespace

TEST_CASE("Dunavant's NACA 65-(12A10)10 cascade: design incidence, deviation, loss, D (Ch. 6 vs the book's statements)",
          "[turbolenza][correlations][dunavant]")
{
    const double o_s = dunavant_throat_over_pitch();
    std::printf("[correlations] NACA 65-(12A10)10, sigma 1, stagger 33 deg: geometric throat o/s = %.5f (BladeGen)\n", o_s);
    REQUIRE(o_s > 0.4);
    REQUIRE(o_s < 0.8);

    CorrelationOptions two_d;  // the planar cascade correlations: Eq. 6-37 (K1 = 0.004, K2 = 1)
    two_d.calibration = LossCalibration::TwoDimensionalCascade;

    // QN1: beta'_1 = 45 deg (W_theta1 = W_m1), a stationary planar cascade; QN2: W_m2 = W_m1, r2 = r1
    // (AVDR 1 assumed: Dunavant's AVDR is not in the book)
    CascadeInletState in;
    in.radius = 1.0;
    in.wm = 100.0;
    in.wtheta = 100.0;
    const DischargeEstimate out{1.0, 100.0};
    const Gas air;  // gamma 1.4

    // the design point is a property of the cascade: the hand values (independent arithmetic from the
    // equations as read): theta = 30.15937, kappa1 = 48.07968, alpha* = 12.55747, i* = -2.52221,
    // beta1* = 45.55747, (delta*_0)_10 = 1.02677, m = 0.22394, delta* = 7.78075, beta2* = 25.70107,
    // W_max* / W1 = 1.28101 (6-35 = 6-39 here), D_eq* = 1.64851, D* = 0.41141, omega* = 0.012418 (6-37)
    for (const double M : {0.3, 0.65}) {
        in.mach = M;
        const CascadePerformance p = evaluate(dunavant_geometry(o_s), in, out, air, two_d);
        print(M < 0.5 ? "M 0.30" : "M 0.65", p);
        CHECK(p.kappa1_deg == Catch::Approx(48.07968).margin(1e-5));
        CHECK(p.design_attack_deg == Catch::Approx(12.55747).margin(1e-5));
        CHECK(p.design_incidence_deg == Catch::Approx(-2.52221).margin(1e-5));
        CHECK(p.design_deviation_deg == Catch::Approx(7.78075).margin(1e-5));
        CHECK(p.design_beta2_deg == Catch::Approx(25.70107).margin(1e-5));
        CHECK(p.wmax_ratio_design == Catch::Approx(1.28101).margin(1e-5));
        CHECK(p.d_eq_design == Catch::Approx(1.64851).margin(1e-5));
        CHECK(p.d_factor_design == Catch::Approx(0.41141).margin(1e-5));
        CHECK(p.losses.design == Catch::Approx(0.012418).margin(1e-6));
        // the operating point: i = 45 - kappa1 = -3.07968, slope (6-76) 0.083624, delta = 7.73413,
        // beta2 = 25.65445, D = 0.39931, D_eq = 1.62985 (i < i*: no alpha term)
        CHECK(p.incidence_deg == Catch::Approx(-3.07968).margin(1e-5));
        CHECK(p.deviation_deg == Catch::Approx(7.73413).margin(1e-5));
        CHECK(p.beta2 == Catch::Approx(25.65445 * kDeg).margin(1e-7));
        CHECK(p.d_factor == Catch::Approx(0.39931).margin(1e-5));
        CHECK(p.d_eq == Catch::Approx(1.62985).margin(1e-5));
        CHECK(p.k_w == 1.0);
        // stall, Mach and loss (hand values, the same arithmetic; agreed by the independent review):
        // M 0.65: i_c -9.8073, i_s 4.8055, i_m -2.5222, choke beta1 34.978, M'_c 0.7421, omega 0.012490;
        // M 0.30: i_c -10.6973, i_s 5.7006, no choke angle, M'_c 0.7190, omega 0.012475
        CHECK(p.minimum_loss_incidence_deg == Catch::Approx(-2.5222).margin(1e-4));
        CHECK(p.losses.minimum == Catch::Approx(p.losses.design).margin(1e-15));  // i_m = i*
        if (M > 0.5) {
            CHECK(p.negative_stall_incidence_deg == Catch::Approx(-9.8073).margin(1e-4));
            CHECK(p.positive_stall_incidence_deg == Catch::Approx(4.8055).margin(1e-4));
            CHECK(p.choke_beta1_deg == Catch::Approx(34.978).margin(1e-3));
            CHECK(p.critical_mach == Catch::Approx(0.7421).margin(1e-4));
            CHECK(p.loss == Catch::Approx(0.012490).margin(1e-6));
            CHECK_FALSE(p.choke_limit_active);
        } else {
            CHECK(p.negative_stall_incidence_deg == Catch::Approx(-10.6973).margin(1e-4));
            CHECK(p.positive_stall_incidence_deg == Catch::Approx(5.7006).margin(1e-4));
            CHECK(std::isnan(p.choke_beta1_deg));
            CHECK(p.critical_mach == Catch::Approx(0.7190).margin(1e-4));
            CHECK(p.loss == Catch::Approx(0.012475).margin(1e-6));
        }

        // the book: alpha = 12 deg is "within 0.5 deg of the optimum" (p. 91).  Eq. 6-8 gives
        // alpha* = 12.557 deg: 0.557 deg away (measured), just outside the stated 0.5 deg; the bound
        // records that.  The measured beta_out = 25.0 deg (Fig. 5-6): Eq. 6-77 gives 25.654 deg,
        // 0.65 deg high (measured).
        CHECK(std::abs(p.design_attack_deg - 12.0) < 0.6);
        CHECK(std::abs(p.beta2_deg - 25.0) < 0.7);
    }

    // the Lieblein design incidence (6-12 .. 6-18) for comparison
    CorrelationOptions lie = two_d;
    lie.incidence = DesignIncidenceModel::Lieblein;
    in.mach = 0.3;
    const CascadePerformance pl = evaluate(dunavant_geometry(o_s), in, out, air, lie);
    print("M 0.30 Lieblein", pl);
    CHECK(std::abs(pl.design_attack_deg - 12.0) < 3.0);  // a different, less subjective optimum (p. 124)
}
