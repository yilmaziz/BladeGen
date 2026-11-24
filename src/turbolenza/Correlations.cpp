// The empirical cascade performance models of Aungier (2003) Ch. 6 (see Correlations.hpp for the
// sources, the conventions, the refusals and the two corrected misprints).  Every function cites the
// equation it implements; every choice the book does not make is marked "ours".
#include "turbolenza/Correlations.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <numbers>
#include <string>

namespace BladeGen::Turbolenza {

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kDeg = kPi / 180.0;

[[nodiscard]] double rad(double deg) noexcept { return deg * kDeg; }
[[nodiscard]] double deg(double r) noexcept { return r / kDeg; }
[[nodiscard]] double nan() noexcept { return std::numeric_limits<double>::quiet_NaN(); }

[[noreturn]] void refuse(const std::string& what) { throw CorrelationError(what); }

void require_finite(double x, const char* name, const char* eq)
{
    if (!std::isfinite(x)) refuse(std::format("{}: {} = {} is not finite", eq, name, x));
}

void require_positive(double x, const char* name, const char* eq)
{
    if (!(std::isfinite(x) && x > 0.0)) refuse(std::format("{}: {} = {} must be finite and > 0", eq, name, x));
}

/// |angle| < 90 deg (a flow or blade angle measured from the axial direction).
void require_angle(double a_deg, const char* name, const char* eq)
{
    if (!(std::isfinite(a_deg) && std::abs(a_deg) < 90.0))
        refuse(std::format("{}: {} = {} deg must lie in (-90, 90) deg", eq, name, a_deg));
}

void require_thickness(double tb_c, const char* eq)
{
    if (!(std::isfinite(tb_c) && tb_c >= 0.0 && tb_c < 1.0))
        refuse(std::format("{}: t_b/c = {} must lie in [0, 1)", eq, tb_c));
}

/// beta1 >= 0 where the correlation raises beta1 to a non-integer power (6-13, 6-15, 6-20, 6-76).
void require_beta1_nonnegative(double beta1_deg, const char* eq)
{
    if (!(std::isfinite(beta1_deg) && beta1_deg >= 0.0 && beta1_deg < 90.0))
        refuse(std::format("{}: the inlet flow angle beta1 = {} deg (positive-camber cascade frame) must lie in [0, 90) deg",
                           eq, beta1_deg));
}

void require_gamma(double gamma)
{
    if (!(std::isfinite(gamma) && gamma > 1.0))
        refuse(std::format("perfect gas: the ratio of specific heats gamma = {} must be > 1", gamma));
}

}  // namespace

// ================================================================================================
// Ch. 4: blade angles; Eq. 6-7
// ================================================================================================

BladeAngles blade_angles(Camberline camberline, double stagger_deg, double camber_deg, double a_c)
{
    require_angle(stagger_deg, "stagger gamma", "Eqs. 4-32, 4-33");
    if (!(std::isfinite(camber_deg) && camber_deg >= 0.0 && camber_deg < 90.0))
        refuse(std::format("Ch. 4 blade angles: the camber theta = {} deg (positive-camber frame) must lie in [0, 90) deg",
                           camber_deg));
    if (!std::isfinite(a_c)) refuse(std::format("Ch. 4 blade angles: a/c = {} is not finite", a_c));
    BladeAngles b;
    switch (camberline) {
    case Camberline::Naca65:
    case Camberline::CircularArc:
        // the circular arc and the NACA 65 equivalent circular arc: a/c = 0.5 (p. 120, p. 121)
        if (a_c != 0.5)
            refuse(std::format("Eq. 4-35: a circular-arc (or NACA 65 equivalent circular-arc) camberline has a/c = 0.5, "
                               "got {}",
                               a_c));
        b.chi1_deg = b.chi2_deg = 0.5 * camber_deg;  // [4-35]
        break;
    case Camberline::ParabolicArc:
    case Camberline::NacaA4K6: {
        // the parabolic arc (A4K6 as its equivalent parabolic arc, p. 66): 0.25 < a/c < 0.75 [4-20]
        if (!(a_c > 0.25 && a_c < 0.75))
            refuse(std::format("Eq. 4-20: a parabolic-arc camberline needs 0.25 < a/c < 0.75, got {}", a_c));
        if (camber_deg == 0.0) break;  // b = 0: chi1 = chi2 = 0 (the limit of Eq. 4-19)
        const double t4 = 4.0 * std::tan(rad(camber_deg));
        const double b_c = (std::sqrt(1.0 + t4 * t4 * (a_c - a_c * a_c - 3.0 / 16.0)) - 1.0) / t4;  // [4-19]
        b.chi1_deg = deg(std::atan(4.0 * b_c / (4.0 * a_c - 1.0)));                                // [4-17]
        b.chi2_deg = deg(std::atan(4.0 * b_c / (3.0 - 4.0 * a_c)));                                // [4-18]
        break;
    }
    }
    b.kappa1_deg = b.chi1_deg + stagger_deg;  // [4-32]
    b.kappa2_deg = stagger_deg - b.chi2_deg;  // [4-33]
    return b;
}

double naca_camber_from_cl0_deg(double cl0, double a_c)
{
    require_finite(cl0, "C_l0", "Eq. 6-7");
    require_positive(a_c, "a/c", "Eq. 6-7");
    return 4.0 * deg(std::atan(0.05515 * cl0 / a_c));  // [6-7] = [4-29]
}

double profile_shape_factor(ProfileFamily profile) noexcept
{
    switch (profile) {  // p. 121
    case ProfileFamily::Naca: return 1.0;
    case ProfileFamily::C4: return 1.1;
    case ProfileFamily::DoubleCircularArc: return 0.7;
    }
    return 1.0;
}

// ================================================================================================
// Sec. 6.2: design incidence
// ================================================================================================

double incidence_thickness_factor(double tb_c)
{
    require_thickness(tb_c, "Eqs. 6-10, 6-11");
    const double q = 0.28 / (0.1 + std::pow(tb_c, 0.3));  // [6-11]
    return std::pow(10.0 * tb_c, q);                      // [6-10]
}

double design_attack_angle_deg(double ksh, double kti, double camber_deg, double a_c, double solidity)
{
    require_positive(solidity, "solidity sigma", "Eq. 6-8");
    require_positive(a_c, "a/c", "Eq. 6-8");
    require_finite(camber_deg, "theta", "Eq. 6-8");
    const double e = 0.65 - 0.002 * camber_deg;                                                 // [6-9]
    return (3.6 * ksh * kti + 0.3532 * camber_deg * std::pow(a_c, 0.25)) * std::pow(solidity, e);  // [6-8]
}

double lieblein_zero_camber_incidence_deg(double beta1_deg, double solidity)
{
    require_beta1_nonnegative(beta1_deg, "Eq. 6-13");
    require_positive(solidity, "solidity sigma", "Eq. 6-13");
    const double s3 = solidity * solidity * solidity;
    const double p = 0.914 + s3 / 160.0;  // [6-14]
    return std::pow(beta1_deg, p) / (5.0 + 46.0 * std::exp(-2.3 * solidity)) -
           0.1 * s3 * std::exp((beta1_deg - 70.0) / 4.0);  // [6-13]
}

double lieblein_incidence_slope(double beta1_deg, double solidity)
{
    require_beta1_nonnegative(beta1_deg, "Eq. 6-15");
    require_positive(solidity, "solidity sigma", "Eq. 6-15");
    return 0.025 * solidity - 0.06 - std::pow(beta1_deg / 90.0, 1.0 + 1.2 * solidity) / (1.5 + 0.43 * solidity);  // [6-15]
}

LiebleinIncidence lieblein_design_incidence(double ksh, double kti, double stagger_deg, double camber_deg, double a_c,
                                            double kappa1_deg, double solidity)
{
    require_angle(stagger_deg, "stagger gamma", "Eq. 6-16");
    require_finite(camber_deg, "theta", "Eq. 6-16");
    require_finite(kappa1_deg, "kappa1", "Eq. 6-18");
    LiebleinIncidence r;
    r.pseudo_kappa1_deg = std::abs(stagger_deg) + 0.5 * camber_deg;  // [6-16]
    // [6-12] with beta1 = kappa1_bar + i*_bar: a fixed point (p. 123: "an iterative solution is needed").
    // Ours: plain substitution from i*_bar = 0 to |change| <= 1e-12 deg, at most 200 steps (the map's
    // slope K_sh K_t,i d(i0)/d(beta1) + theta dn/d(beta1) is ~0.1 in the figures' range).
    double i_bar = 0.0;
    for (int k = 1; k <= 200; ++k) {
        const double beta1 = r.pseudo_kappa1_deg + i_bar;
        const double next = ksh * kti * lieblein_zero_camber_incidence_deg(beta1, solidity) +
                            lieblein_incidence_slope(beta1, solidity) * camber_deg;  // [6-12]
        const double change = std::abs(next - i_bar);
        i_bar = next;
        r.iterations = k;
        if (change <= 1e-12) {
            r.pseudo_incidence_deg = i_bar;
            r.attack_deg = r.pseudo_kappa1_deg + i_bar - std::abs(stagger_deg) + (a_c - 0.5) * camber_deg;  // [6-17]
            r.incidence_deg = r.attack_deg + stagger_deg - kappa1_deg;                                     // [6-18]
            return r;
        }
    }
    refuse(std::format("Eq. 6-12: the design-incidence iteration on beta1 did not converge in 200 steps (kappa1_bar = {} "
                       "deg, sigma = {}, theta = {} deg)",
                       r.pseudo_kappa1_deg, solidity, camber_deg));
}

// ================================================================================================
// Sec. 6.3: design deviation
// ================================================================================================

double lieblein_zero_camber_deviation_deg(double beta1_deg, double solidity)
{
    require_beta1_nonnegative(beta1_deg, "Eq. 6-20");
    require_positive(solidity, "solidity sigma", "Eq. 6-20");
    return 0.01 * solidity * beta1_deg +
           (0.74 * std::pow(solidity, 1.9) + 3.0 * solidity) * std::pow(beta1_deg / 90.0, 1.67 + 1.09 * solidity);  // [6-20]
}

double lieblein_deviation_slope(Camberline camberline, double beta1_deg, double solidity)
{
    require_finite(beta1_deg, "beta1", "Eqs. 6-21 .. 6-24");
    require_positive(solidity, "solidity sigma", "Eq. 6-21");
    const double x = beta1_deg / 100.0;
    double m10 = 0.0;
    switch (camberline) {
    case Camberline::Naca65: m10 = 0.17 - 0.0333 * x + 0.333 * x * x; break;                          // [6-22]
    case Camberline::CircularArc: m10 = 0.249 + 0.074 * x - 0.132 * x * x + 0.316 * x * x * x; break;  // [6-23]
    case Camberline::ParabolicArc:
    case Camberline::NacaA4K6:
        refuse("Eqs. 6-21 .. 6-23: Lieblein's deviation slope exists for the NACA 65 and circular-arc camberlines only "
               "(p. 126); the parabolic arc and A4K6 use Eq. 6-26");
    }
    const double b = 0.9625 - 0.17 * x - 0.85 * x * x * x;  // [6-24]
    return m10 / std::pow(solidity, b);                     // [6-21]
}

double deviation_thickness_factor(double tb_c)
{
    require_thickness(tb_c, "Eq. 6-25");
    return 6.25 * tb_c + 37.5 * tb_c * tb_c;  // [6-25]
}

double lieblein_design_deviation_deg(Camberline camberline, double ksh, double ktd, double beta1_deg, double camber_deg,
                                     double solidity)
{
    require_finite(camber_deg, "theta", "Eq. 6-19");
    return ksh * ktd * lieblein_zero_camber_deviation_deg(beta1_deg, solidity) +
           lieblein_deviation_slope(camberline, beta1_deg, solidity) * camber_deg;  // [6-19]
}

double howell_design_deviation_deg(double ksh, double ktd, double a_c, double kappa2_deg, double camber_deg,
                                   double beta1_deg, double solidity)
{
    require_positive(solidity, "solidity sigma", "Eq. 6-26");
    require_finite(a_c, "a/c", "Eq. 6-26");
    require_finite(kappa2_deg, "kappa2", "Eq. 6-26");
    require_finite(camber_deg, "theta", "Eq. 6-26");
    const double rs = std::sqrt(solidity);
    const double denom = 1.0 - 0.002 * camber_deg / rs;
    if (!(denom > 0.0))
        refuse(std::format("Eq. 6-26: 1 - 0.002 theta/sqrt(sigma) = {} must be > 0 (theta = {} deg, sigma = {})", denom,
                           camber_deg, solidity));
    const double corr = ksh * ktd - 1.0;
    // (delta*_0)_10 only when its factor is not zero (ours: a 10 % NACA profile then needs no beta1 >= 0)
    const double d0 = corr != 0.0 ? lieblein_zero_camber_deviation_deg(beta1_deg, solidity) : 0.0;
    return (0.92 * a_c * a_c + 0.002 * kappa2_deg) / denom * camber_deg / rs + corr * d0;  // [6-26]
}

// ================================================================================================
// Sec. 6.4: diffusion factors and the design loss
// ================================================================================================

double velocity_ratio_w2_w1(double beta1_deg, double beta2_deg, double wm_ratio)
{
    require_angle(beta1_deg, "beta1", "W = W_m / cos(beta)");
    require_angle(beta2_deg, "beta2", "W = W_m / cos(beta)");
    require_positive(wm_ratio, "W_m2/W_m1", "W = W_m / cos(beta)");
    return wm_ratio * std::cos(rad(beta1_deg)) / std::cos(rad(beta2_deg));
}

double diffusion_factor(double beta1_deg, double beta2_deg, double solidity, double wm_ratio)
{
    require_positive(solidity, "solidity sigma", "Eq. 6-31");
    const double w2_w1 = velocity_ratio_w2_w1(beta1_deg, beta2_deg, wm_ratio);
    // (W_theta1 - W_theta2)/W1 = sin(beta1) - (W_m2/W_m1) cos(beta1) tan(beta2)
    const double dwt = std::sin(rad(beta1_deg)) - wm_ratio * std::cos(rad(beta1_deg)) * std::tan(rad(beta2_deg));
    return 1.0 - w2_w1 + dwt / (2.0 * solidity);  // [6-31]
}

double diffusion_factor_from_velocities(double w_max, double w1, double w2)
{
    require_positive(w1, "W1", "Eq. 6-29");
    require_finite(w_max, "W_max", "Eq. 6-29");
    require_finite(w2, "W2", "Eq. 6-29");
    return (w_max - w2) / w1;  // [6-29]
}

double equivalent_diffusion_from_velocities(double w_max, double w2)
{
    require_positive(w2, "W2", "Eq. 6-30");
    require_finite(w_max, "W_max", "Eq. 6-30");
    return w_max / w2;  // [6-30]
}

double lieblein_wmax_ratio_design(double beta1_deg, double beta2_deg, double solidity)
{
    require_angle(beta1_deg, "beta1*", "Eq. 6-35");
    require_angle(beta2_deg, "beta2*", "Eq. 6-35");
    require_positive(solidity, "solidity sigma", "Eq. 6-35");
    const double c1 = std::cos(rad(beta1_deg));
    return 1.12 + 0.61 * c1 * c1 / solidity * (std::tan(rad(beta1_deg)) - std::tan(rad(beta2_deg)));  // [6-35]
}

double wmax_ratio(double beta1_deg, double beta2_deg, double solidity, double wm_ratio, double radius_ratio,
                  double blade_speed_ratio, double incidence_excess_deg, double alpha)
{
    require_angle(beta1_deg, "beta1", "Eq. 6-39");
    require_angle(beta2_deg, "beta2", "Eq. 6-39");
    require_positive(solidity, "solidity sigma", "Eq. 6-39");
    require_positive(wm_ratio, "W_m2/W_m1", "Eq. 6-39");
    require_positive(radius_ratio, "r2/r1", "Eq. 6-39");
    require_finite(blade_speed_ratio, "U1/W_m1", "Eq. 6-39");
    require_finite(alpha, "alpha", "Eq. 6-39");
    if (!(std::isfinite(incidence_excess_deg) && incidence_excess_deg >= 0.0))
        refuse(std::format("Eq. 6-38/6-39: the term alpha (i - i*)^1.43 is stated for i >= i* (p. 132); i - i* = {} deg",
                           incidence_excess_deg));
    const double c1 = std::cos(rad(beta1_deg));
    // (r1 C_theta1 - r2 C_theta2)/(r1 C_m1) with C_theta = W_m tan(beta) + omega r, C_m = W_m, U1 = omega r1:
    //   = tan b1 - (r2/r1)(W_m2/W_m1) tan b2 + (U1/W_m1)(1 - (r2/r1)^2)
    const double circulation = std::tan(rad(beta1_deg)) - radius_ratio * wm_ratio * std::tan(rad(beta2_deg)) +
                               blade_speed_ratio * (1.0 - radius_ratio * radius_ratio);
    return 1.12 + 0.61 * c1 * c1 / solidity * circulation + alpha * std::pow(incidence_excess_deg, 1.43);  // [6-39]
}

double lieblein_alpha(ProfileFamily profile)
{
    switch (profile) {  // p. 132
    case ProfileFamily::Naca: return 0.0117;
    case ProfileFamily::C4: return 0.007;
    case ProfileFamily::DoubleCircularArc: break;
    }
    refuse("Eq. 6-38: alpha is given for NACA 65 (0.0117) and C.4 (0.007) blades only (p. 132); set "
           "CorrelationOptions::lieblein_alpha for a double-circular-arc profile");
}

double equivalent_diffusion(double wmax_r, double beta1_deg, double beta2_deg, double wm_ratio)
{
    require_finite(wmax_r, "W_max/W1", "Eq. 6-36");
    return wmax_r / velocity_ratio_w2_w1(beta1_deg, beta2_deg, wm_ratio);  // [6-30], [6-36]
}

double design_loss_d_factor(double d_star, double beta2_deg, double solidity, double k1, double k2)
{
    require_angle(beta2_deg, "beta2*", "Eq. 6-45");
    require_positive(solidity, "solidity sigma", "Eq. 6-45");
    require_finite(k1, "K1", "Eq. 6-45");
    require_finite(k2, "K2", "Eq. 6-45");
    if (!(std::isfinite(d_star) && d_star <= 0.6))
        refuse(std::format("Eq. 6-34/6-45: D* = {} exceeds the diffusion limit 0.6, where the correlation \"should not be "
                           "applied\" (p. 130)",
                           d_star));
    const double d4 = d_star * d_star * d_star * d_star;
    return 2.0 * solidity / std::cos(rad(beta2_deg)) * k1 * (k2 + 3.5 * d_star + 37.0 * d4);  // [6-45]
}

double design_loss_equivalent_diffusion(double d_eq_star, double beta2_deg, double solidity, double w1_w2, double k1,
                                        double k2)
{
    require_angle(beta2_deg, "beta2*", "Eq. 6-46");
    require_positive(solidity, "solidity sigma", "Eq. 6-46");
    require_positive(w1_w2, "W1*/W2*", "Eq. 6-46");
    require_finite(d_eq_star, "D_eq*", "Eq. 6-46");
    require_finite(k1, "K1", "Eq. 6-46");
    require_finite(k2, "K2", "Eq. 6-46");
    const double x = d_eq_star - 1.0;
    const double x2 = x * x, x8 = x2 * x2 * x2 * x2;
    return 2.0 * solidity / (std::cos(rad(beta2_deg)) * w1_w2 * w1_w2) * k1 * (k2 + 3.1 * x2 + 0.4 * x8);  // [6-46]
}

double wake_loss_coefficient(double theta_w_over_c, double solidity, double beta2_deg, double w2_w1)
{
    require_finite(theta_w_over_c, "theta_w/c", "Eq. 6-27");
    require_positive(solidity, "solidity sigma", "Eq. 6-27");
    require_angle(beta2_deg, "beta2", "Eq. 6-27");
    require_finite(w2_w1, "W2/W1", "Eq. 6-27");
    return 2.0 * theta_w_over_c * solidity / std::cos(rad(beta2_deg)) * w2_w1 * w2_w1;  // [6-27]
}

double howell_mean_angle_deg(double beta1_deg, double beta2_deg)
{
    require_angle(beta1_deg, "beta1", "Eq. 6-41");
    require_angle(beta2_deg, "beta2", "Eq. 6-41");
    return deg(std::atan(0.5 * (std::tan(rad(beta1_deg)) + std::tan(rad(beta2_deg)))));  // [6-41]
}

double howell_loss_from_drag(double cd, double solidity, double beta1_deg, double beta2_deg)
{
    require_finite(cd, "C_D", "Eq. 6-40");
    require_positive(solidity, "solidity sigma", "Eq. 6-40");
    const double cm = std::cos(rad(howell_mean_angle_deg(beta1_deg, beta2_deg)));
    const double c1 = std::cos(rad(beta1_deg));
    return cd * solidity * c1 * c1 / (cm * cm * cm);  // [6-40]
}

double howell_endwall_drag(double pitch_over_height)
{
    require_positive(pitch_over_height, "s/h", "Eq. 6-42");
    return 0.02 * pitch_over_height;  // [6-42]
}

double howell_lift_coefficient(double beta1_deg, double beta2_deg, double solidity)
{
    require_positive(solidity, "solidity sigma", "Eq. 6-44");
    const double cm = std::cos(rad(howell_mean_angle_deg(beta1_deg, beta2_deg)));
    return 2.0 * cm * (std::tan(rad(beta1_deg)) - std::tan(rad(beta2_deg))) / solidity;  // [6-44]
}

double howell_secondary_drag(double lift_coefficient)
{
    require_finite(lift_coefficient, "C_L", "Eq. 6-43");
    return 0.18 * lift_coefficient * lift_coefficient;  // [6-43]
}

double aspect_ratio_k2(double pitch_over_height, double beta2_deg, double k_re, double k1)
{
    require_positive(pitch_over_height, "s/h", "Eq. 6-99/6-100");
    require_angle(beta2_deg, "beta2", "Eq. 6-99/6-100");
    require_finite(k_re, "K_Re", "Eq. 6-100");
    require_positive(k1, "K1", "Eq. 6-100");
    return 1.0 + pitch_over_height * std::cos(rad(beta2_deg)) + 0.004 * k_re / k1;  // [6-99] (k_re = 0), [6-100]
}

double reynolds_factor(double reynolds)
{
    require_positive(reynolds, "Re_c", "Eqs. 6-101, 6-102");
    if (!(reynolds > 1.0)) refuse(std::format("Eq. 6-102: Re_c = {} must be > 1 (log Re_c > 0)", reynolds));
    constexpr double ref = 2.5e5;
    if (reynolds < ref) return std::sqrt(ref / reynolds) - 1.0;    // [6-101] laminar
    return std::pow(std::log(ref) / std::log(reynolds), 2.58) - 1.0;  // [6-102] turbulent (both 0 at 2.5e5)
}

// ================================================================================================
// Sec. 6.5: stall incidence
// ================================================================================================

double negative_stall_attack_offset_deg(double beta1c_deg, double camber_deg)
{
    require_finite(camber_deg, "theta", "Eq. 6-47");
    if (!(std::isfinite(beta1c_deg) && beta1c_deg >= 20.0))
        refuse(std::format("Eq. 6-47: beta_1c = {} deg; Aungier limits beta_1c >= 20 deg (p. 135; the equation is "
                           "singular at 0)",
                           beta1c_deg));
    return -9.0 + (1.0 - std::pow(30.0 / beta1c_deg, 0.48)) * camber_deg / 4.176;  // [6-47]
}

double positive_stall_attack_offset_deg(double beta1s_deg, double camber_deg)
{
    require_finite(beta1s_deg, "beta_1s", "Eq. 6-48");
    require_finite(camber_deg, "theta", "Eq. 6-48");
    return 10.3 + (2.92 - beta1s_deg / 15.6) * camber_deg / 8.2;  // [6-48]
}

StallRanges stall_ranges(double design_attack_deg, double stagger_deg, double camber_deg)
{
    require_finite(design_attack_deg, "alpha*", "Eqs. 6-47, 6-48");
    require_angle(stagger_deg, "stagger gamma", "Eqs. 6-47, 6-48");
    require_finite(camber_deg, "theta", "Eqs. 6-47, 6-48");
    StallRanges r;
    // positive stall: beta_1s = gamma + alpha* + (alpha_s - alpha*)(beta_1s), linear in beta_1s [6-48]
    r.beta1s_deg = (stagger_deg + design_attack_deg + 10.3 + 2.92 * camber_deg / 8.2) /
                   (1.0 + camber_deg / (15.6 * 8.2));
    r.positive_deg = positive_stall_attack_offset_deg(r.beta1s_deg, camber_deg);  // R_s [6-50]
    // negative stall: beta_1c = gamma + alpha* + (alpha_c - alpha*)(max(beta_1c, 20)) [6-47], p. 135
    // Ours: plain substitution from beta_1c = gamma + alpha* - 9 to |change| <= 1e-12 deg (the map's
    // slope 0.48 (30/b)^0.48 theta/(4.176 b) is < 0.5 for b >= 20 deg, theta < 90 deg).
    double b = stagger_deg + design_attack_deg - 9.0;
    for (int k = 1; k <= 500; ++k) {
        const double used = std::max(b, 20.0);
        const double next = stagger_deg + design_attack_deg + negative_stall_attack_offset_deg(used, camber_deg);
        const double change = std::abs(next - b);
        b = next;
        if (change <= 1e-12) {
            r.beta1c_deg = b;
            r.beta1c_limited = b < 20.0;
            r.negative_deg = -negative_stall_attack_offset_deg(std::max(b, 20.0), camber_deg);  // R_c [6-49]
            return r;
        }
    }
    refuse(std::format("Eq. 6-47: the negative-stall iteration on beta_1c did not converge (alpha* = {} deg, gamma = {} "
                       "deg, theta = {} deg)",
                       design_attack_deg, stagger_deg, camber_deg));
}

// ================================================================================================
// Sec. 6.6: Mach number effects (perfect gas where an equation of state is needed)
// ================================================================================================

double mach_negative_stall_incidence_deg(double design_incidence_deg, double range_negative_deg, double mach1)
{
    require_finite(design_incidence_deg, "i*", "Eq. 6-51");
    require_finite(range_negative_deg, "R_c", "Eq. 6-51");
    require_positive(mach1, "M'_1", "Eq. 6-51");
    return design_incidence_deg - range_negative_deg / (1.0 + 0.5 * mach1 * mach1 * mach1);  // [6-51]
}

double mach_positive_stall_incidence_deg(double design_incidence_deg, double range_positive_deg, double mach1, double ksh)
{
    require_finite(design_incidence_deg, "i*", "Eq. 6-52");
    require_finite(range_positive_deg, "R_s", "Eq. 6-52");
    require_positive(mach1, "M'_1", "Eq. 6-52");
    require_positive(ksh, "K_sh", "Eq. 6-52");
    const double km = std::min(ksh, 1.0) * mach1;  // "with the constraint that K_sh <= 1" (p. 136)
    return design_incidence_deg + range_positive_deg / (1.0 + 0.5 * km * km * km);  // [6-52]
}

double mass_flux_ratio_to_sonic(double mach, double gamma)
{
    require_gamma(gamma);
    if (!(std::isfinite(mach) && mach >= 0.0)) refuse(std::format("isentropic mass flux: M = {} must be >= 0", mach));
    return mach * std::pow((gamma + 1.0) / (2.0 + (gamma - 1.0) * mach * mach), (gamma + 1.0) / (2.0 * (gamma - 1.0)));
}

double sonic_velocity_ratio(double mach, double gamma)
{
    require_gamma(gamma);
    require_positive(mach, "M", "W*/W");
    return std::sqrt((2.0 + (gamma - 1.0) * mach * mach) / ((gamma + 1.0) * mach * mach));
}

double choke_inlet_angle_deg(double throat_over_pitch, double mach1, double gamma)
{
    require_positive(throat_over_pitch, "o/s", "Eq. 6-53");
    require_positive(mach1, "M'_1", "Eq. 6-53");
    // rho1 W1 s cos(beta_choke) = o rho* W*  =>  cos(beta_choke) = (o/s) / [rho1 W1 / (rho* W*)]  [6-53]
    const double c = throat_over_pitch / mass_flux_ratio_to_sonic(mach1, gamma);
    if (c >= 1.0) return nan();  // no inlet angle chokes the passage
    return deg(std::acos(c));
}

double minimum_loss_incidence_deg(double ic_deg, double is_deg, double range_negative_deg, double range_positive_deg)
{
    require_finite(ic_deg, "i_c", "Eq. 6-55");
    require_finite(is_deg, "i_s", "Eq. 6-55");
    if (!(range_negative_deg + range_positive_deg > 0.0))
        refuse(std::format("Eq. 6-55: R_c + R_s = {} deg must be > 0", range_negative_deg + range_positive_deg));
    return ic_deg + (is_deg - ic_deg) * range_negative_deg / (range_negative_deg + range_positive_deg);  // [6-55]
}

double critical_mach_number(double mach1, double sonic_ratio, double wmax_ratio_design)
{
    require_positive(mach1, "M'_1", "Eq. 6-59");
    require_positive(sonic_ratio, "W*/W1", "Eq. 6-59");
    require_positive(wmax_ratio_design, "W_max/W1", "Eq. 6-59");
    return mach1 * sonic_ratio / wmax_ratio_design;  // [6-59]: M'_c = M'_1 W*/W_max
}

double minimum_loss(double design_loss, double im_deg, double design_incidence_deg, double range_positive_deg, double ksh,
                    double mach1, double critical_mach, double sonic_ratio)
{
    require_finite(design_loss, "omega*", "Eq. 6-58");
    require_finite(im_deg, "i_m", "Eq. 6-58");
    require_finite(design_incidence_deg, "i*", "Eq. 6-58");
    require_positive(range_positive_deg, "R_s", "Eq. 6-58");
    const double d = (im_deg - design_incidence_deg) / range_positive_deg;
    double w = design_loss * (1.0 + d * d);  // [6-58]
    if (mach1 > critical_mach) {
        require_positive(critical_mach, "M'_c", "Eq. 6-60");
        // [6-60]: M' is M'_1 limited to 1 (p. 138: "the limit M' <= 1 is imposed"); K_sh as defined (the
        // K_sh <= 1 constraint is stated for Eq. 6-52 only; ours: not applied here)
        const double m = std::min(mach1, 1.0);
        const double e = (m / critical_mach - 1.0) * sonic_ratio;
        if (e > 0.0) w += ksh * e * e;
    }
    return w;
}

double normalized_incidence(double i_deg, double im_deg, double ic_deg, double is_deg)
{
    require_finite(i_deg, "i", "Eqs. 6-78, 6-79");
    if (!(is_deg > im_deg && im_deg > ic_deg))
        refuse(std::format("Eqs. 6-78, 6-79: need i_c < i_m < i_s, got i_c = {}, i_m = {}, i_s = {} deg (no low-loss "
                           "incidence range)",
                           ic_deg, im_deg, is_deg));
    if (i_deg >= im_deg) return (i_deg - im_deg) / (is_deg - im_deg);  // [6-78]
    return (i_deg - im_deg) / (im_deg - ic_deg);                       // [6-79]
}

double off_design_loss_factor(double xi) noexcept
{
    if (xi < -2.0) return 5.0 - 4.0 * (xi + 2.0);  // [6-81]
    if (xi > 1.0) return 2.0 + 2.0 * (xi - 1.0);   // [6-82]
    return 1.0 + xi * xi;                          // [6-80] (= Eqs. 6-56, 6-57 inside [i_c, i_s])
}

// ================================================================================================
// Sec. 6.8: off-design deviation; wake blockage
// ================================================================================================

double deviation_incidence_slope(double beta1_deg, double solidity)
{
    require_beta1_nonnegative(beta1_deg, "Eq. 6-76");
    require_positive(solidity, "solidity sigma", "Eq. 6-76");
    const double s4 = solidity * solidity * solidity * solidity;
    return (1.0 + (solidity + 0.25 * s4) * std::pow(beta1_deg / 53.0, 2.5)) / std::exp(3.1 * solidity);  // [6-76]
}

double off_design_deviation_deg(double design_deviation_deg, double slope, double i_deg, double design_incidence_deg,
                                double wm_ratio, bool inlet_guide_vane)
{
    require_finite(design_deviation_deg, "delta*", "Eq. 6-77");
    require_finite(slope, "(d delta/d i)*", "Eq. 6-77");
    require_finite(i_deg, "i", "Eq. 6-77");
    require_finite(design_incidence_deg, "i*", "Eq. 6-77");
    require_positive(wm_ratio, "W_m2/W_m1", "Eq. 6-77");
    const double incidence_term = inlet_guide_vane ? 0.0 : slope * (i_deg - design_incidence_deg);  // p. 144
    return design_deviation_deg + incidence_term + 10.0 * (1.0 - wm_ratio);                         // [6-77]
}

double wake_blockage(double d_eq)
{
    require_finite(d_eq, "D_eq", "Eqs. 6-83, 6-84");
    if (d_eq <= 2.0) return 0.0;              // [6-83]
    return 1.0 - std::pow(2.0 / d_eq, 0.9);  // [6-84]
}

// ================================================================================================
// Sec. 6.7: shock loss (Swan 1961), perfect gas
// ================================================================================================

ShockGeometry shock_expansion_geometry(double camber_deg, double tb_c, double stagger_deg, double solidity)
{
    if (!(std::isfinite(camber_deg) && camber_deg >= 0.0 && camber_deg < 180.0))
        refuse(std::format("Eq. 6-65: theta = {} deg must lie in [0, 180) deg", camber_deg));
    require_thickness(tb_c, "Eq. 6-65");
    require_angle(stagger_deg, "stagger gamma", "Eq. 6-68");
    require_positive(solidity, "solidity sigma", "Eq. 6-69");
    ShockGeometry g;
    g.suction_camber_deg = 4.0 * deg(std::atan(std::tan(rad(camber_deg) / 4.0) + tb_c));  // [6-65]
    g.psi_deg = 90.0 - 0.5 * g.suction_camber_deg - stagger_deg;                             // [6-68]
    if (!(g.psi_deg > 0.0 && g.psi_deg < 180.0))
        refuse(std::format("Eq. 6-68: psi = 90 - theta_u/2 - gamma = {} deg must lie in (0, 180) deg (beta1 = 90 - psi "
                           "must be < 90 deg, Eq. 6-67)",
                           g.psi_deg));
    if (g.suction_camber_deg == 0.0) {  // a flat suction surface: R_u infinite, phi = 0
        g.radius_over_pitch = std::numeric_limits<double>::infinity();
        g.expansion_deg = 0.0;
        return g;
    }
    // c/(2 R_u) = sin(theta_u/2) (Eq. 4-7; 6-61 is printed 2R_u/c = sin(theta_u/2)): R_u/s = sigma/(2 sin(theta_u/2))
    g.radius_over_pitch = solidity / (2.0 * std::sin(rad(g.suction_camber_deg) / 2.0));
    const double p = rad(g.psi_deg);
    g.expansion_deg = deg(std::atan2(std::cos(p), std::sin(p) + g.radius_over_pitch));  // [6-69] (divided by s)
    return g;
}

double prandtl_meyer_angle(double mach, double gamma)
{
    require_gamma(gamma);
    if (!(std::isfinite(mach) && mach >= 1.0)) refuse(std::format("Prandtl-Meyer angle: M = {} must be >= 1", mach));
    const double g = (gamma + 1.0) / (gamma - 1.0);
    const double m2 = mach * mach - 1.0;
    return std::sqrt(g) * std::atan(std::sqrt(m2 / g)) - std::atan(std::sqrt(m2));  // closed form of Eq. 6-70
}

double prandtl_meyer_mach(double nu, double gamma)
{
    require_gamma(gamma);
    const double nu_max = 0.5 * kPi * (std::sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0);
    if (!(std::isfinite(nu) && nu >= 0.0 && nu < nu_max))
        refuse(std::format("Prandtl-Meyer: nu = {} rad must lie in [0, nu_max = {}) rad", nu, nu_max));
    // ours: bisection on M (nu(M) increases strictly) to adjacent doubles
    double lo = 1.0, hi = 2.0;
    while (prandtl_meyer_angle(hi, gamma) < nu) {
        lo = hi;
        hi *= 2.0;
        if (!(hi < 1e12)) refuse(std::format("Prandtl-Meyer: no Mach number below 1e12 has nu = {} rad", nu));
    }
    for (int k = 0; k < 200; ++k) {
        const double mid = 0.5 * (lo + hi);
        if (!(mid > lo && mid < hi)) break;
        (prandtl_meyer_angle(mid, gamma) < nu ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double normal_shock_total_pressure_ratio(double mach, double gamma)
{
    require_gamma(gamma);
    if (!(std::isfinite(mach) && mach >= 1.0)) refuse(std::format("normal shock: M = {} must be >= 1", mach));
    // Eqs. 6-72 .. 6-75 solved in closed form for the perfect gas (the Rankine-Hugoniot relations)
    const double m2 = mach * mach;
    const double a = std::pow((gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0), gamma / (gamma - 1.0));
    const double b = std::pow((gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0)), 1.0 / (gamma - 1.0));
    return a * b;
}

ShockLoss shock_loss(double mach1, double camber_deg, double tb_c, double stagger_deg, double solidity, double gamma)
{
    require_positive(mach1, "M'_1", "Sec. 6.7");
    require_gamma(gamma);
    ShockLoss s;
    s.geometry = shock_expansion_geometry(camber_deg, tb_c, stagger_deg, solidity);
    if (mach1 <= 1.0) {  // no bow shock
        s.surface_mach = s.shock_mach = mach1;
        return s;
    }
    s.surface_mach = prandtl_meyer_mach(prandtl_meyer_angle(mach1, gamma) + rad(s.geometry.expansion_deg), gamma);  // [6-70]
    s.shock_mach = std::sqrt(mach1 * s.surface_mach);                                                             // [6-71]
    s.total_pressure_ratio = normal_shock_total_pressure_ratio(s.shock_mach, gamma);
    // omega_s = (P't1 - P't2)/(P't1 - P1) (the loss coefficient of Eq. 7-49), P1/P't1 isentropic at M'_1
    const double p_pt = std::pow(1.0 + 0.5 * (gamma - 1.0) * mach1 * mach1, -gamma / (gamma - 1.0));
    s.loss = (1.0 - s.total_pressure_ratio) / (1.0 - p_pt);
    return s;
}

// ================================================================================================
// Sec. 6.9, 6.10: tip clearance and shroud seal
// ================================================================================================

TipClearanceLoss tip_clearance_loss(const TipClearanceInput& in)
{
    const char* eq = "Eqs. 6-85 .. 6-89";
    require_positive(in.clearance, "delta_c", eq);
    require_positive(in.tip_radius, "r_tip", eq);
    if (in.blade_count < 1) refuse(std::format("{}: the blade count Z = {} must be >= 1", eq, in.blade_count));
    require_positive(in.chord, "c", eq);
    if (!(std::isfinite(in.stagger_rad) && std::abs(in.stagger_rad) < 0.5 * kPi))
        refuse(std::format("{}: the stagger {} rad must lie in (-pi/2, pi/2)", eq, in.stagger_rad));
    if (in.row_number < 1) refuse(std::format("Eq. 6-87: the row number N_row = {} must be >= 1", in.row_number));
    require_positive(in.r1, "r1", eq);
    require_positive(in.r2, "r2", eq);
    require_positive(in.rho1, "rho1", eq);
    require_positive(in.rho2, "rho2", eq);
    require_finite(in.cm1, "C_m1", eq);
    require_finite(in.cm2, "C_m2", eq);
    require_finite(in.ctheta1, "C_theta1", eq);
    require_finite(in.ctheta2, "C_theta2", eq);
    require_positive(in.mean_density, "rho_bar", eq);
    require_positive(in.mass_flow, "m_dot", eq);
    TipClearanceLoss r;
    r.torque = kPi * in.clearance * (in.r1 * in.rho1 * in.cm1 + in.r2 * in.rho2 * in.cm2) *
               (in.r2 * in.ctheta2 - in.r1 * in.ctheta1);  // [6-85]
    const double axial_chord = in.chord * std::cos(in.stagger_rad);
    // [6-86] with |tau|: the pressure difference across the blade is a magnitude (ours: a stator,
    // r2 C_theta2 < r1 C_theta1, would otherwise give Delta P < 0 and an imaginary U_c in Eq. 6-87)
    r.pressure_difference = std::abs(r.torque) / (in.blade_count * in.tip_radius * in.clearance * axial_chord);
    r.leakage_velocity = 0.816 * std::sqrt(2.0 * r.pressure_difference / in.mean_density) /
                         std::pow(static_cast<double>(in.row_number), 0.2);                               // [6-87]
    r.leakage_mass_flow = in.mean_density * r.leakage_velocity * in.blade_count * in.clearance * axial_chord;  // [6-88]
    r.total_pressure_loss = r.pressure_difference * r.leakage_mass_flow / in.mass_flow;                        // [6-89]
    return r;
}

double seal_contraction_ratio(double clearance_over_fin_thickness)
{
    if (!(std::isfinite(clearance_over_fin_thickness) && clearance_over_fin_thickness >= 0.0))
        refuse(std::format("Eq. 6-92: delta_c/t = {} must be >= 0", clearance_over_fin_thickness));
    return 1.0 - 1.0 / (3.0 + std::pow(54.3 / (1.0 + 100.0 * clearance_over_fin_thickness), 3.45));  // [6-92]
}

double seal_throttling_coefficient(int fins, double pressure_ratio)
{
    if (fins < 1) refuse(std::format("Eq. 6-93: the number of fins N = {} must be >= 1", fins));
    if (!(std::isfinite(pressure_ratio) && pressure_ratio >= 0.0 && pressure_ratio < 1.0))
        refuse(std::format("Eq. 6-93: the seal pressure ratio P_R = {} (lower / higher) must lie in [0, 1)", pressure_ratio));
    const double n = static_cast<double>(fins);
    // Aungier (2000) Eq. 4-35, p. 64: 2.143 [ln(N) - 1.464] / (N - 4.322); the 2003 Eq. 6-93 is
    // printed 2.143 ln(N) - 1.464 (negative for N = 2..4, a pole at 4.322).  N is an integer, so
    // N - 4.322 != 0.
    return 2.143 * (std::log(n) - 1.464) / (n - 4.322) * std::pow(1.0 - pressure_ratio, 0.375 * pressure_ratio);
}

double seal_carryover_coefficient(int fins, double clearance_over_fin_pitch)
{
    if (fins < 1) refuse(std::format("Eq. 6-94: the number of fins N = {} must be >= 1", fins));
    if (!(std::isfinite(clearance_over_fin_pitch) && clearance_over_fin_pitch >= 0.0))
        refuse(std::format("Eq. 6-94: delta_c/p = {} must be >= 0", clearance_over_fin_pitch));
    const double n = static_cast<double>(fins);
    double x1 = 0.0, x2 = 0.0;
    if (fins <= 12) {
        x1 = 15.1 - 0.05255 * std::exp(0.507 * (12.0 - n));  // [6-95] = 2000 [4-36]
        x2 = 1.058 + 0.0218 * n;                             // [6-96] = 2000 [4-37]
    } else {
        x1 = 13.15 + 0.1625 * n;  // Aungier (2000) Eq. 4-36, N > 12 (2003 p. 149 refers there)
        x2 = 1.32;                // Aungier (2000) Eq. 4-37, N > 12
    }
    if (clearance_over_fin_pitch > x2 - 1.0)
        refuse(std::format("Eq. 6-94: delta_c/p = {} exceeds X2 - 1 = {} (N = {}), beyond the range of Egli's model "
                           "(p. 149; 2000 Eq. 4-38)",
                           clearance_over_fin_pitch, x2 - 1.0, fins));
    const double r = clearance_over_fin_pitch;
    return 1.0 + x1 * (r - x2 * std::log(1.0 + r)) / (1.0 - x2);  // [6-94]
}

ShroudSealLoss shroud_seal_loss(const ShroudSealInput& in)
{
    const char* eq = "Eqs. 6-90, 6-91";
    require_positive(in.seal_radius, "r_seal", eq);
    require_positive(in.clearance, "delta_c", eq);
    require_positive(in.fin_thickness, "t", eq);
    require_positive(in.fin_pitch, "p", eq);
    require_positive(in.upstream_density, "rho", eq);
    require_positive(in.upstream_temperature, "T", eq);
    require_positive(in.gas_constant, "R", eq);
    if (!(std::isfinite(in.pressure_difference) && in.pressure_difference >= 0.0))
        refuse(std::format("{}: Delta P = {} Pa must be >= 0", eq, in.pressure_difference));
    require_positive(in.mass_flow, "m_dot", eq);
    ShroudSealLoss r;
    r.contraction_ratio = seal_contraction_ratio(in.clearance / in.fin_thickness);
    r.throttling_coefficient = seal_throttling_coefficient(in.fins, in.pressure_ratio);
    r.carryover_coefficient = seal_carryover_coefficient(in.fins, in.clearance / in.fin_pitch);
    // rho and T on the upstream (higher-pressure) side: Aungier (2000) p. 63
    r.leakage_mass_flow = 2.0 * kPi * in.seal_radius * in.clearance * r.throttling_coefficient *
                          r.carryover_coefficient * r.contraction_ratio * in.upstream_density *
                          std::sqrt(in.gas_constant * in.upstream_temperature);  // [6-91]
    r.total_pressure_loss = in.pressure_difference * r.leakage_mass_flow / in.mass_flow;  // [6-90]
    return r;
}

// ================================================================================================
// The row-level leakage call (contract Sec. 4.1 B) and the p. 147 distribution
// ================================================================================================

RowLeakage endwall_leakage(const RowWallData& row, const WallState& qn1, const WallState& qn2, double mass_flow,
                           const Gas& gas)
{
    if (row.row_sign != 1 && row.row_sign != -1)
        refuse(std::format("RowWallData: row_sign = {} must be +1 or -1", row.row_sign));
    RowLeakage r;
    if (!row.shrouded) {
        TipClearanceInput t;
        t.clearance = row.clearance;
        t.tip_radius = row.gap_radius;
        t.blade_count = row.blade_count;
        t.chord = row.chord;
        t.stagger_rad = row.stagger;
        t.row_number = row.row_number;
        t.r1 = qn1.radius;
        t.r2 = qn2.radius;
        t.rho1 = qn1.density;
        t.rho2 = qn2.density;
        t.cm1 = qn1.cm;
        t.cm2 = qn2.cm;
        const double sd = static_cast<double>(row.row_sign);  // cascade signs (contract Sec. 2.3)
        t.ctheta1 = sd * qn1.ctheta;
        t.ctheta2 = sd * qn2.ctheta;
        t.mean_density = 0.5 * (qn1.density + qn2.density);  // rho_bar of 6-87, 6-88: ours (undefined in the book)
        t.mass_flow = mass_flow;
        r.clearance = tip_clearance_loss(t);
        r.total_pressure_loss = r.clearance.total_pressure_loss;
        r.leakage_mass_flow = r.clearance.leakage_mass_flow;
        return r;
    }
    require_gamma(gas.gamma);
    require_positive(qn1.pressure, "P1", "Eq. 6-90");
    require_positive(qn2.pressure, "P2", "Eq. 6-90");
    const WallState& high = qn2.pressure > qn1.pressure ? qn2 : qn1;
    ShroudSealInput si;
    si.seal_radius = row.gap_radius;
    si.clearance = row.clearance;
    si.fins = row.fins;
    si.fin_thickness = row.fin_thickness;
    si.fin_pitch = row.fin_pitch;
    si.pressure_ratio = std::min(qn1.pressure, qn2.pressure) / std::max(qn1.pressure, qn2.pressure);
    si.upstream_density = high.density;  // Aungier (2000) p. 63: the higher-pressure side of the seal
    si.upstream_temperature = high.temperature;
    si.gas_constant = gas.gas_constant;
    si.pressure_difference = std::abs(qn2.pressure - qn1.pressure);  // "across the blade row" (6-90)
    si.mass_flow = mass_flow;
    r.seal = shroud_seal_loss(si);
    r.total_pressure_loss = r.seal.total_pressure_loss;
    r.leakage_mass_flow = r.seal.leakage_mass_flow;
    r.seal_flow_direction = qn2.pressure > qn1.pressure ? 1 : -1;
    return r;
}

double leakage_loss_at(double row_total_pressure_loss, double mass_fraction, EndWall gap_side)
{
    require_finite(row_total_pressure_loss, "Delta P_t", "p. 147 distribution");
    if (!(std::isfinite(mass_fraction) && mass_fraction >= 0.0 && mass_fraction <= 1.0))
        refuse(std::format("p. 147 distribution: the mass fraction {} must lie in [0, 1]", mass_fraction));
    const double from_opposite = gap_side == EndWall::Shroud ? mass_fraction : 1.0 - mass_fraction;
    return 2.0 * row_total_pressure_loss * from_opposite;  // linear, 0 at the opposite wall, mean = Delta P_t
}

double loss_coefficient_from_pressure(double total_pressure_loss, double relative_total_pressure1, double static_pressure1)
{
    require_finite(total_pressure_loss, "Delta P't", "Eq. 7-49");
    require_finite(relative_total_pressure1, "P't1", "Eq. 7-49");
    require_finite(static_pressure1, "P1", "Eq. 7-49");
    const double q = relative_total_pressure1 - static_pressure1;
    if (!(q > 0.0)) refuse(std::format("Eq. 7-49: P't1 - P1 = {} Pa must be > 0", q));
    return total_pressure_loss / q;  // [7-49]
}

// ================================================================================================
// evaluate: Sec. 7.3's use of the models on one stream surface
// ================================================================================================

namespace {

void note_range(CascadePerformance& r, bool outside, const std::string& what)
{
    if (!outside) return;
    r.extrapolated = true;
    r.notes.push_back(what);
}

}  // namespace

CascadePerformance evaluate(const CascadeGeometry& geo, const InletState& in, const DischargeEstimate& out,
                            const Gas& gas, const CorrelationOptions& opt)
{
    // ---- refusals of the input
    if (geo.row_sign != 1 && geo.row_sign != -1)
        refuse(std::format("CascadeGeometry: row_sign = {} must be +1 or -1", geo.row_sign));
    require_positive(geo.chord, "chord c", "CascadeGeometry");
    require_positive(geo.solidity, "solidity sigma", "CascadeGeometry");
    require_finite(geo.stagger, "stagger", "CascadeGeometry");
    require_finite(geo.camber, "camber", "CascadeGeometry");
    require_thickness(geo.thickness_ratio, "CascadeGeometry");
    require_positive(geo.throat_over_pitch, "throat o/s", "CascadeGeometry (Eq. 6-53)");
    require_positive(in.radius, "r1", "InletState");
    require_positive(in.wm, "W_m1", "InletState");
    require_finite(in.wtheta, "W_theta1", "InletState");
    require_finite(in.omega, "omega", "InletState");
    require_positive(in.mach, "M'_1", "InletState");
    require_positive(out.radius, "r2", "DischargeEstimate");
    require_positive(out.wm, "W_m2", "DischargeEstimate");
    require_gamma(gas.gamma);

    CascadePerformance r;

    // ---- the cascade frame: s_row (contract Sec. 0.2), then Aungier's p. 128 flip for a guide vane
    // of negative camber; a non-IGV row of negative camber after s_row is refused
    const bool igv = geo.kind == BladeRowKind::InletGuideVane;
    int s = geo.row_sign;
    if (static_cast<double>(s) * geo.camber < 0.0) {
        if (!igv)
            refuse(std::format("Sec. 6.3, p. 128: the correlations hold for positive camber; camber {} deg with row_sign "
                               "{} is negative in the cascade frame (only an inlet guide vane is flipped)",
                               deg(geo.camber), geo.row_sign));
        s = -s;
        r.notes.push_back("inlet guide vane of negative camber: all angles sign-changed (p. 128)");
    }
    r.frame_sign = s;
    const double sd = static_cast<double>(s);
    const double gamma_d = sd * deg(geo.stagger);
    const double theta = sd * deg(geo.camber);
    const double beta1 = sd * deg(std::atan2(in.wtheta, in.wm));  // beta'_1 = atan(W_theta1 / W_m1)
    const double wm_ratio = out.wm / in.wm, radius_ratio = out.radius / in.radius;
    const double u_ratio = sd * in.omega * in.radius / in.wm;  // U1 / W_m1: a tangential velocity, sign-changed too
    const double sigma = geo.solidity, tb = geo.thickness_ratio, ac = geo.max_camber_position;
    require_angle(beta1, "beta1", "InletState");
    // Re_c = rho1 W1 c / mu1 (ours, contract F12)
    if (in.density > 0.0 && in.viscosity > 0.0)
        r.reynolds = in.density * std::hypot(in.wm, in.wtheta) * geo.chord / in.viscosity;

    // ---- blade angles (Ch. 4)
    const BladeAngles ba = blade_angles(geo.camberline, gamma_d, theta, ac);
    r.chi1_deg = ba.chi1_deg;
    r.chi2_deg = ba.chi2_deg;
    r.kappa1_deg = ba.kappa1_deg;
    r.kappa2_deg = ba.kappa2_deg;
    const double ksh = profile_shape_factor(geo.profile);
    const double kti = incidence_thickness_factor(tb);
    const double ktd = deviation_thickness_factor(tb);

    // ---- data ranges of the figures (notes, not refusals)
    note_range(r, sigma < 0.4 || sigma > 2.0,
               std::format("solidity {} outside 0.4 .. 2.0, the data of Figs. 6-4, 6-5 (Eqs. 6-13, 6-15)", sigma));
    note_range(r, tb > 0.12, std::format("t_b/c = {} above 0.12, the data of Figs. 6-3, 6-9 (Eqs. 6-10, 6-25)", tb));

    // ---- design incidence (Sec. 6.2)
    if (opt.incidence == DesignIncidenceModel::AngleOfAttack) {
        r.design_attack_deg = design_attack_angle_deg(ksh, kti, theta, ac, sigma);  // [6-8]
        r.design_incidence_deg = r.design_attack_deg + gamma_d - ba.kappa1_deg;     // [6-18]
    } else {
        const LiebleinIncidence li = lieblein_design_incidence(ksh, kti, gamma_d, theta, ac, ba.kappa1_deg, sigma);
        r.design_attack_deg = li.attack_deg;
        r.design_incidence_deg = li.incidence_deg;
    }
    r.design_beta1_deg = ba.kappa1_deg + r.design_incidence_deg;  // [6-2] at i*
    require_angle(r.design_beta1_deg, "beta1* = kappa1 + i*", "Sec. 6.2");
    note_range(r, r.design_beta1_deg < 0.0 || r.design_beta1_deg > 70.0,
               std::format("design inlet angle beta1* = {} deg outside 0 .. 70 deg, the data of Figs. 6-4 .. 6-8",
                           r.design_beta1_deg));

    // ---- design deviation (Sec. 6.3): Howell's Eq. 6-26 for guide vanes and the parabolic / A4K6 lines
    const bool howell = igv || geo.camberline == Camberline::ParabolicArc || geo.camberline == Camberline::NacaA4K6;
    r.design_deviation_deg =
        howell ? howell_design_deviation_deg(ksh, ktd, ac, ba.kappa2_deg, theta, r.design_beta1_deg, sigma)  // [6-26]
               : lieblein_design_deviation_deg(geo.camberline, ksh, ktd, r.design_beta1_deg, theta, sigma);  // [6-19]
    r.design_beta2_deg = ba.kappa2_deg + r.design_deviation_deg;  // [6-3] at i*
    require_angle(r.design_beta2_deg, "beta2* = kappa2 + delta*", "Sec. 6.3");

    // ---- design diffusion and loss (Sec. 6.4).  Ours: the design state uses the current W_m2/W_m1,
    // r2/r1 and U1/W_m1, and W_max of Klapproth's Eq. 6-39 at i = i* (p. 133: "for both design and
    // off-design operation"); for a planar stationary cascade with W_m2 = W_m1 it is Eq. 6-35.
    const double b1s = r.design_beta1_deg, b2s = r.design_beta2_deg;
    r.wmax_ratio_design = wmax_ratio(b1s, b2s, sigma, wm_ratio, radius_ratio, u_ratio, 0.0, 0.0);
    r.d_eq_design = equivalent_diffusion(r.wmax_ratio_design, b1s, b2s, wm_ratio);  // [6-36]
    r.d_factor_design = diffusion_factor(b1s, b2s, sigma, wm_ratio);               // [6-31] at i*
    const double w1_w2_design = 1.0 / velocity_ratio_w2_w1(b1s, b2s, wm_ratio);
    double k1 = 0.0, k2 = 0.0;
    switch (opt.calibration) {
    case LossCalibration::TwoDimensionalCascade:
        k1 = opt.loss == ProfileLossModel::DFactor ? 0.0035 : 0.004;  // Eq. 6-34 / Eq. 6-37
        k2 = 1.0;
        break;
    case LossCalibration::Constants:
        require_positive(opt.k1, "K1", "CorrelationOptions (Eq. 6-45/6-46)");
        require_finite(opt.k2, "K2", "CorrelationOptions (Eq. 6-45/6-46)");
        k1 = opt.k1;
        k2 = opt.k2;
        break;
    case LossCalibration::AspectRatio: {
        if (opt.loss != ProfileLossModel::EquivalentDiffusion)
            refuse("Sec. 6.11: the calibration K1 = 0.0073 with K2 of Eq. 6-99 belongs to Eq. 6-46 (p. 150), not to the "
                   "D-factor Eq. 6-45");
        require_positive(geo.pitch_over_height, "pitch over height s/h", "CascadeGeometry (Eq. 6-99)");
        k1 = 0.0073;  // p. 150
        const double k_re = opt.reynolds_correction ? reynolds_factor(r.reynolds) : 0.0;  // [6-101], [6-102]
        k2 = aspect_ratio_k2(geo.pitch_over_height, b2s, k_re, k1);                       // [6-99], [6-100]
        break;
    }
    }
    if (opt.reynolds_correction && opt.calibration != LossCalibration::AspectRatio)
        refuse("Eq. 6-100: the Reynolds correction is part of K2 of the Sec. 6.11 calibration (LossCalibration::"
               "AspectRatio)");
    r.losses.k1 = k1;
    r.losses.k2 = k2;
    r.losses.design = opt.loss == ProfileLossModel::DFactor
                          ? design_loss_d_factor(r.d_factor_design, b2s, sigma, k1, k2)                       // [6-45]
                          : design_loss_equivalent_diffusion(r.d_eq_design, b2s, sigma, w1_w2_design, k1, k2);  // [6-46]
    if (r.d_eq_design > 2.0) {
        r.design_diffusion_limit_exceeded = true;
        r.notes.push_back(std::format("D_eq* = {} exceeds Lieblein's diffusion limit 2.0 at design incidence (p. 132)",
                                      r.d_eq_design));
    }
    if (r.d_factor_design > 0.6) {  // a refusal with the D-factor model (design_loss_d_factor), a flag here
        r.design_diffusion_limit_exceeded = true;
        r.notes.push_back(std::format("D* = {} exceeds the diffusion limit 0.6 at design incidence (p. 130)",
                                      r.d_factor_design));
    }

    // ---- stall incidence, low speed (Sec. 6.5)
    const StallRanges sr = stall_ranges(r.design_attack_deg, gamma_d, theta);
    r.stall_range_negative_deg = sr.negative_deg;
    r.stall_range_positive_deg = sr.positive_deg;
    r.beta1_negative_stall_deg = sr.beta1c_deg;
    r.beta1_positive_stall_deg = sr.beta1s_deg;
    if (sr.beta1c_limited)
        r.notes.push_back(std::format("Eq. 6-47 evaluated at beta_1c = 20 deg (solved {} deg; Aungier's limit, p. 135)",
                                      sr.beta1c_deg));
    if (!(sr.negative_deg > 0.0 && sr.positive_deg > 0.0))
        refuse(std::format("Eqs. 6-49, 6-50: the low-speed stall ranges R_c = {}, R_s = {} deg must both be > 0",
                           sr.negative_deg, sr.positive_deg));

    // ---- Mach number effects (Sec. 6.6)
    const double M = in.mach;
    double ic = mach_negative_stall_incidence_deg(r.design_incidence_deg, sr.negative_deg, M);       // [6-51]
    const double is = mach_positive_stall_incidence_deg(r.design_incidence_deg, sr.positive_deg, M, ksh);  // [6-52]
    if (ksh > 1.0) r.notes.push_back(std::format("Eq. 6-52 evaluated with K_sh = 1 (K_sh = {}; p. 136)", ksh));
    r.choke_beta1_deg = choke_inlet_angle_deg(geo.throat_over_pitch, M, gas.gamma);  // [6-53]
    if (std::isfinite(r.choke_beta1_deg)) {
        const double limit = r.choke_beta1_deg - ba.kappa1_deg + 1.0;  // [6-54]
        if (ic < limit) {
            r.notes.push_back(std::format("choke limit Eq. 6-54 active: i_c raised from {} to {} deg", ic, limit));
            r.choke_limit_active = true;
            ic = limit;
        }
    }
    r.negative_stall_incidence_deg = ic;
    r.positive_stall_incidence_deg = is;
    if (!(is > ic))
        refuse(std::format("Sec. 6.6: no low-loss incidence range is left, i_c = {} >= i_s = {} deg (M'_1 = {}, choke "
                           "beta1 = {} deg)",
                           ic, is, M, r.choke_beta1_deg));
    r.minimum_loss_incidence_deg = minimum_loss_incidence_deg(ic, is, sr.negative_deg, sr.positive_deg);  // [6-55]
    const double w_star_w1 = sonic_velocity_ratio(M, gas.gamma);
    r.critical_mach = critical_mach_number(M, w_star_w1, r.wmax_ratio_design);  // [6-59]
    r.losses.minimum = minimum_loss(r.losses.design, r.minimum_loss_incidence_deg, r.design_incidence_deg,
                                    sr.positive_deg, ksh, M, r.critical_mach, w_star_w1);  // [6-58], [6-60]
    if (M > r.critical_mach && M > 1.0) r.notes.push_back("Eq. 6-60 evaluated with M' = 1 (M'_1 > 1; p. 138)");

    // ---- shock loss (Sec. 6.7)
    r.losses.shock = shock_loss(M, theta, tb, gamma_d, sigma, gas.gamma).loss;

    // ---- the operating point (Sec. 6.8)
    r.beta1_deg = beta1;
    r.incidence_deg = beta1 - ba.kappa1_deg;  // [6-2]
    r.xi = normalized_incidence(r.incidence_deg, r.minimum_loss_incidence_deg, ic, is);  // [6-78], [6-79]
    const double factor = off_design_loss_factor(r.xi);                                  // [6-80] .. [6-82]
    r.losses.incidence = r.losses.minimum * (factor - 1.0);
    r.losses.total = r.losses.shock + r.losses.minimum * factor;
    r.loss = r.losses.total;
    r.positive_stall_margin_deg = is - r.incidence_deg;
    r.negative_stall_margin_deg = r.incidence_deg - ic;

    // ---- off-design deviation: the slope at the design inlet angle beta1* (ours: Eq. 6-76's beta1 is
    // read as beta1*, the slope being the one "at the design incidence angle", p. 142)
    r.deviation_slope = igv ? 0.0 : deviation_incidence_slope(r.design_beta1_deg, sigma);  // [6-76]
    r.deviation_deg = off_design_deviation_deg(r.design_deviation_deg, r.deviation_slope, r.incidence_deg,
                                               r.design_incidence_deg, wm_ratio, igv);  // [6-77]
    r.beta2_deg = ba.kappa2_deg + r.deviation_deg;                                         // [6-3]
    require_angle(r.beta2_deg, "beta2 = kappa2 + delta", "Eq. 6-77");
    r.beta2 = sd * rad(r.beta2_deg);

    // ---- off-design diffusion and the wake blockage (p. 145)
    const double excess = r.incidence_deg - r.design_incidence_deg;
    double alpha = 0.0;
    if (excess > 0.0) alpha = opt.lieblein_alpha ? *opt.lieblein_alpha : lieblein_alpha(geo.profile);
    else if (excess < 0.0)
        r.notes.push_back("i < i*: the term alpha (i - i*)^1.43 of Eq. 6-39 (stated for i >= i*, p. 132) is omitted "
                          "(ours)");
    r.wmax_ratio = wmax_ratio(beta1, r.beta2_deg, sigma, wm_ratio, radius_ratio, u_ratio, std::max(excess, 0.0),
                              alpha);                                                    // [6-39]
    r.d_eq = equivalent_diffusion(r.wmax_ratio, beta1, r.beta2_deg, wm_ratio);        // [6-30]
    r.d_factor = diffusion_factor(beta1, r.beta2_deg, sigma, wm_ratio);               // [6-31]
    r.wake_blockage = wake_blockage(r.d_eq);                                             // [6-83], [6-84]
    r.k_w = 1.0 - r.wake_blockage;                                                       // [7-6]
    return r;
}

}  // namespace BladeGen::Turbolenza
