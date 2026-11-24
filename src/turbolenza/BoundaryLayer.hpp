/// @file turbolenza/BoundaryLayer.hpp
/// @brief The blade-surface boundary layer of Aungier (2003) Sec. 5.6: an integral method marched
///        along one blade surface from given edge conditions, its wake, and the profile loss.
///
/// Not a v1 port: TURBOlenza v1 has no boundary layer.  The method is Aungier's, read from the
/// page images of R. H. Aungier, Axial-Flow Compressors (ASME Press, 2003), Sec. 5.6, printed
/// pp. 107-113 (Eqs. 5-106 to 5-141), and of R. H. Aungier, Centrifugal Compressors (ASME Press,
/// 2000), Ch. 15, printed pp. 273-279 (Eqs. 15-1 to 15-35, the per-station procedures on p. 276
/// laminar and pp. 278-279 turbulent).  Planar cascade: stream-sheet thickness b = 1.
///
///  * momentum integral [5-106 = 15-1]:  d(rho_e W_e^2 theta)/dx + delta* rho_e W_e dW_e/dx = tau_w,
///    c_f = tau_w / (rho_e W_e^2 / 2) [5-109]; x the arc length along the surface;
///  * laminar, Gruschwitz (1950) [5-110 .. 5-125 = 15-4 .. 15-22]:
///      theta/delta' = 37/315 - Lambda/945 - Lambda^2/9072                          [5-115 = 15-14]
///      delta_E/delta' = (798048 - 4656 Lambda - 758 Lambda^2 - 7 Lambda^3) / 4324320 [5-116 = 15-15]
///      delta_W/delta' = 3/10 - Lambda/120 + F W_e^2 / (2 cp T_e)                    [5-117 = 15-16]
///      F = 0.232912 - 0.831483 (L/100) + 0.650584 (L/100)^2 + 17.8063 (L/100)^3      [5-118 = 15-17]
///        (as printed in both books and kept; the Pr = 1 Crocco-Busemann integral int (1-f)(1-f^2) deta
///        of the profile gives the same first three coefficients and 1.61875 for the cubic one, a
///        difference of up to 0.028 in F at |Lambda| = 12: test_boundary_layer.cpp)
///      delta_h = W_e^2 delta_E / (2 cp T_e)  (adiabatic wall, Pr = 1)               [5-119 = 15-18]
///      delta* = delta_h + delta_W                                                   [5-120 = 15-19]
///      b0 = rho_e / rho_w = T't / T_e                                               [5-121 = 15-21]
///      K = Lambda (theta/delta')^2 = b0 rho_e theta^2 / mu  dW_e/dx                  [5-122 = 15-10]
///      K = (37/315 - Lambda/945 - Lambda^2/9072)^2 Lambda                            [5-124 = 15-22]
///      -12 <= Lambda <= 12, Lambda = -12 laminar separation                         [5-125]
///      1/2 c_f = mu / (rho_e W_e delta') (2 + Lambda/6)  Centrifugal 2000 Eq. 15-20 (the profile's
///              wall slope C1 = 2 + Lambda/6 of Eq. 5-114), or
///      1/2 c_f = mu / (rho_e W_e delta') (1 + Lambda/6)  Axial 2003 Eq. 5-123 as printed
///      (delta - delta*) = delta' (7/10 + Lambda/120)                               [5-138 = 15-34]
///  * transition: Re_theta = rho_e W_e theta / mu > 250 [5-126 = 15-23]; theta and (delta - delta*)
///    carried over [5-136, 5-137 = 15-32, 15-33];
///  * turbulent, Head (1958, 1968) / Green (1968) [5-127 .. 5-135 = 15-24 .. 15-31]:
///      d(rho_e W_e (delta - delta*))/dx = rho_e W_e E                               [5-127 = 15-24]
///      H1 = (delta - delta*)/theta, H = delta*/theta                                [5-128, 5-129]
///      H_k = 1 + [0.9 / (H1 - 3.3)]^0.75                                            [5-132 = 15-28]
///      E = 0.025 (H_k - 1)                                                          [5-133 = 15-27]
///      H = (H_k + 1) T't / T_e - 1                                                  [5-131 = 15-30]
///      c_f,inc = 0.246 exp(-1.561 H_k) Re_theta^-0.268  (Ludwieg and Tillmann 1950) [5-134 = 15-29]
///      c_f = c_f,inc H_k (H_k + 1) / [2 H_k + H (H_k - 1)]                         [5-135 = 15-31]
///      H_k >= 2.4 separation; H_k limited to 2.4 to continue (2003 p. 111, 2000 p. 279);
///  * streamwise derivatives: the backward difference of Eq. 5-94, (u_i - u_{i-1}) / dx (2003 p. 113),
///    the station equations solved iteratively before the next station (2003 p. 113);
///  * loss (Lieblein and Roudebush 1956) [5-139 .. 5-141]:
///      omega = (cos beta_in / cos beta_out)^2 (2 Theta + Delta*^2) / (1 - Delta*)^2,
///      Theta = sum theta / (S cos beta_out), Delta* = sum delta* / (S cos beta_out), the sums over
///      both surfaces at the trailing edge, S the pitch.
///
/// Our choices (Aungier does not give them; each documented where it is made):
///  * Lambda from K: K(Lambda) is strictly increasing on [-12, 12) with K'(12) = 0 exactly
///    (K' = g (37/315 - 3 Lambda/945 - 5 Lambda^2/9072), whose root is Lambda = 12), so the root in
///    [-12, 12] is unique; it is found by bisection to adjacent doubles; K outside [K(-12), K(12)]
///    gives Lambda = -12 or 12 (Eq. 5-125);
///  * the laminar station equation theta = F(theta) (p. 276 steps 2-7) is solved by bracketing and
///    bisection instead of plain substitution, which does not converge at the first station
///    (|F'| -> 1 there, since the upstream momentum flux is 0); the bracket starts from the "safe"
///    guess of Exercise 5.6: at the station after the first, the flat-plate value of the method
///    theta = sqrt(2 g(0) C1(0) mu dx / (rho_e W_e)); further on, the upstream momentum flux held
///    (theta_{k-1} rho_{k-1} W_{k-1}^2 / (rho_k W_k^2)), which accounts for the change of the edge
///    conditions only;
///  * the turbulent station (pp. 278-279 steps 1-7) by nested bisection: for a theta the entrainment
///    equation has one root (delta - delta*), then the momentum equation is bracketed from the
///    upstream momentum flux held (plain substitution diverges where dx / theta is large, right
///    after transition on a coarse grid); no root: BoundaryLayerError;
///  * the first station: theta = 0, so delta' = delta* = 0; its K = 0 gives Lambda = 0 and the
///    reported H is the profile's H at Lambda = 0 (the limit of delta*/theta), c_f is NaN;
///  * the transition station: the first station whose laminar Re_theta exceeds the threshold; its
///    reported values are the turbulent ones of the matched theta and (delta - delta*) (as Fig. 5-15's
///    transition dots, which sit on the turbulent level);
///  * the H_k limit acts in the closures (E, H, c_f); the marched theta and (delta - delta*) are not
///    reset; H1 <= 3.3 (H_k unbounded) is H_k = the limit;
///  * laminar separation (Lambda = -12 reached) is reported; the laminar march goes on with Lambda
///    held at -12 (Gruschwitz's limit, Eq. 5-125) until transition;
///  * the wake (not in Aungier): Eqs. 5-106 and 5-127 with tau_w = 0 along the exit periodic line,
///    started from the sums of both surfaces' theta and (delta - delta*) at the last blade station,
///    the turbulent closures and the same H_k limit.
#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// The boundary-layer analysis refused its input or failed (no value is clipped or repaired).
class BoundaryLayerError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// Sutherland's law (F. M. White, Viscous Fluid Flow, McGraw-Hill, Ch. 1, the Sutherland law for
/// the viscosity of gases): mu = mu_ref (T / T_ref)^1.5 (T_ref + S) / (T + S).  Air with the
/// user's constants mu_ref = 1.716e-5 Pa s, T_ref = 273.15 K, S = 110.4 K (the book was not opened
/// for this port: the equation and the constants are not verified against White's page).
/// Aungier does not say at which temperature mu is taken: here the edge static temperature T_e.
struct SutherlandLaw
{
    double mu_ref = 1.716e-5;   ///< Pa s
    double t_ref = 273.15;      ///< K
    double s = 110.4;           ///< K
};
/// mu(T); BoundaryLayerError unless T and the constants are finite and positive.
[[nodiscard]] double sutherland_viscosity(double temperature, const SutherlandLaw& law = {});

enum class LaminarCfForm { Centrifugal2000, Axial2003AsPrinted };  ///< "centrifugal_2000", "axial_2003_as_printed"
[[nodiscard]] const char* laminar_cf_form_name(LaminarCfForm f) noexcept;

/// The state of one station.
enum class BlRegime { Laminar, LaminarSeparated, Transition, Turbulent, Separated };
[[nodiscard]] const char* bl_regime_name(BlRegime r) noexcept;  ///< "laminar", "laminar_separated", ...

inline constexpr double kGruschwitzLambdaLimit = 12.0;   ///< Eq. 5-125

struct BoundaryLayerSettings
{
    LaminarCfForm laminar_cf_form = LaminarCfForm::Centrifugal2000;
    double transition_re_theta = 250.0;  ///< Eq. 5-126
    double separation_hk = 2.4;          ///< 2003 p. 111
};

/// The edge conditions of one surface at its stations, in SI.  x strictly increasing.
struct EdgeConditions
{
    std::vector<double> x;               ///< arc length from the LE point [m]
    std::vector<double> W, rho, T, mu;   ///< edge speed [m/s], density, static temperature [K], viscosity
    double cp = 0.0;                     ///< J/(kg K)
    double total_temperature = 0.0;      ///< T't, the relative total temperature [K]
};

/// The layer at the stations of one surface (or of the wake).
struct BoundaryLayerProfile
{
    std::vector<double> x, W, rho, T, mu;
    std::vector<double> theta, delta_star, delta, H, Hk, cf, Lambda, Re_theta;
    std::vector<double> mass_deficit;    ///< rho_e W_e delta* [kg/(m s)]
    std::vector<BlRegime> regime;
    int transition_index = -1;           ///< the transition station, -1 none
    int laminar_separation_index = -1;   ///< first station with Lambda = -12, -1 none
    int separation_index = -1;           ///< first turbulent station with H_k >= the limit, -1 none
};

// --- the closures (public for the tests) -----------------------------------------------------

/// theta/delta' of the Gruschwitz profile [5-115].
[[nodiscard]] double gruschwitz_theta_ratio(double lambda) noexcept;
/// K(Lambda) [5-124].
[[nodiscard]] double gruschwitz_k(double lambda) noexcept;
/// Lambda of K on [-12, 12] (unique; bisection to adjacent doubles; clamped outside [K(-12), K(12)]).
[[nodiscard]] double gruschwitz_lambda(double k);
/// The wall slope C1 of the profile in the chosen c_f form: 2 + Lambda/6 or 1 + Lambda/6.
[[nodiscard]] double laminar_wall_factor(double lambda, LaminarCfForm form) noexcept;

/// The thicknesses of a laminar station with momentum thickness theta > 0 and shape factor Lambda.
struct LaminarThicknesses
{
    double delta_prime = 0.0, delta_star = 0.0, delta_minus_delta_star = 0.0, H = 0.0;
};
/// [5-115 .. 5-120, 5-138]; mach2 = W_e^2 / (2 cp T_e).  theta = 0 gives zero thicknesses and the
/// profile's H (the limit delta*/theta).
[[nodiscard]] LaminarThicknesses laminar_thicknesses(double lambda, double theta, double mach2) noexcept;

/// H_k of H1 [5-132]; H1 <= 3.3 is +infinity.
[[nodiscard]] double head_hk(double h1) noexcept;
/// H1 of H_k (the inverse of 5-132), H_k > 1.
[[nodiscard]] double head_h1(double hk) noexcept;
/// E [5-133].
[[nodiscard]] double head_entrainment(double hk) noexcept;
/// H of H_k and T't / T_e [5-131].
[[nodiscard]] double green_h(double hk, double total_to_static) noexcept;
/// c_f of H_k, H and Re_theta [5-134, 5-135].
[[nodiscard]] double green_cf(double hk, double h, double re_theta) noexcept;

// --- the march ---------------------------------------------------------------------------------

/// One surface from its first station (theta = 0) to its last (BoundaryLayerError for invalid
/// edge conditions or a station equation without a solution).
[[nodiscard]] BoundaryLayerProfile march_surface(const EdgeConditions& edge, const BoundaryLayerSettings& settings);

/// The wake along the stations of `edge` (station 0 = the last blade station), from theta_0,
/// (delta - delta*)_0 and delta*_0 there (station 0 holds them as given), tau_w = 0 (H_k limited
/// as on the blade).
[[nodiscard]] BoundaryLayerProfile march_wake(const EdgeConditions& edge, double theta0, double delta_minus_delta_star0,
                                              double delta_star0, const BoundaryLayerSettings& settings);

/// Eqs. 5-139 to 5-141.
struct ProfileLoss
{
    double omega = 0.0, Theta = 0.0, Delta_star = 0.0;
};
[[nodiscard]] ProfileLoss profile_loss(double theta_sum, double delta_star_sum, double pitch, double beta_in,
                                       double beta_out);

}  // namespace BladeGen::Turbolenza
