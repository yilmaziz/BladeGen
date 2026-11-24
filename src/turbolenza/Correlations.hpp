/// @file turbolenza/Correlations.hpp
/// @brief The empirical cascade performance models of Aungier (2003) Ch. 6: design incidence and
///        deviation, design loss and diffusion factors, stall incidence, Mach number effects, bow
///        shock loss, off-design incidence / deviation / loss, tip clearance and shroud seal leakage
///        losses and the wake blockage K_W, for a blade row on one stream surface.
///
/// Not a v1 port (TURBOlenza v1 has no empirical models).  Every equation is read from the page
/// images of R. H. Aungier, Axial-Flow Compressors (ASME Press, 2003), Ch. 6 "Empirical Performance
/// Models Based On Two-Dimensional Cascade Tests", printed pp. 117-152, Eqs. 6-1 .. 6-102; the
/// blade angles from Ch. 4, printed pp. 66-76 (Eqs. 4-7, 4-17 .. 4-20, 4-29, 4-32 .. 4-35); the
/// use by the through-flow from Ch. 7 Sec. 7.3, pp. 162-164 (Eqs. 7-43 .. 7-50) and Eq. 7-6 p. 157
/// (K_W = 1 - B_wake).  Two equations are taken from R. H. Aungier, Centrifugal Compressors (ASME
/// Press, 2000), Ch. 4, printed p. 64, because the 2003 print is wrong or refers there:
///   * the seal throttling coefficient: 2000 Eq. 4-35  C_t = 2.143 [ln N - 1.464] / (N - 4.322) ..
///     (2003 Eq. 6-93 is printed "2.143 ln(N) - 1.464", which is negative for N = 2, 3, 4 and has a
///     pole at N = 4.322; the bracketed 2000 form has a removable singularity there,
///     ln 4.322 = 1.4638, and follows Egli's curves, C_t ~ 1/sqrt(N), Fig. 4-8);
///   * the carryover coefficient for N > 12: 2000 Eqs. 4-36, 4-37 (2003 p. 149 refers to them).
/// One 2003 equation is printed wrong and is used in the form of Ch. 4: Eqs. 6-61 and 6-63 read
/// "2R/c = sin(theta/2)"; the circular arc through both blade ends has c/2 = R sin(theta/2)
/// (Eq. 4-7, p. 66), i.e. c/(2R) = sin(theta/2), and that is what is used.
///
/// ## Conventions (aligned with the q3d data contract between Chs. 5, 6, 7, 8 and 12)
///
///  * Input frame: the BOOK's through-flow frame of Ch. 7 (contract Sec. 0.1): angles measured from
///    the meridional direction on the stream surface, signed with +theta; W_theta = C_theta - omega r
///    (Eq. 7-44), tan beta' = W_theta / W_m (Eq. 7-50); omega signed (TURBOlenza / BladeGen map onto
///    it with omega_book = -2 pi rpm / 60, contract Sec. 0.3).  Ch. 6 calls the reference direction
///    "axial" in its text (p. 119) and "meridional" in its nomenclature (p. 118): meridional is used.
///  * The row sign s_row = CascadeGeometry::row_sign (+1 or -1, contract Sec. 0.2): every signed input
///    (stagger, camber, beta'_1, W_theta1, omega) is multiplied by s_row at entry, which gives the
///    positive-camber CASCADE frame the correlations are written in (theta = kappa1 - kappa2 > 0,
///    Eq. 6-6; turning beta1 - beta2 > 0, Eq. 6-4; Fig. 6-1).  A row whose camber is still negative
///    after s_row is refused, unless it is an inlet guide vane: then Aungier's own rule of p. 128 is
///    applied on top (change the signs of all blade and flow angles, apply the correlations, change
///    them back).  CascadePerformance::frame_sign is the total factor used (s_row, or -s_row for such a
///    guide vane): cascade-frame angle = frame_sign * book-frame angle.  The output beta2 is in the
///    book frame (beta'_2 = frame_sign * beta2_cascade); every _deg field is in the cascade frame.
///  * Units: radians in the structs; DEGREES in every function and field named _deg: the correlations
///    are dimensional in degrees (p. 135: "alpha, beta1 and theta are all expressed in degrees";
///    p. 144 for Eq. 6-77; implied by the figures elsewhere; (i - i*)^1.43 of Eq. 6-38 in degrees as
///    in Lieblein's form, not stated on the page).  SI otherwise (m, Pa, kg/m^3, kg/s, K).
///  * Definitions: i = beta1 - kappa1 (6-2), delta = beta2 - kappa2 (6-3), alpha = beta1 - gamma
///    (6-1), sigma = c / s (6-5), theta = kappa1 - kappa2 (6-6), kappa1 = chi1 + gamma,
///    kappa2 = gamma - chi2 (4-32, 4-33).
///  * Where: inlet = the quasi-normal QN1 upstream of the row, discharge = QN2 downstream (Ch. 7
///    places no quasi-normal inside a row, p. 155).  The through-flow iteration (Sec. 7.3, p. 163):
///    the caller gives its current W_m2 and r2 (DischargeEstimate; start W_m2 = W_m1); evaluate()
///    returns beta'_2, the loss coefficient and K_W; the caller updates W_m2 and calls again.  Inside
///    evaluate() nothing iterates between deviation and diffusion: delta (6-77) does not depend on D_eq,
///    and D_eq (6-39 with C_theta2 = W_m2 tan beta'_2 + r2 omega, 7-50) follows from delta.
///  * Loss coefficient: omega_bar = Delta P't / (P't1 - P1), the relative total pressure loss over the
///    relative inlet dynamic head at QN1 on the stream surface (Eq. 7-49, P't2 = P't2,id -
///    omega_bar (P't1 - P1); Ch. 2 Eq. 2-68).  The correlations' omega_bar (Eqs. 6-46, 6-58 .. 6-82) is
///    used as that.  The shock loss (a normal-shock Delta P't) is converted with the same P't1 - P1 of
///    the perfect gas at M'_1 (ours: the book is silent on the conversion).  The row-level leakage
///    losses are pressures [Pa]; leakage_loss_at() spreads them across the annulus and
///    loss_coefficient_from_pressure() converts them (ours, contract Sec. 2.3).
///  * Wake blockage: K_W = 1 - B_wake (Eq. 7-6), for Eq. 7-5 at QN2.
///  * Chord Reynolds number (Eqs. 6-101, 6-102): Re_c = rho1 W1 c / mu1 at QN1 (ours: the book names
///    Re_c, pp. 118, 150, without a formula; contract F12).  Needed only with the Reynolds correction.
///  * Gas: a thermally and calorically perfect gas (Gas::gamma, Gas::gas_constant) where Ch. 6 needs an
///    equation of state: choke (6-53, rho* W* "for the local inlet relative total conditions", p. 136),
///    W* of Eqs. 6-59 / 6-60, the bow shock (6-70 .. 6-75), the seal (6-91).  Aungier uses the general
///    equation of state of Ch. 2 and says the perfect gas gives the closed forms (p. 141); a real-gas
///    version is not implemented.
///
/// ## Refuse, never repair
///
/// Inputs that make a correlation undefined (non-positive solidity, chord, Mach number, W_m1, W_m2
/// or radii; |angles| >= 90 deg; a non-IGV row of negative camber after s_row; a negative thickness; a/c outside Eq. 4-20 for a parabolic arc; the
/// D-factor loss model beyond D* = 0.6, where Aungier says Eq. 6-34 "should not be applied"; a
/// seal with delta_c/p > X2 - 1, outside Egli's model; a Prandtl-Meyer expansion beyond nu_max; no
/// low-loss incidence range left, i_c >= i_s) throw CorrelationError naming the equation and range.
/// Inputs inside the equations' domain but outside the data range of the book's figures (beta1
/// 0..70 deg, Figs. 6-4 .. 6-8; solidity 0.4 .. 2.0, Figs. 6-4, 6-5; t_b/c <= 0.12, Fig. 6-9) are
/// evaluated and reported in CascadePerformance::notes with extrapolated = true.  The limits that
/// Aungier himself imposes are applied and reported as notes: beta_1c >= 20 deg in Eq. 6-47
/// (p. 135), K_sh <= 1 in Eq. 6-52 (p. 136), the choke limit Eq. 6-54 on i_c, M' <= 1 in Eq. 6-60
/// (p. 138).  Each choice that is ours is marked "ours" where it is made.
#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// A correlation refused its input (no value is clipped or repaired).
class CorrelationError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

// ================================================================================================
// The public data contract: geometry, inlet, options, result
// ================================================================================================

/// The blade profile (thickness distribution) family: sets K_sh (p. 121) and Lieblein's alpha of
/// Eq. 6-38 (p. 132).
enum class ProfileFamily {
    Naca,               ///< NACA 65- (and 63-) series thickness: K_sh = 1.0, alpha = 0.0117
    C4,                 ///< British C.4: K_sh = 1.1, alpha = 0.007
    DoubleCircularArc,  ///< DCA: K_sh = 0.7; alpha of Eq. 6-38 not given (CorrelationOptions::lieblein_alpha)
};

/// The camberline: sets the blade angles (Ch. 4) and the deviation model (Sec. 6.3).
enum class Camberline {
    Naca65,        ///< NACA 65-series (A10) mean line as its equivalent circular arc (p. 120); a/c = 0.5
    CircularArc,   ///< circular arc (also DCA); a/c = 0.5
    ParabolicArc,  ///< parabolic arc, 0.25 < a/c < 0.75 (Eq. 4-20)
    NacaA4K6,      ///< NACA A4K6 guide-vane mean line as an equivalent parabolic arc (p. 66)
};

/// The blade row: an inlet guide vane uses Howell's deviation Eq. 6-26 (p. 128) and no incidence
/// term in Eq. 6-77 (p. 144).  Rotor and stator are treated alike (a rotor through InletState::omega,
/// the blade speed of Eq. 6-39).
enum class BladeRowKind { Rotor, Stator, InletGuideVane };

/// The cascade on one stream surface (the BladeGen section where the stream surface cuts the blade;
/// how it is taken there is the caller's, Ex. 7.1).  Angles in radians, BOOK frame; row_sign maps
/// them to the cascade frame (see the file comment).
struct CascadeGeometry
{
    ProfileFamily profile = ProfileFamily::Naca;
    Camberline camberline = Camberline::Naca65;
    BladeRowKind kind = BladeRowKind::Stator;
    int row_sign = 1;                  ///< s_row = +1 or -1: cascade-frame angle = s_row * book-frame angle
    double chord = 0.0;                ///< c [m] (for Re_c)
    double solidity = 0.0;             ///< sigma = c / s (Eq. 6-5)
    double stagger = 0.0;              ///< gamma [rad]: chord line from the meridional direction (Fig. 6-1)
    double camber = 0.0;               ///< theta = kappa1 - kappa2 [rad] (equivalent circular-arc camber for NACA 65; Eq. 6-7)
    double thickness_ratio = 0.10;     ///< t_b / c
    double max_camber_position = 0.5;  ///< a / c (Fig. 4-6); 0.5 for Naca65 and CircularArc
    double throat_over_pitch = 0.0;    ///< o / s: minimum passage opening over pitch (Fig. 6-1), for the choke limit Eq. 6-53
    double pitch_over_height = 0.0;    ///< s_bar / h_bar: mean pitch over mean blade height of the row (Eqs. 6-99, 6-100)
};

/// The flow at QN1 (upstream of the row) on the stream surface, book frame, relative to this row.
struct InletState
{
    double radius = 0.0;      ///< r1 [m]
    double wm = 0.0;          ///< W_m1 = C_m1 [m/s] (> 0)
    double wtheta = 0.0;      ///< W_theta1 = C_theta1 - omega r1 [m/s] (Eq. 7-44); beta'_1 = atan(W_theta1 / W_m1)
    double omega = 0.0;       ///< omega [rad/s], signed (book: W_theta = C_theta - omega r); 0 for a stator
    double mach = 0.0;        ///< M'_1 = W1 / a1, the relative Mach number
    double density = 0.0;     ///< rho1 [kg/m^3] (for Re_c only)
    double viscosity = 0.0;   ///< mu1 [Pa s] (for Re_c only)
};

/// The through-flow's current estimate at QN2 on the same stream surface (Sec. 7.3, p. 163).
struct DischargeEstimate
{
    double radius = 0.0;  ///< r2 [m]
    double wm = 0.0;      ///< W_m2 [m/s] (> 0); start: W_m2 = W_m1
};

/// The perfect gas used where Ch. 6 needs an equation of state (see the file comment).
struct Gas
{
    double gamma = 1.4;            ///< ratio of specific heats
    double gas_constant = 287.05;  ///< R [J/(kg K)] (Eq. 6-91 only)
};

/// Aungier's two design-incidence models (Sec. 6.2); he prefers the angle of attack (p. 125).
enum class DesignIncidenceModel {
    AngleOfAttack,  ///< Eqs. 6-8, 6-9, 6-18
    Lieblein,       ///< Eqs. 6-12 .. 6-17, 6-18 (iterated on beta1)
};

/// The design profile-loss correlation (Sec. 6.4).
enum class ProfileLossModel {
    DFactor,              ///< Eq. 6-45 (D*, Eq. 6-31; valid for D* <= 0.6)
    EquivalentDiffusion,  ///< Eq. 6-46 (D_eq*, Eq. 6-36 with W_max of Eq. 6-39); Aungier's choice (p. 134)
};

/// K1 and K2 of Eqs. 6-45 / 6-46.
enum class LossCalibration {
    TwoDimensionalCascade,  ///< K1 = 0.0035, K2 = 1 (Eq. 6-34) or K1 = 0.004, K2 = 1 (Eq. 6-37)
    Constants,              ///< CorrelationOptions::k1, k2 (p. 149: K1 = 0.004, K2 ~ 4 "rather good")
    AspectRatio,            ///< Sec. 6.11: K1 = 0.0073 (p. 150) and K2 of Eq. 6-99 (6-100 with the Reynolds correction); Eq. 6-46 only
};

struct CorrelationOptions
{
    DesignIncidenceModel incidence = DesignIncidenceModel::AngleOfAttack;
    ProfileLossModel loss = ProfileLossModel::EquivalentDiffusion;
    LossCalibration calibration = LossCalibration::AspectRatio;
    double k1 = 0.004;                       ///< with LossCalibration::Constants
    double k2 = 4.0;                         ///< with LossCalibration::Constants
    bool reynolds_correction = false;        ///< Eqs. 6-100 .. 6-102 ("not normally used", p. 150); needs AspectRatio
    std::optional<double> lieblein_alpha;    ///< overrides alpha of Eq. 6-38/6-39 (required for DoubleCircularArc at i > i*)
};

/// The loss coefficient and its parts (Eqs. 6-46, 6-58 .. 6-60, 6-78 .. 6-82).
struct LossBreakdown
{
    double design = 0.0;        ///< omega*: design (minimum-loss) profile loss at i*, Eq. 6-45 / 6-46
    double minimum = 0.0;       ///< omega_m: minimum loss after the Mach number correction, Eq. 6-58 / 6-60
    double incidence = 0.0;     ///< omega_m * (curve factor - 1): the off-design part of Eqs. 6-80 .. 6-82
    double shock = 0.0;         ///< omega_s: bow / passage shock loss, Sec. 6.7 (0 for M'_1 <= 1)
    double total = 0.0;         ///< omega = omega_s + omega_m * factor (Eqs. 6-80 .. 6-82)
    double k1 = 0.0, k2 = 0.0;  ///< the K1, K2 used in Eq. 6-45 / 6-46
};

/// What the through-flow needs from one blade row on one stream surface (Sec. 7.3), and the
/// intermediate values.  _deg fields: degrees, positive-camber cascade frame (beta = frame_sign *
/// beta_book).
struct CascadePerformance
{
    double beta2 = 0.0;          ///< discharge relative flow angle beta'_2 [rad], BOOK frame (QN2; Eq. 7-50)
    double loss = 0.0;           ///< omega, the relative total pressure loss coefficient (Eq. 7-49) = losses.total
    double wake_blockage = 0.0;  ///< B_wake (Eqs. 6-83, 6-84)
    double k_w = 1.0;            ///< K_W = 1 - B_wake (Eq. 7-6), the wake blockage factor of Eq. 7-5
    LossBreakdown losses;
    int frame_sign = 1;          ///< cascade-frame angle = frame_sign * book-frame angle (s_row; -s_row for a flipped IGV)
    double reynolds = 0.0;       ///< Re_c = rho1 W1 c / mu1 (ours; 0 when density or viscosity is not given)

    // blade angles and the operating point (cascade frame, degrees)
    double chi1_deg = 0.0, chi2_deg = 0.0;      ///< camberline angles with the chord (Ch. 4)
    double kappa1_deg = 0.0, kappa2_deg = 0.0;  ///< blade angles (4-32, 4-33)
    double beta1_deg = 0.0, beta2_deg = 0.0;    ///< inlet / discharge flow angles
    double incidence_deg = 0.0;                 ///< i = beta1 - kappa1 (6-2)
    double deviation_deg = 0.0;                 ///< delta = beta2 - kappa2 (6-3, 6-77)

    // design point (cascade frame, degrees)
    double design_attack_deg = 0.0;     ///< alpha* (6-8 or 6-17)
    double design_incidence_deg = 0.0;  ///< i* (6-18)
    double design_deviation_deg = 0.0;  ///< delta* (6-19 or 6-26)
    double design_beta1_deg = 0.0;      ///< beta1* = kappa1 + i*
    double design_beta2_deg = 0.0;      ///< beta2* = kappa2 + delta*
    double deviation_slope = 0.0;       ///< (d delta / d i)* (6-76), 0 for an inlet guide vane

    // diffusion
    double d_factor = 0.0;              ///< D at the operating point (6-31)
    double d_factor_design = 0.0;       ///< D* (6-31 at i*, delta*)
    double wmax_ratio = 0.0;            ///< W_max / W_1 at the operating point (6-39)
    double wmax_ratio_design = 0.0;     ///< (W_max / W_1)* (6-39 at i*, delta*)
    double d_eq = 0.0;                  ///< D_eq at the operating point (6-30 with 6-39)
    double d_eq_design = 0.0;           ///< D_eq* (6-36 with 6-39)

    // stall and Mach number effects (cascade frame, degrees)
    double stall_range_negative_deg = 0.0;  ///< R_c = alpha* - alpha_c (6-47, 6-49), low speed
    double stall_range_positive_deg = 0.0;  ///< R_s = alpha_s - alpha* (6-48, 6-50), low speed
    double beta1_negative_stall_deg = 0.0;  ///< beta_1c (solved with Eq. 6-47, before the 20 deg limit)
    double beta1_positive_stall_deg = 0.0;  ///< beta_1s (solved with Eq. 6-48)
    double negative_stall_incidence_deg = 0.0;  ///< i_c (6-51, then the choke limit 6-54)
    double positive_stall_incidence_deg = 0.0;  ///< i_s (6-52)
    double minimum_loss_incidence_deg = 0.0;    ///< i_m (6-55)
    double choke_beta1_deg = 0.0;               ///< beta_1,choke (6-53); NaN when the passage cannot choke (o/s rho*W*/(rho1 W1) >= 1)
    double positive_stall_margin_deg = 0.0;     ///< i_s - i
    double negative_stall_margin_deg = 0.0;     ///< i - i_c
    double xi = 0.0;                            ///< normalized incidence (6-78, 6-79)
    double critical_mach = 0.0;                 ///< M'_c (6-59)

    bool design_diffusion_limit_exceeded = false;  ///< D_eq* > 2 (p. 132) or D* > 0.6 (p. 130) at design incidence
    bool choke_limit_active = false;               ///< i_c raised to the choke limit Eq. 6-54
    bool extrapolated = false;       ///< an input lies outside the data range of the book's figures (see notes)
    std::vector<std::string> notes;  ///< every limit applied, every range left, every choice made at run time
};

/// The cascade performance of Sec. 7.3: exit relative angle, loss (with parts), stall margins and
/// K_W.  Throws CorrelationError (see the file comment).
[[nodiscard]] CascadePerformance evaluate(const CascadeGeometry& geometry, const InletState& inlet,
                                          const DischargeEstimate& discharge, const Gas& gas = {},
                                          const CorrelationOptions& options = {});

// ================================================================================================
// The single correlations (degrees, positive-camber cascade frame)
// ================================================================================================

/// The blade angles of a camberline (Ch. 4).
struct BladeAngles
{
    double chi1_deg = 0.0, chi2_deg = 0.0;      ///< camberline angles with the chord at LE and TE
    double kappa1_deg = 0.0, kappa2_deg = 0.0;  ///< blade angles with the axial direction
};

/// chi1 = chi2 = theta/2 (circular arc, NACA 65 equivalent arc, Eq. 4-35) or the parabolic arc
/// (Eqs. 4-17 .. 4-19, 0.25 < a/c < 0.75 Eq. 4-20); kappa1 = chi1 + gamma, kappa2 = gamma - chi2
/// (4-32, 4-33).  theta in [0, 90) deg.
[[nodiscard]] BladeAngles blade_angles(Camberline camberline, double stagger_deg, double camber_deg, double a_c);

/// The equivalent circular-arc camber of a NACA 65 (or A4K6) mean line from its lift coefficient:
/// tan(theta/4) = 0.05515 C_l0 / (a/c) (Eq. 6-7 = Eq. 4-29).
[[nodiscard]] double naca_camber_from_cl0_deg(double cl0, double a_c);

/// K_sh: 1.0 NACA, 1.1 C.4, 0.7 DCA (p. 121).
[[nodiscard]] double profile_shape_factor(ProfileFamily profile) noexcept;

/// K_t,i = (10 t_b/c)^q, q = 0.28 / [0.1 + (t_b/c)^0.3] (Eqs. 6-10, 6-11).
[[nodiscard]] double incidence_thickness_factor(double tb_c);

/// alpha* = [3.6 K_sh K_t,i + 0.3532 theta (a/c)^0.25] sigma^e, e = 0.65 - 0.002 theta (6-8, 6-9).
[[nodiscard]] double design_attack_angle_deg(double ksh, double kti, double camber_deg, double a_c, double solidity);

/// (i*_0)_10 = beta1^p / (5 + 46 exp(-2.3 sigma)) - 0.1 sigma^3 exp[(beta1 - 70)/4],
/// p = 0.914 + sigma^3/160 (6-13, 6-14).  beta1 >= 0.
[[nodiscard]] double lieblein_zero_camber_incidence_deg(double beta1_deg, double solidity);

/// n = 0.025 sigma - 0.06 - (beta1/90)^(1 + 1.2 sigma) / (1.5 + 0.43 sigma) (6-15).  beta1 >= 0.
[[nodiscard]] double lieblein_incidence_slope(double beta1_deg, double solidity);

/// Lieblein's design incidence with Aungier's pseudo-blade extension (6-12, 6-16 .. 6-18).
struct LiebleinIncidence
{
    double pseudo_kappa1_deg = 0.0;     ///< kappa1_bar = |gamma| + theta/2 (6-16)
    double pseudo_incidence_deg = 0.0;  ///< i*_bar = K_sh K_t,i (i*_0)_10 + n theta at beta1 = kappa1_bar + i*_bar (6-12)
    double attack_deg = 0.0;            ///< alpha* = kappa1_bar + i*_bar - |gamma| + (a/c - 0.5) theta (6-17)
    double incidence_deg = 0.0;         ///< i* = alpha* + gamma - kappa1 (6-18)
    int iterations = 0;
};
[[nodiscard]] LiebleinIncidence lieblein_design_incidence(double ksh, double kti, double stagger_deg, double camber_deg,
                                                          double a_c, double kappa1_deg, double solidity);

/// (delta*_0)_10 = 0.01 sigma beta1 + [0.74 sigma^1.9 + 3 sigma] (beta1/90)^(1.67 + 1.09 sigma) (6-20).
[[nodiscard]] double lieblein_zero_camber_deviation_deg(double beta1_deg, double solidity);

/// m = m_1.0 / sigma^b (6-21); m_1.0 NACA 65 (6-22) or circular arc (6-23), x = beta1/100;
/// b = 0.9625 - 0.17 x - 0.85 x^3 (6-24).  Only Naca65 and CircularArc camberlines (p. 126).
[[nodiscard]] double lieblein_deviation_slope(Camberline camberline, double beta1_deg, double solidity);

/// K_t,delta = 6.25 (t_b/c) + 37.5 (t_b/c)^2 (6-25).
[[nodiscard]] double deviation_thickness_factor(double tb_c);

/// delta* = K_sh K_t,delta (delta*_0)_10 + m theta (6-19).
[[nodiscard]] double lieblein_design_deviation_deg(Camberline camberline, double ksh, double ktd, double beta1_deg,
                                                   double camber_deg, double solidity);

/// Howell's design deviation modified by Lieblein's corrections (6-26):
/// delta* = [0.92 (a/c)^2 + 0.002 kappa2] / [1 - 0.002 theta/sqrt(sigma)] theta/sqrt(sigma)
///          + (K_sh K_t,delta - 1)(delta*_0)_10.
[[nodiscard]] double howell_design_deviation_deg(double ksh, double ktd, double a_c, double kappa2_deg, double camber_deg,
                                                 double beta1_deg, double solidity);

/// The velocity ratios of a stream surface: W2/W1 = (W_m2/W_m1) cos(beta1) / cos(beta2) (from
/// W = W_m / cos beta).
[[nodiscard]] double velocity_ratio_w2_w1(double beta1_deg, double beta2_deg, double wm_ratio);

/// D = 1 - W2/W1 + (W_theta1 - W_theta2) / (2 sigma W1) (6-31).
[[nodiscard]] double diffusion_factor(double beta1_deg, double beta2_deg, double solidity, double wm_ratio);

/// D ~ (W_max - W2)/W1 (6-29) and D_eq ~ W_max/W2 (6-30) from velocities.
[[nodiscard]] double diffusion_factor_from_velocities(double w_max, double w1, double w2);
[[nodiscard]] double equivalent_diffusion_from_velocities(double w_max, double w2);

/// Lieblein's design (W_max/W1)* = 1.12 + 0.61 cos^2(beta1)/sigma [tan beta1 - tan beta2] (6-35).
[[nodiscard]] double lieblein_wmax_ratio_design(double beta1_deg, double beta2_deg, double solidity);

/// Klapproth's general form (6-39):  W_max/W1 = 1.12 + 0.61 cos^2(beta1)/sigma
///   (r1 C_theta1 - r2 C_theta2)/(r1 C_m1) + alpha (i - i*)^1.43,
/// written with W_theta = W_m tan(beta), C_theta = W_theta + omega r, C_m = W_m:
///   (r1 C_theta1 - r2 C_theta2)/(r1 C_m1) = tan b1 - (r2/r1)(W_m2/W_m1) tan b2 + (U1/W_m1)(1 - (r2/r1)^2).
/// incidence_excess_deg = i - i* >= 0 (6-38: "for i >= i*").
[[nodiscard]] double wmax_ratio(double beta1_deg, double beta2_deg, double solidity, double wm_ratio, double radius_ratio,
                                double blade_speed_ratio, double incidence_excess_deg, double alpha);

/// alpha of Eq. 6-38: 0.0117 NACA 65, 0.007 C.4 (p. 132); not given for DCA (throws).
[[nodiscard]] double lieblein_alpha(ProfileFamily profile);

/// D_eq = (W_max/W1)(W1/W2) = (W_max/W1) cos(beta2) W_m1 / (cos(beta1) W_m2) (6-30, 6-36).
[[nodiscard]] double equivalent_diffusion(double wmax_ratio, double beta1_deg, double beta2_deg, double wm_ratio);

/// omega* from omega* cos(beta2*)/(2 sigma) = K1 [K2 + 3.5 D* + 37 D*^4] (6-45; 6-34 with K1 = 0.0035,
/// K2 = 1).  D* <= 0.6 (p. 130).
[[nodiscard]] double design_loss_d_factor(double d_star, double beta2_deg, double solidity, double k1, double k2);

/// omega* from omega* cos(beta2*)/(2 sigma) (W1*/W2*)^2 = K1 [K2 + 3.1 (D_eq* - 1)^2 + 0.4 (D_eq* - 1)^8]
/// (6-46; 6-37 with K1 = 0.004, K2 = 1).
[[nodiscard]] double design_loss_equivalent_diffusion(double d_eq_star, double beta2_deg, double solidity, double w1_w2,
                                                      double k1, double k2);

/// omega = 2 (theta_w/c) sigma / cos(beta2) (W2/W1)^2 (6-27).
[[nodiscard]] double wake_loss_coefficient(double theta_w_over_c, double solidity, double beta2_deg, double w2_w1);

/// Howell (6-40 .. 6-44): omega = C_D sigma cos^2(beta1)/cos^3(beta_m), tan beta_m = (tan b1 + tan b2)/2;
/// C_Da = 0.02 s/h; C_Ds = 0.18 C_L^2; C_L = 2 cos(beta_m)(tan b1 - tan b2)/sigma.
[[nodiscard]] double howell_mean_angle_deg(double beta1_deg, double beta2_deg);
[[nodiscard]] double howell_loss_from_drag(double cd, double solidity, double beta1_deg, double beta2_deg);
[[nodiscard]] double howell_endwall_drag(double pitch_over_height);
[[nodiscard]] double howell_lift_coefficient(double beta1_deg, double beta2_deg, double solidity);
[[nodiscard]] double howell_secondary_drag(double lift_coefficient);

/// K2 = 1 + (s/h) cos(beta2) (6-99); + 0.004 K_Re / K1 (6-100).
[[nodiscard]] double aspect_ratio_k2(double pitch_over_height, double beta2_deg, double k_re, double k1);

/// K_Re = sqrt(2.5e5/Re_c) - 1 for Re_c < 2.5e5 (6-101), [log(2.5e5)/log(Re_c)]^2.58 - 1 above (6-102).
[[nodiscard]] double reynolds_factor(double reynolds);

/// The low-speed stall ranges of Herrig et al. as fitted by Aungier (6-47, 6-48), solved with
/// beta_1c = gamma + alpha_c and beta_1s = gamma + alpha_s (p. 135: "a simple iterative solution").
struct StallRanges
{
    double negative_deg = 0.0;        ///< R_c = alpha* - alpha_c (6-49)
    double positive_deg = 0.0;        ///< R_s = alpha_s - alpha* (6-50)
    double beta1c_deg = 0.0;          ///< the solved beta_1c
    double beta1s_deg = 0.0;          ///< the solved beta_1s (closed form: Eq. 6-48 is linear in beta_1s)
    bool beta1c_limited = false;      ///< Eq. 6-47 was evaluated at beta_1c = 20 deg (Aungier's limit, p. 135)
};
/// alpha_c - alpha* = -9 + [1 - (30/beta_1c)^0.48] theta/4.176 (6-47), beta_1c >= 20 deg;
/// alpha_s - alpha* = 10.3 + [2.92 - beta_1s/15.6] theta/8.2 (6-48).
[[nodiscard]] double negative_stall_attack_offset_deg(double beta1c_deg, double camber_deg);
[[nodiscard]] double positive_stall_attack_offset_deg(double beta1s_deg, double camber_deg);
[[nodiscard]] StallRanges stall_ranges(double design_attack_deg, double stagger_deg, double camber_deg);

/// i_c = i* - R_c / [1 + 0.5 M'^3] (6-51); i_s = i* + R_s / [1 + 0.5 (K_sh M')^3], K_sh <= 1 (6-52).
[[nodiscard]] double mach_negative_stall_incidence_deg(double design_incidence_deg, double range_negative_deg, double mach1);
[[nodiscard]] double mach_positive_stall_incidence_deg(double design_incidence_deg, double range_positive_deg, double mach1,
                                                       double ksh);

/// Perfect gas: rho1 W1 / (rho* W*) = M [(gamma + 1)/(2 + (gamma - 1) M^2)]^((gamma+1)/(2(gamma-1)))
/// (the isentropic mass-flux ratio of the relative total conditions).
[[nodiscard]] double mass_flux_ratio_to_sonic(double mach, double gamma);

/// W*/W1 = sqrt[(2 + (gamma - 1) M^2) / ((gamma + 1) M^2)] (perfect gas, same relative total enthalpy).
[[nodiscard]] double sonic_velocity_ratio(double mach, double gamma);

/// The inlet angle at choke: rho1 W1 s cos(beta_1,choke) = o rho* W* (6-53); NaN when
/// (o/s) rho* W* / (rho1 W1) >= 1 (no inlet angle chokes the passage).
[[nodiscard]] double choke_inlet_angle_deg(double throat_over_pitch, double mach1, double gamma);

/// i_m = i_c + (i_s - i_c) R_c / (R_c + R_s) (6-55).
[[nodiscard]] double minimum_loss_incidence_deg(double ic_deg, double is_deg, double range_negative_deg,
                                                double range_positive_deg);

/// M'_c = M'_1 W* / W_max = M'_1 (W*/W1) / (W_max/W1) (6-59).
[[nodiscard]] double critical_mach_number(double mach1, double sonic_ratio, double wmax_ratio_design);

/// omega_m = omega* [1 + (i_m - i*)^2 / R_s^2] (6-58), + K_sh [(M'/M'_c - 1) W*/W1]^2 when M'_1 > M'_c
/// with M' <= 1 (6-60).
[[nodiscard]] double minimum_loss(double design_loss, double im_deg, double design_incidence_deg, double range_positive_deg,
                                  double ksh, double mach1, double critical_mach, double sonic_ratio);

/// xi = (i - i_m)/(i_s - i_m) for i >= i_m (6-78), (i - i_m)/(i_m - i_c) for i < i_m (6-79).
[[nodiscard]] double normalized_incidence(double i_deg, double im_deg, double ic_deg, double is_deg);

/// The off-design loss factor: 1 + xi^2 for -2 <= xi <= 1 (6-80), 5 - 4 (xi + 2) for xi < -2 (6-81),
/// 2 + 2 (xi - 1) for xi > 1 (6-82); omega = omega_s + omega_m * factor.
[[nodiscard]] double off_design_loss_factor(double xi) noexcept;

/// (d delta/d i)* = [1 + (sigma + 0.25 sigma^4)(beta1/53)^2.5] / exp(3.1 sigma) (6-76).
[[nodiscard]] double deviation_incidence_slope(double beta1_deg, double solidity);

/// delta = delta* + (d delta/d i)* (i - i*) + 10 (1 - W_m2/W_m1) (6-77), without the incidence term for
/// an inlet guide vane (p. 144).
[[nodiscard]] double off_design_deviation_deg(double design_deviation_deg, double slope, double i_deg,
                                              double design_incidence_deg, double wm_ratio, bool inlet_guide_vane);

/// B_wake = 0 for D_eq <= 2 (6-83), 1 - (2/D_eq)^0.9 above (6-84).
[[nodiscard]] double wake_blockage(double d_eq);

// --- Sec. 6.7: shock loss (Swan 1961, perfect gas) ----------------------------------------------

/// The suction-surface expansion of Fig. 6-17 for a sharp-LE circular-arc suction surface:
/// tan(theta_u/4) = tan(theta/4) + t_b/c (6-65); c/(2 R_u) = sin(theta_u/2) (6-61 corrected by Eq. 4-7);
/// psi = 90 - theta_u/2 - gamma (6-67, 6-68); tan(phi) = s cos(psi) / (s sin(psi) + R_u) (6-69).
struct ShockGeometry
{
    double suction_camber_deg = 0.0;  ///< theta_u
    double psi_deg = 0.0;             ///< psi
    double radius_over_pitch = 0.0;   ///< R_u / s (infinite for theta_u = 0)
    double expansion_deg = 0.0;       ///< phi
};
[[nodiscard]] ShockGeometry shock_expansion_geometry(double camber_deg, double tb_c, double stagger_deg, double solidity);

/// The Prandtl-Meyer angle nu(M) [rad] of a perfect gas (the closed form of Eq. 6-70) and its inverse.
[[nodiscard]] double prandtl_meyer_angle(double mach, double gamma);
[[nodiscard]] double prandtl_meyer_mach(double nu, double gamma);

/// The total pressure ratio P_t,out / P_t,in across a normal shock (perfect gas; Eqs. 6-72 .. 6-75).
[[nodiscard]] double normal_shock_total_pressure_ratio(double mach, double gamma);

/// Sec. 6.7: M'_S from nu(M'_S) = nu(M'_1) + phi (6-70), M'_in = sqrt(M'_1 M'_S) (6-71), the normal-shock
/// total pressure loss at M'_in, as omega_s = (P't1 - P't2)/(P't1 - P1).  0 for M'_1 <= 1.
struct ShockLoss
{
    ShockGeometry geometry;
    double surface_mach = 0.0;  ///< M'_S
    double shock_mach = 0.0;    ///< M'_in
    double total_pressure_ratio = 1.0;
    double loss = 0.0;          ///< omega_s
};
[[nodiscard]] ShockLoss shock_loss(double mach1, double camber_deg, double tb_c, double stagger_deg, double solidity,
                                   double gamma);

// --- Sec. 6.9: tip clearance loss (blade-row level) -------------------------------------------

struct TipClearanceInput
{
    double clearance = 0.0;      ///< delta_c [m]
    double tip_radius = 0.0;     ///< r_tip [m]
    int blade_count = 0;         ///< Z
    double chord = 0.0;          ///< c [m] at the tip
    double stagger_rad = 0.0;    ///< gamma [rad] at the tip (cos gamma: the axial chord)
    int row_number = 1;          ///< N_row, sequential through the compressor (6-87)
    double r1 = 0.0, r2 = 0.0;              ///< inlet / exit radius of the clearance stream surface [m]
    double rho1 = 0.0, rho2 = 0.0;          ///< densities there [kg/m^3]
    double cm1 = 0.0, cm2 = 0.0;            ///< meridional velocities there [m/s]
    double ctheta1 = 0.0, ctheta2 = 0.0;    ///< absolute tangential velocities there [m/s]
    double mean_density = 0.0;   ///< rho_bar (6-87, 6-88)
    double mass_flow = 0.0;      ///< m_dot of the blade row [kg/s] (6-89)
};
struct TipClearanceLoss
{
    double torque = 0.0;               ///< tau, the torque on the gap flow [N m] (6-85)
    double pressure_difference = 0.0;  ///< Delta P across the blade at the gap [Pa] (6-86)
    double leakage_velocity = 0.0;     ///< U_c [m/s] (6-87)
    double leakage_mass_flow = 0.0;    ///< m_dot_c [kg/s] (6-88)
    double total_pressure_loss = 0.0;  ///< Delta P_t [Pa], blade-row average (6-89)
};
/// tau = pi delta_c [(r rho C_m)_1 + (r rho C_m)_2][r2 C_theta2 - r1 C_theta1] (6-85),
/// Delta P = tau / (Z r_tip delta_c c cos gamma) (6-86), U_c = 0.816 sqrt(2 Delta P / rho_bar) / N_row^0.2
/// (6-87), m_dot_c = rho_bar U_c Z delta_c c cos gamma (6-88), Delta P_t = Delta P m_dot_c / m_dot (6-89).
/// |tau| is used (the loss does not depend on the sense of the work; ours).  The through-flow
/// distributes Delta P_t linearly across the annulus, zero at the opposite wall (p. 147): leakage_loss_at().
[[nodiscard]] TipClearanceLoss tip_clearance_loss(const TipClearanceInput& in);

// --- Sec. 6.10: shroud seal leakage loss (blade-row level, Egli 1935) --------------------------

/// C_r = 1 - 1 / (3 + [54.3 / (1 + 100 delta_c/t)]^3.45) (6-92 = 2000 Eq. 4-34).
[[nodiscard]] double seal_contraction_ratio(double clearance_over_fin_thickness);
/// C_t = 2.143 [ln N - 1.464] / (N - 4.322) [1 - P_R]^(0.375 P_R) (2000 Eq. 4-35; 2003 Eq. 6-93 misprinted),
/// P_R = lower / higher pressure, 0 <= P_R < 1.
[[nodiscard]] double seal_throttling_coefficient(int fins, double pressure_ratio);
/// C_c = 1 + X1 [delta_c/p - X2 ln(1 + delta_c/p)] / (1 - X2) (6-94), delta_c/p <= X2 - 1;
/// X1 = 15.1 - 0.05255 exp[0.507 (12 - N)], X2 = 1.058 + 0.0218 N for N <= 12 (6-95, 6-96);
/// X1 = 13.15 + 0.1625 N, X2 = 1.32 for N > 12 (2000 Eqs. 4-36, 4-37).
[[nodiscard]] double seal_carryover_coefficient(int fins, double clearance_over_fin_pitch);

struct ShroudSealInput
{
    double seal_radius = 0.0;          ///< r_seal [m]
    double clearance = 0.0;            ///< delta_c [m]
    int fins = 0;                      ///< N
    double fin_thickness = 0.0;        ///< t [m] (Fig. 6-22)
    double fin_pitch = 0.0;            ///< p [m] (Fig. 6-22)
    double pressure_ratio = 0.0;       ///< P_R = lower / higher pressure across the seal
    double upstream_density = 0.0;     ///< rho on the higher-pressure side [kg/m^3] (2000 p. 63)
    double upstream_temperature = 0.0; ///< T there [K]
    double gas_constant = 0.0;         ///< R [J/(kg K)]
    double pressure_difference = 0.0;  ///< Delta P across the blade row [Pa] (6-90)
    double mass_flow = 0.0;            ///< m_dot of the blade row [kg/s]
};
struct ShroudSealLoss
{
    double contraction_ratio = 0.0, throttling_coefficient = 0.0, carryover_coefficient = 0.0;
    double leakage_mass_flow = 0.0;    ///< m_dot_leak = 2 pi r_seal delta_c C_t C_c C_r rho sqrt(R T) (6-91)
    double total_pressure_loss = 0.0;  ///< Delta P_t = Delta P m_dot_leak / m_dot (6-90)
};
[[nodiscard]] ShroudSealLoss shroud_seal_loss(const ShroudSealInput& in);

// --- the row-level call for the through-flow (contract Sec. 4.1 B) ------------------------------

enum class EndWall { Hub, Shroud };

/// One end-wall state of the through-flow (the wall stream surface on the gap side), book frame.
struct WallState
{
    double radius = 0.0;       ///< r [m]
    double density = 0.0;      ///< rho [kg/m^3]
    double cm = 0.0;           ///< C_m [m/s]
    double ctheta = 0.0;       ///< C_theta [m/s], absolute, book sign
    double pressure = 0.0;     ///< static P [Pa]
    double temperature = 0.0;  ///< static T [K]
};

/// The row's end-wall gap: an unshrouded tip / hub clearance (Sec. 6.9) or a shroud seal (Sec. 6.10).
struct RowWallData
{
    bool shrouded = false;               ///< false: tip clearance Eqs. 6-85 .. 6-89; true: seal Eqs. 6-90 .. 6-96
    EndWall gap_side = EndWall::Shroud;  ///< the wall of the gap (rotor tip: shroud; cantilevered stator: hub)
    int row_sign = 1;                    ///< s_row (C_theta enters 6-85 in cascade signs)
    double clearance = 0.0;              ///< delta_c [m]
    double gap_radius = 0.0;             ///< r_tip (6-86) or r_seal (6-91) [m]
    int blade_count = 0;                 ///< Z
    double chord = 0.0;                  ///< c at the gap [m]
    double stagger = 0.0;                ///< gamma at the gap [rad]
    int row_number = 1;                  ///< N_row (6-87)
    int fins = 0;                        ///< N (seal)
    double fin_thickness = 0.0;          ///< t [m] (seal)
    double fin_pitch = 0.0;              ///< p [m] (seal)
};

struct RowLeakage
{
    double total_pressure_loss = 0.0;  ///< Delta P_t [Pa], the row average (6-89 or 6-90)
    double leakage_mass_flow = 0.0;    ///< m_dot_c (6-88) or m_dot_leak (6-91) [kg/s], >= 0
    int seal_flow_direction = 0;       ///< seal: +1 from the discharge (QN2) to the inlet (QN1) when P2 > P1, -1 otherwise; 0 for a tip gap
    TipClearanceLoss clearance;        ///< the tip-gap details (when !shrouded)
    ShroudSealLoss seal;               ///< the seal details (when shrouded)
};

/// The leakage loss of one row from the gap-side wall states at QN1 and QN2 and the row mass flow.
/// Ours where the book is silent (contract Sec. 2.3): tip gap rho_bar = (rho1 + rho2)/2; seal:
/// P_R = min(P1, P2)/max(P1, P2), Delta P = |P2 - P1|, rho and T of the higher-pressure side
/// (Aungier 2000 p. 63: "evaluated on the upstream (higher pressure) side of the seal").
[[nodiscard]] RowLeakage endwall_leakage(const RowWallData& row, const WallState& qn1, const WallState& qn2,
                                         double mass_flow, const Gas& gas = {});

/// The leakage loss on the stream surface at mass fraction psi = m_dot(y)/m_dot(y_s) (Eq. 7-54) from the
/// hub: linear in psi, zero at the wall opposite the gap, mass-weighted mean equal to the row value
/// (p. 147): 2 Delta P_t psi (gap at the shroud) or 2 Delta P_t (1 - psi) (gap at the hub).  Linear in
/// mass fraction and mass-weighted are ours (the book: "a linear distribution across the annulus").
[[nodiscard]] double leakage_loss_at(double row_total_pressure_loss, double mass_fraction, EndWall gap_side);

/// omega_bar = Delta P't / (P't1 - P1) (Eq. 7-49), for adding a pressure loss to the stream surface's
/// coefficient (ours: the book does not state the conversion).
[[nodiscard]] double loss_coefficient_from_pressure(double total_pressure_loss, double relative_total_pressure1,
                                                    double static_pressure1);

}  // namespace BladeGen::Turbolenza
