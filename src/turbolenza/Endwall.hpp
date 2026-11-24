/// @file turbolenza/Endwall.hpp
/// @brief The end-wall (hub and shroud) boundary layers of Aungier (2003) Ch. 8, marched along the
///        end walls through the quasi-normal stations of the meridional through-flow (Ch. 7), and
///        the end-wall blockage factor K_B of Eq. 7-5.
///
/// Not a v1 port.  Every equation is Aungier's, read from the page images of R. H. Aungier,
/// Axial-Flow Compressors (ASME Press, 2003), Ch. 8, printed pp. 175-198 (Eqs. 8-1 to 8-77), with
/// Ch. 3 Sec. 3.6 (pp. 54-56, Eqs. 3-44 to 3-60) and Ch. 7 (p. 157 Eq. 7-5, p. 164 Sec. 7.4,
/// p. 166 Fig. 7-3).  Every choice Aungier does not make is marked "ours".
///
/// ============================== CONVENTIONS (read first) ==============================
///  * SI units throughout; angles in radians; rotation speeds in rad/s.  BOOK SIGNS: omega > 0
///    turns toward +theta, W_theta = C_theta - omega r [8-67]; tangential velocities signed with
///    +theta.  TURBOlenza's frame maps to it with theta_book = theta_T, C_theta,book = C_theta,T and
///    omega_book = -2 pi rpm / 60 (the q3d data contract, Sec. 0.3); every equation is then used
///    verbatim.
///  * The analysis is a pure one-station step (endwall_step: the layer at k-1 and the edges at k-1
///    and k -> the layer at k), because Fig. 7-3 (p. 166) calls the end-wall analysis at each
///    quasi-normal inside the flow-field loop and may repeat it; endwall_inlet starts a wall,
///    endwall_blockage gives B and K_B of a station, analyze_endwalls is the convenience march.
///  * Stations are the through-flow quasi-normals, k = 0 .. N-1 from inlet to exit.  Each station
///    carries one EndwallWallStation per wall (hub, shroud): the point where the quasi-normal meets
///    the end-wall contour, and the inviscid edge state of the end-wall stream surface there.
///  * m: meridional arc length ALONG THAT WALL (each wall has its own), strictly increasing in k.
///  * phi: the wall contour angle with the axial direction, sin(phi) = dr/dm along the wall.
///  * Edge velocities are given in the ABSOLUTE frame: Vm = C_m = W_m = V_m [8-3] (> 0) and
///    Ctheta = C_theta.  The module forms the velocity relative to the wall itself,
///    V_theta = C_theta - r omega_w [8-4], and relative to a blade row, W_theta = C_theta - r omega.
///  * wall_omega at station k is the rotation speed of the wall SEGMENT from station k-1 to
///    station k (rotor drum, rotating shroud band, stationary casing ...).  At station 0 it is the
///    frame of the inlet data.  A wall rotation change is therefore always located at a station,
///    and the layer reported at station k is in the frame of the segment ending at k.
///  * A blade row lies between its inlet station (leading edge) and discharge station (trailing
///    edge).  In the Ch. 7 through-flow a row is one step, station k-1 -> k ("the meridional
///    integration step is always across the blade row", p. 303): endwall_step with the row's data.
///    analyze_endwalls also accepts stations inside a row; with the Ch. 8 model the row then keeps
///    one nu over all its intervals (ours).  The two walls of a row have their own data.
///  * The through-flow mass flow mdot is the compressor mass flow of Eqs. 8-72 and 8-76.
///  * Output per station and wall: the layer (EndwallLayer, in the frame of the segment ending at
///    that station); per station: B [8-76] and K_B = 1 - B [8-77] for Eq. 7-5 (Ch. 7 p. 164: K_B is
///    an area correction on the mass flow with the hub and shroud stream surfaces kept on the
///    walls; the displacement thicknesses delta1* are reported for the alternative of Sec. 7.4,
///    repositioning the end-wall stream surfaces).  Per row and wall: the blade force defect.
///  * Ch. 8 sends nothing else back to the through-flow: the blade force defects stay inside 8-6
///    and 8-7, end-wall losses and work are left to the empirical models of Ch. 6 (p. 175).  The
///    "limits and numerical damping" of the coupled iteration (p. 195) are not given by the book
///    and are the through-flow's to choose (none here).
///  * Reynolds numbers: c_f [8-36] takes rho_e V_e theta / mu of the streamwise theta, raised to 250
///    (p. 186); the layer is kept turbulent by rho_e V_e theta11 / mu >= 250 on the marched state
///    [8-73] (p. 194).  Both limits are applied, each where the book states it.
/// ========================================================================================
///
/// The equations (wall frame: V relative to the wall, omega_w the wall's rotation speed; e = edge):
///  * entrainment [8-5 = 3-51]:   d/dm[r rho_e V_me (delta - delta1*)] = r rho_e V_e E
///  * meridional momentum [8-6 = 3-52]:
///      d/dm[r rho_e V_me^2 theta11] + delta1* r rho_e V_me dV_me/dm
///        - rho_e V_theta_e sin(phi) [V_theta_e (delta2* + theta22) + 2 omega_w r delta2*] = r [tau_mw + f_me nu_m]
///  * tangential momentum [8-7 = 3-53, with 8-9 = 8-10]:
///      d/dm[r^2 rho_e V_me V_theta_e theta12] + r delta1* rho_e V_me [r dV_theta_e/dm + sin(phi)(V_theta_e + 2 r omega_w)]
///        = r^2 [tau_theta_w + f_theta_e nu_theta],   r dV_theta_e/dm + sin(phi)(V_theta_e + 2 r omega_w) = d(r C_theta_e)/dm
///  * the edge blade forces [8-8, 8-9]:
///      f_me = rho_e V_me dV_me/dm + dP_e/dm - sin(phi)/r rho_e (V_theta_e + omega_w r)^2,
///      f_theta_e = rho_e V_me / r d(r C_theta_e)/dm
///    (8-5 .. 8-7 are the integral forms of Ch. 3, 3-51 .. 3-53, with V for W and omega_w for omega;
///    they are used as printed.  Not used here but flagged: the differential forms printed on p. 55
///    lack r inside d/dm of continuity 3-46 and carry -d tau/dy in 3-47, 3-48 where 3-32, 3-36 have
///    +d tau/dy; the integral forms are consistent with r restored and +d tau/dy.  3-52 / 8-6 drop a
///    centrifugal density-defect term sin(phi) omega_w^2 r^2 int (rho_e - rho) dy, zero for the
///    constant density of Sec. 8.3; 8-69 .. 8-71 are exact for constant density only);
///  * the thicknesses [8-11 .. 8-17] (defined as integrals over the layer, V relative to the wall);
///  * power-law profiles [8-18, 8-19]: V_m/V_me = (y/delta)^n, V_theta/V_theta_e = (y/delta)^m, constant
///    density, which give [8-20 .. 8-26]:
///      n = theta11/(delta - delta1* - 2 theta11), m = theta12 (n+1)^2 / [delta - theta12 (n+1)],
///      H1 = delta1*/theta11 = 2n + 1, delta - delta1* = 2 H1 theta11/(H1 - 1), delta1*/delta = n/(n+1),
///      H2 = delta2*/theta22 = 2m + 1, delta2*/delta = m/(m+1);
///  * the extended tangential profile for m < 0.05 [8-28]:
///      V_theta/V_theta_e = (y/delta)^0.05 + 0.1705 (1 - 20m)(1 - y/delta)^2 (y/delta)^0.1, which gives
///      delta2*/delta = 20m/21 [8-29], theta22/delta = 0.95m - 1.684m^2 [8-30],
///      m = 0.05 + [theta12/delta - 0.05/((n+1)(n+1.05))] (n+1.1)(n+2.1)(n+3.1)/6.82 [8-31],
///      H2 = 20 / [21 (0.95 - 1.684m)] [8-32]
///      (8-29 and 8-31 are the exact integrals of 8-28; 8-30 as printed is not: the integral is
///      -4.38e-4 + 0.9502 m - 1.5126 m^2 (test_endwall.cpp); kept as printed);
///  * limits (p. 184): H1 <= 2.4 and H2 <= 2.4, "the boundary layer profile separation limit"
///    (Fig. 8-4 draws its separation-limit line at H2 = 2.2; the text value 2.4 is used);
///  * wall shear (Sec. 8.4): c_f = tau_w / (rho_e V_e^2 / 2) [8-33], directed along the edge streamline,
///      tau_mw = c_f rho_e V_e V_me / 2 [8-34], tau_theta_w = c_f rho_e V_e V_theta_e / 2 [8-35],
///      c_f = 0.246 exp(-1.561 H) (rho_e V_e theta / mu)^-0.268 (Ludwieg and Tillmann 1950) [8-36],
///      the same correlation as BoundaryLayer's Eq. 5-134 (without 5-135's compressibility factor);
///      the streamwise thicknesses (V_s = V_m cos(beta_e) + V_theta sin(beta_e) [8-37], [8-38, 8-39]):
///      delta* = delta1* cos^2 + delta2* sin^2 [8-40],
///      theta = (theta11 + delta1*) cos^4 + (theta22 + delta2*) sin^4 + 2 (theta12 + delta1*) sin^2 cos^2 - delta* [8-41],
///      all thicknesses on the right of 8-40 and 8-41 limited to >= 0, H = delta*/theta <= 2.4 and
///      rho_e V_e theta / mu >= 250 in 8-36 (p. 186);
///  * entrainment: E = 0.025 (H - 1) [8-43] (default: blade rows, the axial-compressor application,
///    p. 187), or E = 0.025 (H1 - 1) [8-42] (blade-free passages); the form of BoundaryLayer's
///    head_entrainment [5-133], called with H or H1 (ours: with H after its 2.4 limit, which p. 186
///    states for 8-36 only; H1 is the limited one in any case);
///  * blade force defects (Sec. 8.5), nu_m = nu_theta = nu (Balsa and Mellor 1975, p. 189), constant
///    across a blade row, nu = nu0 + dnu [8-56]:
///      g = s cos(gamma) (averages of inlet and discharge) [8-2],
///      nu0 = (0.12 g + delta_c/2)(8 theta11/g)^3 / [1 + (8 theta11/g)^3], theta11 averaged [8-44],
///      the collateral-flow condition at the discharge V_theta_e theta12 - W_theta_e theta11 = 0 [8-45, 8-46]
///      (printed without the edge subscript e, as are V_theta in 8-62 and 8-66; edge values meant),
///      xi = [rbar/r + Wbar_theta_e W_theta_e / (Vbar_me V_me)] rbar f_theta_e dm / (r rho_e V_me) [8-53]
///      (from 8-47, 8-48 with f_me = -Wbar_theta_e f_theta_e / Vbar_me [8-49 .. 8-51]),
///      F = (1/2) (8 theta11/g)^3 / [1 + (8 theta11/g)^3] [8-55],
///      dnu -> [dnu - (V_theta_e theta12 - W_theta_e theta11)/xi] F, updated on every iteration [8-57]
///      (its fixed point is 8-54: dnu = -(V_theta_e theta12 - W_theta_e theta11)_{dnu=0} F / xi);
///  * shroud seal leakage (Sec. 8.6, Fig. 8-7; mdot_leak > 0 from the row discharge to its inlet,
///    from Ch. 6 Sec. 6.10, supplied by the caller):
///      inlet:     Delta+ = mdot_leak / (2 pi r rho_e V_me)_in [8-59], (delta - delta1*)+ = (delta - delta1*)_in + Delta+ [8-58],
///                 theta11+ = theta11_in + Delta+ [8-60],
///                 (V_theta theta12)+ = (V_theta theta12)_in + Delta+ V_theta_in - (Delta+ + |Delta+|) U_leak / 2 [8-62],
///      U_leak = r omega_c / 2 + V_theta - C_theta, the leakage's tangential velocity relative to the shroud [8-61],
///      discharge: Delta- = mdot_leak / (2 pi r rho_e V_me)_out [8-63], (delta - delta1*)_out = (delta - delta1*)- - Delta- [8-64],
///                 theta11_out = theta11- - Delta- [8-65],
///                 (V_theta theta12)_out = (V_theta theta12)- - Delta- V_theta_out + (Delta- - |Delta-|) U_leak / 2 [8-66];
///    the layer across a shrouded row is computed with the wall rotating at the blade speed (p. 193);
///  * jump conditions (Sec. 8.7): W_theta = C_theta - omega r [8-67] gives
///      C_theta_e theta12 = W_theta_e theta12' [8-68], C_theta_e delta2* = W_theta_e delta2*' [8-69],
///      C_theta_e^2 theta22 = W_theta_e^2 theta22' + omega r W_theta_e delta2*' [8-70],
///      C_theta_e^2 theta22 (1 + H2) = W_theta_e^2 theta22' (1 + H2') + 2 omega r W_theta_e delta2*' [8-71];
///    8-68 is the quantity under d/dm in 8-7, so the march carries (delta - delta1*), theta11 and
///    q12 = V_theta_e theta12 across a wall rotation change without special logic (p. 193-194).  The
///    bracket of 8-6 (8-71) is NOT carried: it is rebuilt from the profile in the new frame, whose m
///    follows from theta12' = q12 / V_theta_e' by 8-21 or 8-31; the profile family does not keep 8-69
///    and 8-71 exactly (test_endwall.cpp measures -2.9 % and -3.0 % at a stationary-to-rotating hub);
///  * procedure (Sec. 8.8): turbulent layers only; at the inlet H1 = 1.4, m = 1/7 (H2 = 1.286),
///    2 pi r rho_e V_me delta1* = B mdot / 2 for each wall [8-72], B = 0.02 just upstream of the first
///    row; rho_e V_e theta11 / mu >= 250 throughout [8-73]; per station: initial theta11 and theta12 from
///    d/dm[r rho_e V_me^2 theta11] = 0 [8-74] and d/dm[r^2 rho_e V_me V_theta_e theta12] = 0 [8-75] with H1,
///    H2 held, then steps 2-6 (closures, integrate 8-5 .. 8-7, n and delta from 8-20 and 8-24 with
///    theta12 <= 0.99 delta/(n+1), m from 8-21, 8-22 .. 8-26) until theta11, H1 and nu converge;
///    leakage corrections at the row inlet (step 1) and discharge (step 7);
///      B = sum over hub and shroud of 2 pi r rho_e V_me delta1* / mdot [8-76], K_B = 1 - B [8-77].
///
/// Our choices (Aungier does not give them):
///  * the integration of 8-5 .. 8-7 between stations k-1 and k: the fluxes under d/dm differenced, every
///    other term by the trapezoidal rule over the two stations (second order; test_endwall.cpp
///    measures the order), with d/dm of an edge quantity the difference quotient over the interval:
///      delta1* r rho_e V_me dV_me/dm   -> mean of (delta1* r rho_e V_me) times (V_me,k - V_me,k-1),
///      r delta1* rho_e V_me d(r C_theta_e)/dm -> the same mean times (r C_theta_e)_k - (r C_theta_e)_k-1;
///    the station equations are solved by the substitution of steps 2-6, started from 8-74 and 8-75
///    with H1 held (8-23 gives delta - delta1*); H2 is not held: q12 is the marched variable, so 8-75
///    already fixes theta12 and with it m (holding H2 too would overdetermine the start); the start
///    does not change the converged station;
///  * the blade forces of 8-8 and 8-9 are interval quantities (they jump at a leading or trailing
///    edge): f_me and f_theta_e of an interval from the interval means of rho_e, V_me, r, sin(phi),
///    C_theta_e^2 and the difference quotients of V_me, P_e and r C_theta_e; they enter 8-6 and 8-7
///    only through the defect terms r f_me nu and r^2 f_theta_e nu, nu = 0 outside blade rows;
///  * the meridional force of 8-6 is 8-8 (the book's definition); 8-49 is used only where the book
///    uses it, in xi (8-50 .. 8-53).  So 8-57's fixed point is 8-54 up to that difference (measured);
///  * in 8-53, "rbar f_theta_e dm" is the sum over the row's intervals of r f_theta_e dm, rbar =
///    (r_in + r_out)/2, Wbar and Vbar the means of the row's inlet and discharge values, W relative
///    to the row (for a one-step row: the step's two stations and the step length, the q3d contract's
///    reading); xi = 0 (a row without turning) gives dnu = 0;  8-47, 8-48 and 8-50 are not coded
///    (their A and B, B not the blockage, are the nu-free parts used only to derive xi);
///  * |dnu| <= [0.05 + 0.2 (H1 - 1)] delta at the row discharge, Aungier's Eq. 12-50 (2003 p. 306,
///    stated there for a local blade force that "may be extremely small, or even zero") applied to
///    the Ch. 8 row correction: without it dnu ~ 1/xi grows without bound as the turning -> 0 and
///    dnu f_theta_e stays finite (a force defect without a blade force), while xi = 0 gives dnu = 0;
///    with it the defect terms vanish continuously with the blade force;
///  * the row iteration: the row's intervals are marched with nu, then theta11bar = mean of the
///    row-inlet (after the inlet leakage correction) and discharge (before the discharge
///    correction) theta11, nu0 (8-44), F (8-55) and dnu (8-57) are updated, until |nu change| <=
///    tolerance * g;  nu starts from nu0 of the inlet theta11 and dnu = 0;
///  * the caller passes s bar and gamma bar of the row at that wall (8-2's averages of inlet and
///    discharge values); delta_c = 0 for an attached blade; a shrouded wall with delta_c != 0 is
///    refused (its leakage is the seal's, Sec. 8.6; the book is silent on delta_c there);
///  * the limits act in the closures; the marched (delta - delta1*), theta11 and q12 are not reset:
///    H1 limited -> n = 0.7, delta1* = H1 theta11, delta = (delta - delta1*) + delta1* (as BoundaryLayer's
///    delta* = H theta with H_k limited); theta12 limited by the step-4 precaution -> m from the limited
///    value; m limited at (2.4 - 1)/2 and, ours, at m >= -0.85 (H2 >= 0.40, the end of Fig. 8-4's
///    extended curve; `extended_m_min`); 8-41 gets theta11, delta1*, the (precaution-limited) theta12
///    and the profile's delta2* and theta22; H = 2.4 when theta <= 0 after the >= 0 limits of p. 186;
///  * 8-73 (rho_e V_e theta11 / mu >= 250) acts on the state: a station (or the inlet) below it is
///    scaled, (delta - delta1*), theta11 and q12 by the same factor (H1 and m kept);
///  * V_theta_e -> 0 with q12 != 0: theta12 = q12 / V_theta_e -> +-infinity.  Toward the side where
///    theta12 > 0 the H2 limit holds m at 0.7; toward the other, 8-31 sends m -> -infinity and
///    V_theta_e delta2* (8-69's product) to a finite value, so without a second limit the 8-6 bracket
///    jumps at V_theta_e = 0 (0.0569 against 0 in the review's case).  With m >= -0.85 the profile
///    stays finite on both sides, the bracket, tau_theta_w and the streamwise terms (all carrying
///    V_theta_e) tend to 0 from both sides, and the layer is continuous through V_theta_e = 0
///    (tested).  The finite physical limit of 8-69 is not reproduced: the power-law family cannot hold
///    a tangential defect at zero edge velocity, and the book's own H2 limit already drops it on one
///    side.  At V_theta_e = 0 exactly, m is held from the last station where it was defined (the
///    inlet's to start with) and theta12 is reported as NaN when q12 != 0;
///  * leakage: the station records keep the layer outside the row (the uncorrected inlet layer, the
///    corrected discharge layer); the "+" and "-" layers are in EndwallLeakage;  U_leak at the row
///    inlet and discharge from 8-61 with that station's r and edge values (the shroud frame is the
///    row's);  a shrouded row must have its walls' segments rotating at the row speed (refused otherwise);
///  * K_B_geometric = 1 - (delta1*_hub + delta1*_shroud) / quasi-normal length (Ch. 7 p. 157's wording,
///    no correction for the quasi-normal's lean), a diagnostic beside K_B of 8-77.
#pragma once

#include <limits>
#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

/// The end-wall analysis refused its input or failed (no value is clipped or repaired).
class EndwallError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

// ================================ input ================================================

/// One end wall at one station: the wall point and the edge state of the end-wall stream surface.
struct EndwallWallStation
{
    // the wall
    double z = std::numeric_limits<double>::quiet_NaN();    ///< axial coordinate [m] (the record only)
    double r = std::numeric_limits<double>::quiet_NaN();    ///< radius [m], > 0
    double phi = std::numeric_limits<double>::quiet_NaN();  ///< contour angle with the axis [rad], sin(phi) = dr/dm
    double m = std::numeric_limits<double>::quiet_NaN();    ///< arc length along this wall [m], increasing in k
    double wall_omega = 0.0;  ///< omega_w of the wall segment from station k-1 to k [rad/s] (station 0: inlet frame)
    // the edge (absolute frame)
    double Vm = std::numeric_limits<double>::quiet_NaN();      ///< C_m = W_m = V_m [m/s], > 0  [8-3]
    double Ctheta = std::numeric_limits<double>::quiet_NaN();  ///< C_theta [m/s], book sign
    double rho = std::numeric_limits<double>::quiet_NaN();     ///< static density [kg/m^3], > 0
    double P = std::numeric_limits<double>::quiet_NaN();       ///< static pressure [Pa] (8-8, inside blade rows)
    double mu = std::numeric_limits<double>::quiet_NaN();      ///< viscosity [Pa s], > 0 (e.g. sutherland_viscosity(T_e))
};

/// One station of the march: both walls.
struct EndwallStation
{
    EndwallWallStation hub, shroud;
    double quasi_normal_length = std::numeric_limits<double>::quiet_NaN();  ///< hub-to-shroud [m] (K_B_geometric only)
};

/// One wall of a blade row, for a step (or the steps) spanning that row.
struct EndwallRowWall
{
    double omega = 0.0;    ///< blade rotation speed [rad/s], book sign (0 for a stator)
    int blade_count = 0;   ///< Z (not used by Ch. 8; for 12-48)
    double pitch = std::numeric_limits<double>::quiet_NaN();    ///< s bar at this wall [m], > 0 (8-2)
    double stagger = std::numeric_limits<double>::quiet_NaN();  ///< gamma bar at this wall [rad], |gamma| < pi/2 (8-2)
    double clearance = 0.0;     ///< delta_c [m] >= 0: tip / hub gap; 0 where the blade is attached or shrouded
    bool shrouded = false;      ///< a shroud band with a seal at this wall (Sec. 8.6)
    double seal_leakage = 0.0;  ///< mdot_leak [kg/s], > 0 from the row discharge to its inlet (Ch. 6); 0 unless shrouded
};

/// A blade row of the march: its stations and its two walls.
struct EndwallBladeRow
{
    int station_in = -1;         ///< leading-edge station
    int station_out = -1;        ///< trailing-edge station, > station_in
    EndwallRowWall hub, shroud;  ///< the same omega on both
};

/// The inlet layer (Sec. 8.8): blockage B shared equally by the two walls [8-72], H1 and m.
struct EndwallInlet
{
    double blockage = 0.02;  ///< B at station 0 (p. 194: 0.02 just upstream of the first row)
    double H1 = 1.4;         ///< p. 194
    double m = 1.0 / 7.0;    ///< p. 194 (H2 = 1.286)
};

enum class EndwallEntrainment { FreeStreamH, MeridionalH1 };  ///< E of H [8-43] or of H1 [8-42]
[[nodiscard]] const char* endwall_entrainment_name(EndwallEntrainment e) noexcept;  ///< "free_stream_H", "meridional_H1"

/// The blade force defect model inside blade rows.  The march keeps one nu per interval, so a model
/// that varies nu through the passage (Aungier 2003 Ch. 12, Eqs. 12-44 to 12-50, for the quasi-3D
/// in-passage march; 12-46 has unbalanced brackets as printed) is one more case here.  Only Ch. 8's
/// row model is implemented.
enum class EndwallForceDefect { None, RowChapter8 };  ///< nu = 0; Sec. 8.5 (one nu across the row)
[[nodiscard]] const char* endwall_force_defect_name(EndwallForceDefect f) noexcept;  ///< "none", "row_chapter8"

struct EndwallSettings
{
    EndwallEntrainment entrainment = EndwallEntrainment::FreeStreamH;
    EndwallForceDefect force_defect = EndwallForceDefect::RowChapter8;  ///< Sec. 8.5; None: nu = 0 everywhere
    bool force_defect_correction = true;  ///< 8-57; false: nu = nu0
    double shape_factor_limit = 2.4;      ///< H1, H2 (p. 184) and H (p. 186)
    double re_theta_min = 250.0;          ///< 8-36 (p. 186) and 8-73
    double theta12_precaution = 0.99;     ///< theta12 <= 0.99 delta/(n+1) (p. 195 step 4)
    double extended_m_min = -0.85;        ///< ours: m >= this (H2 >= 0.40, the end of Fig. 8-4's curve)
    double tolerance = 1e-12;             ///< ours: relative change of the station state, and of nu/g
    int max_iterations = 200;             ///< ours: per station and per row
};

struct EndwallInput
{
    std::vector<EndwallStation> stations;  ///< >= 2
    std::vector<EndwallBladeRow> rows;     ///< inlet-to-exit order, not overlapping (a row may start where the last ended)
    double mass_flow = std::numeric_limits<double>::quiet_NaN();  ///< compressor mdot [kg/s], > 0 (8-72, 8-76)
    double compressor_omega = 0.0;  ///< omega_c [rad/s], book sign (8-61)
    EndwallInlet inlet;
    EndwallSettings settings;
};

// ================================ output ===============================================

/// The marched, frame-invariant state of one layer (Sec. 8.7).
struct EndwallState
{
    double Q = 0.0;        ///< delta - delta1* [m]
    double theta11 = 0.0;  ///< [m]
    double q12 = 0.0;      ///< V_theta_e theta12 = C_theta_e theta12 = W_theta_e theta12' [m^2/s] (8-68)
};

/// The layer on one wall at one station, in the frame of the wall segment ending there.  It is
/// also what the next step continues from (Q, theta11, q12 and m).
struct EndwallLayer
{
    double Q = 0.0, theta11 = 0.0, q12 = 0.0;        ///< the state
    double wall_omega = 0.0;                          ///< the frame [rad/s]
    double Vtheta_e = 0.0, V_e = 0.0, beta_e = 0.0;   ///< edge relative to the wall [8-4], |V_e|, atan2(V_theta_e, V_me)
    double delta = 0.0, delta1 = 0.0, delta2 = 0.0;   ///< delta, delta1*, delta2*
    double theta12 = 0.0, theta22 = 0.0;              ///< theta12 = q12 / V_theta_e (the profile's when V_theta_e = 0)
    double n = 0.0, m = 0.0, H1 = 0.0, H2 = 0.0;      ///< profile exponents and shape factors (after the limits)
    double delta_star = 0.0, theta = 0.0, H = 0.0;    ///< streamwise [8-40, 8-41], H after its limit
    double Re_theta = 0.0, Re_theta11 = 0.0;          ///< rho_e V_e theta / mu (before its limit), rho_e V_e theta11 / mu
    double cf = 0.0, E = 0.0, tau_mw = 0.0, tau_thetaw = 0.0;  ///< [8-36], [8-42/43], [8-34], [8-35]
    double nu = 0.0, dnu = 0.0;  ///< the force defect of the interval ending here (0 outside rows) [8-56], [8-57]
    bool H1_limited = false, H2_limited = false, theta12_limited = false, H_limited = false;
    bool re_theta_limited = false;    ///< 8-36's Re_theta raised to the limit
    bool re_theta11_limited = false;  ///< 8-73 scaled the state at this station
    int iterations = 0;               ///< substitution passes at this station (0 at the inlet)
};

/// B and K_B of one station [8-76, 8-77].
struct EndwallBlockage
{
    double B_hub = 0.0, B_shroud = 0.0;  ///< each wall's 2 pi r rho_e V_me delta1* / mdot
    double B = 0.0;                      ///< [8-76]
    double K_B = 1.0;                    ///< [8-77], the blockage factor of Eq. 7-5
};

struct EndwallStationResult
{
    EndwallLayer hub, shroud;
    EndwallBlockage blockage;
    double K_B_geometric = std::numeric_limits<double>::quiet_NaN();  ///< ours, diagnostic (NaN without the length)
};

/// The seal-leakage corrections of one shrouded row wall (Sec. 8.6).
struct EndwallLeakage
{
    double delta_plus = 0.0, delta_minus = 0.0;  ///< [8-59], [8-63] [m]
    double U_leak_in = 0.0, U_leak_out = 0.0;    ///< [8-61] at the inlet and discharge stations [m/s]
    EndwallState in, plus, minus, out;          ///< before / after the inlet correction, before / after the discharge one
};

/// The blade force defect of one row wall (Sec. 8.5).
struct EndwallRowWallResult
{
    double g = 0.0;             ///< staggered spacing [8-2] [m]
    double theta11_mean = 0.0;  ///< theta11 bar
    double nu0 = 0.0, dnu = 0.0, nu = 0.0;  ///< [8-44], [8-57], [8-56]; nu is what the reported layers used
    double F = 0.0;                  ///< [8-55]
    double dnu_bound = 0.0;          ///< [0.05 + 0.2 (H1 - 1)] delta at the discharge (12-50, ours here)
    bool dnu_limited = false;        ///< 8-57's dnu exceeded the bound
    double xi = 0.0;                 ///< [8-53] [m/s]
    double collateral_defect = 0.0;  ///< V_theta_e theta12 - W_theta_e theta11 at the discharge [8-45] [m^2/s]
    int iterations = 0;
    bool shrouded = false;
    EndwallLeakage leakage;  ///< zeros unless shrouded
};

struct EndwallRowResult
{
    EndwallRowWallResult hub, shroud;
};

struct EndwallResult
{
    std::vector<EndwallStationResult> stations;
    std::vector<EndwallRowResult> rows;
};

// ================================ the analysis =========================================

/// The inlet layer of one wall at station 0 (Sec. 8.8): delta1* from 8-72, 8-22 .. 8-26, then 8-73.
[[nodiscard]] EndwallLayer endwall_inlet(const EndwallWallStation& wall, const EndwallInlet& inlet, double mass_flow,
                                         const EndwallSettings& settings);

/// The result of one step on one wall.
struct EndwallStepResult
{
    EndwallLayer layer;  ///< the layer at station k (after the discharge leakage correction)
    bool row = false;    ///< the step spanned a blade row
    EndwallRowWallResult row_result;
};

/// One station on one wall (steps 1-7 of p. 194-195): from the layer at k-1 (`prev`, any frame) and
/// the walls at k-1 and k to the layer at k.  `row` non-null: the step spans that blade row (k-1 its
/// leading edge, k its trailing edge): force defect (Sec. 8.5) and, if shrouded, the seal leakage
/// (Sec. 8.6).  Pure: the same arguments give the same result.  EndwallError for inadmissible input
/// or no converged solution.
[[nodiscard]] EndwallStepResult endwall_step(const EndwallLayer& prev, const EndwallWallStation& prev_wall,
                                             const EndwallWallStation& wall, const EndwallRowWall* row,
                                             double compressor_omega, const EndwallSettings& settings);

/// B [8-76] and K_B [8-77] of a station.
[[nodiscard]] EndwallBlockage endwall_blockage(const EndwallLayer& hub, const EndwallWallStation& hub_wall,
                                               const EndwallLayer& shroud, const EndwallWallStation& shroud_wall,
                                               double mass_flow);

/// The convenience march of both walls through all stations (endwall_inlet, then endwall_step; a
/// row with stations inside is marched as one row with one nu, ours).
[[nodiscard]] EndwallResult analyze_endwalls(const EndwallInput& input);

// ================================ the closures (public for the tests) ==================

/// The thicknesses of the profiles 8-18 and 8-19 (8-28 for m < 0.05) of thickness delta.
struct EndwallProfile
{
    double delta = 0.0, n = 0.0, m = 0.0;
    double delta1 = 0.0, theta11 = 0.0, H1 = 0.0;  ///< [8-24], [8-22]
    double delta2 = 0.0, theta22 = 0.0, H2 = 0.0;  ///< [8-26, 8-25] or [8-29, 8-30, 8-32]
    double theta12 = 0.0;                          ///< 8-21 or 8-31 solved for theta12
};
[[nodiscard]] EndwallProfile endwall_profile(double delta, double n, double m);

/// m of theta12 [8-21] (theta12 >= the value of m = 0.05) or [8-31] (below); +infinity when
/// theta12 >= delta/(n+1) (8-21's pole).
[[nodiscard]] double endwall_m_of_theta12(double theta12, double delta, double n) noexcept;

/// H2 of m: 2m + 1 [8-25] for m >= 0.05, 20 / [21 (0.95 - 1.684 m)] [8-32] below.
[[nodiscard]] double endwall_H2(double m) noexcept;

/// delta* [8-40] and theta [8-41], each thickness on the right limited to >= 0 (p. 186).
struct EndwallStreamwise
{
    double delta_star = 0.0, theta = 0.0;
};
[[nodiscard]] EndwallStreamwise endwall_streamwise(double delta1, double delta2, double theta11, double theta12,
                                                   double theta22, double beta) noexcept;

/// c_f of Ludwieg and Tillmann [8-36].
[[nodiscard]] double ludwieg_tillmann_cf(double H, double re_theta) noexcept;

/// The inlet state: delta1* (from 8-72), H1 and m -> [8-22 .. 8-26]; q12 = V_theta_e theta12.
[[nodiscard]] EndwallState endwall_inlet_state(double delta1, double H1, double m, double Vtheta_e);

/// The closures of a state on a wall station in the frame omega_w (steps 4-5 and Sec. 8.4);
/// held_m is m when V_theta_e = 0.  EndwallError unless Q, theta11 > 0 and the edge is admissible.
[[nodiscard]] EndwallLayer endwall_closures(const EndwallState& state, const EndwallWallStation& wall, double omega_w,
                                            double held_m, const EndwallSettings& settings);

/// Staggered spacing g = s cos(gamma) [8-2].
[[nodiscard]] double staggered_spacing(double pitch, double stagger) noexcept;
/// x^3 / (1 + x^3), x = 8 theta11 / g: the shape of 8-44 and 8-55 (and of Figs. 8-5, 8-6 with x = 8 delta* / g).
[[nodiscard]] double force_defect_shape(double x) noexcept;
/// nu0 [8-44].
[[nodiscard]] double base_force_defect(double g, double clearance, double theta11_mean) noexcept;
/// F [8-55].
[[nodiscard]] double force_defect_factor(double g, double theta11_mean) noexcept;

/// Delta = mdot_leak / (2 pi r rho_e V_me) [8-59, 8-63].
[[nodiscard]] double leakage_thickness(double mdot_leak, double r, double rho, double Vm) noexcept;
/// U_leak = r omega_c / 2 + V_theta - C_theta [8-61].
[[nodiscard]] double leakage_velocity(double r, double omega_c, double Vtheta, double Ctheta) noexcept;
/// The row-inlet correction [8-58, 8-60, 8-62].
[[nodiscard]] EndwallState seal_leakage_inlet(const EndwallState& in, double delta_plus, double Vtheta_in,
                                              double U_leak) noexcept;
/// The row-discharge correction [8-64, 8-65, 8-66].
[[nodiscard]] EndwallState seal_leakage_discharge(const EndwallState& minus, double delta_minus, double Vtheta_out,
                                                  double U_leak) noexcept;

}  // namespace BladeGen::Turbolenza
