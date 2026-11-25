/// @file turbolenza/Throughflow.hpp
/// @brief The meridional through-flow analysis of Aungier, Axial-Flow Compressors (ASME Press,
///        2003), Chapter 7 (pp. 153-173): streamline curvature on straight quasi-normals with
///        blade rows between them, in the four model levels of Sec. 7.6.  This header holds the
///        case, the plug-in interfaces (BladeRowModel and RowLossModel for Chapter 6, EndwallModel
///        for Chapter 8) and the solution structs, with their conventions stated here once.
///
/// CONVENTIONS (binding for every implementation of the interfaces and every user of the result)
///
///  * Units: SI throughout (m, s, kg, Pa, K, J/kg, J/(kg K), rad).  Angles are radians; no field
///    of this header is in degrees (Chapter 6's correlations are in degrees: its module converts).
///  * Coordinates (Figs. 7-1, 7-2): z axial (flow toward +z), r radius, theta the polar angle.  The
///    meridional coordinate m runs along a stream surface (Eq. 7-58, m = 0 on the first
///    quasi-normal).  phi is measured from +z toward +r (Eq. 7-2: sin phi = dr/dm); kappa_m =
///    -dphi/dm (Eq. 7-9 = 3-26), so a surface turning toward +r has kappa_m < 0.  The normal n is
///    the meridional direction turned by +90 deg.
///  * Quasi-normals (Eq. 7-1): straight, from the hub point (z_h, r_h) to the shroud point (z_s, r_s);
///    tan(lambda) = (z_h - z_s) / (r_s - r_h): lambda = 0 radial, lambda > 0 leans upstream at the
///    shroud.  y is the distance from the hub, 0 <= y <= y_s, d(z, r)/dy = (-sin lambda, cos lambda);
///    epsilon = phi - lambda (Eq. 7-3).
///  * Stream surfaces j = 0 (hub contour) .. n_ss - 1 (shroud contour); surface j carries the mass
///    fraction psi_j = j / (n_ss - 1) of the K_B-free mass function of Eq. 7-54 (equal stream tubes,
///    "typically", p. 160).  The end-wall surfaces stay on the contours (p. 164).
///  * Rotation and signs, Aungier's frame: omega [rad/s] is a row's angular speed toward +theta;
///    W_theta = C_theta - omega r (Eqs. 3-1, 7-44); tan(alpha) = C_theta / C_m, tan(beta') = W_theta /
///    W_m (p. 158), both signed toward +theta.  BladeGen's and TURBOlenza's theta IS this theta, and
///    omega = -2 pi rpm / 60: a positive rpm (a BladeGen compressor rotor) has omega < 0 and, with
///    axial inflow, beta'_1 > 0.  The equations of Chs. 3, 7 and 8 are used unchanged with a signed
///    omega (they are invariant under theta -> -theta with omega, C_theta, W_theta and the angles).
///    Metal angles kappa_1, kappa_2 and the stagger gamma (BladeSectionGeometry) are in the same frame.
///    Chapter 6's correlations hold for positive camber only (p. 128): each row carries s_row = +1 or
///    -1, and the Chapter 6 module multiplies the signed angles and swirl velocities by it on entry
///    and its returned beta'_2 by it on exit.
///  * Thermodynamics (Ch. 2): a thermally and calorically perfect gas (Eqs. 2-22, 2-25, 2-28; ours),
///    h = cp T, s = cp ln(T / T_ref) - R ln(P / P_ref).  H total enthalpy, I = H - omega r C_theta the
///    rothalpy (Eqs. 3-13, 7-41), H' = I + (r omega)^2 / 2 (Eq. 7-37), h = H' - W^2 / 2 (Eq. 7-38).
///    The viscosity follows Sutherland's law at the local static temperature (the project's choice;
///    ThroughflowGas).
///  * Loss coefficient (Eq. 7-49, the form of Eq. 2-68 with the inlet velocity pressure): P't2 =
///    P't2,id - omegabar (P't1 - P1), P't2,id the relative total pressure at (H'2, s1), P't1 and P1 the
///    relative total and static pressure on the same stream surface at the inlet quasi-normal.  Every
///    loss handed to the through-flow (BladeRowDischarge, RowLossResult) is in, or converted to, this
///    form at the interface.
///  * Two blockages, two names: K_B_endwall, one number per quasi-normal OUTSIDE the integral of
///    Eq. 7-5 (Ch. 8 Eq. 8-77: K_B = 1 - B), and the pointwise factor INSIDE it: Ch. 7's K_W =
///    1 - B_wake (Eq. 7-6; Ch. 6 Eqs. 6-83, 6-84) times Ch. 12's blade-metal blockage K_B_blade
///    (Eq. 12-7; 1 in this version).  K_W applies on the discharge quasi-normal of its row only, 1
///    elsewhere (ours: the book does not say how far a wake reaches).  mdot = 2 pi K_B_endwall
///    int_0^ys K_W K_B_blade r rho W_m cos(epsilon) dy (Eq. 7-5).
///  * Stream-sheet thickness (Ch. 7 defines none).  The sheets are the tubes between neighbouring
///    surfaces, sheet (j, j + 1) carrying dpsi = 1 / (n_ss - 1).  Per SHEET at its mean surface (the
///    y-midpoint, ours): b_sheet = dmdot_sheet / (2 pi r rho W_m) with the actual sheet mass flow
///    dmdot_sheet = mdot / (n_ss - 1) and the mean of the two surfaces' r rho W_m, beside its geometric
///    thickness cos(epsilon) (y_(j+1) - y_j).  Per SURFACE j, the same two thicknesses of a sheet of
///    dpsi centred on the surface: b_j = cos(epsilon_j) D_j, D_j the sheets' y-width interpolated
///    (second order) to the surface, the area term of Eq. 7-32; and b_mass_j = dpsi mdot / (2 pi r
///    K_B_endwall K_W K_B_blade rho W_m) (Eq. 12-39 with Ch. 7's blockages).  Values are given at the
///    quasi-normals; inside a blade row the consumer interpolates linearly in m between the two
///    quasi-normals bounding it (Aungier p. 287).
///
/// CALL ORDER (Figs. 7-3, 7-4).  The quasi-normals are solved in order i = 0 .. n_qn - 1 in every
/// pass.  On quasi-normal i the flow field is iterated (Sec. 7.2 with the linking of Sec. 7.3).  In
/// every iteration, for a blade row between quasi-normals i - 1 and i: BladeRowModel::evaluate once
/// per stream surface with the latest discharge estimate, then RowLossModel::evaluate once (if the
/// row has one); after the flow field of the iteration, EndwallModel::evaluate once (Fig. 7-3,
/// "Compute end-wall boundary layer data", inside the outer loop and before the repositioning).
/// Repeated calls for the same quasi-normal replace the earlier ones.  EndwallModel::begin_pass
/// starts every pass.  The full normal equilibrium model makes many passes (the outer loop); the
/// simplified models one.
///
/// NOT IN THIS VERSION, kept open for Chapter 12 (pp. 287-307):
///  * quasi-normals inside blade rows (QuasiNormal::inside_blade_row, BladeRowSpec::qn_outlet >
///    qn_inlet + 1): there Eqs. 7-19 to 7-21 are invalid (they assume Eq. 7-7) and Ch. 12's 12-17 to
///    12-19 apply with an imposed tan(beta')(m) per surface; validate() refuses them.  The quasi-normal
///    solver is written for a tangential mode per quasi-normal (C_theta known: Eqs. 7-16 to 7-18;
///    beta' known outside a row: 7-19 to 7-21; absolute alpha known at the inlet: the same with
///    omega = 0, p. 162), the pointwise blockage K_B_blade inside the integral, and a choke test
///    choice (ThroughflowOptions::choke_test: Eq. 7-39; Ch. 12's Eq. 12-37 is declared but refused in this
///    version, see ChokeTest);
///  * the damping constants are parameters (B*, the Mach cap of Eq. 7-57 / 12-34, the curvature
///    damping D of Eq. 12-38); the in-passage damping Eq. 12-35 is not implemented;
///  * annulus sizing (Sec. 7.7) and a real-gas equation of state (Ch. 2).
#pragma once

#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza::Throughflow {

/// Refusal of a through-flow case or state (invalid input, non-physical state, reverse flow).
class ThroughflowError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

// ─────────────────────────────────────────────────────────────────────────────
// Gas
// ─────────────────────────────────────────────────────────────────────────────

/// Thermally and calorically perfect gas (Aungier Eqs. 2-22, 2-25, 2-28; Aungier himself warns
/// against limiting an analysis to it, p. 25: ours for this version).  The enthalpy reference is
/// h = cp T; the entropy reference s = 0 at (T_ref, P_ref).  The viscosity: Sutherland's law
/// mu = mu_ref (T / T_mu)^1.5 (T_mu + S) / (T + S) (BoundaryLayer.hpp sutherland_viscosity, air:
/// 1.716e-5 Pa s, 273.15 K, 110.4 K) at the local static temperature, or a constant.
struct ThroughflowGas
{
    double gamma = 1.4;
    double gas_constant = 287.05;      ///< R [J/(kg K)]
    double t_ref = 288.15, p_ref = 101325.0;
    bool sutherland = true;            ///< false: the constant viscosity below
    double mu_ref = 1.716e-5, t_mu = 273.15, s_mu = 110.4;  ///< Sutherland's constants
    double viscosity = 1.8e-5;         ///< constant mu [Pa s] when !sutherland

    void validate() const;  ///< throws ThroughflowError
    [[nodiscard]] double cp() const noexcept { return gamma * gas_constant / (gamma - 1.0); }
    /// Eq. 2-25.
    [[nodiscard]] double entropy(double T, double P) const;
    /// The pressure at (h, s): T = h / cp, P = P_ref exp((cp ln(T / T_ref) - s) / R) (Eq. 2-25
    /// solved for P).  Throws ThroughflowError when h <= 0.
    [[nodiscard]] double pressure(double h, double s) const;
    /// mu(T) [Pa s].  Throws ThroughflowError for T <= 0.
    [[nodiscard]] double mu(double T) const;
};

/// The static state at (h, s).
struct GasState
{
    double T = 0.0, P = 0.0, rho = 0.0, a = 0.0;
};
/// T = h / cp, P from (h, s), rho = P / (R T), a = sqrt(gamma R T) (Eq. 2-28 = 2-26 for the
/// perfect gas).  Throws ThroughflowError when h <= 0 or not finite.
[[nodiscard]] GasState gas_state(const ThroughflowGas& gas, double h, double s);

// ─────────────────────────────────────────────────────────────────────────────
// Blade rows (Sec. 7.3; Chapter 6 plugs in through BladeRowModel and RowLossModel)
// ─────────────────────────────────────────────────────────────────────────────

/// The blade geometry of one section (Exercise 7.1: chord, t_b / c, camber, stagger, kappa_1,
/// kappa_2, a / c; and the throat o of Eq. 6-53), in the book's frame (not yet multiplied by s_row).
struct BladeSectionGeometry
{
    double span = 0.0;               ///< row span fraction (0 hub, 1 shroud; see BladeRowSpec)
    double chord = 0.0;              ///< c [m]
    double stagger = 0.0;            ///< gamma [rad]
    double camber = 0.0;             ///< theta = kappa_1 - kappa_2 [rad]
    double kappa1 = 0.0;             ///< inlet metal angle [rad]
    double kappa2 = 0.0;             ///< exit metal angle [rad]
    double thickness_ratio = 0.0;    ///< t_b / c
    double max_camber_position = 0.5;    ///< a / c
    double max_thickness_position = 0.5; ///< x of t_b / c over c
    double throat = 0.0;             ///< o [m] (Eq. 6-53); 0: not known
};

/// What the through-flow hands a blade-row model for one stream surface j of one row in one
/// iteration of the discharge quasi-normal (Sec. 7.3, p. 163: the inlet relative conditions from the
/// known upstream flow and the current estimate of the discharge, which the empirical models need).
/// Inlet = the quasi-normal upstream of the row (subscript 1), discharge = the one downstream
/// (subscript 2), both outside the blade.  Relative quantities are in the row's frame (omega).
struct BladeRowInlet
{
    int row = 0;                     ///< index of the row in ThroughflowCase::rows
    int row_number = 1;              ///< N_row of Eq. 6-87: 1, 2, ... in the flow direction
    int stream_surface = 0;          ///< j
    int n_stream_surfaces = 0;
    double mass_fraction = 0.0;      ///< psi_j
    double span = 0.0;               ///< row span fraction of this stream surface (BladeRowSpec)
    double omega = 0.0;              ///< row speed [rad/s] (book sign)
    int n_blades = 0;
    int s_row = 1;                   ///< the Chapter 6 cascade sign (BladeRowSpec::s_row)
    bool is_igv = false;
    std::string profile_family;
    double blade_height = 0.0;       ///< mean of the LE and TE line lengths [m] (ours; Eqs. 6-99, 6-100)
    double mass_flow = 0.0;          ///< the case's mdot
    // stream surface at the two quasi-normals
    double z1 = 0.0, r1 = 0.0, z2 = 0.0, r2 = 0.0;
    double pitch = 0.0;              ///< 2 pi r_mean / n_blades, r_mean = (r1 + r2) / 2 (ours)
    double solidity = 0.0;           ///< chord / pitch
    BladeSectionGeometry geometry;   ///< interpolated to `span` (BladeRowSpec), book signs
    // inlet (known, upstream quasi-normal)
    double beta1 = 0.0;              ///< beta'_1 = atan2(W_theta1, W_m1) [rad]
    double alpha1 = 0.0;             ///< absolute inlet angle atan2(C_theta1, C_m1) [rad]
    double W1 = 0.0, Wm1 = 0.0, Wtheta1 = 0.0, Ctheta1 = 0.0;
    double mach1_rel = 0.0;          ///< M'_1 = W1 / a1
    double mach1_abs = 0.0;          ///< C1 / a1
    double reynolds = 0.0;           ///< Re_c = rho1 W1 c / mu(T1) (ours: the book names Re_c only)
    double rho1 = 0.0, P1 = 0.0, T1 = 0.0, h1 = 0.0, a1 = 0.0, mu1 = 0.0;
    double Hrel1 = 0.0, Pt1_rel = 0.0, Tt1_rel = 0.0;  ///< H'_1, P'_t1, T'_t1 (Eq. 7-46)
    double s1 = 0.0, I1 = 0.0;
    // discharge estimate (iterated; W_m2 = W_m1 on the first call, p. 163)
    double Wm2 = 0.0;
    double rho2 = 0.0;               ///< density of the current discharge estimate (rho1 on the first call)
};

/// What a blade-row model returns for one stream surface.
struct BladeRowDischarge
{
    double beta2 = 0.0;              ///< discharge relative flow angle beta'_2 [rad], book sign
    double loss_coefficient = 0.0;   ///< omegabar of Eq. 7-49 (CONVENTIONS)
    double wake_blockage = 1.0;      ///< K_W = 1 - B_wake (Eq. 7-6) on the discharge quasi-normal
};

/// The empirical blade-row model per stream surface (Chapter 6 implements it).  evaluate is called
/// once per stream surface per iteration of the discharge quasi-normal (CALL ORDER); it must be a
/// pure function of its input.
class BladeRowModel
{
public:
    virtual ~BladeRowModel() = default;
    [[nodiscard]] virtual BladeRowDischarge evaluate(const BladeRowInlet& in) const = 0;
};

/// A specified blade row for tests and design studies: beta'_2, omegabar and K_W as functions of
/// the row span fraction (linear interpolation; the span points must cover [0, 1]).
class SpecifiedBladeRow final : public BladeRowModel
{
public:
    /// Throws ThroughflowError unless span is strictly increasing from 0 to 1 and the profiles
    /// have its length (K_W may be empty: 1 everywhere), with 0 < K_W <= 1 and |beta2| < pi / 2.
    SpecifiedBladeRow(std::vector<double> span, std::vector<double> beta2, std::vector<double> loss,
                      std::vector<double> wake_blockage = {});
    [[nodiscard]] BladeRowDischarge evaluate(const BladeRowInlet& in) const override;

private:
    std::vector<double> span_, beta2_, loss_, kw_;
};

/// The state of one end-wall stream surface at one quasi-normal (the wall stream surfaces are the
/// boundary-layer edges).
struct WallState
{
    double z = 0.0, r = 0.0;
    double m = 0.0;                  ///< meridional distance along the wall from quasi-normal 0 (Eq. 7-58)
    double phi = 0.0;                ///< the wall's geometric angle (Eq. 7-2), every model level
    double Cm = 0.0;                 ///< = W_m = V_me (Eq. 8-3)
    double Ctheta = 0.0;             ///< absolute
    double Wtheta = 0.0;             ///< relative to the frame named where the state is handed over
    double rho = 0.0, P = 0.0, T = 0.0, a = 0.0, mu = 0.0;
    double Pt = 0.0, Tt = 0.0;       ///< absolute total state
};

/// What the through-flow hands a row-level loss model (Ch. 6 Secs. 6.9, 6.10, Eqs. 6-85 to 6-96): the
/// end-wall states at both quasi-normals of the row, W_theta in the ROW's frame at both (the
/// discharge at the current estimate, computed with the per-surface losses only), the mass flow, the
/// row data.
struct RowLossInput
{
    int row = 0;
    int row_number = 1;              ///< N_row of Eq. 6-87
    double omega = 0.0;
    int n_blades = 0;
    int s_row = 1;
    bool shrouded = false;
    double hub_clearance = 0.0, tip_clearance = 0.0;
    double mass_flow = 0.0;
    double blade_height = 0.0;
    BladeSectionGeometry hub_section, shroud_section;
    WallState hub1, shroud1, hub2, shroud2;  ///< 1: inlet quasi-normal, 2: discharge quasi-normal
};

/// Row-level total-pressure losses (Ch. 6 Eqs. 6-89, 6-90) and the seal leakage (6-91).
struct RowLossResult
{
    enum class Side
    {
        Hub,
        Shroud,
    };
    double clearance_dPt = 0.0;      ///< Delta P_t of the clearance [Pa], whole row
    Side clearance_side = Side::Shroud;
    double seal_dPt = 0.0;           ///< Delta P_t of the shroud seal [Pa], whole row
    Side seal_side = Side::Hub;
    double leakage_mass_flow = 0.0;  ///< mdot_leak [kg/s], + from discharge to inlet (forwarded to Ch. 8)
};

/// The row-level loss model (Chapter 6 Secs. 6.9, 6.10 implement it).  Its Delta P_t are spread over
/// the stream surfaces by distribute_row_loss and added to omegabar of Eq. 7-49 as Delta P_t,j /
/// (P't1 - P1)_j (ours: the book does not say how the distributed Delta P_t enters Eq. 7-49).
class RowLossModel
{
public:
    virtual ~RowLossModel() = default;
    [[nodiscard]] virtual RowLossResult evaluate(const RowLossInput& in) const = 0;
};

/// The linear distribution of a row-level Delta P_t over the mass fractions psi_j (p. 147: "a linear
/// distribution across the annulus", zero at the wall opposite the gap; ours: linear in psi, its
/// mass-weighted mean, int_0^1 dPt_j dpsi, equals dPt): dPt_j = 2 dPt psi_j with the gap at the
/// shroud, 2 dPt (1 - psi_j) at the hub.
[[nodiscard]] std::vector<double> distribute_row_loss(double dPt, RowLossResult::Side side,
                                                      const std::vector<double>& psi);

/// One blade row between quasi-normals qn_inlet and qn_inlet + 1.  The section geometry is given
/// at row span fractions; a stream surface's span fraction is the mean of its fractions along the
/// LE and TE lines, where the straight segment between its points on the two quasi-normals
/// crosses them (ours, Exercise 7.1; the wall surfaces are 0 and 1), and its geometry the linear
/// interpolation in span.
struct BladeRowSpec
{
    std::string name = "row";
    int qn_inlet = 0;
    /// The discharge quasi-normal; -1: qn_inlet + 1.  This version (Ch. 7) requires qn_inlet + 1:
    /// quasi-normals inside a blade row are Ch. 12's (NOT IN THIS VERSION, see the file comment).
    int qn_outlet = -1;
    double omega = 0.0;              ///< [rad/s], book sign, 0 for a stator
    int n_blades = 0;
    int s_row = 1;                   ///< Ch. 6 cascade sign, +1 or -1 (CONVENTIONS)
    bool is_igv = false;
    bool shrouded = false;
    std::string profile_family = "NACA65";
    double hub_clearance = 0.0;      ///< [m] (for the Ch. 6 row-level and Ch. 8 models)
    double tip_clearance = 0.0;      ///< [m]
    double le_hub_z = 0.0, le_hub_r = 0.0, le_shroud_z = 0.0, le_shroud_r = 0.0;
    double te_hub_z = 0.0, te_hub_r = 0.0, te_shroud_z = 0.0, te_shroud_r = 0.0;
    std::vector<BladeSectionGeometry> sections;  ///< spans strictly increasing, covering [0, 1]
    std::shared_ptr<const BladeRowModel> model;
    std::shared_ptr<const RowLossModel> row_loss;  ///< null: none
};

// ─────────────────────────────────────────────────────────────────────────────
// End walls (Sec. 7.4, Fig. 7-3; Chapter 8 plugs in through EndwallModel)
// ─────────────────────────────────────────────────────────────────────────────

/// The edge state of one wall at one quasi-normal (W_theta in the frame of the quasi-normal).
struct EndwallEdge : WallState
{
    /// Rotation speed of the wall segment between the previous quasi-normal and this one
    /// (QuasiNormal::hub_wall_omega / shroud_wall_omega).
    double wall_omega = 0.0;
    double Vtheta = 0.0;             ///< C_theta - r wall_omega, the wall frame (Eq. 8-4)
};

/// A blade row's data at one wall, for a boundary-layer step across the row.
struct EndwallRowWall
{
    BladeSectionGeometry section;    ///< the section at this wall (span 0 or 1)
    double pitch_mean = 0.0;         ///< s-bar: 2 pi (r_1 + r_2) / 2 / Z on the wall
    double stagger_mean = 0.0;       ///< gamma-bar: gamma at this wall
    double clearance = 0.0;          ///< delta_c at this wall
    double r_inlet = 0.0, r_discharge = 0.0;
    double Wtheta_inlet = 0.0, Wtheta_discharge = 0.0;  ///< in the row frame (Eq. 8-67)
    double delta_m = 0.0;            ///< wall arc length from the inlet to the discharge quasi-normal
    double delta_rCtheta = 0.0;      ///< discharge minus inlet r C_theta on the wall surface
};

/// A blade row bounding a quasi-normal.  For the row upstream (the step that spans it ends here)
/// all fields are set; for the row downstream only the row data and the inlet fields.
struct EndwallRow
{
    int row = -1;                    ///< index in ThroughflowCase::rows, -1: none
    int row_number = 0;              ///< N_row (1, 2, ... in the flow direction)
    double omega = 0.0;
    int n_blades = 0;
    bool shrouded = false;
    double leakage_mass_flow = 0.0;  ///< mdot_leak from the row-level model (Eq. 6-91), + discharge -> inlet
    EndwallRowWall hub, shroud;
};

/// What the through-flow hands an end-wall model on one quasi-normal.
struct EndwallStation
{
    int qn = 0;
    int n_qn = 0;
    int pass = 0;                    ///< outer pass (0 for the simplified models)
    double mass_flow = 0.0;          ///< the case's mdot (Ch. 8 Eq. 8-76)
    double compressor_omega = 0.0;   ///< omega_c of Eq. 8-61 (ThroughflowCase::compressor_omega)
    double y_s = 0.0;                ///< quasi-normal length
    double lambda = 0.0;
    double frame_omega = 0.0;        ///< omega of the frame the quasi-normal was solved in
    int row_upstream = -1;           ///< blade row between qn - 1 and qn (-1: none)
    int row_downstream = -1;         ///< blade row between qn and qn + 1 (-1: none)
    EndwallRow upstream, downstream; ///< the data of those rows (row = -1 when none)
    EndwallEdge hub, shroud;
};

struct EndwallResult
{
    double blockage_factor = 1.0;    ///< K_B_endwall of the quasi-normal (Eq. 7-5), 0 < K_B <= 1
};

/// The end-wall blockage model (Chapter 8 implements it, marching along both walls).  K_B is an area
/// correction only (p. 164): the end-wall stream surfaces stay on the contours; no blade-force defect
/// and no end-wall loss enters the through-flow (Ch. 8 p. 175: no end-wall boundary-layer model is
/// accurate enough for losses; p. 195: K_B "for use in the meridional through-flow analysis").  It may
/// keep state between calls (the march) but must give the same K_B for repeated calls on one
/// quasi-normal with the same input (Fig. 7-3 repeats them).
class EndwallModel
{
public:
    virtual ~EndwallModel() = default;
    /// Called before quasi-normal 0 of every pass.
    virtual void begin_pass(int pass, int n_qn) { (void)pass, (void)n_qn; }
    [[nodiscard]] virtual EndwallResult evaluate(const EndwallStation& st) = 0;
};

/// K_B constant on every quasi-normal (default 1: no blockage).
class ConstantEndwall final : public EndwallModel
{
public:
    explicit ConstantEndwall(double blockage_factor = 1.0);  ///< throws unless 0 < K_B <= 1
    [[nodiscard]] EndwallResult evaluate(const EndwallStation& st) override;

private:
    double kb_ = 1.0;
};

// ─────────────────────────────────────────────────────────────────────────────
// The case
// ─────────────────────────────────────────────────────────────────────────────

/// Sec. 7.6, the four aerodynamic models.  phi and epsilon of the mass equations (Eqs. 7-5, 7-32,
/// 7-54) are always the geometric ones; the model level sets what the normal-momentum equation uses.
enum class ThroughflowModel
{
    FullNormalEquilibrium,   ///< Sec. 7.5: phi, kappa_m, epsilon from the stream surfaces; outer loop with Novak damping
    Approximate,             ///< phi and kappa_m linear in y between the end-wall values (p. 168)
    /// epsilon = kappa_m = 0 in Eqs. 7-12 to 7-19 (p. 168) and in Eq. 7-29 (ours: the book's list stops
    /// at 7-19); the stream surfaces' geometric phi from Eq. 7-70 (p. 172) for the mass equations
    SimpleNonIsentropic,
    SimpleIsentropic,        ///< also the entropy gradient along the quasi-normal dropped (p. 168)
};
[[nodiscard]] std::string to_string(ThroughflowModel m);

/// The choke test on a quasi-normal while W_m(0) is sought.
enum class ChokeTest
{
    MeanMeridionalMach,      ///< Eq. 7-39: (1 / y_s) int M_m dy >= 1 (Chapter 7)
    /// Eq. 12-37: d mdot / d W_m(0) <= 0 (Chapter 12).  REFUSED in this version (validate()): on a
    /// strongly non-uniform profile (e.g. a 45-degree cone near choke) mdot(W_m(0)) still rises when one
    /// node reaches the limiting speed (h -> 0, M_m -> infinity), so its maximum over valid states is a
    /// degenerate state, not a choke.  Ch. 12 needs it with a bounded definition of its own.
    MassFlowDerivative,
};

/// A straight quasi-normal (Eq. 7-1) and the end-wall slope and curvature at its two ends (the
/// Approximate model).  A wall value not given is computed from the quasi-normals' wall points:
/// m by Eq. 7-58, phi by Eq. 7-2 and kappa_m by Eq. 7-9 with Eq. 7-67 inside and Eq. 7-70 at the
/// ends (ours; exact for straight walls, second order otherwise).
struct QuasiNormal
{
    double z_hub = 0.0, r_hub = 0.0, z_shroud = 0.0, r_shroud = 0.0;
    std::optional<double> phi_hub, kappa_hub, phi_shroud, kappa_shroud;
    /// Rotation speed of the hub / shroud wall segment between the previous quasi-normal and this one
    /// [rad/s] (for the end-wall model; Ch. 7 itself does not use it).
    double hub_wall_omega = 0.0, shroud_wall_omega = 0.0;
    /// Inside a blade row (Ch. 12).  Chapter 7 places every quasi-normal outside the rows (p. 155):
    /// validate() refuses true in this version.
    bool inside_blade_row = false;
};

/// Inlet boundary conditions on quasi-normal 0 (p. 162), given at fractions y / y_s of its length
/// (ours) and interpolated linearly to the stream surfaces.
struct InletProfile
{
    enum class Swirl
    {
        Ctheta,  ///< swirl_value = C_theta [m/s]
        Alpha,   ///< swirl_value = alpha [rad], tan(alpha) = C_theta / C_m
    };
    std::vector<double> span;        ///< strictly increasing from 0 to 1
    std::vector<double> Pt, Tt;      ///< absolute total pressure [Pa] and temperature [K]
    Swirl swirl = Swirl::Ctheta;
    std::vector<double> swirl_value;

    /// Uniform Pt, Tt and swirl.
    [[nodiscard]] static InletProfile uniform(double Pt, double Tt, Swirl swirl = Swirl::Ctheta, double value = 0.0);
};

struct ThroughflowOptions
{
    int n_stream_surfaces = 11;
    double bstar = 8.0;                 ///< B* of Eq. 7-57 ("typically about 8", p. 165; Ch. 12: 8 to 16)
    bool adjust_bstar = true;           ///< Eq. 7-59 when the position error grows (p. 167)
    /// Eq. 7-59 is applied when the position error exceeds (1 + this) times the smallest error since
    /// the last change of B* (ours: "if that error is increasing", p. 167, read as Ch. 12's "grow
    /// markedly", p. 299; BstarControl).
    double bstar_growth_threshold = 0.05;
    double max_mach_damping = 0.95;     ///< the M_m <= 0.95 constraint on Eq. 7-57 (p. 165; Ch. 12 uses 0.9)
    /// D of Ch. 12's direct curvature damping, Eq. 12-38: kappa_m -> (kappa_m,i + D kappa_m,i-1) / (1 + D)
    /// over the outer iterations i (full model).  0: Chapter 7 (no curvature damping, the default).
    double curvature_damping = 0.0;
    double min_one_minus_mm2 = 0.1;     ///< the constraint of Eq. 7-30
    ChokeTest choke_test = ChokeTest::MeanMeridionalMach;
    /// f4 of Eq. 7-34 by the stream-tube form Eqs. 7-31 to 7-33, 1 / W_m = rho dA / dmdot (true:
    /// Aungier's, p. 160, which avoids the singularity at W_m = 0), or f2 + f3 / W_m of Eq. 7-15 as it
    /// stands (false).  The tubes are the sheets between the surfaces; their rho dA / dmdot is carried
    /// to the surfaces at second order (ours, see the .cpp).
    bool stream_tube_form = true;
    double flow_tolerance = 1e-10;      ///< quasi-normal iteration: max |dW_m| / max |W_m| (ours)
    double position_tolerance = 1e-7;   ///< outer loop: max |y_new - y| / y_s (ours)
    double mass_tolerance = 1e-12;      ///< |mdot(W_m0) - mdot| / mdot in the W_m(0) iteration (ours)
    int max_flow_iterations = 400;
    int max_outer_iterations = 400;
    /// Repositionings before a choke is taken as real: the full model ends on a choke only after
    /// this many (Fig. 7-3 shows 1; "it may be desirable to require more", p. 166).
    int choke_repositionings = 1;
};

struct ThroughflowCase
{
    ThroughflowGas gas;
    ThroughflowModel model = ThroughflowModel::FullNormalEquilibrium;
    double mass_flow = 0.0;              ///< [kg/s]
    double compressor_omega = 0.0;       ///< omega_c, the compressor's speed (Eq. 8-61; for the end-wall model)
    std::vector<QuasiNormal> quasi_normals;  ///< upstream to downstream, at least 3
    InletProfile inlet;
    std::vector<BladeRowSpec> rows;      ///< in the flow direction, at most one between two quasi-normals
    std::shared_ptr<EndwallModel> endwall;   ///< null: ConstantEndwall(1)
    ThroughflowOptions options;

    /// Throws ThroughflowError with the reason.
    void validate() const;
};

// ─────────────────────────────────────────────────────────────────────────────
// The solution
// ─────────────────────────────────────────────────────────────────────────────

/// One stream surface on one quasi-normal.
struct ThroughflowNode
{
    double z = 0.0, r = 0.0, y = 0.0, m = 0.0;
    double mass_fraction = 0.0;      ///< psi_j
    double phi = 0.0, epsilon = 0.0; ///< geometric (the mass equations)
    double kappa_m = 0.0;            ///< as used by the normal-momentum equation of the model level
    double Wm = 0.0;                 ///< = C_m
    double Ctheta = 0.0, Wtheta = 0.0;               ///< W_theta in the frame of the quasi-normal
    double C = 0.0, W = 0.0;
    double alpha = 0.0;              ///< atan2(C_theta, C_m)
    double beta = 0.0;               ///< beta' = atan2(W_theta, W_m) in the frame of the quasi-normal
    double mach = 0.0, mach_rel = 0.0, mach_m = 0.0, mach_theta = 0.0;  ///< C/a, W/a, W_m/a, C_theta/a (Eq. 7-28)
    double Pt = 0.0, Tt = 0.0;       ///< absolute total
    double Pt_rel = 0.0, Tt_rel = 0.0;  ///< total in the frame of the quasi-normal
    double P = 0.0, T = 0.0, rho = 0.0, a = 0.0, h = 0.0, mu = 0.0;
    double s = 0.0;                  ///< entropy (ThroughflowGas reference)
    double H = 0.0, I = 0.0, Hrel = 0.0;  ///< total enthalpy; rothalpy and H' in the frame of the quasi-normal
    double rCtheta = 0.0;
    double K_W = 1.0;                ///< wake blockage (inside the integral)
    double K_B_blade = 1.0;          ///< blade-metal blockage (Ch. 12, Eq. 12-7): 1 in this version
    double b = 0.0;                  ///< stream-sheet thickness at the surface, the area term of Eq. 7-32 (CONVENTIONS)
    double b_mass = 0.0;             ///< Eq. 12-39 with Ch. 7's blockages (CONVENTIONS)
    double loss_coefficient = 0.0;   ///< omegabar of the row upstream (profile + row-level; 0 when none)
    double row_loss_coefficient = 0.0;   ///< its row-level part (RowLossModel)
};

/// One stream sheet (j, j + 1) at its mean surface on one quasi-normal (CONVENTIONS).
struct ThroughflowSheetNode
{
    double z = 0.0, r = 0.0;
    double b_geometric = 0.0;        ///< cos(epsilon) (y_(j+1) - y_j)
    double b = 0.0;                  ///< dmdot_sheet / (2 pi r rho W_m) with the actual sheet mass flow
};

struct ThroughflowStation
{
    int index = 0;
    double lambda = 0.0, y_s = 0.0;
    double frame_omega = 0.0;
    int row_upstream = -1;
    double K_B_endwall = 1.0;
    double mass_flow = 0.0;          ///< Eq. 7-5 on the solution
    double mass_error = 0.0;         ///< (mass_flow - case mdot) / case mdot
    double mean_mach_m = 0.0;        ///< Eq. 7-39
    bool choked = false;
    bool mach_limit_active = false;  ///< the Eq. 7-30 limit was applied on some node
    /// Nodes with M_m > 1 (by more than 1e-6): outside the validity of Eq. 7-29 (the Eq. 7-30 limit
    /// bounds the term only).  A choke on such a profile says so in ThroughflowResult::message.
    int supersonic_nodes = 0;
    bool converged = false;
    int iterations = 0;
    std::vector<ThroughflowNode> nodes;
    std::vector<ThroughflowSheetNode> sheets;  ///< n_ss - 1
};

/// One stream surface through all quasi-normals: what the blade-to-blade solver takes.
struct ThroughflowStreamSurface
{
    double mass_fraction = 0.0;
    std::vector<double> z, r, m, b, b_mass;  ///< b: Eq. 7-32 area term; b_mass: Eq. 12-39
};

/// One stream sheet (j, j + 1) through all quasi-normals.
struct ThroughflowStreamSheet
{
    double mass_fraction = 0.0;      ///< of its mean, (psi_j + psi_(j+1)) / 2
    std::vector<double> z, r, m, b, b_geometric;
};

/// One outer iteration (full model) or the single pass (simplified models).
struct ThroughflowIteration
{
    int pass = 0;
    double position_error = 0.0;     ///< max over all nodes of |y_calc - y| / y_s (Sec. 7.4, before damping)
    double max_mass_error = 0.0;     ///< max |mass_error| over the quasi-normals
    double bstar = 0.0;              ///< B* used in this pass's repositioning
    double min_damping = 1.0;        ///< smallest F of Eq. 7-57 applied
    int flow_iterations = 0;         ///< sum over the quasi-normals
    bool choked = false;
};

struct ThroughflowResult
{
    enum class Status
    {
        Converged,
        Choked,         ///< ended on a choke (Figs. 7-3, 7-4)
        NotConverged,   ///< an iteration limit was reached
    };
    Status status = Status::NotConverged;
    ThroughflowModel model = ThroughflowModel::FullNormalEquilibrium;
    int choked_station = -1;
    std::string message;
    std::vector<ThroughflowStation> stations;          ///< the quasi-normals solved in the last pass
    std::vector<ThroughflowStreamSurface> stream_surfaces;
    std::vector<ThroughflowStreamSheet> stream_sheets;
    std::vector<ThroughflowIteration> history;
};
[[nodiscard]] std::string to_string(ThroughflowResult::Status s);

/// Runs the analysis (Figs. 7-3, 7-4).  Throws ThroughflowError for an invalid case, a
/// non-physical state, reverse or zero meridional flow on a quasi-normal, or crossing stream
/// surfaces.  A choke is a reported state, not an exception.
[[nodiscard]] ThroughflowResult solve_throughflow(const ThroughflowCase& c);

/// Novak's damping factor, Eq. 7-57: 1 / F = 1 + (1 - M_m^2) y_s^2 / (B* dm^2), with M_m limited to
/// mach_cap (p. 165: M_m <= 0.95).  Throws ThroughflowError unless y_s, dm, B* > 0.
[[nodiscard]] double novak_damping(double mach_m, double y_s, double dm, double bstar, double mach_cap = 0.95);

/// Eq. 7-59, B* -> (B* + 1) / 2, with our trigger (ThroughflowOptions::bstar_growth_threshold): the
/// error of a pass is compared with the smallest error since the last change of B*, so a steady
/// growth below the threshold per pass also triggers it.
struct BstarControl
{
    double bstar = 8.0;
    double threshold = 0.05;
    double min_error = std::numeric_limits<double>::infinity();
    /// Feeds the error of a pass; returns the B* to use for its repositioning.
    double update(double error);
};

/// The fraction t of the crossing a + t (b - a) of the line a -> b with the line through p1, p2
/// (ours, BladeRowSpec).  Throws when they are parallel.
[[nodiscard]] double line_crossing_fraction(double az, double ar, double bz, double br, double z1, double r1,
                                            double z2, double r2);

}  // namespace BladeGen::Turbolenza::Throughflow
