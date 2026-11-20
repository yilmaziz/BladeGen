/// @file profile/ProfileParams.hpp
/// @brief Parameters of one 2-D blade profile in the unstaggered (m', theta) plane.
///
/// Frame: the chord lies on the m' axis, the leading-edge point P_LE = (0, 0),
/// the trailing-edge point P_TE = (1, 0); lengths are fractions of the chord,
/// and the suction side (SS) is the side of larger theta.
///
/// The construction (DLR-style, ported from the blade2d Python engine):
///   * characteristic points at each edge: an ellipse-like edge of semi-axes
///     a = stretch_ab * r (along the edge axis V) and b = r, offset by `asym`
///     normal to V; it fixes P_LE_SS / P_LE_PS (and P_TE_SS / P_TE_PS);
///   * the suction side: cubic B-spline, 5 control points
///     [P_LE_SS, P_LE_SS + h_LE V_LE_SS, CP2, P_TE_SS + h_TE V_TE_SS, P_TE_SS];
///   * the pressure side: cubic B-spline whose interior control points sit at a
///     prescribed normal thickness t(xi) from the suction side;
///   * LE and TE connectors: two-piece quintic B-splines through P_LE (P_TE) that
///     match position, tangent and curvature of the neighbouring sides, so the
///     closed profile is C2 (see EdgeConnector.hpp).
#pragma once

#include <string>
#include <vector>

namespace BladeGen::Profile {

/// How the LE / TE connector meets the suction and pressure sides.
enum class ConnectorMode
{
    /// Default. Curvature-continuous (G2) at both junctions: position, unit tangent
    /// and curvature equal those of the neighbouring side.  Among all such curves it
    /// is the fairest: a degree-5 B-spline with 4 interior knots through the edge
    /// point that minimises the jerk energy  int |C'''|^2  (math/FairCurve.hpp).
    /// The free end speed is tension * (chord length of the connector).
    Fair,
    /// Parametric C2: the connector copies the neighbours' first and second
    /// derivatives (scaled by its parameter span), so the closed loop is one C2
    /// curve in a single parameter.  Exact, but the copied tangential acceleration
    /// of the cubic sides can make the curvature of a small edge spike.
    C2Parametric,
    /// Geometric G2 exactly as the blade2d Python engine (profile_c2_connect.py):
    /// two-span quintic, |r'| = alpha L, zero tangential second derivative.
    G2Legacy,
};

[[nodiscard]] std::string to_string(ConnectorMode m);
[[nodiscard]] ConnectorMode connector_mode_from_string(const std::string& s);

/// Leading- or trailing-edge family.
struct EdgeParams
{
    double radius = 0.012;      ///< b: edge radius (fraction of chord)
    double stretch_ab = 1.2;    ///< a / b: edge stretch along the edge axis V
    double asym = 0.0;          ///< offset of the edge axis normal to V (signed)
    double beta_deg = 15.0;     ///< angle of the edge axis V from the m' axis [deg]
    double d_beta_ss_deg = 20.0;///< wedge half-angle on the suction side [deg]
    double d_beta_ps_deg = 15.0;///< wedge half-angle on the pressure side [deg]
    double tension = 1.0;       ///< connector end speed |r'| / chord length of the connector (shape knob)
};

/// Suction-side B-spline family.
struct SuctionSideParams
{
    double handle_le_frac = 0.16;  ///< LE handle length / |P_TE_SS - P_LE_SS|
    double ctrl2_mprime = 0.40;    ///< m' of the free control point CP2
    double ctrl2_theta = 0.15;     ///< theta of the free control point CP2
    double handle_te_frac = 0.29;  ///< TE handle length / |P_TE_SS - P_LE_SS|
};

/// Thickness-distribution family: a piecewise power law in SS arc length xi.
///   xi <= xi_max: t = t_le + (t_max - t_le) * (xi / xi_max)^shape_le
///   xi >  xi_max: t = t_te + (t_max - t_te) * ((1 - xi) / (1 - xi_max))^shape_te
/// xi_max is the SS arc-length station whose m' equals mprime_max.
struct ThicknessParams
{
    double t_max = 0.08;        ///< peak normal thickness (fraction of chord)
    double mprime_max = 0.30;   ///< chordwise position of the peak (m')
    double t_le = 0.03;         ///< thickness at the SS start (LE side)
    double t_te = 0.006;        ///< thickness at the SS end (TE side)
    double shape_le = 1.0;      ///< exponent of the forward branch (>1 flatter at LE)
    double shape_te = 1.5;      ///< exponent of the aft branch
    bool enforce_c1 = false;    ///< blade2d option: derive shape_te so |dt/dxi| is equal on both sides of the peak (not C1)
};

/// Pressure-side coupling family.
struct PressureSideParams
{
    double handle_le_frac = 0.14;  ///< LE handle length / |P_TE_PS - P_LE_PS|
    double handle_te_frac = 0.18;  ///< TE handle length / |P_TE_PS - P_LE_PS|
    /// m' stations of the free PS control points (thickness applied there).
    std::vector<double> mprime_ctrl = {0.20, 0.45, 0.80};
};

/// Complete profile parameter set (one blade section).
struct ProfileParams
{
    EdgeParams le{};
    EdgeParams te{0.004, 1.2, 0.0, 9.0, 1.0, 1.0, 1.0};
    SuctionSideParams ss{};
    ThicknessParams thickness{};
    PressureSideParams ps{};
    ConnectorMode connector = ConnectorMode::Fair;

    /// NACA 65-(12)10, BladeGen's fairness-weighted inverse fit (Fair connectors,
    /// RMS 0.0104 % chord, no geometry errors).  Default profile.
    [[nodiscard]] static ProfileParams naca65_1210();
    /// The blade2d fit of the same section (naca65_1210_blade_params.json) with the
    /// blade2d connector (G2Legacy, alpha 0.28 / 0.20): the Python engine's geometry.
    /// Its TE connector is S-shaped (found by Profile::check_quality).
    [[nodiscard]] static ProfileParams naca65_1210_blade2d();
    /// The CDA starting point used by blade2d's fit_cda_c2.py (SPEC x0).
    [[nodiscard]] static ProfileParams cda_start();
};

/// One named scalar of a ProfileParams, for generic access (fitting, bindings,
/// spanwise interpolation, notebooks).  `family` groups the knobs.
struct ParamInfo
{
    std::string name;     ///< e.g. "le.radius"
    std::string family;   ///< "Leading edge", "Trailing edge", "Suction side", ...
    std::string label;    ///< human-readable
    double lo;            ///< suggested lower bound
    double hi;            ///< suggested upper bound
};

/// The scalar parameters in a fixed order (PS m' stations are excluded).
[[nodiscard]] const std::vector<ParamInfo>& param_table();
[[nodiscard]] double get_param(const ProfileParams& p, const std::string& name);
void set_param(ProfileParams& p, const std::string& name, double value);
[[nodiscard]] std::vector<double> to_vector(const ProfileParams& p);
void from_vector(ProfileParams& p, const std::vector<double>& v);

/// Linear blend a + w (b - a) of every scalar (connector mode and PS stations from a;
/// PS stations blended too when both have the same count).
[[nodiscard]] ProfileParams blend(const ProfileParams& a, const ProfileParams& b, double w);

}  // namespace BladeGen::Profile
