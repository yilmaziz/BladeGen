/// @file turbolenza/ThroughflowBladeGen.hpp
/// @brief The BladeGen adapter of the through-flow analysis (Throughflow.hpp): a BladeGen design's
///        meridional channel, quasi-normals and blade row.
///
/// The channel (blade/FlowPath.hpp): hub and casing are degree-5 B-splines in u in [0, 3] (u in
/// [0, 1] the inlet extension, u = 1 the LE point, u = 2 the TE point, [2, 3] the outlet
/// extension), lengths in mm; the through-flow is in m (length_scale = 1e-3).
///
/// Automatic quasi-normals (ours), upstream to downstream:
///   the far inlet, n_inlet stations in the inlet extension, one just before the LE line, one just
///   after the TE line, n_exit stations in the outlet extension, the far exit;
///  * Radial (the default): radial lines z = const.  The far inlet / exit is at the larger / smaller
///    of the two walls' end z (both walls must exist there); before the LE at
///    min(z_LE,hub, z_LE,casing) - gap, after the TE at max(z_TE,hub, z_TE,casing) + gap, with gap =
///    edge_gap (z_TE,hub - z_LE,hub); the extension stations equally spaced in z between them;
///  * EqualParameter (when the channel needs otherwise, e.g. a strongly leaning edge line or a
///    radial flow): the straight lines hub(u) -> casing(u) at equal u (the LE / TE lines are of this
///    kind), u = 0, ..., 1 - edge_gap, 2 + edge_gap, ..., 3.
/// The end-wall slope and curvature at the quasi-normal ends (Approximate model) are exact, from
/// the wall splines: phi = atan2(r', z'), kappa_m = -dphi/dm = -(z' r'' - r' z'') / (z'^2 + r'^2)^(3/2)
/// (Eq. 7-9 with Eq. 7-2).
///
/// The blade row (one per BladeGen design) lies between the two edge quasi-normals.  Its sections
/// are the design's layers placed on their stream surfaces (Design::place_layer): span = the layer
/// span (BladeGen's layer at span s has its LE / TE at the fraction s of the LE / TE lines, which is
/// the row span fraction of Throughflow.hpp), chord = the physical chord, stagger and metal angles
/// from the placed camber line (LayerResult::gamma, beta_le, beta_te), t_b / c and its position from
/// the camber line, a / c the chordwise position of the largest camber-line ordinate (ours).
///
/// Sign convention (stated here once).  BladeGen's theta is the
/// through-flow's theta (Throughflow.hpp CONVENTIONS): every BladeGen angle is used as it is, and
/// omega = -2 pi rpm / 60.  A positive rpm turns the default BladeGen rotor in the compressor sense:
/// with axial inflow W_theta = -omega r > 0, so beta'_1 > 0 like its LE metal angle kappa_1, and
/// TURBOlenza's blade-to-blade frame (positive relative inlet angle for rpm > 0) needs no sign
/// change.  The Chapter 6 cascade sign is s_row = +1 when the mid-span layer has kappa_1 > kappa_2,
/// else -1.
/// The hub wall segment across the rotor (from the quasi-normal before the LE to the one after the TE)
/// rotates at omega, all other wall segments are at rest (QuasiNormal::hub_wall_omega of the discharge
/// quasi-normal; for the end-wall model only, ours).  compressor_omega = omega (Eq. 8-61's omega_c) is
/// returned for ThroughflowCase::compressor_omega.  The throat o of every section (Eq. 6-53) is the
/// cascade module's geometric throat (Cascade.hpp closest_points: suction side against the neighbour's
/// pressure side) of the placed section in the conformal plane, times r_LE c' (ours, approximate where r
/// varies across the passage).
#pragma once

#include "turbolenza/Throughflow.hpp"

#include "design/BladeDesign.hpp"

#include <vector>

namespace BladeGen::Turbolenza::Throughflow {

struct BladeGenThroughflowOptions
{
    enum class Layout
    {
        Radial,
        EqualParameter,
    };
    Layout layout = Layout::Radial;
    /// Stations strictly inside the inlet / outlet extension (>= 0).  Novak's damping (Eq. 7-57) falls
    /// with (y_s / dm)^2: on the default channel (y_s ~ 63 mm, 6 mm extensions) one station each
    /// converges the full model in ~700 passes, four stations in > 3000 (measured, see the report).
    int n_inlet = 1;
    int n_exit = 1;
    double edge_gap = 0.02;      ///< see the file comment (> 0)
    double length_scale = 1e-3;  ///< BladeGen length unit -> m
    double rpm = 0.0;            ///< > 0: the compressor sense (file comment)
};

struct BladeGenThroughflowGeometry
{
    std::vector<QuasiNormal> quasi_normals;
    BladeRowSpec row;            ///< everything but the model
    int qn_before_le = 0, qn_after_te = 0;
    double compressor_omega = 0.0;  ///< omega of the row (ThroughflowCase::compressor_omega)
};

/// The quasi-normals and the blade row of a design.  Throws std::invalid_argument (or what BladeGen
/// throws) for an invalid design, a layer that cannot be built, or stations that do not fit.
[[nodiscard]] BladeGenThroughflowGeometry throughflow_geometry(const Design::BladeDesign& d, Design::Evaluator& ev,
                                                               const BladeGenThroughflowOptions& opt = {});

/// The end-wall slope phi and curvature kappa_m (Eqs. 7-2, 7-9) of a wall spline at parameter u.
void wall_slope_curvature(const Math::BSpline2& wall, double u, double& phi, double& kappa);

}  // namespace BladeGen::Turbolenza::Throughflow
