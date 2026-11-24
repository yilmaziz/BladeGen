/// @file turbolenza/ThroughflowSheet.hpp
/// @brief A Chapter 7 through-flow stream surface as the blade-to-blade solver's stream sheet r(m), b(m), and the
///        blade-to-blade case on it: Aungier (Axial-Flow Compressors, ASME Press, 2003), Ch. 12 p. 287, "the
///        simple one-way approach".  Frames: ThroughflowCoupling.hpp (FRAMES).
///
/// THE RULES (throughflow_surface_sheet, solver_stream_sheet, blade_to_blade_on_surface).
/// Ch. 12 p. 287, the simple one-way approach: the stream-surface geometry from the through-flow at the
/// quasi-normals bounding the row, "assuming a linear variation between those stations".
///  * The surface j of the through-flow's last pass: z, r, m (Eq. 7-58, m = 0 on QN 0) at every quasi-normal,
///    and the stream-sheet thickness b at the surface: SheetThickness::Geometric (default) the stream-tube
///    spacing normal to the streamline, cos(epsilon) times the tubes' y-width, the area term of Eq. 7-32
///    (ThroughflowStreamSurface::b); SheetThickness::Mass the thickness that carries the tube's actual mass,
///    Eq. 12-39 with Ch. 7's blockages (ThroughflowStreamSurface::b_mass).  Both in metres; the solver uses
///    b(m) / b(m_in) only (StreamSheet.hpp).
///  * r(m) and b(m) linear in m between neighbouring quasi-normals (p. 287 inside the row; ours outside it),
///    constant beyond the first and last quasi-normal (ours: a cylindrical extension).
///  * The LE and TE on the surface: where the straight segment between the surface's points on QN1 and QN2
///    crosses the row's LE and TE lines (Throughflow.hpp line_crossing_fraction), at m linear between them.
///  * The solver's sheet: m from the LE crossing divided by the blade's normalizing length (StreamSurfaceBlade::
///    chord), r divided by it, b in metres; a dense table (step table_step_chords, the cubic spline of
///    StreamSheet.hpp follows the piecewise-linear r to the table's resolution); r_ref the blade's own (so the
///    pitch r_ref 2 pi / Z of make_cascade stays the blade's), or r_LE / L for a blade without a sheet.
///  * The blade: the caller's StreamSurfaceBlade (e.g. BladeGenSource's stream_surface_blade of the layer at
///    the surface's span, or a section by hand) with its sheet replaced: the section is not re-placed on the
///    through-flow surface (ours; the difference of r_LE is reported).
///  * The inlet of the case: the through-flow's absolute state on the surface at the grid's nominal far inlet,
///    m_LE - upstream_extension_chords L: P_t, T_t, r C_theta and W_m linear in m between the quasi-normals of
///    the blade-free stretch that ends at the row's QN1 (rC_theta, H and s are conserved there, Eqs. 7-40 ..
///    7-42), C_theta = r C_theta / r, the static state from H; QN 0's state before QN 0, QN1's beyond QN1 (ours:
///    Ch. 7 knows nothing between QN1 and the LE).  A far inlet inside an upstream row is refused.  A rotor (omega != 0) runs as RowMotion::Rotating with rpm = -60 omega / (2 pi), a
///    stator as RowMotion::Stationary with the absolute (= relative) Mach number and angle.
///  * A row with s_row = -1 runs mirrored (theta -> -theta): rpm and the inlet angle sign-changed; the caller's
///    blade is then the mirrored section (suction side +y, Cascade.hpp); book angle = s_row * the case's angle.
#pragma once

#include "turbolenza/Case.hpp"
#include "turbolenza/Throughflow.hpp"

#include <memory>
#include <vector>

namespace BladeGen::Turbolenza::Throughflow {

enum class SheetThickness
{
    Geometric,  ///< the stream-tube spacing normal to the streamline (Eq. 7-32's area term)
    Mass,       ///< Eq. 12-39 with Ch. 7's blockages
};

/// One through-flow stream surface through a blade row (see the file comment).  Lengths in metres.
struct ThroughflowSurfaceSheet
{
    int row = -1, surface = -1;
    int qn_inlet = 0, qn_outlet = 0;
    int qn_free = 0;                   ///< the first quasi-normal of the blade-free stretch that ends at qn_inlet
    int s_row = 1;
    double omega = 0.0;                ///< the row's [rad/s], book sign
    double mass_fraction = 0.0;
    /// The surface's row span fraction as the through-flow takes it for the row model (BladeRowSpec): the mean
    /// of the fractions along the LE and TE lines where the segment QN1 -> QN2 crosses them; 0 / 1 on the walls.
    /// The span of the blade section to put on the sheet (e.g. the BladeGen layer).
    double span = 0.0;
    SheetThickness thickness = SheetThickness::Geometric;
    std::vector<double> m, z, r, b;    ///< at every quasi-normal (b by `thickness`)
    /// the absolute flow on the surface at every quasi-normal (for the inlet of the case)
    std::vector<double> Pt, Tt, rCtheta, Wm;
    double m_le = 0.0, m_te = 0.0;     ///< m of the LE and TE crossings
    double z_le = 0.0, r_le = 0.0, z_te = 0.0, r_te = 0.0;
    double beta2 = 0.0;                ///< the through-flow's beta' on QN2 (the row model's discharge angle), book sign

    /// r(m) and b(m) by the rule of the file comment (any m).
    [[nodiscard]] double radius(double m) const;
    [[nodiscard]] double thickness_at(double m) const;
};

/// The surface `surface` (0 hub .. n_ss - 1 shroud) of `result` through row `row` of `c`.  ThroughflowError for an
/// index out of range, a result that does not reach the row's discharge quasi-normal, or a surface whose segment
/// between the row's quasi-normals misses the LE or TE line.
[[nodiscard]] ThroughflowSurfaceSheet throughflow_surface_sheet(const ThroughflowCase& c, const ThroughflowResult& result,
                                                                int row, int surface,
                                                                SheetThickness thickness = SheetThickness::Geometric);

/// The solver's StreamSheet of a surface sheet: M = (m - m_LE) / L on [M_lo, M_hi] (chords), r / L, b [m], r_ref
/// [chords]; stations every `step` chords from M_lo (and M_hi).  std::invalid_argument for L <= 0, r_ref <= 0,
/// M_lo >= M_hi or step <= 0.
[[nodiscard]] std::shared_ptr<const StreamSheet> solver_stream_sheet(const ThroughflowSurfaceSheet& s, double L,
                                                                     double r_ref, double M_lo, double M_hi,
                                                                     double step = 0.0025);

/// The absolute state on the surface at m (see the file comment).
struct SurfaceInletState
{
    double m = 0.0, r = 0.0;
    double Pt = 0.0, Tt = 0.0, Ctheta = 0.0, Wm = 0.0;
    double T = 0.0, a = 0.0;
    double mach = 0.0;   ///< C / a
    double alpha = 0.0;  ///< atan2(C_theta, C_m) [rad], book sign
};
/// ThroughflowError when m lies upstream of the blade-free stretch (inside an upstream row) or the static enthalpy
/// is not positive.
[[nodiscard]] SurfaceInletState surface_inlet_state(const ThroughflowSurfaceSheet& s, const ThroughflowGas& gas, double m);

struct BladeToBladeOptions
{
    SheetThickness thickness = SheetThickness::Geometric;
    double table_step_chords = 0.0025;  ///< the solver sheet's station spacing
    double margin_chords = 1.0;         ///< the sheet reaches this far beyond the grid's far inlet and exit
};

/// A blade-to-blade case on a through-flow stream surface.
struct BladeToBladeOnSurface
{
    ThroughflowSurfaceSheet sheet;
    StreamSurfaceBlade blade;           ///< the caller's, with the through-flow's sheet
    CaseSettings settings;              ///< the caller's, with the stream surface, the row's motion and the inlet
    SurfaceInletState inlet;            ///< at the nominal far inlet (book frame)
    int frame_sign = 1;                 ///< s_row: book angle = frame_sign * the case's angle
    double rpm = 0.0;                   ///< the case's (frame_sign * -60 omega / (2 pi))
    double r_le_ratio = 1.0;            ///< the through-flow's r_LE over the blade's own (r_ref L): the placement mismatch
};

/// See the file comment.  `blade` needs n_blades > 0 and chord > 0; a blade with a sheet keeps its r_ref and
/// solidity, one without gets r_ref = r_LE / L and solidity = L Z / (2 pi r_LE).  std::invalid_argument (or
/// ThroughflowError) for what does not fit.
[[nodiscard]] BladeToBladeOnSurface blade_to_blade_on_surface(const ThroughflowCase& c, const ThroughflowResult& result,
                                                              int row, int surface, StreamSurfaceBlade blade,
                                                              CaseSettings base, const BladeToBladeOptions& opt = {});

}  // namespace BladeGen::Turbolenza::Throughflow
