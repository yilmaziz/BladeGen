/// @file turbolenza/BladeGenSource.hpp
/// @brief The BladeGen adapter: one design layer's profile -> TURBOlenza's BladeSection.
///
/// BladeGen builds every layer profile in TURBOlenza's blade frame (chord 1, LE point (0, 0), TE
/// point (1, 0), suction side +y; design/BladeDesign.hpp), so no transformation is needed: per
/// side, LE -> TE, BladeGen's own curve (Naca65Section or Profile2D, sampled densely at exact
/// curve points with their exact curvature) is approximated by one clamped quintic B-spline with
/// simple interior knots, the representation the solver reads (Curve.hpp, v1's C2 surfaces).
///
/// The fit (least squares, end points fixed).  The samples P_k get the parameter u_k = s_k / S,
/// s_k the chord-length sum of the dense sample (equal to the arc length to ~1e-9 at the sample
/// densities used) and S the side's length, which becomes the curve's reference_length.  For a
/// knot vector the control points c_1 .. c_(N-2) minimise sum_k |X(u_k) - P_k|^2 with
/// c_0 = (0, 0) and c_(N-1) = (1, 0) fixed (so X(0) and X(1) are the end points exactly): the
/// normal equations B^T B c = B^T r are banded and solved by a sparse LDL^T factorisation.  The
/// knots start cosine-clustered toward both ends; every span whose samples miss the position
/// tolerance (|X(u_k) - P_k|) or the curvature tolerance (|kappa_fit - kappa_k| <= tol *
/// max(|kappa_k|, 1)) is halved, and the fit is repeated (adaptive knot insertion; least squares
/// with fixed end points: Piegl & Tiller, The NURBS Book, 2nd ed., 1997, Sec. 9.4.1).  So the knots
/// cluster where BladeGen's curve needs them (the LE and TE fairs and their joins).  A span is
/// halved only while both halves keep 2 (p + 1) samples (and stay wider than 1e-8 in u): the fit
/// stays well posed, and its curvature cannot wiggle between samples.  The position tolerance must
/// be met (else std::runtime_error).  The curvature tolerance is a target: BladeGen's pieces join
/// G2 (curvature continuous) but not G3, and at a join with a jump in dkappa/ds (e.g. 1.4e5 / c^2 at
/// the start of the Kovach TE fair) a C4 spline with simple knots rounds the kink over a knot span;
/// the report gives the remaining curvature error and the overshoot (the wiggle measure).
///
/// Blade angles (Cascade.hpp: kappa_le = stagger + chi_le, kappa_te = stagger - chi_te):
///  * a NACA 65 layer: the equivalent circular arc, tan(camber / 4) = 0.1103 cl0 and
///    chi_le = chi_te = camber / 2 (Aungier 2003 Eq. 4-6, as v1; Naca::equivalent_arc_angles);
///  * a Custom (or fitted NACA 65 + Custom) layer: BladeGen's own camber-line angles at the LE and
///    TE, the ones its designer reports as the metal angles beta_LE, beta_TE minus the stagger
///    (Section2D::camber, profile/Camber.hpp): the camber line is the locus of the inscribed-circle
///    centres, and its end slope is that of a quadratic least-squares fit of the centres within
///    12 % chord of the edge (the direction of the mean line where the inscribed circles end).
///    chi_le = its LE angle (> 0 toward +y), chi_te = minus its TE angle (> 0 toward -y).  (The
///    custom edge-axis angles le.beta_deg / te.beta_deg are construction parameters, not the
///    camber line: for the default naca65_1210 profile they differ from it by 6.8 and 2.8 deg.)
/// LE / TE radius and centre: the NACA 65 metadata (nose circle, Kovach TE circle; 0 for the sharp
/// Emery TE) or the custom edge radius b with the centre of the edge ellipse P_mid.
#pragma once

#include "turbolenza/Cascade.hpp"

#include "design/BladeDesign.hpp"

#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// One side of a profile, LE point -> TE point, as exact curve points.
struct SideSample
{
    std::vector<Point> xy;
    std::vector<double> kappa;  ///< signed curvature of the LE -> TE curve (> 0 turning counterclockwise)
};

struct SideFitOptions
{
    double tolerance = 3e-8;            ///< max |X(u_k) - P_k| at the fit samples [chord]
    double curvature_tolerance = 1e-3;  ///< target max |kappa_fit - kappa_k| / max(|kappa_k|, 1) (see the file comment)
    int initial_spans = 24;             ///< cosine-clustered starting knots
    int max_control_points = 3000;      ///< the refinement stops (std::runtime_error) beyond this
    int n_edge = 6000;                  ///< dense sample: points per edge piece
    int n_body = 20000;                 ///< dense sample: points on the side body
    int n_parameters = 201;             ///< BladeSection::*_parameters (cosine-clustered in u)
};

struct SideFitReport
{
    int control_points = 0;
    int refinements = 0;                ///< knot-insertion rounds
    int samples = 0;
    double arc_length = 0.0;            ///< S = reference_length
    double max_residual = 0.0;          ///< max |X(u_k) - P_k| at the fit samples
    double max_curvature_error = 0.0;   ///< max |kappa_fit(u_k) - kappa_k|
    double max_relative_curvature_error = 0.0;  ///< max of the above / max(|kappa_k|, 1)
    /// The wiggle measure: how far the fit's curvature leaves BladeGen's own curvature range over
    /// its knot span and one span width on either side, / max(|kappa_k|, 1).  0: no overshoot.
    double max_curvature_overshoot = 0.0;
};

struct SectionFitReport
{
    SideFitReport suction, pressure;
    std::string chi_source;             ///< how chi_le / chi_te were found
};

/// The clamped quintic fit of one dense side sample (see the file comment).  Throws
/// std::invalid_argument when the sample is too short or does not run from (0, 0) to (1, 0)
/// (within 1e-9), std::runtime_error when the tolerances need more than max_control_points.
[[nodiscard]] ParametricBladeCurve fit_side(const SideSample& sample, const SideFitOptions& opt = {},
                                            SideFitReport* report = nullptr);

/// BladeGen's own side of a layer's profile (NACA 65 section or custom Profile2D, as the design
/// uses it), n_edge points per edge piece and n_body on the side body.
[[nodiscard]] SideSample sample_layer_side(const Design::Layer& layer, Design::Evaluator& ev, Profile::Side side,
                                           int n_edge, int n_body);

/// u_k = (1 - cos(pi k / (n - 1))) / 2, k = 0 .. n - 1: the sampled surface parameters.
[[nodiscard]] std::vector<double> surface_parameters(int n);

/// The layer's blade section for TURBOlenza (see the file comment).  Throws what BladeGen throws
/// for an invalid profile and what fit_side throws.
[[nodiscard]] BladeSection section_from_layer(const Design::Layer& layer, Design::Evaluator& ev,
                                              SectionFitReport* report = nullptr, const SideFitOptions& opt = {});

// ── the blade on its stream surface (not in v1) ──────────────────────────────
//
// stream_surface_blade: the layer's section as BladeGen places it on its stream surface (the user's
// decisions 1 and 2), in the solver's frame (StreamSheet.hpp, Cascade.hpp StreamSurfaceBlade):
//  * the placement: BladeGen's own (Design::Evaluator::evaluate of the design with the window's copy of
//    the layer at its index): m' = m0 + c' (x cos g - y sin g), theta = t0 + c' (x sin g + y cos g)
//    (Design::Placement; stagger g, conformal chord c', the placed LE point (m0, t0), re-stagger and
//    stacking included); a mirrored design is refused (its suction side faces -theta);
//  * m(m'): the layer streamline X(u) = (x, r)(u) of BladeGen's FlowPath (the exact degree-5 B-spline,
//    hub and casing blended at equal u), m = int |X_u| du and m' = int |X_u| / r du (Gauss-Legendre,
//    8 points on 3 x 512 table intervals of u, then on the remainder), u(m') by Newton; m from the
//    placed LE point.  BladeGen's own StreamSurface tables the same integrals (5 points, 3 x 128
//    intervals, Hermite interpolation): the two maps agree to its interpolation error (reported);
//  * the solver frame: M = m - m_LE, Y = r_LE (theta - t0), r_ref = r_LE; L = |(M, Y)(TE)| the
//    normalizing length (the chord in the (m, r_LE theta) plane); the chord line's angle
//    gamma_my = atan2(Y_TE, M_TE) is the cascade's stagger, and the section the placed blade turned
//    back by it and divided by L (LE (0, 0), TE (1, 0)), fitted per side by fit_side (the dense
//    BladeGen sample mapped point by point, its curvature by the chain rule of the map: dM/ds = r dm'/ds,
//    d2M/ds2 = r d2m'/ds2 + r sin(phi) (dm'/ds)^2, Y'' = r_LE theta'');
//  * the blade angles: the conformal map keeps angles, so the TRUE angles on the stream surface are
//    BladeGen's kappa_LE = g + chi_LE, kappa_TE = g - chi_TE (chi as section_from_layer); the section
//    carries chi_LE - (gamma_my - g), chi_TE + (gamma_my - g), so the cascade's kappa are the true ones;
//  * the pitch: r_LE 2 pi / N (N = design.n_blades), solidity = L / pitch;
//  * the sheet (decision 1, our model): r(m) from the streamline; b(m) = |dX/ds . n|, X(u, s) = (1 - s)
//    hub(u) + s casing(u) (BladeGen's streamline family), n the unit normal of the layer's streamline, at
//    the layer's span s: the normal spacing of the neighbouring stream surfaces, so the sheet carries a
//    fixed fraction of the channel's flow (Aungier takes b from a hub-to-shroud analysis; we take the
//    channel's own family).  Exact for a cylinder (b = casing - hub radius, constant).  Beyond the
//    channel's ends every streamline continues straight along its end tangent (BladeGen's convention),
//    the family with it (X = X_end + (u - u_end) X_u,end): the table reaches coverage_chords beyond the
//    LE and TE, step table_step_chords, stations uniform in m.
struct StreamSurfaceOptions
{
    double length_unit = 1e-3;        ///< BladeGen's length unit in metres (the designer works in mm)
    double coverage_chords = 6.0;     ///< the sheet's table reaches this far beyond the LE and the TE
    double table_step_chords = 0.0025; ///< its station spacing (the spline meets the streamline r to 2e-10: test)
    SideFitOptions fit;
};

struct StreamSurfaceReport
{
    SectionFitReport fit;
    double map_max_difference = 0.0;  ///< max |m(m') - BladeGen's StreamSurface::m_of_mprime| over the blade [chords]
    double b_le = 0.0, b_te = 0.0;    ///< b / b(LE point) at the TE (b_le = 1)
};

/// The (window's copy of) layer k of `design` on its stream surface (see above).  Throws
/// std::invalid_argument for a mirrored design, a layer BladeGen cannot place, or a section that
/// leaves the channel; what fit_side throws.
[[nodiscard]] StreamSurfaceBlade stream_surface_blade(const Design::BladeDesign& design, int k,
                                                      const Design::Layer& layer, Design::Evaluator& ev,
                                                      const StreamSurfaceOptions& opt = {},
                                                      StreamSurfaceReport* report = nullptr);

}  // namespace BladeGen::Turbolenza
