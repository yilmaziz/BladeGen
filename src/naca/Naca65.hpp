/// @file naca/Naca65.hpp
/// @brief Parametric NACA 65-series blade section (A10 mean line), C2 design geometry.
///
/// Migrated from TURBOlenza (turbolenza/src/turbolenza/b2b/geometry.py, the frozen C2
/// construction of 2026-11-14), the blade geometry of the BladeGen -> ENIGMAgrid ->
/// TURBOlenza chain.  Frame and conventions are TURBOlenza's: chord 1, LE point (0, 0),
/// TE point (1, 0), camber toward +y, suction side = +y side.
///
/// Kept from TURBOlenza (sources: Aungier 2003 Secs. 4.1, 4.2, Table 4-1; Emery et al.
/// 1958; Kovach & Sandercock 1961):
///  * Table 4-1 (NACA 65-010 thickness, a = 1.0 mean line), with the 2.5 % misprint
///    corrected (3.480 -> 3.148 %); thickness scaled by t_max_over_c / 0.10;
///  * the closed-form A10 mean line y/c = -(cl0 / 4 pi) [(1-x) ln(1-x) + x ln x];
///  * the thickness station = arc length along the camber line / its length, s / S
///    (straight end pieces of XI_MIN at both ends, as TURBOlenza);
///  * the thickness reconstruction "quintic_natural_sqrt": the interpolating quintic
///    spline in z = sqrt(station) with h''' = h'''' = 0 at both ends (C4, exact at every
///    station); the legacy PCHIP (C1, curvature sawtooth) is not migrated;
///  * surfaces X = X_c +- (t/2) n_c;
///  * the LE nose circle R_LE = 0.00687 (t/0.10)^2 through (0, 0), centred on the mean-line
///    tangent at x/c = 0.005, and the C2 leading edge "monotone_curvature_fair": per side a
///    Whewell-form fair r(phi), r = R_LE for phi <= 20 deg, then a cubic B-spline with
///    non-decreasing control radii (curvature falls monotonically from 1/R_LE: no
///    overshoot, no inflection, no S-bend), G2 onto the surface at an anchor, least
///    int (dkappa/ds)^2 ds (reweighted LSEI); the anchor 0.5 R_LE past the earliest anchor
///    at which such a fair exists (scan + bisection); LE_FAIRNESS_LIMIT reported;
///  * the two profile variants: "emery_base" (the table to a sharp TE) and
///    "kovach_sandercock" (the table up to station 0.60, then a straight taper tangent to a
///    TE circle of radius 0.008 c centred on the camber line and through (1, 0));
///  * the equivalent circular-arc blade angles tan(camber / 4) = 0.1103 cl0.
///
/// BladeGen's own choices (instead of TURBOlenza's kernel smoothing + least-squares C2 refit):
///  * Kovach start (station 0.60): the table / taper corner (0.08 deg) is blended C2 by a
///    quintic smoothstep over station 0.60 +- KOVACH_BLEND (deviation ~1e-6 c);
///  * Kovach TE: the pressure surface is concave near the TE, so no monotone fair reaches
///    it; the TE is a G2 fair connector (minimum jerk, math/FairCurve.hpp) from each
///    surface's tangency point with the TE circle through (1, 0), with the circle's
///    tangent and curvature 1/R_TE there;
///  * Emery TE: each surface ends with a C2 closing piece from x/c = 0.995 to (1, 0) along
///    TURBOlenza's closing direction: a sharp TE (wedge), C2 up to it.
#pragma once

#include "math/BSpline.hpp"
#include "profile/Profile2D.hpp"
#include "profile/Quality.hpp"

#include <array>
#include <string>
#include <vector>

namespace BladeGen::Naca {

using Math::Vec2;

enum class Variant { KovachSandercock, EmeryBase };

[[nodiscard]] std::string to_string(Variant v);
[[nodiscard]] Variant variant_from_string(const std::string& s);

struct Naca65Params
{
    double cl0 = 1.2;           ///< design lift coefficient of the A10 mean line (the "12" in 65-(12A10)10 / 10)
    double t_max_over_c = 0.10; ///< nominal thickness ratio (the last "10")
    Variant variant = Variant::KovachSandercock;

    void validate() const;
    [[nodiscard]] std::string designation() const;  ///< e.g. "NACA 65-(12A10)10"
};

/// Table 4-1 (Aungier 2003, p. 63), fractions of chord: x, mean line (Cl0 = 1), 65-010
/// thickness as printed, and corrected.
struct Table41
{
    std::array<double, 26> x, camber, thickness_printed, thickness;
};
[[nodiscard]] const Table41& table41();

/// The A10 mean line and its first three derivatives at x in (0, 1).
[[nodiscard]] std::array<double, 4> a10(double x, double cl0);

/// Equivalent circular arc: camber angle, chi_LE, chi_TE [rad] (Aungier Eq. 4-6).
[[nodiscard]] std::array<double, 3> equivalent_arc_angles(double cl0);

/// Whewell-form edge fair: radius of curvature r against the turning angle phi.
struct WhewellFair
{
    Vec2 apex{0.0, 0.0};
    double apex_angle = 0.0;  ///< tangent angle at the apex
    double sense = 1.0;       ///< +1 counter-clockwise turning, -1 clockwise
    double radius = 0.0;      ///< R of the edge circle
    double min_arc = 0.0;     ///< circle kept up to this turning angle
    double anchor_angle = 0.0;
    std::vector<double> control; ///< non-decreasing control radii, first ones = R
    double fairness = 0.0;    ///< R^3 int (dkappa/ds)^2 ds
    double anchor_s = 0.0;    ///< source arc length of the anchor
    double earliest_s = 0.0;  ///< earliest source arc length with a fair

    [[nodiscard]] double r(double phi) const;
    [[nodiscard]] double r_slope(double phi) const;
    [[nodiscard]] Vec2 point(double phi) const;
    [[nodiscard]] Vec2 tangent(double phi) const;
};

struct Naca65Metadata
{
    std::string designation;
    double le_radius = 0.0;
    Vec2 le_center{0.0, 0.0};
    double te_radius = 0.0;        ///< 0 for the sharp Emery TE
    Vec2 te_center{1.0, 0.0};
    double camber_angle = 0.0, chi_le = 0.0, chi_te = 0.0;
    double max_thickness = 0.0, max_thickness_station = 0.0;
    double le_fairness[2] = {0, 0};      ///< suction, pressure
    double le_anchor_x[2] = {0, 0};
    double le_anchor_angle_deg[2] = {0, 0};
    double te_start_x[2] = {0, 0};       ///< where each surface hands over to the TE piece
    double te_circle_deviation = 0.0;    ///< Kovach: max | |X - centre| - R_TE | along the TE fair
    double camber_length = 0.0;
};

class Naca65Section
{
public:
    explicit Naca65Section(Naca65Params p);

    [[nodiscard]] const Naca65Params& params() const noexcept { return p_; }
    [[nodiscard]] const Naca65Metadata& metadata() const noexcept { return meta_; }

    // ── laws ─────────────────────────────────────────────────────────────────
    [[nodiscard]] double station(double x) const;          ///< s / S of the camber point at x
    [[nodiscard]] double x_of_station(double station) const;
    [[nodiscard]] double half_thickness(double station) const;  ///< incl. the Kovach taper
    [[nodiscard]] double table_half_thickness(double station) const; ///< the C4 table law only
    [[nodiscard]] Vec2 camber_point(double x) const;
    /// Surface point and its first two x-derivatives (side +1 suction, -1 pressure).
    [[nodiscard]] std::array<Vec2, 3> surface(double x, double side) const;

    // ── geometry ─────────────────────────────────────────────────────────────
    [[nodiscard]] const WhewellFair& le_fair(Profile::Side s) const { return le_[static_cast<int>(s)]; }
    /// One side from the LE point to the TE point: LE fair, table surface, TE piece.
    /// Pieces as in Profile2D (LeadingEdge, SuctionSide / PressureSide, TrailingEdge).
    [[nodiscard]] Profile::CurveSample sample_side(Profile::Side side, int n_edge = 60, int n_body = 240) const;
    /// Closed loop from the LE point over the suction side, the TE point and back
    /// (last point = first), the same layout as Profile2D::sample_loop.
    [[nodiscard]] Profile::CurveSample sample_loop(int n_edge = 40, int n_body = 160) const;
    [[nodiscard]] bool sharp_te() const noexcept { return p_.variant == Variant::EmeryBase; }

private:
    struct Point3 { Vec2 p, d1, d2; };
    void build();
    void build_le(int k);
    void build_te();
    [[nodiscard]] Point3 te_piece(int k, double u) const;   ///< u in [0, 1]: surface end -> TE point
    [[nodiscard]] double arc_station_raw(double x) const;

    Naca65Params p_;
    Naca65Metadata meta_;
    Math::BSpline<1> H_;                 // half thickness of the 65-010 base form in z = sqrt(station)
    std::vector<double> arc_x_, arc_s_;  // camber arc-length table
    double S_ = 0.0, end_piece_ = 0.0;
    // Kovach taper
    double sig_tan_ = 1.0, h_tan_ = 0.0, h_start_ = 0.0;
    WhewellFair le_[2];
    double x_body_start_[2] = {0, 0}, x_body_end_[2] = {1, 1};
    Math::BSpline2 te_conn_;  // Kovach: SS end -> TE point -> PS end
    double te_umid_ = 0.5;
    Math::BSpline2 te_close_[2];  // Emery: closing pieces
};

/// Profile quality of a NACA 65 section (the same checks as Profile2D's; a sharp Emery TE
/// is a design corner, reported as information).
[[nodiscard]] Profile::QualityReport check_quality(const Naca65Section& s, const Profile::QualityOptions& opt = {});

}  // namespace BladeGen::Naca
