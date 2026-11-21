/// @file profile/Camber.hpp
/// @brief Camber line and thickness of a closed profile: the inscribed-circle (medial
///        axis) construction, independent of how the profile was built.
///
/// For a suction-side point p with unit inward normal n, the inscribed circle tangent to
/// the suction side at p has its centre c = p + rho n where the distance from c to the
/// pressure side equals rho:
///
///     f(rho) = dist(p + rho n, PS) - rho = 0,     f(0) = dist(p, PS) > 0,
///
/// f is strictly decreasing for rho below the local radius of curvature, so the root is
/// bracketed by doubling and found by bisection.  The camber line is the locus of the
/// centres c (closed by the LE and TE points), the thickness is 2 rho, and the camber
/// (metal) angle is the slope angle of that locus, smoothed by a least-squares cubic
/// B-spline y(x) pinned at the LE and TE points.  This is the classical definition of
/// the mean line of a blade section as the locus of the centres of the inscribed circles
/// (e.g. S. L. Dixon and C. A. Hall, "Fluid Mechanics and Thermodynamics of
/// Turbomachinery", 7th ed., 2014, Sec. 3.3; Aungier 2003, Sec. 4.1, where the NACA
/// surfaces are built from the mean line by the same normal offsets).  Used here for
/// NACA 65 and Custom sections alike, so both are measured the same way.
#pragma once

#include "math/Vec.hpp"

#include <vector>

namespace BladeGen::Profile {

using Math::Vec2;

struct CamberLine
{
    std::vector<Vec2> point;        ///< n points, cosine-spaced in x, LE point to TE point (smoothed locus)
    std::vector<double> thickness;  ///< 2 rho at every point (0 at the LE / TE points), chord units
    std::vector<double> angle_deg;  ///< slope angle atan2(dy, dx) of the camber line [deg]
    double max_thickness = 0.0;
    double max_thickness_x = 0.0;   ///< chordwise position (x) of the maximum
};

/// Camber line of a profile given its two sides as polylines from the LE point to the
/// TE point (suction side = +y side, x from the LE to the TE).  n SS stations,
/// cosine-clustered in SS arc length.  Throws std::invalid_argument on degenerate input.
[[nodiscard]] CamberLine camber_line(const std::vector<Vec2>& ss, const std::vector<Vec2>& ps, int n = 101);

/// Distance from q to a polyline.
[[nodiscard]] double polyline_distance(const Vec2& q, const std::vector<Vec2>& poly);

/// Signed area of a closed polygon (counter-clockwise > 0).
[[nodiscard]] double polygon_area(const std::vector<Vec2>& loop);

}  // namespace BladeGen::Profile
