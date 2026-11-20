/// @file fit/SectionFit.hpp
/// @brief Inverse design of a blade section on its stream surface, and of a whole
///        blade from published sections (e.g. NASA Rotor 67).
///
/// A published closed section in conformal coordinates (m', theta) is brought to
/// the unstaggered frame (Fit::normalize_contour), fitted (Fit::fit_profile), and
/// placed back (Blade::SectionPlacement from the frame: LE point, chord, stagger,
/// mirror).  The placement reproduces the frame exactly, so the fitted section
/// in 3-D differs from the published one only by the profile fit.
#pragma once

#include "blade/Blade.hpp"
#include "fit/InverseFit.hpp"

#include <vector>

namespace BladeGen::Fit {

using Math::Vec3;

struct SectionFitResult
{
    Normalized frame;
    FitResult fit;
    Blade::BladeStation station;       ///< span, fitted profile, placement, streamline
    std::vector<Vec3> published_xyz;   ///< the published loop [m]
    std::vector<Vec3> fitted_xyz;      ///< the fitted section loop [m]
    double rms_3d = 0.0, max_3d = 0.0; ///< distance published point -> fitted surface curve [m]
};

/// Placement that maps the unstaggered frame onto `frame` on stream line `sl`.
[[nodiscard]] Blade::SectionPlacement placement_from_frame(const Normalized& frame, const Blade::MeridionalCurve& sl);

/// Fit one section given as a closed loop in (m', theta) of stream line `sl`.
[[nodiscard]] SectionFitResult fit_section(const std::vector<Vec2>& conformal_loop, const Blade::MeridionalCurve& sl,
                                           double span, const Profile::ProfileParams& start,
                                           const FitOptions& opt = {});

/// Distances of points to a closed 3-D polyline.
[[nodiscard]] std::vector<double> polyline_distances(const std::vector<Vec3>& pts, const std::vector<Vec3>& poly);

}  // namespace BladeGen::Fit
