/// @file profile/Quality.hpp
/// @brief Geometry quality checks of a profile: curvature along both sides from
///        the LE point to the TE point, continuity, self-intersection, thickness.
///
/// The checks never stop the construction: a profile with problems is still
/// built and returned, and the report lists every problem with its severity and
/// location, so a plot can mark it (ERROR in red, WARNING in orange).
#pragma once

#include "profile/Profile2D.hpp"

#include <string>
#include <vector>

namespace BladeGen::Profile {

enum class Severity : int { Info = 0, Warning = 1, Error = 2 };

[[nodiscard]] std::string to_string(Severity s);

struct Issue
{
    Severity severity = Severity::Warning;
    std::string code;     ///< machine-readable, e.g. "ss_inflection"
    std::string message;  ///< human-readable
    Vec2 where{0.0, 0.0}; ///< location in (m', theta)
    int side = -1;        ///< 0 = suction, 1 = pressure, -1 = whole profile
    double s = 0.0;       ///< arc-length fraction from the LE point along that side
};

struct QualityOptions
{
    /// curvature-extremum detection ignores ripples below this fraction of the
    /// side's curvature range
    double wiggle_rel = 0.02;
    /// more curvature extrema than this on one side body is a warning
    int max_extrema_body = 3;
    double junction_kappa_rel = 1e-6;  ///< |dkappa| / max|kappa| allowed at a junction
    double junction_angle = 1e-7;      ///< tangent angle allowed at a junction [rad]
    double min_thickness_frac = 0.0;   ///< measured thickness below this (chord) is an error
    int n_edge = 80;                   ///< sampling of every half connector
    int n_body = 400;                  ///< sampling of each side body
};

struct QualityReport
{
    std::vector<Issue> issues;
    CurveSample ss;            ///< suction side, LE point -> TE point
    CurveSample ps;            ///< pressure side, LE point -> TE point
    std::array<Junction, 4> junctions{};
    double max_curvature = 0.0;
    double min_thickness = 0.0;
    int ss_inflections = 0;
    int ps_inflections = 0;
    int ss_extrema = 0;
    int ps_extrema = 0;
    int self_intersections = 0;

    [[nodiscard]] int count(Severity s) const;
    [[nodiscard]] bool ok() const { return count(Severity::Error) == 0; }
    [[nodiscard]] std::string summary() const;
};

[[nodiscard]] QualityReport check_quality(const Profile2D& prof, const QualityOptions& opt = {});

/// Number of crossings between non-adjacent segments of a closed polyline
/// (the last point must repeat the first).
[[nodiscard]] int count_self_intersections(const std::vector<Vec2>& closed);

}  // namespace BladeGen::Profile
