/// @file fit/InverseFit.hpp
/// @brief Inverse design: profile parameters that reproduce published coordinates.
///
/// Port of blade2d's tools/profile_inverse_fit.py (with the residual kernel of
/// CDA/fit_cda_c2.py), in C++:
///   * target points are in the unstaggered frame (LE (0,0), TE (1,0), SS above);
///     normalize_contour() brings a closed published contour there;
///   * residuals are side aware: suction-side points are measured against the
///     upper path (LE connector + SS + TE connector), pressure-side points against
///     the lower path (LE connector + PS + TE connector), so a mid-chord SS point
///     cannot snap onto the PS;
///   * bounded Levenberg–Marquardt (math/LeastSquares.hpp) in several stages;
///   * optional fairness penalty: concave dips of the LE / TE connector curvature
///     below its junction values (the S-shaped edges found by Profile/Quality) are
///     penalised, so a fit cannot hide an edge defect.
#pragma once

#include "profile/Profile2D.hpp"

#include <string>
#include <vector>

namespace BladeGen::Fit {

using Math::Vec2;

/// One fitted scalar: profile parameter name (Profile::param_table), start, bounds.
struct FitParam
{
    std::string name;
    double x0;
    double lo;
    double hi;
};

/// The blade2d fit_cda_c2.py SPEC shape parameters (24) with its bounds.
[[nodiscard]] std::vector<FitParam> default_fit_params(const Profile::ProfileParams& start);

struct FitOptions
{
    std::vector<FitParam> params;   ///< empty: default_fit_params(start)
    int stages = 3;                 ///< LM restarts from the previous result
    int max_iterations = 150;       ///< per stage
    double fairness_weight = 0.0;   ///< weight of the connector-dip penalty (0: off)
    int n_edge = 120;               ///< samples per connector
    int n_body = 320;               ///< samples per side
};

struct FitResult
{
    Profile::ProfileParams params;
    double rms = 0.0;               ///< RMS point distance (chord fractions)
    double max = 0.0;               ///< max point distance (chord fractions)
    double rms_ss = 0.0, rms_ps = 0.0;
    std::vector<double> d_ss, d_ps; ///< distance of every target point
    std::vector<std::string> at_bound;
    int evaluations = 0;
    double seconds = 0.0;
    std::string status;
};

/// Side-aware distances of the target points to a profile (chord fractions).
void side_distances(const Profile::Profile2D& p, const std::vector<Vec2>& ss, const std::vector<Vec2>& ps,
                    std::vector<double>& d_ss, std::vector<double>& d_ps, int n_edge = 120, int n_body = 320);

/// Fit the profile to suction- and pressure-side points (unstaggered frame).
[[nodiscard]] FitResult fit_profile(const std::vector<Vec2>& ss, const std::vector<Vec2>& ps,
                                    const Profile::ProfileParams& start, const FitOptions& opt = {});

/// A closed published contour brought to the unstaggered frame.
struct Normalized
{
    std::vector<Vec2> ss, ps;  ///< target points, LE -> TE, SS above
    Vec2 le, te;               ///< LE / TE points in the input frame
    double chord = 0.0;        ///< |TE - LE| in the input frame
    double stagger = 0.0;      ///< angle of TE - LE from the input x axis [rad]
    bool mirrored = false;     ///< the input had its convex side toward -y

    [[nodiscard]] Vec2 to_input(const Vec2& xieta) const;  ///< back to the input frame
};

/// LE and TE are the two contour points farthest apart (the longest chord); the
/// contour is split there into its two sides; the side with the larger mean eta
/// becomes the suction side (mirroring the frame when needed so that the mean
/// camber is positive).  `points` is a closed loop (last point may repeat the first).
[[nodiscard]] Normalized normalize_contour(const std::vector<Vec2>& points);

}  // namespace BladeGen::Fit
