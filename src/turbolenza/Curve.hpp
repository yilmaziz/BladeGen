/// @file turbolenza/Curve.hpp
/// @brief The authoritative C2 blade surfaces: clamped B-spline curves X(u) = (m(u), y(u)).
///
/// Port of v1 geometry.py ParametricBladeCurve and SurfaceCurves (TURBOlenza v1, the C2
/// surface representation frozen 2026-11-14).  One blade surface is a parametric curve from the
/// LE point (u = 0) to the TE point (u = 1): a clamped B-spline of degree p (5 in v1) with simple
/// interior knots, so C^(p-1) at every knot.  m is the curve's first coordinate (the axial
/// direction in the cascade frame, the chord in the blade frame), y the second.
///
/// Evaluation: the B-spline basis and its derivatives by the standard recurrence (Piegl &
/// Tiller, The NURBS Book, 2nd ed., 1997, Algorithms A2.1 and A2.3), the span of u found as
/// t_l <= u < t_{l+1} with the last non-degenerate span at u = 1 (the convention of
/// scipy.interpolate.BSpline, which v1 uses).  At the clamped ends X(0) and X(1) are the end
/// control points exactly (v1 ParametricBladeCurve.derivative).
#pragma once

#include "turbolenza/Field.hpp"

#include <array>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 CURVE_PARAMETER_SLACK: u this far outside [0, 1] is round-off and is clipped.
inline constexpr double kCurveParameterSlack = 1e-12;
/// v1 PROJECTION_ITERATIONS: Newton steps of ParametricBladeCurve::project.
inline constexpr int kProjectionIterations = 8;

class ParametricBladeCurve
{
public:
    ParametricBladeCurve() = default;
    /// Throws std::invalid_argument unless the knots are clamped on [0, 1] with simple interior
    /// knots and n_knots = n_control + degree + 1 (v1 __post_init__).
    ParametricBladeCurve(std::vector<double> knots, std::vector<Point> control, double reference_length, int degree);

    [[nodiscard]] const std::vector<double>& knots() const noexcept { return t_; }
    [[nodiscard]] const std::vector<Point>& control_points() const noexcept { return c_; }
    [[nodiscard]] double reference_length() const noexcept { return L_; }
    [[nodiscard]] int degree() const noexcept { return p_; }
    /// The simple interior knots: the joins between polynomial pieces.
    [[nodiscard]] std::vector<double> interior_knots() const;

    /// d^order X / du^order at u, order 0 .. degree.  u outside [0, 1] by more than
    /// kCurveParameterSlack throws std::invalid_argument (v1 _curve_parameter).
    [[nodiscard]] Point derivative(double u, int order) const;
    [[nodiscard]] Point evaluate(double u) const { return derivative(u, 0); }
    /// X, X', ..., X^(order) at u in one basis evaluation.
    [[nodiscard]] std::vector<Point> derivatives(double u, int order) const;
    [[nodiscard]] Point unit_tangent(double u) const;
    /// Signed curvature (v1 signed_curvature: > 0 turning counterclockwise as u grows).
    [[nodiscard]] double curvature(double u) const;

    /// (u, signed normal distance) of the nearest curve point to q: Newton steps from u_guess
    /// (or from the nearest of 4097 samples when u_guess < 0), clipped to [0, 1]; the distance
    /// is positive on the curve's left (v1 ParametricBladeCurve.project).
    [[nodiscard]] std::array<double, 2> project(const Point& q, double u_guess = -1.0) const;

    /// Every length multiplied by factor; u unchanged.
    [[nodiscard]] ParametricBladeCurve scaled(double factor) const;
    /// Turned counterclockwise by angle [rad] about the origin (v1 rotate_by_stagger).
    [[nodiscard]] ParametricBladeCurve rotated(double angle) const;

private:
    [[nodiscard]] int span(double u) const;

    std::vector<double> t_;
    std::vector<Point> c_;
    double L_ = 1.0;
    int p_ = 0;
};

/// v1 SurfaceCurves: the authoritative surfaces of a cascade in one length unit, in the cascade
/// frame and (_unstaggered) in the blade frame, LE -> TE.  The sampled surfaces of the Cascade
/// are these curves at suction_parameters and pressure_parameters.
struct SurfaceCurves
{
    ParametricBladeCurve suction, pressure;
    ParametricBladeCurve suction_unstaggered, pressure_unstaggered;
    std::vector<double> suction_parameters, pressure_parameters;

    [[nodiscard]] SurfaceCurves scaled(double factor) const;
};

/// v1 geometry.signed_curvature for one point.
[[nodiscard]] double signed_curvature(const Point& d1, const Point& d2);
/// v1 SLOPE_MIN_COS = sqrt(machine epsilon): |cos(beta)| below this, y(m) is singular.
[[nodiscard]] double slope_min_cos() noexcept;
/// v1 geometry.wall_slopes for one point: (dy/dm, d2y/dm2), NaN where |m_u| < min_cos |X_u|.
[[nodiscard]] std::array<double, 2> wall_slopes(const Point& d1, const Point& d2, double min_cos);
/// Blade-frame point turned counterclockwise by stagger about the LE (v1 rotate_by_stagger):
/// X = x cos - y sin, Y = x sin + y cos.
[[nodiscard]] Point rotate_by_stagger(const Point& p, double stagger);

}  // namespace BladeGen::Turbolenza
