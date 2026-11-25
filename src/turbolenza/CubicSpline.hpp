/// @file turbolenza/CubicSpline.hpp
/// @brief A planar cubic interpolating spline with not-a-knot end conditions: the port of
///        scipy.interpolate.CubicSpline(x, points) as v1 geometry._chord_length_spline uses it.
///
/// Internal to TURBOlenza (the geometric throat, Cascade.cpp); declared here so that the tests
/// can check it.  scipy (1.x, _cubic.py) builds the slopes s_i at the knots from the tridiagonal
/// system of C2 continuity,
///     dx_i s_{i-1} + 2 (dx_{i-1} + dx_i) s_i + dx_{i-1} s_{i+1} = 3 (dx_i d_{i-1} + dx_{i-1} d_i),
/// d_i = (y_{i+1} - y_i) / dx_i, closed by the not-a-knot rows (third derivative continuous at
/// x_1 and x_{n-2}), and solves it with LAPACK gtsv (Gaussian elimination with partial pivoting,
/// ported here as is).  n = 2 gives the line and n = 3 the parabola through the points (scipy's
/// own special cases).  The pieces are PPoly coefficients c0 s^3 + c1 s^2 + c2 s + c3,
/// s = x - x_i, evaluated as scipy's evaluate_poly1 does (ascending powers, not Horner), on the
/// interval x_i <= x < x_{i+1} (the last one for x = x_{n-1}; extrapolated outside).
#pragma once

#include "turbolenza/Field.hpp"

#include <array>
#include <vector>

namespace BladeGen::Turbolenza::detail {

class CubicSpline
{
public:
    /// x strictly increasing, one point per x, n >= 2; std::invalid_argument otherwise.
    CubicSpline(std::vector<double> x, const std::vector<Point>& y);

    /// The nu-th derivative (0 .. 3) at x.
    [[nodiscard]] Point operator()(double x, int nu = 0) const;
    [[nodiscard]] const std::vector<double>& x() const noexcept { return x_; }

private:
    std::vector<double> x_;
    std::vector<std::array<Point, 4>> c_;  // per interval: c0 (s^3) .. c3 (s^0)
};

}  // namespace BladeGen::Turbolenza::detail
