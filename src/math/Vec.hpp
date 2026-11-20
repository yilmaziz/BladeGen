/// @file math/Vec.hpp
/// @brief Small fixed-size vectors (Eigen) and 2-D helpers used across BladeGen.
#pragma once

#include <Eigen/Dense>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace BladeGen::Math {

template <int D>
using Vec = Eigen::Matrix<double, D, 1>;
using Vec2 = Vec<2>;
using Vec3 = Vec<3>;

inline constexpr double kPi = std::numbers::pi;

[[nodiscard]] inline double deg2rad(double d) noexcept { return d * kPi / 180.0; }
[[nodiscard]] inline double rad2deg(double r) noexcept { return r * 180.0 / kPi; }

/// z-component of the 2-D cross product a x b.
[[nodiscard]] inline double cross2(const Vec2& a, const Vec2& b) noexcept
{
    return a.x() * b.y() - a.y() * b.x();
}

/// Left (counter-clockwise) normal of v: (-v_y, v_x).
[[nodiscard]] inline Vec2 left_normal(const Vec2& v) noexcept { return {-v.y(), v.x()}; }

/// Right (clockwise) normal of v: (v_y, -v_x).
[[nodiscard]] inline Vec2 right_normal(const Vec2& v) noexcept { return {v.y(), -v.x()}; }

/// Unit vector; throws on a (near) zero vector rather than inventing a direction.
template <int D>
[[nodiscard]] Vec<D> unit(const Vec<D>& v)
{
    const double n = v.norm();
    if (!(n > 1e-300)) throw std::invalid_argument("unit(): zero-length vector");
    return v / n;
}

/// Signed curvature of a planar parametric curve from its first two derivatives.
[[nodiscard]] inline double signed_curvature(const Vec2& d1, const Vec2& d2)
{
    const double s2 = d1.squaredNorm();
    if (!(s2 > 0.0)) throw std::invalid_argument("signed_curvature: zero parametric speed");
    return cross2(d1, d2) / (s2 * std::sqrt(s2));
}

}  // namespace BladeGen::Math
