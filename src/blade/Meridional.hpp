/// @file blade/Meridional.hpp
/// @brief Meridional streamlines (x axial, r radial), their arc length m and the
///        conformal coordinate m' = \int dm / r.
///
/// The conformal map (m', theta) of a stream surface of revolution is angle
/// preserving: ds^2 = dm^2 + r^2 dtheta^2 = r^2 (dm'^2 + dtheta^2), so a blade
/// profile drawn in (m', theta) has the true flow angles (e.g. J. D. Denton,
/// "Loss mechanisms in turbomachines", and the blade-to-blade plane of
/// ANSYS BladeGen / ref. Aungier, "Axial-Flow Compressors", 2003, ch. 7).
///
/// Beyond its ends a streamline continues along its end tangent (straight line in
/// the meridional plane); m' there follows from r(m) = r_end + r'_end (m - m_end):
///     m' = m'_end + ln(r / r_end) / r'_end      (r'_end != 0),
///     m' = m'_end + (m - m_end) / r_end          (r'_end  = 0).
/// The blade layer uses this to build sections slightly outside hub and casing.
#pragma once

#include "math/BSpline.hpp"

#include <vector>

namespace BladeGen::Blade {

using Math::Vec2;

class MeridionalCurve
{
public:
    MeridionalCurve() = default;
    /// Any B-spline in (x, r); r must stay > 0 on it.
    explicit MeridionalCurve(Math::BSpline2 curve);
    /// C2 cubic (not-a-knot) through (x, r) points, chord-length parameters.
    [[nodiscard]] static MeridionalCurve through_points(const std::vector<Vec2>& xr);
    /// Straight segment from a to b (as a cubic, so it blends with other curves).
    [[nodiscard]] static MeridionalCurve line(const Vec2& a, const Vec2& b);

    [[nodiscard]] const Math::BSpline2& curve() const noexcept { return c_; }
    [[nodiscard]] double length() const noexcept { return m_.back(); }        ///< total m
    [[nodiscard]] double mprime_length() const noexcept { return mp_.back(); } ///< total m'

    /// Point at meridional distance m from the start (extrapolated outside [0, length]).
    [[nodiscard]] Vec2 at_m(double m) const;
    /// Unit tangent d(x, r)/dm at m.
    [[nodiscard]] Vec2 tangent_at_m(double m) const;
    /// m' at meridional distance m, and the inverse.
    [[nodiscard]] double mprime_of_m(double m) const;
    [[nodiscard]] double m_of_mprime(double mprime) const;
    [[nodiscard]] Vec2 at_mprime(double mprime) const { return at_m(m_of_mprime(mprime)); }

    /// Parameter of the curve at meridional distance m in [0, length].
    [[nodiscard]] double t_of_m(double m) const;

    /// Sampled (x, r) at n points equally spaced in m.
    [[nodiscard]] std::vector<Vec2> sample(int n) const;

private:
    void tabulate();
    [[nodiscard]] double integrate_m(double ta, double tb) const;
    [[nodiscard]] double integrate_mp(double ta, double tb) const;

    Math::BSpline2 c_;
    std::vector<double> t_, m_, mp_;  // table: parameter, m, m'
};

/// Hub and shroud; streamline(span) interpolates them pointwise at equal fractions
/// of their meridional lengths (span 0 = hub, 1 = shroud; outside [0, 1] it
/// extrapolates linearly, which is how sections beyond the endwalls are made).
class MeridionalChannel
{
public:
    MeridionalChannel() = default;
    MeridionalChannel(MeridionalCurve hub, MeridionalCurve shroud, int n_samples = 161);

    [[nodiscard]] const MeridionalCurve& hub() const noexcept { return hub_; }
    [[nodiscard]] const MeridionalCurve& shroud() const noexcept { return shroud_; }
    [[nodiscard]] MeridionalCurve streamline(double span) const;

    /// A simple axial compressor channel: straight hub and casing lines between x0 and x1.
    [[nodiscard]] static MeridionalChannel straight(double x0, double x1, double r_hub0, double r_hub1,
                                                    double r_tip0, double r_tip1);

private:
    MeridionalCurve hub_, shroud_;
    int n_ = 161;
};

}  // namespace BladeGen::Blade
