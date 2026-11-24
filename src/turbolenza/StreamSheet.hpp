/// @file turbolenza/StreamSheet.hpp
/// @brief The axisymmetric stream sheet of the blade-to-blade solution: the radius r(m) and the
///        sheet thickness b(m) along the meridional coordinate m (Aungier 2003 Sec. 5.1, Fig. 5-2;
///        Sec. 5.3, Eqs. 5-9 to 5-26).  Not in v1, which solves the planar cascade (r = b = 1).
///
/// The equations that carry the sheet (Aungier 2003, read from the page images, printed pp. 81-89):
///  * S = r (theta1 - theta0)                                         [5-9]
///  * tan(beta) = r d theta / dm |eta                                  [5-3]
///  * sin(phi) = dr/dm                                                 [5-15]
///  * mdot psi_eta = S rho b W_m,  mdot psi_m = -rho b (W_theta - W_m tan(beta))   [5-16, 5-17]
///  * A psi_eta,eta - 2 B psi_eta,m + C psi_m,m + D psi_eta + E psi_m = 2 S omega sin(phi)   [5-21]
///    A = mdot / (S b rho cos^2 beta), B = mdot tan(beta) / (b rho), C = mdot S / (b rho)  [5-22 .. 5-24]
///  * h = I + (r omega)^2 / 2 - W^2 / 2                                [5-51]
///
/// The solver's frame (GridData.hpp) stays the planar one: a grid node is (m, y) with
/// y = r_ref theta, r_ref the radius of the sheet at the blade's LE point.  So the pitch is the
/// constant y-distance r_ref 2 pi / N, every grid operation (station cuts, passage checks, eta map)
/// is the planar one, and the sheet enters through the scale sigma(m) = r(m) / r_ref:
/// S = sigma (y1 - y0) and tan(beta) = sigma dy/dm.  A planar cascade has no sheet (sigma = 1,
/// b = 1, sin(phi) = 0): every module is then v1's, bit for bit.
///
/// Representation (ours): station tables (m_k, r_k, b_k), strictly increasing m, interpolated by
/// the not-a-knot cubic spline of CubicSpline.hpp, (r, b) as one planar point per m, so
/// sin(phi) = dr/dm is the derivative of the r interpolant itself.  Lengths in the solver's unit
/// (chords); b in any unit (the solver uses b(m) / b(m_in), Aungier's b normalized to 1 at the far
/// inlet).  A sheet whose r and b are both constant (a cylinder) returns them exactly (sigma = 1,
/// sin(phi) = 0): the cylinder is the planar cascade of the same pitch, to the last bit.
/// Evaluation outside [m_front, m_back] is refused (std::out_of_range): no extrapolation.
///
/// The side lines (Aungier p. 86: beta linear in m from the blade angle to the far flow angle,
/// theta0(m) by "simple integration of Eq. 5-3"): dy/dm = tan(beta(m)) / sigma(m), so
/// y(m) = y_j + int_{m_j}^{m} tan(beta(s)) / sigma(s) ds.  The closed form of v1 (mean_tan) holds
/// for sigma = 1 only; with a sheet the integral is numerical (ours): Gauss-Legendre, 8 points per
/// piece, the pieces split at the table stations and at most kSideLinePiece long.
#pragma once

#include "turbolenza/CubicSpline.hpp"

#include <functional>
#include <optional>
#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

/// The longest piece of the side-line quadrature [solver length unit].
inline constexpr double kSideLinePiece = 0.02;

class StreamSheet
{
public:
    /// m strictly increasing, at least 2 stations; r > 0 and b > 0 finite; r_ref > 0 finite.
    /// std::invalid_argument otherwise.
    StreamSheet(std::vector<double> m, std::vector<double> r, std::vector<double> b, double r_ref);

    [[nodiscard]] double radius(double m) const;      ///< r(m) [solver length unit]
    [[nodiscard]] double sin_phi(double m) const;     ///< dr/dm [5-15]
    [[nodiscard]] double thickness(double m) const;   ///< b(m) [any unit]
    [[nodiscard]] double scale(double m) const { return radius(m) / r_ref_; }  ///< sigma = r / r_ref
    [[nodiscard]] double r_ref() const noexcept { return r_ref_; }
    [[nodiscard]] bool cylinder() const noexcept { return cylinder_; }  ///< r and b constant
    /// A cylinder whose radius is r_ref: sigma = 1 and b / b_in = 1 exactly, the planar cascade itself.
    [[nodiscard]] bool unit_scale() const noexcept { return cylinder_ && r_.front() == r_ref_; }
    [[nodiscard]] double m_front() const noexcept { return m_.front(); }
    [[nodiscard]] double m_back() const noexcept { return m_.back(); }
    [[nodiscard]] const std::vector<double>& m() const noexcept { return m_; }
    [[nodiscard]] const std::vector<double>& r() const noexcept { return r_; }
    [[nodiscard]] const std::vector<double>& b() const noexcept { return b_; }

    /// int_a^b f(s) ds: Gauss-Legendre, 8 points per piece, pieces split at the table stations and
    /// at most kSideLinePiece long (a > b gives the negative integral).
    [[nodiscard]] double integrate(double a, double b, const std::function<double(double)>& f) const;
    /// y(m) - y_j of the side line through (m_j, y_j) whose angle runs linearly in m from
    /// beta_j at m_j to beta_far at m_far: int_{m_j}^{m} tan(beta(s)) / sigma(s) ds.
    [[nodiscard]] double side_line_rise(double m_j, double m, double beta_j, double beta_far, double m_far) const;

private:
    void check(double m) const;
    std::vector<double> m_, r_, b_;
    double r_ref_ = 1.0;
    bool cylinder_ = false;
    std::optional<detail::CubicSpline> spline_;  // (r, b) per m; empty for a cylinder
};

}  // namespace BladeGen::Turbolenza
