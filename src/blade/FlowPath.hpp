/// @file blade/FlowPath.hpp
/// @brief Initial meridional configuration of an axial blade row: hub and casing from
///        eight cylindrical coordinates, layer streamlines, and the conformal map
///        (m', theta) -> (x, r, theta) of each layer's surface of revolution.
///
/// Input (axial compressor, image of the "initial meridional configuration"):
///
///     casing  (z_LE, R_LE) -------------- (z_TE, R_TE)
///               |          blade            |
///     hub     (z_LE, R_LE) -------------- (z_TE, R_TE)
///
/// z is the axial coordinate (the rotation axis, BladeGen's x), R the radius.  Each wall
/// between its LE and TE points follows one of two laws:
///
///  * Linear: the straight line through the LE and TE points, continued straight
///    upstream and downstream (no kink anywhere);
///  * BellMehta: the fifth-order contraction contour of Bell & Mehta (1988),
///        R(xi) = R_LE + (R_TE - R_LE) (6 xi^5 - 15 xi^4 + 10 xi^3),  xi = (z - z_LE) / (z_TE - z_LE),
///    whose first and second derivatives vanish at both ends, so the wall continues as
///    a cylinder (R = const) upstream and downstream with continuous slope and
///    curvature (G2).  J. H. Bell and R. D. Mehta, "Contraction design for small
///    low-speed wind tunnels", NASA CR-182747, 1988.
///
/// The inlet / outlet extensions are fractions of each wall's axial LE -> TE length.
///
/// Exact representation.  Each wall is a degree-5 B-spline in u in [0, 3]: u in [0, 1]
/// the inlet extension, [1, 2] the bladed part (u = 1 the LE point, u = 2 the TE point),
/// [2, 3] the outlet extension; z is linear in u on every piece.  The Bell & Mehta
/// polynomial is the Bezier curve with R control values R_LE + dR [0, 0, 0, 1, 1, 1]
/// (its Bernstein coefficients), so the wall is exact, not fitted.  Hub and casing share
/// the knot vector, so a layer streamline at span fraction s is the B-spline whose
/// control points are (1 - s) P_hub + s P_casing: the exact pointwise blend at equal u,
/// which puts every layer's LE / TE on the straight LE / TE lines between the walls.
/// Because the extensions scale with each wall's axial length, the blend keeps G1 (and,
/// for Bell & Mehta / Linear, G2) at the piece joins.
///
/// Conformal map (Drela & Youngren, MISES user's guide 1989; Vos et al., Computers and
/// Fluids 292 (2025) 106588, Eq. 2): m' = int dm / r measured from the layer's LE point.
/// A profile drawn in (m', theta) maps to (x(m'), r(m'), theta), angle preserving.
///
///  * cylindrical wrapping: on a layer of constant radius r0 (straight cylindrical
///    walls) the map is closed form, x = x_LE + r0 m', r = r0: the (m', theta) plane is
///    the unrolled cylinder (r0 m', r0 theta);
///  * mapping onto a surface of revolution: otherwise x(m'), r(m') follow the layer
///    streamline r(z); here they come from a table of (m', x, r) at Gauss-integrated
///    m' with cubic Hermite interpolation using the exact slopes d(x, r)/dm' = r t.
#pragma once

#include "blade/Meridional.hpp"
#include "math/BSpline.hpp"
#include "math/Vec.hpp"

#include <string>
#include <vector>

namespace BladeGen::Blade {

using Math::Vec2;
using Math::Vec3;

enum class WallLaw { Linear, BellMehta };

[[nodiscard]] std::string to_string(WallLaw w);
[[nodiscard]] WallLaw wall_law_from_string(const std::string& s);

/// The Bell & Mehta fifth-order contour shape 6 xi^5 - 15 xi^4 + 10 xi^3 and its first
/// two derivatives (xi in [0, 1]).
[[nodiscard]] double bell_mehta(double xi);
[[nodiscard]] double bell_mehta_d1(double xi);
[[nodiscard]] double bell_mehta_d2(double xi);

/// LE and TE points of one wall in cylindrical coordinates (z axial, r radius).
struct WallPoints
{
    double z_le = -15.0, r_le = 40.0;
    double z_te = 15.0, r_te = 60.0;
};

struct MeridionalConfig
{
    WallPoints hub{-15.0, 40.0, 15.0, 43.0};
    WallPoints casing{-11.0, 106.0, 10.0, 102.0};
    WallLaw hub_law = WallLaw::BellMehta;
    WallLaw casing_law = WallLaw::BellMehta;
    double inlet_extension = 0.5;   ///< upstream length / the wall's axial LE -> TE length
    double outlet_extension = 0.5;  ///< downstream length / the wall's axial LE -> TE length

    /// Throws std::invalid_argument: z_TE > z_LE on both walls, r > 0 everywhere,
    /// casing above the hub, extensions in [0.05, 5].
    void validate() const;
};

/// Exact degree-5 B-spline of one wall in u in [0, 3] (see the file comment).
[[nodiscard]] Math::BSpline2 wall_spline(const WallPoints& w, WallLaw law, double inlet_ext, double outlet_ext);

/// The conformal map of one layer (surface of revolution), m' measured from the LE point.
class StreamSurface
{
public:
    StreamSurface() = default;
    StreamSurface(const Math::BSpline2& streamline, double span);

    [[nodiscard]] double span() const noexcept { return span_; }
    [[nodiscard]] bool cylindrical() const noexcept { return cylindrical_; }
    [[nodiscard]] double radius() const noexcept { return r0_; }      ///< cylindrical radius (when cylindrical)
    [[nodiscard]] double mprime_inlet() const noexcept { return mp_.front(); }  ///< m' of the inlet end (< 0)
    [[nodiscard]] double mprime_te() const noexcept { return mp_te_; }          ///< m' of the TE point
    [[nodiscard]] double mprime_outlet() const noexcept { return mp_.back(); }
    [[nodiscard]] double m_te() const noexcept { return m_te_; }                ///< meridional LE -> TE length
    [[nodiscard]] Vec2 le_point() const noexcept { return le_; }
    [[nodiscard]] Vec2 te_point() const noexcept { return te_; }

    /// (x, r) at m'; beyond the ends straight along the end tangent.
    [[nodiscard]] Vec2 xr(double mprime) const;
    /// m' -> meridional arc length m (from the LE point).
    [[nodiscard]] double m_of_mprime(double mprime) const;
    /// (x, y, z) = (x, r sin theta, r cos theta) of (m', theta [rad]).
    [[nodiscard]] Vec3 cartesian(double mprime, double theta) const;
    /// (x, r) samples of the whole streamline (inlet end -> outlet end).
    [[nodiscard]] std::vector<Vec2> polyline(int n_per_piece = 40) const;
    /// Signed curvature of the streamline in the (x, r) plane at m'.
    [[nodiscard]] double curvature(double mprime) const;

private:
    Math::BSpline2 c_;
    double span_ = 0.0;
    bool cylindrical_ = false;
    double r0_ = 0.0, x_le_ = 0.0;
    std::vector<double> u_, mp_, m_;
    std::vector<Vec2> xr_, dxr_;  // (x, r) and d(x, r)/dm' at the table points
    double mp_te_ = 0.0, m_te_ = 0.0;
    Vec2 le_{0.0, 0.0}, te_{0.0, 0.0};
};

class FlowPath
{
public:
    explicit FlowPath(MeridionalConfig cfg);

    [[nodiscard]] const MeridionalConfig& config() const noexcept { return cfg_; }
    [[nodiscard]] const Math::BSpline2& hub() const noexcept { return hub_; }
    [[nodiscard]] const Math::BSpline2& casing() const noexcept { return casing_; }
    /// Exact layer streamline at span fraction s (0 hub, 1 casing), u in [0, 3].
    [[nodiscard]] Math::BSpline2 streamline(double span) const;
    [[nodiscard]] StreamSurface surface(double span) const { return {streamline(span), span}; }
    /// The streamline as a MeridionalCurve (arc-length tables; for Blade / OCCT).
    [[nodiscard]] MeridionalCurve streamline_curve(double span) const { return MeridionalCurve(streamline(span)); }
    [[nodiscard]] Vec2 le_point(double span) const;
    [[nodiscard]] Vec2 te_point(double span) const;

private:
    MeridionalConfig cfg_;
    Math::BSpline2 hub_, casing_;
};

}  // namespace BladeGen::Blade
