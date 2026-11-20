/// @file profile/Profile2D.hpp
/// @brief Closed C2 blade profile in the unstaggered (m', theta) plane.
///
/// The closed curve is four B-spline pieces joined in one global parameter t:
///
///     piece            curve            global t                       from -> to
///     LeadingEdge      quintic, 2 spans [0, D_le]                      P_LE_PS -> P_LE -> P_LE_SS
///     SuctionSide      cubic            [D_le, D_le + 1]               P_LE_SS -> P_TE_SS
///     TrailingEdge     quintic, 2 spans [D_le + 1, D_le + 1 + D_te]    P_TE_SS -> P_TE -> P_TE_PS
///     PressureSide     cubic, reversed  [D_le + 1 + D_te, period]      P_TE_PS -> P_LE_PS
///
/// so the loop runs clockwise (LE up, SS left->right, TE down, PS right->left)
/// and is periodic with period = D_le + D_te + 2.  In ConnectorMode::C2 the
/// connectors copy the neighbours' derivatives (dC/dt, d2C/dt2), so the loop is
/// C2 in t everywhere; in G2Legacy it is curvature-continuous (G2) only.
#pragma once

#include "math/BSpline.hpp"
#include "profile/ProfileParams.hpp"

#include <array>
#include <string>
#include <vector>

namespace BladeGen::Profile {

using Math::Vec2;

/// Characteristic points of one edge (blade2d geometry/characteristicpoints.py).
struct EdgeGeometry
{
    Vec2 p_edge, p_asym, p_mid, p_ss, p_ps, v_axis, n_axis;
    double a = 0.0, b = 0.0;
};

/// Unit vector at angle beta from the m' axis, pointing right (x_sign = +1, LE)
/// or left (x_sign = -1, TE): (x_sign cos beta, sin beta).
[[nodiscard]] Vec2 unit_from_beta(double beta_deg, double x_sign);

[[nodiscard]] EdgeGeometry build_edge_geometry(const Vec2& p_edge, const EdgeParams& e, double x_sign);

/// The piecewise power thickness law at SS arc-length fraction xi (see ThicknessParams).
[[nodiscard]] double thickness_law(double xi, const ThicknessParams& t, double xi_max);

/// Throws std::invalid_argument when a parameter is outside its domain.
void validate(const ProfileParams& p);

enum class Piece : int { LeadingEdge = 0, SuctionSide = 1, TrailingEdge = 2, PressureSide = 3 };
enum class Side : int { Suction = 0, Pressure = 1 };

[[nodiscard]] std::string to_string(Piece p);

/// Junction report: continuity of the closed curve at a piece boundary.
struct Junction
{
    std::string name;
    double t = 0.0;              ///< global parameter
    double gap = 0.0;            ///< |C(t-) - C(t+)|
    double tangent_angle = 0.0;  ///< angle between the unit tangents [rad]
    double curvature_jump = 0.0; ///< |kappa(t-) - kappa(t+)|
    double d1_jump = 0.0;        ///< |C'(t-) - C'(t+)| / |C'|      (parametric C1)
    double d2_jump = 0.0;        ///< |C''(t-) - C''(t+)| / |C'|^2  (parametric C2)
};

/// Points along a curve with their arc length and curvature.
struct CurveSample
{
    std::vector<Vec2> xy;
    std::vector<double> t;       ///< global loop parameter
    std::vector<double> s;       ///< arc length from the first point
    std::vector<double> kappa;   ///< convexity-signed curvature (> 0: bulges outward)
    std::vector<int> piece;      ///< Piece of every point
};

/// Measured normal thickness: rays from the suction side along its inward normal.
struct ThicknessSample
{
    std::vector<double> xi;       ///< SS arc-length fraction
    std::vector<double> mprime;   ///< m' of the SS point
    std::vector<double> measured; ///< distance to the pressure side (NaN when no hit)
    std::vector<double> law;      ///< thickness law t(xi)
};

class Profile2D
{
public:
    explicit Profile2D(ProfileParams params);

    [[nodiscard]] const ProfileParams& params() const noexcept { return params_; }

    // ── characteristic points ────────────────────────────────────────────────
    [[nodiscard]] const EdgeGeometry& le_geometry() const noexcept { return le_; }
    [[nodiscard]] const EdgeGeometry& te_geometry() const noexcept { return te_; }
    [[nodiscard]] Vec2 p_le() const { return {0.0, 0.0}; }
    [[nodiscard]] Vec2 p_te() const { return {1.0, 0.0}; }
    [[nodiscard]] const Vec2& p_le_ss() const noexcept { return le_.p_ss; }
    [[nodiscard]] const Vec2& p_le_ps() const noexcept { return le_.p_ps; }
    [[nodiscard]] const Vec2& p_te_ss() const noexcept { return te_.p_ss; }
    [[nodiscard]] const Vec2& p_te_ps() const noexcept { return te_.p_ps; }

    // ── the four pieces, each in its own parameter ───────────────────────────
    [[nodiscard]] const Math::BSpline2& suction_side() const noexcept { return ss_; }   ///< LE -> TE, u in [0,1]
    [[nodiscard]] const Math::BSpline2& pressure_side() const noexcept { return ps_; }  ///< LE -> TE, u in [0,1]
    [[nodiscard]] const Math::BSpline2& le_connector() const noexcept { return le_arc_; } ///< P_LE_PS -> P_LE_SS
    [[nodiscard]] const Math::BSpline2& te_connector() const noexcept { return te_arc_; } ///< P_TE_SS -> P_TE_PS
    [[nodiscard]] double le_u_mid() const noexcept { return le_umid_; }
    [[nodiscard]] double te_u_mid() const noexcept { return te_umid_; }

    // thickness coupling, for plots and checks
    [[nodiscard]] double ss_length() const noexcept { return ss_len_; }
    [[nodiscard]] double xi_max() const noexcept { return xi_max_; }
    [[nodiscard]] const std::vector<double>& xi_ctrl() const noexcept { return xi_ctrl_; }
    [[nodiscard]] const std::vector<Vec2>& ps_free_points() const noexcept { return ps_free_; }

    // ── the closed curve in the global parameter ─────────────────────────────
    [[nodiscard]] double period() const noexcept { return breaks_[4]; }
    /// Global t at the piece starts: [0] LE arc, [1] SS, [2] TE arc, [3] PS, [4] = period.
    [[nodiscard]] const std::array<double, 5>& breaks() const noexcept { return breaks_; }
    [[nodiscard]] double t_le_point() const noexcept { return breaks_[0] + le_umid_ * dle_; }
    [[nodiscard]] double t_te_point() const noexcept { return breaks_[2] + te_umid_ * dte_; }
    [[nodiscard]] Piece piece_at(double t) const;
    /// C(t), C'(t), C''(t), C'''(t) in the global parameter (t wraps modulo the period).
    /// `from_left` picks the piece ending at t when t is a break.
    [[nodiscard]] std::array<Vec2, 4> derivs(double t, bool from_left = false) const;
    [[nodiscard]] Vec2 point(double t) const { return derivs(t)[0]; }
    /// Convexity-signed curvature (> 0 where the profile bulges outward).
    [[nodiscard]] double curvature(double t, bool from_left = false) const;

    [[nodiscard]] std::array<Junction, 4> junctions() const;

    // ── sampling ─────────────────────────────────────────────────────────────
    /// One side from the LE point P_LE to the TE point P_TE: half LE connector,
    /// SS (or PS), half TE connector.  n_edge points per half connector, n_body
    /// on the side (cosine-clustered toward both ends).
    [[nodiscard]] CurveSample sample_side(Side side, int n_edge = 60, int n_body = 240) const;
    /// The closed loop from P_LE over the suction side, P_TE and back over the
    /// pressure side; the last point repeats the first.  Every profile built with
    /// the same counts has the same point layout (corresponding points: LE, TE and
    /// the four junctions sit at fixed indices), which the 3-D skinning relies on.
    [[nodiscard]] CurveSample sample_loop(int n_edge = 40, int n_body = 160) const;
    /// Measured normal thickness at n SS stations, with the thickness law.
    [[nodiscard]] ThicknessSample thickness(int n = 200) const;
    /// m' -> suction-side parameter u (first crossing; clamped to the SS range).
    [[nodiscard]] double ss_u_at_mprime(double mprime) const;

private:
    void build();
    [[nodiscard]] Math::BSpline2 connector(const Vec2& p0, const Vec2& pm, const Vec2& p1,
                                           const Vec2& d0, const Vec2& dd0, const Vec2& d1,
                                           const Vec2& dd1, double& u_mid) const;
    [[nodiscard]] Math::BSpline2 fair_connector(const Vec2& p0, const Vec2& pm, const Vec2& p1,
                                                const Vec2& t0, double k0, const Vec2& t1, double k1,
                                                double speed, double& u_mid) const;
    void append_range(CurveSample& out, double ta, double tb, int n, bool cluster_start,
                      bool cluster_end, bool skip_first) const;

    ProfileParams params_;
    EdgeGeometry le_{}, te_{};
    Math::BSpline2 ss_, ps_, le_arc_, te_arc_;
    double le_umid_ = 0.5, te_umid_ = 0.5;
    double dle_ = 0.0, dte_ = 0.0;
    std::array<double, 5> breaks_{};
    double ss_len_ = 0.0, xi_max_ = 0.0;
    std::vector<double> xi_ctrl_;
    std::vector<Vec2> ps_free_;
};

}  // namespace BladeGen::Profile
