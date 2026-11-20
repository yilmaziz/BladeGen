#include "profile/Profile2D.hpp"

#include "math/FairCurve.hpp"
#include "math/Interpolation.hpp"

#include <tuple>
#include "math/Numerics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BladeGen::Profile {

using Math::BSpline2;
using Math::Constraint;
using Math::cross2;
using Math::deg2rad;
using Math::left_normal;
using Math::right_normal;
using Math::unit;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

void require(bool ok, const std::string& msg)
{
    if (!ok) throw std::invalid_argument("ProfileParams: " + msg);
}

bool finite_all(const ProfileParams& p)
{
    for (double v : to_vector(p))
        if (!std::isfinite(v)) return false;
    for (double v : p.ps.mprime_ctrl)
        if (!std::isfinite(v)) return false;
    return true;
}

}  // namespace

std::string to_string(Piece p)
{
    switch (p) {
        case Piece::LeadingEdge: return "LE";
        case Piece::SuctionSide: return "SS";
        case Piece::TrailingEdge: return "TE";
        case Piece::PressureSide: return "PS";
    }
    return "?";
}

Vec2 unit_from_beta(double beta_deg, double x_sign)
{
    const double b = deg2rad(beta_deg);
    return unit<2>(Vec2(x_sign * std::cos(b), std::sin(b)));
}

EdgeGeometry build_edge_geometry(const Vec2& p_edge, const EdgeParams& e, double x_sign)
{
    EdgeGeometry g;
    g.p_edge = p_edge;
    g.b = e.radius;
    g.a = e.stretch_ab * e.radius;
    g.v_axis = unit_from_beta(e.beta_deg, x_sign);
    g.n_axis = left_normal(g.v_axis);
    g.p_asym = p_edge + e.asym * g.n_axis;
    g.p_mid = g.p_asym + g.a * g.v_axis;
    const Vec2 p1 = g.p_mid + g.b * g.n_axis;
    const Vec2 p2 = g.p_mid - g.b * g.n_axis;
    // the higher-theta point is on the suction side
    if (p1.y() >= p2.y()) { g.p_ss = p1; g.p_ps = p2; }
    else { g.p_ss = p2; g.p_ps = p1; }
    return g;
}

double thickness_law(double xi, const ThicknessParams& t, double xi_max)
{
    xi = std::clamp(xi, 0.0, 1.0);
    xi_max = std::clamp(xi_max, 1e-6, 1.0 - 1e-6);
    const double shp_le = std::max(t.shape_le, 1e-6);
    double shp_te = std::max(t.shape_te, 1e-6);
    if (t.enforce_c1) {
        // dt/dxi continuous at the peak:
        //   shape_le (t_max - t_le) / xi_max = shape_te (t_max - t_te) / (1 - xi_max)
        const double num = shp_le * (t.t_max - t.t_le) * (1.0 - xi_max);
        const double den = (t.t_max - t.t_te) * xi_max;
        shp_te = std::max(den > 1e-14 ? num / den : t.shape_te, 1e-6);
    }
    if (xi <= xi_max) return t.t_le + (t.t_max - t.t_le) * std::pow(xi / xi_max, shp_le);
    return t.t_te + (t.t_max - t.t_te) * std::pow((1.0 - xi) / (1.0 - xi_max), shp_te);
}

void validate(const ProfileParams& p)
{
    require(finite_all(p), "non-finite value");
    for (const auto* e : {&p.le, &p.te}) {
        const std::string w = (e == &p.le) ? "le." : "te.";
        require(e->radius > 0.0, w + "radius must be > 0");
        require(e->stretch_ab > 0.0, w + "stretch_ab must be > 0");
        require(e->tension > 0.0, w + "tension must be > 0");
        require(std::abs(e->beta_deg) < 89.0, w + "beta_deg must be within (-89, 89)");
    }
    require(p.ss.handle_le_frac > 0.0 && p.ss.handle_te_frac > 0.0, "ss handles must be > 0");
    require(p.ps.handle_le_frac > 0.0 && p.ps.handle_te_frac > 0.0, "ps handles must be > 0");
    require(p.thickness.shape_le > 0.0 && p.thickness.shape_te > 0.0, "thickness exponents must be > 0");
    require(p.thickness.mprime_max > 0.0 && p.thickness.mprime_max < 1.0, "thickness.mprime_max must be in (0, 1)");
    require(!p.ps.mprime_ctrl.empty(), "ps.mprime_ctrl needs at least one station");
    for (std::size_t i = 0; i < p.ps.mprime_ctrl.size(); ++i) {
        require(p.ps.mprime_ctrl[i] > 0.0 && p.ps.mprime_ctrl[i] < 1.0, "ps.mprime_ctrl stations must be in (0, 1)");
        if (i > 0) require(p.ps.mprime_ctrl[i] > p.ps.mprime_ctrl[i - 1], "ps.mprime_ctrl must increase");
    }
}

// ─────────────────────────────────────────────────────────────────────────────

Profile2D::Profile2D(ProfileParams params) : params_(std::move(params))
{
    validate(params_);
    build();
}

double Profile2D::ss_u_at_mprime(double mprime) const
{
    const double x0 = ss_.eval(0.0).x(), x1 = ss_.eval(1.0).x();
    const double lo = std::min(x0, x1), hi = std::max(x0, x1);
    const double m = std::clamp(mprime, lo, hi);
    auto f = [&](double u) { return ss_.eval(u).x() - m; };
    const double u = Math::first_root(f, 0.0, 1.0, 256);
    if (std::isfinite(u)) return u;
    double best = 0.0, err = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= 2048; ++i) {
        const double uu = i / 2048.0, e = std::abs(f(uu));
        if (e < err) { err = e; best = uu; }
    }
    return best;
}

BSpline2 Profile2D::connector(const Vec2& p0, const Vec2& pm, const Vec2& p1, const Vec2& d0,
                              const Vec2& dd0, const Vec2& d1, const Vec2& dd1, double& u_mid) const
{
    const double L0 = (pm - p0).norm(), L1 = (p1 - pm).norm();
    if (!(L0 > 1e-12) || !(L1 > 1e-12))
        throw std::invalid_argument("Profile2D: degenerate edge (edge point coincides with a side end)");
    u_mid = L0 / (L0 + L1);
    // Two-span quintic, knots {0^6, u_mid, 1^6}: 7 control points for 7 conditions
    // (position at 0, u_mid, 1; first and second derivatives at both ends).
    std::vector<double> U = {0, 0, 0, 0, 0, 0, u_mid, 1, 1, 1, 1, 1, 1};
    std::vector<Constraint<2>> c = {
        {0.0, 0, p0}, {0.0, 1, d0}, {0.0, 2, dd0}, {u_mid, 0, pm},
        {1.0, 0, p1}, {1.0, 1, d1}, {1.0, 2, dd1},
    };
    return Math::solve_bspline<2>(5, std::move(U), c);
}

BSpline2 Profile2D::fair_connector(const Vec2& p0, const Vec2& pm, const Vec2& p1, const Vec2& t0,
                                   double k0, const Vec2& t1, double k1, double speed, double& u_mid) const
{
    const double L0 = (pm - p0).norm(), L1 = (p1 - pm).norm();
    if (!(L0 > 1e-12) || !(L1 > 1e-12))
        throw std::invalid_argument("Profile2D: degenerate edge (edge point coincides with a side end)");
    u_mid = L0 / (L0 + L1);
    // degree 5, 4 uniform interior knots: 10 control points (20 unknowns) for 12
    // conditions; the 8 remaining degrees of freedom minimise the jerk energy
    std::vector<double> U = {0, 0, 0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1, 1, 1};
    const Vec2 ex(1.0, 0.0), ey(0.0, 1.0);
    std::vector<Math::LinearCondition> c;
    for (const auto& [u, q] : {std::pair{0.0, p0}, std::pair{u_mid, pm}, std::pair{1.0, p1}}) {
        c.push_back({u, 0, ex, q.x()});
        c.push_back({u, 0, ey, q.y()});
    }
    for (const auto& [u, t, k] : {std::tuple{0.0, t0, k0}, std::tuple{1.0, t1, k1}}) {
        c.push_back({u, 1, ex, speed * t.x()});  // C' = speed * t
        c.push_back({u, 1, ey, speed * t.y()});
        c.push_back({u, 2, left_normal(t), k * speed * speed});  // n . C'' = kappa |C'|^2
    }
    return Math::fair_bspline(5, std::move(U), c, 3);
}

void Profile2D::build()
{
    const auto& p = params_;
    le_ = build_edge_geometry(p_le(), p.le, +1.0);
    te_ = build_edge_geometry(p_te(), p.te, -1.0);

    // ── suction side ─────────────────────────────────────────────────────────
    const Vec2 v_le_ss = unit_from_beta(p.le.beta_deg + p.le.d_beta_ss_deg, +1.0);
    const Vec2 v_te_ss = unit_from_beta(p.te.beta_deg + p.te.d_beta_ss_deg, -1.0);
    const double span = (te_.p_ss - le_.p_ss).norm();
    if (!(span > 1e-12)) throw std::invalid_argument("Profile2D: P_LE_SS and P_TE_SS coincide");
    ss_ = BSpline2::clamped(3, {le_.p_ss, le_.p_ss + p.ss.handle_le_frac * span * v_le_ss,
                                Vec2(p.ss.ctrl2_mprime, p.ss.ctrl2_theta),
                                te_.p_ss + p.ss.handle_te_frac * span * v_te_ss, te_.p_ss});

    // ── thickness coupling on SS arc length ──────────────────────────────────
    ss_len_ = Math::arc_length(ss_, 0.0, 1.0);
    auto xi_of_u = [&](double u) { return Math::arc_length(ss_, 0.0, u) / ss_len_; };
    xi_max_ = xi_of_u(ss_u_at_mprime(p.thickness.mprime_max));

    // ── pressure side ────────────────────────────────────────────────────────
    const Vec2 v_le_ps = unit_from_beta(p.le.beta_deg - p.le.d_beta_ps_deg, +1.0);
    const Vec2 v_te_ps = unit_from_beta(p.te.beta_deg - p.te.d_beta_ps_deg, -1.0);
    const double span_ps = (te_.p_ps - le_.p_ps).norm();
    if (!(span_ps > 1e-12)) throw std::invalid_argument("Profile2D: P_LE_PS and P_TE_PS coincide");
    std::vector<Vec2> ctrl = {le_.p_ps, le_.p_ps + p.ps.handle_le_frac * span_ps * v_le_ps};
    xi_ctrl_.clear();
    ps_free_.clear();
    for (double mp : p.ps.mprime_ctrl) {
        const double u = ss_u_at_mprime(mp);
        const double xi = std::clamp(xi_of_u(u), 1e-6, 1.0 - 1e-6);
        const auto d = ss_.derivs(u, 1);
        const Vec2 n = right_normal(unit<2>(d[1]));  // toward the pressure side
        const Vec2 q = d[0] + thickness_law(xi, p.thickness, xi_max_) * n;
        xi_ctrl_.push_back(xi);
        ps_free_.push_back(q);
        ctrl.push_back(q);
    }
    ctrl.push_back(te_.p_ps + p.ps.handle_te_frac * span_ps * v_te_ps);
    ctrl.push_back(te_.p_ps);
    ps_ = BSpline2::clamped(3, std::move(ctrl));

    // ── edge connectors ──────────────────────────────────────────────────────
    const auto S0 = ss_.derivs(0.0, 2), S1 = ss_.derivs(1.0, 2);
    const auto P0 = ps_.derivs(0.0, 2), P1 = ps_.derivs(1.0, 2);
    for (const auto* d : {&S0, &S1, &P0, &P1})
        if (!((*d)[1].norm() > 1e-12)) throw std::invalid_argument("Profile2D: zero end speed of a side");

    const double Lle = (p_le() - le_.p_ps).norm() + (le_.p_ss - p_le()).norm();
    const double Lte = (p_te() - te_.p_ss).norm() + (te_.p_ps - p_te()).norm();
    // parameter spans of the connectors in the global parameter: the tension is the
    // connector speed relative to its chord length, measured against the
    // geometric mean of the neighbours' speeds.
    dle_ = p.le.tension * Lle / std::sqrt(P0[1].norm() * S0[1].norm());
    dte_ = p.te.tension * Lte / std::sqrt(S1[1].norm() * P1[1].norm());

    const double kS0 = Math::signed_curvature(S0[1], S0[2]), kS1 = Math::signed_curvature(S1[1], S1[2]);
    const double kP0 = Math::signed_curvature(P0[1], P0[2]), kP1 = Math::signed_curvature(P1[1], P1[2]);
    if (p.connector == ConnectorMode::Fair) {
        // G2 at both junctions (unit tangent and curvature of the neighbour; the
        // reversed PS has the opposite tangent and curvature), the edge point
        // interpolated, minimum jerk energy for the remaining freedom.
        le_arc_ = fair_connector(le_.p_ps, p_le(), le_.p_ss, -unit<2>(P0[1]), -kP0, unit<2>(S0[1]), kS0,
                                 p.le.tension * Lle, le_umid_);
        te_arc_ = fair_connector(te_.p_ss, p_te(), te_.p_ps, unit<2>(S1[1]), kS1, -unit<2>(P1[1]), -kP1,
                                 p.te.tension * Lte, te_umid_);
    } else if (p.connector == ConnectorMode::C2Parametric) {
        // LE: w = 0 meets the reversed PS (dC/dt = -P'(0), d2C/dt2 = P''(0)),
        //     w = 1 meets the SS (S'(0), S''(0)); dC/dw = D dC/dt, d2C/dw2 = D^2 d2C/dt2.
        le_arc_ = connector(le_.p_ps, p_le(), le_.p_ss, -dle_ * P0[1], dle_ * dle_ * P0[2],
                            dle_ * S0[1], dle_ * dle_ * S0[2], le_umid_);
        te_arc_ = connector(te_.p_ss, p_te(), te_.p_ps, dte_ * S1[1], dte_ * dte_ * S1[2],
                            -dte_ * P1[1], dte_ * dte_ * P1[2], te_umid_);
    } else {
        // blade2d profile_c2_connect.py: |r'| = alpha * L_total along the matched unit
        // tangent, r'' = kappa (alpha L)^2 n (no tangential part).
        auto ends = [](const Vec2& t_dir, double kappa, double speed) {
            return std::pair<Vec2, Vec2>{speed * t_dir, kappa * speed * speed * left_normal(t_dir)};
        };
        const auto [a0, aa0] = ends(-unit<2>(P0[1]), -kP0, p.le.tension * Lle);
        const auto [a1, aa1] = ends(unit<2>(S0[1]), kS0, p.le.tension * Lle);
        le_arc_ = connector(le_.p_ps, p_le(), le_.p_ss, a0, aa0, a1, aa1, le_umid_);
        const auto [b0, bb0] = ends(unit<2>(S1[1]), kS1, p.te.tension * Lte);
        const auto [b1, bb1] = ends(-unit<2>(P1[1]), -kP1, p.te.tension * Lte);
        te_arc_ = connector(te_.p_ss, p_te(), te_.p_ps, b0, bb0, b1, bb1, te_umid_);
    }
    breaks_ = {0.0, dle_, dle_ + 1.0, dle_ + 1.0 + dte_, dle_ + dte_ + 2.0};
}

Piece Profile2D::piece_at(double t) const
{
    const double T = period();
    t = std::fmod(t, T);
    if (t < 0.0) t += T;
    for (int k = 3; k >= 0; --k)
        if (t >= breaks_[static_cast<std::size_t>(k)]) return static_cast<Piece>(k);
    return Piece::LeadingEdge;
}

std::array<Vec2, 4> Profile2D::derivs(double t, bool from_left) const
{
    const double T = period();
    if (t < 0.0 || t > T) {
        t = std::fmod(t, T);
        if (t < 0.0) t += T;
    }
    if (from_left && t == 0.0) t = T;
    int k = 3;
    for (int j = 0; j < 4; ++j) {
        const double a = breaks_[static_cast<std::size_t>(j)], b = breaks_[static_cast<std::size_t>(j + 1)];
        const bool in = from_left ? (t > a && t <= b) : (t >= a && t < b);
        if (in) { k = j; break; }
    }
    if (t >= T) k = 3;
    const double a = breaks_[static_cast<std::size_t>(k)];
    const double L = breaks_[static_cast<std::size_t>(k + 1)] - a;
    const BSpline2* c = nullptr;
    switch (k) {
        case 0: c = &le_arc_; break;
        case 1: c = &ss_; break;
        case 2: c = &te_arc_; break;
        default: c = &ps_; break;
    }
    const bool rev = (k == 3);
    const double w = std::clamp((t - a) / L, 0.0, 1.0);
    const auto d = c->derivs(rev ? 1.0 - w : w, 3);
    const double s = rev ? -1.0 / L : 1.0 / L;
    return {d[0], s * d[1], s * s * d[2], s * s * s * d[3]};
}

double Profile2D::curvature(double t, bool from_left) const
{
    const auto d = derivs(t, from_left);
    return -Math::signed_curvature(d[1], d[2]);  // loop is clockwise: convex = negative signed
}

std::array<Junction, 4> Profile2D::junctions() const
{
    const std::array<std::string, 4> names = {"PS -> LE (P_LE_PS)", "LE -> SS (P_LE_SS)",
                                              "SS -> TE (P_TE_SS)", "TE -> PS (P_TE_PS)"};
    std::array<Junction, 4> out;
    for (std::size_t k = 0; k < 4; ++k) {
        const double t = breaks_[k];
        const auto L = derivs(t, true), R = derivs(t, false);
        Junction& j = out[k];
        j.name = names[k];
        j.t = t;
        j.gap = (L[0] - R[0]).norm();
        const Vec2 a = unit<2>(L[1]), b = unit<2>(R[1]);
        j.tangent_angle = std::atan2(std::abs(cross2(a, b)), a.dot(b));
        j.curvature_jump = std::abs(Math::signed_curvature(L[1], L[2]) - Math::signed_curvature(R[1], R[2]));
        const double sp = R[1].norm();
        j.d1_jump = (L[1] - R[1]).norm() / sp;
        j.d2_jump = (L[2] - R[2]).norm() / (sp * sp);
    }
    return out;
}

void Profile2D::append_range(CurveSample& out, double ta, double tb, int n, bool cluster_start,
                             bool cluster_end, bool skip_first) const
{
    if (n < 2) throw std::invalid_argument("Profile2D sampling: need >= 2 points per piece");
    for (int i = skip_first ? 1 : 0; i < n; ++i) {
        double s = static_cast<double>(i) / (n - 1);
        if (cluster_start && cluster_end) s = 0.5 * (1.0 - std::cos(Math::kPi * s));
        else if (cluster_start) s = 1.0 - std::cos(0.5 * Math::kPi * s);
        else if (cluster_end) s = std::sin(0.5 * Math::kPi * s);
        const double t = ta + (tb - ta) * s;
        // at a break, evaluate on the piece that lies inside [min(ta,tb), max(ta,tb)]:
        // the upper end of the range is reached from the left
        const bool from_left = (i == 0) ? (ta > tb) : (i == n - 1 ? (tb > ta) : false);
        const auto d = derivs(t, from_left);
        out.xy.push_back(d[0]);
        out.t.push_back(t);
        out.kappa.push_back(-Math::signed_curvature(d[1], d[2]));
        out.piece.push_back(static_cast<int>(piece_at(from_left ? t - 1e-12 * period() : t)));
    }
}

static void finish_arclength(CurveSample& c)
{
    c.s.assign(c.xy.size(), 0.0);
    for (std::size_t i = 1; i < c.xy.size(); ++i) c.s[i] = c.s[i - 1] + (c.xy[i] - c.xy[i - 1]).norm();
}

CurveSample Profile2D::sample_side(Side side, int n_edge, int n_body) const
{
    CurveSample c;
    const auto& B = breaks_;
    if (side == Side::Suction) {
        append_range(c, t_le_point(), B[1], n_edge, false, true, false);
        append_range(c, B[1], B[2], n_body, true, true, true);
        append_range(c, B[2], t_te_point(), n_edge, true, false, true);
    } else {
        append_range(c, t_le_point(), B[0], n_edge, false, true, false);
        append_range(c, B[4], B[3], n_body, true, true, true);
        append_range(c, B[3], t_te_point(), n_edge, true, false, true);
    }
    finish_arclength(c);
    return c;
}

CurveSample Profile2D::sample_loop(int n_edge, int n_body) const
{
    CurveSample c;
    const auto& B = breaks_;
    append_range(c, t_le_point(), B[1], n_edge, false, true, false);
    append_range(c, B[1], B[2], n_body, true, true, true);
    append_range(c, B[2], t_te_point(), n_edge, true, false, true);
    append_range(c, t_te_point(), B[3], n_edge, false, true, true);
    append_range(c, B[3], B[4], n_body, true, true, true);
    append_range(c, B[0], t_le_point(), n_edge, true, false, true);
    finish_arclength(c);
    return c;
}

ThicknessSample Profile2D::thickness(int n) const
{
    if (n < 2) throw std::invalid_argument("Profile2D::thickness: need n >= 2");
    // pressure-side polyline
    std::vector<Vec2> ps;
    const int m = 800;
    for (int i = 0; i <= m; ++i) ps.push_back(ps_.eval(static_cast<double>(i) / m));
    ThicknessSample out;
    for (int i = 0; i < n; ++i) {
        const double u = static_cast<double>(i) / (n - 1);
        const auto d = ss_.derivs(u, 1);
        const Vec2 nrm = right_normal(unit<2>(d[1]));
        const double xi = Math::arc_length(ss_, 0.0, u) / ss_len_;
        double best = kNaN;
        for (std::size_t k = 0; k + 1 < ps.size(); ++k) {
            const Vec2 A = ps[k], AB = ps[k + 1] - ps[k], AP = A - d[0];
            const double den = cross2(nrm, AB);
            if (std::abs(den) < 1e-14) continue;
            const double tt = cross2(AP, AB) / den;   // along the ray
            const double uu = cross2(AP, nrm) / den;  // along the segment
            if (tt > 1e-9 && uu >= 0.0 && uu <= 1.0 && !(tt >= best)) best = tt;
        }
        out.xi.push_back(xi);
        out.mprime.push_back(d[0].x());
        out.measured.push_back(best);
        out.law.push_back(thickness_law(xi, params_.thickness, xi_max_));
    }
    return out;
}

}  // namespace BladeGen::Profile
