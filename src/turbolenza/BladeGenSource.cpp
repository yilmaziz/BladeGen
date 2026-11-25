// The BladeGen adapter: a design layer's profile -> TURBOlenza's BladeSection (see BladeGenSource.hpp).
#include "turbolenza/BladeGenSource.hpp"

#include "math/Vec.hpp"
#include "naca/Naca65.hpp"
#include "profile/Profile2D.hpp"
#include "blade/FlowPath.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace BladeGen::Turbolenza {

namespace {

constexpr int kDegree = 5;
constexpr double kEndTolerance = 1e-9;  // BladeGen's sides start at (0, 0) and end at (1, 0) exactly
constexpr double kMinSpan = 1e-8;       // the narrowest knot span the refinement makes (in u)

std::size_t sz(int k) { return static_cast<std::size_t>(k); }

std::string sci(double v)
{
    char b[32];
    std::snprintf(b, sizeof b, "%.3e", v);
    return b;
}

/// The p + 1 non-zero basis functions N_(l-p) .. N_l at u on span l (The NURBS Book, A2.2).
void basis_functions(const std::vector<double>& t, int l, double u, int p, double* N)
{
    double left[kDegree + 1], right[kDegree + 1];
    N[0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = u - t[sz(l + 1 - j)];
        right[j] = t[sz(l + j)] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            const double temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
}

/// The basis functions on span l and their first two derivatives: D[k][a] = d^k N_(l-p+a) / du^k
/// (The NURBS Book, A2.3, for order 2; fixed-size arrays: the fit evaluates every sample in every
/// round).
void basis_derivatives(const std::vector<double>& t, int l, double u, double D[3][kDegree + 1])
{
    constexpr int p = kDegree;
    double ndu[p + 1][p + 1], left[p + 1], right[p + 1], a[2][p + 1];
    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = u - t[sz(l + 1 - j)];
        right[j] = t[sz(l + j)] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    for (int j = 0; j <= p; ++j) D[0][j] = ndu[j][p];
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= 2; ++k) {
            double d = 0.0;
            const int rk = r - k, pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            const int j1 = rk >= -1 ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            D[k][r] = d;
            std::swap(s1, s2);
        }
    }
    for (int j = 0; j <= p; ++j) {
        D[1][j] *= p;
        D[2][j] *= p * (p - 1);
    }
}

/// The span l with t_l <= u < t_(l+1) (the last span at u = 1), as ParametricBladeCurve.
int find_span(const std::vector<double>& t, int n_ctrl, double u)
{
    if (u >= t[sz(n_ctrl)]) return n_ctrl - 1;
    const auto it = std::upper_bound(t.begin() + kDegree, t.begin() + n_ctrl + 1, u);
    return static_cast<int>(it - t.begin()) - 1;
}

/// Least-squares control points for the clamped quintic with the given interior knots.
ParametricBladeCurve solve_fit(const std::vector<Point>& P, const std::vector<double>& u, const std::vector<double>& interior,
                               double length)
{
    std::vector<double> t(sz(kDegree + 1), 0.0);
    t.insert(t.end(), interior.begin(), interior.end());
    t.insert(t.end(), sz(kDegree + 1), 1.0);
    const int n_ctrl = static_cast<int>(t.size()) - kDegree - 1;
    const int n_free = n_ctrl - 2;  // c_0 = (0, 0) and c_(n-1) = (1, 0) are fixed
    std::vector<Eigen::Triplet<double>> trip;
    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(n_free, 2);
    std::vector<double> M(sz(n_free * (2 * kDegree + 1)), 0.0);  // band of B^T B, row a, offset b - a + p
    double N[kDegree + 1];
    for (std::size_t k = 0; k < P.size(); ++k) {
        const int l = find_span(t, n_ctrl, u[k]);
        basis_functions(t, l, u[k], kDegree, N);
        // the fixed end control points move to the right-hand side
        double rm = P[k].m, ry = P[k].y;
        for (int a = 0; a <= kDegree; ++a) {
            const int ia = l - kDegree + a;
            if (ia == n_ctrl - 1) rm -= N[a] * 1.0;  // c_(n-1) = (1, 0); c_0 = (0, 0) contributes nothing
        }
        for (int a = 0; a <= kDegree; ++a) {
            const int ia = l - kDegree + a - 1;  // free index
            if (ia < 0 || ia >= n_free) continue;
            rhs(ia, 0) += N[a] * rm;
            rhs(ia, 1) += N[a] * ry;
            for (int b = 0; b <= kDegree; ++b) {
                const int ib = l - kDegree + b - 1;
                if (ib < 0 || ib >= n_free) continue;
                M[sz(ia * (2 * kDegree + 1) + (ib - ia + kDegree))] += N[a] * N[b];
            }
        }
    }
    for (int a = 0; a < n_free; ++a)
        for (int o = -kDegree; o <= kDegree; ++o) {
            const int b = a + o;
            if (b < 0 || b >= n_free) continue;
            const double v = M[sz(a * (2 * kDegree + 1) + (o + kDegree))];
            if (v != 0.0) trip.emplace_back(a, b, v);
        }
    Eigen::SparseMatrix<double> A(n_free, n_free);
    A.setFromTriplets(trip.begin(), trip.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(A);
    if (ldlt.info() != Eigen::Success)
        throw std::runtime_error("fit_side: the least-squares system is singular (a knot span without samples)");
    const Eigen::MatrixXd c = ldlt.solve(rhs);
    if (ldlt.info() != Eigen::Success || !c.allFinite()) throw std::runtime_error("fit_side: the least-squares solve failed");
    std::vector<Point> ctrl(sz(n_ctrl));
    ctrl.front() = {0.0, 0.0};
    ctrl.back() = {1.0, 0.0};
    for (int a = 0; a < n_free; ++a) ctrl[sz(a + 1)] = {c(a, 0), c(a, 1)};
    return {std::move(t), std::move(ctrl), length, kDegree};
}

SideSample to_side_sample(const Profile::CurveSample& c, Profile::Side side)
{
    // CurveSample::kappa is convexity-signed (> 0 bulging outward).  LE -> TE the suction side
    // (outward = +y, on the left) bulges outward when it turns clockwise, the pressure side
    // (outward on the right) when it turns counterclockwise.
    const double sign = side == Profile::Side::Suction ? -1.0 : 1.0;
    SideSample s;
    s.xy.reserve(c.xy.size());
    for (std::size_t k = 0; k < c.xy.size(); ++k) {
        s.xy.push_back({c.xy[k].x(), c.xy[k].y()});
        s.kappa.push_back(sign * c.kappa[k]);
    }
    return s;
}

}  // namespace

std::vector<double> surface_parameters(int n)
{
    if (n < 2) throw std::invalid_argument("surface_parameters: need n >= 2");
    std::vector<double> u(sz(n));
    for (int k = 0; k < n; ++k) u[sz(k)] = 0.5 * (1.0 - std::cos(Math::kPi * k / (n - 1)));
    u.front() = 0.0;
    u.back() = 1.0;
    return u;
}

ParametricBladeCurve fit_side(const SideSample& s, const SideFitOptions& opt, SideFitReport* report)
{
    const std::size_t n = s.xy.size();
    if (n < 100 || s.kappa.size() != n) throw std::invalid_argument("fit_side: need >= 100 samples with their curvature");
    if (std::hypot(s.xy.front().m, s.xy.front().y) > kEndTolerance || std::hypot(s.xy.back().m - 1.0, s.xy.back().y) > kEndTolerance)
        throw std::invalid_argument("fit_side: the side must run from the LE point (0, 0) to the TE point (1, 0)");
    if (!(opt.tolerance > 0.0) || !(opt.curvature_tolerance > 0.0) || opt.initial_spans < 2)
        throw std::invalid_argument("fit_side: tolerances must be > 0 and initial_spans >= 2");
    // chord-length parameter (exact curve points, dense: equal to the arc length to round-off)
    std::vector<Point> P(s.xy);
    P.front() = {0.0, 0.0};
    P.back() = {1.0, 0.0};
    std::vector<double> u(n, 0.0);
    for (std::size_t k = 1; k < n; ++k) u[k] = u[k - 1] + std::hypot(P[k].m - P[k - 1].m, P[k].y - P[k - 1].y);
    const double S = u.back();
    if (!(S > 1.0)) throw std::invalid_argument("fit_side: a side from (0, 0) to (1, 0) is longer than the chord");
    for (auto& v : u) v /= S;
    u.back() = 1.0;
    // starting knots, cosine-clustered toward both ends
    std::vector<double> interior;
    for (int i = 1; i < opt.initial_spans; ++i) interior.push_back(0.5 * (1.0 - std::cos(Math::kPi * i / opt.initial_spans)));
    auto first_at = [&](double v) { return static_cast<std::size_t>(std::lower_bound(u.begin(), u.end(), v) - u.begin()); };
    std::vector<double> kfit(n), err(n);
    int rounds = 0;
    for (;;) {
        ParametricBladeCurve curve = solve_fit(P, u, interior, S);
        std::vector<double> bounds = {0.0};
        bounds.insert(bounds.end(), interior.begin(), interior.end());
        bounds.push_back(1.0);
        const std::size_t n_span = bounds.size() - 1;
        {
            const auto& t = curve.knots();
            const auto& C = curve.control_points();
            const int n_ctrl = static_cast<int>(C.size());
            double D[3][kDegree + 1];
            for (std::size_t k = 0; k < n; ++k) {
                const int l = find_span(t, n_ctrl, u[k]);
                basis_derivatives(t, l, u[k], D);
                Point d[3];
                for (int q = 0; q < 3; ++q) {
                    double m = 0.0, y = 0.0;
                    for (int a = 0; a <= kDegree; ++a) {
                        m += D[q][a] * C[sz(l - kDegree + a)].m;
                        y += D[q][a] * C[sz(l - kDegree + a)].y;
                    }
                    d[q] = {m, y};
                }
                err[k] = std::hypot(d[0].m - P[k].m, d[0].y - P[k].y);
                kfit[k] = signed_curvature(d[1], d[2]);
            }
        }
        // per span: refine when a sample misses the position or the curvature tolerance; the
        // overshoot: the fit's curvature outside BladeGen's own range over the span and one span
        // width on either side (a wiggle), relative to max(|kappa|, 1)
        double res = 0.0, kerr = 0.0, krel = 0.0, over = 0.0;
        std::vector<double> next;
        bool split_any = false, missed = false;
        for (std::size_t a = 0; a < n_span; ++a) {
            const double lo = bounds[a], hi = bounds[a + 1], w = hi - lo;
            const std::size_t i0 = first_at(lo), i1 = a + 1 == n_span ? n : first_at(hi);
            const std::size_t w0 = first_at(lo - w), w1 = std::min(n, first_at(hi + w) + 1);
            double klo = s.kappa[w0], khi = s.kappa[w0];
            for (std::size_t k = w0; k < w1; ++k) {
                klo = std::min(klo, s.kappa[k]);
                khi = std::max(khi, s.kappa[k]);
            }
            bool refine = false, pos_miss = false;
            for (std::size_t k = i0; k < i1; ++k) {
                const double scale = std::max(std::abs(s.kappa[k]), 1.0);
                const double ke = std::abs(kfit[k] - s.kappa[k]);
                res = std::max(res, err[k]);
                kerr = std::max(kerr, ke);
                krel = std::max(krel, ke / scale);
                over = std::max(over, std::max({0.0, klo - kfit[k], kfit[k] - khi}) / scale);
                if (err[k] > opt.tolerance) pos_miss = true;
                if (err[k] > opt.tolerance || ke > opt.curvature_tolerance * scale) refine = true;
            }
            // a span is halved only while both halves keep enough samples for a well-posed fit
            // and stay wider than kMinSpan
            const double mid = 0.5 * (lo + hi);
            const std::size_t im = first_at(mid);
            const bool can = 0.5 * w >= kMinSpan && im - i0 >= sz(2 * (kDegree + 1)) && i1 - im >= sz(2 * (kDegree + 1));
            if (refine && can) {
                next.push_back(mid);
                split_any = true;
            } else if (pos_miss) {
                missed = true;
            }
            if (a + 1 < n_span) next.push_back(hi);
        }
        if (!split_any) {
            if (missed || res > opt.tolerance)
                throw std::runtime_error("fit_side: the position tolerance " + sci(opt.tolerance) +
                                         " cannot be met with the sample density (residual " + sci(res) + ")");
            if (report) {
                report->control_points = static_cast<int>(curve.control_points().size());
                report->refinements = rounds;
                report->samples = static_cast<int>(n);
                report->arc_length = S;
                report->max_residual = res;
                report->max_curvature_error = kerr;
                report->max_relative_curvature_error = krel;
                report->max_curvature_overshoot = over;
            }
            return curve;
        }
        std::sort(next.begin(), next.end());
        interior = std::move(next);
        if (static_cast<int>(interior.size()) + kDegree + 1 > opt.max_control_points)
            throw std::runtime_error("fit_side: more than " + std::to_string(opt.max_control_points) +
                                     " control points needed (residual " + sci(res) + ")");
        ++rounds;
    }
}

SideSample sample_layer_side(const Design::Layer& layer, Design::Evaluator& ev, Profile::Side side, int n_edge, int n_body)
{
    if (layer.uses_custom()) {
        const Profile::Profile2D prof(layer.custom);
        return to_side_sample(prof.sample_side(side, n_edge, n_body), side);
    }
    const auto naca = ev.naca(layer.naca);
    return to_side_sample(naca->sample_side(side, n_edge, n_body), side);
}

namespace {

/// The section's designation, blade angles chi_le / chi_te and edge metadata (see the file comment).
void section_metadata(const Design::Layer& layer, Design::Evaluator& ev, BladeSection& b, SectionFitReport& rep)
{
    if (layer.uses_custom()) {
        const Profile::Profile2D prof(layer.custom);
        b.designation = layer.method == Design::Method::Custom ? "Custom" : "NACA 65 + Custom (fitted)";
        // BladeGen's own LE / TE metal angles of the section (Section2D::camber, cached)
        const auto sec = ev.section(layer, Design::MeshOptions{});
        b.chi_le = Math::deg2rad(sec->camber.angle_deg.front());
        b.chi_te = -Math::deg2rad(sec->camber.angle_deg.back());
        b.le_radius = layer.custom.le.radius;
        b.te_radius = layer.custom.te.radius;
        const auto& le = prof.le_geometry().p_mid;
        const auto& te = prof.te_geometry().p_mid;
        b.le_center = {le.x(), le.y()};
        b.te_center = {te.x(), te.y()};
        rep.chi_source = "custom profile: BladeGen's camber-line end slopes (inscribed-circle centres, Section2D::camber)";
    } else {
        const auto naca = ev.naca(layer.naca);
        const auto& m = naca->metadata();
        b.designation = m.designation + " (" + Naca::to_string(layer.naca.variant) + ")";
        const auto arc = Naca::equivalent_arc_angles(layer.naca.cl0);
        b.chi_le = arc[1];
        b.chi_te = arc[2];
        b.le_radius = m.le_radius;
        b.te_radius = m.te_radius;
        b.le_center = {m.le_center.x(), m.le_center.y()};
        b.te_center = {m.te_center.x(), m.te_center.y()};
        rep.chi_source = "NACA 65: equivalent circular arc, tan(camber / 4) = 0.1103 cl0 (Aungier Eq. 4-6)";
    }
}

}  // namespace

BladeSection section_from_layer(const Design::Layer& layer, Design::Evaluator& ev, SectionFitReport* report, const SideFitOptions& opt)
{
    BladeSection b;
    SectionFitReport rep;
    // a side whose features are too fine for the default sample (e.g. the Kovach blend of a thin
    // NACA 65 at station 0.60) is sampled again, three times denser, once
    auto fit = [&](Profile::Side side, SideFitReport& r) {
        try {
            return fit_side(sample_layer_side(layer, ev, side, opt.n_edge, opt.n_body), opt, &r);
        } catch (const std::runtime_error&) {
            return fit_side(sample_layer_side(layer, ev, side, 3 * opt.n_edge, 3 * opt.n_body), opt, &r);
        }
    };
    b.suction = fit(Profile::Side::Suction, rep.suction);
    b.pressure = fit(Profile::Side::Pressure, rep.pressure);
    b.suction_parameters = surface_parameters(opt.n_parameters);
    b.pressure_parameters = b.suction_parameters;
    section_metadata(layer, ev, b, rep);
    if (report) *report = rep;
    return b;
}

// ── the blade on its stream surface (not in v1; see BladeGenSource.hpp) ───────

namespace {

// 8-point Gauss-Legendre on [-1, 1] (Abramowitz & Stegun, Table 25.4).
constexpr double kGx[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
constexpr double kGw[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
constexpr int kTable = 3 * 512;  // table intervals of u in [0, 3]

/// One layer streamline of BladeGen's FlowPath with its exact arc lengths m(u), m'(u) (from u = 0)
/// and the stream-sheet thickness b of the streamline family (see BladeGenSource.hpp).
class Streamline
{
public:
    Streamline(const Blade::FlowPath& fp, double span) : X_(fp.streamline(span)), hub_(fp.hub()), cas_(fp.casing())
    {
        h_ = 3.0 / kTable;
        m_.assign(sz(kTable + 1), 0.0);
        mp_.assign(sz(kTable + 1), 0.0);
        for (int k = 0; k < kTable; ++k) {
            const auto [dm, dmp] = arcs(h_ * k, h_ * (k + 1));
            m_[sz(k + 1)] = m_[sz(k)] + dm;
            mp_[sz(k + 1)] = mp_[sz(k)] + dmp;
        }
    }
    [[nodiscard]] double m_of_u(double u) const { return at(u, true); }
    [[nodiscard]] double mp_of_u(double u) const { return at(u, false); }
    /// u(m') and u(m) inside the channel, by Newton from the table (std::invalid_argument outside it).
    [[nodiscard]] double u_of_mp(double mp) const { return invert(mp, false); }
    [[nodiscard]] double u_of_m(double m) const { return invert(m, true); }
    [[nodiscard]] Math::Vec2 xr(double u) const { return X_.eval(u); }
    /// sin(phi) = dr/dm at u.
    [[nodiscard]] double sin_phi(double u) const
    {
        const auto d = X_.derivs(u, 1);
        return d[1].y() / d[1].norm();
    }

    struct Local
    {
        double r = 0.0, sin_phi = 0.0, b = 0.0;
    };
    /// r, sin(phi) = dr/dm and b at the arc length m from u = 0; beyond the channel's ends the
    /// streamline and its family continue straight along their end tangents.
    [[nodiscard]] Local at_m(double m) const
    {
        if (m >= 0.0 && m <= m_.back()) return local(u_of_m(m), 0.0);
        const bool before = m < 0.0;
        return local(before ? 0.0 : 3.0, before ? m : m - m_.back());
    }

private:
    [[nodiscard]] std::pair<double, double> arcs(double a, double b) const
    {
        const double c = 0.5 * (a + b), w = 0.5 * (b - a);
        double dm = 0.0, dmp = 0.0;
        for (int g = 0; g < 4; ++g)
            for (double sgn : {-1.0, 1.0}) {
                const auto d = X_.derivs(c + sgn * w * kGx[g], 1);
                const double v = d[1].norm();
                dm += kGw[g] * w * v;
                dmp += kGw[g] * w * v / d[0].y();
            }
        return {dm, dmp};
    }
    [[nodiscard]] double at(double u, bool arc) const
    {
        const int k = std::clamp(static_cast<int>(std::floor(u / h_)), 0, kTable - 1);
        const auto [dm, dmp] = arcs(h_ * k, u);
        return arc ? m_[sz(k)] + dm : mp_[sz(k)] + dmp;
    }
    [[nodiscard]] double invert(double target, bool arc) const
    {
        const std::vector<double>& t = arc ? m_ : mp_;
        if (!(target >= t.front() && target <= t.back()))
            throw std::invalid_argument("the blade leaves the BladeGen channel (its m' lies outside the layer streamline)");
        const auto it = std::upper_bound(t.begin(), t.end(), target);
        const int k = std::clamp(static_cast<int>(it - t.begin()) - 1, 0, kTable - 1);
        double u = h_ * (k + (target - t[sz(k)]) / (t[sz(k + 1)] - t[sz(k)]));
        for (int n = 0; n < 20; ++n) {
            const auto d = X_.derivs(u, 1);
            const double slope = arc ? d[1].norm() : d[1].norm() / d[0].y();
            const double step = (at(u, arc) - target) / slope;
            u -= step;
            if (std::abs(step) <= 1e-15 * std::max(1.0, std::abs(u))) break;
        }
        return u;
    }
    [[nodiscard]] Local local(double u, double beyond) const
    {
        const auto d = X_.derivs(u, 1);
        const Math::Vec2 t = d[1] / d[1].norm();
        const Math::Vec2 n(-t.y(), t.x());
        const auto dh = hub_.derivs(u, 1), dc = cas_.derivs(u, 1);
        Math::Vec2 ds = dc[0] - dh[0];  // dX/ds of X(u, s) = (1 - s) hub(u) + s casing(u)
        if (beyond != 0.0) ds += (beyond / d[1].norm()) * (dc[1] - dh[1]);
        Local l;
        l.r = d[0].y() + beyond * t.y();
        l.sin_phi = t.y();
        l.b = std::abs(ds.dot(n));
        return l;
    }

    Math::BSpline2 X_, hub_, cas_;
    double h_ = 0.0;
    std::vector<double> m_, mp_;
};

}  // namespace

double mapped_curvature(const Point& t, double kappa, double gamma, double c, double r, double dr_dmprime, double r_le,
                        double L)
{
    const double cg = std::cos(gamma), sg = std::sin(gamma);
    const Point n{-t.y, t.m};  // the chord-frame unit normal (left of t)
    // d(m', theta)/ds and d2(m', theta)/ds2 of the similarity (s the chord-frame arc length)
    const double a1 = c * (t.m * cg - t.y * sg), a2 = c * (t.m * sg + t.y * cg);
    const double kc = kappa * c;
    const double b1 = kc * (n.m * cg - n.y * sg), b2 = kc * (n.m * sg + n.y * cg);
    const double M1 = r * a1, M2 = r * b1 + dr_dmprime * a1 * a1;
    const double Y1 = r_le * a2, Y2 = r_le * b2;
    return L * (M1 * Y2 - Y1 * M2) / std::pow(M1 * M1 + Y1 * Y1, 1.5);
}

StreamSurfaceBlade stream_surface_blade(const Design::BladeDesign& design, int k, const Design::Layer& layer,
                                        Design::Evaluator& ev, const StreamSurfaceOptions& opt, StreamSurfaceReport* report)
{
    if (k < 0 || sz(k) >= design.layers.size())
        throw std::invalid_argument("stream_surface_blade: no layer " + std::to_string(k) + " in the design");
    if (design.mirror)
        throw std::invalid_argument("stream_surface_blade: a mirrored design (suction side toward -theta) is not "
                                    "supported: TURBOlenza's passage lies on the +theta side of the suction surface");
    if (!(opt.length_unit > 0.0 && opt.coverage_chords > 0.0 && opt.table_step_chords > 0.0))
        throw std::invalid_argument("stream_surface_blade: length_unit, coverage_chords and table_step_chords must be > 0");
    Design::BladeDesign d = design;
    d.layers[sz(k)] = layer;
    const Design::DesignResult res = ev.evaluate(d);
    if (!res.errors.empty()) throw std::invalid_argument("BladeGen cannot place the design: " + res.errors.front());
    const Design::LayerResult& placed = res.layers[sz(k)];
    if (!placed.error.empty()) throw std::invalid_argument("BladeGen cannot place the layer: " + placed.error);
    Design::Placement conf;
    conf.gamma = placed.gamma;
    conf.c = placed.chord_mprime;
    conf.m0 = placed.mprime_le;
    conf.t0 = placed.theta_le;
    conf.sign = 1.0;

    const Blade::FlowPath fp(d.meridional);
    const Streamline line(fp, layer.span);
    const double mp_origin = line.mp_of_u(1.0);  // BladeGen measures m' from the layer's LE point (u = 1)
    const double u_le = line.u_of_mp(mp_origin + conf.m0);
    const double m_le = line.m_of_u(u_le);
    const double r_le = line.xr(u_le).y();
    // (x, y) of the chord frame -> (M, Y) of the (m, r_LE theta) plane [BladeGen length unit], and its u
    auto map = [&](double x, double y, double& u) {
        const Math::Vec2 q = conf(Math::Vec2(x, y));
        u = line.u_of_mp(mp_origin + q.x());
        return Math::Vec2(line.m_of_u(u) - m_le, r_le * (q.y() - conf.t0));
    };
    double u_te = 0.0;
    const Math::Vec2 te = map(1.0, 0.0, u_te);
    const double L = te.norm();
    const double g_my = std::atan2(te.y(), te.x());
    const double cg = std::cos(g_my), sg = std::sin(g_my);
    auto turn_back = [&](const Math::Vec2& p) {
        return Point{(p.x() * cg + p.y() * sg) / L, (-p.x() * sg + p.y() * cg) / L};
    };
    const double m_le_bg = placed.surface.m_of_mprime(conf.m0);

    StreamSurfaceReport rep;
    double M_lo = 0.0, M_hi = te.x();
    // one side: BladeGen's dense sample mapped point by point, its curvature by the chain rule
    auto mapped_side = [&](Profile::Side side, int n_edge, int n_body) {
        const SideSample s = sample_layer_side(layer, ev, side, n_edge, n_body);
        SideSample out;
        const std::size_t n = s.xy.size();
        out.xy.resize(n);
        out.kappa.resize(n);
        for (std::size_t q = 0; q < n; ++q) {
            // the unit tangent of the chord-frame sample (central differences of the dense exact points)
            const std::size_t a = q == 0 ? 0 : q - 1, b = q + 1 == n ? q : q + 1;
            Math::Vec2 t(s.xy[b].m - s.xy[a].m, s.xy[b].y - s.xy[a].y);
            t /= t.norm();
            double u = 0.0;
            const Math::Vec2 P = map(s.xy[q].m, s.xy[q].y, u);
            M_lo = std::min(M_lo, P.x());
            M_hi = std::max(M_hi, P.x());
            out.xy[q] = turn_back(P);
            const double r = line.xr(u).y();
            out.kappa[q] = mapped_curvature({t.x(), t.y()}, s.kappa[q], conf.gamma, conf.c, r, r * line.sin_phi(u), r_le, L);
            // the same point by BladeGen's own conformal map (its StreamSurface table)
            const double bg = placed.surface.m_of_mprime(conf(Math::Vec2(s.xy[q].m, s.xy[q].y)).x()) - m_le_bg;
            rep.map_max_difference = std::max(rep.map_max_difference, std::abs(bg - P.x()) / L);
        }
        return out;
    };
    StreamSurfaceBlade blade;
    BladeSection& b = blade.section;
    auto fit = [&](Profile::Side side, SideFitReport& r) {
        SideSample& kept = side == Profile::Side::Suction ? rep.suction_sample : rep.pressure_sample;
        try {
            SideSample s = mapped_side(side, opt.fit.n_edge, opt.fit.n_body);
            if (opt.keep_samples) kept = s;
            return fit_side(s, opt.fit, &r);
        } catch (const std::runtime_error&) {
            SideSample s = mapped_side(side, 3 * opt.fit.n_edge, 3 * opt.fit.n_body);
            if (opt.keep_samples) kept = s;
            return fit_side(s, opt.fit, &r);
        }
    };
    b.suction = fit(Profile::Side::Suction, rep.fit.suction);
    b.pressure = fit(Profile::Side::Pressure, rep.fit.pressure);
    b.suction_parameters = surface_parameters(opt.fit.n_parameters);
    b.pressure_parameters = b.suction_parameters;
    section_metadata(layer, ev, b, rep.fit);
    // the true blade angles kappa = g + chi_le, g - chi_te (the conformal map keeps angles), measured
    // from the (m, y) chord line, whose angle is g_my
    const double turn = g_my - conf.gamma;
    b.chi_le = b.chi_le - turn;
    b.chi_te = b.chi_te + turn;
    double u_edge = 0.0;
    b.le_center = turn_back(map(b.le_center.m, b.le_center.y, u_edge));
    b.te_center = turn_back(map(b.te_center.m, b.te_center.y, u_edge));
    const double scale = conf.c * r_le / L;  // the edge radii to first order (metadata for reports)
    b.le_radius = b.le_radius * scale;
    b.te_radius = b.te_radius * scale;
    b.designation = b.designation + " on the " + layer.name + " stream surface";

    // the sheet: stations uniform in m, coverage_chords beyond the blade's ends
    const double lo = M_lo - opt.coverage_chords * L, hi = M_hi + opt.coverage_chords * L;
    const int n = static_cast<int>(std::ceil((hi - lo) / (opt.table_step_chords * L))) + 1;
    std::vector<double> ms(sz(n)), rs(sz(n)), bs(sz(n));
    for (int q = 0; q < n; ++q) {
        const double M = lo + (hi - lo) * static_cast<double>(q) / static_cast<double>(n - 1);
        const Streamline::Local l = line.at_m(m_le + M);
        ms[sz(q)] = M / L;
        rs[sz(q)] = l.r / L;
        bs[sz(q)] = l.b;
    }
    rep.b_le = 1.0;
    rep.b_te = line.at_m(m_le + te.x()).b / line.at_m(m_le).b;
    blade.sheet = std::make_shared<const StreamSheet>(std::move(ms), std::move(rs), std::move(bs), r_le / L);
    blade.stagger_deg = Math::rad2deg(g_my);
    blade.n_blades = d.n_blades;
    blade.solidity = 1.0 / ((r_le / L) * (2.0 * Math::kPi / static_cast<double>(d.n_blades)));
    blade.chord = L * opt.length_unit;
    blade.source = "BladeGen layer " + layer.name;
    blade.layer_stagger_deg = Math::rad2deg(conf.gamma);
    blade.conformal_chord = conf.c;
    blade.theta_le_deg = Math::rad2deg(conf.t0);
    blade.r_le = r_le * opt.length_unit;
    blade.r_te = line.xr(u_te).y() * opt.length_unit;
    blade.max_fit_residual = std::max(rep.fit.suction.max_residual, rep.fit.pressure.max_residual);
    if (report) *report = rep;
    return blade;
}

}  // namespace BladeGen::Turbolenza
