// The BladeGen adapter: a design layer's profile -> TURBOlenza's BladeSection (see BladeGenSource.hpp).
#include "turbolenza/BladeGenSource.hpp"

#include "math/Vec.hpp"
#include "naca/Naca65.hpp"
#include "profile/Profile2D.hpp"

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
    if (report) *report = rep;
    return b;
}

}  // namespace BladeGen::Turbolenza
