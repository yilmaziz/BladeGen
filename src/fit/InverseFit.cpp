#include "fit/InverseFit.hpp"

#include "math/LeastSquares.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BladeGen::Fit {

using Profile::Piece;
using Profile::Profile2D;
using Profile::ProfileParams;

std::vector<FitParam> default_fit_params(const ProfileParams& s)
{
    // blade2d CDA/fit_cda_c2.py SPEC (bounds), starting at `s` (clamped into the box)
    std::vector<FitParam> v = {
        {"le.radius", 0, 0.002, 0.040},       {"le.stretch_ab", 0, 1.0, 6.0},
        {"le.asym", 0, -0.010, 0.010},        {"le.beta_deg", 0, -5.0, 35.0},
        {"le.d_beta_ss_deg", 0, 2.0, 35.0},   {"le.d_beta_ps_deg", 0, 2.0, 35.0},
        {"te.radius", 0, 0.001, 0.020},       {"te.stretch_ab", 0, 0.2, 4.0},
        {"te.asym", 0, -0.005, 0.005},        {"te.beta_deg", 0, -10.0, 25.0},
        {"te.d_beta_ss_deg", 0, 0.0, 25.0},   {"te.d_beta_ps_deg", 0, 0.0, 25.0},
        {"ss.handle_le_frac", 0, 0.02, 0.50}, {"ss.ctrl2_mprime", 0, 0.15, 0.70},
        {"ss.ctrl2_theta", 0, 0.0, 0.35},     {"ss.handle_te_frac", 0, 0.02, 0.60},
        {"thickness.t_max", 0, 0.010, 0.150}, {"thickness.mprime_max", 0, 0.10, 0.70},
        {"thickness.t_le", 0, 0.002, 0.080},  {"thickness.t_te", 0, 0.001, 0.030},
        {"thickness.shape_le", 0, 0.3, 3.0},  {"thickness.shape_te", 0, 0.3, 4.0},
        {"ps.handle_le_frac", 0, 0.02, 0.50}, {"ps.handle_te_frac", 0, 0.02, 0.50},
    };
    for (auto& p : v) {
        const double x = Profile::get_param(s, p.name);
        p.x0 = std::clamp(x, p.lo, p.hi);
    }
    return v;
}

namespace {

double seg_dist(const Vec2& p, const Vec2& a, const Vec2& b)
{
    const Vec2 ab = b - a;
    const double L2 = ab.squaredNorm();
    const double t = L2 > 1e-300 ? std::clamp((p - a).dot(ab) / L2, 0.0, 1.0) : 0.0;
    return (p - (a + t * ab)).norm();
}

std::vector<double> poly_dist(const std::vector<Vec2>& pts, const std::vector<Vec2>& poly)
{
    std::vector<double> d(pts.size(), std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t k = 0; k + 1 < poly.size(); ++k) d[i] = std::min(d[i], seg_dist(pts[i], poly[k], poly[k + 1]));
    return d;
}

std::vector<Vec2> sample(const Math::BSpline2& c, int n, bool reverse = false)
{
    std::vector<Vec2> out;
    for (int i = 0; i < n; ++i) {
        const double u = static_cast<double>(i) / (n - 1);
        out.push_back(c.eval(reverse ? 1.0 - u : u));
    }
    return out;
}

/// Concave dips of a connector's curvature below its junction values, normalised.
void fairness_residuals(const Profile2D& p, int k, double w, std::vector<double>& out)
{
    const double ta = p.breaks()[static_cast<std::size_t>(k)], tb = p.breaks()[static_cast<std::size_t>(k + 1)];
    const int n = 40;
    std::vector<double> kap(n);
    for (int i = 0; i < n; ++i) kap[static_cast<std::size_t>(i)] = p.curvature(ta + (tb - ta) * i / (n - 1), i == n - 1);
    const double kend = std::min(kap.front(), kap.back());
    const double scale = std::max(1.0, *std::max_element(kap.begin(), kap.end()));
    for (double v : kap) out.push_back(w * std::max(0.0, std::min(0.0, kend) - v) / scale);
}

}  // namespace

void side_distances(const Profile2D& p, const std::vector<Vec2>& ss, const std::vector<Vec2>& ps,
                    std::vector<double>& d_ss, std::vector<double>& d_ps, int n_edge, int n_body)
{
    // upper: LE connector (P_LE_PS -> P_LE_SS) + SS + TE connector (P_TE_SS -> P_TE_PS)
    // lower: LE connector reversed + PS + TE connector reversed
    const auto le = sample(p.le_connector(), n_edge), te = sample(p.te_connector(), n_edge);
    const auto sb = sample(p.suction_side(), n_body), pb = sample(p.pressure_side(), n_body);
    std::vector<Vec2> upper(le);
    upper.insert(upper.end(), sb.begin() + 1, sb.end());
    upper.insert(upper.end(), te.begin() + 1, te.end());
    std::vector<Vec2> lower(le.rbegin(), le.rend());
    lower.insert(lower.end(), pb.begin() + 1, pb.end());
    lower.insert(lower.end(), te.rbegin() + 1, te.rend());
    d_ss = poly_dist(ss, upper);
    d_ps = poly_dist(ps, lower);
}

FitResult fit_profile(const std::vector<Vec2>& ss, const std::vector<Vec2>& ps, const ProfileParams& start,
                      const FitOptions& opt)
{
    if (ss.size() < 3 || ps.size() < 3) throw std::invalid_argument("fit_profile: need at least 3 points per side");
    const auto spec = opt.params.empty() ? default_fit_params(start) : opt.params;
    const auto n = static_cast<Eigen::Index>(spec.size());
    Eigen::VectorXd x0(n), lb(n), ub(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        const auto& s = spec[static_cast<std::size_t>(i)];
        x0[i] = s.x0;
        lb[i] = s.lo;
        ub[i] = s.hi;
    }
    auto params_of = [&](const Eigen::VectorXd& x) {
        ProfileParams p = start;
        for (Eigen::Index i = 0; i < n; ++i) Profile::set_param(p, spec[static_cast<std::size_t>(i)].name, x[i]);
        return p;
    };
    const auto m = static_cast<Eigen::Index>(ss.size() + ps.size());
    const Eigen::Index n_fair = opt.fairness_weight > 0.0 ? 80 : 0;
    std::vector<double> dss, dps, fair;
    auto residuals = [&](const Eigen::VectorXd& x) {
        Eigen::VectorXd r(m + n_fair);
        try {
            const Profile2D p(params_of(x));
            side_distances(p, ss, ps, dss, dps, opt.n_edge, opt.n_body);
            Eigen::Index k = 0;
            for (double d : dss) r[k++] = d;
            for (double d : dps) r[k++] = d;
            if (n_fair) {
                fair.clear();
                fairness_residuals(p, 0, opt.fairness_weight, fair);
                fairness_residuals(p, 2, opt.fairness_weight, fair);
                for (double f : fair) r[k++] = f;
            }
            if (!r.allFinite()) throw std::runtime_error("non-finite");
        } catch (const std::exception&) {
            r.setConstant(0.5);  // an invalid parameter set: a large, flat residual
        }
        return r;
    };
    const auto t0 = std::chrono::steady_clock::now();
    FitResult res;
    Math::LMOptions lo;
    lo.max_iterations = opt.max_iterations;
    lo.fd_step = 1e-6;
    Eigen::VectorXd x = x0;
    for (int s = 0; s < std::max(1, opt.stages); ++s) {
        const auto r = Math::levenberg_marquardt(residuals, x, lb, ub, lo);
        x = r.x;
        res.evaluations += r.evaluations;
        res.status = r.status;
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    res.params = params_of(x);
    const Profile2D p(res.params);
    side_distances(p, ss, ps, res.d_ss, res.d_ps, 4 * opt.n_edge, 4 * opt.n_body);
    double s2 = 0.0, sss = 0.0, sps = 0.0;
    for (double d : res.d_ss) { s2 += d * d; sss += d * d; res.max = std::max(res.max, d); }
    for (double d : res.d_ps) { s2 += d * d; sps += d * d; res.max = std::max(res.max, d); }
    res.rms = std::sqrt(s2 / static_cast<double>(res.d_ss.size() + res.d_ps.size()));
    res.rms_ss = std::sqrt(sss / static_cast<double>(res.d_ss.size()));
    res.rms_ps = std::sqrt(sps / static_cast<double>(res.d_ps.size()));
    for (Eigen::Index i = 0; i < n; ++i)
        if (std::min(x[i] - lb[i], ub[i] - x[i]) < 1e-6 * (ub[i] - lb[i]))
            res.at_bound.push_back(spec[static_cast<std::size_t>(i)].name);
    return res;
}

Vec2 Normalized::to_input(const Vec2& q) const
{
    const Vec2 e(q.x(), mirrored ? -q.y() : q.y());
    const double c = std::cos(stagger), s = std::sin(stagger);
    return le + chord * Vec2(c * e.x() - s * e.y(), s * e.x() + c * e.y());
}

Normalized normalize_contour(const std::vector<Vec2>& in)
{
    std::vector<Vec2> pts(in);
    if (pts.size() > 1 && (pts.front() - pts.back()).norm() < 1e-14) pts.pop_back();
    if (pts.size() < 8) throw std::invalid_argument("normalize_contour: need at least 8 points");
    // longest chord: alternate farthest-point searches
    std::size_t ile = 0, ite = 0;
    for (std::size_t i = 1; i < pts.size(); ++i) {
        if (pts[i].x() < pts[ile].x()) ile = i;
        if (pts[i].x() > pts[ite].x()) ite = i;
    }
    for (int it = 0; it < 8; ++it) {
        std::size_t a = ile, b = ite;
        for (std::size_t i = 0; i < pts.size(); ++i)
            if ((pts[i] - pts[ite]).norm() > (pts[a] - pts[ite]).norm()) a = i;
        for (std::size_t i = 0; i < pts.size(); ++i)
            if ((pts[i] - pts[a]).norm() > (pts[b] - pts[a]).norm()) b = i;
        if (a == ile && b == ite) break;
        ile = a;
        ite = b;
    }
    Normalized out;
    out.le = pts[ile];
    out.te = pts[ite];
    const Vec2 d = out.te - out.le;
    out.chord = d.norm();
    out.stagger = std::atan2(d.y(), d.x());
    const double c = std::cos(out.stagger), s = std::sin(out.stagger);
    auto to_local = [&](const Vec2& p) {
        const Vec2 q = (p - out.le) / out.chord;
        return Vec2(c * q.x() + s * q.y(), -s * q.x() + c * q.y());
    };
    // the two sides from LE to TE
    const std::size_t n = pts.size();
    std::vector<Vec2> a, b;
    for (std::size_t k = 0;; ++k) {
        const std::size_t i = (ile + k) % n;
        a.push_back(to_local(pts[i]));
        if (i == ite) break;
    }
    for (std::size_t k = 0;; ++k) {
        const std::size_t i = (ile + n - k) % n;
        b.push_back(to_local(pts[i]));
        if (i == ite) break;
    }
    auto mean_eta = [](const std::vector<Vec2>& v) {
        double sum = 0.0;
        for (const auto& q : v) sum += q.y();
        return sum / static_cast<double>(v.size());
    };
    const double camber = 0.5 * (mean_eta(a) + mean_eta(b));
    out.mirrored = camber < 0.0;
    if (out.mirrored) {
        for (auto& q : a) q.y() = -q.y();
        for (auto& q : b) q.y() = -q.y();
    }
    if (mean_eta(a) >= mean_eta(b)) { out.ss = a; out.ps = b; }
    else { out.ss = b; out.ps = a; }
    return out;
}

}  // namespace BladeGen::Fit
