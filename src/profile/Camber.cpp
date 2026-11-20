#include "profile/Camber.hpp"

#include "math/BSpline.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BladeGen::Profile {

double polyline_distance(const Vec2& q, const std::vector<Vec2>& poly)
{
    if (poly.empty()) throw std::invalid_argument("polyline_distance: empty polyline");
    double best = (q - poly.front()).squaredNorm();
    for (std::size_t i = 1; i < poly.size(); ++i) {
        const Vec2 a = poly[i - 1], d = poly[i] - a;
        const double l2 = d.squaredNorm();
        double t = l2 > 0.0 ? (q - a).dot(d) / l2 : 0.0;
        t = std::clamp(t, 0.0, 1.0);
        best = std::min(best, (q - (a + t * d)).squaredNorm());
    }
    return std::sqrt(best);
}

double polygon_area(const std::vector<Vec2>& loop)
{
    double a = 0.0;
    for (std::size_t i = 0; i < loop.size(); ++i) a += Math::cross2(loop[i], loop[(i + 1) % loop.size()]);
    return 0.5 * a;
}

CamberLine camber_line(const std::vector<Vec2>& ss, const std::vector<Vec2>& ps, int n)
{
    if (ss.size() < 4 || ps.size() < 4) throw std::invalid_argument("camber_line: each side needs at least 4 points");
    if (n < 5) throw std::invalid_argument("camber_line: n >= 5");
    for (const auto* side : {&ss, &ps})
        for (const auto& p : *side)
            if (!p.allFinite()) throw std::invalid_argument("camber_line: non-finite point");
    // SS arc length
    std::vector<double> s(ss.size(), 0.0);
    for (std::size_t i = 1; i < ss.size(); ++i) s[i] = s[i - 1] + (ss[i] - ss[i - 1]).norm();
    const double L = s.back();
    if (!(L > 0.0)) throw std::invalid_argument("camber_line: zero-length suction side");

    auto at = [&](double sq, Vec2& p, Vec2& t) {
        const auto it = std::upper_bound(s.begin(), s.end(), sq);
        const auto k = std::clamp<std::size_t>(static_cast<std::size_t>(it - s.begin()), 1, s.size() - 1);
        const double w = (sq - s[k - 1]) / std::max(s[k] - s[k - 1], 1e-300);
        p = ss[k - 1] + w * (ss[k] - ss[k - 1]);
        // vertex tangents by central differences, interpolated along the segment
        auto vt = [&](std::size_t i) {
            const std::size_t a = i == 0 ? 0 : i - 1, b = std::min(i + 1, ss.size() - 1);
            return Math::unit<2>(ss[b] - ss[a]);
        };
        t = Math::unit<2>((1.0 - w) * vt(k - 1) + w * vt(k));
    };

    const double gap_tol = 1e-9 * L;
    const bool open_le = (ss.front() - ps.front()).norm() > gap_tol, open_te = (ss.back() - ps.back()).norm() > gap_tol;
    CamberLine c;
    c.point.push_back(ss.front());
    c.thickness.push_back(0.0);
    for (int k = 0; k < n; ++k) {
        // stations clustered toward both ends, strictly inside (0, 1)
        const double sig = 0.5 * (1.0 - std::cos(Math::kPi * (k + 1.0) / (n + 1.0)));
        Vec2 p, t;
        at(sig * L, p, t);
        const Vec2 nrm = Math::right_normal(t);  // SS runs LE -> TE with the profile below it
        auto f = [&](double rho) { return polyline_distance(p + rho * nrm, ps) - rho; };
        const double d0 = polyline_distance(p, ps);
        if (!(d0 > 1e-12)) continue;  // touching sides: no inscribed circle
        double lo = 0.0, hi = 0.5 * d0;
        int grow = 0;
        while (f(hi) > 0.0 && grow++ < 40) {
            lo = hi;
            hi *= 2.0;
        }
        if (grow > 40) continue;
        for (int it = 0; it < 60 && hi - lo > 1e-13 * std::max(1.0, hi); ++it) {
            const double mid = 0.5 * (lo + hi);
            (f(mid) > 0.0 ? lo : hi) = mid;
        }
        const double rho = 0.5 * (lo + hi);
        const Vec2 ctr = p + rho * nrm;
        if (ctr.x() <= c.point.back().x()) continue;  // keep the locus monotone in x
        // a circle that touches the pressure side only at one of its end points (open edge)
        // is not an inscribed circle of the section
        // an inscribed circle stays inside the section: it must not cut the suction side
        // (it does on a nose whose radius of curvature is below rho)
        if (polyline_distance(ctr, ss) < rho * (1.0 - 2e-3)) continue;  // (tolerance: polyline sag)
        // the nose region (circles nearly through the LE / TE point) is degenerate for the
        // medial axis: those centres are left out and the spline continues to the edge point
        if ((ctr - 0.5 * (ss.front() + ps.front())).norm() < 1.5 * rho) continue;
        if ((ctr - 0.5 * (ss.back() + ps.back())).norm() < 1.5 * rho) continue;
        if (open_le && (ctr - ps.front()).norm() <= rho * (1.0 + 1e-9)) continue;
        if (open_te && (ctr - ps.back()).norm() <= rho * (1.0 + 1e-9)) continue;
        c.point.push_back(ctr);
        c.thickness.push_back(2.0 * rho);
    }
    if (c.point.size() < 4) throw std::invalid_argument("camber_line: no inscribed circles found (sides crossed?)");
    // The locus is smoothed by a least-squares cubic B-spline y(x) pinned at the LE and TE
    // points (knots clustered toward both ends), so the camber angle is the derivative of
    // a C2 curve instead of a difference quotient of the bisection results.
    const Vec2 pa = 0.5 * (ss.front() + ps.front()), pb = 0.5 * (ss.back() + ps.back());  // open edges: mid points
    const double x0 = pa.x(), x1 = pb.x();
    if (!(x1 > x0)) throw std::invalid_argument("camber_line: the TE point must lie at larger x than the LE point");
    std::vector<double> rx, ry, rt;
    for (std::size_t i = 1; i < c.point.size(); ++i)
        if (c.point[i].x() > x0 && c.point[i].x() < x1) {
            rx.push_back((c.point[i].x() - x0) / (x1 - x0));
            ry.push_back(c.point[i].y());
            rt.push_back(c.thickness[i]);
        }
    constexpr int p = 3;
    const int nc = std::clamp(static_cast<int>(rx.size()) / 4, 6, 18);
    std::vector<double> U(static_cast<std::size_t>(p + 1), 0.0);
    for (int k = 1; k < nc - p; ++k) U.push_back(0.5 * (1.0 - std::cos(Math::kPi * k / (nc - p))));
    for (int k = 0; k <= p; ++k) U.push_back(1.0);
    auto basis = [&](double t, int nd, int& span) {
        span = Math::find_span(nc - 1, p, t, U);
        return Math::ders_basis_funs(span, t, p, nd, U);
    };
    const int nf = nc - 2;  // free control values (the first and the last are the LE / TE y)
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(rx.size()), nf);
    Eigen::VectorXd rhs(static_cast<Eigen::Index>(rx.size()));
    for (std::size_t i = 0; i < rx.size(); ++i) {
        int span = 0;
        const auto N = basis(rx[i], 0, span);
        double r = ry[i];
        for (int j = 0; j <= p; ++j) {
            const int idx = span - p + j;
            const double w = N[0][static_cast<std::size_t>(j)];
            if (idx == 0) r -= w * pa.y();
            else if (idx == nc - 1) r -= w * pb.y();
            else A(static_cast<Eigen::Index>(i), idx - 1) += w;
        }
        rhs(static_cast<Eigen::Index>(i)) = r;
    }
    // a light second-difference penalty keeps the system well posed where data are sparse
    Eigen::MatrixXd Q = A.transpose() * A;
    Eigen::VectorXd g = A.transpose() * rhs;
    const double lam = 1e-9 * Q.trace() / nf;
    for (int j = 0; j < nf; ++j) Q(j, j) += lam;
    const Eigen::VectorXd cy = Q.ldlt().solve(g);
    std::vector<double> P(static_cast<std::size_t>(nc));
    P.front() = pa.y();
    P.back() = pb.y();
    for (int j = 0; j < nf; ++j) P[static_cast<std::size_t>(j + 1)] = cy(j);

    CamberLine out;
    const int m = std::max(n, 21);
    for (int k = 0; k < m; ++k) {
        const double t = 0.5 * (1.0 - std::cos(Math::kPi * k / (m - 1)));
        int span = 0;
        const auto N = basis(t, 1, span);
        double y = 0.0, dy = 0.0;
        for (int j = 0; j <= p; ++j) {
            const double v = P[static_cast<std::size_t>(span - p + j)];
            y += N[0][static_cast<std::size_t>(j)] * v;
            dy += N[1][static_cast<std::size_t>(j)] * v;
        }
        out.point.emplace_back(x0 + t * (x1 - x0), y);
        out.angle_deg.push_back(Math::rad2deg(std::atan(dy / (x1 - x0))));
        // thickness: the measured inscribed diameters, 0 at the LE / TE points
        double th = 0.0;
        if (t > rx.front() && t < rx.back()) {
            const auto it = std::upper_bound(rx.begin(), rx.end(), t);
            const auto q = static_cast<std::size_t>(it - rx.begin());
            const double w = (t - rx[q - 1]) / (rx[q] - rx[q - 1]);
            th = rt[q - 1] + w * (rt[q] - rt[q - 1]);
        } else if (t <= rx.front()) {
            th = rt.front() * t / rx.front();
        } else {
            th = rt.back() * (1.0 - t) / (1.0 - rx.back());
        }
        out.thickness.push_back(k == 0 || k == m - 1 ? 0.0 : th);
    }
    c = std::move(out);
    const auto it = std::max_element(c.thickness.begin(), c.thickness.end());
    c.max_thickness = *it;
    c.max_thickness_x = c.point[static_cast<std::size_t>(it - c.thickness.begin())].x();
    return c;
}

}  // namespace BladeGen::Profile
