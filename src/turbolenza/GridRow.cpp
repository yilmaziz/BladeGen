// The reference blade of a C2 cascade as grid.py reads it: port of v1 grid._turning_points,
// _Cut, _Row.narrowest, _CurveRow and _check_walls (grid.py sections 3, 4, 7; Aungier 2003
// Sec. 5.3, Fig. 5-5).  See GridDetail.hpp.
#include "turbolenza/GridDetail.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace BladeGen::Turbolenza::detail {

namespace {

std::size_t sz(int k) { return static_cast<std::size_t>(k); }

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

double horner(const std::vector<double>& c, double s)
{
    double r = 0.0;
    for (std::size_t k = c.size(); k-- > 0;) r = r * s + c[k];
    return r;
}

/// One root of a polynomial whose values at the ends of [x0, x1] differ in sign (f0 != 0).
double bisect_root(const std::vector<double>& c, double x0, double x1, double f0)
{
    double lo = x0, hi = x1, flo = f0;
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (!(mid > lo && mid < hi)) break;
        const double fm = horner(c, mid);
        if (fm == 0.0) return mid;
        if ((fm < 0.0) == (flo < 0.0)) {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

std::vector<double> roots_on(std::vector<double> c, double a, double b)
{
    int d = static_cast<int>(c.size()) - 1;
    while (d >= 0 && c[sz(d)] == 0.0) --d;
    if (d <= 0) return {};
    c.resize(sz(d + 1));
    if (d == 1) {
        const double r = -c[0] / c[1];
        return (a <= r && r <= b) ? std::vector<double>{r} : std::vector<double>{};
    }
    std::vector<double> dc(sz(d));
    for (int k = 0; k < d; ++k) dc[sz(k)] = static_cast<double>(k + 1) * c[sz(k + 1)];
    std::vector<double> pts{a};
    for (double r : roots_on(dc, a, b))
        if (r > pts.back() && r < b) pts.push_back(r);
    pts.push_back(b);
    std::vector<double> out;
    for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
        const double f0 = horner(c, pts[i]), f1 = horner(c, pts[i + 1]);
        if (f0 == 0.0) {
            out.push_back(pts[i]);
            continue;
        }
        if (f1 != 0.0 && ((f0 < 0.0) != (f1 < 0.0))) out.push_back(bisect_root(c, pts[i], pts[i + 1], f0));
    }
    if (horner(c, b) == 0.0) out.push_back(b);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

double component(const Point& p, int axis) { return axis == 0 ? p.m : p.y; }

}  // namespace

// ── NumPy / Python semantics ─────────────────────────────────────────────────

double np_mod(double a, double b)
{
    // numpy npy_divmod: fmod, moved into the divisor's sign
    double mod = std::fmod(a, b);
    if (mod != 0.0) {
        if ((b < 0.0) != (mod < 0.0)) mod += b;
    } else {
        mod = std::copysign(0.0, b);
    }
    return mod;
}

long long py_round(double x) { return static_cast<long long>(std::nearbyint(x)); }

std::vector<double> linspace(double a, double b, int n)
{
    std::vector<double> out(sz(std::max(n, 0)));
    if (n == 1) out[0] = a;
    if (n < 2) return out;
    const double step = (b - a) / static_cast<double>(n - 1);
    for (int k = 0; k < n; ++k) out[sz(k)] = static_cast<double>(k) * step + a;
    out.back() = b;
    return out;
}

std::vector<std::size_t> stable_argsort(const std::vector<double>& v)
{
    std::vector<std::size_t> order(v.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    // NumPy sorts NaN last
    std::stable_sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) {
        const double a = v[i], b = v[j];
        if (std::isnan(a)) return false;
        if (std::isnan(b)) return true;
        return a < b;
    });
    return order;
}

std::vector<double> unique(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

std::string fmt_g(double x, int precision)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", precision, x);
    return buf;
}

std::string fmt_e(double x, int precision)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*e", precision, x);
    return buf;
}

// ── polynomials and turning points ───────────────────────────────────────────

std::vector<double> polynomial_roots(const std::vector<double>& c, double h)
{
    if (!(h >= 0.0)) throw std::invalid_argument("polynomial_roots: the interval [0, h] needs h >= 0");
    return roots_on(c, 0.0, h);
}

std::vector<double> turning_points(const ParametricBladeCurve& curve, int axis)
{
    // v1 _turning_points: PPoly.from_spline(BSpline(...).derivative()).roots(extrapolate=False)
    const int p = curve.degree();
    const auto& t = curve.knots();
    const int n = static_cast<int>(curve.control_points().size());
    std::vector<double> roots;
    double previous_end = 0.0;  // the previous span polynomial at its right end
    bool have_previous = false;
    for (int l = p; l < n; ++l) {
        const double a = t[sz(l)], b = t[sz(l + 1)];
        if (!(b > a)) continue;
        const auto d = curve.derivatives(a, p);
        std::vector<double> c(sz(p));
        double factorial = 1.0;
        for (int k = 0; k < p; ++k) {
            if (k > 0) factorial *= static_cast<double>(k);
            c[sz(k)] = component(d[sz(k + 1)], axis) / factorial;
        }
        // scipy's discontinuity root: a sign change across the join of two pieces
        if (have_previous && ((previous_end < 0.0 && c[0] > 0.0) || (previous_end > 0.0 && c[0] < 0.0)))
            roots.push_back(a);
        for (double s : polynomial_roots(c, b - a)) {
            const double r = a + s;
            if (a <= r && r <= b) roots.push_back(r);
        }
        previous_end = horner(c, b - a);
        have_previous = true;
    }
    std::vector<double> kept;
    for (double r : roots)
        if (std::isfinite(r) && r > 0.0 && r < 1.0) kept.push_back(r);
    kept = unique(std::move(kept));
    if (kept.size() > 1) {
        std::vector<double> merged{kept[0]};
        for (std::size_t k = 1; k < kept.size(); ++k)
            if (kept[k] - kept[k - 1] > kTurnMerge) merged.push_back(kept[k]);
        kept = std::move(merged);
    }
    if (kept.empty()) return kept;
    std::vector<double> lo, hi, start;
    for (std::size_t k = 0; k < kept.size(); ++k) {
        const double left = kept[k] - (k == 0 ? 0.0 : kept[k - 1]);
        const double right = (k + 1 == kept.size() ? 1.0 : kept[k + 1]) - kept[k];
        const double half = std::min(kTurnBracket, 0.5 * std::min(left, right));
        const double a = kept[k] - half, b = kept[k] + half;
        const double sa = component(curve.derivative(a, 1), axis);
        if (np_sign(sa) * np_sign(component(curve.derivative(b, 1), axis)) < 0.0) {
            lo.push_back(a);
            hi.push_back(b);
            start.push_back(sa);
        }
    }
    std::vector<double> out(lo.size());
    for (std::size_t k = 0; k < lo.size(); ++k) {
        double a = lo[k], b = hi[k], s = start[k];
        for (int it = 0; it < kBisections; ++it) {
            const double mid = 0.5 * (a + b);
            const double value = component(curve.derivative(mid, 1), axis);
            if (np_sign(value) == np_sign(s)) {
                a = mid;
                s = value;
            } else {
                b = mid;
            }
        }
        out[k] = 0.5 * (a + b);
    }
    return out;
}

// ── the cut ──────────────────────────────────────────────────────────────────

int Cut::top_segment(std::size_t i) const
{
    const auto& s = segment[i];
    return count[i] >= 1 ? s[sz(count[i] - 1)] : s.back();
}

// ── CurveRow ─────────────────────────────────────────────────────────────────

CurveRow::CurveRow(const PassageBoundaries& b)
{
    owner = b.blade_curves;
    if (!owner) throw GridError("the passage boundaries carry no blade curves: v2 grids are built from the C2 curves");
    curves = {&owner->suction, &owner->pressure};
    touch = kTouch * b.chord;
    lengths = {curves[0]->reference_length(), curves[1]->reference_length()};
    total = lengths[0] + lengths[1];
    for (int index = 0; index < 2; ++index) {
        std::vector<double> splits{0.0};
        for (int axis = 0; axis < 2; ++axis)
            for (double u : turning_points(*curves[sz(index)], axis)) splits.push_back(u);
        splits.push_back(1.0);
        splits = unique(std::move(splits));
        for (std::size_t k = 0; k + 1 < splits.size(); ++k) {
            curve_of.push_back(index);
            u0.push_back(splits[k]);
            u1.push_back(splits[k + 1]);
        }
    }
    const int nb = n_branches();
    std::vector<double> ends_m, ends_tau;
    for (int k = 0; k < nb; ++k) {
        m0.push_back(point(k, u0[sz(k)]).m);
        m1.push_back(point(k, u1[sz(k)]).m);
        lo_m.push_back(std::min(m0.back(), m1.back()));
        hi_m.push_back(std::max(m0.back(), m1.back()));
    }
    for (int k = 0; k < nb; ++k) {
        ends_m.push_back(m0[sz(k)]);
        ends_tau.push_back(tau(k, u0[sz(k)]));
    }
    for (int k = 0; k < nb; ++k) {
        ends_m.push_back(m1[sz(k)]);
        ends_tau.push_back(tau(k, u1[sz(k)]));
    }
    const auto order = stable_argsort(ends_m);
    tau_tip = ends_tau[order.front()];
    tau_tail = ends_tau[order.back()];
    m_tip = ends_m[order.front()];
    m_tail = ends_m[order.back()];
    // the branch end behind each extreme (ends_m: the u0 ends of the branches, then their u1 ends)
    auto branch_end = [&](std::size_t e, int& branch, double& u) {
        const std::size_t n = sz(nb);
        branch = static_cast<int>(e < n ? e : e - n);
        u = e < n ? u0[e] : u1[e - n];
    };
    branch_end(order.front(), tip_branch, tip_u);
    branch_end(order.back(), tail_branch, tail_u);
    tip_point = point(tip_branch, tip_u);
    tail_point = point(tail_branch, tail_u);
    const double span = np_mod(tau_tail - tau_tip, total);
    for (int k = 0; k < nb; ++k) {
        const double middle = tau(k, 0.5 * (u0[sz(k)] + u1[sz(k)]));
        lower_contour.push_back(np_mod(middle - tau_tip, total) < span);
    }
}

double CurveRow::tau(int branch, double u) const
{
    return curve_of[sz(branch)] == 0 ? (1.0 - u) * lengths[0] : lengths[0] + u * lengths[1];
}

Point CurveRow::point(int branch, double u) const { return curves[sz(curve_of[sz(branch)])]->evaluate(u); }

double CurveRow::solve(int branch, double m) const
{
    double lo = u0[sz(branch)], hi = u1[sz(branch)];
    double start = m0[sz(branch)] - m;  // m(u) - m at the branch start; its sign stays on the lo side
    for (int it = 0; it < kBisections; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double value = point(branch, mid).m - m;
        if (np_sign(value) == np_sign(start)) {
            lo = mid;
            start = value;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

Cut CurveRow::cut(const std::vector<double>& m) const
{
    const int nb = n_branches();
    Cut out;
    out.m = m;
    out.touch = touch;
    bool odd = false;
    for (double station : m) {
        std::vector<double> height(sz(nb), kInf), u(sz(nb), kNaN);
        int count = 0;
        for (int k = 0; k < nb; ++k) {
            if (lo_m[sz(k)] <= station && station < hi_m[sz(k)]) {
                u[sz(k)] = solve(k, station);
                height[sz(k)] = point(k, u[sz(k)]).y;
                ++count;
            }
        }
        const auto order = stable_argsort(height);
        std::vector<int> seg(sz(nb));
        std::vector<double> h(sz(nb)), uu(sz(nb));
        for (int k = 0; k < nb; ++k) {
            seg[sz(k)] = static_cast<int>(order[sz(k)]);
            h[sz(k)] = height[order[sz(k)]];
            uu[sz(k)] = u[order[sz(k)]];
        }
        if (count % 2) odd = true;
        // metal intervals, bottom to top; numpy's sum: the first term plus the rest in order
        std::vector<double> pairs;
        for (int j = 0; 2 * j + 1 < nb; ++j) {
            const double top = std::isfinite(h[sz(2 * j + 1)]) ? h[sz(2 * j + 1)] : 0.0;
            const double bottom = std::isfinite(h[sz(2 * j)]) ? h[sz(2 * j)] : 0.0;
            pairs.push_back(j < count / 2 ? top - bottom : 0.0);
        }
        double metal = 0.0;
        if (!pairs.empty()) {
            double rest = 0.0;
            for (std::size_t j = 1; j < pairs.size(); ++j) rest += pairs[j];
            metal = pairs.size() > 1 ? pairs[0] + rest : pairs[0];
        }
        out.count.push_back(count);
        out.segment.push_back(std::move(seg));
        out.height.push_back(std::move(h));
        out.u.push_back(std::move(uu));
        out.metal.push_back(metal);
    }
    if (odd) throw GridError("a station line crosses the blade's curves an odd number of times: the outline is not closed");
    return out;
}

namespace {

void check_walls(const WallSide& lower, const WallSide& upper, const std::array<std::vector<double>, 2>& splits,
                 double tolerance, int first_blade_i)
{
    const char* names[2] = {"lower", "upper"};
    const WallSide* walls[2] = {&lower, &upper};
    for (int w = 0; w < 2; ++w) {
        const WallSide& wall = *walls[w];
        for (std::size_t k = 0; k + 1 < wall.contour.size(); ++k)
            if (!(wall.contour[k + 1] - wall.contour[k] >= -tolerance))
                throw GridError(std::string("the ") + names[w] + " wall runs backward along its contour between blade stations " +
                                std::to_string(first_blade_i + static_cast<int>(k)) + " and " +
                                std::to_string(first_blade_i + static_cast<int>(k) + 1));
        for (std::size_t k = 0; k + 1 < wall.surface.size(); ++k) {
            if (wall.surface[k + 1] == wall.surface[k]) continue;
            bool passed = false;
            for (double split : splits[sz(w)])
                if (wall.contour[k] - tolerance <= split && split <= wall.contour[k + 1] + tolerance) passed = true;
            if (!passed) {
                auto name = [](Surface s) { return s == Surface::Suction ? "suction" : "pressure"; };
                throw GridError(std::string("the ") + names[w] + " wall changes from the " + name(wall.surface[k]) +
                                " to the " + name(wall.surface[k + 1]) + " surface between blade stations " +
                                std::to_string(first_blade_i + static_cast<int>(k)) + " and " +
                                std::to_string(first_blade_i + static_cast<int>(k) + 1) +
                                " without passing geometry.py's LE or TE point");
            }
        }
    }
}

}  // namespace

WallIdentity CurveRow::identity(const Cut& cut, int first, int last, bool touching) const
{
    WallIdentity out;
    std::vector<int> top, bottom;
    std::vector<double> u_top, u_bottom;
    // touching: the stations between the two touching stations are read from the cut (as v1)
    const int cut_first = touching ? first + 1 : first, cut_last = touching ? last - 1 : last;
    if (touching) {
        // the tip entry first (its station `first`); the cut's rows follow
        out.stations.push_back(first);
        for (WallSide* w : {&out.lower, &out.upper}) {
            w->surface.push_back(surface_of(tip_branch));
            w->parameter.push_back(tip_u);
            w->u.push_back(tip_u);
            w->hidden.push_back(0);
        }
        out.lower.contour.push_back(np_mod(tau_tip - tau(tip_branch, tip_u), total));
        out.upper.contour.push_back(np_mod(tau(tip_branch, tip_u) - tau_tip, total));
    }
    for (int i = cut_first; i <= cut_last; ++i) {
        const std::size_t r = sz(i);
        const int c = cut.count[r];
        top.push_back(cut.top_segment(r));
        bottom.push_back(cut.segment[r][0]);
        u_top.push_back(c >= 1 ? cut.u[r][sz(c - 1)] : cut.u[r].back());
        u_bottom.push_back(cut.u[r][0]);
    }
    for (std::size_t k = 0; k < top.size(); ++k)
        if (lower_contour[sz(top[k])] || !lower_contour[sz(bottom[k])])
            throw GridError("a wall node is not on its contour: the top of a cut must lie on the blade's upper contour "
                            "and the bottom on its lower contour");
    for (std::size_t k = 0; k < top.size(); ++k) {
        const std::size_t r = sz(cut_first) + k;
        int on_lower = 0, off_lower = 0;
        for (int j = 0; j < cut.count[r]; ++j) (lower_contour[sz(cut.segment[r][sz(j)])] ? on_lower : off_lower) += 1;
        out.stations.push_back(cut_first + static_cast<int>(k));
        out.lower.surface.push_back(surface_of(top[k]));
        out.lower.parameter.push_back(u_top[k]);
        out.lower.u.push_back(u_top[k]);
        out.lower.contour.push_back(np_mod(tau_tip - tau(top[k], u_top[k]), total));
        out.lower.hidden.push_back(off_lower - 1);
        out.upper.surface.push_back(surface_of(bottom[k]));
        out.upper.parameter.push_back(u_bottom[k]);
        out.upper.u.push_back(u_bottom[k]);
        out.upper.contour.push_back(np_mod(tau(bottom[k], u_bottom[k]) - tau_tip, total));
        out.upper.hidden.push_back(on_lower - 1);
    }
    if (touching) {
        // the tail entry last (its station `last`): the lower wall reaches the tail end along the
        // upper contour (decreasing tau from the tip), the upper wall along the lower contour
        out.stations.push_back(last);
        for (WallSide* w : {&out.lower, &out.upper}) {
            w->surface.push_back(surface_of(tail_branch));
            w->parameter.push_back(tail_u);
            w->u.push_back(tail_u);
            w->hidden.push_back(0);
        }
        out.lower.contour.push_back(np_mod(tau_tip - tau(tail_branch, tail_u), total));
        out.upper.contour.push_back(np_mod(tau(tail_branch, tail_u) - tau_tip, total));
    }
    if (touching && out.stations.size() >= 3) {
        // the touching nodes' labels: each wall keeps the surface it has on the station next to the
        // touching one (normally lower suction, upper pressure), as if the extreme point were the
        // shared end of both surfaces; so the unchanged Kutta condition and post-processing read
        // one suction and one pressure node there.  u and parameter stay the extreme point's u on
        // the curve it lies on (Grid::touching_surface, touching_u)
        const std::size_t last_k = out.stations.size() - 1;
        for (WallSide* w : {&out.lower, &out.upper}) {
            w->surface.front() = w->surface[1];
            w->surface.back() = w->surface[last_k - 1];
        }
    }
    // geometry.py's LE point (tau = L_s) and TE point (tau = 0 = total) on each contour
    std::array<std::vector<double>, 2> splits;
    const double span = np_mod(tau_tail - tau_tip, total);
    for (double t : {lengths[0], 0.0}) {
        if (np_mod(t - tau_tip, total) <= span)
            splits[1].push_back(np_mod(t - tau_tip, total));
        else
            splits[0].push_back(np_mod(tau_tip - t, total));
    }
    check_walls(out.lower, out.upper, splits, touch, first);
    return out;
}

const CurveRow::Pieces& CurveRow::pieces() const
{
    if (!pieces_.empty()) return pieces_.front();
    std::vector<double> vertices;
    vertices.insert(vertices.end(), m0.begin(), m0.end());
    vertices.insert(vertices.end(), m1.begin(), m1.end());
    for (const auto* curve : curves)
        for (double u : curve->interior_knots()) vertices.push_back(curve->evaluate(u).m);
    vertices = unique(std::move(vertices));
    Pieces p;
    std::vector<double> middle;
    for (std::size_t k = 0; k + 1 < vertices.size(); ++k) {
        p.lo.push_back(vertices[k]);
        p.hi.push_back(vertices[k + 1]);
        middle.push_back(0.5 * (vertices[k] + vertices[k + 1]));
    }
    const Cut c = cut(middle);
    for (std::size_t k = 0; k < middle.size(); ++k) {
        p.top.push_back(c.top_segment(k));
        p.bottom.push_back(c.segment[k][0]);
    }
    pieces_.push_back(std::move(p));
    return pieces_.front();
}

std::array<double, 2> CurveRow::widest_cut() const
{
    const Pieces& p = pieces();
    const std::size_t n = p.lo.size();
    std::vector<double> span(2 * n), ends(2 * n);
    for (std::size_t k = 0; k < 2 * n; ++k) {
        const std::size_t q = k % n;
        ends[k] = k < n ? p.lo[q] : p.hi[q];
        span[k] = smooth(p.top[q], ends[k]) - smooth(p.bottom[q], ends[k]);
    }
    const auto order = stable_argsort(span);
    return {span[order.back()], ends[order.back()]};
}

std::array<double, 2> CurveRow::narrowest(double pitch) const
{
    const Pieces& P = pieces();
    std::vector<std::array<double, 2>> found;  // (width, m) of every evaluation round
    struct Rows
    {
        std::vector<int> piece;
        std::vector<std::vector<double>> m, upper, lower, width;
    };
    auto evaluate = [&](Rows& r) {
        r.upper.assign(r.m.size(), {});
        r.lower.assign(r.m.size(), {});
        r.width.assign(r.m.size(), {});
        double best = 0.0, best_m = 0.0;
        bool first = true;
        for (std::size_t i = 0; i < r.m.size(); ++i) {
            const int piece = r.piece[i];
            for (double m : r.m[i]) {
                const double up = smooth(P.top[sz(piece)], m), low = smooth(P.bottom[sz(piece)], m);
                const double w = pitch - (up - low);
                r.upper[i].push_back(up);
                r.lower[i].push_back(low);
                r.width[i].push_back(w);
                // np.argsort(width, axis=None, kind="stable")[0]: the first smallest, NaN last
                if (!std::isnan(w) && (first || w < best)) {
                    best = w;
                    best_m = m;
                    first = false;
                }
            }
        }
        if (first && !r.m.empty() && !r.m[0].empty()) {  // every width NaN
            best = r.width[0][0];
            best_m = r.m[0][0];
        }
        found.push_back({best, best_m});
    };
    auto doubtful = [&](const Rows& r, std::size_t i, std::size_t c) {
        const double bound = pitch - std::max(r.upper[i][c], r.upper[i][c + 1]) + std::min(r.lower[i][c], r.lower[i][c + 1]);
        return !(bound > touch);
    };
    Rows ends;
    for (std::size_t k = 0; k < P.lo.size(); ++k) {
        ends.piece.push_back(static_cast<int>(k));
        ends.m.push_back({P.lo[k], P.hi[k]});
    }
    evaluate(ends);
    std::vector<int> doubt;
    for (std::size_t i = 0; i < ends.m.size(); ++i)
        if (doubtful(ends, i, 0)) doubt.push_back(ends.piece[i]);
    if (!doubt.empty()) {
        const auto fraction = linspace(0.0, 1.0, kSamples + 1);
        Rows cells;
        for (int piece : doubt) {
            cells.piece.push_back(piece);
            std::vector<double> m;
            for (double f : fraction) m.push_back(P.lo[sz(piece)] + f * (P.hi[sz(piece)] - P.lo[sz(piece)]));
            cells.m.push_back(std::move(m));
        }
        evaluate(cells);
        Rows zoom;
        std::vector<double> a, b;
        for (std::size_t i = 0; i < cells.m.size(); ++i)
            for (std::size_t c = 0; c + 1 < cells.m[i].size(); ++c)
                if (doubtful(cells, i, c)) {
                    zoom.piece.push_back(cells.piece[i]);
                    a.push_back(cells.m[i][c]);
                    b.push_back(cells.m[i][c + 1]);
                }
        for (int round = 0; round < (zoom.piece.empty() ? 0 : kZooms); ++round) {
            zoom.m.assign(a.size(), {});
            for (std::size_t i = 0; i < a.size(); ++i)
                for (double f : fraction) zoom.m[i].push_back(a[i] + f * (b[i] - a[i]));
            evaluate(zoom);
            for (std::size_t i = 0; i < a.size(); ++i) {
                const auto order = stable_argsort(zoom.width[i]);
                const int k = static_cast<int>(order[0]);
                a[i] = zoom.m[i][sz(std::max(k - 1, 0))];
                b[i] = zoom.m[i][sz(std::min(k + 1, kSamples))];
            }
        }
    }
    std::vector<double> widths;
    for (const auto& f : found) widths.push_back(f[0]);
    return found[stable_argsort(widths)[0]];
}

}  // namespace BladeGen::Turbolenza::detail
