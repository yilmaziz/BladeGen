// Port of v1 metrics.py (see Metrics.hpp): build_metrics, _wall_derivatives, _curve_slopes,
// _metric_fields, _uniform_spacing, _runs, _d_dm, _d_dm_run.
#include "turbolenza/Metrics.hpp"

#include "turbolenza/Detail.hpp"
#include "turbolenza/StreamSheet.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>
#include <tuple>

namespace BladeGen::Turbolenza {

using detail::sz;

double singular_cos() noexcept
{
    static const double value = std::sqrt(kEpsilon);
    return value;
}

namespace {

// v1 _uniform_spacing: the constant step of strictly increasing, uniformly spaced values.
double uniform_spacing(const std::vector<double>& values, const char* name, double scale)
{
    const std::size_t n = values.size();
    for (std::size_t k = 0; k + 1 < n; ++k)
        if (!(values[k + 1] - values[k] > 0.0))
            throw MetricsError(std::format("{} is not strictly increasing at index {} ({:.17g} -> {:.17g})", name, k,
                                           values[k], values[k + 1]));
    const double step = (values[n - 1] - values[0]) / static_cast<double>(n - 1);
    double worst = -1.0;
    std::size_t at = 0;
    for (std::size_t k = 0; k + 1 < n; ++k) {
        const double off = std::abs((values[k + 1] - values[k]) - step);
        if (off > worst) {
            worst = off;
            at = k;
        }
    }
    if (!(worst <= kRoundOff * scale))
        throw MetricsError(std::format("{} is not uniformly spaced: step {} is {:.17g}, the mean step {:.17g}", name, at,
                                       values[at + 1] - values[at], step));
    return step;
}

// v1 _d_dm_run on one run [start, stop) of rows, column c of a row-major array with `stride` columns.
void d_dm_run(const std::vector<double>& f, std::vector<double>& d, int stride, int start, int stop, double dm)
{
    const double h2 = 2.0 * dm;
    for (int c = 0; c < stride; ++c) {
        auto at = [&](int i) { return f[sz(i) * sz(stride) + sz(c)]; };
        auto out = [&](int i) -> double& { return d[sz(i) * sz(stride) + sz(c)]; };
        for (int i = start + 1; i < stop - 1; ++i) out(i) = (at(i + 1) - at(i - 1)) / h2;
        out(start) = (-3.0 * at(start) + 4.0 * at(start + 1) - at(start + 2)) / h2;
        out(stop - 1) = (3.0 * at(stop - 1) - 4.0 * at(stop - 2) + at(stop - 3)) / h2;
    }
}

struct WallDerivatives
{
    std::vector<Point> d1_lower, d1_upper, d2_lower, d2_upper;
};

// v1 _wall_derivatives: the curve derivatives of the C2 walls, once order 0 is checked to be the
// wall nodes.
WallDerivatives wall_derivatives(const Grid& g)
{
    const int f = g.first_blade_i, l = g.last_blade_i;
    const int n_blade = l - f + 1;
    const int J = g.n_eta() - 1;
    double scale = std::max({std::abs(g.m.front()), std::abs(g.m.back()), g.m.back() - g.m.front(), g.pitch()});
    for (double y : g.node_y.v) scale = std::max(scale, std::abs(y));
    for (int w = 0; w < 2; ++w) {
        const char* name = w == 0 ? "lower" : "upper";
        const std::vector<Point>& points = w == 0 ? g.wall_d[0].lower : g.wall_d[0].upper;
        const int j = w == 0 ? 0 : J;
        double off = 0.0;
        if (static_cast<int>(points.size()) != n_blade) {
            off = std::numeric_limits<double>::infinity();
        } else {
            for (int b = 0; b < n_blade; ++b) {
                off = std::max(off, std::abs(points[sz(b)].m - g.node_m(f + b, j)));
                off = std::max(off, std::abs(points[sz(b)].y - g.node_y(f + b, j)));
                if (!std::isfinite(points[sz(b)].m) || !std::isfinite(points[sz(b)].y))
                    off = std::numeric_limits<double>::infinity();
            }
        }
        if (!(off <= kRoundOff * scale))
            throw MetricsError(std::format("the {} wall's curve points are up to {:.3g} off its wall nodes: the curve "
                                           "derivatives would not belong to the nodes",
                                           name, off));
    }
    return {g.wall_d[1].lower, g.wall_d[1].upper, g.wall_d[2].lower, g.wall_d[2].upper};
}

struct CurveSlopes
{
    std::vector<double> tan0, tan1, dtan0, dtan1;
};

// v1 _curve_slopes: tan = y_u / m_u, d tan / dm = (y_uu m_u - y_u m_uu) / m_u^3.
CurveSlopes curve_slopes(const WallDerivatives& w, const std::vector<double>& m, int f, int l)
{
    const std::size_t n_blade = sz(l - f + 1);
    if (w.d1_lower.size() != n_blade || w.d1_upper.size() != n_blade || w.d2_lower.size() != n_blade ||
        w.d2_upper.size() != n_blade)
        throw MetricsError(std::format("the C2 wall derivatives must be four (n_blade, 2) = ({}, 2) arrays, got shapes "
                                       "[({}, 2), ({}, 2), ({}, 2), ({}, 2)]",
                                       n_blade, w.d1_lower.size(), w.d1_upper.size(), w.d2_lower.size(),
                                       w.d2_upper.size()));
    for (const auto* a : {&w.d1_lower, &w.d1_upper, &w.d2_lower, &w.d2_upper})
        for (const Point& p : *a)
            if (!std::isfinite(p.m) || !std::isfinite(p.y)) throw MetricsError("the C2 wall derivatives are not finite");
    CurveSlopes out;
    const double sc = singular_cos();
    for (int wall = 0; wall < 2; ++wall) {
        const char* name = wall == 0 ? "lower" : "upper";
        const std::vector<Point>& d1 = wall == 0 ? w.d1_lower : w.d1_upper;
        const std::vector<Point>& d2 = wall == 0 ? w.d2_lower : w.d2_upper;
        std::vector<double> tan(n_blade), dtan(n_blade);
        for (std::size_t k = 0; k < n_blade; ++k) {
            const double m_u = d1[k].m, y_u = d1[k].y, m_uu = d2[k].m, y_uu = d2[k].y;
            const double speed = std::hypot(m_u, y_u);
            if (!(std::abs(m_u) > sc * speed)) {
                const int i = f + static_cast<int>(k);
                const double cosine = speed > 0.0 ? std::abs(m_u) / speed : 0.0;
                throw MetricsError(std::format(
                    "the {} blade wall is numerically vertical at blade station {} (m = {:.10g}): |dm/du| = {:.3g} "
                    "against |dX/du| = {:.3g}, so cos(beta) = {:.3g} is not above {:.3g}. The curve is smooth there, "
                    "but m is no coordinate along it: tan(beta) = y_u / m_u and d tan(beta)/dm = (y_uu m_u - y_u m_uu) "
                    "/ m_u^3 do not exist",
                    name, i, m[sz(i)], std::abs(m_u), speed, cosine, sc));
            }
            tan[k] = y_u / m_u;
            dtan[k] = (y_uu * m_u - y_u * m_uu) / std::pow(m_u, 3.0);
        }
        if (wall == 0) {
            out.tan0 = std::move(tan);
            out.dtan0 = std::move(dtan);
        } else {
            out.tan1 = std::move(tan);
            out.dtan1 = std::move(dtan);
        }
    }
    return out;
}

}  // namespace

std::vector<MRun> m_runs(int n_m, int first_blade_i, int last_blade_i)
{
    if (!(0 <= first_blade_i && first_blade_i <= last_blade_i && last_blade_i < n_m))
        throw MetricsError(
            std::format("blade stations {} to {} do not fit {} stations", first_blade_i, last_blade_i, n_m));
    const std::vector<MRun> runs = {{"upstream periodic", 0, first_blade_i},
                                    {"blade", first_blade_i, last_blade_i + 1},
                                    {"downstream periodic", last_blade_i + 1, n_m}};
    for (const MRun& r : runs) {
        if (r.stop - r.start < 3) {
            const std::string name = r.name;
            const char* cure = name == "blade"             ? "raise n_passage"
                               : name == "upstream periodic" ? "lengthen upstream_extension_chords or raise n_passage"
                                                             : "lengthen downstream_extension_chords or raise n_passage";
            throw MetricsError(std::format("the {} run holds {} station(s), and a second-order m derivative inside one "
                                           "run needs three: {}",
                                           r.name, r.stop - r.start, cure));
        }
    }
    return runs;
}

std::vector<double> d_dm_runs(const std::vector<double>& values, int stride, double dm, const std::vector<MRun>& runs)
{
    std::vector<double> d(values.size(), 0.0);
    for (const MRun& r : runs) d_dm_run(values, d, stride, r.start, r.stop, dm);
    return d;
}

GridMetrics build_metrics(std::shared_ptr<const Grid> grid_ptr)
{
    if (!grid_ptr) throw MetricsError("build_metrics takes a Grid, got None");
    const Grid& g = *grid_ptr;
    const int n_m = g.n_m(), n_eta = g.n_eta();
    if (g.node_m.n_m != n_m || g.node_m.n_eta != n_eta || !g.node_m.same_shape(g.node_y))
        throw MetricsError(std::format("m ({},), eta ({},) and nodes ({}, {}, 2) do not fit: nodes must be (n_m, n_eta, 2)",
                                       n_m, n_eta, g.node_m.n_m, g.node_m.n_eta));
    if (n_eta < 2)
        throw MetricsError(std::format("eta has {} node(s); a passage needs at least its two walls", n_eta));
    for (int i = 0; i < n_m; ++i)
        if (!std::isfinite(g.m[sz(i)])) throw MetricsError(std::format("m is not finite at index ({},)", i));
    for (int j = 0; j < n_eta; ++j)
        if (!std::isfinite(g.eta[sz(j)])) throw MetricsError(std::format("eta is not finite at index ({},)", j));
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j) {
            if (!std::isfinite(g.node_m(i, j)))
                throw MetricsError(std::format("nodes is not finite at index ({}, {}, 0)", i, j));
            if (!std::isfinite(g.node_y(i, j)))
                throw MetricsError(std::format("nodes is not finite at index ({}, {}, 1)", i, j));
        }
    const double pitch = g.pitch();
    if (!std::isfinite(pitch) || !(pitch > 0.0))
        throw MetricsError(std::format("the pitch {:.17g} must be finite and positive", pitch));

    // C2 (v2 always): the curve derivatives at the wall nodes, checked first (v1 build_metrics order)
    if (!(g.first_blade_i >= 0 && g.first_blade_i <= g.last_blade_i && g.last_blade_i < n_m))
        throw MetricsError(std::format("blade stations {} to {} do not fit {} stations", g.first_blade_i,
                                       g.last_blade_i, n_m));
    const WallDerivatives walls = wall_derivatives(g);

    const std::vector<double>& m = g.m;
    const std::vector<double>& eta = g.eta;
    const double m_scale = std::max({std::abs(m.front()), std::abs(m.back()), m.back() - m.front()});
    const double dm = uniform_spacing(m, "m", m_scale);
    const double deta = uniform_spacing(eta, "eta", 1.0);
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (!(std::abs(g.node_m(i, j) - m[sz(i)]) <= kRoundOff * m_scale))
                throw MetricsError(
                    std::format("the nodes of station {} do not lie on its line x = m = {:.17g}", i, m[sz(i)]));
    const int f = g.first_blade_i, l = g.last_blade_i;
    const std::vector<MRun> runs = m_runs(n_m, f, l);

    const int J = n_eta - 1;
    std::vector<double> y0(sz(n_m)), y1(sz(n_m)), S(sz(n_m));
    double scale = pitch;
    for (double y : g.node_y.v) scale = std::max(scale, std::abs(y));
    for (int i = 0; i < n_m; ++i) {
        y0[sz(i)] = g.node_y(i, 0);
        y1[sz(i)] = g.node_y(i, J);
        S[sz(i)] = y1[sz(i)] - y0[sz(i)];
    }
    for (int i = 0; i < n_m; ++i)
        if (!(S[sz(i)] > kRoundOff * scale))
            throw MetricsError(std::format("S = y1 - y0 = {:.6g} at station {} (m = {:.10g}) is not above round-off "
                                           "({:.3g}): the map from (m, eta) collapses or reverses there",
                                           S[sz(i)], i, m[sz(i)], kRoundOff * scale));
    {
        double worst = -1.0;
        int wi = 0, wj = 0;
        bool bad = false;
        for (int i = 0; i < n_m; ++i)
            for (int j = 0; j < n_eta; ++j) {
                const double mismatch = std::abs(y0[sz(i)] + eta[sz(j)] * S[sz(i)] - g.node_y(i, j));
                if (!(mismatch <= kRoundOff * scale)) bad = true;
                if (mismatch > worst || std::isnan(mismatch)) {
                    worst = mismatch;
                    wi = i;
                    wj = j;
                }
            }
        if (bad)
            throw MetricsError(std::format("the node ({}, {}) is {:.3g} off y0 + eta S: the nodes are not the linear eta "
                                           "map between the boundaries",
                                           wi, wj, worst));
    }
    for (const MRun& r : runs) {
        if (std::string(r.name) == "blade") continue;
        double worst = -1.0;
        int at = r.start;
        bool bad = false;
        for (int i = r.start; i < r.stop; ++i) {
            const double drift = std::abs(S[sz(i)] - pitch);
            if (!(drift <= kRoundOff * scale)) bad = true;
            if (drift > worst) {
                worst = drift;
                at = i;
            }
        }
        if (bad)
            throw MetricsError(std::format("S = {:.17g} at {} station {} is not the pitch {:.17g}: the side lines must "
                                           "be one pitch apart",
                                           S[sz(at)], r.name, at, pitch));
    }

    // not in v1: the stream sheet at the stations (planar: the neutral factors, never applied)
    const StreamSheet* sheet = g.boundaries.sheet.get();
    std::vector<double> sigma(sz(n_m), 1.0), radius(sz(n_m), 1.0), thickness(sz(n_m), 1.0), sin_phi(sz(n_m), 0.0);
    if (sheet) {
        try {
            const double b_in = sheet->thickness(m[0]);
            for (int i = 0; i < n_m; ++i) {
                radius[sz(i)] = sheet->radius(m[sz(i)]);
                sigma[sz(i)] = radius[sz(i)] / sheet->r_ref();
                thickness[sz(i)] = sheet->thickness(m[sz(i)]) / b_in;  // Aungier's b, 1 at the far inlet
                sin_phi[sz(i)] = sheet->sin_phi(m[sz(i)]);             // Eq. 5-15
            }
        } catch (const std::out_of_range& error) {
            throw MetricsError(std::string("the grid leaves its stream sheet: ") + error.what());
        }
    }

    std::vector<double> tan0 = d_dm_runs(y0, 1, dm, runs), tan1 = d_dm_runs(y1, 1, dm, runs);
    if (sheet)  // Eq. 5-3: tan(beta) = r d theta / dm = sigma dy/dm
        for (int i = 0; i < n_m; ++i) {
            tan0[sz(i)] = sigma[sz(i)] * tan0[sz(i)];
            tan1[sz(i)] = sigma[sz(i)] * tan1[sz(i)];
        }
    const CurveSlopes exact = curve_slopes(walls, m, f, l);
    for (int i = f; i <= l; ++i) {
        tan0[sz(i)] = exact.tan0[sz(i - f)];
        tan1[sz(i)] = exact.tan1[sz(i - f)];
        if (sheet) {
            tan0[sz(i)] = sigma[sz(i)] * exact.tan0[sz(i - f)];
            tan1[sz(i)] = sigma[sz(i)] * exact.tan1[sz(i - f)];
        }
    }
    std::vector<double> dtan0 = d_dm_runs(tan0, 1, dm, runs), dtan1 = d_dm_runs(tan1, 1, dm, runs);
    for (int i = f; i <= l; ++i) {
        dtan0[sz(i)] = exact.dtan0[sz(i - f)];
        dtan1[sz(i)] = exact.dtan1[sz(i - f)];
        if (sheet) {  // d(sigma y')/dm = sigma' y' + sigma y'', sigma' = sin(phi) / r_ref
            const double dsigma = sin_phi[sz(i)] / sheet->r_ref();
            dtan0[sz(i)] = dsigma * exact.tan0[sz(i - f)] + sigma[sz(i)] * exact.dtan0[sz(i - f)];
            dtan1[sz(i)] = dsigma * exact.tan1[sz(i - f)] + sigma[sz(i)] * exact.dtan1[sz(i - f)];
        }
    }
    std::vector<double> dS_dm(sz(n_m));
    for (int i = 0; i < n_m; ++i) dS_dm[sz(i)] = tan1[sz(i)] - tan0[sz(i)];
    if (sheet)  // S = r (theta1 - theta0) [Eq. 5-9]: dS/dm = sin(phi) Theta + r (theta1' - theta0')
        for (int i = 0; i < n_m; ++i) {
            dS_dm[sz(i)] = sin_phi[sz(i)] * (S[sz(i)] / sheet->r_ref()) + (tan1[sz(i)] - tan0[sz(i)]);
            S[sz(i)] = sigma[sz(i)] * S[sz(i)];
        }

    GridMetrics out;
    out.grid = grid_ptr;
    out.tan_beta = Field(n_m, n_eta);
    out.beta = Field(n_m, n_eta);
    out.sin_beta = Field(n_m, n_eta);
    out.cos_beta = Field(n_m, n_eta);
    out.dtan_beta_dm = Field(n_m, n_eta);
    const double sc = singular_cos();
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j) {
            const double e = eta[sz(j)];
            const double t = (1.0 - e) * tan0[sz(i)] + e * tan1[sz(i)];
            out.tan_beta(i, j) = t;
            out.cos_beta(i, j) = 1.0 / std::hypot(1.0, t);
        }
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j) {
            const double c = out.cos_beta(i, j);
            if (!(c > sc))
                throw MetricsError(std::format("cos(beta) = {:.3g} at node ({}, {}), beta = {:.10g} deg: the grid line "
                                               "there is numerically vertical (tan(beta) = {:.3g}) and 1/cos(beta) "
                                               "singular",
                                               c, i, j, detail::degrees(std::atan(out.tan_beta(i, j))),
                                               out.tan_beta(i, j)));
        }
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j) {
            const double e = eta[sz(j)];
            const double t = out.tan_beta(i, j);
            out.dtan_beta_dm(i, j) = (1.0 - e) * dtan0[sz(i)] + e * dtan1[sz(i)];
            out.beta(i, j) = std::atan(t);
            out.sin_beta(i, j) = t * out.cos_beta(i, j);
        }
    out.dS_dm = std::move(dS_dm);
    out.radius_ratio = std::move(sigma);
    out.radius = std::move(radius);
    out.thickness = std::move(thickness);
    out.sin_phi = std::move(sin_phi);
    out.on_sheet = sheet != nullptr;
    out.S = std::move(S);
    out.tan_beta0 = std::move(tan0);
    out.tan_beta1 = std::move(tan1);
    out.dtan_beta0_dm = std::move(dtan0);
    out.dtan_beta1_dm = std::move(dtan1);

    for (const auto& [name, spacing, expected] :
         {std::tuple<const char*, double, double>{"dm", g.dm(), dm}, {"deta", g.deta(), deta}})
        if (!(std::abs(spacing - expected) <= kRoundOff * std::abs(expected)))
            throw MetricsError(
                std::format("grid.{} = {:.17g} is not the spacing of its stations, {:.17g}", name, spacing, expected));
    return out;
}

}  // namespace BladeGen::Turbolenza
