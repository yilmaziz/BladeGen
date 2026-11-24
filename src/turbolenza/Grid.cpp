// Grid generation for one blade passage, the C2 path of v1 grid.py sections 1-8: the passage
// boundaries (make_passage_boundaries), the uniform stations, the blade-passage topology, the
// wall nodes, the Grid invariants, rebuild_exit and the admissible n_passage (Aungier 2003
// Secs. 5.2, 5.3, Figs. 5-3, 5-5).  The standard levels (section 9) are in GridStandard.cpp, the
// blade row in GridRow.cpp.
//
// The EdgeNodes::Touching EXPERIMENT (not in v1; see Grid.hpp) lives here too: the side lines
// anchored at the nose tip / tail end (make_passage_boundaries), the touching grid
// (grid_on_touching_stations, touching_lattice), the relaxed Grid invariants of its two touching
// stations (finish_grid) and its admissibility rule (validate_uniform_grid_topology and friends).
// The Inside path (v1) is unchanged: every touching branch is guarded by the layout.
#include "turbolenza/Grid.hpp"

#include "turbolenza/GridDetail.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace BladeGen::Turbolenza {

namespace {

using detail::CurveRow;
using detail::Cut;
using detail::fmt_g;

std::size_t sz(int k) { return static_cast<std::size_t>(k); }

constexpr double kPi = 3.141592653589793;

std::string repr(double x)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", x);
    return buf;
}

// ── input checks (v1 _check_angle, _check_blade_angle, _positive_number, _whole_number) ──

void check_angle(double angle, const char* name)
{
    if (!(std::isfinite(angle) && std::abs(angle) < 0.5 * kPi))
        throw std::invalid_argument(std::string(name) + " must be a finite angle strictly between -pi/2 and pi/2 radians, got " +
                                    repr(angle) + "; convert an angle in degrees with math.radians");
}

void check_blade_angle(double kappa, const char* name, const char* side)
{
    if (!(std::abs(kappa) < 0.5 * kPi))
        throw std::invalid_argument(std::string(name) + " = " + fmt_g(kappa * (180.0 / kPi)) + " deg: the " + side +
                                    " side boundary cannot turn to a blade angle at or beyond 90 deg from axial");
}

double positive_number(double value, const std::string& name)
{
    if (!(std::isfinite(value) && value > 0.0))
        throw std::invalid_argument(name + " must be a positive number of chords, got " + repr(value));
    return value;
}

int whole_number(std::optional<int> value, const char* name, int minimum)
{
    if (!value || *value < minimum)
        throw std::invalid_argument(std::string(name) + " must be a whole number of at least " + std::to_string(minimum) +
                                    ", got " + (value ? std::to_string(*value) : std::string("None")));
    return *value;
}

// ── side lines (v1 _samples, _sample_fraction, _inlet_side, _exit_side) ─────

std::vector<double> samples(double start, double stop, int n_passage, double blade_length)
{
    const long long n_cells = std::max(1LL, detail::py_round((stop - start) * static_cast<double>(n_passage) / blade_length));
    return detail::linspace(start, stop, static_cast<int>(n_cells) + 1);
}

std::vector<double> sample_fraction(std::size_t n_samples, bool junction_last)
{
    std::vector<double> out(n_samples);
    const double n = static_cast<double>(n_samples) - 1.0;
    for (std::size_t k = 0; k < n_samples; ++k) {
        const double kk = static_cast<double>(k);
        out[k] = junction_last ? (n - kk) / n : kk / n;
    }
    return out;
}

/// The side line sampled at m: beta and the lower boundary points.
void side_samples(const std::vector<double>& m, bool junction_last, const Point& junction, double junction_angle,
                  double far_angle, std::vector<double>& beta, std::vector<Point>& lower)
{
    const auto fraction = sample_fraction(m.size(), junction_last);
    beta.resize(m.size());
    lower.resize(m.size());
    for (std::size_t k = 0; k < m.size(); ++k) {
        const auto by = detail::side_line(m[k], fraction[k], junction, junction_angle, far_angle);
        beta[k] = by[0];
        lower[k] = {m[k], by[1]};
    }
}

std::vector<Point> translated(const std::vector<Point>& p, const Point& t)
{
    std::vector<Point> out(p.size());
    for (std::size_t k = 0; k < p.size(); ++k) out[k] = {p[k].m + t.m, p[k].y + t.y};
    return out;
}

// The station fraction is measured from the side line's own junction (its last / first sample).
// Inside (v1) that junction is the LE / TE point, whose m is b.m_le / b.m_te exactly (linspace
// ends there exactly), so v1's fraction (m - m_le) / (m_in - m_le) is unchanged bit for bit;
// Touching it is the nose tip / tail end, so the lines and rebuild_exit keep that anchor.
double inlet_side(const PassageBoundaries& b, double m)
{
    const Point& junction = b.upstream_lower_boundary.back();
    const double fraction = (m - junction.m) / (b.m_in - junction.m);
    return detail::side_line(m, fraction, junction, b.upstream_beta.back(), b.inlet_angle)[1];
}

double exit_side(const PassageBoundaries& b, double m)
{
    const Point& junction = b.downstream_lower_boundary.front();
    const double fraction = (m - junction.m) / (b.m_out - junction.m);
    return detail::side_line(m, fraction, junction, b.downstream_beta.front(), b.exit_angle)[1];
}

/// The touching stations' assigned edge direction (not in v1; GridData.hpp, Grid::wall_d): order
/// k = 1 or 2 of the side line at its blade end, parametrized by arc length.  exit: the TE station.
Point touching_edge_derivative(const PassageBoundaries& b, bool exit, int k)
{
    const double kappa = exit ? b.kappa_te : b.kappa_le;
    const double far = exit ? b.exit_angle : b.inlet_angle;
    const double m_far = exit ? b.m_out : b.m_in;
    const double m_junction = exit ? b.downstream_lower_boundary.front().m : b.upstream_lower_boundary.back().m;
    const double c = std::cos(kappa), s = std::sin(kappa);
    if (k == 1) return {c, s};
    const double dbeta_ds = c * (far - kappa) / (m_far - m_junction);
    return {-dbeta_ds * s, dbeta_ds * c};
}

// ── section 8 helpers ────────────────────────────────────────────────────────

EdgeStationTopology edge_topology(const CurveRow& row, const PassageBoundaries& b, const Cut& cut, std::size_t k,
                                  const char* edge, int station)
{
    EdgeStationTopology t;
    t.edge = edge;
    t.station = station;
    t.m = cut.m[k];
    t.crossings = cut.count[k];
    for (int j = 0; j < t.crossings; ++j)
        (row.surface_of(cut.segment[k][sz(j)]) == Surface::Suction ? t.suction : t.pressure) += 1;
    const auto ext = row.extremes();
    const bool le = std::string(edge) == "LE";
    t.branch_point = le ? b.m_le : b.m_te;
    t.extreme = le ? ext[0] : ext[1];
    t.on_branch_point = std::abs(t.m - t.branch_point) <= row.touch;
    return t;
}

/// v1 _lattice_topology: (LE, TE) of the uniform grid of n_passage, or empty when make_grid would
/// refuse its stations.
std::optional<std::array<EdgeStationTopology, 2>> lattice_topology(const PassageBoundaries& b, const CurveRow& row,
                                                                   int n_passage)
{
    const auto m = detail::stations(b, n_passage);
    const Cut cut = row.cut(m);
    std::vector<int> blade;
    for (std::size_t i = 0; i < m.size(); ++i)
        if (cut.cuts_blade(i)) blade.push_back(static_cast<int>(i));
    if (blade.empty() || static_cast<int>(blade.size()) != blade.back() - blade.front() + 1 || blade.front() == 0 ||
        blade.back() == static_cast<int>(m.size()) - 1)
        return std::nullopt;
    return std::array<EdgeStationTopology, 2>{edge_topology(row, b, cut, sz(blade.front()), "LE", blade.front()),
                                              edge_topology(row, b, cut, sz(blade.back()), "TE", blade.back())};
}

/// v1 _ceiling.
std::optional<int> ceiling_of(const PassageBoundaries& b, const CurveRow& row)
{
    const auto ext = row.extremes();
    std::vector<double> windows;
    for (double w : {b.m_le - ext[0], ext[1] - b.m_te})
        if (w > row.touch) windows.push_back(w);
    if (windows.empty()) return std::nullopt;
    const double width = b.m_out - b.m_in, length = b.m_te - b.m_le;
    auto spacing = [&](long long n) {
        return width / static_cast<double>(std::max(1LL, detail::py_round(width / (length / static_cast<double>(n)))));
    };
    std::optional<int> found;
    for (double w : windows) {
        long long n = std::max(1LL, static_cast<long long>(length / w) - 2);
        while (spacing(n) > w) ++n;
        while (n > 1 && spacing(n - 1) <= w) --n;
        found = found ? std::min(*found, static_cast<int>(n)) : static_cast<int>(n);
    }
    return found;
}

// ── EdgeNodes::Touching (not in v1) ──────────────────────────────────────────

/// The stations next to the touching ones on the touching lattice of n_passage cells
/// (m_tip + dm and m_tail - dm), read as v1 reads an edge station ("LE" and "TE").  Empty when
/// n_passage < 2 (no station between the touching ones).
std::optional<std::array<EdgeStationTopology, 2>> touching_topology(const PassageBoundaries& b, const CurveRow& row,
                                                                    int n_passage, int first = 0)
{
    if (n_passage < 2) return std::nullopt;
    const double dm = (row.m_tail - row.m_tip) / static_cast<double>(n_passage);
    const Cut cut = row.cut({row.m_tip + dm, row.m_tail - dm});
    return std::array<EdgeStationTopology, 2>{edge_topology(row, b, cut, 0, "LE", first + 1),
                                              edge_topology(row, b, cut, 1, "TE", first + n_passage - 1)};
}

bool touching_admissible(const PassageBoundaries& b, const CurveRow& row, int n_passage)
{
    const auto t = touching_topology(b, row, n_passage);
    return t && (*t)[0].valid() && (*t)[1].valid();
}

/// The touching ceiling: the first n_passage with dm = (m_tail - m_tip) / n_passage <= the wider
/// of v1's windows (m_le - nose tip, tail end - m_te); from there on the station next to a
/// touching one lies in its window.  Empty without a window above TOUCH.
std::optional<int> touching_ceiling(const PassageBoundaries& b, const CurveRow& row)
{
    double width = 0.0;
    for (double w : {b.m_le - row.m_tip, row.m_tail - b.m_te})
        if (w > row.touch) width = std::max(width, w);
    if (!(width > 0.0)) return std::nullopt;
    const double length = row.m_tail - row.m_tip;
    long long n = std::max(1LL, static_cast<long long>(length / width) - 2);
    while (length / static_cast<double>(n) > width) ++n;
    while (n > 1 && length / static_cast<double>(n - 1) <= width) --n;
    return static_cast<int>(n);
}

std::string describe_edge(const EdgeStationTopology& t, bool touching = false)
{
    const bool le = t.edge == "LE";
    const char* where = touching ? (le ? "Station next to the nose tip" : "Station next to the tail end") : le ? "First" : "Last";
    const char* point = le ? "LE point" : "TE point";
    const char* extreme = le ? "nose tip" : "tail end";
    const bool inside = le ? (t.extreme < t.m && t.m < t.branch_point) : (t.branch_point < t.m && t.m < t.extreme);
    const char* interval = le ? "the nose-return interval between the nose tip and the LE point"
                              : "the tail-return interval between the TE point and the tail end";
    std::string text = std::string(where) + (touching ? " i = " : " blade station i = ") + std::to_string(t.station) +
                       ": m = " + fmt_g(t.m) + "; " +
                       extreme + " m = " + fmt_g(t.extreme) + "; " + point +
                       " (suction / pressure branch point) m = " + fmt_g(t.branch_point) + ". ";
    if (t.on_branch_point)
        return text + "The station line passes through the " + point +
               " itself, where the two surfaces meet, so its wall nodes are not one suction and one pressure node.";
    if (inside) {
        const char* surface = t.suction >= t.pressure ? "suction" : "pressure";
        return text + "The station lies inside " + interval + ", so its line crosses the " + surface + " surface twice (" +
               std::to_string(t.suction) + " suction, " + std::to_string(t.pressure) +
               " pressure crossings) instead of producing one suction and one pressure wall node.";
    }
    return text + "Its line crosses the blade " + std::to_string(t.crossings) + " times (" + std::to_string(t.suction) +
           " suction, " + std::to_string(t.pressure) + " pressure) instead of once on each surface.";
}

/// The touching layout's admissibility (not in v1; the counterpart of v1
/// validate_uniform_grid_topology): stations f and l on the nose tip and the tail end (to TOUCH),
/// a blade station between them, and the stations f + 1 and l - 1 each cutting one suction and
/// one pressure surface away from the LE / TE point.  GridTopologyError otherwise (rejects, never
/// repairs), with the nearest admissible n_passage of the touching lattice.
void validate_touching(const PassageBoundaries& b, const CurveRow& row, const std::vector<double>& m, int f, int l,
                       std::optional<int> n_passage, std::optional<double> stagger)
{
    const std::string named = stagger ? " (stagger " + fmt_g(*stagger * (180.0 / kPi)) + " deg)" : "";
    const std::string head = (n_passage ? "n_passage = " + std::to_string(*n_passage) + " is not admissible for this cascade"
                                        : std::string("this touching grid is not admissible for its cascade")) +
                             named + " in the touching layout (EXPERIMENT, Aungier 2003 Fig. 5-5)";
    const int n = static_cast<int>(m.size());
    if (!(1 <= f && l - f >= 2 && l <= n - 2))
        throw GridTopologyError(head + ": touching stations " + std::to_string(f) + " and " + std::to_string(l) + " of " +
                                    std::to_string(n) + " need a blade station between them and periodic stations on both sides",
                                "LE and TE", n_passage, std::nullopt, std::nullopt, std::nullopt, std::nullopt);
    const double off_tip = std::abs(m[sz(f)] - row.m_tip), off_tail = std::abs(m[sz(l)] - row.m_tail);
    if (!(off_tip <= row.touch && off_tail <= row.touch))
        throw GridTopologyError(head + ": the touching stations must lie on the nose tip m = " + fmt_g(row.m_tip, 17) +
                                    " and the tail end m = " + fmt_g(row.m_tail, 17) + ", not at m = " + fmt_g(m[sz(f)], 17) +
                                    " and " + fmt_g(m[sz(l)], 17),
                                off_tip <= row.touch ? "TE" : (off_tail <= row.touch ? "LE" : "LE and TE"), n_passage,
                                std::nullopt, std::nullopt, std::nullopt, std::nullopt);
    const Cut cut = row.cut({m[sz(f + 1)], m[sz(l - 1)]});
    const EdgeStationTopology le = edge_topology(row, b, cut, 0, "LE", f + 1);
    const EdgeStationTopology te = edge_topology(row, b, cut, 1, "TE", l - 1);
    if (le.valid() && te.valid()) return;
    const auto ceiling = touching_ceiling(b, row);
    std::optional<int> finest, lower, higher;
    if (ceiling)
        for (int k = *ceiling - 1; k > 1; --k)
            if (touching_admissible(b, row, k)) {
                finest = k;
                break;
            }
    if (n_passage) {
        const int np = *n_passage;
        if (ceiling && np >= *ceiling) {
            lower = finest;
        } else {
            for (int k = np - 1; k > 1; --k)
                if (touching_admissible(b, row, k)) {
                    lower = k;
                    break;
                }
            const int top = ceiling ? *ceiling : 4 * np + 100;
            for (int k = np + 1; k < top; ++k)
                if (touching_admissible(b, row, k)) {
                    higher = k;
                    break;
                }
        }
    }
    std::string edge;
    std::string text = head + ": the stations next to the touching ones (first_blade_i + 1 and last_blade_i - 1) must "
                              "each cut one suction and one pressure surface.";
    for (const EdgeStationTopology* t : {&le, &te})
        if (!t->valid()) {
            edge += (edge.empty() ? "" : " and ") + t->edge;
            text += "\n" + describe_edge(*t, true);
        }
    auto opt = [](const std::optional<int>& v) { return v ? std::to_string(*v) : std::string("none"); };
    text += "\nChoose another admissible n_passage; none is substituted.";
    if (n_passage) text += " Nearest admissible: lower " + opt(lower) + ", higher " + opt(higher) + ".";
    if (ceiling)
        text += " Finest admissible n_passage: " + opt(finest) + "; every n_passage from " + std::to_string(*ceiling) +
                " on fails (dm = (m_tail - m_tip) / n_passage no wider than the wider LE / TE window).";
    throw GridTopologyError(text, edge, n_passage, lower, higher, finest, ceiling);
}

PassageBoundaries boundaries_for(const Cascade& cascade, const GridSettings& settings, double inlet_angle)
{
    return make_passage_boundaries(cascade, settings, inlet_angle, std::nullopt);
}

}  // namespace

// ── detail: shared with GridStandard.cpp and the tests ───────────────────────

namespace detail {

double mean_tan(double a, double b)
{
    // v1 _mean_tan: ln(cos a / cos b) / (b - a) = q log1p(x) / x, x = d q, d = b - a,
    // q = tan(b) sinc(d) - sin(d / 2) sinc(d / 2), with numpy's sinc(z) = sin(pi z) / (pi z)
    auto sinc = [](double z) {
        const double y = kPi * (z == 0.0 ? 1.0e-20 : z);
        return std::sin(y) / y;
    };
    const double d = b - a;
    const double q = std::tan(b) * sinc(d / kPi) - std::sin(0.5 * d) * sinc(0.5 * d / kPi);
    const double x = d * q;
    const double ratio = x != 0.0 ? std::log1p(x) / x : 1.0;
    return q * ratio;
}

std::array<double, 2> side_line(double m, double fraction, const Point& junction, double junction_angle, double far_angle)
{
    // v1 _side_line: beta linear in m from the nearer end, y the exact primitive of tan(beta)
    // (Aungier 2003 Sec. 5.3, p. 86: the far-station grid lines along the flow angle, Eq. 5-33)
    const double turn = far_angle - junction_angle;
    const double beta = fraction <= 0.5 ? junction_angle + fraction * turn : far_angle - (1.0 - fraction) * turn;
    return {beta, junction.y + (m - junction.m) * mean_tan(beta, junction_angle)};
}

std::vector<double> stations(const PassageBoundaries& b, int n_passage)
{
    const double dm_target = (b.m_te - b.m_le) / static_cast<double>(n_passage);
    return linspace(b.m_in, b.m_out, static_cast<int>(std::max(1LL, py_round((b.m_out - b.m_in) / dm_target))) + 1);
}

void finish_grid(Grid& g, const CurveRow& row) { finish_grid(g, row, row.cut(g.m)); }

void finish_grid(Grid& g, const CurveRow& row, const Cut& cut)
{
    // v1 Grid.__post_init__
    const int n = g.n_m(), ne = g.n_eta(), f = g.first_blade_i, l = g.last_blade_i;
    const bool touching = g.touching();
    bool topology = 1 <= f && f <= l && l < n - 1 && (!touching || l - f >= 2);
    for (int i = 0; topology && i < n; ++i) {
        // RELAXED for the two touching stations (Touching only): v1's TOUCH rule makes a station
        // line whose cut of the blade is no longer than TOUCH (a line through the nose tip or the
        // tail end) a periodic station; here those two lines are blade stations by definition
        if (g.touching_station(i)) continue;
        topology = cut.cuts_blade(sz(i)) == (i >= f && i <= l);
    }
    if (!topology)
        throw GridError("stations " + std::to_string(f) + " to " + std::to_string(l) +
                        " as blade stations break the blade-passage topology: exactly the stations whose line cuts the "
                        "blade must be blade stations, with periodic stations on both sides");
    int shut = 0, first_shut = -1;
    for (int i = 0; i < n; ++i)
        if (!(g.node_y(i, ne - 1) - g.node_y(i, 0) > 0.0)) {
            if (first_shut < 0) first_shut = i;
            ++shut;
        }
    if (shut)
        throw GridError("the passage width is not positive at " + std::to_string(shut) + " station(s), first at m = " +
                        fmt_g(g.m[sz(first_shut)]));
    const double ty = g.boundaries.periodic_translation.y;
    if (touching) {
        // the touching stations lie on the nose tip and the tail end (to TOUCH)
        const double off_tip = std::abs(g.m[sz(f)] - row.m_tip), off_tail = std::abs(g.m[sz(l)] - row.m_tail);
        if (!(off_tip <= row.touch && off_tail <= row.touch))
            throw GridError("the touching stations " + std::to_string(f) + " (m = " + fmt_g(g.m[sz(f)], 17) + ") and " +
                            std::to_string(l) + " (m = " + fmt_g(g.m[sz(l)], 17) + ") must lie on the nose tip m = " +
                            fmt_g(row.m_tip, 17) + " and the tail end m = " + fmt_g(row.m_tail, 17) + " (to TOUCH)");
    }
    for (int i = f; i <= l; ++i) {
        if (g.touching_station(i)) {
            // RELAXED (Touching only): the two wall nodes of a touching station are one point of the
            // blade, the upper node that point moved up one pitch (zero blade thickness)
            const double gap = std::abs(g.node_y(i, 0) - (g.node_y(i, ne - 1) - ty));
            if (!(gap <= row.touch))
                throw GridError("the wall nodes of the touching station " + std::to_string(i) + " (m = " +
                                fmt_g(g.m[sz(i)]) + ") are " + fmt_g(gap, 3) +
                                " apart (modulo the pitch): both must be the one blade point the station touches");
            continue;
        }
        if (!(g.node_y(i, 0) - (g.node_y(i, ne - 1) - ty) > 0.0))
            throw GridError("the wall nodes of blade station " + std::to_string(i) + " (m = " + fmt_g(g.m[sz(i)]) +
                            ") are not two distinct points of the blade's cut: the lower node, its top, must lie above "
                            "the upper node moved down one pitch, its bottom");
    }
    g.wall_identity = row.identity(cut, f, l, touching);
    const SurfaceCurves& curves = *g.boundaries.blade_curves;
    const Point shift = g.boundaries.periodic_translation;
    if (touching) {
        g.touching_surface = {row.surface_of(row.tip_branch), row.surface_of(row.tail_branch)};
        g.touching_u = {row.tip_u, row.tail_u};
    }
    for (int k = 0; k < 3; ++k) {
        auto side = [&](const WallSide& w, bool upper) {
            std::vector<Point> out(w.u.size());
            for (std::size_t b = 0; b < w.u.size(); ++b) {
                const auto& curve = w.surface[b] == Surface::Suction ? curves.suction : curves.pressure;
                out[b] = curve.derivative(w.u[b], k);
                if (k == 0 && upper) out[b] = {out[b].m + shift.m, out[b].y + shift.y};
            }
            return out;
        };
        g.wall_d[k].lower = side(g.wall_identity.lower, false);
        g.wall_d[k].upper = side(g.wall_identity.upper, true);
        if (touching) {
            // the touching stations (GridData.hpp, Grid::wall_d): k = 0 the touching point on the
            // curve it lies on; k = 1, 2 the mesh's edge direction, the side line at its end
            for (int e = 0; e < 2; ++e) {
                const std::size_t b = e == 0 ? 0 : sz(l - f);
                Point d;
                if (k == 0) {
                    const auto& curve = g.touching_surface[sz(e)] == Surface::Suction ? curves.suction : curves.pressure;
                    d = curve.derivative(g.touching_u[sz(e)], 0);
                } else {
                    d = touching_edge_derivative(g.boundaries, e == 1, k);
                }
                g.wall_d[k].lower[b] = d;
                g.wall_d[k].upper[b] = k == 0 ? Point{d.m + shift.m, d.y + shift.y} : d;
            }
        }
        if (k == 0) {
            // v1 _check_on_curves
            for (int s = 0; s < 2; ++s) {
                const auto& pts = s == 0 ? g.wall_d[0].lower : g.wall_d[0].upper;
                for (int i = f; i <= l; ++i) {
                    const Point& c = pts[sz(i - f)];
                    const double off =
                        std::max(std::abs(c.m - g.m[sz(i)]), std::abs(c.y - g.node_y(i, s == 0 ? 0 : ne - 1)));
                    if (!(off <= row.touch))
                        throw GridError(std::string("the ") + (s == 0 ? "lower" : "upper") + " wall node of blade station " +
                                        std::to_string(i) + " (m = " + fmt_g(g.m[sz(i)]) + ") is " + fmt_g(off, 3) +
                                        " off its authoritative curve point: a C2 wall node must be X(u) of the curve and "
                                        "u that wall_identity records");
                }
            }
        }
    }
}

Grid grid_on_stations(const PassageBoundaries& b, const CurveRow& row, const std::vector<double>& m, int n_eta,
                      const Cut* given)
{
    // v1 _grid_on_stations (sections 3 to 5): the wall nodes on the passage-facing sides of the cut
    // (Aungier 2003 Fig. 5-3) and y = (1 - eta) y_lower + eta y_upper (Aungier 2003 Eq. 5-2)
    const Cut owned = given ? Cut{} : row.cut(m);
    const Cut& cut = given ? *given : owned;
    std::vector<int> blade;
    for (std::size_t i = 0; i < m.size(); ++i)
        if (cut.cuts_blade(i)) blade.push_back(static_cast<int>(i));
    if (blade.empty())
        throw GridError("no station cuts the blade: dm = " + fmt_g((m[1] - m[0]) / b.chord, 4) +
                        " chords lets the blade fall between two stations; raise n_passage");
    const int first = blade.front(), last = blade.back(), n = static_cast<int>(m.size());
    if (static_cast<int>(blade.size()) != last - first + 1)
        throw GridError("the stations that cut the blade, " + std::to_string(first) + " to " + std::to_string(last) +
                        ", are not one run");
    if (first == 0 || last == n - 1) throw GridError("the blade reaches a far station: lengthen the extensions");

    const double ty = b.periodic_translation.y;
    std::vector<double> lower(m.size()), upper(m.size());
    for (int i = 0; i < n; ++i) {
        if (i < first)
            lower[sz(i)] = inlet_side(b, m[sz(i)]);
        else if (i > last)
            lower[sz(i)] = exit_side(b, m[sz(i)]);
        else
            lower[sz(i)] = cut.height[sz(i)][sz(cut.count[sz(i)] - 1)];
        upper[sz(i)] = lower[sz(i)] + ty;
        if (i >= first && i <= last) upper[sz(i)] = cut.height[sz(i)][0] + ty;
    }
    Grid g;
    g.boundaries = b;
    g.m = m;
    g.eta = linspace(0.0, 1.0, n_eta);
    g.node_m = Field(n, n_eta);
    g.node_y = Field(n, n_eta);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n_eta; ++j) {
            const double eta = g.eta[sz(j)];
            g.node_m(i, j) = m[sz(i)];
            g.node_y(i, j) = (1.0 - eta) * lower[sz(i)] + eta * upper[sz(i)];
        }
    g.first_blade_i = first;
    g.last_blade_i = last;
    finish_grid(g, row, cut);
    return g;
}

Grid grid_on_touching_stations(const PassageBoundaries& b, const CurveRow& row, const std::vector<double>& m, int n_eta,
                               int first, int last, const Cut* given)
{
    // the touching layout (not in v1; Aungier 2003 Fig. 5-5 "outside", 2000 Fig. 13-3): stations
    // first and last through the nose tip and the tail end, both wall nodes on that point; the
    // stations between as v1 builds a blade station (top / bottom of the cut, Fig. 5-3); the
    // periodic stations on the side lines, which end at the same two points; Eq. 5-2 across
    const Cut owned = given ? Cut{} : row.cut(m);
    const Cut& cut = given ? *given : owned;
    const int n = static_cast<int>(m.size());
    if (!(1 <= first && last - first >= 2 && last <= n - 2))
        throw GridError("touching stations " + std::to_string(first) + " and " + std::to_string(last) + " of " +
                        std::to_string(n) + " stations: they need a blade station between them and periodic stations on "
                        "both sides (n_passage >= 2, an extension of at least one cell)");
    for (int i = first + 1; i < last; ++i)
        if (!cut.cuts_blade(sz(i)))
            throw GridError("station " + std::to_string(i) + " (m = " + fmt_g(m[sz(i)]) +
                            ") between the touching stations does not cut the blade");
    const double ty = b.periodic_translation.y;
    std::vector<double> lower(m.size()), upper(m.size());
    for (int i = 0; i < n; ++i) {
        if (i < first)
            lower[sz(i)] = inlet_side(b, m[sz(i)]);
        else if (i > last)
            lower[sz(i)] = exit_side(b, m[sz(i)]);
        else if (i == first)
            lower[sz(i)] = row.tip_point.y;
        else if (i == last)
            lower[sz(i)] = row.tail_point.y;
        else
            lower[sz(i)] = cut.height[sz(i)][sz(cut.count[sz(i)] - 1)];
        upper[sz(i)] = lower[sz(i)] + ty;
        if (i > first && i < last) upper[sz(i)] = cut.height[sz(i)][0] + ty;
    }
    Grid g;
    g.boundaries = b;
    g.m = m;
    g.eta = linspace(0.0, 1.0, n_eta);
    g.node_m = Field(n, n_eta);
    g.node_y = Field(n, n_eta);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n_eta; ++j) {
            const double eta = g.eta[sz(j)];
            g.node_m(i, j) = m[sz(i)];
            g.node_y(i, j) = (1.0 - eta) * lower[sz(i)] + eta * upper[sz(i)];
        }
    g.first_blade_i = first;
    g.last_blade_i = last;
    g.edge_nodes = EdgeNodes::Touching;
    finish_grid(g, row, cut);
    return g;
}

TouchingLattice touching_lattice(const CurveRow& row, double m_le, double m_te, double chord, int n_cells,
                                 double upstream_chords, double downstream_chords)
{
    // the formulas of v1 _candidate with m_first = m_tip, m_last = m_tail (extensions measured in
    // chords from the LE / TE point, rounded to whole cells)
    TouchingLattice t;
    t.dm = (row.m_tail - row.m_tip) / static_cast<double>(n_cells);
    t.k_up = py_round((row.m_tip - (m_le - upstream_chords * chord)) / t.dm);
    t.k_down = py_round(((m_te + downstream_chords * chord) - row.m_tail) / t.dm);
    t.upstream_chords = (m_le - (row.m_tip - static_cast<double>(t.k_up) * t.dm)) / chord;
    t.downstream_chords = ((row.m_tail + static_cast<double>(t.k_down) * t.dm) - m_te) / chord;
    return t;
}

void check_passage(const Cascade& cascade, const PassageBoundaries& b, const CurveRow& row)
{
    // v1 _check_passage (section 4 on the curves, section 7)
    const double pitch = b.periodic_translation.y;
    const std::string blade = (cascade.designation.empty() ? std::string("the blade") : cascade.designation) +
                              " (trailing edge: not given) at stagger " + fmt_g(cascade.stagger * (180.0 / kPi)) +
                              " deg and solidity " + fmt_g(b.chord / pitch, 10);
    const auto widest = row.widest_cut();
    if (!(widest[0] < pitch))
        throw GridError("the cut of the curves (at a passage-check piece end) of " + blade + " spans " +
                        fmt_g(widest[0] / b.chord) + " chords at m = " + fmt_g(widest[1]) +
                        ", at least the pitch: the blades overlap there and the passage closes");
    const auto narrow = row.narrowest(pitch);
    if (!(narrow[0] > row.touch))
        throw GridError("the smooth walls close the passage of " + blade + ": their narrowest width over the blade is " +
                        detail::fmt_e(narrow[0] / b.chord, 6) + " chords, at m = " + fmt_g(narrow[1]) + " (m - m_le = " +
                        fmt_g((narrow[1] - b.m_le) / b.chord) + " chords), not above TOUCH = 1e-12 chords. The check "
                        "covers the whole blade, not the stations, so no n_passage or n_eta builds this cascade; the "
                        "straight-segment outlines that geometry.py checks for overlap can still leave a gap");
}

}  // namespace detail

// ── public API ───────────────────────────────────────────────────────────────

PassageBoundaries make_passage_boundaries(const Cascade& cascade, const GridSettings& settings, double inlet_angle,
                                          std::optional<double> exit_angle)
{
    // v1 make_passage_boundaries (section 1; Aungier 2003 Sec. 5.3, p. 86)
    check_angle(inlet_angle, "inlet_angle");
    double exit = cascade.kappa_te;
    if (exit_angle) {
        check_angle(*exit_angle, "exit_angle");
        exit = *exit_angle;
    }
    check_blade_angle(cascade.kappa_le, "kappa_le", "inlet");
    check_blade_angle(cascade.kappa_te, "kappa_te", "exit");
    const double upstream = positive_number(settings.upstream_extension_chords, "upstream_extension_chords");
    const double downstream = positive_number(settings.downstream_extension_chords, "downstream_extension_chords");
    const int n_passage = whole_number(settings.n_passage, "n_passage", 1);
    if (!cascade.active_curves) throw std::invalid_argument("make_passage_boundaries: the cascade has no C2 curves");

    const CascadeCoordinates& a = cascade.active;
    const Point le = a.leading_edge_point, te = a.trailing_edge_point;
    const double m_le = le.m, m_te = te.m;
    if (!(m_te > m_le))
        throw std::invalid_argument("the TE (m = " + repr(m_te) + ") must lie downstream of the LE (m = " + repr(m_le) + ")");
    PassageBoundaries b;
    b.chord = a.chord;
    b.pitch = a.pitch;
    b.periodic_translation = a.periodic_translation;
    b.inlet_angle = inlet_angle;
    b.exit_angle = exit;
    b.m_in = m_le - upstream * a.chord;
    b.m_le = m_le;
    b.m_te = m_te;
    b.m_out = m_te + downstream * a.chord;
    b.kappa_le = cascade.kappa_le;
    b.kappa_te = cascade.kappa_te;
    b.lower_blade_boundary = a.suction_surface;
    b.upper_blade_boundary = translated(a.pressure_surface, a.periodic_translation);
    b.blade_curves = cascade.active_curves;
    // the side lines' junctions: v1 the LE / TE point; Touching (not in v1) the nose tip / tail end
    // (the side lines end at the touching nodes, Aungier 2003 Fig. 5-5), with the same blade angles
    Point j_in = le, j_out = te;
    if (settings.edge_nodes == EdgeNodes::Touching) {
        const CurveRow row(b);
        j_in = row.tip_point;
        j_out = row.tail_point;
        if (!(b.m_in < j_in.m && j_out.m < b.m_out))
            throw std::invalid_argument("the extensions end inside the blade's m range (m_in " + repr(b.m_in) +
                                        ", nose tip " + repr(j_in.m) + "; tail end " + repr(j_out.m) + ", m_out " +
                                        repr(b.m_out) + "): the touching side lines need m_in < m_tip and m_tail < m_out");
    }
    b.upstream_m = samples(b.m_in, j_in.m, n_passage, m_te - m_le);
    side_samples(b.upstream_m, true, j_in, cascade.kappa_le, inlet_angle, b.upstream_beta, b.upstream_lower_boundary);
    b.upstream_upper_boundary = translated(b.upstream_lower_boundary, a.periodic_translation);
    b.downstream_m = samples(j_out.m, b.m_out, n_passage, m_te - m_le);
    side_samples(b.downstream_m, false, j_out, cascade.kappa_te, exit, b.downstream_beta, b.downstream_lower_boundary);
    b.downstream_upper_boundary = translated(b.downstream_lower_boundary, a.periodic_translation);
    return b;
}

Grid make_grid(const Cascade& cascade, const GridSettings& settings, double inlet_angle, std::optional<double> exit_angle,
               bool require_ss_ps)
{
    // v1 make_grid
    if (settings.level) {
        if (settings.n_passage || settings.n_eta)
            throw std::invalid_argument("GridSettings with level = '" + *settings.level +
                                        "' take no n_passage or n_eta: the level sets the resolution; give either a "
                                        "level or n_passage and n_eta (the explicit mode)");
        return build_grid(cascade, *settings.level, inlet_angle, exit_angle,
                          {settings.upstream_extension_chords, settings.downstream_extension_chords}, settings.edge_nodes)
            .grid;
    }
    const PassageBoundaries b = make_passage_boundaries(cascade, settings, inlet_angle, exit_angle);
    const int n_eta = whole_number(settings.n_eta, "n_eta", 2);
    const CurveRow row(b);
    detail::check_passage(cascade, b, row);
    if (settings.edge_nodes == EdgeNodes::Touching) {
        // not in v1: n_passage cells between the nose tip and the tail end, the extensions in
        // whole cells; refused (GridTopologyError) unless the stations next to the touching ones
        // each cut one suction and one pressure surface, whatever require_ss_ps says (the layout
        // has no meaning otherwise)
        const int n = *settings.n_passage;
        if (n < 2)
            throw GridError("n_passage = " + std::to_string(n) + " in the touching layout leaves no blade station between "
                            "the nose tip and the tail end: n_passage >= 2");
        const detail::TouchingLattice lattice =
            detail::touching_lattice(row, b.m_le, b.m_te, b.chord, n, settings.upstream_extension_chords,
                                     settings.downstream_extension_chords);
        if (lattice.k_up < 1 || lattice.k_down < 1)
            throw GridError("the blade reaches a far station: lengthen the extensions");
        if (!touching_admissible(b, row, n)) {
            // the refusal with its message, before any grid is built (the lattice's stations)
            std::vector<double> m(static_cast<std::size_t>(lattice.k_up + n + lattice.k_down + 1));
            for (std::size_t k = 0; k < m.size(); ++k)
                m[k] = row.m_tip + (static_cast<double>(k) - static_cast<double>(lattice.k_up)) * lattice.dm;
            validate_touching(b, row, m, static_cast<int>(lattice.k_up), static_cast<int>(lattice.k_up) + n,
                              settings.n_passage, cascade.stagger);
        }
        GridSettings actual = settings;
        actual.upstream_extension_chords = lattice.upstream_chords;
        actual.downstream_extension_chords = lattice.downstream_chords;
        const PassageBoundaries ba = make_passage_boundaries(cascade, actual, inlet_angle, exit_angle);
        const int stations = static_cast<int>(lattice.k_up + lattice.k_down) + n + 1;
        Grid g = detail::grid_on_touching_stations(ba, row, detail::linspace(ba.m_in, ba.m_out, stations), n_eta,
                                                   static_cast<int>(lattice.k_up), static_cast<int>(lattice.k_up) + n);
        validate_uniform_grid_topology(g, settings.n_passage, cascade.stagger);
        return g;
    }
    Grid g = detail::grid_on_stations(b, row, detail::stations(b, *settings.n_passage), n_eta);
    if (require_ss_ps) validate_uniform_grid_topology(g, settings.n_passage, cascade.stagger);
    return g;
}

PassageBoundaries rebuild_exit(const PassageBoundaries& item, double exit_angle)
{
    // v1 rebuild_exit (PassageBoundaries)
    check_angle(exit_angle, "exit_angle");
    PassageBoundaries b = item;
    std::vector<double> beta;
    std::vector<Point> lower;
    side_samples(item.downstream_m, false, item.downstream_lower_boundary.front(), item.downstream_beta.front(), exit_angle,
                 beta, lower);
    b.exit_angle = exit_angle;
    b.downstream_beta = std::move(beta);
    b.downstream_upper_boundary = translated(lower, item.periodic_translation);
    b.downstream_lower_boundary = std::move(lower);
    return b;
}

Grid rebuild_exit(const Grid& item, double exit_angle)
{
    // v1 rebuild_exit (Grid): only the exit stations' nodes move; the Grid is checked again
    Grid g = item;
    g.boundaries = rebuild_exit(item.boundaries, exit_angle);
    const double ty = g.boundaries.periodic_translation.y;
    for (int i = item.last_blade_i + 1; i < g.n_m(); ++i) {
        const double lower = exit_side(g.boundaries, g.m[sz(i)]);
        const double upper = lower + ty;
        for (int j = 0; j < g.n_eta(); ++j) {
            const double eta = g.eta[sz(j)];
            g.node_y(i, j) = (1.0 - eta) * lower + eta * upper;
        }
    }
    detail::finish_grid(g, CurveRow(g.boundaries));
    return g;
}

std::array<EdgeStationTopology, 2> wall_station_topology(const Grid& grid)
{
    const CurveRow row(grid.boundaries);
    const Cut cut = row.cut({grid.m[sz(grid.first_blade_i)], grid.m[sz(grid.last_blade_i)]});
    return {edge_topology(row, grid.boundaries, cut, 0, "LE", grid.first_blade_i),
            edge_topology(row, grid.boundaries, cut, 1, "TE", grid.last_blade_i)};
}

void validate_uniform_grid_topology(const Grid& grid) { validate_uniform_grid_topology(grid, std::nullopt, std::nullopt); }

void validate_uniform_grid_topology(const Grid& grid, std::optional<int> n_passage, std::optional<double> stagger)
{
    // v1 validate_uniform_grid_topology (section 8; Aungier 2003 pp. 89-90, Fig. 5-5): rejects, never repairs
    if (grid.touching()) {
        validate_touching(grid.boundaries, CurveRow(grid.boundaries), grid.m, grid.first_blade_i, grid.last_blade_i,
                          n_passage, stagger);
        return;
    }
    const auto topology = wall_station_topology(grid);
    const EdgeStationTopology& le = topology[0];
    const EdgeStationTopology& te = topology[1];
    if (le.valid() && te.valid()) return;
    const PassageBoundaries& b = grid.boundaries;
    const CurveRow row(b);
    const auto ceiling = ceiling_of(b, row);
    auto admissible = [&](int n) {
        const auto t = lattice_topology(b, row, n);
        return t && (*t)[0].valid() && (*t)[1].valid();
    };
    std::optional<int> finest, lower, higher;
    if (ceiling)
        for (int n = *ceiling - 1; n > 0; --n)
            if (admissible(n)) {
                finest = n;
                break;
            }
    if (n_passage) {
        const int np = whole_number(n_passage, "n_passage", 1);
        if (ceiling && np >= *ceiling) {
            lower = finest;
        } else {
            for (int n = np - 1; n > 0; --n)
                if (admissible(n)) {
                    lower = n;
                    break;
                }
            const int top = ceiling ? *ceiling : 4 * np + 100;
            for (int n = np + 1; n < top; ++n)
                if (admissible(n)) {
                    higher = n;
                    break;
                }
        }
    }
    std::vector<const EdgeStationTopology*> bad;
    for (const auto* t : {&le, &te})
        if (!t->valid()) bad.push_back(t);
    std::string edge;
    for (const auto* t : bad) edge += (edge.empty() ? "" : " and ") + t->edge;
    const std::string named = stagger ? " (stagger " + fmt_g(*stagger * (180.0 / kPi)) + " deg)" : "";
    const std::string head = n_passage ? "n_passage = " + std::to_string(*n_passage) + " is not admissible for this cascade" + named
                                       : "this uniform grid is not admissible for its cascade" + named;
    std::string text = head + ": the first and the last blade station must each cut one suction and one pressure surface.";
    for (const auto* t : bad) text += "\n" + describe_edge(*t);
    auto opt = [](const std::optional<int>& v) { return v ? std::to_string(*v) : std::string("none"); };
    std::string choice = "Choose another admissible n_passage; none is substituted.";
    if (n_passage) choice += " Nearest admissible: lower " + opt(lower) + ", higher " + opt(higher) + ".";
    if (ceiling) {
        const auto ext = row.extremes();
        // the widest window above TOUCH (v1: max over the windows dict, first one on a tie)
        const double w_le = b.m_le - ext[0], w_te = ext[1] - b.m_te;
        std::string name = "LE window m_le - nose tip";
        double width = w_le;
        if (!(w_le > row.touch) || (w_te > row.touch && w_te > w_le)) {
            name = "TE window tail end - m_te";
            width = w_te;
        }
        choice += " Finest admissible n_passage: " + opt(finest) + "; every n_passage from " + std::to_string(*ceiling) +
                  " on fails (dm <= " + fmt_g(width, 4) + ", the " + name +
                  "). Coarser values can fail too: which ones pass depends on where the uniform stations fall "
                  "(find_valid_n_passage lists them).";
    } else {
        choice += " This blade has no LE or TE window, so no refinement ceiling exists.";
    }
    text += "\n" + choice;
    throw GridTopologyError(text, edge, n_passage, lower, higher, finest, ceiling);
}

std::vector<int> find_valid_n_passage(const Cascade& cascade, const GridSettings& settings, double inlet_angle,
                                      const std::vector<int>& n_values)
{
    const PassageBoundaries b = boundaries_for(cascade, settings, inlet_angle);
    const CurveRow row(b);
    detail::check_passage(cascade, b, row);
    std::vector<int> out;
    const bool touching = settings.edge_nodes == EdgeNodes::Touching;
    for (int value : n_values) {
        const int n = whole_number(value, "n_passage", 1);
        if (touching) {
            if (n < 2) continue;
            const auto lattice = detail::touching_lattice(row, b.m_le, b.m_te, b.chord, n, settings.upstream_extension_chords,
                                                          settings.downstream_extension_chords);
            if (lattice.k_up >= 1 && lattice.k_down >= 1 && touching_admissible(b, row, n)) out.push_back(n);
            continue;
        }
        const auto t = lattice_topology(b, row, n);
        if (t && (*t)[0].valid() && (*t)[1].valid()) out.push_back(n);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::optional<int> n_passage_ceiling(const Cascade& cascade, const GridSettings& settings, double inlet_angle)
{
    const PassageBoundaries b = boundaries_for(cascade, settings, inlet_angle);
    if (settings.edge_nodes == EdgeNodes::Touching) return touching_ceiling(b, CurveRow(b));
    return ceiling_of(b, CurveRow(b));
}

}  // namespace BladeGen::Turbolenza
