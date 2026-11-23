// The geometry-driven standard grids: port of v1 grid.py section 9 (passage_transitions,
// build_grid, the candidate search with its rejection reasons, validate_standard_grid,
// StandardGridInfo).  Constant dm and deta (Aungier 2003 Sec. 5.3, p. 87); the lattice puts the
// cell faces at the middles of the LE and TE windows, so the first and the last blade station
// each cut one suction and one pressure surface (Aungier 2003 pp. 89-90, Fig. 5-5).
//
// Not ported: v1's timing (time.perf_counter), geometry_id (a sha256 fingerprint used in the
// messages; the cascade designation names the geometry here), build_grid_family (the window
// builds the levels it needs one by one), and the last check of validate_standard_grid,
// "metrics with positive S" (metrics.py; in v2 the metrics come after the grid).
#include "turbolenza/Grid.hpp"

#include "turbolenza/GridDetail.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <stdexcept>

namespace BladeGen::Turbolenza {

namespace {

using detail::CurveRow;
using detail::Cut;
using detail::fmt_e;
using detail::fmt_g;

std::size_t sz(int k) { return static_cast<std::size_t>(k); }

std::string fmt_f(double x, int precision)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", precision, x);
    return buf;
}

std::string fmt_signed_e(double x, int precision)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%+.*e", precision, x);
    return buf;
}

/// v1 _Rejected: a candidate fails one condition.
struct Rejected
{
    std::string reason, detail;
};

struct Interrogation
{
    const Cascade* cascade = nullptr;
    PassageBoundaries template_b;
    std::shared_ptr<CurveRow> row;
    PassageTransitions transitions;
    std::string geometry;  ///< names the geometry in the messages (v1: geometry_id)
};

std::string geometry_name(const Cascade& cascade)
{
    return "'" + (cascade.designation.empty() ? std::string("the blade") : cascade.designation) + "'";
}

std::string level_name(const std::string& level)
{
    for (const char* name : kGridLevels)
        if (level == name) return level;
    throw std::invalid_argument("level must be one of very_coarse, coarse, medium, fine, very_fine, got '" + level + "'");
}

double positive(double value, const char* name)
{
    if (!(std::isfinite(value) && value > 0.0))
        throw std::invalid_argument(std::string(name) + " must be a positive number of chords, got " + fmt_g(value, 17));
    return value;
}

/// v1 _transitions.
PassageTransitions transitions_of(const PassageBoundaries& b, const CurveRow& row, const std::string& geometry)
{
    std::vector<double> ends = row.lo_m;
    ends.insert(ends.end(), row.hi_m.begin(), row.hi_m.end());
    ends = detail::unique(std::move(ends));
    PassageTransitions tr;
    std::vector<bool> valid;
    for (std::size_t i = 0; i + 1 < ends.size(); ++i) {
        const double middle = 0.5 * (ends[i] + ends[i + 1]);
        int ss = 0, ps = 0;
        for (int k = 0; k < row.n_branches(); ++k)
            if (row.lo_m[sz(k)] <= middle && middle < row.hi_m[sz(k)]) (row.curve_of[sz(k)] == 0 ? ss : ps) += 1;
        tr.coverage.push_back({ends[i], ends[i + 1], ss, ps});
        valid.push_back(ss == 1 && ps == 1);
    }
    const double centre = 0.5 * (b.m_le + b.m_te);
    const int k = static_cast<int>(std::upper_bound(ends.begin(), ends.end(), centre) - ends.begin()) - 1;
    if (!(0 <= k && k < static_cast<int>(valid.size()) && valid[sz(k)])) {
        std::string found = "None";
        if (0 <= k && k < static_cast<int>(tr.coverage.size())) {
            const auto& c = tr.coverage[sz(k)];
            found = "(" + fmt_g(c.lo, 17) + ", " + fmt_g(c.hi, 17) + ", " + std::to_string(c.suction) + ", " +
                    std::to_string(c.pressure) + ")";
        }
        throw GridGenerationError("geometry " + geometry +
                                      " could not produce an admissible grid at any level: the line through the middle of "
                                      "the chord does not cross the blade once on the suction and once on the pressure "
                                      "surface (" + found + "): no SS + PS-valid passage interval",
                                  std::nullopt, {}, {{kRejectBranch, 1}});
    }
    int first = k, last = k;
    while (first > 0 && valid[sz(first - 1)]) --first;
    while (last < static_cast<int>(valid.size()) - 1 && valid[sz(last + 1)]) ++last;
    const auto ext = row.extremes();
    tr.m_tip = ext[0];
    tr.m_le = b.m_le;
    tr.m_start = ends[sz(first)];
    tr.m_end = ends[sz(last + 1)];
    tr.m_te = b.m_te;
    tr.m_tail = ext[1];
    tr.chord = b.chord;
    return tr;
}

/// v1 _interrogate: what every level shares.
Interrogation interrogate(const Cascade& cascade)
{
    Interrogation it;
    it.cascade = &cascade;
    it.geometry = geometry_name(cascade);
    GridSettings s;
    s.n_passage = 1;
    s.n_eta = 2;
    it.template_b = make_passage_boundaries(cascade, s, 0.0, 0.0);
    it.row = std::make_shared<CurveRow>(it.template_b);
    it.transitions = transitions_of(it.template_b, *it.row, it.geometry);
    try {
        detail::check_passage(cascade, it.template_b, *it.row);
    } catch (const GridError& error) {
        throw GridGenerationError("geometry " + it.geometry + " could not produce an admissible grid at any level: " +
                                      error.what(),
                                  std::nullopt, {}, {{kRejectOther, 1}});
    }
    return it;
}

/// Suction and pressure crossings of station i of the cut.
std::array<int, 2> crossings(const CurveRow& row, const Cut& cut, std::size_t i)
{
    std::array<int, 2> out{0, 0};
    for (int j = 0; j < cut.count[i]; ++j) out[row.surface_of(cut.segment[i][sz(j)]) == Surface::Suction ? 0 : 1] += 1;
    return out;
}

/// v1 _classify.
void classify(const CurveRow& row, const PassageTransitions& tr, const Cut& cut, int first, int last)
{
    for (int i = first - 1; i >= 0; --i)
        if (cut.count[sz(i)] > 0)
            throw Rejected{kRejectLE, "station " + std::to_string(i) + " ahead of the first SS + PS station meets the blade (" +
                                          std::to_string(cut.count[sz(i)]) + " crossings) at (m - m_start)/c = " +
                                          fmt_signed_e((cut.m[sz(i)] - tr.m_start) / tr.chord, 3) + ": the LE window is " +
                                          fmt_e(tr.le_window() / tr.chord, 3) + " c"};
    for (int i = last + 1; i < static_cast<int>(cut.m.size()); ++i)
        if (cut.count[sz(i)] > 0)
            throw Rejected{kRejectTE, "station " + std::to_string(i) + " behind the last SS + PS station meets the blade (" +
                                          std::to_string(cut.count[sz(i)]) + " crossings) at (m - m_end)/c = " +
                                          fmt_signed_e((cut.m[sz(i)] - tr.m_end) / tr.chord, 3) + ": the TE window is " +
                                          fmt_e(tr.te_window() / tr.chord, 3) + " c"};
    for (int i = first; i <= last; ++i) {
        const auto c = crossings(row, cut, sz(i));
        if (!(cut.count[sz(i)] == 2 && c[0] == 1 && c[1] == 1))
            throw Rejected{kRejectBranch, "blade station " + std::to_string(i) + " crosses the blade " +
                                              std::to_string(cut.count[sz(i)]) + " times (" + std::to_string(c[0]) +
                                              " suction, " + std::to_string(c[1]) + " pressure)"};
    }
    for (int i = first; i <= last; ++i) {
        const double m = cut.m[sz(i)];
        if (std::abs(m - tr.m_le) <= row.touch || std::abs(m - tr.m_te) <= row.touch)
            throw Rejected{kRejectBranch, "blade station " + std::to_string(i) + " passes through the LE or TE point"};
    }
}

/// v1 _standard_checks (without "metrics with positive S", see the file comment).
std::vector<std::string> standard_checks(const Grid& grid, const CurveRow& row, const PassageTransitions& tr, const Cut& cut)
{
    const double chord = grid.chord();
    const auto& m = grid.m;
    const auto& eta = grid.eta;
    const int f = grid.first_blade_i, l = grid.last_blade_i, n = grid.n_m(), ne = grid.n_eta();
    std::vector<std::string> checks;
    auto need = [&](bool ok, const char* name, const char* reason, const std::string& detail) {
        if (!ok) throw Rejected{reason, std::string(name) + ": " + detail};
        checks.emplace_back(name);
    };
    bool increasing = true;
    double spread = 0.0, eta_spread = 0.0;
    const double dm = grid.dm(), deta = grid.deta();
    for (int i = 0; i + 1 < n; ++i) {
        const double step = m[sz(i + 1)] - m[sz(i)];
        increasing = increasing && step > 0.0;
        spread = std::max(spread, std::abs(step - dm));
    }
    for (int j = 0; j + 1 < ne; ++j) eta_spread = std::max(eta_spread, std::abs((eta[sz(j + 1)] - eta[sz(j)]) - deta));
    need(increasing, "distinct stations", kRejectOther, "two stations coincide or run backward");
    need(spread <= kUniformSpacingTolerance * chord, "uniform dm", kRejectOther,
         "spacing varies by " + fmt_g(spread / chord, 3) + " c");
    need(eta.front() == 0.0 && eta.back() == 1.0 && eta_spread <= kUniformSpacingTolerance, "uniform deta", kRejectOther,
         "eta spacing varies by " + fmt_g(eta_spread, 3));
    need(1 <= f && f <= l && l <= n - 2, "f <= l inside the far stations", kRejectOther,
         "f = " + std::to_string(f) + ", l = " + std::to_string(l) + ", n_m = " + std::to_string(n));
    need(!tr.valid_cut(m[sz(f - 1)]) && cut.count[sz(f - 1)] == 0, "f - 1 outside the blade", kRejectLE,
         "station " + std::to_string(f - 1) + " meets the blade");
    need(!tr.valid_cut(m[sz(l + 1)]) && cut.count[sz(l + 1)] == 0, "l + 1 outside the blade", kRejectTE,
         "station " + std::to_string(l + 1) + " meets the blade");
    bool inside = true, one_each = true;
    for (int i = f; i <= l; ++i) {
        inside = inside && tr.valid_cut(m[sz(i)]);
        const auto c = crossings(row, cut, sz(i));
        one_each = one_each && cut.count[sz(i)] == 2 && c[0] == 1 && c[1] == 1;
    }
    need(inside, "every blade station inside the SS + PS-valid interval", kRejectBranch, "a blade station lies outside it");
    need(one_each, "one SS + one PS crossing per blade station", kRejectBranch, "a blade station has other crossings");
    const auto& id = grid.wall_identity;
    const std::set<Surface> lower(id.lower.surface.begin(), id.lower.surface.end());
    const std::set<Surface> upper(id.upper.surface.begin(), id.upper.surface.end());
    auto names = [](const std::set<Surface>& s) {
        std::string out = "[";
        for (Surface x : s) out += std::string(out.size() > 1 ? ", " : "") + (x == Surface::Suction ? "'suction'" : "'pressure'");
        return out + "]";
    };
    need(lower.size() == 1 && upper.size() == 1 && lower != upper, "wall identity constant", kRejectBranch,
         "lower wall on " + names(lower) + ", upper wall on " + names(upper));
    bool hidden = false;
    for (int h : id.lower.hidden) hidden = hidden || h != 0;
    for (int h : id.upper.hidden) hidden = hidden || h != 0;
    need(!hidden, "no hidden crossings", kRejectBranch, "a wall node hides further crossings of its contour");
    bool finite = true;
    for (std::size_t k = 0; k < grid.node_m.size(); ++k)
        finite = finite && std::isfinite(grid.node_m.v[k]) && std::isfinite(grid.node_y.v[k]);
    need(finite, "finite nodes", kRejectOther, "a node is not finite");
    bool positive_width = true;
    double min_width = 0.0;
    for (int i = 0; i < n; ++i) {
        const double w = grid.node_y(i, ne - 1) - grid.node_y(i, 0);
        positive_width = positive_width && w > row.touch;
        min_width = i == 0 ? w : std::min(min_width, w);
    }
    need(positive_width, "positive width", kRejectOther, "width " + fmt_g(min_width / chord, 3) + " c");
    return checks;
}

/// v1 _info.
StandardGridInfo info_of(const Interrogation& it, const std::string& level, const Grid& grid, int n_cells, double target,
                         std::vector<std::string> checks)
{
    const PassageTransitions& tr = it.transitions;
    const double chord = tr.chord;
    const PassageBoundaries& b = grid.boundaries;
    const int f = grid.first_blade_i, l = grid.last_blade_i, ne = grid.n_eta();
    const Point shift = b.periodic_translation;
    const Point le = b.lower_blade_boundary.front(), te = b.lower_blade_boundary.back();
    const Point axis{te.m - le.m, te.y - le.y};
    auto points = [&](int i) {
        const Point lower{grid.node_m(i, 0), grid.node_y(i, 0)};
        const Point upper{grid.node_m(i, ne - 1) - shift.m, grid.node_y(i, ne - 1) - shift.y};
        const bool suction = grid.wall_identity.lower.surface[sz(i - f)] == Surface::Suction;
        return std::array<Point, 2>{suction ? lower : upper, suction ? upper : lower};
    };
    auto x_c = [&](const Point& p) {
        return ((p.m - le.m) * axis.m + (p.y - le.y) * axis.y) / (axis.m * axis.m + axis.y * axis.y);
    };
    const auto first = points(f), last = points(l);
    StandardGridInfo info;
    info.level = level;
    info.target_dm_over_chord = target / chord;
    info.n_blade_cells = n_cells;
    info.dm = grid.dm();
    info.dm_over_chord = grid.dm() / chord;
    info.deta = grid.deta();
    info.n_m = grid.n_m();
    info.n_eta = grid.n_eta();
    info.total_nodes = grid.n_m() * grid.n_eta();
    info.first_blade_i = f;
    info.last_blade_i = l;
    info.n_blade_stations = l - f + 1;
    info.m_first = grid.m[sz(f)];
    info.m_last = grid.m[sz(l)];
    info.first_ss_wall_point = first[0];
    info.first_ps_wall_point = first[1];
    info.last_ss_wall_point = last[0];
    info.last_ps_wall_point = last[1];
    info.first_ss_x_c = x_c(first[0]);
    info.first_ps_x_c = x_c(first[1]);
    info.last_ss_x_c = x_c(last[0]);
    info.last_ps_x_c = x_c(last[1]);
    info.le_gap_over_chord = (grid.m[sz(f)] - tr.m_start) / chord;
    info.te_gap_over_chord = (tr.m_end - grid.m[sz(l)]) / chord;
    info.m_first_minus_m_le_over_chord = (grid.m[sz(f)] - tr.m_le) / chord;
    info.m_te_minus_m_last_over_chord = (tr.m_te - grid.m[sz(l)]) / chord;
    info.upstream_extension_chords = (tr.m_le - grid.m.front()) / chord;
    info.downstream_extension_chords = (grid.m.back() - tr.m_te) / chord;
    info.cell_aspect = b.pitch * grid.deta() / grid.dm();
    info.admissibility = "admissible";
    info.checks = std::move(checks);
    return info;
}

/// v1 _candidate: one candidate spacing; Rejected unless every condition holds.
StandardGrid candidate(const Interrogation& it, const std::string& level, int n_cells, double target, double inlet_angle,
                       std::optional<double> exit_angle, double up, double down)
{
    const PassageTransitions& tr = it.transitions;
    const CurveRow& row = *it.row;
    const double chord = tr.chord;
    const double dm = tr.lattice_length() / static_cast<double>(n_cells);
    const double m_first = tr.le_face() + 0.5 * dm, m_last = tr.te_face() - 0.5 * dm;
    for (int e = 0; e < 2; ++e) {
        const char* reason = e == 0 ? kRejectLE : kRejectTE;
        const char* edge = e == 0 ? "LE" : "TE";
        const double window = e == 0 ? tr.le_window() : tr.te_window();
        const double clearance = 0.5 * (dm - window);  // of the transition stations from the window, on both sides
        if (!(clearance >= kTransitionClearance * dm))
            throw Rejected{reason, std::string("the ") + edge + " window " + fmt_e(window / chord, 3) +
                                       " c leaves the stations beside it " + fmt_e(clearance / chord, 3) +
                                       " c of clearance, below " + fmt_g(kTransitionClearance) + " dm (" +
                                       fmt_e(kTransitionClearance * dm / chord, 3) + " c): dm must exceed the window by " +
                                       fmt_g(2 * kTransitionClearance) + " dm"};
    }
    const long long k_up = detail::py_round((m_first - (tr.m_le - up * chord)) / dm);
    const long long k_down = detail::py_round(((tr.m_te + down * chord) - m_last) / dm);
    const double upstream = (tr.m_le - (m_first - static_cast<double>(k_up) * dm)) / chord;
    const double downstream = ((m_last + static_cast<double>(k_down) * dm) - tr.m_te) / chord;
    for (int e = 0; e < 2; ++e) {
        const double actual = e == 0 ? upstream : downstream, nominal = e == 0 ? up : down;
        if (!(std::abs(actual - nominal) <= kExtensionToleranceChords))
            throw Rejected{kRejectExtension, std::string(e == 0 ? "upstream" : "downstream") + " extension " +
                                                 fmt_f(actual, 4) + " c, nominal " + fmt_g(nominal) + " c, tolerance " +
                                                 fmt_g(kExtensionToleranceChords) + " c"};
    }
    if (k_up < 1 || k_down < 1) throw Rejected{kRejectExtension, "an extension holds no whole cell"};
    const long long eta_cells = detail::py_round(it.template_b.pitch / (kCellAspect * dm));
    if (eta_cells < kMinEtaCells)
        throw Rejected{kRejectEta, std::to_string(eta_cells) + " eta cells (pitch / (" + fmt_g(kCellAspect) +
                                       " dm)), fewer than " + std::to_string(kMinEtaCells)};
    const int first = static_cast<int>(k_up), last = static_cast<int>(k_up) + n_cells - 1;

    GridSettings settings;
    settings.n_passage = n_cells;
    settings.n_eta = static_cast<int>(eta_cells) + 1;
    settings.upstream_extension_chords = upstream;
    settings.downstream_extension_chords = downstream;
    const PassageBoundaries b = make_passage_boundaries(*it.cascade, settings, inlet_angle, exit_angle);
    const auto m = detail::linspace(b.m_in, b.m_out, static_cast<int>(k_up + n_cells + k_down));
    const Cut cut = row.cut(m);
    classify(row, tr, cut, first, last);
    StandardGrid out;
    try {
        out.grid = detail::grid_on_stations(b, row, m, static_cast<int>(eta_cells) + 1, &cut);
    } catch (const GridError& error) {
        std::string what = error.what();
        throw Rejected{kRejectOther, what.substr(0, what.find('\n'))};
    }
    if (out.grid.first_blade_i != first || out.grid.last_blade_i != last)
        throw Rejected{kRejectOther, "blade run " + std::to_string(out.grid.first_blade_i) + " .. " +
                                         std::to_string(out.grid.last_blade_i) + ", expected " + std::to_string(first) +
                                         " .. " + std::to_string(last)};
    auto checks = standard_checks(out.grid, row, tr, cut);
    out.info = info_of(it, level, out.grid, n_cells, target, std::move(checks));
    return out;
}

/// v1 _build_level: the candidate search of one level.
StandardGrid build_level(const Interrogation& it, const std::string& level, double inlet_angle,
                         std::optional<double> exit_angle, std::array<double, 2> extension_chords)
{
    const PassageTransitions& tr = it.transitions;
    const double up = positive(extension_chords[0], "upstream extension_chords");
    const double down = positive(extension_chords[1], "downstream extension_chords");
    const double target = level_target_dm_over_chord(level) * tr.chord;
    const double low = target * (1.0 - kCandidateBand), high = target * (1.0 + kCandidateBand);
    const double length = tr.lattice_length();
    std::vector<int> ranked;
    for (long long n = std::max(1LL, static_cast<long long>(std::ceil(length / high)));
         n <= static_cast<long long>(std::floor(length / low)); ++n)
        ranked.push_back(static_cast<int>(n));
    auto key = [&](int n) { return std::abs(std::log(length / static_cast<double>(n) / target)); };
    std::stable_sort(ranked.begin(), ranked.end(), [&](int a, int b) {
        const double ka = key(a), kb = key(b);
        return ka < kb || (ka == kb && a < b);
    });
    std::vector<CandidateRecord> records;
    const std::array<double, 2> band{low / tr.chord, high / tr.chord};
    int rank = 0;
    for (int n_cells : ranked) {
        ++rank;
        try {
            StandardGrid out = candidate(it, level, n_cells, target, inlet_angle, exit_angle, up, down);
            records.push_back({n_cells, out.info.dm_over_chord, rank, true, "", "accepted"});
            out.candidates = std::move(records);
            out.target_dm_over_chord = target / tr.chord;
            out.band = band;
            return out;
        } catch (const Rejected& r) {
            records.push_back({n_cells, length / static_cast<double>(n_cells) / tr.chord, rank, false, r.reason, r.detail});
        }
    }
    std::vector<std::pair<std::string, int>> reasons;
    std::string tried;
    for (const auto& r : records) {
        auto at = std::find_if(reasons.begin(), reasons.end(), [&](const auto& p) { return p.first == r.reason; });
        if (at == reasons.end())
            reasons.emplace_back(r.reason, 1);
        else
            ++at->second;
        tried += (tried.empty() ? "" : "; ") + std::to_string(r.n_cells) + " cells (dm/c " + fmt_f(r.dm_over_chord, 5) +
                 "): " + r.reason + ", " + r.detail;
    }
    if (tried.empty()) tried = "none";
    throw GridGenerationError("geometry " + it.geometry + " could not produce an admissible " + level + " grid: target dm/c " +
                                  fmt_f(target / tr.chord, 5) + ", band " + fmt_f(band[0], 5) + " .. " + fmt_f(band[1], 5) +
                                  ", " + std::to_string(records.size()) + " candidate(s) tried: " + tried,
                              level, records, reasons);
}

}  // namespace

double level_target_dm_over_chord(const std::string& level)
{
    for (std::size_t k = 0; k < kGridLevels.size(); ++k)
        if (level == kGridLevels[k])
            return kMediumDmOverChord * std::pow(std::sqrt(2.0), 2.0 - static_cast<double>(k));
    (void)level_name(level);  // throws for an unknown level
    return 0.0;
}

PassageTransitions passage_transitions(const Cascade& cascade) { return interrogate(cascade).transitions; }

StandardGrid build_grid(const Cascade& cascade, const std::string& level, double inlet_angle, std::optional<double> exit_angle)
{
    return build_grid(cascade, level, inlet_angle, exit_angle, {kNominalExtensionChords, kNominalExtensionChords});
}

StandardGrid build_grid(const Cascade& cascade, const std::string& level, double inlet_angle, std::optional<double> exit_angle,
                        std::array<double, 2> extension_chords)
{
    const std::string name = level_name(level);
    const Interrogation it = interrogate(cascade);
    return build_level(it, name, inlet_angle, exit_angle, extension_chords);
}

std::vector<std::string> validate_standard_grid(const Grid& grid)
{
    const CurveRow row(grid.boundaries);
    const PassageTransitions tr = transitions_of(grid.boundaries, row, "of this grid");
    try {
        return standard_checks(grid, row, tr, row.cut(grid.m));
    } catch (const Rejected& r) {
        throw GridError(r.reason + ": " + r.detail);
    }
}

}  // namespace BladeGen::Turbolenza
