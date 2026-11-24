// Cascade generation: port of the cascade part of v1 geometry.py (make_cascade without the NACA 65
// construction, _check_cascade_inputs, place_cascade, find_geometric_throat / closest_points and
// their helpers, signed_area, _paths_cross, and the checks of validate_cascade that apply to a C2
// cascade built from given curves; Aungier 2003 Secs. 4.9, 4.10).
//
// Ported checks of validate_cascade: _check_shapes, _check_frame, _check_surfaces,
// _check_periodicity, _check_throat, _check_angles (its kappa / camber consistency) and
// _check_curves (parameters, samples on their curves, the cascade-frame curves the blade-frame
// ones turned by the stagger).  Not ported, because they test v1's own NACA 65 construction:
// _check_thickness (against metadata max_thickness_over_c), the chi_le = chi_te = camber / 2 rule
// of the equivalent circular arc, _check_edges (the Kovach TE circle and the table LE radius),
// _check_scaling (v2 is chord-normalized: active == normalized), the quintic degree of
// _check_curves and _check_leading_edge (the curvature of v1's LE fair).
#include "turbolenza/Cascade.hpp"

#include "turbolenza/CubicSpline.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <stdexcept>

namespace BladeGen::Turbolenza {

namespace {

using detail::CubicSpline;

constexpr double kPi = 3.141592653589793;

// v1 geometry.py constants
constexpr double CANDIDATE_MARGIN = 1e-3;     // chords: local minima this close to the best are refined
constexpr int MAX_THROAT_ITERATIONS = 200;
constexpr double THROAT_PARAMETER_TOL = 1e-13;  // chords of spline parameter: refinement converged
constexpr double ARMIJO = 1e-4;               // sufficient-decrease constant of the line search
constexpr double MIN_STEP_FRACTION = 1e-12;
constexpr double NEWTON_ZONE = 1e-6;          // chords: Newton steps this small are taken without line search
constexpr double EXACT_TOL = 1e-12;           // chords: quantities that hold by construction
constexpr double THROAT_TOL = 1e-9;           // chords: suction/pressure throat against the whole-blade one

std::string g(double x, int precision = 6)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", precision, x);
    return buf;
}

std::string repr(double x) { return g(x, 17); }

std::vector<Point> rotated(const std::vector<Point>& p, double stagger)
{
    std::vector<Point> out(p.size());
    for (std::size_t k = 0; k < p.size(); ++k) out[k] = rotate_by_stagger(p[k], stagger);
    return out;
}

std::vector<Point> translated(const std::vector<Point>& p, const Point& t)
{
    std::vector<Point> out(p.size());
    for (std::size_t k = 0; k < p.size(); ++k) out[k] = {p[k].m + t.m, p[k].y + t.y};
    return out;
}

/// v1 closed_loop: TE -> suction -> LE -> pressure -> TE.
std::vector<Point> closed_loop(const std::vector<Point>& suction, const std::vector<Point>& pressure)
{
    std::vector<Point> loop(suction.rbegin(), suction.rend());
    loop.insert(loop.end(), pressure.begin() + 1, pressure.end());
    return loop;
}

double turn(const Point& o, const Point& a, const Point& b)
{
    return (a.m - o.m) * (b.y - o.y) - (a.y - o.y) * (b.m - o.m);
}

/// v1 _paths_cross: any segment of path_a properly crosses a segment of path_b.
bool paths_cross(const std::vector<Point>& a, const std::vector<Point>& b)
{
    for (std::size_t i = 0; i + 1 < a.size(); ++i)
        for (std::size_t j = 0; j + 1 < b.size(); ++j) {
            const bool straddle_b = turn(a[i], a[i + 1], b[j]) * turn(a[i], a[i + 1], b[j + 1]) < 0.0;
            if (!straddle_b) continue;
            if (turn(b[j], b[j + 1], a[i]) * turn(b[j], b[j + 1], a[i + 1]) < 0.0) return true;
        }
    return false;
}

/// v1 signed_area: shoelace area of a closed loop (last point = first).
double signed_area(const std::vector<Point>& loop)
{
    double s = 0.0;
    for (std::size_t k = 0; k + 1 < loop.size(); ++k) s += loop[k].m * loop[k + 1].y - loop[k + 1].m * loop[k].y;
    return 0.5 * s;
}

bool close(double a, double b, double tol = EXACT_TOL) { return std::abs(a - b) <= tol; }
bool close(const Point& a, const Point& b, double tol = EXACT_TOL)
{
    return std::max(std::abs(a.m - b.m), std::abs(a.y - b.y)) <= tol;
}

/// v1 _chord_length_spline: cubic spline through points, parametrized by cumulative chord length.
CubicSpline chord_length_spline(const std::vector<Point>& points)
{
    std::vector<double> x{0.0};
    for (std::size_t k = 0; k + 1 < points.size(); ++k)
        x.push_back(x.back() + std::hypot(points[k + 1].m - points[k].m, points[k + 1].y - points[k].y));
    return {std::move(x), points};
}

struct VertexDistances
{
    std::vector<double> distance, along;
    std::vector<std::size_t> segment;
};

/// v1 _vertex_distances: distance from each point to a polyline, the nearest segment, the fraction along it.
VertexDistances vertex_distances(const std::vector<Point>& points, const std::vector<Point>& polyline)
{
    VertexDistances out;
    for (const Point& p : points) {
        double best = std::numeric_limits<double>::infinity(), best_along = 0.0;
        std::size_t best_segment = 0;
        bool first = true;
        for (std::size_t j = 0; j + 1 < polyline.size(); ++j) {
            const double sm = polyline[j + 1].m - polyline[j].m, sy = polyline[j + 1].y - polyline[j].y;
            const double om = p.m - polyline[j].m, oy = p.y - polyline[j].y;
            const double along = std::clamp((om * sm + oy * sy) / (sm * sm + sy * sy), 0.0, 1.0);
            const double d = std::hypot(om - along * sm, oy - along * sy);
            if (first || d < best) {  // np.argmin: the first smallest
                best = d;
                best_segment = j;
                best_along = along;
                first = false;
            }
        }
        out.distance.push_back(best);
        out.segment.push_back(best_segment);
        out.along.push_back(best_along);
    }
    return out;
}

std::vector<std::size_t> local_minima(const std::vector<double>& v, double ceiling)
{
    std::vector<std::size_t> out;
    const double inf = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < v.size(); ++i) {
        const double left = i == 0 ? inf : v[i - 1], right = i + 1 == v.size() ? inf : v[i + 1];
        if (v[i] <= left && v[i] <= right && v[i] <= ceiling) out.push_back(i);
    }
    return out;
}

/// v1 _refine_closest: projected Newton with an Armijo line search on f(u, v) = |a(u) - b(v)|^2 / 2.
std::array<double, 2> refine_closest(const CubicSpline& A, const CubicSpline& B, std::array<double, 2> p)
{
    const std::array<double, 2> upper{A.x().back(), B.x().back()};
    auto value = [&](const std::array<double, 2>& q) {
        const Point a = A(q[0]), b = B(q[1]);
        const double dm = a.m - b.m, dy = a.y - b.y;
        return 0.5 * (dm * dm + dy * dy);
    };
    auto clip = [&](double x, int k) { return std::clamp(x, 0.0, upper[static_cast<std::size_t>(k)]); };
    auto dot = [](const Point& a, const Point& b) { return a.m * b.m + a.y * b.y; };
    double f = value(p);
    for (int it = 0; it < MAX_THROAT_ITERATIONS; ++it) {
        const Point a = A(p[0], 0), a1 = A(p[0], 1), a2 = A(p[0], 2);
        const Point b = B(p[1], 0), b1 = B(p[1], 1), b2 = B(p[1], 2);
        const Point gap{a.m - b.m, a.y - b.y};
        const std::array<double, 2> grad{dot(gap, a1), -dot(gap, b1)};
        const double h00 = dot(a1, a1) + dot(gap, a2), h01 = -dot(a1, b1), h11 = dot(b1, b1) - dot(gap, b2);
        std::array<bool, 2> free{};
        for (int k = 0; k < 2; ++k) {
            const auto K = static_cast<std::size_t>(k);
            free[K] = !((p[K] <= 0.0 && grad[K] > 0.0) || (p[K] >= upper[K] && grad[K] < 0.0));
        }
        if (!free[0] && !free[1]) return p;
        std::array<double, 2> step{0.0, 0.0};
        bool newton = false;
        if (free[0] && free[1]) {
            // eigvalsh of the 2 x 2 block: its smaller eigenvalue
            const double mean = 0.5 * (h00 + h11), half = 0.5 * (h00 - h11);
            newton = mean - std::sqrt(half * half + h01 * h01) > 0.0;
            if (newton) {
                // LAPACK gesv: partial pivoting
                double a00 = h00, a01 = h01, a10 = h01, a11 = h11, r0 = -grad[0], r1 = -grad[1];
                bool swapped = std::abs(a10) > std::abs(a00);
                if (swapped) {
                    std::swap(a00, a10);
                    std::swap(a01, a11);
                    std::swap(r0, r1);
                }
                const double l = a10 / a00;
                const double u11 = a11 - l * a01;
                const double y1 = r1 - l * r0;
                step[1] = y1 / u11;
                step[0] = (r0 - a01 * step[1]) / a00;
            } else {
                step = {-grad[0], -grad[1]};
            }
        } else {
            const int k = free[0] ? 0 : 1;
            const auto K = static_cast<std::size_t>(k);
            const double h = k == 0 ? h00 : h11;
            newton = h > 0.0;
            step[K] = newton ? -grad[K] / h : -grad[K];
        }
        if (newton && std::max(std::abs(step[0]), std::abs(step[1])) <= NEWTON_ZONE) {
            const std::array<double, 2> trial{clip(p[0] + step[0], 0), clip(p[1] + step[1], 1)};
            if (std::max(std::abs(trial[0] - p[0]), std::abs(trial[1] - p[1])) <= THROAT_PARAMETER_TOL) return trial;
            p = trial;
            f = value(trial);
            continue;
        }
        double fraction = 1.0;
        std::array<double, 2> trial{};
        double f_trial = 0.0;
        while (true) {
            trial = {clip(p[0] + fraction * step[0], 0), clip(p[1] + fraction * step[1], 1)};
            f_trial = value(trial);
            if (f_trial <= f + ARMIJO * (grad[0] * (trial[0] - p[0]) + grad[1] * (trial[1] - p[1])) ||
                fraction < MIN_STEP_FRACTION)
                break;
            fraction *= 0.5;
        }
        if (std::max(std::abs(trial[0] - p[0]), std::abs(trial[1] - p[1])) <= THROAT_PARAMETER_TOL) return trial;
        if (f_trial > f) return p;  // no step lowers f: stationary to round-off
        p = trial;
        f = f_trial;
    }
    throw GeometryError("geometric throat refinement did not converge");
}

}  // namespace

ClosestPoints closest_points(const std::vector<Point>& curve_a, const std::vector<Point>& curve_b)
{
    // v1 closest_points / find_geometric_throat (Aungier 2003 Sec. 4.9): candidate starts at the local minima of the vertex-to-polyline distance
    // along either curve (within CANDIDATE_MARGIN of the smallest), each refined on the splines
    if (curve_a.size() < 2 || curve_b.size() < 2) throw std::invalid_argument("closest_points: needs two polylines");
    const CubicSpline A = chord_length_spline(curve_a), B = chord_length_spline(curve_b);
    const auto& ka = A.x();
    const auto& kb = B.x();
    const VertexDistances da = vertex_distances(curve_a, curve_b);
    const VertexDistances db = vertex_distances(curve_b, curve_a);
    const double ceiling = std::min(*std::min_element(da.distance.begin(), da.distance.end()),
                                    *std::min_element(db.distance.begin(), db.distance.end())) +
                           CANDIDATE_MARGIN;
    std::vector<std::array<double, 2>> starts;
    for (std::size_t i : local_minima(da.distance, ceiling)) {
        const std::size_t s = da.segment[i];
        starts.push_back({ka[i], kb[s] + da.along[i] * (kb[s + 1] - kb[s])});
    }
    for (std::size_t j : local_minima(db.distance, ceiling)) {
        const std::size_t s = db.segment[j];
        starts.push_back({ka[s] + db.along[j] * (ka[s + 1] - ka[s]), kb[j]});
    }
    std::optional<ClosestPoints> best;
    for (const auto& start : starts) {
        const auto uv = refine_closest(A, B, start);
        const Point pa = A(uv[0]), pb = B(uv[1]);
        const double distance = std::hypot(pa.m - pb.m, pa.y - pb.y);
        if (!best || distance < best->distance) best = ClosestPoints{distance, pa, pb};
    }
    return *best;  // the smallest vertex distance is always a candidate
}

namespace {

void check_cascade_inputs(double solidity, double stagger_deg)
{
    if (!(std::isfinite(solidity) && solidity > 0.0))
        throw std::invalid_argument("solidity must be a positive finite number (sigma = c / s), got " + repr(solidity));
    if (!(std::isfinite(stagger_deg) && -90.0 < stagger_deg && stagger_deg < 90.0))
        throw std::invalid_argument("stagger_deg must lie strictly between -90 and 90 degrees, got " + repr(stagger_deg));
}

bool finite(const Point& p) { return std::isfinite(p.m) && std::isfinite(p.y); }

/// The C2-relevant checks of v1 validate_cascade (see the file comment); every problem listed.
void validate(const Cascade& c, double solidity, double stagger_deg)
{
    const CascadeCoordinates& n = c.active;
    std::vector<std::string> problems;
    // _check_shapes
    const std::pair<const std::vector<Point>*, const char*> curves[] = {
        {&n.suction_surface_unstaggered, "suction_surface_unstaggered"},
        {&n.pressure_surface_unstaggered, "pressure_surface_unstaggered"},
        {&n.suction_surface, "suction_surface"},
        {&n.pressure_surface, "pressure_surface"}};
    for (const auto& [pts, name] : curves) {
        bool ok = pts->size() >= 3;
        for (const Point& p : *pts) ok = ok && finite(p);
        if (!ok) problems.push_back(std::string("active.") + name + " is not a finite point or (n, 2) curve");
    }
    const std::pair<Point, const char*> points[] = {{n.leading_edge_point, "leading_edge_point"},
                                                    {n.trailing_edge_point, "trailing_edge_point"},
                                                    {n.leading_edge_center, "leading_edge_center"},
                                                    {n.trailing_edge_center, "trailing_edge_center"},
                                                    {n.periodic_translation, "periodic_translation"},
                                                    {n.throat_point_reference, "throat_point_reference"},
                                                    {n.throat_point_neighbor, "throat_point_neighbor"}};
    for (const auto& [p, name] : points)
        if (!finite(p)) problems.push_back(std::string("active.") + name + " is not a finite point or (n, 2) curve");
    for (const auto& [x, name] : {std::pair<double, const char*>{n.chord, "chord"}, {n.pitch, "pitch"},
                                  {n.leading_edge_radius, "leading_edge_radius"},
                                  {n.trailing_edge_radius, "trailing_edge_radius"}, {n.throat_opening, "throat_opening"}})
        if (!std::isfinite(x)) problems.push_back(std::string("active.") + name + " is not finite");
    if (problems.empty()) {
        // _check_frame
        const Point te_expected{std::cos(c.stagger), std::sin(c.stagger)};
        if (!close(c.stagger, stagger_deg * (kPi / 180.0))) problems.emplace_back("stagger is not radians(spec.stagger_deg)");
        if (!close(n.chord, 1.0)) problems.push_back("normalized chord is " + repr(n.chord) + ", not 1");
        if (!close(n.pitch, 1.0 / solidity)) problems.push_back("pitch " + repr(n.pitch) + " is not 1 / solidity");
        if (!close(n.periodic_translation, Point{0.0, n.pitch})) problems.emplace_back("periodic_translation is not (0, pitch)");
        if (!(close(n.leading_edge_point, Point{0.0, 0.0}) && close(n.trailing_edge_point, te_expected)))
            problems.emplace_back("LE is not at the origin or TE is not at (cos(stagger), sin(stagger))");
        for (int s = 0; s < 2; ++s) {
            const char* name = s == 0 ? "suction_surface" : "pressure_surface";
            const auto& un = s == 0 ? n.suction_surface_unstaggered : n.pressure_surface_unstaggered;
            const auto& st = s == 0 ? n.suction_surface : n.pressure_surface;
            if (!(close(un.front(), Point{0.0, 0.0}) && close(un.back(), Point{1.0, 0.0})))
                problems.push_back(std::string(name) + "_unstaggered does not run from (0, 0) to (1, 0)");
            if (!(close(st.front(), n.leading_edge_point) && close(st.back(), n.trailing_edge_point)))
                problems.push_back(std::string(name) + " does not run from the LE point to the TE point");
            bool turned = st.size() == un.size();
            for (std::size_t k = 0; turned && k < st.size(); ++k) turned = close(st[k], rotate_by_stagger(un[k], c.stagger));
            if (!turned) problems.push_back(std::string(name) + " is not " + name + "_unstaggered turned by the stagger");
        }
        // _check_surfaces
        for (int s = 0; s < 2; ++s) {
            const auto& pts = s == 0 ? n.suction_surface_unstaggered : n.pressure_surface_unstaggered;
            double shortest = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k + 1 < pts.size(); ++k)
                shortest = std::min(shortest, std::hypot(pts[k + 1].m - pts[k].m, pts[k + 1].y - pts[k].y));
            if (shortest <= EXACT_TOL) problems.push_back(std::string(s == 0 ? "suction" : "pressure") + " surface repeats a point");
        }
        const auto blade = closed_loop(n.suction_surface_unstaggered, n.pressure_surface_unstaggered);
        if (paths_cross(blade, blade)) problems.emplace_back("suction and pressure surfaces cross");
        if (signed_area(blade) <= 0.0) problems.emplace_back("the suction surface is not on the +y side of the blade frame");
        // _check_periodicity
        const auto loop = closed_loop(n.suction_surface, n.pressure_surface);
        if (paths_cross(loop, translated(loop, n.periodic_translation))) problems.emplace_back("neighboring blades overlap");
        // _check_throat
        if (!(0.0 < n.throat_opening && n.throat_opening < n.pitch))
            problems.push_back("throat opening " + repr(n.throat_opening) + " is not between 0 and the pitch");
        if (!close(std::hypot(n.throat_point_reference.m - n.throat_point_neighbor.m,
                              n.throat_point_reference.y - n.throat_point_neighbor.y),
                   n.throat_opening))
            problems.emplace_back("the throat points are not one throat opening apart");
        if (!(close(c.throat_over_pitch, n.throat_opening / n.pitch) && close(c.throat_over_chord, n.throat_opening / n.chord)))
            problems.emplace_back("throat ratios do not match the throat opening");
        const double whole_blades = closest_points(loop, translated(loop, n.periodic_translation)).distance;
        if (whole_blades < n.throat_opening - THROAT_TOL)
            problems.push_back("the blades come closer (" + g(whole_blades, 9) + ") than the suction-to-pressure throat (" +
                               g(n.throat_opening, 9) + ") away from those surfaces");
        // _check_angles (kappa and camber consistency)
        if (!(close(c.camber_angle, c.kappa_le - c.kappa_te) && close(c.kappa_le, c.stagger + c.chi_le) &&
              close(c.kappa_te, c.stagger - c.chi_te)))
            problems.emplace_back("kappa_le, kappa_te, stagger and camber_angle are inconsistent");
        // _check_curves
        const SurfaceCurves& sc = *c.active_curves;
        for (int s = 0; s < 2; ++s) {
            const char* name = s == 0 ? "suction" : "pressure";
            const auto& u = s == 0 ? sc.suction_parameters : sc.pressure_parameters;
            const auto& st = s == 0 ? n.suction_surface : n.pressure_surface;
            const auto& un = s == 0 ? n.suction_surface_unstaggered : n.pressure_surface_unstaggered;
            bool ok = !u.empty() && u.front() == 0.0 && u.back() == 1.0 && u.size() == st.size();
            for (std::size_t k = 0; ok && k + 1 < u.size(); ++k) ok = u[k + 1] > u[k];
            if (!ok) {
                problems.push_back(std::string("active ") + name + " parameters do not run from 0 to 1, one per surface point");
                continue;
            }
            const ParametricBladeCurve& cc = s == 0 ? sc.suction : sc.pressure;
            const ParametricBladeCurve& cu = s == 0 ? sc.suction_unstaggered : sc.pressure_unstaggered;
            for (int f = 0; f < 2; ++f) {
                const auto& curve = f == 0 ? cc : cu;
                const auto& pts = f == 0 ? st : un;
                bool on = true;
                for (std::size_t k = 0; on && k < u.size(); ++k) on = close(curve.evaluate(u[k]), pts[k]);
                if (!on)
                    problems.push_back(std::string("active.") + name + "_surface" + (f == 0 ? "" : "_unstaggered") +
                                       " is not its curve sampled at its parameters");
            }
            bool turned = cc.control_points().size() == cu.control_points().size() && cc.knots() == cu.knots();
            for (std::size_t k = 0; turned && k < cc.control_points().size(); ++k)
                turned = close(cc.control_points()[k], rotate_by_stagger(cu.control_points()[k], c.stagger));
            if (!turned) problems.push_back(std::string("active ") + name + " curve is not the blade-frame curve turned by the stagger");
        }
    }
    if (!problems.empty()) {
        std::string text = "invalid cascade geometry:";
        for (const auto& p : problems) text += "\n- " + p;
        throw GeometryError(text);
    }
}

}  // namespace

Cascade make_cascade(const BladeSection& section, double solidity, double stagger_deg)
{
    // v1 make_cascade: the inputs, the one degree-to-radian conversion, the sampled blade-frame
    // surfaces (v1 _naca65_c2_profile: the C2 curves at their parameters), place_cascade, the
    // curves turned by the stagger, the blade angles and the checks
    check_cascade_inputs(solidity, stagger_deg);
    const double stagger = stagger_deg * (kPi / 180.0);
    if (section.suction.degree() < 1 || section.pressure.degree() < 1)
        throw std::invalid_argument("make_cascade: the blade section has no suction or pressure curve");
    if (!(std::isfinite(section.chi_le) && std::isfinite(section.chi_te)))
        throw std::invalid_argument("make_cascade: chi_le and chi_te must be finite");
    std::vector<Point> suction_un, pressure_un;
    for (double u : section.suction_parameters) suction_un.push_back(section.suction.evaluate(u));
    for (double u : section.pressure_parameters) pressure_un.push_back(section.pressure.evaluate(u));
    if (suction_un.size() < 2 || pressure_un.size() < 2)
        throw GeometryError("invalid cascade geometry:\n- the surfaces need at least two sample parameters each");

    // v1 place_cascade (Aungier 2003 Sec. 4.10: the blade turned by the stagger about its LE)
    const double pitch = 1.0 / solidity;
    const Point translation{0.0, pitch};
    const auto suction = rotated(suction_un, stagger);
    const auto pressure = rotated(pressure_un, stagger);
    const auto loop = closed_loop(suction, pressure);
    if (paths_cross(loop, translated(loop, translation)))
        throw GeometryError("neighboring blades overlap at solidity " + g(solidity, 17) + " and stagger " +
                            g(stagger * (180.0 / kPi)) + " deg");
    const ClosestPoints throat = closest_points(suction, translated(pressure, translation));

    Cascade c;
    c.designation = section.designation;
    c.solidity = solidity;
    c.stagger = stagger;
    c.chi_le = section.chi_le;
    c.chi_te = section.chi_te;
    c.camber_angle = section.chi_le + section.chi_te;
    c.kappa_le = stagger + section.chi_le;
    c.kappa_te = stagger - section.chi_te;
    CascadeCoordinates& a = c.active;
    a.chord = std::hypot(suction_un.back().m - suction_un.front().m, suction_un.back().y - suction_un.front().y);
    a.pitch = pitch;
    a.suction_surface_unstaggered = suction_un;
    a.pressure_surface_unstaggered = pressure_un;
    a.suction_surface = suction;
    a.pressure_surface = pressure;
    a.leading_edge_point = suction.front();
    a.trailing_edge_point = suction.back();
    a.leading_edge_center = rotate_by_stagger(section.le_center, stagger);
    a.leading_edge_radius = section.le_radius;
    a.trailing_edge_center = rotate_by_stagger(section.te_center, stagger);
    a.trailing_edge_radius = section.te_radius;
    a.periodic_translation = translation;
    a.throat_opening = throat.distance;
    a.throat_point_reference = throat.a;
    a.throat_point_neighbor = throat.b;
    c.throat_over_pitch = a.throat_opening / a.pitch;
    c.throat_over_chord = a.throat_opening / a.chord;
    c.active_curves = std::make_shared<const SurfaceCurves>(
        SurfaceCurves{section.suction.rotated(stagger), section.pressure.rotated(stagger), section.suction,
                      section.pressure, section.suction_parameters, section.pressure_parameters});
    validate(c, solidity, stagger_deg);
    return c;
}

Cascade make_cascade(const StreamSurfaceBlade& blade)
{
    // not in v1: the planar cascade of the turned-back section, on its stream sheet
    if (!blade.sheet) throw std::invalid_argument("make_cascade: the stream-surface blade has no stream sheet");
    if (blade.n_blades < 1) throw std::invalid_argument("make_cascade: the stream-surface blade needs its number of blades");
    if (!(std::isfinite(blade.chord) && blade.chord > 0.0))
        throw std::invalid_argument("make_cascade: the stream-surface blade's chord must be a positive length");
    Cascade c = make_cascade(blade.section, blade.solidity, blade.stagger_deg);
    const double pitch = blade.sheet->r_ref() * (2.0 * kPi / static_cast<double>(blade.n_blades));
    if (!(std::abs(c.active.pitch - pitch) <= 1e-12 * pitch))
        throw GeometryError("invalid cascade geometry:\n- the pitch 1 / solidity = " + repr(c.active.pitch) +
                            " is not r_ref 2 pi / N = " + repr(pitch) + " of the stream sheet");
    double lo = 0.0, hi = 0.0;
    for (const auto* pts : {&c.active.suction_surface, &c.active.pressure_surface})
        for (const Point& p : *pts) {
            lo = std::min(lo, p.m);
            hi = std::max(hi, p.m);
        }
    if (!(blade.sheet->m_front() < lo && hi < blade.sheet->m_back()))
        throw GeometryError("invalid cascade geometry:\n- the stream sheet (m " + g(blade.sheet->m_front()) + " .. " +
                            g(blade.sheet->m_back()) + ") does not cover the blade (m " + g(lo) + " .. " + g(hi) + ")");
    c.sheet = blade.sheet;
    return c;
}

}  // namespace BladeGen::Turbolenza
