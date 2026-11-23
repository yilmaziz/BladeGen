// Grid.hpp without v1 data: closed forms and defining equations of the grid generation (side
// lines, the linear-angle law and its primitive, uniform stations, turning points, station cuts),
// the refusals, Aungier's edge rule on every standard level, and rebuild_exit.  Tests marked
// [v1ref] use the v1 reference cascade as a realistic blade but compare with no v1 output.
#include "CascadeFixture.hpp"

#include "turbolenza/Grid.hpp"
#include "turbolenza/GridDetail.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <set>

namespace T = BladeGen::Turbolenza;
namespace D = BladeGen::Turbolenza::detail;
using fixture::kPi;

namespace {

double deg(double d) { return d * kPi / 180.0; }

std::size_t z(int i) { return static_cast<std::size_t>(i); }

const T::Cascade& synthetic_cascade()
{
    static const T::Cascade c = T::make_cascade(fixture::synthetic_section(), 1.25, 30.0);
    return c;
}

T::GridSettings explicit_settings(int n_passage, int n_eta)
{
    T::GridSettings s;
    s.n_passage = n_passage;
    s.n_eta = n_eta;
    return s;
}

/// Gauss-Legendre (5 points) composite integral of f over [a, b] in n pieces.
template <class F>
double integrate(F f, double a, double b, int n)
{
    const double x[5] = {0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640};
    const double w[5] = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891,
                         0.2369268850561891};
    double s = 0.0;
    const double h = (b - a) / n;
    for (int k = 0; k < n; ++k) {
        const double c = a + (k + 0.5) * h;
        for (int q = 0; q < 5; ++q) s += w[q] * f(c + 0.5 * h * x[q]) * 0.5 * h;
    }
    return s;
}

/// Every invariant a v2 grid must show whatever built it (v1 grid.py section 6).
void check_invariants(const T::Grid& g)
{
    const int n = g.n_m(), ne = g.n_eta();
    const double pitch = g.boundaries.periodic_translation.y;
    // constant dm and deta (Aungier p. 87)
    for (int i = 0; i + 1 < n; ++i) CHECK(std::abs((g.m[z(i + 1)] - g.m[z(i)]) - g.dm()) <= 1e-13 * g.chord());
    CHECK(g.eta.front() == 0.0);
    CHECK(g.eta.back() == 1.0);
    for (int i = 0; i < n; ++i) {
        CHECK(g.node_m(i, 0) == g.m[z(i)]);
        // periodic stations: the upper side line is the lower one moved one pitch, exactly
        if (g.region(i) != T::Region::Passage) CHECK(g.node_y(i, ne - 1) == g.node_y(i, 0) + pitch);
        // the nodes are the linear map of Eq. 5-2 between the two wall values
        for (int j = 0; j < ne; ++j)
            CHECK(g.node_y(i, j) == (1.0 - g.eta[z(j)]) * g.node_y(i, 0) + g.eta[z(j)] * g.node_y(i, ne - 1));
    }
    // the wall nodes are curve points at their recorded u
    for (int b = 0; b <= g.last_blade_i - g.first_blade_i; ++b) {
        const int i = g.first_blade_i + b;
        CHECK(std::abs(g.wall_d[0].lower[z(b)].m - g.m[z(i)]) <= 1e-12);
        CHECK(g.wall_d[0].lower[z(b)].y == g.node_y(i, 0));
        CHECK(std::abs(g.wall_d[0].upper[z(b)].y - g.node_y(i, ne - 1)) <= 1e-15);
    }
}

}  // namespace

TEST_CASE("the mean of tan over [a, b] is the exact primitive of the linear-angle law", "[turbolenza][grid]")
{
    // T(a, b) = ln(cos a / cos b) / (b - a); equal angles give tan(b)
    for (double b : {-1.2, -0.3, 0.0, 0.31, 0.84, 1.4}) {
        CHECK(D::mean_tan(b, b) == Catch::Approx(std::tan(b)).epsilon(1e-15));
        for (double a : {-1.0, -0.2, 0.1, 0.5, 1.3}) {
            const double closed = std::log(std::cos(a) / std::cos(b)) / (b - a);
            CHECK(D::mean_tan(a, b) == Catch::Approx(closed).epsilon(1e-12));
            // the mean of tan by quadrature
            const double mean = integrate([](double x) { return std::tan(x); }, std::min(a, b), std::max(a, b), 200) /
                                std::abs(b - a);
            CHECK(D::mean_tan(a, b) == Catch::Approx(mean).epsilon(1e-12));
        }
        // no cancellation as a -> b: mean = tan(b) - d sec^2(b) / 2 + O(d^2), d = b - a
        for (double d : {1e-4, 1e-8, 1e-12}) {
            const double sec2 = 1.0 / (std::cos(b) * std::cos(b));
            CHECK(std::abs(D::mean_tan(b - d, b) - (std::tan(b) - 0.5 * d * sec2)) <= 100.0 * d * d + 1e-14);
        }
    }
}

TEST_CASE("side lines: the linear-angle law, its primitive, the far-station angle, one pitch apart", "[turbolenza][grid]")
{
    const T::Cascade& c = synthetic_cascade();
    const double inlet = deg(52.0), exit = deg(-5.0);
    const T::PassageBoundaries b = T::make_passage_boundaries(c, explicit_settings(30, 11), inlet, exit);
    // the far stations lie along the flow angle; the junctions take the blade angles
    CHECK(b.upstream_beta.front() == inlet);
    CHECK(b.upstream_beta.back() == c.kappa_le);
    CHECK(b.downstream_beta.front() == c.kappa_te);
    CHECK(b.downstream_beta.back() == exit);
    CHECK(b.upstream_m.front() == b.m_in);
    CHECK(b.upstream_m.back() == b.m_le);
    CHECK(b.downstream_m.front() == b.m_te);
    CHECK(b.downstream_m.back() == b.m_out);
    CHECK(b.m_in == b.m_le - 1.0);
    CHECK(b.m_out == b.m_te + 1.0);
    // the lines end on the exact LE and TE points
    CHECK(b.upstream_lower_boundary.back().y == c.active.leading_edge_point.y);
    CHECK(b.downstream_lower_boundary.front().y == c.active.trailing_edge_point.y);
    struct Side
    {
        const std::vector<double>& m;
        const std::vector<double>& beta;
        const std::vector<T::Point>& lower;
        const std::vector<T::Point>& upper;
        double m_j, y_j, kappa, m_far, far;
    };
    const Side sides[2] = {
        {b.upstream_m, b.upstream_beta, b.upstream_lower_boundary, b.upstream_upper_boundary, b.m_le,
         c.active.leading_edge_point.y, c.kappa_le, b.m_in, inlet},
        {b.downstream_m, b.downstream_beta, b.downstream_lower_boundary, b.downstream_upper_boundary, b.m_te,
         c.active.trailing_edge_point.y, c.kappa_te, b.m_out, exit}};
    for (const Side& s : sides) {
        for (std::size_t k = 0; k < s.m.size(); ++k) {
            // beta linear in m (Aungier 2003 Sec. 5.3, p. 86)
            const double f = (s.m[k] - s.m_j) / (s.m_far - s.m_j);
            CHECK(s.beta[k] == Catch::Approx(s.kappa + f * (s.far - s.kappa)).margin(1e-15));
            // y the integral of tan(beta(m)) from the junction
            auto slope = [&](double m) { return std::tan(s.kappa + (m - s.m_j) / (s.m_far - s.m_j) * (s.far - s.kappa)); };
            const double rise = s.m[k] == s.m_j ? 0.0
                                                : integrate(slope, std::min(s.m_j, s.m[k]), std::max(s.m_j, s.m[k]), 64) *
                                                      (s.m[k] > s.m_j ? 1.0 : -1.0);
            CHECK(s.lower[k].y == Catch::Approx(s.y_j + rise).margin(1e-13));
            // the upper line is the lower one moved one pitch, never integrated on its own
            CHECK(s.upper[k].m == s.lower[k].m);
            CHECK(s.upper[k].y == s.lower[k].y + c.active.pitch);
        }
        // dy/dm at the far station is tan(flow angle) (central difference of the closed form)
        const double h = 1e-5 * std::abs(s.m_far - s.m_j);
        auto y = [&](double m) {
            return D::side_line(m, (m - s.m_j) / (s.m_far - s.m_j), {s.m_j, s.y_j}, s.kappa, s.far)[1];
        };
        CHECK((y(s.m_far + h) - y(s.m_far - h)) / (2.0 * h) == Catch::Approx(std::tan(s.far)).epsilon(1e-8));
        // at the junction the slope is tan(kappa)
        CHECK((y(s.m_j + h) - y(s.m_j - h)) / (2.0 * h) == Catch::Approx(std::tan(s.kappa)).epsilon(1e-8));
    }
}

TEST_CASE("polynomial_roots isolates every real root of odd multiplicity", "[turbolenza][grid]")
{
    // (s - r1)(s - r2)(s - r3)(s - r4) in ascending coefficients
    auto poly = [](std::vector<double> r) {
        std::vector<double> c{1.0};
        for (double x : r) {
            std::vector<double> n(c.size() + 1, 0.0);
            for (std::size_t k = 0; k < c.size(); ++k) {
                n[k + 1] += c[k];
                n[k] -= x * c[k];
            }
            c = n;
        }
        return c;
    };
    const auto a = D::polynomial_roots(poly({0.1, 0.35, 0.3500001, 0.8}), 1.0);
    REQUIRE(a.size() == 4);
    CHECK(a[0] == Catch::Approx(0.1).margin(1e-14));
    CHECK(a[1] == Catch::Approx(0.35).margin(1e-11));
    CHECK(a[2] == Catch::Approx(0.3500001).margin(1e-11));
    CHECK(a[3] == Catch::Approx(0.8).margin(1e-14));
    // roots outside [0, h] are not reported; ends are
    const auto b = D::polynomial_roots(poly({-0.5, 0.0, 0.25, 2.0}), 0.25);
    REQUIRE(b.size() == 2);
    CHECK(b[0] == 0.0);
    CHECK(b[1] == Catch::Approx(0.25).margin(1e-15));
    // a triple root changes sign: found; a quartic without real roots: none
    const auto c = D::polynomial_roots(poly({0.3, 0.3, 0.3, 0.9}), 1.0);
    REQUIRE(!c.empty());
    CHECK(c.front() == Catch::Approx(0.3).margin(1e-5));
    CHECK(D::polynomial_roots({1.0, 0.0, 2.0, 0.0, 1.0}, 1.0).empty());  // (1 + s^2)^2
    CHECK(D::polynomial_roots({0.0, 0.0, 0.0}, 1.0).empty());            // identically zero (PPoly: NaN, dropped)
    CHECK(D::polynomial_roots({3.0}, 1.0).empty());
    const auto e = D::polynomial_roots({-0.5, 2.0}, 1.0);
    REQUIRE(e.size() == 1);
    CHECK(e[0] == 0.25);
}

TEST_CASE("turning points of a B-spline with interior knots agree with a dense sign scan", "[turbolenza][grid]")
{
    // a degree-5 clamped B-spline with five interior knots: m(u) and y(u) both turn several times
    std::vector<double> t{0, 0, 0, 0, 0, 0, 0.15, 0.3, 0.45, 0.7, 0.85, 1, 1, 1, 1, 1, 1};
    std::vector<T::Point> P{{0.0, 0.0}, {-0.2, 0.3}, {0.5, 0.6}, {0.1, -0.4}, {0.9, 0.2},
                            {0.4, 0.8},  {1.2, -0.1}, {0.8, 0.5}, {1.3, 0.3}, {1.0, 0.0}, {1.5, 0.2}};
    const T::ParametricBladeCurve curve(t, P, 1.0, 5);
    for (int axis = 0; axis < 2; ++axis) {
        INFO("axis " << axis);
        const auto tp = D::turning_points(curve, axis);
        // the sign changes of the derivative on a fine scan, each bisected
        auto d = [&](double u) { return axis == 0 ? curve.derivative(u, 1).m : curve.derivative(u, 1).y; };
        std::vector<double> scan;
        const int N = 20000;
        for (int k = 0; k < N; ++k) {
            double a = static_cast<double>(k) / N, b = static_cast<double>(k + 1) / N;
            if ((d(a) < 0.0) != (d(b) < 0.0)) {
                for (int it = 0; it < 60; ++it) {
                    const double mid = 0.5 * (a + b);
                    ((d(mid) < 0.0) == (d(a) < 0.0) ? a : b) = mid;
                }
                scan.push_back(0.5 * (a + b));
            }
        }
        REQUIRE(tp.size() == scan.size());
        REQUIRE(!tp.empty());  // m turns once, y three times
        for (std::size_t k = 0; k < tp.size(); ++k) {
            CHECK(tp[k] == Catch::Approx(scan[k]).margin(1e-14));
            CHECK(d(tp[k] - 1e-9) * d(tp[k] + 1e-9) < 0.0);
        }
    }
}

TEST_CASE("station cuts: every crossing is a curve point on the station line, branches monotone", "[turbolenza][grid]")
{
    const T::Cascade& c = synthetic_cascade();
    const auto b = T::make_passage_boundaries(c, explicit_settings(20, 5), deg(45.0), std::nullopt);
    const D::CurveRow row(b);
    // branches: m and y strictly monotone between their ends
    for (int k = 0; k < row.n_branches(); ++k) {
        double pm = row.point(k, row.u0[z(k)]).m, py = row.point(k, row.u0[z(k)]).y;
        const double sm = row.m1[z(k)] > row.m0[z(k)] ? 1.0 : -1.0;
        const double ey = row.point(k, row.u1[z(k)]).y, sy = ey > py ? 1.0 : -1.0;
        for (int q = 1; q <= 200; ++q) {
            const auto p = row.point(k, row.u0[z(k)] + (row.u1[z(k)] - row.u0[z(k)]) * q / 200.0);
            CHECK(sm * (p.m - pm) > 0.0);
            CHECK(sy * (p.y - py) > 0.0);
            pm = p.m;
            py = p.y;
        }
    }
    const auto m = D::linspace(row.m_tip + 1e-6, row.m_tail - 1e-6, 157);
    const D::Cut cut = row.cut(m);
    for (std::size_t i = 0; i < m.size(); ++i) {
        CHECK(cut.count[i] % 2 == 0);
        CHECK(cut.count[i] >= 2);
        for (int j = 0; j < cut.count[i]; ++j) {
            const auto p = row.point(cut.segment[i][z(j)], cut.u[i][z(j)]);
            CHECK(std::abs(p.m - m[i]) <= 1e-15);
            CHECK(p.y == cut.height[i][z(j)]);
            if (j > 0) CHECK(cut.height[i][z(j)] >= cut.height[i][z(j - 1)]);
        }
        CHECK(cut.metal[i] > 0.0);
        CHECK(cut.metal[i] < b.pitch);
    }
    // a line ahead of the nose tip or behind the tail meets nothing
    const D::Cut out = row.cut({row.m_tip - 1e-9, row.m_tail + 1e-9});
    CHECK(out.count[0] == 0);
    CHECK(out.count[1] == 0);
}

TEST_CASE("an explicit grid of the synthetic cascade keeps every invariant", "[turbolenza][grid]")
{
    const T::Cascade& c = synthetic_cascade();
    for (int n : {12, 25, 40}) {
        INFO("n_passage " << n);
        const T::Grid g = T::make_grid(c, explicit_settings(n, 9), deg(48.0), std::nullopt);
        check_invariants(g);
        CHECK(g.m.front() == g.boundaries.m_in);
        CHECK(g.m.back() == g.boundaries.m_out);
        // dm near (m_te - m_le) / n_passage (v1 section 2)
        const double target = (g.m_te() - g.m_le()) / n;
        CHECK(std::abs(g.dm() - target) <= target / (2.0 * (g.n_m() - 1)) + 1e-15);
        // blade stations: exactly the stations whose line meets the blade with metal
        const D::CurveRow row(g.boundaries);
        const D::Cut cut = row.cut(g.m);
        for (int i = 0; i < g.n_m(); ++i) CHECK(cut.cuts_blade(z(i)) == (g.region(i) == T::Region::Passage));
    }
}

TEST_CASE("make_grid refuses: short extensions, bad settings, bad angles", "[turbolenza][grid]")
{
    const T::Cascade& c = synthetic_cascade();
    // an inlet extension shorter than the nose: station 0 cuts the blade
    const auto ext = D::CurveRow(T::make_passage_boundaries(c, explicit_settings(30, 9), 0.5, std::nullopt)).extremes();
    REQUIRE(c.active.leading_edge_point.m - ext[0] > 1e-4);  // the nose reaches ahead of the LE point
    auto s = explicit_settings(30, 9);
    s.upstream_extension_chords = 0.5 * (c.active.leading_edge_point.m - ext[0]);
    CHECK_THROWS_AS(T::make_grid(c, s, 0.5, std::nullopt), T::GridError);
    s = explicit_settings(30, 9);
    s.upstream_extension_chords = 0.0;
    CHECK_THROWS_AS(T::make_grid(c, s, 0.5, std::nullopt), std::invalid_argument);
    CHECK_THROWS_AS(T::make_grid(c, explicit_settings(30, 1), 0.5, std::nullopt), std::invalid_argument);
    CHECK_THROWS_AS(T::make_grid(c, explicit_settings(0, 9), 0.5, std::nullopt), std::invalid_argument);
    T::GridSettings missing;
    missing.n_eta = 9;
    CHECK_THROWS_AS(T::make_grid(c, missing, 0.5, std::nullopt), std::invalid_argument);
    CHECK_THROWS_AS(T::make_grid(c, explicit_settings(30, 9), 0.5 * kPi, std::nullopt), std::invalid_argument);
    CHECK_THROWS_AS(T::make_grid(c, explicit_settings(30, 9), 0.5, -0.5 * kPi), std::invalid_argument);
    CHECK_THROWS_AS(T::make_grid(c, explicit_settings(30, 9), std::nan(""), std::nullopt), std::invalid_argument);
    // a grid so coarse that no station meets the blade
    CHECK_THROWS_AS(T::make_grid(c, explicit_settings(1, 9), 0.5, std::nullopt), T::GridError);
    // a level together with n_passage, an unknown level
    T::GridSettings both;
    both.level = "coarse";
    both.n_passage = 20;
    CHECK_THROWS_AS(T::make_grid(c, both, 0.5, std::nullopt), std::invalid_argument);
    T::GridSettings unknown;
    unknown.level = "ultra";
    CHECK_THROWS_AS(T::make_grid(c, unknown, 0.5, std::nullopt), std::invalid_argument);
    CHECK_THROWS_AS(T::rebuild_exit(T::make_grid(c, explicit_settings(30, 9), 0.5, std::nullopt), 2.0), std::invalid_argument);
}

TEST_CASE("the Grid refuses nodes that break its invariants, whoever built them", "[turbolenza][grid]")
{
    const T::Cascade& c = synthetic_cascade();
    const T::Grid good = T::make_grid(c, explicit_settings(30, 9), 0.5, std::nullopt);
    const D::CurveRow row(good.boundaries);
    {
        T::Grid g = good;  // blade run shifted by one station
        g.first_blade_i += 1;
        CHECK_THROWS_AS(D::finish_grid(g, row), T::GridError);
    }
    {
        T::Grid g = good;  // a wall node off its curve point
        g.node_y(g.first_blade_i + 2, 0) += 1e-9;
        CHECK_THROWS_AS(D::finish_grid(g, row), T::GridError);
    }
    {
        T::Grid g = good;  // a closed station
        const int i = g.first_blade_i - 1;
        g.node_y(i, g.n_eta() - 1) = g.node_y(i, 0);
        CHECK_THROWS_AS(D::finish_grid(g, row), T::GridError);
    }
    T::Grid g = good;
    CHECK_NOTHROW(D::finish_grid(g, row));
}

TEST_CASE("rebuild_exit moves only the exit stations and equals a fresh build", "[turbolenza][grid]")
{
    const T::Cascade& c = synthetic_cascade();
    const auto s = explicit_settings(28, 11);
    const T::Grid g = T::make_grid(c, s, deg(50.0), std::nullopt);
    const T::Grid e = T::rebuild_exit(g, deg(12.0));
    const T::Grid fresh = T::make_grid(c, s, deg(50.0), deg(12.0));
    CHECK(e.boundaries.exit_angle == deg(12.0));
    CHECK(e.boundaries.downstream_beta.back() == deg(12.0));
    for (int i = 0; i < g.n_m(); ++i)
        for (int j = 0; j < g.n_eta(); ++j) {
            if (i <= g.last_blade_i) CHECK(e.node_y(i, j) == g.node_y(i, j));
            CHECK(e.node_y(i, j) == fresh.node_y(i, j));
        }
    CHECK(e.boundaries.downstream_lower_boundary.size() == fresh.boundaries.downstream_lower_boundary.size());
    for (std::size_t k = 0; k < e.boundaries.downstream_lower_boundary.size(); ++k)
        CHECK(e.boundaries.downstream_lower_boundary[k].y == fresh.boundaries.downstream_lower_boundary[k].y);
    check_invariants(e);
}

TEST_CASE("GridSettings with a level routes make_grid to the standard grid", "[turbolenza][grid]")
{
    const T::Cascade& c = synthetic_cascade();
    T::GridSettings s;
    s.level = "coarse";
    const T::Grid g = T::make_grid(c, s, deg(45.0), std::nullopt, true);
    const T::StandardGrid sg = T::build_grid(c, "coarse", deg(45.0), std::nullopt);
    CHECK(g.m == sg.grid.m);
    CHECK(g.node_y.v == sg.grid.node_y.v);
    CHECK(T::level_target_dm_over_chord("medium") == 0.02);
    CHECK(T::level_target_dm_over_chord("very_coarse") == Catch::Approx(0.04).epsilon(1e-15));
    CHECK(T::level_target_dm_over_chord("very_fine") == Catch::Approx(0.01).epsilon(1e-15));
    CHECK_THROWS_AS(T::level_target_dm_over_chord("x"), std::invalid_argument);
}

TEST_CASE("standard levels: Aungier's edge rule and the section 9 invariants on every level", "[turbolenza][grid]")
{
    // the synthetic cascade, and the v1 reference cascade when its data is there
    std::vector<const T::Cascade*> cascades{&synthetic_cascade()};
    if (const auto* R = fixture::reference()) cascades.push_back(&R->cascade);
    for (const T::Cascade* c : cascades) {
        INFO(c->designation);
        const T::PassageTransitions tr = T::passage_transitions(*c);
        double previous_dm = 1e9, previous_deta = 1e9;
        for (const char* level : T::kGridLevels) {
            INFO(level);
            const T::StandardGrid sg = T::build_grid(*c, level, deg(45.0), std::nullopt);
            const T::Grid& g = sg.grid;
            check_invariants(g);
            // the first and the last blade station each cut one suction and one pressure surface
            // (Aungier 2003 pp. 89-90, Fig. 5-5), off the LE and TE points
            const auto edges = T::wall_station_topology(g);
            for (const auto& e : edges) {
                CHECK(e.valid());
                CHECK(e.crossings == 2);
                CHECK(e.suction == 1);
                CHECK(e.pressure == 1);
            }
            CHECK_NOTHROW(T::validate_uniform_grid_topology(g));
            for (std::size_t b = 0; b < g.wall_identity.lower.surface.size(); ++b) {
                CHECK(g.wall_identity.lower.surface[b] == T::Surface::Suction);
                CHECK(g.wall_identity.upper.surface[b] == T::Surface::Pressure);
            }
            // f - 1 and l + 1 meet no blade; the blade stations lie inside the SS + PS interval
            const D::Cut cut = D::CurveRow(g.boundaries).cut(g.m);
            CHECK(cut.count[z(g.first_blade_i - 1)] == 0);
            CHECK(cut.count[z(g.last_blade_i + 1)] == 0);
            CHECK(tr.valid_cut(g.m_first_blade()));
            CHECK(tr.valid_cut(g.m_last_blade()));
            // the lattice: cell faces at the window middles, dm within the band of the target
            CHECK(g.m_first_blade() == Catch::Approx(tr.le_face() + 0.5 * g.dm()).margin(1e-13));
            CHECK(g.m_last_blade() == Catch::Approx(tr.te_face() - 0.5 * g.dm()).margin(1e-13));
            const double target = T::level_target_dm_over_chord(level);
            CHECK(std::abs(g.dm() / target - 1.0) <= T::kCandidateBand + 1e-12);
            // eta from CELL_ASPECT, extensions within tolerance, both falling from level to level
            CHECK(g.n_eta() - 1 == std::llround(g.pitch() / (T::kCellAspect * g.dm())));
            CHECK(std::abs(sg.info.upstream_extension_chords - 1.0) <= T::kExtensionToleranceChords);
            CHECK(std::abs(sg.info.downstream_extension_chords - 1.0) <= T::kExtensionToleranceChords);
            CHECK(g.dm() < previous_dm);
            CHECK(g.deta() < previous_deta);
            previous_dm = g.dm();
            previous_deta = g.deta();
            CHECK(T::validate_standard_grid(g) == sg.info.checks);
            CHECK(sg.info.checks.size() == 12);
        }
    }
}

TEST_CASE("an inadmissible explicit n_passage is refused with the nearest admissible values", "[turbolenza][grid][v1ref]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    const T::Cascade& c = R->cascade;
    const double inlet = deg(45.0);
    std::vector<int> range;
    for (int n = 10; n <= 80; ++n) range.push_back(n);
    const auto valid = T::find_valid_n_passage(c, explicit_settings(1, 2), inlet, range);
    REQUIRE(!valid.empty());
    const std::set<int> ok(valid.begin(), valid.end());
    int bad = -1;
    for (int n = valid.front() + 1; n < valid.back(); ++n)
        if (!ok.count(n)) {
            bad = n;
            break;
        }
    REQUIRE(bad > 0);
    // without require_ss_ps the grid builds (the tracked-contour walls describe it) ...
    CHECK_NOTHROW(T::make_grid(c, explicit_settings(bad, 11), inlet, std::nullopt));
    // ... with it, GridTopologyError naming the nearest admissible n_passage on either side
    try {
        (void)T::make_grid(c, explicit_settings(bad, 11), inlet, std::nullopt, true);
        FAIL("an inadmissible n_passage was accepted");
    } catch (const T::GridTopologyError& e) {
        INFO(e.what());
        CHECK(e.n_passage == bad);
        CHECK(!e.edge.empty());
        const int lower = *std::prev(ok.lower_bound(bad));
        const int higher = *ok.upper_bound(bad);
        CHECK(e.nearest_lower == lower);
        CHECK(e.nearest_higher == higher);
        CHECK(std::string(e.what()).find("Nearest admissible: lower " + std::to_string(lower) + ", higher " +
                                         std::to_string(higher)) != std::string::npos);
        REQUIRE(e.ceiling);
        REQUIRE(e.finest);
        CHECK(*e.finest < *e.ceiling);
        CHECK(e.ceiling == T::n_passage_ceiling(c, explicit_settings(1, 2), inlet));
        // the nearest admissible values do build with require_ss_ps
        for (int n : {lower, higher}) CHECK_NOTHROW(T::make_grid(c, explicit_settings(n, 11), inlet, std::nullopt, true));
    }
    // every n_passage from the ceiling on fails (dm <= the LE window)
    const auto ceiling = T::n_passage_ceiling(c, explicit_settings(1, 2), inlet);
    REQUIRE(ceiling);
    for (int n : {*ceiling, *ceiling + 7}) CHECK_THROWS_AS(T::make_grid(c, explicit_settings(n, 5), inlet, std::nullopt, true),
                                                           T::GridTopologyError);
}

TEST_CASE("the v1 blade at a solidity so high that the blades overlap is refused", "[turbolenza][grid][v1ref]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    const auto s = fixture::v1_section(R->family.at("cascade"));
    CHECK_THROWS_AS(T::make_cascade(s, 30.0, 33.0), T::GeometryError);
    CHECK_THROWS_AS(T::make_cascade(s, 30.0, 0.0), T::GeometryError);
}
