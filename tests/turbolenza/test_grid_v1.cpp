// Grid.hpp against TURBOlenza v1: the explicit 45 x 21 grid, the five standard levels with their
// info and rebuilt exits, and the passage transitions of the reference cascade (NACA
// 65-(12A10)10, C2, solidity 1, stagger 33 deg, inlet 45 deg, exit kappa_te).
#include "CascadeFixture.hpp"

#include "turbolenza/Grid.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <sstream>
#include <string>

namespace T = BladeGen::Turbolenza;

namespace {

double max_abs(const std::vector<double>& a, const std::vector<double>& b)
{
    REQUIRE(a.size() == b.size());
    double d = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) d = std::max(d, std::abs(a[k] - b[k]));
    return d;
}

double max_abs(const std::vector<T::Point>& a, const std::vector<T::Point>& b)
{
    REQUIRE(a.size() == b.size());
    double d = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) d = std::max({d, std::abs(a[k].m - b[k].m), std::abs(a[k].y - b[k].y)});
    return d;
}

/// The largest differences of a port grid against a v1 grid, by kind.
struct GridDiff
{
    double m = 0, eta = 0, nodes = 0, wall_u = 0, wall_contour = 0, boundaries = 0, wall_d0 = 0, wall_d1 = 0, wall_d2 = 0;
    void merge(const GridDiff& o)
    {
        m = std::max(m, o.m);
        eta = std::max(eta, o.eta);
        nodes = std::max(nodes, o.nodes);
        wall_u = std::max(wall_u, o.wall_u);
        wall_contour = std::max(wall_contour, o.wall_contour);
        boundaries = std::max(boundaries, o.boundaries);
        wall_d0 = std::max(wall_d0, o.wall_d0);
        wall_d1 = std::max(wall_d1, o.wall_d1);
        wall_d2 = std::max(wall_d2, o.wall_d2);
    }
};

GridDiff compare(const T::Grid& g, const ref::json& v, const ref::json& vb)
{
    const T::Grid r = ref::grid(v);
    // integers and names exactly
    REQUIRE(g.n_m() == r.n_m());
    REQUIRE(g.n_eta() == r.n_eta());
    CHECK(g.first_blade_i == r.first_blade_i);
    CHECK(g.last_blade_i == r.last_blade_i);
    CHECK(g.wall_identity.stations == r.wall_identity.stations);
    CHECK(g.wall_identity.lower.surface == r.wall_identity.lower.surface);
    CHECK(g.wall_identity.upper.surface == r.wall_identity.upper.surface);
    CHECK(g.wall_identity.lower.hidden == r.wall_identity.lower.hidden);
    CHECK(g.wall_identity.upper.hidden == r.wall_identity.upper.hidden);
    GridDiff d;
    d.m = max_abs(g.m, r.m);
    d.eta = max_abs(g.eta, r.eta);
    d.nodes = std::max(max_abs(g.node_m.v, r.node_m.v), max_abs(g.node_y.v, r.node_y.v));
    d.wall_u = std::max({max_abs(g.wall_identity.lower.u, r.wall_identity.lower.u),
                         max_abs(g.wall_identity.upper.u, r.wall_identity.upper.u),
                         max_abs(g.wall_identity.lower.parameter, r.wall_identity.lower.parameter),
                         max_abs(g.wall_identity.upper.parameter, r.wall_identity.upper.parameter)});
    d.wall_contour = std::max(max_abs(g.wall_identity.lower.contour, r.wall_identity.lower.contour),
                              max_abs(g.wall_identity.upper.contour, r.wall_identity.upper.contour));
    // wall_d[k]: relative to the largest |X^(k)| over the wall nodes (as test_curve.cpp)
    double* wd[3] = {&d.wall_d0, &d.wall_d1, &d.wall_d2};
    for (int k = 0; k < 3; ++k) {
        double scale = 0.0;
        for (const auto* w : {&r.wall_d[k].lower, &r.wall_d[k].upper})
            for (const auto& p : *w) scale = std::max(scale, std::hypot(p.m, p.y));
        *wd[k] = std::max(max_abs(g.wall_d[k].lower, r.wall_d[k].lower), max_abs(g.wall_d[k].upper, r.wall_d[k].upper)) /
                 scale;
    }
    // boundaries
    const auto& b = g.boundaries;
    CHECK(b.chord == ref::num(vb.at("chord")));
    CHECK(b.pitch == ref::num(vb.at("pitch")));
    CHECK(b.inlet_angle == ref::num(vb.at("inlet_angle")));
    CHECK(b.exit_angle == ref::num(vb.at("exit_angle")));
    double bd = 0.0;
    for (const char* name : {"m_in", "m_le", "m_te", "m_out"}) {
        const double x = std::string(name) == "m_in" ? b.m_in
                         : std::string(name) == "m_le" ? b.m_le
                         : std::string(name) == "m_te" ? b.m_te
                                                         : b.m_out;
        bd = std::max(bd, std::abs(x - ref::num(vb.at(name))));
    }
    bd = std::max(bd, max_abs(b.lower_blade_boundary, ref::points(vb.at("lower_blade_boundary"))));
    bd = std::max(bd, max_abs(b.upper_blade_boundary, ref::points(vb.at("upper_blade_boundary"))));
    bd = std::max(bd, max_abs(b.upstream_m, ref::vec(vb.at("upstream_m"))));
    bd = std::max(bd, max_abs(b.upstream_beta, ref::vec(vb.at("upstream_beta"))));
    bd = std::max(bd, max_abs(b.upstream_lower_boundary, ref::points(vb.at("upstream_lower_boundary"))));
    bd = std::max(bd, max_abs(b.upstream_upper_boundary, ref::points(vb.at("upstream_upper_boundary"))));
    bd = std::max(bd, max_abs(b.downstream_m, ref::vec(vb.at("downstream_m"))));
    bd = std::max(bd, max_abs(b.downstream_beta, ref::vec(vb.at("downstream_beta"))));
    bd = std::max(bd, max_abs(b.downstream_lower_boundary, ref::points(vb.at("downstream_lower_boundary"))));
    bd = std::max(bd, max_abs(b.downstream_upper_boundary, ref::points(vb.at("downstream_upper_boundary"))));
    d.boundaries = bd;
    CHECK(b.kappa_le == b.upstream_beta.back());
    CHECK(b.kappa_te == b.downstream_beta.front());
    REQUIRE(b.blade_curves);
    return d;
}

std::string text(const char* what, const GridDiff& d)
{
    std::ostringstream o;
    o << what << ": m " << d.m << ", eta " << d.eta << ", nodes " << d.nodes << ", wall u " << d.wall_u << ", contour "
      << d.wall_contour << ", boundaries " << d.boundaries << ", wall_d " << d.wall_d0 << " " << d.wall_d1 << " "
      << d.wall_d2;
    return o.str();
}

}  // namespace

TEST_CASE("the explicit 45 x 21 grid matches v1", "[turbolenza][grid][v1ref]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    const ref::json& fam = R->family;
    T::GridSettings s;
    s.n_passage = 45;
    s.n_eta = 21;
    const double beta1 = ref::num(fam.at("inlet_angle")), beta2 = ref::num(fam.at("exit_angle"));
    const T::Grid g = T::make_grid(R->cascade, s, beta1, beta2);
    const GridDiff d = compare(g, fam.at("explicit_45x21"), fam.at("explicit_45x21_boundaries"));
    INFO(text("explicit 45 x 21 vs v1", d));
    // measured (MSVC 19.51, v1 @ 67599a6 with numpy 2.4.4): m, eta and boundaries 0 (bit for bit),
    // nodes 1.1e-16, wall u 5.6e-17, contour 0; wall_d relative to the largest |X^(k)| at the wall
    // nodes: 6.4e-17, 2.7e-14, 1.7e-10 (the basis recurrence against scipy's de Boor on the short
    // nose spans, as test_curve.cpp).  Bounds about 10x the measurement.
    CHECK(d.m <= 1e-15);
    CHECK(d.eta == 0.0);
    CHECK(d.nodes <= 2e-15);
    CHECK(d.wall_u <= 1e-15);
    CHECK(d.wall_contour <= 1e-15);
    CHECK(d.boundaries <= 1e-15);
    CHECK(d.wall_d0 <= 1e-15);
    CHECK(d.wall_d1 <= 3e-13);
    CHECK(d.wall_d2 <= 2e-9);
}

TEST_CASE("the five standard levels match v1: grid, info and rebuilt exit", "[turbolenza][grid][v1ref]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    const ref::json& fam = R->family;
    const double beta1 = ref::num(fam.at("inlet_angle")), beta2 = ref::num(fam.at("exit_angle"));
    GridDiff all;
    double info_d = 0.0, rebuilt = 0.0;
    for (const char* level : T::kGridLevels) {
        INFO(level);
        const ref::json& v = fam.at("levels").at(level);
        const T::StandardGrid sg = T::build_grid(R->cascade, level, beta1, beta2);
        all.merge(compare(sg.grid, v, v.at("boundaries")));
        // info
        const ref::json& vi = v.at("info");
        const T::StandardGridInfo& i = sg.info;
        CHECK(i.level == level);
        CHECK(i.n_blade_cells == vi.at("n_blade_cells").get<int>());
        CHECK(i.n_m == vi.at("n_m").get<int>());
        CHECK(i.n_eta == vi.at("n_eta").get<int>());
        CHECK(i.total_nodes == vi.at("total_nodes").get<int>());
        CHECK(i.first_blade_i == vi.at("first_blade_i").get<int>());
        CHECK(i.last_blade_i == vi.at("last_blade_i").get<int>());
        CHECK(i.n_blade_stations == vi.at("n_blade_stations").get<int>());
        CHECK(i.admissibility == vi.at("admissibility").get<std::string>());
        CHECK(i.target_dm_over_chord == ref::num(vi.at("target_dm_over_chord")));
        // v1's checks, less its last ("metrics with positive S": the metrics are not the grid's)
        auto checks = vi.at("checks").get<std::vector<std::string>>();
        REQUIRE(checks.back() == "metrics with positive S");
        checks.pop_back();
        CHECK(i.checks == checks);
        const std::pair<double, const char*> numbers[] = {
            {i.dm, "dm"},
            {i.dm_over_chord, "dm_over_chord"},
            {i.deta, "deta"},
            {i.m_first, "m_first"},
            {i.m_last, "m_last"},
            {i.first_ss_x_c, "first_ss_x_c"},
            {i.first_ps_x_c, "first_ps_x_c"},
            {i.last_ss_x_c, "last_ss_x_c"},
            {i.last_ps_x_c, "last_ps_x_c"},
            {i.le_gap_over_chord, "le_gap_over_chord"},
            {i.te_gap_over_chord, "te_gap_over_chord"},
            {i.m_first_minus_m_le_over_chord, "m_first_minus_m_le_over_chord"},
            {i.m_te_minus_m_last_over_chord, "m_te_minus_m_last_over_chord"},
            {i.upstream_extension_chords, "upstream_extension_chords"},
            {i.downstream_extension_chords, "downstream_extension_chords"},
            {i.cell_aspect, "cell_aspect"}};
        for (const auto& [x, name] : numbers) info_d = std::max(info_d, std::abs(x - ref::num(vi.at(name))));
        const std::pair<T::Point, const char*> wall_points[] = {{i.first_ss_wall_point, "first_ss_wall_point"},
                                                                {i.first_ps_wall_point, "first_ps_wall_point"},
                                                                {i.last_ss_wall_point, "last_ss_wall_point"},
                                                                {i.last_ps_wall_point, "last_ps_wall_point"}};
        for (const auto& [p, name] : wall_points) {
            const auto& q = vi.at(name);  // a tuple in v1: a JSON list
            info_d = std::max({info_d, std::abs(p.m - ref::num(q[0])), std::abs(p.y - ref::num(q[1]))});
        }
        CHECK(sg.candidates.back().accepted);
        // the exit turned to 20 deg
        const T::Grid e = T::rebuild_exit(sg.grid, 20.0 * fixture::kPi / 180.0);
        const auto nodes = ref::vec(v.at("rebuilt_exit_20deg"));
        REQUIRE(nodes.size() == 2 * e.node_m.size());
        for (std::size_t k = 0; k < e.node_m.size(); ++k)
            rebuilt = std::max({rebuilt, std::abs(e.node_m.v[k] - nodes[2 * k]), std::abs(e.node_y.v[k] - nodes[2 * k + 1])});
    }
    INFO(text("standard levels vs v1", all));
    INFO("info " << info_d << ", rebuilt exit " << rebuilt);
    // measured (MSVC 19.51, five levels): m, eta, boundaries and info 0 (bit for bit), nodes 2.2e-16,
    // wall u 4.4e-16, contour 4.4e-16, wall_d relative 1.9e-16, 2.2e-13, 2.5e-11, rebuilt exit
    // nodes 2.2e-16.  Bounds about 10x the measurement.
    CHECK(all.m <= 1e-15);
    CHECK(all.eta == 0.0);
    CHECK(all.nodes <= 3e-15);
    CHECK(all.wall_u <= 5e-15);
    CHECK(all.wall_contour <= 5e-15);
    CHECK(all.boundaries <= 1e-15);
    CHECK(all.wall_d0 <= 2e-15);
    CHECK(all.wall_d1 <= 3e-12);
    CHECK(all.wall_d2 <= 3e-10);
    CHECK(info_d <= 1e-15);
    CHECK(rebuilt <= 3e-15);
}

TEST_CASE("the passage transitions match v1", "[turbolenza][grid][v1ref]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    const ref::json& v = R->family.at("transitions");
    const T::PassageTransitions t = T::passage_transitions(R->cascade);
    const double d = std::max({std::abs(t.m_tip - ref::num(v.at("m_tip"))), std::abs(t.m_le - ref::num(v.at("m_le"))),
                               std::abs(t.m_start - ref::num(v.at("m_start"))), std::abs(t.m_end - ref::num(v.at("m_end"))),
                               std::abs(t.m_te - ref::num(v.at("m_te"))), std::abs(t.m_tail - ref::num(v.at("m_tail")))});
    const auto& cov = v.at("coverage");
    REQUIRE(t.coverage.size() == cov.size());
    double dc = 0.0;
    for (std::size_t k = 0; k < cov.size(); ++k) {
        dc = std::max({dc, std::abs(t.coverage[k].lo - ref::num(cov[k][0])), std::abs(t.coverage[k].hi - ref::num(cov[k][1]))});
        CHECK(t.coverage[k].suction == cov[k][2].get<int>());
        CHECK(t.coverage[k].pressure == cov[k][3].get<int>());
    }
    INFO("transitions vs v1: " << d << ", coverage " << dc);
    // measured: 8.7e-19 (the nose tip; the landmarks exactly), coverage 1.1e-16
    CHECK(d <= 1e-17);
    CHECK(dc <= 1e-15);
}
