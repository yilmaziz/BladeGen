// Metrics (port of v1 metrics.py): v1's metrics of the very_coarse grid of pass 0 and of the
// rebuilt grid of pass 1, and the closed forms of a straight sheared channel; the refusals.
#include "SolverRef.hpp"

#include "turbolenza/Metrics.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>

using namespace BladeGen::Turbolenza;

TEST_CASE("metrics: v1 fields of the dumped grids", "[turbolenza][metrics][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    for (int k = 0; k < 2; ++k) {
        const ref::json& p = d.at("passes").at(k);
        const auto g = sref::grid(p.at("grid"));
        const GridMetrics m = build_metrics(g);
        const GridMetrics v1 = sref::metrics(p.at("metrics"), g);
        const std::string tag = " pass " + std::to_string(k);
        // measured: S, tan_beta0/1, dS_dm, tan_beta, dtan_beta_dm 0 (bit for bit); cos_beta,
        // sin_beta, beta <= 1.1e-16 (hypot / arctan of MSVC vs NumPy); bound 1e-15
        CHECK(sref::measure("metrics S" + tag, ref::rel_diff(m.S, v1.S)) <= 1e-15);
        CHECK(sref::measure("metrics tan_beta0" + tag, ref::rel_diff(m.tan_beta0, v1.tan_beta0)) <= 1e-15);
        CHECK(sref::measure("metrics tan_beta1" + tag, ref::rel_diff(m.tan_beta1, v1.tan_beta1)) <= 1e-15);
        CHECK(sref::measure("metrics dS_dm" + tag, ref::rel_diff(m.dS_dm, v1.dS_dm)) <= 1e-15);
        CHECK(sref::measure("metrics dtan_beta0_dm" + tag, ref::rel_diff(m.dtan_beta0_dm, v1.dtan_beta0_dm)) <= 1e-15);
        CHECK(sref::measure("metrics dtan_beta1_dm" + tag, ref::rel_diff(m.dtan_beta1_dm, v1.dtan_beta1_dm)) <= 1e-15);
        CHECK(sref::measure("metrics tan_beta" + tag, ref::rel_diff(m.tan_beta, v1.tan_beta)) <= 1e-15);
        CHECK(sref::measure("metrics cos_beta" + tag, ref::rel_diff(m.cos_beta, v1.cos_beta)) <= 1e-15);
        CHECK(sref::measure("metrics sin_beta" + tag, ref::rel_diff(m.sin_beta, v1.sin_beta)) <= 1e-15);
        CHECK(sref::measure("metrics beta" + tag, ref::rel_diff(m.beta, v1.beta)) <= 1e-15);
        CHECK(sref::measure("metrics dtan_beta_dm" + tag, ref::rel_diff(m.dtan_beta_dm, v1.dtan_beta_dm)) <= 1e-15);
    }
}

TEST_CASE("metrics: a straight sheared channel has its exact slopes", "[turbolenza][metrics]")
{
    const double angle = 0.4;
    const auto g = std::make_shared<const Grid>(sref::channel(31, 7, 8, 20, angle, 0.8));
    const GridMetrics m = build_metrics(g);
    const double t = std::tan(angle);
    for (int i = 0; i < g->n_m(); ++i) {
        CHECK(std::abs(m.S[static_cast<std::size_t>(i)] - 0.8) <= 1e-14);
        // the side-line slopes are differences of a linear function: exact to round-off
        CHECK(std::abs(m.tan_beta0[static_cast<std::size_t>(i)] - t) <= 1e-12);
        CHECK(std::abs(m.tan_beta1[static_cast<std::size_t>(i)] - t) <= 1e-12);
        CHECK(std::abs(m.dtan_beta0_dm[static_cast<std::size_t>(i)]) <= 1e-10);
        for (int j = 0; j < g->n_eta(); ++j) {
            CHECK(std::abs(m.cos_beta(i, j) - std::cos(angle)) <= 1e-13);
            CHECK(std::abs(m.beta(i, j) - angle) <= 1e-12);
        }
    }
    // the blade walls take the curve's own slope y_u / m_u: exactly tan(angle) * L / L
    CHECK(m.tan_beta0[10] == (t * (g->m_te() - g->m_le())) / (g->m_te() - g->m_le()));
}

TEST_CASE("metrics: refusals", "[turbolenza][metrics]")
{
    // a run of fewer than three stations
    CHECK_THROWS_AS(build_metrics(std::make_shared<const Grid>(sref::channel(31, 7, 2, 20, 0.2))), MetricsError);
    CHECK_THROWS_AS(build_metrics(std::make_shared<const Grid>(sref::channel(31, 7, 8, 28, 0.2))), MetricsError);
    // a node off the linear eta map
    Grid g = sref::channel(31, 7, 8, 20, 0.2);
    g.node_y(4, 3) += 1e-6;
    CHECK_THROWS_AS(build_metrics(std::make_shared<const Grid>(g)), MetricsError);
    // non-uniform stations
    Grid h = sref::channel(31, 7, 8, 20, 0.2);
    h.m[5] += 1e-4;
    for (int j = 0; j < 7; ++j) h.node_m(5, j) = h.m[5];
    CHECK_THROWS_AS(build_metrics(std::make_shared<const Grid>(h)), MetricsError);
    // a wall curve point that is not its node
    Grid w = sref::channel(31, 7, 8, 20, 0.2);
    w.wall_d[0].lower[3].y += 1e-9;
    CHECK_THROWS_AS(build_metrics(std::make_shared<const Grid>(w)), MetricsError);
    // a numerically vertical wall (m_u = 0)
    Grid v = sref::channel(31, 7, 8, 20, 0.2);
    v.wall_d[1].upper[4] = {0.0, 1.0};
    CHECK_THROWS_AS(build_metrics(std::make_shared<const Grid>(v)), MetricsError);
    // S = 0: a collapsed passage
    Grid s = sref::channel(31, 7, 8, 20, 0.0);
    for (int j = 0; j < 7; ++j) s.node_y(12, j) = 0.0;
    s.wall_d[0].upper[4].y = 0.0;
    CHECK_THROWS_AS(build_metrics(std::make_shared<const Grid>(s)), MetricsError);
}
