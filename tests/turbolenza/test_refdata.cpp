// The v1 reference data reads back into the port's grid contract (GridData.hpp).
#include "RefData.hpp"

#include "turbolenza/Case.hpp"  // the contract headers compile on their own

#include <catch2/catch_test_macros.hpp>

TEST_CASE("v1 reference grid loads into turbolenza::Grid", "[turbolenza][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(fam, "grid_family");
    const auto g = ref::grid(fam.at("levels").at("very_coarse"));
    REQUIRE(g.n_m() == static_cast<int>(g.node_m.n_m));
    REQUIRE(g.n_eta() == g.node_m.n_eta);
    CHECK(g.first_blade_i >= 1);
    CHECK(g.last_blade_i < g.n_m() - 1);
    CHECK(g.wall_d[0].lower.size() == static_cast<std::size_t>(g.last_blade_i - g.first_blade_i + 1));
    CHECK(g.wall_identity.lower.surface.size() == g.wall_d[0].lower.size());
    for (int i = 0; i < g.n_m(); ++i) CHECK(g.node_m(i, 0) == g.m[static_cast<std::size_t>(i)]);
    // wall node k = 0 of the curve derivative is the node itself (v1 Grid._check_on_curves)
    for (int b = 0; b <= g.last_blade_i - g.first_blade_i; ++b)
        CHECK(std::abs(g.wall_d[0].lower[static_cast<std::size_t>(b)].y - g.node_y(g.first_blade_i + b, 0)) <= 1e-12);
}
