// BoundaryConditions (port of v1 boundary_conditions.py): the classification and relations of
// the v1 grid, the eliminations 4/3, -1/3, the corner identity, the ghost rule, the refusals.
#include "SolverRef.hpp"

#include "turbolenza/BoundaryConditions.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

using namespace BladeGen::Turbolenza;

TEST_CASE("boundary conditions: the v1 relations of the very_coarse grid", "[turbolenza][boundary][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    const ref::json& p = d.at("passes").at(0);
    const BoundaryConditions b = build_boundary_conditions(sref::grid(p.at("grid")));
    const ref::json& v1 = p.at("boundary_conditions");
    const ref::json& rel = v1.at("relations");
    REQUIRE(rel.size() == b.relations.size());
    for (std::size_t k = 0; k < rel.size(); ++k) {
        const Relation& r = b.relations[k];
        CHECK(rel[k].at("kind").get<std::string>() == relation_kind_name(r.kind));
        CHECK(ref::num(rel[k].at("rhs")) == r.rhs);
        REQUIRE(rel[k].at("terms").size() == r.terms.size());
        for (std::size_t t = 0; t < r.terms.size(); ++t) {
            CHECK(rel[k].at("terms")[t].at("node")[0].get<int>() == r.terms[t].node.i);
            CHECK(rel[k].at("terms")[t].at("node")[1].get<int>() == r.terms[t].node.j);
            CHECK(ref::num(rel[k].at("terms")[t].at("coefficient")) == r.terms[t].coefficient);
        }
    }
    const std::pair<const char*, const std::vector<char>*> masks[] = {
        {"lower_wall", &b.lower_wall}, {"upper_wall", &b.upper_wall}, {"periodic_lower", &b.periodic_lower},
        {"periodic_upper", &b.periodic_upper}, {"upstream", &b.upstream}, {"downstream", &b.downstream}};
    for (const auto& [name, mask] : masks) {
        const ref::json& data = v1.at(name).at("data");
        REQUIRE(data.size() == mask->size());
        for (std::size_t q = 0; q < mask->size(); ++q) CHECK(data[q].get<bool>() == ((*mask)[q] != 0));
    }
}

TEST_CASE("boundary conditions: eliminations, corners and the ghost rule", "[turbolenza][boundary]")
{
    const BoundaryConditions b = build_boundary_conditions(std::make_shared<const Grid>(sref::channel(31, 7, 8, 20, 0.2)));
    const Elimination u = b.upstream_elimination(3);
    REQUIRE(u.terms.size() == 2);
    CHECK(u.terms[0].coefficient == 4.0 / 3.0);
    CHECK(u.terms[1].coefficient == -1.0 / 3.0);
    CHECK(u.constant == 0.0);
    const Elimination dn = b.downstream_elimination(6);
    CHECK(dn.node == Node{30, 6});
    CHECK(dn.terms[0].node == Node{29, 6});
    // the corner identity: U(0, J) - U(0, 0) = -3 P(0) + 4 P(1) - P(2) for any psi
    Field psi(31, 7);
    for (std::size_t q = 0; q < psi.size(); ++q) psi.v[q] = std::sin(0.37 * static_cast<double>(q));
    auto P = [&](int i) { return b.periodic_relation(i).residual(psi); };
    double U0 = 0.0, UJ = 0.0;
    for (const Relation& r : b.relations)
        if (r.kind == RelationKind::Upstream) {
            if (*r.target == Node{0, 0}) U0 = r.residual(psi);
            if (*r.target == Node{0, 6}) UJ = r.residual(psi);
        }
    CHECK(std::abs((UJ - U0) - (-3.0 * P(0) + 4.0 * P(1) - P(2))) <= 1e-14);
    // the ghosts: psi[i, -1] = psi[i, J-1] - 1, psi[i, n_eta] = psi[i, 1] + 1
    CHECK(b.periodic_ghost(4, -1).value(psi) == psi(4, 5) - 1.0);
    CHECK(b.periodic_ghost(4, 7).value(psi) == psi(4, 1) + 1.0);
    CHECK_THROWS_AS(b.periodic_ghost(10, 7), BoundaryConditionError);  // a blade station
    CHECK_THROWS_AS(b.periodic_ghost(4, 3), BoundaryConditionError);
    CHECK(b.interior(5, 3));
    CHECK(!b.interior(0, 3));
    CHECK(b.periodic_line(5, 6));
    CHECK(!b.periodic_line(0, 6));
}

TEST_CASE("boundary conditions: refusals", "[turbolenza][boundary]")
{
    CHECK_THROWS_AS(build_boundary_conditions(std::make_shared<const Grid>(sref::channel(31, 2, 8, 20, 0.2))),
                    BoundaryConditionError);
    CHECK_THROWS_AS(build_boundary_conditions(std::make_shared<const Grid>(sref::channel(31, 7, 2, 20, 0.2))),
                    BoundaryConditionError);
    CHECK_THROWS_AS(build_boundary_conditions(std::make_shared<const Grid>(sref::channel(31, 7, 8, 28, 0.2))),
                    BoundaryConditionError);
    Grid g = sref::channel(31, 7, 8, 20, 0.2);
    g.wall_identity.stations[2] += 1;
    CHECK_THROWS_AS(build_boundary_conditions(std::make_shared<const Grid>(g)), BoundaryConditionError);
    CHECK_THROWS_AS(Relation(RelationKind::Periodic, {Term{{1, 1}, 1.0}, Term{{1, 1}, -1.0}}, 1.0, std::nullopt),
                    BoundaryConditionError);
    CHECK_THROWS_AS(Relation(RelationKind::LowerWall, {Term{{1, 1}, 0.0}}, 0.0, Node{1, 1}), BoundaryConditionError);
}
