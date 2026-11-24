// MatrixAssembly (port of v1 matrix_assembly.py): F and Q of v1 pass 0 and pass 1 entry by entry
// (COO, as sorted triplets) from v1's coefficients and metrics; the numbering; uniform flow in a
// straight channel solves the assembled system; the refusals.
#include "SolverRef.hpp"

#include "turbolenza/MatrixAssembly.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <tuple>

using namespace BladeGen::Turbolenza;

namespace {

using Triplet = std::tuple<int, int, double>;

std::vector<Triplet> triplets(const SparseMatrix& F)
{
    std::vector<Triplet> t;
    for (int c = 0; c < F.outerSize(); ++c)
        for (SparseMatrix::InnerIterator it(F, c); it; ++it)
            t.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
    std::sort(t.begin(), t.end());
    return t;
}

std::vector<double> multiply(const SparseMatrix& F, const std::vector<double>& x)
{
    std::vector<double> y(static_cast<std::size_t>(F.rows()), 0.0);
    for (int c = 0; c < F.outerSize(); ++c)
        for (SparseMatrix::InnerIterator it(F, c); it; ++it)
            y[static_cast<std::size_t>(it.row())] += it.value() * x[static_cast<std::size_t>(c)];
    return y;
}

}  // namespace

TEST_CASE("matrix assembly: v1's F and Q entry by entry", "[turbolenza][assembly][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    for (int k = 0; k < 2; ++k) {
        const ref::json& p = d.at("passes").at(k);
        const auto g = sref::grid(p.at("grid"));
        const GridMetrics metrics = sref::metrics(p.at("metrics"), g);
        const BoundaryConditions bcs = build_boundary_conditions(g);
        const AssemblyResult a = assemble_from_coefficients(sref::coefficients(p.at("coefficients")), metrics, bcs);
        const ref::json& v1 = p.at("assembly");
        const ref::json& M = v1.at("matrix");
        std::vector<Triplet> ref_t;
        for (std::size_t q = 0; q < M.at("row").size(); ++q)
            ref_t.emplace_back(M.at("row")[q].get<int>(), M.at("col")[q].get<int>(), ref::num(M.at("val")[q]));
        std::sort(ref_t.begin(), ref_t.end());
        const std::vector<Triplet> t = triplets(a.matrix);
        const std::string tag = " pass " + std::to_string(k);
        REQUIRE(a.matrix.rows() == M.at("sparse_shape")[0].get<int>());
        REQUIRE(t.size() == ref_t.size());  // the same pattern, exact zeros dropped as v1
        double worst = 0.0, scale = 0.0;
        for (std::size_t q = 0; q < t.size(); ++q) {
            REQUIRE(std::get<0>(t[q]) == std::get<0>(ref_t[q]));
            REQUIRE(std::get<1>(t[q]) == std::get<1>(ref_t[q]));
            worst = std::max(worst, std::abs(std::get<2>(t[q]) - std::get<2>(ref_t[q])));
            scale = std::max(scale, std::abs(std::get<2>(ref_t[q])));
        }
        // measured 0 (bit for bit) from v1's own coefficients; bound 1e-15
        CHECK(sref::measure("assembly F entries vs v1" + tag, worst / scale) <= 1e-15);
        CHECK(sref::measure("assembly Q vs v1" + tag, ref::rel_diff(a.rhs, ref::vec(v1.at("rhs")))) <= 1e-15);
        const std::vector<double> active = ref::vec(v1.at("active_nodes"));
        REQUIRE(active.size() == a.active_nodes.size());
        for (std::size_t r = 0; r < active.size(); ++r) CHECK(static_cast<int>(active[r]) == a.active_nodes[r]);
        const ref::json& kinds = v1.at("row_kind");
        for (std::size_t r = 0; r < a.row_kind.size(); ++r)
            CHECK(kinds[r].get<std::string>() == row_kind_name(a.row_kind[r]));
        CHECK(v1.at("substitutions").get<int>() == a.substitutions);
        CHECK(v1.at("eliminations").size() == a.eliminations.size());
        CHECK(v1.at("implied").size() == a.implied.size());
    }
}

TEST_CASE("matrix assembly: uniform flow in a straight channel solves F Psi = Q", "[turbolenza][assembly]")
{
    for (const double angle : {0.0, 0.35}) {
        const auto g = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, angle, 0.8));
        const GridMetrics m = build_metrics(g);
        const BoundaryConditions b = build_boundary_conditions(g);
        const AssemblyResult a = assemble_from_coefficients(build_coefficients(m, Field(31, 9, 1.1), 7.0, 0.1), m, b);
        CHECK(a.n_unknowns() == (31 - 2) * 9);  // the far stations eliminated
        CHECK(a.implied.size() == 2);
        std::vector<double> x(static_cast<std::size_t>(a.n_unknowns()));
        for (std::size_t r = 0; r < x.size(); ++r) x[r] = g->eta[static_cast<std::size_t>(a.active_nodes[r] % 9)];
        const std::vector<double> Fx = multiply(a.matrix, x);
        double worst = 0.0;
        for (std::size_t r = 0; r < x.size(); ++r) worst = std::max(worst, std::abs(Fx[r] - a.rhs[r]));
        CHECK(worst <= 1e-12);
        const Field psi = reconstruct(a, x);
        for (int j = 0; j < 9; ++j) CHECK(std::abs(psi(0, j) - g->eta[static_cast<std::size_t>(j)]) <= 1e-15);
        for (double r : implied_residuals(a, psi)) CHECK(std::abs(r) <= 1e-15);
    }
    CHECK(full_index(3, 4, 10, 7) == 25);
    CHECK(node_of(25, 10, 7) == Node{3, 4});
    CHECK_THROWS_AS(full_index(-1, 0, 10, 7), MatrixAssemblyError);
    CHECK_THROWS_AS(node_of(70, 10, 7), MatrixAssemblyError);
}

TEST_CASE("matrix assembly: refusals", "[turbolenza][assembly]")
{
    const auto g = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, 0.2));
    const GridMetrics m = build_metrics(g);
    const BoundaryConditions b = build_boundary_conditions(g);
    const CoefficientFields c = build_coefficients(m, Field(31, 9, 1.0), 1.0, 0.1);
    const InteriorStencil s = build_interior_stencil(c, m);
    std::vector<PeriodicPDERow> rows = build_periodic_pde_rows(c, m, b);
    // a periodic row missing: its node has no equation and the jump has two free images
    std::vector<PeriodicPDERow> fewer(rows.begin() + 1, rows.end());
    CHECK_THROWS_AS(assemble_system(s, fewer, b), MatrixAssemblyError);
    // a periodic row twice: one node owns two equations
    std::vector<PeriodicPDERow> twice = rows;
    twice.push_back(rows.front());
    CHECK_THROWS_AS(assemble_system(s, twice, b), MatrixAssemblyError);
    // a stencil of another grid
    const InteriorStencil other = interior_stencil(Field(30, 9, 1.0), Field(30, 9, 0.0), Field(30, 9, 1.0),
                                                   Field(30, 9, 0.0), Field(30, 9, 0.0), 0.1, 0.1);
    CHECK_THROWS_AS(assemble_system(other, rows, b), MatrixAssemblyError);
    // a periodic row whose target is not on a periodic line
    std::vector<PeriodicPDERow> moved = rows;
    moved.front().target = Node{10, 8};
    CHECK_THROWS_AS(assemble_system(s, moved, b), MatrixAssemblyError);
}
