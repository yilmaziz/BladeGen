// LinearSolver (port of v1 linear_solver.py): psi of v1 pass 0 / 1 from the assembled system with
// SparseLU and PARDISO (1 and 4 threads); a hand system with a known solution; singular and
// non-finite systems refused.
#include "SolverRef.hpp"

#include "turbolenza/LinearSolver.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <vector>

using namespace BladeGen::Turbolenza;

namespace {

SparseMatrix dense(const std::vector<std::vector<double>>& a)
{
    std::vector<Eigen::Triplet<double, int>> t;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[i].size(); ++j)
            if (a[i][j] != 0.0) t.emplace_back(static_cast<int>(i), static_cast<int>(j), a[i][j]);
    SparseMatrix F(static_cast<int>(a.size()), static_cast<int>(a.size()));
    F.setFromTriplets(t.begin(), t.end());
    return F;
}

std::vector<LinearSolverKind> kinds()
{
    std::vector<LinearSolverKind> k{LinearSolverKind::SparseLU};
    if (pardiso_available()) k.push_back(LinearSolverKind::Pardiso);
    return k;
}

}  // namespace

TEST_CASE("linear solver: psi of v1 passes from the assembled system", "[turbolenza][linear][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    for (int k = 0; k < 2; ++k) {
        const ref::json& p = d.at("passes").at(k);
        const auto g = sref::grid(p.at("grid"));
        const AssemblyResult a = assemble_from_coefficients(sref::coefficients(p.at("coefficients")),
                                                            sref::metrics(p.at("metrics"), g),
                                                            build_boundary_conditions(g));
        const Field v1 = ref::field(p.at("psi"));
        const std::string tag = " pass " + std::to_string(k);
        const LinearSolveResult lu = solve_linear_system(a, LinearSolverKind::SparseLU);
        // SuperLU (v1) against Eigen SparseLU: measured 1.1e-14 / 5.1e-15; bound 1e-13
        CHECK(sref::measure("linear psi SparseLU vs v1" + tag, ref::rel_diff(lu.psi, v1)) <= 1e-13);
        CHECK(sref::measure("linear residual SparseLU" + tag, lu.residual_inf) <= 1e-14);
        CHECK(lu.implied_residual <= 1e-14);
        if (pardiso_available()) {
            for (int threads : {1, 4}) {
                const LinearSolveResult pd = solve_linear_system(a, LinearSolverKind::Pardiso, threads);
                const std::string t2 = tag + " threads " + std::to_string(threads);
                CHECK(sref::measure("linear psi PARDISO vs v1" + t2, ref::rel_diff(pd.psi, v1)) <= 1e-13);
                CHECK(sref::measure("linear psi PARDISO vs SparseLU" + t2, ref::rel_diff(pd.psi, lu.psi)) <= 1e-13);
                CHECK(sref::measure("linear residual PARDISO" + t2, pd.residual_inf) <= 1e-13);
            }
        }
    }
}

TEST_CASE("linear solver: a hand system, and the refusals", "[turbolenza][linear]")
{
    // F x = Q with x = (1, -2, 3) and a zero diagonal (the TE row's centre can be near zero: no row is scaled)
    const SparseMatrix F = dense({{0.0, 2.0, 1.0}, {1.0, 1.0, 0.0}, {4.0, 0.0, -1.0}});
    const std::vector<double> Q{-1.0, -1.0, 1.0};
    for (LinearSolverKind kind : kinds()) {
        const std::vector<double> x = solve(F, Q, kind);
        CHECK(std::abs(x[0] - 1.0) <= 1e-15);
        CHECK(std::abs(x[1] + 2.0) <= 1e-15);
        CHECK(std::abs(x[2] - 3.0) <= 1e-15);
        CHECK(residual_norm(F, x, Q) <= 1e-15);
        CHECK_THROWS_AS(solve(dense({{1.0, 1.0}, {1.0, 1.0}}), {1.0, 2.0}, kind), LinearSolveError);
        CHECK_THROWS_AS(solve(dense({{1.0, 0.0}, {0.0, NAN}}), {1.0, 2.0}, kind), LinearSolveError);
        CHECK_THROWS_AS(solve(F, {1.0, 2.0}, kind), LinearSolveError);
        CHECK_THROWS_AS(solve(F, {1.0, INFINITY, 0.0}, kind), LinearSolveError);
        CHECK_THROWS_AS(solve(SparseMatrix(0, 0), {}, kind), LinearSolveError);
    }
    if (!pardiso_available()) CHECK_THROWS_AS(solve(F, Q, LinearSolverKind::Pardiso), LinearSolveError);
}

TEST_CASE("linear solver: timing of the backends on the largest dumped system", "[turbolenza][linear][v1ref][timing]")
{
    TURBOLENZA_REF_OR_SKIP(d, "fixed_exit_explicit_45x21");
    const ref::json& p = d.at("passes").at(0);
    const auto g = sref::grid(p.at("grid"));
    const AssemblyResult a = assemble_from_coefficients(sref::coefficients(p.at("coefficients")),
                                                        sref::metrics(p.at("metrics"), g), build_boundary_conditions(g));
    const Field v1 = ref::field(p.at("psi"));
    auto time = [&](LinearSolverKind kind, int threads, const char* name) {
        LinearSolveResult r;
        const int repeat = 20;
        const auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < repeat; ++k) r = solve_linear_system(a, kind, threads);
        const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / repeat;
        sref::measure(std::string("linear solve time [s] ") + name, s);
        CHECK(sref::measure(std::string("linear psi vs v1 45x21 ") + name, ref::rel_diff(r.psi, v1)) <= 1e-13);
    };
    std::printf("MEASURE system 45x21: %d unknowns, %lld non-zeros\n", a.n_unknowns(),
                static_cast<long long>(a.matrix.nonZeros()));
    time(LinearSolverKind::SparseLU, 1, "SparseLU");
    if (pardiso_available()) {
        time(LinearSolverKind::Pardiso, 1, "PARDISO 1 thread");
        time(LinearSolverKind::Pardiso, 4, "PARDISO 4 threads");
    }
}
