// The general operator on TURBOlenza's own cascade (verification stages 1 and 3 of the plan): the tensor
// K = mdot' kappa J g^ab against Aungier's A, B, C on the standard grids, AC - B^2 = mu^2, and the solve of
// the reference NACA 65-(12A10)10 case (v1's blade, solidity 1, stagger 33 deg, beta1 45 deg) by the
// general operator against the Aungier solver on the same nodes, at M1 0.3 and 0.65, levels very_coarse
// .. medium: psi, surface Mach, beta_out.  Then the same passage on a curvilinear grid (curved stations),
// solved by the general operator and compared with its solution on Aungier's grid (one passage, two meshes:
// the differences shrink with refinement).  Against Aungier the differences are bounded, not shown to
// shrink: on this grid both discretizations are of low order next to the blade edges (the transition
// cells, whose blade thickness the round NACA edges make ~ sqrt(dm); see test_general_operator.cpp, the
// H-grid with a blade, where both converge to the exact solution with different error constants), and
// both discharge angles still move 0.6 .. 1.1 deg per level.  [needs-grid]: v1's blade.
#include "CascadeFixture.hpp"
#include "GeneralFixture.hpp"

#include "turbolenza/Coefficients.hpp"
#include "turbolenza/Detail.hpp"
#include "turbolenza/GeneralAdapter.hpp"
#include "turbolenza/GeneralSolver.hpp"
#include "turbolenza/Grid.hpp"
#include "turbolenza/Metrics.hpp"
#include "turbolenza/Solver.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <format>
#include <memory>

using namespace BladeGen::Turbolenza;
using gfix::measure;

namespace {

std::size_t z(int k) { return static_cast<std::size_t>(k); }

constexpr double kChord = 0.5;  // c_ref [m] of the user's case

SolverConfig aungier_config(double mach, double tol)
{
    SolverConfig c;
    c.flow.inlet_mach = mach;
    c.flow.inlet_angle_deg = 45.0;
    c.iteration.density_tolerance = tol;
    c.iteration.kutta_tolerance = tol;
    c.iteration.max_density_iterations = 4000;
    c.iteration.max_kutta_iterations = 4000;
    c.iteration.density_relaxation = 0.5;
    c.iteration.kutta_relaxation = 0.5;
    return c;
}

GeneralConfig general_config(const SolverConfig& a, double exit_angle, bool kutta)
{
    GeneralConfig g;
    g.flow = a.flow;
    g.flow.exit_angle_guess_deg = detail::degrees(exit_angle);
    g.gas = a.gas;
    g.iteration = a.iteration;
    g.kutta = kutta;
    return g;
}

struct Diff
{
    double psi = 0.0, mach_surface = 0.0, beta_deg = 0.0;
    double mach_surface_inner = 0.0;  ///< blade stations 0.1 c or more from the leading and trailing edge
    int at_i = 0, at_j = 0;           ///< where mach_surface is reached
};

Diff compare(const Solution& a, const GeneralSolution& g)
{
    const Grid& grid = *a.grid;
    Diff d;
    for (int i = 0; i < grid.n_m(); ++i)
        for (int j = 0; j < grid.n_eta(); ++j)
            d.psi = std::max(d.psi, std::abs(a.psi(i, j) - g.psi[z(general_node_of(grid, i, j))]));
    for (int i = grid.first_blade_i; i <= grid.last_blade_i; ++i)
        for (int j : {0, grid.n_eta() - 1}) {
            const double e = std::abs(a.mach(i, j) - g.mach[z(general_node_of(grid, i, j))]);
            if (e > d.mach_surface) {
                d.at_i = i;
                d.at_j = j;
            }
            d.mach_surface = std::max(d.mach_surface, e);
            if (grid.m[z(i)] >= grid.m_le() + 0.1 && grid.m[z(i)] <= grid.m_te() - 0.1) d.mach_surface_inner = std::max(d.mach_surface_inner, e);
        }
    d.beta_deg = detail::degrees(g.exit_angle - a.exit_angle);
    return d;
}


/// Two general solutions on meshes with the same wall nodes (Aungier's grid and its curved version).
Diff compare_walls(const Grid& grid, const GeneralSolution& a, const GeneralSolution& b)
{
    Diff d;
    for (int i = grid.first_blade_i; i <= grid.last_blade_i; ++i)
        for (int j : {0, grid.n_eta() - 1}) {
            const double e = std::abs(a.mach[z(general_node_of(grid, i, j))] - b.mach[z(general_node_of(grid, i, j))]);
            d.mach_surface = std::max(d.mach_surface, e);
            if (grid.m[z(i)] >= grid.m_le() + 0.1 && grid.m[z(i)] <= grid.m_te() - 0.1) d.mach_surface_inner = std::max(d.mach_surface_inner, e);
        }
    d.beta_deg = detail::degrees(b.exit_angle - a.exit_angle);
    return d;
}

/// The same passage on a curvilinear grid: every interior node slides along its own eta line (Catmull-Rom
/// in i at fixed j) by dm' = A sin(pi eta) sin(3 pi s), s = (m - m_in) / (m_out - m_in): curved stations,
/// the walls, the periodic lines, the inlet and the outlet unchanged (A in chords, independent of the level).
GeneralMesh curved_stations(const Grid& grid, double A)
{
    GeneralMesh mesh = mesh_from_grid(grid, kChord);
    MeshBlock& B = mesh.blocks[0];
    const int nm = grid.n_m(), ne = grid.n_eta();
    const double m_in = grid.m.front(), m_out = grid.m.back(), dm = grid.dm();
    for (int j = 1; j < ne - 1; ++j)
        for (int i = 1; i < nm - 1; ++i) {
            const double eta = static_cast<double>(j) / (ne - 1);
            const double s = (grid.m[z(i)] - m_in) / (m_out - m_in);
            const double target = grid.m[z(i)] + A * std::sin(gfix::kPi * eta) * std::sin(3.0 * gfix::kPi * s);
            const int k = std::clamp(static_cast<int>(std::floor((target - m_in) / dm)), 0, nm - 2);
            const double t = (target - grid.m[z(k)]) / dm;
            auto Y = [&](int q) { return grid.node_y(std::clamp(q, 0, nm - 1), j); };
            const double p0 = Y(k - 1), p1 = Y(k), p2 = Y(k + 1), p3 = Y(k + 2);
            const double y = 0.5 * ((2.0 * p1) + (-p0 + p2) * t + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t * t +
                                    (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t * t * t);
            const std::size_t q = B.at(i, j);
            B.x[q] = target;
            B.y[q] = y;
            B.m[q] = kChord * target;
        }
    validate_mesh(mesh);
    return mesh;
}

}  // namespace

TEST_CASE("general cascade: K = mdot' kappa J g reproduces Aungier's A, B, C", "[turbolenza][general][needs-grid]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    double previous_fd = 0.0;
    for (const char* level : {"very_coarse", "coarse", "medium", "fine"}) {
        const auto grid = std::make_shared<const Grid>(build_grid(R->cascade, level, detail::radians(45.0), std::nullopt).grid);
        const GridMetrics metrics = build_metrics(grid);
        Field rho(grid->n_m(), grid->n_eta());
        for (int i = 0; i < grid->n_m(); ++i)
            for (int j = 0; j < grid->n_eta(); ++j) rho(i, j) = 1.1 + 0.05 * std::sin(0.3 * i);
        const double mdot = 2.5;
        const CoefficientFields cf = build_coefficients(metrics, rho, mdot, kChord);
        const GeneralOperator op(mesh_from_grid(*grid, kChord));
        const double dm = grid->dm(), de = grid->deta();
        double exact = 0.0, fd_inner = 0.0, fd_wall = 0.0, fd_edge = 0.0, fd_ext = 0.0, det = 0.0;
        // within 0.1 chord of a blade end the grid lines follow the nose and the tail (slopes that vary on
        // the scale of the edge radius) and kink where the side lines meet the walls: a difference of the
        // nodes cannot match the pointwise slopes Aungier uses there (reported, not bounded)
        int worst_i = 0, worst_j = 0;
        auto near_edge = [&](int i) {
            const double m = grid->m[z(i)];
            return std::abs(m - grid->m_le()) < 0.1 || std::abs(m - grid->m_te()) < 0.1;
        };
        for (int i = 0; i < grid->n_m(); ++i)
            for (int j = 0; j < grid->n_eta(); ++j) {
                const double mu = mdot / rho(i, j), tb = metrics.tan_beta(i, j), S = metrics.S[z(i)];
                // Aungier's flux form times dm de in index space: [[C de / (L dm), -B], [-B, A L dm / de]]
                const double A = cf.A(i, j), B = cf.B(i, j), C = cf.C(i, j);
                const std::array<double, 3> want{C * de / (kChord * dm), -B, A * kChord * dm / de};
                const auto K = contravariant_tensor(dm, 0.0, tb * dm, S * de, mu);
                const NodalMetric& mt = op.metric(general_node_of(*grid, i, j));
                const auto Kfd = contravariant_tensor(mt.x_z, mt.x_e, mt.y_z, mt.y_e, mu);
                const double scale = std::max(std::abs(want[0]), std::abs(want[2]));
                for (int k = 0; k < 3; ++k) {
                    exact = std::max(exact, std::abs(K[z(k)] - want[z(k)]) / scale);
                    const double e = std::abs(Kfd[z(k)] - want[z(k)]) / scale;
                    const bool extension = i < grid->first_blade_i || i > grid->last_blade_i;
                    double& slot = near_edge(i) ? fd_edge : extension ? fd_ext : (j == 0 || j == grid->n_eta() - 1 ? fd_wall : fd_inner);
                    if (&slot == &fd_inner && e > fd_inner) {
                        worst_i = i;
                        worst_j = j;
                    }
                    slot = std::max(slot, e);
                }
                det = std::max(det, std::abs(A * C - B * B - mu * mu) / (mu * mu));
            }
        CHECK(measure(std::format("general K(Aungier metrics) vs A,B,C {}", level), exact) <= 1e-14);  // measured ~2e-16
        CHECK(measure(std::format("general AC - B^2 = mu^2 {}", level), det) <= 1e-13);              // measured ~1e-15
        measure(std::format("general K(node differences) vs A,B,C blade passage off the walls, > 0.1 c from the blade ends {}", level), fd_inner);
        measure(std::format("general K(node differences) vs A,B,C blade walls, > 0.1 c from the blade ends {}", level), fd_wall);
        measure(std::format("general K(node differences) vs A,B,C within 0.1 c of a blade end {}", level), fd_edge);
        measure(std::format("general K(node differences) vs A,B,C on the extensions (> 0.1 c from the blade) {}", level), fd_ext);
        CHECK(fd_ext <= 1e-12);  // measured 9e-15 .. 3e-14: the side lines' nodes give the slopes to round-off
        std::printf("  worst off-wall node (%d, %d), m - m_le = %.4f c, of %d x %d (blade %d .. %d)\n", worst_i, worst_j,
                    grid->m[z(worst_i)] - grid->m_le(), grid->n_m(), grid->n_eta(), grid->first_blade_i, grid->last_blade_i);
        if (previous_fd > 0.0) CHECK(fd_inner < previous_fd);
        previous_fd = fd_inner;
    }
}

TEST_CASE("general cascade: the general operator against the Aungier solver on the same nodes",
          "[turbolenza][general][needs-grid][slow]")
{
    TURBOLENZA_CASCADE_OR_SKIP(R);
    for (double mach : {0.3, 0.65}) {
        Diff first_fixed, last_fixed, first_grids, last_grids;
        int n = 0;
        for (const char* level : {"very_coarse", "coarse", "medium"}) {
            const SolverConfig ac = aungier_config(mach, 1e-10);
            const auto grid0 = std::make_shared<const Grid>(build_grid(R->cascade, level, detail::radians(45.0), std::nullopt).grid);
            // 1. fixed discharge angle (the grid's: the TE blade angle), the same grid
            const auto t0 = std::chrono::steady_clock::now();
            const Solution af = solve_at_exit_angle(grid0, ac, kChord);
            const auto t1 = std::chrono::steady_clock::now();
            const GeneralSolution gf = solve_general(mesh_from_grid(*grid0, kChord), general_config(ac, grid0->exit_angle(), false));
            const auto t2 = std::chrono::steady_clock::now();
            REQUIRE(af.converged);
            REQUIRE(gf.converged);
            CHECK(gf.mass_flow_per_depth == af.mass_flow_per_depth);
            const Diff df = compare(af, gf);
            const std::string tag = std::format("M1 {} {} ({} x {})", mach, level, grid0->n_m(), grid0->n_eta());
            measure("general vs Aungier fixed angle max |dpsi| " + tag, df.psi);
            // psi at the inlet corner (0, 0) of the periodic line: where the two differ most upstream
            measure("Aungier fixed angle psi(0, 0) " + tag, af.psi(0, 0));
            measure("general fixed angle psi(0, 0) " + tag, gf.psi[z(general_node_of(*grid0, 0, 0))]);
            measure("general vs Aungier fixed angle max |dM| on the walls " + tag, df.mach_surface);
            measure("general vs Aungier fixed angle max |dM| on the walls, 0.1 c .. 0.9 c " + tag, df.mach_surface_inner);
            std::printf("  largest wall difference at station %d (blade %d .. %d), j = %d\n", df.at_i, grid0->first_blade_i,
                        grid0->last_blade_i, df.at_j);
            std::printf("  fixed angle %s: Aungier %d passes %.2f s, general %d passes %.2f s, outlet mass-flow error %.2e\n",
                        tag.c_str(), af.iterations, std::chrono::duration<double>(t1 - t0).count(), gf.iterations,
                        std::chrono::duration<double>(t2 - t1).count(), gf.passes.back().outlet_mass_flow_error);
            // 2. the Kutta-coupled solve: Aungier's own loop (exit rebuilt), then the general Kutta loop on
            //    Aungier's final grid (outlet angle only), started at Aungier's beta_out
            const Solution ak = solve_coupled(grid0, ac, kChord);
            REQUIRE(ak.converged);
            const GeneralSolution gk = solve_general(mesh_from_grid(*ak.grid, kChord), general_config(ac, ak.exit_angle, true));
            REQUIRE(gk.converged);
            const Diff dk = compare(ak, gk);
            measure("general vs Aungier Kutta beta_out difference [deg] " + tag, dk.beta_deg);
            measure("general vs Aungier Kutta max |dpsi| " + tag, dk.psi);
            measure("general vs Aungier Kutta max |dM| on the walls " + tag, dk.mach_surface);
            measure("general vs Aungier Kutta max |dM| on the walls, 0.1 c .. 0.9 c " + tag, dk.mach_surface_inner);
            std::printf("  Kutta %s: beta_out Aungier %.6f deg, general %.6f deg (%d passes, eps_K %.1e)\n", tag.c_str(),
                        detail::degrees(ak.exit_angle), detail::degrees(gk.exit_angle), gk.iterations, gk.kutta->normalized);
            // 3. the same passage on a curvilinear grid (curved stations): the general Kutta solve against
            //    Aungier's on his grid, at the wall nodes (the same nodes) and in beta_out
            const GeneralSolution gc = solve_general(curved_stations(*ak.grid, 0.02), general_config(ac, ak.exit_angle, true));
            REQUIRE(gc.converged);
            const Diff dc = compare(ak, gc);
            measure("general (curved stations) vs Aungier Kutta beta_out difference [deg] " + tag, dc.beta_deg);
            measure("general (curved stations) vs Aungier Kutta max |dM| on the walls " + tag, dc.mach_surface);
            measure("general (curved stations) vs Aungier Kutta max |dM| on the walls, 0.1 c .. 0.9 c " + tag, dc.mach_surface_inner);
            const Diff dg = compare_walls(*ak.grid, gk, gc);
            measure("general curved vs general on Aungier's grid beta_out difference [deg] " + tag, dg.beta_deg);
            measure("general curved vs general on Aungier's grid max |dM| on the walls " + tag, dg.mach_surface);
            measure("general curved vs general on Aungier's grid max |dM| on the walls, 0.1 c .. 0.9 c " + tag, dg.mach_surface_inner);
            // regression bounds against Aungier, about 1.5 x measured (at most 1.05 deg, 0.027 in psi,
            // 0.032 in M at 0.1 .. 0.9 c on these levels): both H-grid discretizations are of low order at
            // the edge transition cells (test_general_operator.cpp), with different error constants
            CHECK(std::abs(dk.beta_deg) <= 1.6);
            CHECK(dk.psi <= 0.05);
            CHECK(dk.mach_surface_inner <= 0.08);
            if (n == 0) {
                first_fixed = df;
                first_grids = dg;
            }
            last_fixed = df;
            last_grids = dg;
            ++n;
        }
        // the general operator on two meshes of one passage (Aungier's grid, curved stations) converges to
        // one solution: beta_out and the wall Mach numbers of the two approach each other
        CHECK(std::abs(last_grids.beta_deg) < std::abs(first_grids.beta_deg));
        CHECK(last_grids.mach_surface_inner < first_grids.mach_surface_inner);
        CHECK(last_fixed.psi < first_fixed.psi);
    }
}
