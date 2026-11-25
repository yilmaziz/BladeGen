// solve_general (GeneralSolver.hpp) on cases that need no reference data: a rotating radial stream sheet
// with r(m), b(m) against its closed form (the rotation source, the rothalpy at every node's radius, the
// velocity recovery with r and b, the rotation sign of the contract), the Kutta-coupled solve of a
// synthetic cascade against the Aungier solver on the same nodes, and the orientation of the outlet sums.
#include "CascadeFixture.hpp"
#include "GeneralFixture.hpp"

#include "turbolenza/Detail.hpp"
#include "turbolenza/GeneralAdapter.hpp"
#include "turbolenza/GeneralSolver.hpp"
#include "turbolenza/Grid.hpp"
#include "turbolenza/Solver.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <format>
#include <memory>

using namespace BladeGen::Turbolenza;
using gfix::measure;

namespace {

std::size_t z(int k) { return static_cast<std::size_t>(k); }

/// A radial stream sheet r = r0 + m, b = 1 - 0.5 m, turning at Omega (positive rpm), with no absolute swirl
/// (C_theta = 0): absolute vorticity zero, W_m = mdot' / (rho b r Theta), W_theta = C_theta - omega r =
/// Omega r (omega = -Omega, the contract), h = I + (omega r)^2 / 2 - W^2 / 2 = cp Tt - W_m^2 / 2, rho
/// isentropic to the inlet totals.  The walls are the streamlines theta_w(m) = (Theta / mdot') int_0^m
/// rho b Omega r dm (psi = 0) and theta_w + Theta (psi = 1): the closed form of the whole solve.
struct RotatingSheet
{
    double r0 = 0.2, L = 0.3, Omega = 300.0, Theta = 2.0 * gfix::kPi / 20.0, Wm1 = 80.0;
    PerfectGas gas;
    double Tt = 288.15, Pt = 101325.0;
    double mdot = 0.0;

    double r(double m) const { return r0 + m; }
    double b(double m) const { return 1.0 - 0.5 * m; }
    double rho_t() const { return Pt / (gas.gas_constant * Tt); }
    double rho_of(double W_m) const
    {
        return rho_t() * std::pow(1.0 - W_m * W_m / (2.0 * gas.cp() * Tt), 1.0 / (gas.gamma - 1.0));
    }
    /// rho(m) on the subsonic branch: the fixed point of rho = rho_of(mdot' / (rho b r Theta)).
    double rho(double m) const
    {
        double x = rho_t();
        for (int k = 0; k < 200; ++k) {
            const double next = rho_of(mdot / (x * b(m) * r(m) * Theta));
            if (next == x) break;
            x = next;
        }
        return x;
    }
    double W_m(double m) const { return mdot / (rho(m) * b(m) * r(m) * Theta); }
    double W_theta(double m) const { return Omega * r(m); }
    /// theta_w(m) by Simpson's rule on 4000 intervals per unit m (converged to round-off).
    double theta_w(double m) const
    {
        if (m == 0.0) return 0.0;
        const int n = 2 * std::max(1, static_cast<int>(2000.0 * m));
        const double h = m / n;
        auto f = [&](double s) { return rho(s) * b(s) * Omega * r(s); };
        double sum = f(0.0) + f(m);
        for (int k = 1; k < n; ++k) sum += (k % 2 ? 4.0 : 2.0) * f(k * h);
        return Theta / mdot * sum * h / 3.0;
    }
    RotatingSheet() { mdot = rho_of(Wm1) * Wm1 * b(0.0) * r0 * Theta; }

    GeneralMesh mesh(int ni, int nj) const
    {
        GeneralMesh M;
        M.pitch = Theta;
        MeshBlock B;
        B.name = "rotating radial sheet";
        B.ni = ni;
        B.nj = nj;
        const std::size_t n = static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
        B.x.resize(n);
        B.y.resize(n);
        B.m.resize(n);
        B.r.resize(n);
        B.b.resize(n);
        for (int i = 0; i < ni; ++i) {
            const double m = L * i / (ni - 1), tw = theta_w(m);
            for (int j = 0; j < nj; ++j) {
                const std::size_t k = B.at(i, j);
                B.x[k] = std::log(r(m) / r0);  // m' = int dm / r
                B.y[k] = tw + Theta * j / (nj - 1);
                B.m[k] = m;
                B.r[k] = r(m);
                B.b[k] = b(m);
            }
        }
        M.blocks.push_back(std::move(B));
        FaceSegment lower = gfix::whole(M, 0, BlockFace::JMin, SegmentKind::Blade);
        lower.blade = 0;
        lower.surface = Surface::Suction;
        FaceSegment upper = gfix::whole(M, 0, BlockFace::JMax, SegmentKind::Blade);
        upper.blade = 1;
        upper.surface = Surface::Pressure;
        M.segments = {lower, upper, gfix::whole(M, 0, BlockFace::IMin, SegmentKind::Inlet),
                      gfix::whole(M, 0, BlockFace::IMax, SegmentKind::Outlet)};
        M.tags = {{EdgeKind::Leading, 0, Surface::Suction, 0, 0, 0}, {EdgeKind::Trailing, 0, Surface::Suction, 0, ni - 1, 0},
                  {EdgeKind::Leading, 1, Surface::Pressure, 0, 0, nj - 1},
                  {EdgeKind::Trailing, 1, Surface::Pressure, 0, ni - 1, nj - 1}};
        return M;
    }

    GeneralConfig config(RotationSense sense) const
    {
        GeneralConfig c;
        c.gas = gas;
        c.flow.inlet_absolute_total_pressure = Pt;
        c.flow.inlet_absolute_total_temperature = Tt;
        c.flow.inlet_angle_deg = detail::degrees(std::atan(W_theta(0.0) / W_m(0.0)));
        c.flow.exit_angle_guess_deg = detail::degrees(std::atan(W_theta(L) / W_m(L)));
        c.flow.angular_speed = sense == RotationSense::AgainstTheta ? Omega : -Omega;
        c.rotation = sense;
        c.flow.stream_surface_radius = r0;
        c.flow.inlet_absolute_tangential_velocity = 0.0;
        c.mass_flow_per_depth = mdot;
        c.kutta = false;
        c.iteration.density_tolerance = 1e-12;
        c.iteration.max_density_iterations = 400;
        c.iteration.density_relaxation = 0.7;
        return c;
    }
};

}  // namespace

TEST_CASE("general solver: a rotating radial stream sheet against its closed form", "[turbolenza][general]")
{
    const RotatingSheet S;
    std::vector<double> e_rho, e_wm, e_wt, e_psi;
    for (int N : {8, 16, 32, 64}) {
        const GeneralSolution s = solve_general(S.mesh(2 * N + 1, N + 1), S.config(RotationSense::AgainstTheta));
        REQUIRE(s.converged);
        const GeneralOperator& op = *s.op;
        const int ni = 2 * N + 1;
        std::vector<double> rho(z(ni)), wm(z(ni)), wt(z(ni)), tw(z(ni));
        for (int i = 0; i < ni; ++i) {
            const double m = S.L * i / (ni - 1);
            rho[z(i)] = S.rho(m);
            wm[z(i)] = S.W_m(m);
            wt[z(i)] = S.W_theta(m);
            tw[z(i)] = S.theta_w(m);
        }
        double er = 0.0, em = 0.0, et = 0.0, ep = 0.0;
        bool positive = true;
        for (int g = 0; g < op.n_nodes(); ++g) {
            const std::size_t i = z(op.location(g)[1]);
            er = std::max(er, std::abs(s.density[z(g)] / rho[i] - 1.0));
            em = std::max(em, std::abs(s.W_m[z(g)] / wm[i] - 1.0));
            et = std::max(et, std::abs(s.W_theta[z(g)] / wt[i] - 1.0));
            ep = std::max(ep, std::abs(s.psi[z(g)] - (op.y(g) - tw[i]) / S.Theta));
            positive = positive && s.W_theta[z(g)] > 0.0;
        }
        // positive rpm, axial (swirl-free) absolute inflow: a positive relative angle (the contract)
        CHECK(positive);
        e_rho.push_back(er);
        e_wm.push_back(em);
        e_wt.push_back(et);
        e_psi.push_back(ep);
        measure(std::format("general rotating sheet N={} max rel error rho", N), er);
        measure(std::format("general rotating sheet N={} max rel error W_m", N), em);
        measure(std::format("general rotating sheet N={} max rel error W_theta", N), et);
        measure(std::format("general rotating sheet N={} max error psi", N), ep);
    }
    // measured: orders 1.80, 1.82, 1.91; errors at N = 64 3.0e-5, 3.8e-4, 3.7e-4
    CHECK(measure("general rotating sheet order rho 32->64", gfix::order(e_rho[2], e_rho[3])) > 1.7);
    CHECK(measure("general rotating sheet order W_m 32->64", gfix::order(e_wm[2], e_wm[3])) > 1.7);
    CHECK(measure("general rotating sheet order W_theta 32->64", gfix::order(e_wt[2], e_wt[3])) > 1.7);
    CHECK(e_rho[3] < 4.5e-5);
    CHECK(e_wm[3] < 6e-4);
    CHECK(e_wt[3] < 6e-4);

    // the same row in the book's sign (omega_book = angular_speed): bit for bit the same solve
    const GeneralSolution a = solve_general(S.mesh(17, 9), S.config(RotationSense::AgainstTheta));
    const GeneralSolution b = solve_general(S.mesh(17, 9), S.config(RotationSense::AlongTheta));
    CHECK(a.psi == b.psi);
    CHECK(a.density == b.density);
    // and a rotating row needs its mass flow
    GeneralConfig c = S.config(RotationSense::AgainstTheta);
    c.mass_flow_per_depth.reset();
    c.flow.inlet_mach = 0.3;
    CHECK_THROWS_AS(solve_general(S.mesh(9, 5), c), GeneralOperatorError);
}

TEST_CASE("general solver: the Kutta-coupled synthetic cascade against the Aungier solver", "[turbolenza][general]")
{
    // no reference data: the synthetic blade of the fixture, the Aungier solver's own coupled run, then the
    // general Kutta loop on its final grid started at the grid's first discharge angle (the TE blade angle)
    const Cascade cascade = make_cascade(fixture::synthetic_section(), 1.25, 30.0);
    const double c_ref = 0.5;
    for (const char* level : {"very_coarse", "coarse"}) {
        const auto grid = std::make_shared<const Grid>(build_grid(cascade, level, detail::radians(45.0), std::nullopt).grid);
        SolverConfig ac;
        ac.flow.inlet_mach = 0.5;
        ac.flow.inlet_angle_deg = 45.0;
        ac.iteration.density_tolerance = 1e-9;
        ac.iteration.kutta_tolerance = 1e-9;
        ac.iteration.max_density_iterations = 3000;
        ac.iteration.max_kutta_iterations = 3000;
        ac.iteration.density_relaxation = 0.5;
        ac.iteration.kutta_relaxation = 0.5;
        const Solution ak = solve_coupled(grid, ac, c_ref);
        REQUIRE(ak.converged);
        GeneralConfig gc;
        gc.flow = ac.flow;
        gc.flow.exit_angle_guess_deg = detail::degrees(grid->exit_angle());
        gc.iteration = ac.iteration;
        gc.kutta = true;
        const GeneralSolution gk = solve_general(mesh_from_grid(*ak.grid, c_ref), gc);
        REQUIRE(gk.converged);
        REQUIRE(gk.kutta);
        CHECK(gk.kutta->normalized <= 1e-9);
        CHECK(gk.mass_flow_per_depth == ak.mass_flow_per_depth);
        const Grid& G = *ak.grid;
        double dM = 0.0;
        for (int i = G.first_blade_i; i <= G.last_blade_i; ++i) {
            if (!(G.m[z(i)] >= G.m_le() + 0.1 && G.m[z(i)] <= G.m_te() - 0.1)) continue;
            for (int j : {0, G.n_eta() - 1})
                dM = std::max(dM, std::abs(ak.mach(i, j) - gk.mach[z(general_node_of(G, i, j))]));
        }
        const double dbeta = detail::degrees(gk.exit_angle - ak.exit_angle);
        const std::string tag = std::format("synthetic M1 0.5 {}", level);
        measure("general vs Aungier Kutta beta_out difference [deg] " + tag, dbeta);
        measure("general vs Aungier max |dM| on the walls 0.1 .. 0.9 c " + tag, dM);
        std::printf("  %s: beta_out Aungier %.6f deg (%d passes), general %.6f deg (%d passes, from %.3f deg)\n",
                    tag.c_str(), detail::degrees(ak.exit_angle), ak.iterations, detail::degrees(gk.exit_angle),
                    gk.iterations, detail::degrees(grid->exit_angle()));
        // regression bounds, about 1.5 x measured (beta_out -2.95 / -2.23 deg, mid-chord dM 9.1e-3 / 8.9e-3 on
        // very_coarse / coarse): the two discretizations differ most at the trailing-edge transition cells
        CHECK(std::abs(dbeta) <= 4.5);
        CHECK(dM <= 0.015);
    }
}

TEST_CASE("general solver: outlet sums take every segment toward +theta", "[turbolenza][general]")
{
    // two blocks stacked in theta, the upper one turned by 180 degrees: its outlet face (imin) runs toward
    // -theta.  A uniform W_m = 2, W_theta = 0, rho = b = r = 1: the flux is 2 times the theta extent 1, the
    // mean W_m 2 (a signed sum would give 0 and refuse).
    GeneralMesh M;
    M.pitch = 1.0;
    M.blocks.push_back(gfix::block([](double u, double v) { return std::array<double, 2>{2.0 * u, 0.5 * v}; }, 9, 5));
    M.blocks.push_back(gfix::block([](double u, double v) { return std::array<double, 2>{2.0 * (1.0 - u), 1.0 - 0.5 * v}; }, 9, 5));
    using F = BlockFace;
    using K = SegmentKind;
    M.segments = {gfix::whole(M, 0, F::IMin, K::Dirichlet), gfix::whole(M, 0, F::IMax, K::Outlet),
                  gfix::whole(M, 0, F::JMin, K::Dirichlet), gfix::whole(M, 0, F::JMax, K::Interface),
                  gfix::whole(M, 1, F::IMax, K::Dirichlet), gfix::whole(M, 1, F::IMin, K::Outlet),
                  gfix::whole(M, 1, F::JMin, K::Dirichlet), gfix::whole(M, 1, F::JMax, K::Interface)};
    gfix::connect(M, 3, 7, 0);
    const GeneralOperator op(M);
    const std::size_t n = static_cast<std::size_t>(op.n_nodes());
    const std::vector<double> rho(n, 1.0), Wm(n, 2.0), Wt(n, 0.0);
    CHECK(std::abs(general_boundary_flux(op, SegmentKind::Outlet, rho, Wm, Wt) - 2.0) <= 1e-14);
    GeneralTrailingEdge te;
    te.suction = op.node(0, 4, 0);
    te.pressure = op.node(1, 4, 0);
    te.plus = te.suction;
    te.minus = te.pressure;
    te.passage_width = 1.0;
    te.dm = 0.1;
    te.cos_plus = te.cos_minus = 1.0;
    const GeneralKutta k = general_kutta(op, te, Wm, std::vector<double>(n, 3.0));
    CHECK(std::abs(k.outlet_mean_W_m - 2.0) <= 1e-14);
    CHECK(k.delta_tan == 0.0);  // equal TE speeds
}
