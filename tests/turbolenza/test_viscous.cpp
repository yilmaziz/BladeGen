// The boundary-layer coupling (Viscous.hpp, the transpiration of BoundaryConditions.hpp): zero
// transpiration is v1's system exactly; the transpiration's relations, rows and psi_m; the wall arc
// length against a dense polyline; one way = the inviscid solve bit for bit; the constant wake's
// TE continuity; 1 and 4 threads bit-identical; the Reynolds-number helper; the blade-station fraction
// f of the layers and their files; the coupling defaults.  On the synthetic blade (no v1 data needed).
#include "CascadeFixture.hpp"
#include "SolverRef.hpp"

#include "turbolenza/Case.hpp"
#include "turbolenza/Coefficients.hpp"
#include "turbolenza/LinearSolver.hpp"
#include "turbolenza/MatrixAssembly.hpp"
#include "turbolenza/Viscous.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace BladeGen::Turbolenza;

namespace {

CaseSettings synthetic_settings()
{
    CaseSettings s;
    s.case_name = "synthetic viscous";
    s.chord = 0.127;
    s.solidity = 1.0;
    s.stagger_deg = 33.0;
    s.inlet_flow_angle_deg = 45.0;
    s.inlet_mach = 0.5;
    s.inlet_total_temperature = 335.0;
    s.inlet_total_pressure = 79000.0;
    s.grid_level = std::string("very_coarse");
    s.density_tolerance = 1e-6;
    s.kutta_tolerance = 1e-6;
    s.density_relaxation = 0.5;
    s.kutta_relaxation = 0.5;
    s.max_density_iterations = 3000;
    s.max_kutta_iterations = 3000;
    return s;
}

const Case& synthetic_case()
{
    static const Case c = build_case(synthetic_settings(), fixture::synthetic_section());
    return c;
}

bool same(const SparseMatrix& a, const SparseMatrix& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols() || a.nonZeros() != b.nonZeros()) return false;
    for (int c = 0; c < a.outerSize(); ++c) {
        SparseMatrix::InnerIterator x(a, c), y(b, c);
        for (; x && y; ++x, ++y)
            if (x.row() != y.row() || x.value() != y.value()) return false;
        if (x || y) return false;
    }
    return true;
}

bool same_records(const std::vector<PassRecord>& a, const std::vector<PassRecord>& b)
{
    if (a.size() != b.size()) return false;
    auto eq = [](double x, double y) { return x == y || (std::isnan(x) && std::isnan(y)); };
    for (std::size_t k = 0; k < a.size(); ++k) {
        const PassRecord &p = a[k], &q = b[k];
        if (p.index != q.index || p.final != q.final || p.bl_iteration != q.bl_iteration) return false;
        for (auto [x, y] : {std::pair{p.density_change, q.density_change}, std::pair{p.velocity_change, q.velocity_change},
                            std::pair{p.psi_change, q.psi_change}, std::pair{p.mach_max, q.mach_max},
                            std::pair{p.linear_residual, q.linear_residual}, std::pair{p.exit_angle, q.exit_angle},
                            std::pair{p.kutta_normalized, q.kutta_normalized},
                            std::pair{p.mass_flow_error_max, q.mass_flow_error_max}})
            if (!eq(x, y)) return false;
    }
    return true;
}

}  // namespace

TEST_CASE("viscous: zero transpiration is v1's system exactly", "[turbolenza][viscous]")
{
    const Case& c = synthetic_case();
    auto grid = std::make_shared<const Grid>(c.grid);
    const GridMetrics metrics = build_metrics(grid);
    const Field rho(grid->n_m(), grid->n_eta(), 0.8);
    const CoefficientFields co = build_coefficients(metrics, rho, c.mass_flow_per_depth, c.settings.chord);
    const int n_m = grid->n_m();
    Transpiration zero{std::vector<double>(static_cast<std::size_t>(n_m), 0.0),
                       std::vector<double>(static_cast<std::size_t>(n_m), 1.0),
                       std::vector<double>(static_cast<std::size_t>(n_m), 1.0)};
    const BoundaryConditions b0 = build_boundary_conditions(grid), b1 = build_boundary_conditions(grid, &zero);
    CHECK(!b0.transpiration);
    CHECK(b1.transpiration);
    REQUIRE(b0.relations.size() == b1.relations.size());
    for (std::size_t k = 0; k < b0.relations.size(); ++k) CHECK(b0.relations[k].rhs == b1.relations[k].rhs);
    const AssemblyResult a0 = assemble_from_coefficients(co, metrics, b0), a1 = assemble_from_coefficients(co, metrics, b1);
    CHECK(same(a0.matrix, a1.matrix));
    CHECK(a0.rhs == a1.rhs);
    // and the velocities of one solve: psi_m on the walls is the wall psi's derivative, exactly 0 here
    const Field psi = reconstruct(a0, solve(a0.matrix, a0.rhs));
    const FlowField f0 = velocity_components(psi, rho, co, metrics, b0), f1 = velocity_components(psi, rho, co, metrics, b1);
    CHECK(f0.W.v == f1.W.v);
    CHECK(f0.psi_m.v == f1.psi_m.v);
}

TEST_CASE("viscous: the transpiration's relations, rows and wall derivative", "[turbolenza][viscous]")
{
    const Case& c = synthetic_case();
    auto grid = std::make_shared<const Grid>(c.grid);
    const GridMetrics metrics = build_metrics(grid);
    const Field rho(grid->n_m(), grid->n_eta(), 0.8);
    const CoefficientFields co = build_coefficients(metrics, rho, c.mass_flow_per_depth, c.settings.chord);
    const int n_m = grid->n_m(), J = grid->n_eta() - 1, f = grid->first_blade_i, l = grid->last_blade_i;
    // lower psi falls linearly along the blade, upper psi rises, the exit jump grows
    const double s_lo = -2e-3, s_up = 3e-3;  // per station
    Transpiration t{std::vector<double>(static_cast<std::size_t>(n_m), 0.0),
                    std::vector<double>(static_cast<std::size_t>(n_m), 1.0),
                    std::vector<double>(static_cast<std::size_t>(n_m), 1.0)};
    for (int i = f; i <= l; ++i) {
        t.lower_wall[static_cast<std::size_t>(i)] = s_lo * (i - f);
        t.upper_wall[static_cast<std::size_t>(i)] = 1.0 + s_up * (i - f);
    }
    const double te_flow = t.upper_wall[static_cast<std::size_t>(l)] - t.lower_wall[static_cast<std::size_t>(l)];
    for (int i = l + 1; i < n_m; ++i) t.jump[static_cast<std::size_t>(i)] = te_flow;
    const BoundaryConditions b = build_boundary_conditions(grid, &t);
    for (const Relation& r : b.relations) {
        if (r.kind == RelationKind::LowerWall) CHECK(r.rhs == t.lower_wall[static_cast<std::size_t>(r.target->i)]);
        if (r.kind == RelationKind::UpperWall) CHECK(r.rhs == t.upper_wall[static_cast<std::size_t>(r.target->i)]);
        if (r.kind == RelationKind::Periodic) CHECK(r.rhs == t.jump[static_cast<std::size_t>(r.terms[0].node.i)]);
    }
    CHECK(b.station_flow(l) == te_flow);
    CHECK(b.station_flow(l + 1) == te_flow);
    CHECK(b.periodic_ghost(l + 2, J + 1, b.station_jump(l + 2)).constant == te_flow);
    // the solved psi meets every relation, the far exit's implied jump included (constant wake)
    const AssemblyResult a = assemble_from_coefficients(co, metrics, b);
    const Field psi = reconstruct(a, solve(a.matrix, a.rhs));
    double worst = 0.0;
    for (const Relation& r : b.relations) worst = std::max(worst, std::abs(r.residual(psi)));
    CHECK(sref::measure("viscous transpiration relations residual", worst) <= 1e-12);
    // psi_m on the walls: the run differences of a linear wall psi are its slope exactly
    const FlowField flow = velocity_components(psi, rho, co, metrics, b);
    const double dm = c.settings.chord * grid->dm();
    double off = 0.0;
    for (int i = f; i <= l; ++i) {
        off = std::max(off, std::abs(flow.psi_m(i, 0) - s_lo / dm) / std::abs(s_lo / dm));
        off = std::max(off, std::abs(flow.psi_m(i, J) - s_up / dm) / std::abs(s_up / dm));
    }
    CHECK(sref::measure("viscous psi_m on the walls vs d(psi_wall)/dm", off) <= 1e-9);
    // refusals
    Transpiration bad = t;
    bad.jump.pop_back();
    CHECK_THROWS_AS(build_boundary_conditions(grid, &bad), BoundaryConditionError);
    bad = t;
    bad.lower_wall[static_cast<std::size_t>(f)] = std::nan("");
    CHECK_THROWS_AS(build_boundary_conditions(grid, &bad), BoundaryConditionError);
}

TEST_CASE("viscous: the wall arc length against a dense polyline of the curve", "[turbolenza][viscous]")
{
    const Grid& g = synthetic_case().grid;
    for (int j : {0, g.n_eta() - 1}) {
        const std::vector<double> s = wall_arc_length(g, j);
        const WallSide& side = j == 0 ? g.wall_identity.lower : g.wall_identity.upper;
        const ParametricBladeCurve& curve =
            side.surface[0] == Surface::Suction ? g.boundaries.blade_curves->suction : g.boundaries.blade_curves->pressure;
        double worst = 0.0;
        for (std::size_t b = 0; b < s.size(); b += 3) {
            const int n = 200000;
            double poly = 0.0;
            Point prev = curve.evaluate(0.0);
            for (int k = 1; k <= n; ++k) {
                const Point p = curve.evaluate(side.u[b] * k / n);
                poly += std::hypot(p.m - prev.m, p.y - prev.y);
                prev = p;
            }
            worst = std::max(worst, std::abs(s[b] - poly));
            if (b > 0) CHECK(s[b] > s[b - 1]);
        }
        CHECK(sref::measure(j == 0 ? "viscous arc length lower vs polyline" : "viscous arc length upper vs polyline", worst) <=
              1e-8);
    }
    Grid touching = g;
    touching.edge_nodes = EdgeNodes::Touching;
    CHECK_THROWS_AS(wall_arc_length(touching, 0), BoundaryLayerError);
}

TEST_CASE("viscous: one way is the inviscid solve bit for bit; two way converges; threads", "[turbolenza][viscous][slow]")
{
    const CaseSettings base = synthetic_settings();
    const RunResult inviscid = run_case(build_case(base, fixture::synthetic_section()));
    REQUIRE(inviscid.converged);
    CHECK(!inviscid.viscous);

    CaseSettings one = base;
    one.viscous = true;
    one.bl_coupling = BlCoupling::OneWay;
    const RunResult r1 = run_case(build_case(one, fixture::synthetic_section()));
    INFO(r1.error);
    REQUIRE(r1.converged);
    REQUIRE(r1.viscous);
    CHECK(same_records(r1.records, inviscid.records));
    CHECK(r1.solution->psi.v == inviscid.solution->psi.v);
    CHECK(r1.solution->exit_angle == inviscid.solution->exit_angle);
    CHECK(r1.viscous->history.size() == 1);
    CHECK(r1.viscous->history[0].max_change_m == 1.0);
    CHECK(std::isfinite(r1.viscous->loss_omega));
    CHECK(r1.viscous->loss_omega > 0.0);

    CaseSettings two = base;
    two.viscous = true;
    two.bl_coupling = BlCoupling::TwoWay;  // the default omega_BL and max_bl_iterations
    const RunResult r2 = run_case(build_case(two, fixture::synthetic_section()));
    INFO(r2.error);
    REQUIRE(r2.viscous);
    const ViscousData& v = *r2.viscous;
    REQUIRE(!v.history.empty());
    sref::measure("viscous synthetic two-way BL iterations", static_cast<double>(v.history.size()));
    sref::measure("viscous synthetic two-way beta_out k = 0 [deg]", v.history.front().beta_out * 180.0 / 3.141592653589793);
    sref::measure("viscous synthetic two-way beta_out coupled [deg]", v.history.back().beta_out * 180.0 / 3.141592653589793);
    sref::measure("viscous synthetic two-way loss omega", v.loss_omega);
    CHECK(r2.converged);
    CHECK(v.history.back().max_change_m <= two.bl_tolerance);
    CHECK(v.te_jump_mismatch <= 1e-15);  // the constant wake: continuous at the TE by construction
    // the passes of every BL iteration are in the records, labeled
    std::size_t passes = 0;
    for (const BlIteration& h : v.history) passes += static_cast<std::size_t>(h.inviscid_passes);
    CHECK(passes == r2.records.size());
    CHECK(r2.records.back().bl_iteration == v.history.back().k);

    // the same run on 4 threads: bit-identical (SparseLU)
    if (max_threads() >= 4) {
        CaseSettings four = two;
        four.threads = 4;
        const RunResult r4 = run_case(build_case(four, fixture::synthetic_section()));
        REQUIRE(r4.viscous);
        CHECK(same_records(r4.records, r2.records));
        CHECK(r4.solution->psi.v == r2.solution->psi.v);
        CHECK(r4.viscous->suction.delta_star == v.suction.delta_star);
        CHECK(r4.viscous->pressure_side.theta == v.pressure_side.theta);
        CHECK(r4.viscous->loss_omega == v.loss_omega);
    }

    // the marched wake converges too, and the report and files carry the layer
    CaseSettings marched = two;
    marched.wake_model = WakeModel::Marched;
    const RunResult rm = run_case(build_case(marched, fixture::synthetic_section()));
    INFO(rm.error);
    CHECK(rm.converged);
    REQUIRE(rm.viscous);
    CHECK(!rm.viscous->wake_theta.empty());
    sref::measure("viscous synthetic marched wake TE jump mismatch", rm.viscous->te_jump_mismatch);
    CHECK(format_run_report(rm).find("Boundary layer") != std::string::npos);
}

TEST_CASE("viscous: a run that stops in BL iteration k >= 1 still reports the last layers computed", "[turbolenza][viscous][slow]")
{
    // The two-way run is cancelled at the third pass of BL iteration 1 (the run then ends like a failed
    // solve, with that solve's last completed pass reported).  The layers of iteration 0 must stay in
    // the result, labeled as such: they are the one-way layers of the inviscid solution, bit for bit.
    CaseSettings one = synthetic_settings();
    one.viscous = true;
    one.bl_coupling = BlCoupling::OneWay;
    const RunResult r1 = run_case(build_case(one, fixture::synthetic_section()));
    REQUIRE(r1.converged);
    REQUIRE(r1.viscous);

    CaseSettings two = one;
    two.bl_coupling = BlCoupling::TwoWay;
    std::atomic<bool> flag{false};
    struct Stop final : RunObserver
    {
        std::atomic<bool>* f = nullptr;
        int seen = 0;
        void update(const PassRecord& p) override
        {
            if (p.bl_iteration == 1 && ++seen == 3) f->store(true);
        }
    } stop;
    stop.f = &flag;
    const RunResult r = run_case(build_case(two, fixture::synthetic_section()), &stop, &flag);
    INFO(r.error);
    CHECK(!r.converged);
    CHECK(r.error.starts_with("Cancelled"));
    REQUIRE(r.viscous);
    const ViscousData& v = *r.viscous;
    CHECK(v.layers_iteration == 0);
    CHECK(!v.layers_of_last_solve);       // the reported solution is iteration 1's last pass, not the layers' owner
    REQUIRE(!v.suction.x_c.empty());      // before the fix: empty ("no boundary layer in this run")
    CHECK(v.suction.H == r1.viscous->suction.H);
    CHECK(v.pressure_side.theta == r1.viscous->pressure_side.theta);
    CHECK(v.loss_omega == r1.viscous->loss_omega);  // before the fix: NaN
    CHECK(format_run_report(r).find("the last one computed") != std::string::npos);
}

TEST_CASE("viscous: the total pressure of a reference Reynolds number", "[turbolenza][viscous]")
{
    const double P0 = total_pressure_for_reynolds(550000.0, 0.30, 335.0, 0.127, 1.4, 287.05);
    sref::measure("viscous Dunavant P0 [Pa] for Re 5.5e5 at M 0.30", P0);
    CHECK(std::abs(chord_reynolds_number(0.30, P0, 335.0, 0.127, 1.4, 287.05) / 550000.0 - 1.0) < 1e-14);
    CHECK(std::abs(P0 - 79.1e3) < 0.2e3);  // the case's ~79.1 kPa
    // by hand: T1 = 335 / 1.018, a1 = sqrt(1.4 R T1), mu = Sutherland(T1), rho1 = Re mu / (W1 c), P0 = rho1 R T1 1.018^3.5
    const double T1 = 335.0 / 1.018, W1 = 0.3 * std::sqrt(1.4 * 287.05 * T1), mu = sutherland_viscosity(T1);
    const double hand = 550000.0 * mu / (W1 * 0.127) * 287.05 * T1 * std::pow(1.018, 3.5);
    CHECK(std::abs(P0 / hand - 1.0) < 1e-13);
    sref::measure("viscous Dunavant chord Re at M 0.65", chord_reynolds_number(0.65, P0, 335.0, 0.127, 1.4, 287.05));
    CHECK_THROWS_AS(total_pressure_for_reynolds(-1.0, 0.3, 335.0, 0.127, 1.4, 287.05), CaseError);
    CHECK_THROWS_AS(chord_reynolds_number(0.0, 1e5, 335.0, 0.127, 1.4, 287.05), CaseError);
}

TEST_CASE("viscous: build_case refuses what the coupling cannot run", "[turbolenza][viscous]")
{
    CaseSettings s = synthetic_settings();
    s.viscous = true;
    CHECK_NOTHROW(build_case(s, fixture::synthetic_section()));
    for (auto set : std::vector<void (*)(CaseSettings&)>{
             [](CaseSettings& x) { x.bl_relaxation = 0.0; }, [](CaseSettings& x) { x.bl_relaxation = 1.5; },
             [](CaseSettings& x) { x.bl_tolerance = 0.0; }, [](CaseSettings& x) { x.max_bl_iterations = -1; },
             [](CaseSettings& x) { x.transition_re_theta = -250.0; }, [](CaseSettings& x) { x.separation_hk = 1.0; },
             [](CaseSettings& x) { x.sutherland.s = -1.0; }, [](CaseSettings& x) { x.edge_nodes = EdgeNodes::Touching; }}) {
        CaseSettings bad = s;
        set(bad);
        CHECK_THROWS_AS(build_case(bad, fixture::synthetic_section()), CaseError);
    }
}

TEST_CASE("viscous: the blade-station fraction f", "[turbolenza][viscous]")
{
    const Grid& g = synthetic_case().grid;
    const std::vector<double> f = blade_station_fraction(g);
    const int first = g.first_blade_i, last = g.last_blade_i, n = last - first;
    REQUIRE(f.size() == static_cast<std::size_t>(n + 1));
    // by hand: the stations are uniform in m (GridData.hpp), so f of blade station b is b / (last - first)
    CHECK(f.front() == 0.0);
    CHECK(f.back() == 1.0);
    const int mid = n / 2;
    CHECK(std::abs(f[static_cast<std::size_t>(mid)] - static_cast<double>(mid) / n) <= 1e-14);
    double worst = 0.0;
    for (int b = 0; b <= n; ++b) {
        worst = std::max(worst, std::abs(f[static_cast<std::size_t>(b)] - static_cast<double>(b) / n));
        // and its definition, (m_i - m_first) / (m_last - m_first)
        const double m = g.m[static_cast<std::size_t>(first + b)], m0 = g.m[static_cast<std::size_t>(first)],
                     m1 = g.m[static_cast<std::size_t>(last)];
        if (b < n) CHECK(f[static_cast<std::size_t>(b)] == (m - m0) / (m1 - m0));
    }
    CHECK(sref::measure("viscous station fraction vs b / (last - first)", worst) <= 1e-14);
    // degenerate grids are refused, not repaired
    Grid one = g;
    one.last_blade_i = one.first_blade_i;
    CHECK_THROWS_AS(blade_station_fraction(one), BoundaryLayerError);
    Grid flat = g;
    flat.m[static_cast<std::size_t>(first + 1)] = flat.m[static_cast<std::size_t>(first)];
    CHECK_THROWS_AS(blade_station_fraction(flat), BoundaryLayerError);
    Grid out = g;
    out.last_blade_i = static_cast<int>(out.m.size());
    CHECK_THROWS_AS(blade_station_fraction(out), BoundaryLayerError);
}

TEST_CASE("viscous: the layers and bl_*.csv carry f beside x/c", "[turbolenza][viscous]")
{
    CaseSettings s = synthetic_settings();
    s.viscous = true;
    s.bl_coupling = BlCoupling::OneWay;
    const Case c = build_case(s, fixture::synthetic_section());
    const RunResult r = run_case(c);
    INFO(r.error);
    REQUIRE(r.converged);
    REQUIRE(r.viscous);
    const std::vector<double> f = blade_station_fraction(c.grid);
    const ViscousData& v = *r.viscous;
    for (const BoundaryLayerSurfaceData* d : {&v.suction, &v.pressure_side}) {
        CHECK(d->f == f);
        REQUIRE(d->x_c.size() == f.size());
        // the event abscissae are the stations' own, on both axes
        int tr = -1;
        for (std::size_t k = 0; k < d->regime.size(); ++k)
            if (d->regime[k] == "transition") tr = static_cast<int>(k);
        REQUIRE(tr >= 0);
        CHECK(d->transition_f == f[static_cast<std::size_t>(tr)]);
        CHECK(d->transition_x_c == d->x_c[static_cast<std::size_t>(tr)]);
        CHECK(std::isnan(d->laminar_separation_f) == std::isnan(d->laminar_separation_x_c));
        CHECK(std::isnan(d->separation_f) == std::isnan(d->separation_x_c));
    }
    const auto dir = std::filesystem::path(BLADEGEN_TEST_OUTPUT) / "viscous_station_fraction";
    std::filesystem::remove_all(dir);
    (void)export_run(r, c, dir.string());
    for (const auto& [name, d] : {std::pair{"bl_suction.csv", &v.suction}, std::pair{"bl_pressure.csv", &v.pressure_side}}) {
        std::ifstream in(dir / name);
        REQUIRE(in);
        std::string line;
        std::getline(in, line);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        CHECK(line == "x_c,f,x,W_e,rho_e,T_e,mu,theta,delta_star,delta,H,H_k,c_f,Lambda,Re_theta,regime,mass_deficit");
        std::size_t k = 0;
        while (std::getline(in, line)) {
            std::stringstream row(line);
            std::string x_c, fk, x;
            std::getline(row, x_c, ',');
            std::getline(row, fk, ',');
            std::getline(row, x, ',');
            REQUIRE(k < d->f.size());
            CHECK(std::stod(x_c) == d->x_c[k]);
            CHECK(std::stod(fk) == d->f[k]);
            CHECK(std::stod(x) == d->x[k]);
            ++k;
        }
        CHECK(k == d->f.size());
    }
}

TEST_CASE("viscous: the coupling defaults are omega_BL 0.2 and 60 BL iterations", "[turbolenza][viscous]")
{
    // PORT.md "The Fig. 5-15 study": 0.5 diverges on the Dunavant case, 0.2 converges in 35 iterations
    const CaseSettings s;
    CHECK(s.bl_relaxation == 0.2);
    CHECK(s.max_bl_iterations == 60);
    const ViscousSettings v;
    CHECK(v.relaxation == 0.2);
    CHECK(v.max_iterations == 60);
    CHECK(v.tolerance == s.bl_tolerance);
}
