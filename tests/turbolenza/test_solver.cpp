// Solver (port of v1 solver.py): the pass records and fields of v1 passes 0 and 1 of the coupled
// M1 0.65 case (pass 1 on v1's rebuilt grid, injected), the fixed-angle loop end to end on the
// explicit 45 x 21 grid, thread independence, PARDISO against SparseLU, the observers, cancel,
// errors carrying the last pass, the refusals.  [needs-grid]: the whole coupled runs with the
// grid port's rebuild_exit.
#include "SolverRef.hpp"

#include "turbolenza/Detail.hpp"
#include "turbolenza/Solver.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <chrono>
#include <cmath>
#include <memory>

using namespace BladeGen::Turbolenza;

namespace {

/// The v1 grids of the dumped passes: rebuild_exit at any new angle gives pass 1's grid.
GridOps injected(const ref::json& d)
{
    auto pass1 = std::make_shared<Grid>(ref::grid(d.at("passes").at(1).at("grid")));
    GridOps ops;
    ops.rebuild_exit = [pass1](const Grid&, double) { return *pass1; };
    ops.validate_topology = [](const Grid&) {};
    return ops;
}

struct Counter final : RunObserver
{
    int starts = 0, updates = 0, finals = 0, finishes = 0;
    std::vector<PassRecord> records;
    bool converged = false;
    std::string error;
    void start(const SolveInfo&) override { ++starts; }
    void update(const PassRecord& r) override
    {
        ++updates;
        finals += r.final ? 1 : 0;
        records.push_back(r);
    }
    void finish(bool c, const std::string& e) override
    {
        ++finishes;
        converged = c;
        error = e;
    }
};

struct Thrower final : RunObserver
{
    void update(const PassRecord&) override { throw std::runtime_error("observer failure"); }
};

/// The dumped first grid with v1's C2 curves attached (rebuild_exit needs them; RefData leaves them out).
std::shared_ptr<const Grid> grid_with_curves(const ref::json& d)
{
    Grid g = ref::grid(d.at("first_grid"));
    g.boundaries.blade_curves =
        std::make_shared<const SurfaceCurves>(ref::surface_curves(d.at("cascade").at("active_curves")));
    return std::make_shared<const Grid>(std::move(g));
}

bool same_bits(const Field& a, const Field& b) { return a.same_shape(b) && a.v == b.v; }

bool same_record(const PassRecord& a, const PassRecord& b)
{
    auto eq = [](double x, double y) { return (std::isnan(x) && std::isnan(y)) || x == y; };
    return a.index == b.index && a.final == b.final && eq(a.density_change, b.density_change) &&
           eq(a.velocity_change, b.velocity_change) && eq(a.psi_change, b.psi_change) &&
           eq(a.target_mismatch, b.target_mismatch) && eq(a.mach_max, b.mach_max) &&
           eq(a.linear_residual, b.linear_residual) && eq(a.mass_flow_error_max, b.mass_flow_error_max) &&
           eq(a.exit_angle_next, b.exit_angle_next) && eq(a.kutta_residual, b.kutta_residual) &&
           eq(a.speed_max, b.speed_max) && eq(a.density_min, b.density_min);
}

/// The fixed-angle 45 x 21 case.
struct FixedCase
{
    std::shared_ptr<const Grid> grid;
    SolverConfig config;
    double c_ref = 0.0;
};
FixedCase fixed_case(const ref::json& d)
{
    return {sref::grid(d.at("first_grid")), sref::config(d), ref::num(d.at("settings").at("chord"))};
}

}  // namespace

TEST_CASE("solver: v1 coupled passes 0 and 1 (pass 1 on v1's rebuilt grid)", "[turbolenza][solver][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    SolverConfig cfg = sref::config(d);
    cfg.iteration.max_density_iterations = 2;
    const double c_ref = ref::num(d.at("settings").at("chord"));
    SolveOptions o;
    o.grid_ops = injected(d);
    std::vector<std::string> checks;
    int seen = 0;
    o.inspect = [&](const PassFields& f) {
        const ref::json& p = d.at("passes").at(static_cast<std::size_t>(seen));
        const std::string tag = " pass " + std::to_string(seen);
        // pass 0 from the start state: SuperLU (v1) vs SparseLU round-off only; pass 1 carries
        // pass 0's differences through rho(1) and beta_out(1)
        CHECK(sref::measure("solver density_resident" + tag, ref::rel_diff(f.density_resident, ref::field(p.at("density_resident")))) <=
              1e-13);
        CHECK(sref::measure("solver coefficients A" + tag, ref::rel_diff(f.coefficients.A, ref::field(p.at("coefficients").at("A")))) <=
              1e-13);
        CHECK(sref::measure("solver coefficients D" + tag, ref::rel_diff(f.coefficients.D, ref::field(p.at("coefficients").at("D")))) <=
              1e-11);
        CHECK(sref::measure("solver psi" + tag, ref::rel_diff(f.linear.psi, ref::field(p.at("psi")))) <= 1e-13);
        CHECK(sref::measure("solver W" + tag, ref::rel_diff(f.flow.W, ref::field(p.at("flow").at("W")))) <= 1e-12);
        CHECK(sref::measure("solver mach" + tag, ref::rel_diff(f.mach, ref::field(p.at("mach")))) <= 1e-12);
        CHECK(sref::measure("solver density_next" + tag, ref::rel_diff(f.density_next, ref::field(p.at("density_next")))) <=
              1e-13);
        ++seen;
    };
    try {
        (void)solve_coupled(sref::grid(d.at("first_grid")), cfg, c_ref, std::nullopt, o);
        FAIL("two passes cannot converge to 1e-8");
    } catch (const CoupledNotConvergedError& e) {
        REQUIRE(e.solution);
        REQUIRE(e.solution->passes.size() == 2);
        CHECK(e.solution->iterations == 2);
        CHECK(!e.solution->converged);
        for (int k = 0; k < 2; ++k) {
            std::string worst;
            const double diff = sref::record_diff(e.solution->passes[static_cast<std::size_t>(k)],
                                                  d.at("records").at(static_cast<std::size_t>(k)), &worst);
            // measured 9.4e-14 (pass 0, eps_K) and 3.2e-14 (pass 1); bound 1e-12
            CHECK(sref::measure("solver record pass " + std::to_string(k) + " (" + worst + ")", diff) <= 1e-12);
        }
    }
    CHECK(seen == 2);
}

TEST_CASE("solver: the fixed-angle loop end to end (explicit 45 x 21)", "[turbolenza][solver][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "fixed_exit_explicit_45x21");
    const FixedCase fc = fixed_case(d);
    Counter counter;
    SolveOptions o;
    o.observer = &counter;
    const Solution s = solve_at_exit_angle(fc.grid, fc.config, fc.c_ref, std::nullopt, o);
    const ref::json& v1 = d.at("solution");
    CHECK(s.converged);
    CHECK(s.iterations == v1.at("iterations").get<int>());
    REQUIRE(s.passes.size() == d.at("records").size());
    double worst = 0.0;
    std::string where;
    for (std::size_t k = 0; k < s.passes.size(); ++k) {
        std::string w;
        const double diff = sref::record_diff(s.passes[k], d.at("records").at(k), &w);
        if (diff > worst) {
            worst = diff;
            where = std::to_string(k) + " " + w;
        }
    }
    // measured 1.0e-12 (pass 5, eps_K, absolute); bound 1e-11
    CHECK(sref::measure("solver fixed 45x21 records worst (" + where + ")", worst) <= 1e-11);
    // measured 9.0e-15, 2.9e-14, 4.4e-14
    CHECK(sref::measure("solver fixed 45x21 density", ref::rel_diff(s.density, ref::field(v1.at("density")))) <= 1e-13);
    CHECK(sref::measure("solver fixed 45x21 psi", ref::rel_diff(s.psi, ref::field(v1.at("psi")))) <= 5e-13);
    CHECK(sref::measure("solver fixed 45x21 mach", ref::rel_diff(s.mach, ref::field(v1.at("mach")))) <= 5e-13);
    CHECK(s.passes.back().final);
    CHECK(std::isnan(s.passes.back().tan_correction));  // no correction at a fixed angle
    CHECK(s.passes.back().exit_angle_next == s.passes.back().exit_angle);
    CHECK(s.limitation == kLimitation);
    // every pass reported once, the final evaluation once and marked
    CHECK(counter.starts == 1);
    CHECK(counter.finishes == 1);
    CHECK(counter.converged);
    CHECK(counter.updates == static_cast<int>(s.passes.size()));
    CHECK(counter.finals == 1);
    CHECK(counter.records.back().final);
}

TEST_CASE("solver: 1 and 4 threads give bit-identical results (SparseLU)", "[turbolenza][solver][threads][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "fixed_exit_explicit_45x21");
    const FixedCase fc = fixed_case(d);
    if (max_threads() < 4) SKIP("fewer than 4 hardware threads");
    SolveOptions one, four;
    one.threads = 1;
    four.threads = 4;
    const Solution a = solve_at_exit_angle(fc.grid, fc.config, fc.c_ref, std::nullopt, one);
    const Solution b = solve_at_exit_angle(fc.grid, fc.config, fc.c_ref, std::nullopt, four);
    REQUIRE(a.passes.size() == b.passes.size());
    for (std::size_t k = 0; k < a.passes.size(); ++k) CHECK(same_record(a.passes[k], b.passes[k]));
    CHECK(same_bits(a.psi, b.psi));
    CHECK(same_bits(a.density, b.density));
    CHECK(same_bits(a.mach, b.mach));
    CHECK(same_bits(a.coefficients.D, b.coefficients.D));
    CHECK(same_bits(a.static_state.pressure, b.static_state.pressure));
    CHECK(same_bits(a.flow.W_theta, b.flow.W_theta));
    // the coupled pass too (two passes on the injected v1 grids)
    TURBOLENZA_REF_OR_SKIP(m, "m065_very_coarse_12");
    SolverConfig cfg = sref::config(m);
    cfg.iteration.max_density_iterations = 2;
    one.grid_ops = four.grid_ops = injected(m);
    std::shared_ptr<const Solution> x, y;
    try {
        (void)solve_coupled(sref::grid(m.at("first_grid")), cfg, 0.1, std::nullopt, one);
    } catch (const CoupledNotConvergedError& e) {
        x = e.solution;
    }
    try {
        (void)solve_coupled(sref::grid(m.at("first_grid")), cfg, 0.1, std::nullopt, four);
    } catch (const CoupledNotConvergedError& e) {
        y = e.solution;
    }
    REQUIRE((x && y));
    for (std::size_t k = 0; k < 2; ++k) CHECK(same_record(x->passes[k], y->passes[k]));
    CHECK(same_bits(x->psi, y->psi));
}

TEST_CASE("solver: PARDISO against SparseLU on the full fixed-angle solve", "[turbolenza][solver][pardiso][v1ref]")
{
    if (!pardiso_available()) SKIP("built without MKL PARDISO");
    TURBOLENZA_REF_OR_SKIP(d, "fixed_exit_explicit_45x21");
    const FixedCase fc = fixed_case(d);
    auto run = [&](LinearSolverKind kind, int threads, const char* name) {
        SolveOptions o;
        o.linear_solver = kind;
        o.threads = std::min(threads, max_threads());
        const auto t0 = std::chrono::steady_clock::now();
        Solution s = solve_at_exit_angle(fc.grid, fc.config, fc.c_ref, std::nullopt, o);
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        sref::measure(std::string("solver 45x21 time per pass [s] ") + name, t / static_cast<double>(s.passes.size()));
        return s;
    };
    const Solution lu = run(LinearSolverKind::SparseLU, 1, "SparseLU 1 thread");
    (void)run(LinearSolverKind::SparseLU, 4, "SparseLU 4 threads");
    const Solution p1 = run(LinearSolverKind::Pardiso, 1, "PARDISO 1 thread");
    const Solution p4 = run(LinearSolverKind::Pardiso, 4, "PARDISO 4 threads");
    for (const auto& [name, p] : {std::pair<const char*, const Solution*>{"1 thread", &p1}, {"4 threads", &p4}}) {
        CHECK(p->iterations == lu.iterations);
        CHECK(sref::measure(std::string("solver 45x21 PARDISO vs SparseLU psi ") + name, ref::rel_diff(p->psi, lu.psi)) <=
              1e-13);  // measured 9.8e-15
        CHECK(sref::measure(std::string("solver 45x21 PARDISO vs SparseLU density ") + name,
                            ref::rel_diff(p->density, lu.density)) <= 2e-13);  // measured 1.9e-14
        CHECK(sref::measure(std::string("solver 45x21 PARDISO vs SparseLU mach ") + name, ref::rel_diff(p->mach, lu.mach)) <=
              5e-13);  // measured 5.7e-14
        CHECK(sref::measure(std::string("solver 45x21 PARDISO vs v1 psi ") + name,
                            ref::rel_diff(p->psi, ref::field(d.at("solution").at("psi")))) <= 2e-13);  // measured 2.5e-14
    }
    sref::measure("solver 45x21 PARDISO 1 vs 4 threads psi", ref::rel_diff(p1.psi, p4.psi));
}

TEST_CASE("solver: observers never change the solve; cancel; errors carry the last pass", "[turbolenza][solver]")
{
    // a uniform channel: converges in a few passes at a fixed angle (no v1 data needed)
    const auto g = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, 0.3, 0.8));
    SolverConfig cfg;
    cfg.flow.inlet_mach = 0.4;
    cfg.flow.inlet_angle_deg = BladeGen::Turbolenza::detail::degrees(0.3);
    cfg.iteration.density_tolerance = 1e-12;
    cfg.iteration.max_density_iterations = 60;
    cfg.iteration.kutta_tolerance = 1e-12;
    cfg.iteration.max_kutta_iterations = 60;
    const Solution plain = solve_at_exit_angle(g, cfg, 0.1);
    CHECK(plain.converged);
    // uniform flow: psi = eta and the inlet Mach number everywhere
    for (int i = 0; i < 31; ++i)
        for (int j = 0; j < 9; ++j) {
            CHECK(std::abs(plain.psi(i, j) - g->eta[static_cast<std::size_t>(j)]) <= 1e-12);
            CHECK(std::abs(plain.mach(i, j) - 0.4) <= 1e-12);
        }
    Thrower thrower;
    SolveOptions o;
    o.observer = &thrower;
    const Solution watched = solve_at_exit_angle(g, cfg, 0.1, std::nullopt, o);
    CHECK(same_bits(watched.psi, plain.psi));
    CHECK(watched.passes.size() == plain.passes.size());

    // cancel after pass 1: ends like an unconverged run with the last completed pass
    struct Canceller final : RunObserver
    {
        std::atomic<bool>* flag;
        void update(const PassRecord& r) override
        {
            if (r.index == 1) flag->store(true);
        }
    };
    std::atomic<bool> flag{false};
    Canceller canceller;
    canceller.flag = &flag;
    SolveOptions oc;
    oc.observer = &canceller;
    oc.cancel = &flag;
    try {
        (void)solve_at_exit_angle(g, cfg, 0.1, std::nullopt, oc);
        FAIL("the run was cancelled");
    } catch (const CancelledError& e) {
        CHECK(error_type_name(e) == "Cancelled");
        REQUIRE(e.solution);
        CHECK(e.solution->passes.size() == 2);
        CHECK(!e.solution->converged);
    }

    // a called module fails in pass 1 (the exit rebuild moves a blade node): its type is kept and
    // the last completed pass is handed on
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    Grid moved = ref::grid(d.at("passes").at(1).at("grid"));
    moved.node_y(30, 4) += 1e-3;
    GridOps ops;
    ops.rebuild_exit = [moved](const Grid&, double) { return moved; };
    ops.validate_topology = [](const Grid&) {};
    SolveFailure failure;
    SolveOptions of;
    of.grid_ops = ops;
    of.failure = &failure;
    CHECK_THROWS_AS(solve_coupled(sref::grid(d.at("first_grid")), sref::config(d), 0.1, std::nullopt, of),
                    ExitRebuildError);
    CHECK(failure.failed_pass == 1);
    REQUIRE(failure.solution);
    CHECK(failure.solution->passes.size() == 1);
}

TEST_CASE("solver: refusals", "[turbolenza][solver]")
{
    const auto g = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, 0.3, 0.8));
    SolverConfig cfg;
    cfg.flow.inlet_mach = 0.4;
    cfg.flow.inlet_angle_deg = 17.0;
    cfg.iteration.density_tolerance = 1e-8;
    cfg.iteration.max_density_iterations = 10;
    CHECK_THROWS_AS(solve_at_exit_angle(g, cfg, 0.0), SolverInputError);
    SolveOptions o;
    o.threads = 0;
    CHECK_THROWS_AS(solve_at_exit_angle(g, cfg, 0.1, std::nullopt, o), SolverInputError);
    o.threads = kMaxThreads + 1;
    CHECK_THROWS_AS(solve_at_exit_angle(g, cfg, 0.1, std::nullopt, o), SolverInputError);
    o.threads = 1;
    o.initial_density = Field(3, 3, 1.0);
    CHECK_THROWS_AS(solve_at_exit_angle(g, cfg, 0.1, std::nullopt, o), SolverInputError);
    CHECK_THROWS_AS(solve_at_exit_angle(g, cfg, 0.1, 123.0), SolverInputError);  // contradicts M1
    SolverConfig none = cfg;
    none.flow.inlet_mach.reset();
    CHECK_THROWS_AS(solve_at_exit_angle(g, none, 0.1), SolverInputError);
    SolverConfig rot = cfg;
    rot.flow.angular_speed = 10.0;
    CHECK_THROWS_AS(solve_at_exit_angle(g, rot, 0.1), SolverInputError);
    SolverConfig zero = cfg;
    zero.iteration.max_density_iterations = 0;
    CHECK_THROWS_AS(solve_at_exit_angle(g, zero, 0.1), SolverInputError);
    // a mass flow above the channel's choking flux (about 19.3 kg/(m s) through 0.08 m): the
    // density loop accelerates the flow until M_rel reaches 1, and stops there
    SolverConfig choked = cfg;
    choked.flow.inlet_mach.reset();
    choked.iteration.max_density_iterations = 200;
    try {
        (void)solve_at_exit_angle(g, choked, 0.1, 25.0);
        FAIL("a choked channel has no subsonic solution");
    } catch (const TransonicFlowError& e) {
        REQUIRE(e.solution);
        CHECK(e.solution->passes.back().mach_max >= 1.0);  // the stopping pass is reported
    }
}

TEST_CASE("solver: the coupled M1 0.65 case, 12 passes (v1 records)", "[turbolenza][solver][v1ref][needs-grid]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_very_coarse_12");
    const SolverConfig cfg = sref::config(d);
    try {
        (void)solve_coupled(grid_with_curves(d), cfg, ref::num(d.at("settings").at("chord")));
        FAIL("12 passes do not converge");
    } catch (const CoupledNotConvergedError& e) {
        REQUIRE(e.solution);
        REQUIRE(e.solution->passes.size() == 12);
        double worst = 0.0;
        for (std::size_t k = 0; k < 12; ++k) worst = std::max(worst, sref::record_diff(e.solution->passes[k], d.at("records").at(k)));
        CHECK(sref::measure("solver m065 12 passes records worst", worst) <= 3e-12);  // measured 2.7e-13
    }
}

TEST_CASE("solver: the coupled M1 0.30 case converges as v1 (92 passes)", "[turbolenza][solver][v1ref][needs-grid]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m030_very_coarse_converged");
    const SolverConfig cfg = sref::config(d);
    const Solution s = solve_coupled(grid_with_curves(d), cfg, ref::num(d.at("settings").at("chord")));
    const ref::json& v1 = d.at("solution");
    CHECK(s.converged);
    CHECK(s.iterations == v1.at("iterations").get<int>());
    CHECK(sref::measure("solver m030 exit angle", sref::rel(s.exit_angle, ref::num(v1.at("exit_angle")))) <= 2e-13);  // measured 1.9e-14
    CHECK(sref::measure("solver m030 mach", ref::rel_diff(s.mach, ref::field(v1.at("mach")))) <= 4e-13);  // measured 3.5e-14
    CHECK(sref::measure("solver m030 psi", ref::rel_diff(s.psi, ref::field(v1.at("psi")))) <= 2e-13);  // measured 1.5e-14
}
