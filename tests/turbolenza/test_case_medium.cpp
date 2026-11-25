// The user's production case end to end against v1: NACA 65-(12A10)10 at M1 0.65, beta1 45 deg, the
// medium standard grid, lambda_rho = lambda_K = 0.3, tolerances 1e-5 (scripts/run_full_solver.py of
// 2026-11-24), on v1's frozen C2 blade curves.  Reference: m065_medium_user.json (dump_v1.py --medium).
// Also the same case with PARDISO on 4 threads (agreement to round-off, not bit for bit).
#include "CascadeFixture.hpp"
#include "SolverRef.hpp"

#include "turbolenza/Case.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>

using namespace BladeGen::Turbolenza;

namespace {

CaseSettings user_settings()
{
    CaseSettings c;            // the defaults are the user's run_full_solver.py values
    c.case_name = "m065_medium_user";
    c.chord = 0.5;
    c.grid_level = std::string("medium");
    c.density_tolerance = 1e-5;
    c.kutta_tolerance = 1e-5;
    c.density_relaxation = 0.3;
    c.kutta_relaxation = 0.3;
    c.max_density_iterations = 5000;
    c.max_kutta_iterations = 6000;
    return c;
}

}  // namespace

TEST_CASE("case: the user's M1 0.65 medium case matches v1 end to end", "[turbolenza][case][v1ref][needs-grid][slow]")
{
    TURBOLENZA_REF_OR_SKIP(d, "m065_medium_user");
    const Case c = build_case(user_settings(), fixture::v1_section(d.at("cascade")));
    const auto t0 = std::chrono::steady_clock::now();
    const RunResult r = run_case(c);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const bool v1_converged = d.at("error").is_null();
    CHECK(r.converged == v1_converged);
    REQUIRE(r.records.size() == d.at("records").size());
    double worst = 0.0;
    for (std::size_t k = 0; k < r.records.size(); ++k) worst = std::max(worst, sref::record_diff(r.records[k], d.at("records").at(k)));
    WARN("medium: " << r.records.size() << " passes in " << seconds << " s (" << 1e3 * seconds / static_cast<double>(r.records.size())
                    << " ms/pass, SparseLU, 1 thread); worst record difference " << worst);
    CHECK(sref::measure("case m065 medium records worst", worst) <= 1e-11);  // measured 1.1e-12 (254 passes)
    REQUIRE(r.solution);
    REQUIRE(d.contains("solution"));
    const double beta = r.solution->exit_angle, beta_v1 = ref::num(d.at("solution").at("exit_angle"));
    WARN("medium: beta_out " << beta * 180.0 / 3.141592653589793 << " deg (v1 " << beta_v1 * 180.0 / 3.141592653589793 << ")");
    CHECK(sref::measure("case m065 medium exit angle", sref::rel(beta, beta_v1)) <= 3e-13);  // measured 2.7e-14
    CHECK(sref::measure("case m065 medium surface mach ss", ref::rel_diff(r.solution->suction.mach, ref::vec(d.at("surfaces").at("mach_ss")))) <=
          2e-12);  // measured 1.7e-13
    CHECK(sref::measure("case m065 medium surface mach ps",
                        ref::rel_diff(r.solution->pressure_side.mach, ref::vec(d.at("surfaces").at("mach_ps")))) <= 2e-12);  // measured 1.7e-13

    // PARDISO, 4 threads: the same fixed point to round-off
    if (pardiso_available()) {
        CaseSettings p = user_settings();
        p.linear_solver = LinearSolverKind::Pardiso;
        p.threads = 4;
        const Case cp = build_case(p, fixture::v1_section(d.at("cascade")));
        const auto t1 = std::chrono::steady_clock::now();
        const RunResult rp = run_case(cp);
        const double sp = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
        REQUIRE(rp.solution);
        CHECK(rp.converged == r.converged);
        WARN("medium PARDISO x4: " << rp.records.size() << " passes in " << sp << " s ("
                                   << 1e3 * sp / static_cast<double>(rp.records.size()) << " ms/pass); beta_out "
                                   << rp.solution->exit_angle * 180.0 / 3.141592653589793 << " deg");
        CHECK(sref::measure("case m065 medium PARDISO exit angle vs SparseLU", sref::rel(rp.solution->exit_angle, beta)) <= 3e-13);  // measured 3.0e-14
    }
}
