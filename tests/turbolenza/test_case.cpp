// Case (port of v1 controller.py): build_case from the v1 blade section (the grid port builds the
// grid), run_case end to end (fixed angle 45 x 21, coupled M1 0.65 and 0.30) against v1's records,
// the reports, export_run, the refusals, format_record.
#include "CascadeFixture.hpp"
#include "SolverRef.hpp"

#include "turbolenza/Case.hpp"
#include "turbolenza/Detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <thread>

using namespace BladeGen::Turbolenza;

namespace {

/// CaseSettings of a v1 dump (dump_v1.py settings()).
CaseSettings settings_of(const ref::json& d)
{
    const ref::json& s = d.at("settings");
    CaseSettings c;
    c.case_name = s.at("case_name").get<std::string>();
    c.solver_mode = s.at("solver_mode").get<std::string>() == "fixed_exit" ? SolverMode::FixedExit : SolverMode::CoupledKutta;
    c.chord = ref::num(s.at("chord"));
    c.solidity = ref::num(s.at("solidity"));
    c.stagger_deg = ref::num(s.at("stagger_deg"));
    c.inlet_mach = ref::num(s.at("inlet_mach"));
    c.inlet_flow_angle_deg = ref::num(s.at("inlet_flow_angle_deg"));
    if (s.at("grid_level").is_null()) {
        c.grid_level.reset();
        c.n_passage = s.at("n_passage").get<int>();
        c.n_eta = s.at("n_eta").get<int>();
    } else {
        c.grid_level = s.at("grid_level").get<std::string>();
    }
    c.density_tolerance = ref::num(s.at("density_tolerance"));
    c.max_density_iterations = s.at("max_density_iterations").get<int>();
    c.density_relaxation = ref::num(s.at("density_relaxation"));
    c.kutta_tolerance = ref::num(s.at("kutta_tolerance"));
    c.max_kutta_iterations = s.at("max_kutta_iterations").get<int>();
    c.kutta_relaxation = ref::num(s.at("kutta_relaxation"));
    return c;
}

double grid_diff(const Grid& a, const Grid& b)
{
    REQUIRE(a.node_m.same_shape(b.node_m));
    return std::max(ref::rel_diff(a.node_m, b.node_m), ref::rel_diff(a.node_y, b.node_y));
}

}  // namespace

TEST_CASE("case: format_record is v1's console line", "[turbolenza][case]")
{
    PassRecord r;
    r.index = 34;
    r.density_change = 2.1e-5;
    r.velocity_change = 3.0e-5;
    r.kutta_normalized = 8.4e-4;
    r.exit_angle = BladeGen::Turbolenza::detail::radians(23.21);
    r.mach_max = 0.401;
    r.linear_residual = 6.0e-16;
    CHECK(format_record(r) ==
          "iter    34 | rho 2.1e-05 | W 3.0e-05 | Kutta 8.4e-04 | beta2 23.2100 deg | Mmax 0.401 | lin 6.0e-16");
    r.final = true;
    r.velocity_change = std::numeric_limits<double>::quiet_NaN();
    CHECK(format_record(r).starts_with("iter final | rho 2.1e-05 | W nan |"));
    CHECK(BladeGen::Turbolenza::detail::py_repr(0.1) == "0.1");
    CHECK(BladeGen::Turbolenza::detail::py_repr(1e-05) == "1e-05");
    CHECK(BladeGen::Turbolenza::detail::py_repr(101325.0) == "101325.0");
    CHECK(BladeGen::Turbolenza::detail::py_repr(1e16) == "1e+16");
    CHECK(BladeGen::Turbolenza::detail::py_repr(-2.5e-310) == "-2.5e-310");
    CHECK(max_threads() >= 1);
    CHECK(max_threads() <= kMaxThreads);
}

TEST_CASE("case: refusals before any solve", "[turbolenza][case]")
{
    const BladeSection section = fixture::synthetic_section();
    CaseSettings s;
    s.threads = kMaxThreads + 1;
    CHECK_THROWS_AS(build_case(s, section), CaseError);
    s.threads = 1;
    s.chord = -1.0;
    CHECK_THROWS_AS(build_case(s, section), CaseError);
    s.chord = 0.1;
    s.grid_level = "coarse";
    s.n_passage = 40;
    CHECK_THROWS_AS(build_case(s, section), CaseError);
    s.n_passage.reset();
    s.gamma = 0.9;
    CHECK_THROWS_AS(build_case(s, section), CaseError);
    s.gamma = 1.4;
    s.solver_mode = SolverMode::CoupledKutta;
    s.initial_exit_angle_mode = ExitAngleStart::Manual;
    CHECK_THROWS_AS(build_case(s, section), CaseError);
    if (!pardiso_available()) {
        s.initial_exit_angle_mode = ExitAngleStart::BladeTeAngle;
        s.linear_solver = LinearSolverKind::Pardiso;
        CHECK_THROWS_AS(build_case(s, section), CaseError);
    }
}

TEST_CASE("case: build_case, run_case and export_run (fixed angle 45 x 21)", "[turbolenza][case][v1ref][needs-grid]")
{
    TURBOLENZA_REF_OR_SKIP(d, "fixed_exit_explicit_45x21");
    const CaseSettings s = settings_of(d);
    const Case c = build_case(s, fixture::v1_section(d.at("cascade")));
    CHECK(sref::measure("case 45x21 grid vs v1", grid_diff(c.grid, ref::grid(d.at("first_grid")))) <= 1e-13);
    CHECK(sref::measure("case 45x21 mdot vs v1",
                        sref::rel(c.mass_flow_per_depth, ref::num(d.at("solution").at("mass_flow_per_depth")))) <= 1e-14);
    CHECK(!format_case_summary(c).empty());
    CHECK(format_grid_summary(c).find("n_m") != std::string::npos);
    struct Watch final : RunObserver
    {
        int updates = 0;
        void update(const PassRecord&) override { ++updates; }
    } watch;
    const RunResult r = run_case(c, &watch);
    CHECK(r.error.empty());
    CHECK(r.converged);
    CHECK(r.status == "converged");
    CHECK(r.iterations == d.at("solution").at("iterations").get<int>());
    REQUIRE(r.records.size() == d.at("records").size());
    CHECK(watch.updates == static_cast<int>(r.records.size()));
    double worst = 0.0;
    for (std::size_t k = 0; k < r.records.size(); ++k) worst = std::max(worst, sref::record_diff(r.records[k], d.at("records").at(k)));
    CHECK(sref::measure("case 45x21 records worst", worst) <= 1e-11);  // measured 1.1e-12
    REQUIRE(r.solution);
    REQUIRE(r.summary);
    CHECK(sref::measure("case 45x21 surface mach ss",
                        ref::rel_diff(r.solution->suction.mach, ref::vec(d.at("surfaces").at("mach_ss")))) <= 1e-12);  // measured 8.1e-14
    CHECK(r.runtime.wall_time_s > 0.0);
    sref::measure("case 45x21 wall time [s]", r.runtime.wall_time_s);
    const std::string report = format_run_report(r);
    CHECK(report.find("SOLVER RESULT") != std::string::npos);
    const auto dir = std::filesystem::temp_directory_path() / "turbolenza_test_export";
    std::filesystem::remove_all(dir);
    const std::vector<std::string> files = export_run(r, c, dir.string());
    CHECK(files.size() == 5);
    for (const std::string& f : files) CHECK(std::filesystem::exists(dir / f));
    CHECK_THROWS(export_run(r, c, dir.string()));  // existing files are kept unless asked
    CHECK_NOTHROW(export_run(r, c, dir.string(), true));
    {
        std::ifstream in(dir / "run_summary.json");
        const ref::json summary = ref::json::parse(in);
        CHECK(summary.at("status").get<std::string>() == "converged");
        CHECK(summary.at("results").at("iterations").get<int>() == r.iterations);
        std::ifstream fields(dir / "fields.csv");
        int lines = 0;
        for (std::string line; std::getline(fields, line);) ++lines;
        CHECK(lines == 1 + c.grid.n_m() * c.grid.n_eta());
    }
    std::filesystem::remove_all(dir);
}

TEST_CASE("case: the coupled M1 0.65 and 0.30 cases through build_case and run_case", "[turbolenza][case][v1ref][needs-grid]")
{
    TURBOLENZA_REF_OR_SKIP(d65, "m065_very_coarse_12");
    {
        const Case c = build_case(settings_of(d65), fixture::v1_section(d65.at("cascade")));
        CHECK(sref::measure("case m065 first grid vs v1", grid_diff(c.grid, ref::grid(d65.at("first_grid")))) <= 1e-13);
        REQUIRE(c.standard_info);
        const RunResult r = run_case(c);
        CHECK(r.error.starts_with("CoupledNotConvergedError"));
        REQUIRE(r.records.size() == 12);
        double worst = 0.0;
        for (std::size_t k = 0; k < 12; ++k) worst = std::max(worst, sref::record_diff(r.records[k], d65.at("records").at(k)));
        CHECK(sref::measure("case m065 12 records worst", worst) <= 2e-12);  // measured 1.8e-13
        CHECK(format_run_report(r).find("SOLVER FAILURE") != std::string::npos);
    }
    TURBOLENZA_REF_OR_SKIP(d30, "m030_very_coarse_converged");
    const Case c = build_case(settings_of(d30), fixture::v1_section(d30.at("cascade")));
    const RunResult r = run_case(c);
    CHECK(r.converged);
    CHECK(r.iterations == d30.at("solution").at("iterations").get<int>());
    REQUIRE(r.records.size() == d30.at("records").size());
    double worst = 0.0;
    for (std::size_t k = 0; k < r.records.size(); ++k) worst = std::max(worst, sref::record_diff(r.records[k], d30.at("records").at(k)));
    CHECK(sref::measure("case m030 records worst", worst) <= 5e-12);  // measured 4.9e-13
    REQUIRE(r.solution);
    CHECK(sref::measure("case m030 exit angle", sref::rel(r.solution->exit_angle, ref::num(d30.at("solution").at("exit_angle")))) <=
          3e-13);  // measured 3.0e-14
    CHECK(sref::measure("case m030 surface mach ss",
                        ref::rel_diff(r.solution->suction.mach, ref::vec(d30.at("surfaces").at("mach_ss")))) <= 2e-13);  // measured 1.3e-14
    CHECK(sref::measure("case m030 Dunavant S ss",
                        ref::rel_diff(r.solution->suction.dunavant_S,
                                      ref::vec(d30.at("surface_pressure_coefficient").at("values_ss")))) <= 3e-13);  // measured 2.5e-14
    // cancel: an observer sets the flag after pass 3; the run ends like an unconverged one
    std::atomic<bool> flag{false};
    struct Stop final : RunObserver
    {
        std::atomic<bool>* f;
        void update(const PassRecord& p) override
        {
            if (p.index == 3) f->store(true);
        }
    } stop;
    stop.f = &flag;
    const RunResult cancelled = run_case(c, &stop, &flag);
    CHECK(cancelled.error.starts_with("Cancelled: "));
    CHECK(cancelled.records.size() == 4);
    CHECK(cancelled.solution);
    CHECK(!cancelled.converged);
}
