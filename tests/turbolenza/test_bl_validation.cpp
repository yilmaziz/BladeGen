// The validation run of the boundary layer (hidden, "[.validation]": a report, not a pass/fail test):
// the Dunavant case of Aungier (2003) Fig. 5-15 on v1's frozen NACA 65-(12A10)10 blade, medium grid,
// lambda_rho = lambda_K = 0.3, tolerances 1e-5, chord 0.127 m, T0 335 K, P0 from Re_ref 5.5e5 at M 0.30.
//   (a) one way, both laminar c_f forms;  (b) two way, both wake models.
// Each run prints "VALIDATION" lines (H against the digitized Fig. 5-15, transition and separation,
// beta_out, loss, the coupling history) and writes its files to BLADEGEN_TEST_OUTPUT/bl_validation/<run>.
//   build\<dir>\bin\test_turbolenza.exe "[.validation]"
#include "CascadeFixture.hpp"
#include "SolverRef.hpp"

#include "turbolenza/Case.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace BladeGen::Turbolenza;

namespace {

struct FigPoint
{
    double x_c, H;
    std::string regime;
};

std::vector<FigPoint> figure(const std::string& name)
{
    std::ifstream in(std::string(BLADEGEN_SOURCE_DIR) + "/tests/turbolenza/validation/aungier2003_fig5_15/" + name);
    std::vector<FigPoint> out;
    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string a, b, c;
        std::getline(ss, a, ',');
        std::getline(ss, b, ',');
        std::getline(ss, c, ',');
        if (!c.empty() && c.back() == '\r') c.pop_back();
        out.push_back({std::stod(a), std::stod(b), c});
    }
    return out;
}

// our H at x/c by linear interpolation between stations (NaN outside our stations)
double h_at(const BoundaryLayerSurfaceData& d, double x)
{
    for (std::size_t k = 1; k < d.x_c.size(); ++k)
        if (d.x_c[k - 1] <= x && x <= d.x_c[k]) {
            const double t = (x - d.x_c[k - 1]) / (d.x_c[k] - d.x_c[k - 1]);
            return (1.0 - t) * d.H[k - 1] + t * d.H[k];
        }
    return std::nan("");
}

struct Compare
{
    double max_turbulent = 0.0, mean_turbulent = 0.0, max_laminar = 0.0;
    int n_turbulent = 0, n_laminar = 0;
};

// against the figure's vertices: turbulent (and separated) ones downstream of both transitions, and the laminar ones
Compare compare(const BoundaryLayerSurfaceData& d, const std::vector<FigPoint>& fig)
{
    Compare c;
    double fig_transition = 0.0;
    for (const FigPoint& p : fig)
        if (p.regime == "transition") fig_transition = p.x_c;
    double sum = 0.0;
    for (const FigPoint& p : fig) {
        const double h = h_at(d, p.x_c);
        if (!std::isfinite(h)) continue;
        if ((p.regime == "turbulent" || p.regime == "separated") && p.x_c > std::max(fig_transition, d.transition_x_c)) {
            c.max_turbulent = std::max(c.max_turbulent, std::abs(h - p.H));
            sum += std::abs(h - p.H);
            ++c.n_turbulent;
        } else if (p.regime == "laminar" && (!std::isfinite(d.transition_x_c) || p.x_c < d.transition_x_c)) {
            c.max_laminar = std::max(c.max_laminar, std::abs(h - p.H));
            ++c.n_laminar;
        }
    }
    c.mean_turbulent = c.n_turbulent ? sum / c.n_turbulent : std::nan("");
    return c;
}

void report(const std::string& name, const RunResult& r, const Case& c)
{
    std::printf("VALIDATION ===== %s: status %s; %s\n", name.c_str(), r.status.c_str(), r.error.c_str());
    if (r.chord_reynolds) std::printf("VALIDATION %s chord Re %.6g\n", name.c_str(), *r.chord_reynolds);
    if (!r.viscous) return;
    const ViscousData& v = *r.viscous;
    for (const BlIteration& h : v.history)
        std::printf("VALIDATION %s k %2d passes %4d beta_out %.6f deg change %.3e tr SS %.4f PS %.4f sep SS %.4f PS %.4f "
                    "lamsep SS %.4f PS %.4f loss %.5f\n",
                    name.c_str(), h.k, h.inviscid_passes, h.beta_out * 180.0 / 3.141592653589793, h.max_change_m,
                    h.transition_x_c_ss, h.transition_x_c_ps, h.separation_x_c_ss, h.separation_x_c_ps,
                    h.laminar_separation_x_c_ss, h.laminar_separation_x_c_ps, h.loss);
    std::printf("VALIDATION %s beta_out inviscid %.6f deg, reported %.6f deg; loss omega %.5f (Theta %.5f, Delta* %.5f); "
                "TE jump mismatch %.3e\n",
                name.c_str(), v.inviscid_exit_angle * 180.0 / 3.141592653589793,
                r.solution ? r.solution->exit_angle * 180.0 / 3.141592653589793 : std::nan(""), v.loss_omega, v.loss_Theta,
                v.loss_Delta_star, v.te_jump_mismatch);
    if (v.suction.x_c.empty()) return;
    for (const auto& [label, d, fig] : {std::tuple{"SS", &v.suction, figure("H_suction.csv")},
                                        std::tuple{"PS", &v.pressure_side, figure("H_pressure.csv")}}) {
        const Compare cmp = compare(*d, fig);
        std::printf("VALIDATION %s %s transition x/c %.4f (fig %s), separation x/c %.4f, laminar separation x/c %.4f; "
                    "H vs Fig. 5-15: turbulent max |dH| %.3f mean %.3f (%d pts), laminar max |dH| %.3f (%d pts); "
                    "H at x/c 0.5 %.3f, TE %.3f\n",
                    name.c_str(), label, d->transition_x_c, label == std::string("SS") ? "0.130" : "0.173", d->separation_x_c,
                    d->laminar_separation_x_c, cmp.max_turbulent, cmp.mean_turbulent, cmp.n_turbulent, cmp.max_laminar,
                    cmp.n_laminar, h_at(*d, 0.5), d->H.back());
        for (std::size_t k = 0; k < d->x_c.size(); ++k)
            std::printf("VALIDATION-H %s %s %.5f %.5f %s Re_theta %.1f Lambda %.3f Hk %.3f M %.4f\n", name.c_str(), label,
                        d->x_c[k], d->H[k], d->regime[k].c_str(), d->Re_theta[k], d->Lambda[k], d->H_k[k],
                        d->W_e[k] / std::sqrt(1.4 * 287.05 * d->T_e[k]));
    }
    const auto dir = std::filesystem::path(BLADEGEN_TEST_OUTPUT) / "bl_validation" / name;
    std::filesystem::create_directories(dir);
    (void)export_run(r, c, dir.string(), true);
}

}  // namespace

TEST_CASE("validation: Aungier (2003) Fig. 5-15, the Dunavant case on the medium grid", "[.validation][turbolenza]")
{
    TURBOLENZA_CASCADE_OR_SKIP(ref);
    const BladeSection section = fixture::v1_section(ref->family.at("cascade"));
    CaseSettings s;
    s.case_name = "dunavant_fig5_15";
    s.chord = 0.127;
    s.solidity = 1.0;
    s.stagger_deg = 33.0;
    s.inlet_flow_angle_deg = 45.0;
    s.inlet_mach = 0.65;
    s.inlet_total_temperature = 335.0;
    s.inlet_total_pressure = total_pressure_for_reynolds(550000.0, 0.30, 335.0, 0.127, s.gamma, s.gas_constant, s.sutherland);
    std::printf("VALIDATION P0 = %.6f Pa (Re_ref 5.5e5 at M 0.30, T0 335 K, c 0.127 m); Re at M 0.65 = %.6g\n",
                s.inlet_total_pressure,
                chord_reynolds_number(0.65, s.inlet_total_pressure, 335.0, 0.127, s.gamma, s.gas_constant, s.sutherland));
    s.grid_level = std::string("medium");
    s.density_relaxation = 0.3;
    s.kutta_relaxation = 0.3;
    s.density_tolerance = 1e-5;
    s.kutta_tolerance = 1e-5;
    s.viscous = true;
    struct Run
    {
        const char* name;
        BlCoupling coupling;
        LaminarCfForm form;
        WakeModel wake;
        double omega;
        int max_iterations;
    };
    const Run runs[] = {
        {"one_way_centrifugal", BlCoupling::OneWay, LaminarCfForm::Centrifugal2000, WakeModel::Constant, 0.5, 30},
        {"one_way_axial_printed", BlCoupling::OneWay, LaminarCfForm::Axial2003AsPrinted, WakeModel::Constant, 0.5, 30},
        {"two_way_constant", BlCoupling::TwoWay, LaminarCfForm::Centrifugal2000, WakeModel::Constant, 0.5, 30},
        {"two_way_marched", BlCoupling::TwoWay, LaminarCfForm::Centrifugal2000, WakeModel::Marched, 0.5, 30},
        // a smaller omega_BL and more iterations (settings, not a new damping): measured 2026-11-24, omega 0.5
        // diverges (sonic at the SS TE node in BL iteration 2), 0.2 contracts by about 0.77 per iteration,
        // 0.1 by about 0.88
        {"two_way_constant_omega0.2", BlCoupling::TwoWay, LaminarCfForm::Centrifugal2000, WakeModel::Constant, 0.2, 60},
        {"two_way_marched_omega0.2", BlCoupling::TwoWay, LaminarCfForm::Centrifugal2000, WakeModel::Marched, 0.2, 60},
    };
    for (const Run& run : runs) {
        CaseSettings x = s;
        x.case_name = run.name;
        x.bl_coupling = run.coupling;
        x.laminar_cf_form = run.form;
        x.wake_model = run.wake;
        x.bl_relaxation = run.omega;
        x.max_bl_iterations = run.max_iterations;
        const Case c = build_case(x, section);
        const auto t0 = std::chrono::steady_clock::now();
        const RunResult r = run_case(c);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("VALIDATION %s wall time %.1f s, %zu passes\n", run.name, sec, r.records.size());
        report(run.name, r, c);
        CHECK(r.viscous);
    }
}
