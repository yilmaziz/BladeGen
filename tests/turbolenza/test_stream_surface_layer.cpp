// A BladeGen layer on its stream surface (BladeGenSource.hpp stream_surface_blade, not in v1): the
// placed section in the solver's frame, BladeGen's blade angles, the pitch, r(m) and b(m); a cylindrical
// channel against the planar cascade of the same pitch and chord; and (hidden, [.rotor]) one rotating
// case on the default design's Mid layer.
#include "SolverRef.hpp"

#include "turbolenza/BladeGenSource.hpp"
#include "turbolenza/Case.hpp"
#include "turbolenza/Detail.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>

namespace T = BladeGen::Turbolenza;
namespace D = BladeGen::Design;

namespace {

constexpr double kPi = std::numbers::pi;

T::CaseSettings layer_settings()
{
    T::CaseSettings s;
    s.stream_surface = T::StreamSurfaceKind::Layer;
    s.grid_level = std::string("very_coarse");
    s.density_tolerance = 1e-6;
    s.kutta_tolerance = 1e-6;
    s.density_relaxation = 0.5;
    s.kutta_relaxation = 0.5;
    s.max_density_iterations = 3000;
    s.max_kutta_iterations = 3000;
    s.inlet_total_pressure = 101325.0;
    s.inlet_total_temperature = 288.15;
    return s;
}

}  // namespace

TEST_CASE("a BladeGen layer on its stream surface: the placed section, angles, pitch, r(m), b(m)", "[turbolenza][bladegen][sheet]")
{
    D::BladeDesign design = D::BladeDesign::defaults();
    D::Evaluator ev;
    const int k = 2;  // Mid
    T::StreamSurfaceReport rep;
    const T::StreamSurfaceBlade b = T::stream_surface_blade(design, k, design.layers[2], ev, {}, &rep);
    const D::DesignResult res = ev.evaluate(design);
    const D::LayerResult& L = res.layers[2];
    REQUIRE(L.error.empty());
    // the section: LE (0, 0), TE (1, 0), fitted within the fit's tolerance
    CHECK(std::hypot(b.section.suction.evaluate(0.0).m, b.section.suction.evaluate(0.0).y) <= 1e-15);
    CHECK(std::abs(b.section.suction.evaluate(1.0).m - 1.0) <= 1e-15);
    CHECK(b.max_fit_residual <= 1e-7);
    // our exact map against BladeGen's StreamSurface table (Hermite, 128 intervals per piece)
    CHECK(sref::measure("layer map m(m') vs BladeGen's StreamSurface [chords]", rep.map_max_difference) <= 1e-8);
    // the true blade angles: kappa_le, kappa_te are BladeGen's metal angles beta_LE, beta_TE (custom layers)
    const T::Cascade c = T::make_cascade(b);
    CHECK(std::abs(T::detail::degrees(c.kappa_le) - L.beta_le) <= 1e-9);
    CHECK(std::abs(T::detail::degrees(c.kappa_te) - L.beta_te) <= 1e-9);
    // the pitch r_LE 2 pi / N, the radius at the LE r_LE (sigma = 1 there)
    CHECK(std::abs(b.r_le - 1e-3 * L.r_le) <= 1e-12 * b.r_le);
    // the sheet's cubic spline through its stations (step 0.0025 chords) against the exact streamline
    CHECK(sref::measure("Mid layer sheet r(LE) vs the streamline (relative)",
                        std::abs(b.sheet->radius(0.0) * b.chord / b.r_le - 1.0)) <= 1e-9);
    CHECK(std::abs(b.sheet->scale(0.0) - 1.0) <= 1e-9);
    CHECK(std::abs(1.0 / b.solidity - b.sheet->r_ref() * 2.0 * kPi / design.n_blades) <= 1e-12);
    CHECK(b.n_blades == design.n_blades);
    // r at the TE point (the sheet at the TE's m = cos(stagger) chords) is BladeGen's r_TE
    const double m_te = std::cos(T::detail::radians(b.stagger_deg));
    CHECK(sref::measure("Mid layer sheet r(TE) vs BladeGen r_TE (relative)", std::abs(b.sheet->radius(m_te) * b.chord / (1e-3 * L.r_te) - 1.0)) <= 1e-9);
    // b(m) of the channel family: positive, smooth, b_te / b_le reported
    CHECK(rep.b_te > 0.0);
    sref::measure("Mid layer b_TE / b_LE", rep.b_te);
    sref::measure("Mid layer r_TE / r_LE", b.r_te / b.r_le);
    CHECK_THROWS_AS(T::stream_surface_blade(design, 7, design.layers[2], ev), std::invalid_argument);
    D::BladeDesign mirrored = design;
    mirrored.mirror = true;
    CHECK_THROWS_AS(T::stream_surface_blade(mirrored, k, design.layers[2], ev), std::invalid_argument);
}

TEST_CASE("a cylindrical channel: the layer on its stream surface = the planar cascade of the same pitch and chord", "[turbolenza][bladegen][sheet]")
{
    // straight cylindrical walls: r = r0 on the Mid streamline, b = casing - hub radius, constant
    D::BladeDesign design = D::BladeDesign::naca65_defaults();
    design.meridional.hub = {-15.0, 60.0, 15.0, 60.0};
    design.meridional.casing = {-15.0, 100.0, 15.0, 100.0};
    design.meridional.hub_law = design.meridional.casing_law = BladeGen::Blade::WallLaw::Linear;
    design.n_blades = 30;
    D::Evaluator ev;
    const D::Layer& layer = design.layers[2];
    const T::StreamSurfaceBlade b = T::stream_surface_blade(design, 2, layer, ev);
    const double r0 = 80.0;  // mid span
    for (double m : {-3.0, 0.0, 0.7, 4.0}) {
        CHECK(std::abs(b.sheet->radius(m) * b.chord / 1e-3 - r0) <= 1e-11 * r0);
        CHECK(std::abs(b.sheet->sin_phi(m)) <= 1e-12);
        CHECK(std::abs(b.sheet->thickness(m) - 40.0) <= 1e-10);
    }
    const D::DesignResult res = ev.evaluate(design);
    const D::LayerResult& L = res.layers[2];
    // the planar cascade: chord r0 c', pitch 2 pi r0 / N, BladeGen's stagger
    T::CaseSettings planar = layer_settings();
    planar.stream_surface = T::StreamSurfaceKind::Planar;
    planar.chord = 1e-3 * r0 * L.chord_mprime;
    planar.solidity = L.chord_mprime * design.n_blades / (2.0 * kPi);
    planar.stagger_deg = T::detail::degrees(L.gamma);
    planar.inlet_mach = 0.5;
    planar.inlet_flow_angle_deg = planar.stagger_deg + 10.0;
    T::CaseSettings onsheet = planar;
    onsheet.stream_surface = T::StreamSurfaceKind::Layer;
    const T::Case cp = T::build_case(planar, T::section_from_layer(layer, ev));
    const T::Case cs = T::build_case(onsheet, b);
    CHECK(sref::measure("cylinder: chord vs r0 c' (relative)", std::abs(cs.settings.chord / planar.chord - 1.0)) <= 1e-12);
    CHECK(sref::measure("cylinder: stagger vs BladeGen's [deg]", std::abs(cs.settings.stagger_deg - planar.stagger_deg)) <= 1e-10);
    CHECK(std::abs(cs.settings.solidity / planar.solidity - 1.0) <= 1e-12);
    CHECK(std::abs(cs.cascade.kappa_le - cp.cascade.kappa_le) <= 1e-12);
    CHECK(std::abs(cs.cascade.kappa_te - cp.cascade.kappa_te) <= 1e-12);
    CHECK(std::abs(cs.mass_flow_per_depth / cp.mass_flow_per_depth - 1.0) <= 1e-12);
    const T::RunResult a = T::run_case(cp), r = T::run_case(cs);
    REQUIRE(a.converged);
    REQUIRE(r.converged);
    const double d_beta = std::abs(a.solution->exit_angle - r.solution->exit_angle);
    CHECK(sref::measure("cylinder layer vs planar cascade: beta_out [rad]", d_beta) <= 1e-8);
    double dm = 0.0;
    for (std::size_t q = 0; q < a.solution->suction.mach.size() && q < r.solution->suction.mach.size(); ++q)
        dm = std::max(dm, std::abs(a.solution->suction.mach[q] - r.solution->suction.mach[q]));
    CHECK(a.solution->suction.mach.size() == r.solution->suction.mach.size());
    CHECK(sref::measure("cylinder layer vs planar cascade: SS Mach", dm) <= 1e-7);
}

TEST_CASE("a conical channel: r(m), sin(phi), b(m) in closed form; the mapped curvature by the chain rule", "[turbolenza][bladegen][sheet]")
{
    // straight walls: the hub a cylinder r = 40, the casing a cone 100 -> 120 mm over the same axial length 30 mm
    // (Linear law, continued straight).  BladeGen's family X(u, s) = (1 - s) hub(u) + s casing(u) is then, with
    // tau = u - 1 (the LE line at tau = 0): z = -15 + 30 tau, r = 40 + s (60 + 20 tau), so
    //   X_u = (30, 20 s), m = tau |X_u|, sin(phi) = 20 s / |X_u|, dX/ds = (0, 60 + 20 tau),
    //   b = |dX/ds . n| = 30 (60 + 20 tau) / |X_u|   (n = (-20 s, 30) / |X_u|; |dX/ds| would be 60 + 20 tau)
    D::BladeDesign design = D::BladeDesign::naca65_defaults();
    design.meridional.hub = {-15.0, 40.0, 15.0, 40.0};
    design.meridional.casing = {-15.0, 100.0, 15.0, 120.0};
    design.meridional.hub_law = design.meridional.casing_law = BladeGen::Blade::WallLaw::Linear;
    design.n_blades = 30;
    D::Evaluator ev;
    const D::Layer& layer = design.layers[2];
    const double s = layer.span, speed = std::sqrt(900.0 + 400.0 * s * s);
    T::StreamSurfaceOptions opt;
    opt.keep_samples = true;
    opt.coverage_chords = 2.0;  // the hub line r = 70 + 10 tau of mid span reaches r = 0 about 5 chords upstream
    T::StreamSurfaceReport rep;
    const T::StreamSurfaceBlade b = T::stream_surface_blade(design, 2, layer, ev, opt, &rep);
    const double L = b.chord / 1e-3;  // the normalizing length [mm]
    double wr = 0.0, wb = 0.0, wp = 0.0;
    for (double mc : {-1.9, -1.5, -0.3, 0.0, 0.4, 0.8, 1.7, 2.6}) {
        const double tau = mc * L / speed;
        wr = std::max(wr, std::abs(b.sheet->radius(mc) * L / (40.0 + s * (60.0 + 20.0 * tau)) - 1.0));
        wb = std::max(wb, std::abs(b.sheet->thickness(mc) / (30.0 * (60.0 + 20.0 * tau) / speed) - 1.0));
        wp = std::max(wp, std::abs(b.sheet->sin_phi(mc) * L / L - 20.0 * s / speed));
    }
    CHECK(sref::measure("conical channel: r(m) vs closed form (relative)", wr) <= 1e-10);
    CHECK(sref::measure("conical channel: b(m) vs 30 (60 + 20 tau) / |X_u| (relative)", wb) <= 1e-10);
    CHECK(sref::measure("conical channel: sin(phi) vs 20 s / |X_u|", wp) <= 1e-10);
    // the mapped sides' curvature (the chain rule of the map, d2M/ds2 = r d2m'/ds2 + r sin(phi) (dm'/ds)^2) against
    // the curvature of the circle through three mapped points (independent: the positions only), on the body
    for (const T::SideSample* side : {&rep.suction_sample, &rep.pressure_sample}) {
        REQUIRE(side->xy.size() > 1000);
        const std::size_t n = side->xy.size(), k = 8;
        double worst = 0.0, kmax = 0.0;
        int used = 0;
        for (std::size_t q = n / 10; q < n - n / 10; q += 37) {
            // where BladeGen's curvature is smooth over the three points (not at a fair's join, where a circle
            // through three points straddles a jump in d kappa / ds)
            const double scale = std::max(std::abs(side->kappa[q]), 1.0);
            if (std::abs(side->kappa[q + k] - side->kappa[q - k]) > 1e-3 * scale) continue;
            const T::Point A = side->xy[q - k], B = side->xy[q], C = side->xy[q + k];
            const double ux = B.m - A.m, uy = B.y - A.y, vx = C.m - B.m, vy = C.y - B.y;
            const double cross = ux * vy - uy * vx;
            const double kfd = 2.0 * cross / (std::hypot(ux, uy) * std::hypot(vx, vy) * std::hypot(C.m - A.m, C.y - A.y));
            worst = std::max(worst, std::abs(kfd - side->kappa[q]) / scale);
            ++used;
            kmax = std::max(kmax, std::abs(side->kappa[q]));
        }
        // measured 2026-11-24: 3.7e-3 / 5.2e-4 here and 3.9e-3 / 1.8e-3 on the cylinder (dr/dm' = 0): the floor of
        // the three-point circle against BladeGen's own sample curvature, not the map (mapped_curvature is tested
        // against closed forms below)
        CHECK(sref::measure("conical channel: chain-rule curvature vs three-point circle", worst) <= 1e-2);
        CHECK(kmax > 0.0);
        CHECK(used > 200);
    }
}

TEST_CASE("mapped_curvature: the chain rule of the map against closed forms", "[turbolenza][bladegen][sheet]")
{
    const double g = 0.5, c = 0.4, r = 70.0, sphi = 0.3, r_le = 65.0, L = 30.0;
    // a straight chord-frame line along x (kappa = 0): in (m', theta) the line theta - t0 = tan(g) (m' - m0); mapped,
    // M(m') = m(m'), Y = r_LE tan(g) (m' - m0), so with M' = r, M'' = dr/dm' = r sin(phi), Y'' = 0 per unit m',
    // kappa = -L Y' M'' / (M'^2 + Y'^2)^(3/2), Y' = r_LE tan(g).  Without the dr/dm' term it would be 0.
    {
        const double Yp = r_le * std::tan(g), Mp = r, Mpp = r * sphi;
        const double exact = -L * Yp * Mpp / std::pow(Mp * Mp + Yp * Yp, 1.5);
        const double got = T::mapped_curvature({1.0, 0.0}, 0.0, g, c, r, r * sphi, r_le, L);
        CHECK(std::abs(got / exact - 1.0) <= 1e-13);
        CHECK(exact < 0.0);  // r grows along m: the mapped line bends toward -Y
    }
    // a circle of curvature k0 on a cylinder at r = r_LE (dr/dm' = 0): the map is the similarity of scale c' r / L,
    // so kappa = k0 L / (c' r), whatever the tangent
    for (double ang : {0.0, 0.7, 2.0, -1.2}) {
        const T::Point t{std::cos(ang), std::sin(ang)};
        const double got = T::mapped_curvature(t, 3.0, g, c, r_le, 0.0, r_le, L);
        CHECK(std::abs(got / (3.0 * L / (c * r_le)) - 1.0) <= 1e-13);
    }
    // a general point by finite differences of the map of the osculating circle, on a cone r(m') = r0 exp(sin(phi) m')
    // (dr/dm' = r sin(phi): dm = r dm', dr = sin(phi) dm), m(m') = (r(m') - r0) / sin(phi)
    {
        const double r0 = 60.0, k0 = 2.5, ang = 0.4, m0p = 0.3;
        auto rad = [&](double mp) { return r0 * std::exp(sphi * mp); };
        auto mapped = [&](double s) {  // the circle through the chord-frame origin with tangent ang, curvature k0
            const double x = (std::sin(ang + k0 * s) - std::sin(ang)) / k0, y = (std::cos(ang) - std::cos(ang + k0 * s)) / k0;
            const double mp = m0p + c * (x * std::cos(g) - y * std::sin(g)), th = c * (x * std::sin(g) + y * std::cos(g));
            return T::Point{(rad(mp) - r0) / sphi / L, r_le * th / L};
        };
        const double h = 1e-4;
        const T::Point A = mapped(-h), B = mapped(0.0), C = mapped(h);
        const double d1m = (C.m - A.m) / (2 * h), d1y = (C.y - A.y) / (2 * h);
        const double d2m = (C.m - 2 * B.m + A.m) / (h * h), d2y = (C.y - 2 * B.y + A.y) / (h * h);
        const double fd = (d1m * d2y - d1y * d2m) / std::pow(d1m * d1m + d1y * d1y, 1.5);
        const double rr = rad(m0p);
        const double got = T::mapped_curvature({std::cos(ang), std::sin(ang)}, k0, g, c, rr, rr * sphi, r_le, L);
        CHECK(sref::measure("mapped_curvature vs finite differences of the mapped circle", std::abs(got / fd - 1.0)) <= 1e-6);
    }
}

TEST_CASE("a rotating row on the default design's Mid layer (report)", "[.rotor][turbolenza][bladegen][sheet]")
{
    D::BladeDesign design = D::BladeDesign::defaults();
    D::Evaluator ev;
    const T::StreamSurfaceBlade b = T::stream_surface_blade(design, 2, design.layers[2], ev);
    const T::Cascade c = T::make_cascade(b);
    T::CaseSettings s = layer_settings();
    s.grid_level = std::string("medium");
    s.density_tolerance = 1e-5;
    s.kutta_tolerance = 1e-5;
    s.density_relaxation = 0.3;
    s.kutta_relaxation = 0.3;
    s.max_density_iterations = 6000;
    s.max_kutta_iterations = 6000;
    s.row = T::RowMotion::Rotating;
    s.inlet_absolute_mach = 0.3;
    s.inlet_absolute_flow_angle_deg = 0.0;
    // rpm for an incidence of 2 deg at the nominal far inlet: tan(beta1) = Omega r_in / C_m
    const double r_in = b.chord * b.sheet->radius(-s.upstream_extension_chords);
    const double T1 = s.inlet_total_temperature / (1.0 + 0.2 * 0.09);
    const double C = 0.3 * std::sqrt(1.4 * 287.05 * T1);
    const double beta = c.kappa_le + T::detail::radians(2.0);
    s.rpm = std::round(C * std::tan(beta) / r_in * 60.0 / (2.0 * kPi));
    const T::Case cs = T::build_case(s, b);
    std::printf("%s\n%s\n", T::format_case_summary(cs).c_str(), T::format_grid_summary(cs).c_str());
    const T::RunResult r = T::run_case(cs);
    std::printf("%s\n", T::format_run_report(r).c_str());
    CHECK(r.converged);
    REQUIRE(r.summary);
    REQUIRE(r.summary->rotor);
    const T::RotorSummary& o = *r.summary->rotor;
    // rothalpy: cp (T0_2 - T0_1) = U2 C_theta2 - U1 C_theta1 at every node, so in the mass averages too
    CHECK(std::abs(o.total_enthalpy_rise - o.euler_work) <= 1e-6 * std::abs(o.euler_work));
    CHECK(o.total_pressure_ratio > 1.0);  // a compressor rotor
}
