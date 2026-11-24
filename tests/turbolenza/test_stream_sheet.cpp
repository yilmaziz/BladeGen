// The axisymmetric stream sheet and the rotating row (not in v1): StreamSheet, the side lines on
// a sheet, the metrics and coefficients on a cone (hand-derived), the rotation source of Eq. 5-21
// (sign and magnitude on a manufactured free vortex), the velocity triangle, the Kutta correction
// with r, H' per station, mass conservation on a conical sheet, threads, and the cylinder that is
// the planar cascade to the last bit.
#include "CascadeFixture.hpp"
#include "SolverRef.hpp"

#include "turbolenza/Case.hpp"
#include "turbolenza/Detail.hpp"
#include "turbolenza/Discretization.hpp"
#include "turbolenza/GridDetail.hpp"
#include "turbolenza/Kutta.hpp"
#include "turbolenza/StreamSheet.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>

using namespace BladeGen::Turbolenza;

namespace {

constexpr double kPi = std::numbers::pi;
std::size_t sz(int k) { return static_cast<std::size_t>(k); }

/// A cone r = r0 + m sin(phi) and b = b0 (1 + k m), tabulated on [lo, hi].
std::shared_ptr<const StreamSheet> cone(double r0, double sphi, double r_ref, double b0 = 0.3, double k = 0.0,
                                        double lo = -4.0, double hi = 5.0, int n = 901)
{
    std::vector<double> m(sz(n)), r(sz(n)), b(sz(n));
    for (int q = 0; q < n; ++q) {
        m[sz(q)] = lo + (hi - lo) * q / (n - 1);
        r[sz(q)] = r0 + m[sz(q)] * sphi;
        b[sz(q)] = b0 * (1.0 + k * m[sz(q)]);
    }
    return std::make_shared<const StreamSheet>(m, r, b, r_ref);
}

std::shared_ptr<const StreamSheet> cylinder(double r0, double b0, double lo, double hi)
{
    return std::make_shared<const StreamSheet>(std::vector<double>{lo, 0.5 * (lo + hi), hi}, std::vector<double>(3, r0),
                                               std::vector<double>(3, b0), r0);
}

double simpson(const std::function<double(double)>& f, double a, double b, int n)
{
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int k = 1; k < n; ++k) s += (k % 2 ? 4.0 : 2.0) * f(a + h * k);
    return s * h / 3.0;
}

bool same_bits(const Field& a, const Field& b) { return a.same_shape(b) && a.v == b.v; }

/// The manufactured channel of SolverRef.hpp on a sheet.
std::shared_ptr<const Grid> channel_on(std::shared_ptr<const StreamSheet> sheet, double angle, double pitch = 0.8,
                                       int n_eta = 9)
{
    Grid g = sref::channel(31, n_eta, 8, 20, angle, pitch);
    g.boundaries.sheet = std::move(sheet);
    return std::make_shared<const Grid>(std::move(g));
}

}  // namespace

TEST_CASE("stream sheet: a cone and a cylinder exactly; refusals", "[turbolenza][sheet]")
{
    const auto c = cone(2.0, 0.25, 2.0, 0.3, 0.1);
    for (double m : {-3.7, -1.0, 0.0, 0.123, 2.5, 4.9}) {
        CHECK(std::abs(c->radius(m) - (2.0 + 0.25 * m)) <= 1e-14 * (2.0 + 0.25 * m));
        CHECK(std::abs(c->sin_phi(m) - 0.25) <= 1e-12);  // the cubic spline of linear data is the line
        CHECK(std::abs(c->thickness(m) - 0.3 * (1.0 + 0.1 * m)) <= 1e-14);
        CHECK(std::abs(c->scale(m) - (2.0 + 0.25 * m) / 2.0) <= 1e-14);
    }
    CHECK_FALSE(c->cylinder());
    const auto y = cylinder(1.7, 0.4, -2.0, 3.0);
    CHECK(y->cylinder());
    CHECK(y->radius(0.3) == 1.7);
    CHECK(y->sin_phi(0.3) == 0.0);
    CHECK(y->scale(-1.9) == 1.0);
    CHECK(y->thickness(2.0) == 0.4);
    CHECK_THROWS_AS(c->radius(5.1), std::out_of_range);
    CHECK_THROWS_AS(c->sin_phi(-4.2), std::out_of_range);
    CHECK_THROWS_AS(StreamSheet({0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(StreamSheet({0.0, 1.0}, {1.0, -1.0}, {1.0, 1.0}, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(StreamSheet({0.0, 1.0}, {1.0, 1.0}, {0.0, 1.0}, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(StreamSheet({0.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(StreamSheet({0.0}, {1.0}, {1.0}, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(StreamSheet({0.0, 1.0}, {1.0, std::nan("")}, {1.0, 1.0}, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(StreamSheet({0.0, 1.0}, {1.0}, {1.0, 1.0}, 1.0), std::invalid_argument);
}

TEST_CASE("side lines on a sheet: r d theta / dm = tan(beta) integrated (Eq. 5-3, p. 86)", "[turbolenza][sheet][grid]")
{
    const auto c = cone(2.0, 0.25, 2.0);
    // constant beta: y - y_j = tan(beta) r_ref / sin(phi) ln(r(m) / r(m_j)), closed form
    const double beta = 0.6;
    for (double m : {-3.0, -0.5, 1.7, 4.5}) {
        const double exact = std::tan(beta) * 2.0 / 0.25 * std::log(c->radius(m) / c->radius(0.2));
        const double got = c->side_line_rise(0.2, m, beta, beta, 4.0);
        CHECK(sref::measure("sheet side line const beta vs log closed form", std::abs(got - exact) / std::abs(exact)) <=
              1e-13);  // measured 2026-11-24: <= 2e-15
    }
    // beta linear in m from the junction (0, 0.7) to the far station (-3, 0.2): against a fine Simpson rule
    {
        const double mj = 0.0, mf = -3.0, bj = 0.7, bf = 0.2;
        auto f = [&](double s) { return std::tan(bj + (s - mj) / (mf - mj) * (bf - bj)) / ((2.0 + 0.25 * s) / 2.0); };
        for (double m : {-0.4, -1.5, -3.0}) {
            const double fine = simpson(f, mj, m, 200000);
            const double got = c->side_line_rise(mj, m, bj, bf, mf);
            CHECK(sref::measure("sheet side line linear beta vs Simpson 2e5", std::abs(got - fine) / std::abs(fine)) <=
                  1e-12);  // measured 2026-11-24: <= 3e-15
        }
    }
    // sigma constant != 1 (a cylinder with r_ref != r): v1's closed form over sigma
    {
        const StreamSheet s({-4.0, 0.0, 4.0}, {3.0, 3.0, 3.0}, {1.0, 1.0, 1.0}, 2.0);
        const double planar = (-2.0 - 0.0) * detail::mean_tan(0.5, 0.9);
        const double got = s.side_line_rise(0.0, -2.0, 0.9, 0.1, -4.0);  // beta at m = -2: 0.9 + 0.5 (0.1 - 0.9) = 0.5
        CHECK(std::abs(got - planar / 1.5) <= 1e-14 * std::abs(planar));
    }
    // a cylinder whose r_ref is its radius takes v1's closed form, bit for bit
    {
        const auto y = cylinder(1.3, 0.2, -5.0, 5.0);
        const Point j{0.1, 0.25};
        for (double fr : {0.0, 0.3, 0.5, 0.8, 1.0}) {
            const double m = 0.1 + fr * (-2.9);
            const auto a = detail::side_line(m, fr, j, 0.6, 0.2);
            const auto b = detail::side_line(m, fr, j, 0.6, 0.2, y.get(), -2.8);
            CHECK(a[0] == b[0]);
            CHECK(a[1] == b[1]);
        }
    }
}

TEST_CASE("metrics and coefficients on a cone: hand-derived (Eqs. 5-3, 5-9, 5-15, 5-22 .. 5-26)", "[turbolenza][sheet][coefficients]")
{
    const double r0 = 2.0, sphi = 0.25, t = std::tan(0.3), pitch = 0.8, rho0 = 1.2, mdot = 3.0, c_ref = 0.1;
    // constant b: every field by hand, D and E exact (B, C linear in m, A constant in eta)
    {
        const auto g = channel_on(cone(r0, sphi, r0), 0.3, pitch);
        const GridMetrics mt = build_metrics(g);
        REQUIRE(mt.on_sheet);
        for (int i = 0; i < g->n_m(); ++i) {
            const double m = g->m[sz(i)], sigma = (r0 + sphi * m) / r0;
            CHECK(std::abs(mt.S[sz(i)] - sigma * pitch) <= 1e-14);
            CHECK(std::abs(mt.tan_beta0[sz(i)] - sigma * t) <= 1e-12);
            CHECK(std::abs(mt.tan_beta1[sz(i)] - sigma * t) <= 1e-12);
            CHECK(std::abs(mt.sin_phi[sz(i)] - sphi) <= 1e-12);
            CHECK(mt.thickness[sz(i)] == 1.0);
            CHECK(std::abs(mt.dS_dm[sz(i)] - sphi * pitch / r0) <= 1e-12);  // sin(phi) Theta
        }
        const Field rho(g->n_m(), g->n_eta(), rho0);
        const CoefficientFields co = build_coefficients(mt, rho, mdot, c_ref);
        double worst = 0.0;
        for (int i = 0; i < g->n_m(); ++i)
            for (int j = 0; j < g->n_eta(); ++j) {
                const double sigma = (r0 + sphi * g->m[sz(i)]) / r0, tb = sigma * t;
                const double S = c_ref * sigma * pitch, cos2 = 1.0 / (1.0 + tb * tb);
                const double A = mdot / (S * rho0 * cos2), B = mdot * tb / rho0, C = mdot * S / rho0;
                const double D = -mdot / rho0 * t * sphi / (r0 * c_ref);  // -dB/dm, m physical
                const double E = mdot / rho0 * sphi / r0 * pitch;         // dC/dm
                worst = std::max({worst, std::abs(co.A(i, j) / A - 1.0), std::abs(co.B(i, j) / B - 1.0),
                                  std::abs(co.C(i, j) / C - 1.0), std::abs(co.D(i, j) / D - 1.0),
                                  std::abs(co.E(i, j) / E - 1.0)});
            }
        CHECK(sref::measure("cone coefficients A..E vs hand (relative)", worst) <= 1e-11);  // measured: 1.3e-13
        CHECK_FALSE(co.has_source());
        const CoefficientFields rot = build_coefficients(mt, rho, mdot, c_ref, 50.0);
        REQUIRE(rot.has_source());
        for (int i = 0; i < g->n_m(); ++i) {
            const double S = c_ref * (r0 + sphi * g->m[sz(i)]) / r0 * pitch;
            CHECK(std::abs(rot.Q(i, 3) - 2.0 * S * 50.0 * sphi) <= 1e-13 * std::abs(2.0 * S * 50.0 * sphi));
        }
        CHECK(same_bits(rot.A, co.A));  // the source changes nothing else
    }
    // b varying along m (b / b(m_in) in rho b), A, B, C by hand
    {
        const double k = 0.1;
        const auto g = channel_on(cone(r0, sphi, r0, 0.3, k), 0.3, pitch);
        const GridMetrics mt = build_metrics(g);
        const Field rho(g->n_m(), g->n_eta(), rho0);
        const CoefficientFields co = build_coefficients(mt, rho, mdot, c_ref);
        double worst = 0.0;
        for (int i = 0; i < g->n_m(); ++i) {
            const double b = (1.0 + k * g->m[sz(i)]) / (1.0 + k * g->m.front());
            CHECK(std::abs(mt.thickness[sz(i)] - b) <= 1e-14);
            const double sigma = (r0 + sphi * g->m[sz(i)]) / r0, tb = sigma * t;
            const double S = c_ref * sigma * pitch, cos2 = 1.0 / (1.0 + tb * tb);
            worst = std::max({worst, std::abs(co.A(i, 4) / (mdot / (S * b * rho0 * cos2)) - 1.0),
                              std::abs(co.B(i, 4) / (mdot * tb / (b * rho0)) - 1.0),
                              std::abs(co.C(i, 4) / (mdot * S / (b * rho0)) - 1.0)});
        }
        CHECK(sref::measure("cone, linear b: A, B, C vs hand (relative)", worst) <= 1e-12);
        CHECK(mt.thickness.front() == 1.0);
    }
    // a planar grid takes no rotation
    {
        const auto g = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, 0.3, pitch));
        const GridMetrics mt = build_metrics(g);
        CHECK_FALSE(mt.on_sheet);
        CHECK_THROWS_AS(build_coefficients(mt, Field(31, 9, 1.0), mdot, c_ref, 10.0), CoefficientError);
    }
}

TEST_CASE("the rotation source 2 S omega sin(phi) of Eq. 5-21: sign and magnitude on a free vortex", "[turbolenza][sheet][rotation]")
{
    // A cone (r = 0.1 (2 + 0.3 m*) m, sin(phi) = 0.3), constant rho and b, grid lines along m (tan beta = 0).
    // The absolutely irrotational flow r C_theta = K has W_theta = K / r - omega r; with psi = eta + g(m),
    // Eq. 5-19 gives g'(m) = -(b rho / mdot)(K / r - omega r), and Eq. 5-21 reduces to
    // (C g')' = d/dm[-S (K / r - omega r)] = Theta d/dm(omega r^2) = 2 S omega sin(phi): the source exactly.
    const double c_ref = 0.1, r0 = 2.0, sphi = 0.3, rho0 = 1.2, mdot = 3.0, omega = -400.0, K = 6.0;
    const auto g = channel_on(cone(r0, sphi, r0, 1.0), 0.0, 0.8);
    const GridMetrics mt = build_metrics(g);
    const Field rho(g->n_m(), g->n_eta(), rho0);
    const CoefficientFields with = build_coefficients(mt, rho, mdot, c_ref, omega);
    const CoefficientFields none = build_coefficients(mt, rho, mdot, c_ref, 0.0);
    Field psi(g->n_m(), g->n_eta());
    for (int i = 0; i < g->n_m(); ++i) {
        const double r = c_ref * (r0 + sphi * g->m[sz(i)]);  // dr/dm = sin(phi) in physical m
        const double G = -(rho0 / mdot) * (K * std::log(r) / sphi - omega * r * r / (2.0 * sphi));
        for (int j = 0; j < g->n_eta(); ++j) psi(i, j) = g->eta[sz(j)] + G;
    }
    const InteriorStencil sw = build_interior_stencil(with, mt), sn = build_interior_stencil(none, mt);
    double worst = 0.0, ratio_worst = 0.0;
    for (int i = 1; i < g->n_m() - 1; ++i)
        for (int j = 1; j < g->n_eta() - 1; ++j) {
            auto lhs = [&](const InteriorStencil& s) {
                return s.center(i, j) * psi(i, j) + s.west(i, j) * psi(i - 1, j) + s.east(i, j) * psi(i + 1, j) +
                       s.south(i, j) * psi(i, j - 1) + s.north(i, j) * psi(i, j + 1) +
                       s.mixed(i, j) * (psi(i + 1, j + 1) - psi(i + 1, j - 1) - psi(i - 1, j + 1) + psi(i - 1, j - 1));
            };
            const double q_over_h = with.Q(i, j) / sw.H(i, j);
            worst = std::max(worst, std::abs(lhs(sw) - sw.rhs(i, j)) / std::abs(q_over_h));
            // without the source the exact psi leaves -Q / H: the right sign and size of the source
            ratio_worst = std::max(ratio_worst, std::abs((lhs(sn) - sn.rhs(i, j)) / -q_over_h - 1.0));
            CHECK(sw.rhs(i, j) == -q_over_h);
        }
    CHECK(sref::measure("rotation source: residual with Q / |Q / H|", worst) <= 2e-3);  // truncation, O(dm^2)
    CHECK(sref::measure("rotation source: residual without Q / (-Q / H) - 1", ratio_worst) <= 2e-3);
    // the periodic PDE rows carry it too
    const BoundaryConditions bcs = build_boundary_conditions(g);
    const auto rows_w = build_periodic_pde_rows(with, mt, bcs), rows_n = build_periodic_pde_rows(none, mt, bcs);
    REQUIRE(rows_w.size() == rows_n.size());
    for (std::size_t q = 0; q < rows_w.size(); ++q) {
        const int i = rows_w[q].target.i, J = g->n_eta() - 1;
        CHECK(std::abs(rows_w[q].rhs - (rows_n[q].rhs - with.Q(i, J) / rows_w[q].H)) <= 1e-12 * std::abs(rows_n[q].rhs) + 1e-15);
    }
}

TEST_CASE("velocity triangle: closed form, both signs (TURBOlenza frame, omega = -2 pi rpm / 60)", "[turbolenza][rotation]")
{
    const PerfectGas gas(1.4, 287.05);
    const double P0 = 101325.0, T0 = 288.15, r = 0.2;
    const double T1 = T0 / (1.0 + 0.2 * 0.25), a1 = std::sqrt(1.4 * 287.05 * T1), C = 0.5 * a1;
    const double Om = 2.0 * kPi * 6000.0 / 60.0;
    const VelocityTriangle p = velocity_triangle(6000.0, 0.0, 0.5, P0, T0, gas, r);
    CHECK(std::abs(p.omega + Om) <= 1e-12 * Om);
    CHECK(std::abs(p.U + Om * r) <= 1e-12 * Om * r);           // the blades move toward -theta
    CHECK(std::abs(p.W_theta - Om * r) <= 1e-12 * Om * r);     // W_theta = C_theta - U > 0
    CHECK(std::abs(p.beta - std::atan(Om * r / C)) <= 1e-14);  // beta1 > 0: a compressor rotor
    CHECK(p.beta > 0.0);
    CHECK(std::abs(p.mach_relative - std::hypot(C, Om * r) / a1) <= 1e-14);
    CHECK(std::abs(p.rothalpy - gas.cp() * T0) <= 1e-12 * gas.cp() * T0);  // no swirl: I = cp T0
    const VelocityTriangle n = velocity_triangle(-6000.0, 0.0, 0.5, P0, T0, gas, r);
    CHECK(n.beta == -p.beta);
    CHECK(n.mach_relative == p.mach_relative);
    // swirl: alpha1 = 20 deg toward +theta, against the rotation of rpm > 0
    const VelocityTriangle s = velocity_triangle(6000.0, 20.0, 0.5, P0, T0, gas, r);
    const double ct = C * std::sin(20.0 * kPi / 180.0), cm = C * std::cos(20.0 * kPi / 180.0);
    CHECK(std::abs(s.C_theta - ct) <= 1e-12 * C);
    CHECK(std::abs(s.W_theta - (ct + Om * r)) <= 1e-12 * C);
    CHECK(std::abs(s.beta - std::atan((ct + Om * r) / cm)) <= 1e-14);
    CHECK(std::abs(s.rothalpy - (gas.cp() * T0 + Om * r * ct)) <= 1e-12 * gas.cp() * T0);  // I = cp T0 - U C_theta
    // the relative inlet state the solver closes from it (inlet_state, rotating): the same static state
    FlowConditions f;
    f.inlet_mach = s.mach_relative;
    f.inlet_angle_deg = detail::degrees(s.beta);
    f.inlet_absolute_total_pressure = P0;
    f.inlet_absolute_total_temperature = T0;
    f.angular_speed = s.omega;
    f.stream_surface_radius = r;
    f.inlet_absolute_tangential_velocity = s.C_theta;
    const InletState in = inlet_state(f, gas, 0.05);
    CHECK(sref::measure("triangle vs inlet_state T1 (relative)", std::abs(in.static_temperature / T1 - 1.0)) <= 1e-12);
    CHECK(std::abs(in.meridional_velocity / cm - 1.0) <= 1e-10);
    CHECK(in.rotating);
    CHECK(in.relative_total_pressure > P0);  // W > C here
    CHECK(std::abs(in.mass_flow_per_depth - in.static_density * cm * 0.05) <= 1e-9 * in.mass_flow_per_depth);
    CHECK_THROWS_AS(velocity_triangle(6000.0, 0.0, 0.0, P0, T0, gas, r), InletStateError);
    CHECK_THROWS_AS(velocity_triangle(6000.0, 90.0, 0.5, P0, T0, gas, r), InletStateError);
    CHECK_THROWS_AS(velocity_triangle(6000.0, 0.0, 0.5, P0, T0, gas, 0.0), InletStateError);
    CHECK_THROWS_AS(velocity_triangle(std::nan(""), 0.0, 0.5, P0, T0, gas, r), InletStateError);
}

TEST_CASE("gas state: H' = I + (r_i omega)^2 / 2 per station (Eq. 5-51)", "[turbolenza][rotation][eos]")
{
    const PerfectGas gas;
    const TotalState total(101325.0, 288.15);
    BladeRowReference ref = blade_row_reference(total, gas, -500.0, 0.2, 30.0);
    ref.station_blade_speed = {-100.0, -110.0, -120.0};
    Field W(3, 2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) W(i, j) = 150.0 + 10.0 * i + 5.0 * j;
    const StaticState s = static_state_from_velocity(W, gas, ref);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j) {
            const double U = ref.station_blade_speed[sz(i)];
            CHECK(s.enthalpy(i, j) == (ref.rothalpy + 0.5 * (U * U)) - 0.5 * (W(i, j) * W(i, j)));
        }
    ref.station_blade_speed = {1.0, 2.0};
    CHECK_THROWS_AS(static_state_from_velocity(W, gas, ref), ThermodynamicStateError);
}

TEST_CASE("Kutta correction with r: Eqs. 13-57 (r_te) and 13-59 (r_out)", "[turbolenza][sheet][kutta]")
{
    const auto g = channel_on(cone(2.0, 0.25, 2.0), 0.3, 0.8);
    const GridMetrics mt = build_metrics(g);
    FlowField flow{Field(31, 9), Field(31, 9), Field(31, 9, 100.0), Field(31, 9, 10.0), Field(31, 9), Field(31, 9)};
    const int l = g->last_blade_i;
    flow.W_m(l, 0) = 110.0;  // the TE wall speeds differ
    const KuttaResidual res = kutta_residual(flow, *g);
    const ExitAngleCorrection c = exit_angle_correction(res, flow, mt);
    const double r_te = mt.radius[sz(l)], r_out = mt.radius.back();
    const double t0 = mt.tan_beta0[sz(l)], t1 = mt.tan_beta1[sz(l)];
    const double hand = 2.0 * res.wall_difference * g->dm() * r_te /
                        (mt.S[sz(l)] * (1.0 / std::sqrt(1.0 + t0 * t0) + 1.0 / std::sqrt(1.0 + t1 * t1))) /
                        (r_out * 100.0);
    CHECK(std::abs(c.delta_tan - hand) <= 1e-13 * std::abs(hand));
    CHECK(std::abs(c.radius_ratio - r_te / r_out) <= 1e-15);
    CHECK(c.radius_ratio < 1.0);  // the cone opens downstream
}

namespace {

SolverConfig channel_config()
{
    SolverConfig cfg;
    cfg.flow.inlet_angle_deg = 17.0;
    cfg.iteration.density_tolerance = 1e-11;
    cfg.iteration.max_density_iterations = 300;
    return cfg;
}

}  // namespace

namespace {

/// Eq. 5-17 at every node, rho b S W_m = mdot' psi_eta (relative), and the station diagnostic.
struct MassCheck
{
    double pointwise = 0.0, station = 0.0;
};
MassCheck mass_check(const Solution& s, double mdot)
{
    MassCheck c;
    const Grid& g = *s.grid;
    for (int i = 0; i < g.n_m(); ++i)
        for (int j = 0; j < g.n_eta(); ++j) {
            const double flux = s.density(i, j) * s.metrics.thickness_at(i) * (s.c_ref * s.metrics.S[sz(i)]) * s.flow.W_m(i, j);
            c.pointwise = std::max(c.pointwise, std::abs(flux / (mdot * s.flow.psi_eta(i, j)) - 1.0));
        }
    for (double q : station_mass_flow(s.flow, s.density, s.coefficients, s.metrics))
        c.station = std::max(c.station, std::abs(q / mdot - 1.0));
    return c;
}

}  // namespace

TEST_CASE("mass conservation on a conical stream sheet (Eqs. 5-16, 5-17)", "[turbolenza][sheet][solver]")
{
    // Eq. 5-17 holds node by node by construction; the station integral (trapezoidal rule of the recovered
    // W_m, v1's diagnostic) is exact for psi linear in eta (the planar channel) and otherwise second order
    // in deta: on the cone psi bends across the passage (D = -dB/dm != 0).  v1's own cascades report
    // 1.7e-3 .. 3.4e-3 (its records): that is "the solver's level"
    const double c_ref = 0.1, mdot = 6.0;
    const auto planar = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, 0.3, 0.8));
    const Solution p = solve_at_exit_angle(planar, channel_config(), c_ref, mdot);
    CHECK(sref::measure("planar channel: station mass-flow error", mass_check(p, mdot).station) <= 1e-14);
    for (double k : {0.0, 0.15}) {
        std::vector<double> errors;
        for (int n_eta : {9, 17, 33}) {
            const auto g = channel_on(cone(2.0, 0.25, 2.0, 0.3, k), 0.3, 0.8, n_eta);
            const Solution s = solve_at_exit_angle(g, channel_config(), c_ref, mdot);
            REQUIRE(s.converged);
            const MassCheck m = mass_check(s, mdot);
            CHECK(sref::measure("cone channel: Eq. 5-17 node by node", m.pointwise) <= 1e-13);
            errors.push_back(sref::measure(std::string("cone channel, b ") + (k == 0.0 ? "const" : "linear") +
                                               ", n_eta " + std::to_string(n_eta) + ": station error",
                                           m.station));
            if (n_eta == 9) CHECK(s.flow.W_m(30, 4) < s.flow.W_m(0, 4));  // S b grows: the flow slows down
        }
        CHECK(errors[0] < 1e-3);
        // second order (dm held, so the ratio tends to 4 from below): halving deta divides the error by 2.6 .. 3.2
        CHECK(sref::measure("cone channel: station error ratio n_eta 9 / 17", errors[0] / errors[1]) > 2.2);  // measured 2026-11-24: 2.63, 2.65
        CHECK(sref::measure("cone channel: station error ratio n_eta 17 / 33", errors[1] / errors[2]) > 2.8);  // measured: 3.16, 3.04 (toward 4)
    }
}

namespace {

/// A rotating row on the cone channel (mass flow given; the rothalpy at the far-inlet radius).
SolverConfig rotating_config(const Grid& g, double c_ref, double omega)
{
    SolverConfig cfg = channel_config();
    cfg.flow.angular_speed = omega;
    cfg.flow.stream_surface_radius = c_ref * g.boundaries.sheet->radius(g.m.front());
    cfg.flow.inlet_absolute_tangential_velocity = 20.0;
    return cfg;
}

}  // namespace

TEST_CASE("a rotating row on a cone: converges, mass conserved, 1 and 4 threads bit-identical (SparseLU)", "[turbolenza][sheet][rotation][threads]")
{
    const double c_ref = 0.1, mdot = 6.0, omega = -1500.0;
    const auto g = channel_on(cone(2.0, 0.25, 2.0, 0.3, 0.1), 0.3, 0.8);
    SolveOptions one;
    const Solution a = solve_at_exit_angle(g, rotating_config(*g, c_ref, omega), c_ref, mdot, one);
    REQUIRE(a.converged);
    CHECK(sref::measure("rotating cone channel: Eq. 5-17 node by node", mass_check(a, mdot).pointwise) <= 1e-13);
    {  // the station diagnostic: second order in deta, as without rotation
        const auto g17 = channel_on(cone(2.0, 0.25, 2.0, 0.3, 0.1), 0.3, 0.8, 17);
        const Solution a17 = solve_at_exit_angle(g17, rotating_config(*g17, c_ref, omega), c_ref, mdot, one);
        REQUIRE(a17.converged);
        const double e9 = sref::measure("rotating cone channel, n_eta 9: station error", mass_check(a, mdot).station);
        const double e17 = sref::measure("rotating cone channel, n_eta 17: station error", mass_check(a17, mdot).station);
        CHECK(sref::measure("rotating cone channel: station error ratio 9 / 17", e9 / e17) > 2.2);  // measured 2026-11-24: 2.50
    }
    CHECK(a.coefficients.has_source());
    // H' varies with r: on the cone r grows, so h = I + U^2/2 - W^2/2 at equal W grows downstream
    const Solution st = solve_at_exit_angle(g, channel_config(), c_ref, mdot, one);
    CHECK_FALSE(same_bits(a.psi, st.psi));  // the source and H'(r) act
    if (max_threads() < 4) SKIP("fewer than 4 hardware threads");
    SolveOptions four;
    four.threads = 4;
    const Solution b = solve_at_exit_angle(g, rotating_config(*g, c_ref, omega), c_ref, mdot, four);
    REQUIRE(a.passes.size() == b.passes.size());
    CHECK(same_bits(a.psi, b.psi));
    CHECK(same_bits(a.density, b.density));
    CHECK(same_bits(a.mach, b.mach));
    CHECK(same_bits(a.static_state.pressure, b.static_state.pressure));
    CHECK(same_bits(a.flow.W_theta, b.flow.W_theta));
    // refusals: a rotating row needs its sheet and the far-inlet radius
    const auto planar = std::make_shared<const Grid>(sref::channel(31, 9, 8, 20, 0.3, 0.8));
    CHECK_THROWS_AS(solve_at_exit_angle(planar, rotating_config(*g, c_ref, omega), c_ref, mdot), SolverInputError);
    SolverConfig wrong = rotating_config(*g, c_ref, omega);
    wrong.flow.stream_surface_radius = *wrong.flow.stream_surface_radius * 1.01;
    CHECK_THROWS_AS(solve_at_exit_angle(g, wrong, c_ref, mdot), SolverInputError);
}

namespace {

CaseSettings synthetic_settings()
{
    CaseSettings s;  // test_viscous.cpp's synthetic case (converges in a few hundred passes)
    s.case_name = "synthetic sheet";
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

}  // namespace

namespace {

/// The synthetic blade on a cylinder of the same pitch: r_ref = pitch N / (2 pi).
StreamSurfaceBlade on_cylinder(const CaseSettings& s, int N)
{
    StreamSurfaceBlade b;
    b.section = fixture::synthetic_section();
    b.stagger_deg = s.stagger_deg;
    b.solidity = s.solidity;
    b.chord = s.chord;
    b.n_blades = N;
    const double r = (1.0 / s.solidity) * N / (2.0 * kPi);
    b.sheet = cylinder(r, 0.25, -6.0, 7.0);
    b.source = "test cylinder";
    return b;
}

}  // namespace

TEST_CASE("cylinder: the planar cascade of the same pitch and chord, to the last bit", "[turbolenza][sheet][case]")
{
    CaseSettings s = synthetic_settings();
    const Case planar = build_case(s, fixture::synthetic_section());
    s.stream_surface = StreamSurfaceKind::Layer;
    const StreamSurfaceBlade blade = on_cylinder(s, 30);
    const Case cyl = build_case(s, blade);
    REQUIRE(cyl.cascade.sheet);
    CHECK(std::abs(1.0 / cyl.cascade.solidity - 1.0 / planar.cascade.solidity) <= 1e-15);
    CHECK(same_bits(planar.grid.node_y, cyl.grid.node_y));
    CHECK(planar.mass_flow_per_depth == cyl.mass_flow_per_depth);
    const RunResult a = run_case(planar), b = run_case(cyl);
    REQUIRE(a.converged);
    REQUIRE(b.converged);
    REQUIRE(a.records.size() == b.records.size());
    CHECK(a.solution->exit_angle == b.solution->exit_angle);
    CHECK(same_bits(a.solution->psi, b.solution->psi));
    CHECK(same_bits(a.solution->mach, b.solution->mach));
    CHECK(same_bits(a.solution->pressure, b.solution->pressure));
    sref::measure("cylinder vs planar beta_out difference [rad]", std::abs(a.solution->exit_angle - b.solution->exit_angle));
    CHECK(b.solution->station_radius.size() == static_cast<std::size_t>(b.solution->grid.n_m()));
    CHECK(a.solution->station_radius.empty());
}

TEST_CASE("cylinder on the v1 45 x 21 grid: the planar solve bit for bit", "[turbolenza][sheet][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(d, "fixed_exit_explicit_45x21");
    const auto planar = sref::grid(d.at("first_grid"));
    Grid g = *planar;
    const double r = g.pitch() * 20.0 / (2.0 * kPi);
    g.boundaries.sheet = cylinder(r, 0.05, g.m.front() - 1.0, g.m.back() + 1.0);
    const double r_ref_sigma = g.boundaries.sheet->scale(0.0);
    REQUIRE(r_ref_sigma == 1.0);
    const SolverConfig cfg = sref::config(d);
    const double c_ref = ref::num(d.at("settings").at("chord"));
    const Solution a = solve_at_exit_angle(planar, cfg, c_ref);
    const Solution b = solve_at_exit_angle(std::make_shared<const Grid>(g), cfg, c_ref);
    REQUIRE(a.passes.size() == b.passes.size());
    CHECK(same_bits(a.psi, b.psi));
    CHECK(same_bits(a.density, b.density));
    CHECK(same_bits(a.coefficients.E, b.coefficients.E));
    CHECK(same_bits(a.flow.W_theta, b.flow.W_theta));
}

TEST_CASE("case: stream surface and row refusals; the rotating case's derived inlet", "[turbolenza][sheet][case]")
{
    CaseSettings s = synthetic_settings();
    // a rotating planar cascade has no radius
    CaseSettings rp = s;
    rp.row = RowMotion::Rotating;
    rp.rpm = 3000.0;
    CHECK_THROWS_AS(build_case(rp, fixture::synthetic_section()), CaseError);
    // omega with a stationary row
    CaseSettings om = s;
    om.omega = 10.0;
    CHECK_THROWS_AS(build_case(om, fixture::synthetic_section()), CaseError);
    // a layer setting without its blade, and a blade with a planar setting
    CaseSettings ls = s;
    ls.stream_surface = StreamSurfaceKind::Layer;
    CHECK_THROWS_AS(build_case(ls, fixture::synthetic_section()), CaseError);
    CHECK_THROWS_AS(build_case(s, on_cylinder(s, 30)), CaseError);
    // the boundary layer: on a cylinder yes, on a cone no (Eq. 5-106 has b(x))
    CaseSettings v = ls;
    v.viscous = true;
    CHECK_NOTHROW(build_case(v, on_cylinder(s, 30)));
    StreamSurfaceBlade coned = on_cylinder(s, 30);
    const double r = coned.sheet->r_ref();
    coned.sheet = cone(r, 0.1, r, 0.3, 0.0, -6.0, 7.0);
    CHECK_THROWS_AS(build_case(v, coned), CaseError);
    // the rotating row on the cone: beta1, M1_rel and omega derived at the grid's far inlet
    CaseSettings rot = ls;
    rot.row = RowMotion::Rotating;
    rot.rpm = 3000.0;
    rot.inlet_absolute_mach = 0.35;
    rot.inlet_absolute_flow_angle_deg = 0.0;
    const Case c = build_case(rot, coned);
    REQUIRE(c.triangle);
    const double r_in = rot.chord * coned.sheet->radius(c.grid.m.front());
    CHECK(std::abs(c.triangle->radius - r_in) <= 1e-15 * r_in);
    CHECK(c.settings.omega == c.triangle->omega);
    CHECK(c.settings.omega < 0.0);
    CHECK(std::abs(c.settings.inlet_flow_angle_deg - detail::degrees(c.triangle->beta)) <= 1e-12);
    CHECK(c.triangle->beta > 0.0);
    CHECK(std::abs(c.inlet->angle - c.triangle->beta) <= 1e-15);
    CHECK(std::abs(c.inlet->mach - c.triangle->mach_relative) <= 1e-15);
    rot.viscous = true;
    CHECK_THROWS_AS(build_case(rot, coned), CaseError);
}
