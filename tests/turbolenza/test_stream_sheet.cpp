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
            // d tan(beta)/dm = d(sigma t)/dm = t sin(phi) / r_ref: the blade walls by d sigma (Metrics.cpp), the
            // side lines by the run differences of a linear function
            for (int j = 0; j < g->n_eta(); ++j) CHECK(std::abs(mt.dtan_beta_dm(i, j) - t * sphi / r0) <= 1e-12);
            // d tan(beta)/d eta = tan(beta1) - tan(beta0) = 0 here, while dS/dm = sin(phi) Theta != 0
            CHECK(std::abs(mt.dtan_beta_deta()[sz(i)]) <= 1e-14);  // round-off of tan1 - tan0, both about 0.31 sigma
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

namespace {

/// Kelvin on the unbladed stations: Gamma_i = r_i int C_theta d eta, C_theta = W_theta + omega r_i (the
/// circulation round one pitch over Theta, constant where the absolute flow is irrotational and no blade is
/// crossed: Aungier Eq. 5-11).  The largest relative change of Gamma over stations [from, to), against its
/// value at `from`, and that value.
std::pair<double, double> circulation_drift(const Grid& g, const Field& W_theta, const std::vector<double>& r, double omega,
                                            int from, int to)
{
    std::vector<double> ct(sz(g.n_eta()));
    auto gamma = [&](int i) {
        for (int j = 0; j < g.n_eta(); ++j) ct[sz(j)] = W_theta(i, j) + omega * r[sz(i)];
        return r[sz(i)] * detail::trapezoid(ct.data(), g.eta);
    };
    const double ref = gamma(from);
    double drift = 0.0;
    for (int i = from; i < to; ++i) drift = std::max(drift, std::abs(gamma(i) - ref));
    return {drift / std::abs(ref), ref};
}

}  // namespace

TEST_CASE("a rotating row on a cone: converges, mass conserved, 1 and 4 threads bit-identical (SparseLU)", "[turbolenza][sheet][rotation][threads]")
{
    // omega = -200 rad/s: U up to 50 m/s against W_m about 60 m/s, no recirculation (the reviewer measured the
    // recirculating flow of omega = -1500: converged in eta, not in m)
    const double c_ref = 0.1, mdot = 6.0, omega = -200.0;
    const auto g = channel_on(cone(2.0, 0.25, 2.0, 0.3, 0.1), 0.3, 0.8);
    SolveOptions one;
    const Solution a = solve_at_exit_angle(g, rotating_config(*g, c_ref, omega), c_ref, mdot, one);
    REQUIRE(a.converged);
    CHECK(sref::measure("rotating cone channel: Eq. 5-17 node by node", mass_check(a, mdot).pointwise) <= 1e-13);
    {  // the station diagnostic: second order in deta at fixed dm (the ratio tends to 4 from below), as without rotation
        const auto g17 = channel_on(cone(2.0, 0.25, 2.0, 0.3, 0.1), 0.3, 0.8, 17);
        const Solution a17 = solve_at_exit_angle(g17, rotating_config(*g17, c_ref, omega), c_ref, mdot, one);
        REQUIRE(a17.converged);
        const double e9 = sref::measure("rotating cone channel, n_eta 9: station error", mass_check(a, mdot).station);
        const double e17 = sref::measure("rotating cone channel, n_eta 17: station error", mass_check(a17, mdot).station);
        CHECK(sref::measure("rotating cone channel: station error ratio 9 / 17", e9 / e17) > 2.2);
    }
    {  // Kelvin: r C_theta constant on the unbladed stations (the source 2 S omega sin(phi) makes it so)
        const std::vector<double> r = [&] {
            std::vector<double> v;
            for (int i = 0; i < g->n_m(); ++i) v.push_back(c_ref * g->boundaries.sheet->radius(g->m[sz(i)]));
            return v;
        }();
        const auto in = circulation_drift(*g, a.flow.W_theta, r, omega, 0, g->first_blade_i);
        const auto out = circulation_drift(*g, a.flow.W_theta, r, omega, g->last_blade_i + 1, g->n_m());
        CHECK(sref::measure("rotating cone channel: r C_theta drift, inlet stations", in.first) <= 2e-3);
        CHECK(sref::measure("rotating cone channel: r C_theta drift, exit stations", out.first) <= 2e-3);
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
    // the mass flow per passage: mdot' = rho1 C_m S_in, S_in = r_in Theta = sigma(m_in) pitch (Eq. 5-9; Grid.cpp
    // inlet_passage_width), by hand from the triangle's static state
    {
        const double S_in = rot.chord * (coned.sheet->radius(c.grid.m.front()) / coned.sheet->r_ref()) / coned.solidity;
        const double hand = c.triangle->static_density * c.triangle->C_m * S_in;
        CHECK(sref::measure("rotating cone: mdot' vs rho1 C_m S_in by hand", std::abs(c.mass_flow_per_depth / hand - 1.0)) <=
              1e-12);
    }
    // the same for the stationary row on the cone: isentropic inlet state by hand (M1 0.5, beta1 45 deg)
    {
        const Case st = build_case(ls, coned);
        const double T0 = ls.inlet_total_temperature, P0 = ls.inlet_total_pressure, M = *ls.inlet_mach;
        const double T1 = T0 / (1.0 + 0.2 * M * M), p1 = P0 * std::pow(T1 / T0, 3.5), rho1 = p1 / (287.05 * T1);
        const double W1 = M * std::sqrt(1.4 * 287.05 * T1);
        const double S_in = ls.chord * (coned.sheet->radius(st.grid.m.front()) / coned.sheet->r_ref()) / coned.solidity;
        const double hand = rho1 * W1 * std::cos(detail::radians(ls.inlet_flow_angle_deg)) * S_in;
        CHECK(sref::measure("stationary cone: mdot' vs rho1 W1 cos(beta1) S_in by hand",
                            std::abs(st.mass_flow_per_depth / hand - 1.0)) <= 1e-12);
        CHECK(std::abs(coned.sheet->radius(st.grid.m.front()) / coned.sheet->r_ref() - 1.0) > 0.01);  // sigma_in != 1
    }
    rot.viscous = true;
    CHECK_THROWS_AS(build_case(rot, coned), CaseError);
    // the touching edge-node experiment is planar only
    CaseSettings tch = ls;
    tch.edge_nodes = EdgeNodes::Touching;
    CHECK_THROWS_AS(build_case(tch, coned), CaseError);
    CHECK_NOTHROW(build_case(tch, on_cylinder(s, 30)));  // sigma = 1: the planar cascade
    // extensions beyond the sheet (it ends 6 chords ahead of the LE): stationary in the grid, rotating at the triangle
    CaseSettings far = ls;
    far.upstream_extension_chords = 8.0;
    CHECK_THROWS_AS(build_case(far, coned), CaseError);
    CaseSettings far_rot = rot;
    far_rot.viscous = false;
    far_rot.upstream_extension_chords = 8.0;
    try {
        (void)build_case(far_rot, coned);
        FAIL("an inlet beyond the sheet must be refused");
    } catch (const CaseError& e) {
        CHECK(std::string(e.what()).find("leaves the stream sheet") != std::string::npos);
    }
    // the boundary layer on a cylinder whose radius is not r_ref (sigma != 1): refused, the march's arc length is wrong
    {
        StreamSurfaceBlade sig = on_cylinder(s, 30);
        const double rr = sig.sheet->r_ref();
        sig.sheet = std::make_shared<const StreamSheet>(std::vector<double>{-6.0, 0.0, 7.0}, std::vector<double>(3, 1.1 * rr),
                                                        std::vector<double>(3, 0.25), rr);
        REQUIRE(sig.sheet->cylinder());
        CHECK_FALSE(sig.sheet->unit_scale());
        CHECK_NOTHROW(build_case(ls, sig));
        CHECK_THROWS_AS(build_case(v, sig), CaseError);
    }
}

TEST_CASE("make_cascade of a stream-surface blade and the metrics refuse what does not fit", "[turbolenza][sheet][case]")
{
    const CaseSettings s = synthetic_settings();
    StreamSurfaceBlade b = on_cylinder(s, 30);
    const double r = b.sheet->r_ref();
    b.sheet = cone(r, 0.1, r, 0.3, 0.0, -6.0, 7.0);
    CHECK_NOTHROW(make_cascade(b));
    StreamSurfaceBlade pitch = b;
    pitch.solidity = 1.01 * b.solidity;  // 1 / solidity is not r_ref 2 pi / N
    CHECK_THROWS_AS(make_cascade(pitch), GeometryError);
    StreamSurfaceBlade cover = b;
    cover.sheet = cone(r, 0.1, r, 0.3, 0.0, 0.5, 7.0);  // starts behind the LE
    CHECK_THROWS_AS(make_cascade(cover), GeometryError);
    StreamSurfaceBlade none = b;
    none.sheet.reset();
    CHECK_THROWS_AS(make_cascade(none), std::invalid_argument);
    StreamSurfaceBlade zero = b;
    zero.n_blades = 0;
    CHECK_THROWS_AS(make_cascade(zero), std::invalid_argument);
    StreamSurfaceBlade chord = b;
    chord.chord = 0.0;
    CHECK_THROWS_AS(make_cascade(chord), std::invalid_argument);
    // a grid longer than its sheet (m -1 .. 2 against -0.5 .. 5)
    CHECK_THROWS_AS(build_metrics(channel_on(cone(2.0, 0.25, 2.0, 0.3, 0.0, -0.5, 5.0), 0.3)), MetricsError);
}

TEST_CASE("a rotating cascade on a cone: Kelvin, the absolute frame, the rotor summary (hand values)", "[turbolenza][sheet][rotation][case]")
{
    // the synthetic blade on a cone with sin(phi) = 0.3 (r grows 6 % per chord), N = 30, very coarse grid, rpm for
    // beta1 = kappa_LE at the nominal far inlet, alpha1 = 0, M1 0.3 absolute
    CaseSettings s = synthetic_settings();
    s.stream_surface = StreamSurfaceKind::Layer;
    s.row = RowMotion::Rotating;
    s.inlet_absolute_mach = 0.3;
    s.inlet_absolute_flow_angle_deg = 0.0;
    StreamSurfaceBlade b = on_cylinder(s, 30);
    const double rr = b.sheet->r_ref(), sphi = 0.3;
    b.sheet = cone(rr, sphi, rr, 0.25, 0.0, -6.0, 7.0);
    const Cascade cc = make_cascade(b);
    const double T0 = s.inlet_total_temperature, P0 = s.inlet_total_pressure;
    const double T1 = T0 / (1.0 + 0.2 * 0.09), C1 = 0.3 * std::sqrt(1.4 * 287.05 * T1);
    const double r_in_nominal = s.chord * b.sheet->radius(-s.upstream_extension_chords);
    s.rpm = std::round(C1 * std::tan(cc.kappa_le) / r_in_nominal * 60.0 / (2.0 * kPi));
    const Case c = build_case(s, b);
    REQUIRE(c.triangle);
    const RunResult r = run_case(c);
    REQUIRE(r.converged);
    sref::measure("rotating cone cascade: passes", r.iterations);
    sref::measure("rotating cone cascade: beta_out [deg]", detail::degrees(r.solution->exit_angle));
    const SolutionData& d = *r.solution;
    REQUIRE(d.absolute);
    REQUIRE(r.summary);
    REQUIRE(r.summary->rotor);
    const AbsoluteFrameData& a = *d.absolute;
    const PerfectGas gas(1.4, 287.05);
    const double omega = c.triangle->omega, cp = gas.cp(), ex = 3.5;
    CHECK(std::abs(omega - (-2.0 * kPi * s.rpm / 60.0)) <= 1e-12 * std::abs(omega));
    // every node: r from the sheet, C_theta = W_theta + omega r, P0 = p (1 + C^2 / (2 cp T))^(g / (g - 1)), and the
    // isentropic identity of the inviscid solve P0 / P0_in = (T0 / T0_in)^(g / (g - 1))
    double worst_ct = 0.0, worst_p0 = 0.0, worst_iso = 0.0;
    for (int i = 0; i < d.grid.n_m(); ++i) {
        const double ri = s.chord * b.sheet->radius(d.grid.m[sz(i)]);
        CHECK(std::abs(d.station_radius[sz(i)] - ri) <= 1e-14 * ri);
        for (int j = 0; j < d.grid.n_eta(); ++j) {
            const double ct = d.W_theta(i, j) + omega * ri, C2 = d.W_m(i, j) * d.W_m(i, j) + ct * ct;
            const double T = d.temperature(i, j);
            worst_ct = std::max(worst_ct, std::abs(a.C_theta(i, j) - ct) / std::sqrt(C2));
            worst_p0 = std::max(worst_p0, std::abs(a.total_pressure(i, j) / (d.pressure(i, j) * std::pow(1.0 + C2 / (2.0 * cp * T), ex)) - 1.0));
            worst_iso = std::max(worst_iso, std::abs(a.total_pressure(i, j) / (P0 * std::pow(a.total_temperature(i, j) / T0, ex)) - 1.0));
        }
    }
    CHECK(sref::measure("rotating cone cascade: C_theta vs W_theta + omega r", worst_ct) <= 1e-13);
    CHECK(sref::measure("rotating cone cascade: P0 vs p (1 + C^2/(2 cp T))^3.5", worst_p0) <= 1e-12);
    CHECK(sref::measure("rotating cone cascade: P0/P0_in vs (T0/T0_in)^3.5 (isentropic)", worst_iso) <= 1e-12);
    // the rotor summary against mass averages formed here (weights rho W_m at the far exit)
    const int ie = d.grid.n_m() - 1, ne = d.grid.n_eta();
    std::vector<double> w(sz(ne)), x(sz(ne));
    for (int j = 0; j < ne; ++j) w[sz(j)] = d.density(ie, j) * d.W_m(ie, j);
    const double wsum = detail::trapezoid(w.data(), d.grid.eta);
    auto avg = [&](auto f) {
        for (int j = 0; j < ne; ++j) x[sz(j)] = w[sz(j)] * f(j);
        return detail::trapezoid(x.data(), d.grid.eta) / wsum;
    };
    const double r_out = d.station_radius.back(), U_out = omega * r_out;
    const double P0_out = avg([&](int j) {
        const double ct = d.W_theta(ie, j) + U_out, T = d.temperature(ie, j);
        return d.pressure(ie, j) * std::pow(1.0 + (d.W_m(ie, j) * d.W_m(ie, j) + ct * ct) / (2.0 * cp * T), ex);
    });
    const double ct_out = avg([&](int j) { return d.W_theta(ie, j) + U_out; });
    const double T0_out = avg([&](int j) {
        const double ct = d.W_theta(ie, j) + U_out;
        return d.temperature(ie, j) + (d.W_m(ie, j) * d.W_m(ie, j) + ct * ct) / (2.0 * cp);
    });
    const RotorSummary& o = *r.summary->rotor;
    CHECK(std::abs(o.total_pressure_ratio / (P0_out / P0) - 1.0) <= 1e-12);
    CHECK(std::abs(o.total_temperature_ratio / (T0_out / T0) - 1.0) <= 1e-12);
    CHECK(std::abs(o.C_theta_out / ct_out - 1.0) <= 1e-12);
    const double euler = U_out * ct_out - c.triangle->U * c.triangle->C_theta;
    CHECK(std::abs(o.euler_work / euler - 1.0) <= 1e-12);
    CHECK(std::abs(o.total_enthalpy_rise - cp * (T0_out - T0)) <= 1e-9 * std::abs(o.total_enthalpy_rise));
    CHECK(sref::measure("rotating cone cascade: Euler work vs cp dT0 (relative)", std::abs(o.euler_work / o.total_enthalpy_rise - 1.0)) <=
          1e-10);  // rothalpy: cp T0 = I + U C_theta at every node
    CHECK(o.total_pressure_ratio > 1.0);
    CHECK(o.C_theta_out < 0.0);  // swirl in the direction of rotation (toward -theta for rpm > 0)
    // Kelvin on the unbladed stations of the cascade: r C_theta constant ahead of the blade (zero: alpha1 = 0) and
    // behind it (the work done); the rotation source 2 S omega sin(phi) makes it so on the cone
    const auto in = circulation_drift(d.grid, d.W_theta, d.station_radius, omega, 0, d.grid.first_blade_i);
    const auto out = circulation_drift(d.grid, d.W_theta, d.station_radius, omega, d.grid.last_blade_i + 1, d.grid.n_m());
    CHECK(sref::measure("rotating cone cascade: |r C_theta| at the far inlet / (r U)", std::abs(in.second) / (r_in_nominal * std::abs(c.triangle->U))) <=
          1e-3);
    CHECK(sref::measure("rotating cone cascade: r C_theta drift, exit stations", out.first) <= 2e-3);
    // and at the far inlet its Delta to the exit is the Euler work over omega: U2 C2 - U1 C1 = omega (r2 C2 - r1 C1)
    CHECK(std::abs(omega * (r_out * ct_out) - euler) <= 1e-9 * std::abs(euler));
}

