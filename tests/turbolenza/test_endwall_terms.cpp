// The end-wall boundary layers of Aungier (2003) Ch. 8 (Endwall.hpp), term by term: the march
// against independent integrations of the book's differential equations (von Karman with the
// power-law closure; 8-5, 8-6, 8-7 with the edge forces 8-8, 8-9 in closed form on a conical
// rotating wall and through a blade row), xi of 8-53 against its formula, the continuity of the
// layer through V_theta_e = 0, the force defect as the turning vanishes, and every limit and option.
#include "SolverRef.hpp"

#include "turbolenza/BoundaryLayer.hpp"
#include "turbolenza/Endwall.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numbers>
#include <string>
#include <vector>

using namespace BladeGen::Turbolenza;

namespace {

constexpr double kPi = std::numbers::pi;

double rel(double a, double b) { return std::abs(a - b) / std::max(std::abs(b), 1e-300); }

EndwallWallStation wall(double m, double r, double Vm, double Ct, double rho, double P, double mu, double omega_w = 0.0,
                        double phi = 0.0)
{
    EndwallWallStation w;
    w.z = m * std::cos(phi);
    w.r = r;
    w.phi = phi;
    w.m = m;
    w.wall_omega = omega_w;
    w.Vm = Vm;
    w.Ctheta = Ct;
    w.rho = rho;
    w.P = P;
    w.mu = mu;
    return w;
}

// An edge given in closed form along one straight (conical) wall: r = r0 + m sin(phi) and the
// functions below with their m-derivatives.
struct Edge
{
    double r0, phi, omega_w;
    std::function<double(double)> Vm, dVm, rC, drC, rho, drho, P, dP;
    double mu;
    double r(double m) const { return r0 + m * std::sin(phi); }
    EndwallWallStation at(double m) const { return wall(m, r(m), Vm(m), rC(m) / r(m), rho(m), P(m), mu, omega_w, phi); }
};

using Y3 = std::array<double, 3>;  // (delta - delta1*, theta11, V_theta_e theta12)

// d/dm of the state from the book's 8-5, 8-6, 8-7 in differential form, with the closures of the
// module (tested elsewhere) but every term of the equations written here: the bracket of 8-6, the
// d(r C_theta)/dm term of 8-7 and the edge forces 8-8, 8-9 from the closed-form edge.
Y3 rhs(const Edge& e, double m, const Y3& y, double nu, double held_m, const EndwallSettings& set)
{
    const EndwallWallStation w = e.at(m);
    const EndwallLayer L = endwall_closures({y[0], y[1], y[2]}, w, e.omega_w, held_m, set);
    const double r = e.r(m), dr = std::sin(e.phi), rho = e.rho(m), drho = e.drho(m), V = e.Vm(m), dV = e.dVm(m);
    const double C = e.rC(m) / r, drC = e.drC(m);
    const double fm = rho * V * dV + e.dP(m) - std::sin(e.phi) / r * rho * C * C;  // [8-8], (V_theta + omega_w r) = C
    const double ft = rho * V / r * drC;                                           // [8-9]
    const double a1 = r * rho * V, da1 = dr * rho * V + r * drho * V + r * rho * dV;
    const double a2 = a1 * V, da2 = da1 * V + a1 * dV;
    const double a3 = r * a1, da3 = dr * a1 + r * da1;
    const double bracket =
        rho * L.Vtheta_e * std::sin(e.phi) * (L.Vtheta_e * (L.delta2 + L.theta22) + 2.0 * e.omega_w * r * L.delta2);
    Y3 d;
    d[0] = (r * rho * L.V_e * L.E - y[0] * da1) / a1;                                                       // [8-5]
    d[1] = (r * (L.tau_mw + fm * nu) + bracket - L.delta1 * r * rho * V * dV - y[1] * da2) / a2;             // [8-6]
    d[2] = (r * r * (L.tau_thetaw + ft * nu) - r * L.delta1 * rho * V * drC - y[2] * da3) / a3;              // [8-7]
    return d;
}

// classical RK4 from m0 to m1 in n steps
Y3 rk4(const Edge& e, double m0, double m1, int n, Y3 y, double nu, double held_m, const EndwallSettings& set)
{
    const double h = (m1 - m0) / n;
    auto add = [](const Y3& a, const Y3& b, double k) { return Y3{a[0] + k * b[0], a[1] + k * b[1], a[2] + k * b[2]}; };
    for (int i = 0; i < n; ++i) {
        const double m = m0 + h * i;
        const Y3 k1 = rhs(e, m, y, nu, held_m, set);
        const Y3 k2 = rhs(e, m + 0.5 * h, add(y, k1, 0.5 * h), nu, held_m, set);
        const Y3 k3 = rhs(e, m + 0.5 * h, add(y, k2, 0.5 * h), nu, held_m, set);
        const Y3 k4 = rhs(e, m + h, add(y, k3, h), nu, held_m, set);
        for (int j = 0; j < 3; ++j) y[static_cast<std::size_t>(j)] += h / 6.0 * (k1[static_cast<std::size_t>(j)] +
                                                                                 2.0 * k2[static_cast<std::size_t>(j)] +
                                                                                 2.0 * k3[static_cast<std::size_t>(j)] +
                                                                                 k4[static_cast<std::size_t>(j)]);
    }
    return y;
}

Y3 state3(const EndwallLayer& L) { return {L.Q, L.theta11, L.q12}; }

}  // namespace

TEST_CASE("end wall terms: a decelerating edge is von Karman's equation with the power-law closure", "[turbolenza][endwall]")
{
    // V = 60 (1 - 0.25 m), no swirl, straight annulus: 8-5 and 8-6 reduce to
    //   (delta - delta1*)' = E - (delta - delta1*) V'/V,  theta' = c_f/2 - (2 + H1) theta V'/V
    // with H1 = Q/(Q - 2 theta) [8-23], E = 0.025 (H1 - 1) [8-43, H = H1], c_f of 8-36; integrated here by RK4
    const double rho = 1.2, mu = 1.8e-5, L = 1.0;
    const int N = 401;
    EndwallInput in;
    for (int k = 0; k < N; ++k) {
        const double m = L * k / (N - 1), V = 60.0 * (1.0 - 0.25 * m);
        EndwallStation s;
        s.hub = wall(m, 0.5, V, 0.0, rho, 1e5, mu);
        s.shroud = wall(m, 0.6, V, 0.0, rho, 1e5, mu);
        in.stations.push_back(s);
    }
    in.mass_flow = rho * 60.0 * kPi * (0.36 - 0.25);
    in.inlet.blockage = 1e-9;
    const EndwallResult r = analyze_endwalls(in);
    double Q = r.stations[0].hub.Q, th = r.stations[0].hub.theta11;
    const int n = 20000;
    const double h = L / n;
    auto f = [&](double m, double q, double t, double& dq, double& dt) {
        const double V = 60.0 * (1.0 - 0.25 * m), dV = -15.0;
        const double H1 = q / (q - 2.0 * t);
        const double cf = 0.246 * std::exp(-1.561 * H1) * std::pow(std::max(rho * V * t / mu, 250.0), -0.268);
        dq = 0.025 * (H1 - 1.0) - q * dV / V;
        dt = 0.5 * cf - (2.0 + H1) * t * dV / V;
    };
    for (int i = 0; i < n; ++i) {
        const double m = h * i;
        double q1, t1, q2, t2, q3, t3, q4, t4;
        f(m, Q, th, q1, t1);
        f(m + 0.5 * h, Q + 0.5 * h * q1, th + 0.5 * h * t1, q2, t2);
        f(m + 0.5 * h, Q + 0.5 * h * q2, th + 0.5 * h * t2, q3, t3);
        f(m + h, Q + h * q3, th + h * t3, q4, t4);
        Q += h / 6.0 * (q1 + 2.0 * q2 + 2.0 * q3 + q4);
        th += h / 6.0 * (t1 + 2.0 * t2 + 2.0 * t3 + t4);
    }
    const EndwallLayer& e = r.stations.back().hub;
    CHECK(!e.H1_limited);
    const double d_th = sref::measure("endwall terms decel: theta11(L) march / RK4 - 1", e.theta11 / th - 1.0);
    const double d_q = sref::measure("endwall terms decel: (delta - delta1*)(L) march / RK4 - 1", e.Q / Q - 1.0);
    CHECK(std::abs(d_th) <= 1e-4);  // measured 9.8e-6 (the march's O(dm^2), dm = 2.5e-3)
    CHECK(std::abs(d_q) <= 1e-4);   // measured 8.4e-6
    sref::measure("endwall terms decel: H1(L)", e.H1);
    // and the Sec. 5.6 march (Head / Green, another closure) on the same edge
    EdgeConditions ec;
    for (const EndwallStation& st : in.stations) {
        ec.x.push_back(st.hub.m);
        ec.W.push_back(st.hub.Vm);
    }
    ec.rho.assign(ec.x.size(), rho);
    ec.T.assign(ec.x.size(), 300.0);
    ec.mu.assign(ec.x.size(), mu);
    ec.cp = 1e12;
    ec.total_temperature = 300.0;
    const BoundaryLayerProfile bl = march_surface(ec, {});
    const double d_bl = sref::measure("endwall terms decel: theta11(L) / Sec. 5.6 theta(L) - 1", e.theta11 / bl.theta.back() - 1.0);
    CHECK(std::abs(d_bl) <= 0.07);  // measured 3.5e-2: two closures
}

namespace {

// A rotor on a conical hub (phi = asin 0.4, drum and blades at 300 rad/s) and a stationary conical
// casing: m in [0, 0.1] blade-free, the row in [0.1, 0.2] with 200 stations inside.
struct Conical
{
    Edge hub, shroud;
    EndwallInput in;
};

Conical conical(int n_free, int n_row, double clearance)
{
    Conical c;
    const double phi = std::asin(0.4), omega = 300.0;
    auto Vm = [](double m) { return 120.0 - 100.0 * m; };
    auto dVm = [](double) { return -100.0; };
    auto rho = [](double m) { return 1.2 + 0.5 * m; };
    auto drho = [](double) { return 0.5; };
    auto P = [](double m) { return 1e5 + 2e4 * m + 1e5 * m * m; };
    auto dP = [](double m) { return 2e4 + 2e5 * m; };
    c.hub = {0.30, phi, omega, Vm, dVm, [](double m) { return 15.0 + 40.0 * m + 150.0 * m * m; },
             [](double m) { return 40.0 + 300.0 * m; }, rho, drho, P, dP, 1.8e-5};
    c.shroud = {0.50, phi, 0.0, Vm, dVm, [](double m) { return 25.0 + 60.0 * m + 250.0 * m * m; },
                [](double m) { return 60.0 + 500.0 * m; }, rho, drho, P, dP, 1.8e-5};
    const int N = n_free + n_row + 1;
    for (int k = 0; k < N; ++k) {
        const double m = k <= n_free ? 0.1 * k / n_free : 0.1 + 0.1 * (k - n_free) / n_row;
        EndwallStation s;
        s.hub = c.hub.at(m);
        s.shroud = c.shroud.at(m);
        s.quasi_normal_length = 0.2;
        c.in.stations.push_back(s);
    }
    c.in.mass_flow = 1.2 * 120.0 * kPi * (0.25 - 0.09);
    c.in.compressor_omega = omega;
    c.in.inlet.blockage = 0.1;
    EndwallBladeRow row;
    row.station_in = n_free;
    row.station_out = N - 1;
    const double r_in_h = c.hub.r(0.1), r_out_h = c.hub.r(0.2), r_in_s = c.shroud.r(0.1), r_out_s = c.shroud.r(0.2);
    row.hub = {omega, 40, kPi * (r_in_h + r_out_h) / 40.0, 0.4, 0.0, false, 0.0};
    row.shroud = {omega, 40, kPi * (r_in_s + r_out_s) / 40.0, 0.9, clearance, false, 0.0};
    c.in.rows.push_back(row);
    return c;
}

}  // namespace

TEST_CASE("end wall terms: 8-5 .. 8-9 on a conical rotating wall and through a row, against RK4", "[turbolenza][endwall]")
{
    const Conical c = conical(200, 200, 1e-3);
    const EndwallResult r = analyze_endwalls(c.in);
    const int k_le = 200, k_te = 400;
    for (bool hub : {true, false}) {
        const Edge& e = hub ? c.hub : c.shroud;
        auto L = [&](int k) -> const EndwallLayer& {
            const EndwallStationResult& s = r.stations[static_cast<std::size_t>(k)];
            return hub ? s.hub : s.shroud;
        };
        const EndwallRowWallResult& rw = hub ? r.rows[0].hub : r.rows[0].shroud;
        for (int k = 0; k <= k_te; ++k) {
            CHECK(L(k).Vtheta_e != 0.0);
            CHECK((!L(k).H1_limited && !L(k).H2_limited && !L(k).H_limited && !L(k).re_theta11_limited));
        }
        CHECK(!rw.dnu_limited);
        // the blade-free cone: 8-5 .. 8-7 with nu = 0 (the 8-6 bracket with its Coriolis term, the 8-7 d(rC)/dm term)
        const Y3 y1 = rk4(e, 0.0, 0.1, 4000, state3(L(0)), 0.0, L(0).m, c.in.settings);
        const std::string w = hub ? "hub" : "casing";
        double worst = 0.0;
        for (int j = 0; j < 3; ++j) {
            const double ref = y1[static_cast<std::size_t>(j)], got = state3(L(k_le))[static_cast<std::size_t>(j)];
            worst = std::max(worst, std::abs(got - ref) / std::abs(ref));
        }
        CHECK(sref::measure("endwall terms cone " + w + ": march / RK4 - 1 at the LE (worst of Q, theta11, q12)", worst) <=
              3e-6);  // measured 1.3e-7 hub, 2.4e-7 casing
        // the row: the same with the row's nu and the closed-form f_me (8-8) and f_theta_e (8-9)
        const Y3 y2 = rk4(e, 0.1, 0.2, 4000, state3(L(k_le)), rw.nu, L(k_le).m, c.in.settings);
        worst = 0.0;
        for (int j = 0; j < 3; ++j) {
            const double ref = y2[static_cast<std::size_t>(j)], got = state3(L(k_te))[static_cast<std::size_t>(j)];
            worst = std::max(worst, std::abs(got - ref) / std::abs(ref));
        }
        CHECK(sref::measure("endwall terms row " + w + ": march / RK4 - 1 at the TE (worst of Q, theta11, q12)", worst) <=
              3e-6);  // measured 2.7e-7 hub, 3.6e-7 casing
        // the force-defect terms are not small here (else the check above would not see them)
        const Y3 y0 = rk4(e, 0.1, 0.2, 4000, state3(L(k_le)), 0.0, L(k_le).m, c.in.settings);
        const double effect = sref::measure("endwall terms row " + w + ": theta11(TE) with nu / without - 1",
                                            y2[1] / y0[1] - 1.0);
        CHECK(std::abs(effect) > 1e-2);  // measured 0.21 hub, -0.031 casing
        sref::measure("endwall terms row " + w + ": nu / delta1*(TE)", rw.nu / L(k_te).delta1);
    }
}

TEST_CASE("end wall terms: xi is Eq. 8-53 on a row with changing radius", "[turbolenza][endwall]")
{
    // one step across the row (the Ch. 7 case): r_in != r_out, W_theta_in != W_theta_out
    const Conical c0 = conical(4, 1, 1e-3);
    for (bool hub : {true, false}) {
        const EndwallWallStation& a = hub ? c0.in.stations[4].hub : c0.in.stations[4].shroud;
        const EndwallWallStation& b = hub ? c0.in.stations[5].hub : c0.in.stations[5].shroud;
        const EndwallRowWall& rw = hub ? c0.in.rows[0].hub : c0.in.rows[0].shroud;
        const EndwallResult r = analyze_endwalls(c0.in);
        const EndwallRowWallResult& res = hub ? r.rows[0].hub : r.rows[0].shroud;
        REQUIRE(a.r != b.r);
        // 8-53 with rbar f_theta_e dm = int r f_theta_e dm = rho_mean Vm_mean [(r C)_out - (r C)_in] (8-9 over the step)
        const double rbar = 0.5 * (a.r + b.r);
        const double Wi = a.Ctheta - a.r * rw.omega, Wo = b.Ctheta - b.r * rw.omega;
        const double Wbar = 0.5 * (Wi + Wo), Vbar = 0.5 * (a.Vm + b.Vm);
        const double rf_dm = 0.5 * (a.rho + b.rho) * Vbar * (b.r * b.Ctheta - a.r * a.Ctheta);
        const double xi = (rbar / b.r + Wbar * Wo / (Vbar * b.Vm)) * rf_dm / (b.r * b.rho * b.Vm);
        CHECK(rel(res.xi, xi) < 1e-13);
        CHECK(std::abs(rbar / b.r - 1.0) > 1e-3);  // the first term is not 1
        CHECK(std::abs(Wi / Wo - 1.0) > 1e-2);
        sref::measure(hub ? "endwall terms xi hub" : "endwall terms xi casing", res.xi);
    }
}

TEST_CASE("end wall terms: the layer is continuous through V_theta_e = 0", "[turbolenza][endwall]")
{
    // a hub rotating at 120 rad/s, r = 0.25: V_theta_e = C_theta - 30 at station b; q12 != 0 there
    // (r omega = 30 exactly, so V_theta_e = 0 exactly also where the compiler fuses the multiply-add)
    auto run = [](double Ct1) {
        const EndwallWallStation a = wall(0.0, 0.25, 100.0, 60.0, 1.2, 1e5, 1.8e-5, 120.0);
        const EndwallWallStation b = wall(0.01, 0.25, 100.0, Ct1, 1.2, 1e5, 1.8e-5, 120.0, 0.3);
        const EndwallWallStation c = wall(0.02, 0.25 + 0.01 * std::sin(0.3), 100.0, 30.0, 1.2, 1e5, 1.8e-5, 120.0, 0.3);
        const EndwallLayer L0 = endwall_inlet(a, {}, 20.0, {});
        const EndwallStepResult s1 = endwall_step(L0, a, b, nullptr, 0.0, {});
        const EndwallStepResult s2 = endwall_step(s1.layer, b, c, nullptr, 0.0, {});
        return std::array<EndwallLayer, 2>{s1.layer, s2.layer};
    };
    const auto z = run(30.0);
    REQUIRE(z[0].Vtheta_e == 0.0);
    REQUIRE(z[0].q12 != 0.0);
    CHECK(std::isnan(z[0].theta12));
    double worst = 0.0;
    for (double eps : {1e-3, 1e-6, 1e-9}) {
        const auto lo = run(30.0 - eps), hi = run(30.0 + eps);
        // the 8-6 bracket's products at b vanish from both sides (V_theta_e delta2*, V_theta_e^2 theta22)
        for (const auto& s : {lo, hi}) {
            CHECK(std::abs(s[0].Vtheta_e * s[0].delta2) <= 5.0 * eps * 1e-2);
            CHECK(std::abs(s[0].Vtheta_e * s[0].Vtheta_e * s[0].theta22) <= 5.0 * eps * 1e-2);
        }
        const double jump = std::max({rel(lo[1].theta11, z[1].theta11), rel(hi[1].theta11, z[1].theta11),
                                      rel(lo[1].q12, z[1].q12), rel(hi[1].q12, z[1].q12), rel(lo[0].theta11, z[0].theta11),
                                      rel(hi[0].theta11, z[0].theta11)});
        worst = std::max(worst, jump / eps);
        CHECK(jump <= 1e-2 * eps + 1e-14);  // measured 9.2e-4 eps (was 3.3e-3 for eps = 1e-6 before the m limit)
    }
    sref::measure("endwall terms continuity at V_theta_e = 0: max relative change / dC_theta [1/(m/s)]", worst);
}

TEST_CASE("end wall terms: the force defect vanishes with the blade force", "[turbolenza][endwall]")
{
    // a rotor with almost no turning: xi -> 0; 8-57 alone would give dnu ~ 1/xi; with 12-50's bound dnu f -> 0
    auto run = [](double turn) {
        EndwallInput in;
        const double omega = 800.0;
        for (int k = 0; k < 5; ++k) {
            const double m = 0.025 * k, t = k <= 1 ? 0.0 : 1.0;
            EndwallStation s;
            s.hub = wall(m, 0.30, 150.0, turn * t, 1.2, 1e5, 1.8e-5, (k == 2) ? omega : 0.0);
            s.shroud = wall(m, 0.45, 150.0, turn * t, 1.2, 1e5, 1.8e-5);
            in.stations.push_back(s);
        }
        in.mass_flow = 1.2 * 150.0 * kPi * (0.45 * 0.45 - 0.09);
        in.compressor_omega = omega;
        EndwallBladeRow row;
        row.station_in = 1;
        row.station_out = 2;
        row.hub = {omega, 30, kPi * 0.6 / 30.0, 0.3, 0.0, false, 0.0};
        row.shroud = {omega, 30, kPi * 0.9 / 30.0, 0.9, 1e-3, false, 0.0};
        in.rows.push_back(row);
        return analyze_endwalls(in);
    };
    const EndwallResult r0 = run(0.0);
    CHECK(r0.rows[0].shroud.xi == 0.0);
    CHECK(r0.rows[0].shroud.dnu == 0.0);
    double last = 1.0;
    for (double turn : {10.0, 1.0, 0.1, 0.01, 1e-4}) {
        const EndwallResult r = run(turn);
        const EndwallRowWallResult& w = r.rows[0].shroud;
        CHECK(std::abs(w.dnu) <= w.dnu_bound * (1.0 + 1e-15));
        const double d = rel(r.stations[2].shroud.theta11, r0.stations[2].shroud.theta11);
        sref::measure("endwall terms casing theta11(TE) / no turning - 1 at C_theta " + std::to_string(turn), d);
        CHECK(d <= last);  // the layer tends to the turning-free one
        last = d;
        if (turn <= 0.01) CHECK(w.dnu_limited);
    }
    CHECK(last <= 1e-12);  // measured 9.3e-15 at C_theta = 1e-4 m/s (and ~ turn^2 above)
}

TEST_CASE("end wall terms: the limits and options of the closures", "[turbolenza][endwall]")
{
    const EndwallSettings set;
    const EndwallWallStation w = wall(0.0, 0.3, 100.0, 80.0, 1.2, 1e5, 1.8e-5);
    const double t = 1e-3;
    // H1 <= 2.4 (p. 184): delta - delta1* = 3 theta11 gives H1 = 3 by 8-23
    const EndwallLayer h1 = endwall_closures({3.0 * t, t, 80.0 * 0.1 * t}, w, 0.0, 1.0 / 7.0, set);
    CHECK(h1.H1_limited);
    CHECK(h1.H1 == 2.4);
    CHECK(rel(h1.delta1, 2.4 * t) < 1e-15);
    // a regular state: Q = 7 theta11 (H1 = 1.4, n = 0.2, delta = 8.4 theta11)
    const double delta = 8.4 * t, n = 0.2;
    // H2 <= 2.4: theta12 = 2 delta/(n+1) > 0.99 delta/(n+1) is beyond the step-4 precaution and m = 0.7
    const EndwallLayer h2 = endwall_closures({7.0 * t, t, 80.0 * 2.0 * delta / (n + 1.0)}, w, 0.0, 1.0 / 7.0, set);
    CHECK(rel(h2.H1, 1.4) < 1e-14);
    CHECK(h2.H2_limited);
    CHECK(h2.m == 0.7);
    CHECK(h2.theta12_limited);
    // 8-41 gets the precaution-limited theta12 (p. 195 step 4)
    const EndwallStreamwise sw = endwall_streamwise(h2.delta1, h2.delta2, h2.theta11, 0.99 * h2.delta / (h2.n + 1.0),
                                                    h2.theta22, h2.beta_e);
    CHECK(rel(h2.theta, sw.theta) < 1e-14);
    // the precaution alone (H2 limit far away): m of the limited theta12 by 8-21, 99 (n + 1)
    EndwallSettings wide;
    wide.shape_factor_limit = 1000.0;
    const EndwallLayer pre = endwall_closures({7.0 * t, t, 80.0 * 2.0 * delta / (n + 1.0)}, w, 0.0, 1.0 / 7.0, wide);
    CHECK(pre.theta12_limited);
    CHECK(!pre.H2_limited);
    CHECK(rel(pre.m, 99.0 * (n + 1.0)) < 1e-10);
    // the lower limit (ours): theta12 very negative, m held at -0.85 (H2 = 0.40, the end of Fig. 8-4)
    const EndwallLayer lo = endwall_closures({7.0 * t, t, -80.0 * 5.0 * delta}, w, 0.0, 1.0 / 7.0, set);
    CHECK(lo.H2_limited);
    CHECK(lo.m == -0.85);
    CHECK(rel(lo.H2, endwall_H2(-0.85)) < 1e-15);
    // Re_theta >= 250 in 8-36 (p. 186): a thin layer, rho V theta / mu = 67
    const EndwallLayer re = endwall_closures({7e-4 * 0.1, 1e-5, 0.0}, wall(0.0, 0.3, 100.0, 0.0, 1.2, 1e5, 1.8e-5), 0.0,
                                             1.0 / 7.0, set);
    CHECK(re.re_theta_limited);
    CHECK(re.Re_theta < 250.0);
    CHECK(re.cf == ludwieg_tillmann_cf(re.H, 250.0));
    // entrainment: 8-43 (E of H, default) or 8-42 (E of H1) on a skewed layer, H != H1
    const EndwallLayer sk = endwall_closures({7.0 * t, t, 80.0 * 0.6 * t}, w, 0.0, 1.0 / 7.0, set);
    REQUIRE(std::abs(sk.H - sk.H1) > 1e-3);
    CHECK(sk.E == head_entrainment(sk.H));
    EndwallSettings e42;
    e42.entrainment = EndwallEntrainment::MeridionalH1;
    const EndwallLayer sk42 = endwall_closures({7.0 * t, t, 80.0 * 0.6 * t}, w, 0.0, 1.0 / 7.0, e42);
    CHECK(sk42.E == head_entrainment(sk42.H1));
    CHECK(sk42.E != sk.E);
    CHECK(std::string(endwall_entrainment_name(e42.entrainment)) == "meridional_H1");
    CHECK(std::string(endwall_force_defect_name(EndwallForceDefect::None)) == "none");
    // H <= 2.4 (p. 186): with H1 and H2 limited to 2.4 the streamwise H stays below it on this scan (measured
    // 2.39), so the limit is a guard (theta <= 0 cannot occur after the >= 0 limits for the profile family)
    double h_max = 0.0;
    for (double q : {2.05, 2.5, 3.0, 5.0})
        for (double s12 : {-5.0, -1.0, 0.1, 0.3, 0.6, 2.0})
            for (double ct : {10.0, 80.0, 400.0}) {
                const EndwallLayer x = endwall_closures({q * t, t, ct * s12 * t}, wall(0.0, 0.3, 100.0, ct, 1.2, 1e5, 1.8e-5),
                                                        0.0, 1.0 / 7.0, set);
                h_max = std::max(h_max, x.delta_star / x.theta);
                if (x.H_limited) CHECK(x.H == 2.4);
            }
    sref::measure("endwall terms largest delta*/theta over the scan (H limit 2.4)", h_max);
}

TEST_CASE("end wall terms: 8-73 acts after a step in an accelerating flow", "[turbolenza][endwall]")
{
    EndwallInput in;
    const int N = 101;
    for (int k = 0; k < N; ++k) {
        const double m = 0.2 * k / (N - 1), V = 20.0 + 500.0 * m;
        EndwallStation s;
        s.hub = wall(m, 0.3, V, 0.0, 1.2, 1e5, 1.8e-5);
        s.shroud = wall(m, 0.4, V, 0.0, 1.2, 1e5, 1.8e-5);
        in.stations.push_back(s);
    }
    in.mass_flow = 1.2 * 20.0 * kPi * (0.16 - 0.09);
    in.inlet.blockage = 1e-9;
    const EndwallResult r = analyze_endwalls(in);
    int limited = 0;
    for (std::size_t k = 1; k < r.stations.size(); ++k) {
        const EndwallLayer& L = r.stations[k].hub;
        CHECK(L.Re_theta11 >= 250.0 * (1.0 - 1e-13));
        if (L.re_theta11_limited) {
            ++limited;
            CHECK(rel(L.Re_theta11, 250.0) < 1e-13);
        }
    }
    sref::measure("endwall terms accelerating flow: stations held at Re_theta11 = 250", limited);
    CHECK(limited > 0);
}

TEST_CASE("end wall terms: seal leakage on a rotating shroud band", "[turbolenza][endwall]")
{
    // a shrouded rotor (casing side): the shroud band and its wall segments turn with the blades at 500 rad/s
    EndwallInput in;
    const double omega = 500.0, leak = 0.03;
    for (int k = 0; k < 5; ++k) {
        const double m = 0.03 * k, t = k <= 1 ? 0.0 : 1.0;
        EndwallStation s;
        s.hub = wall(m, 0.30, 120.0 - 8.0 * t, 40.0 * t, 1.2, 1e5 + 5e3 * t, 1.8e-5, (k == 2) ? omega : 0.0);
        s.shroud = wall(m, 0.45, 120.0 - 8.0 * t, 30.0 * t, 1.2, 1e5 + 5e3 * t, 1.8e-5, (k == 2) ? omega : 0.0);
        in.stations.push_back(s);
    }
    in.mass_flow = 1.2 * 120.0 * kPi * (0.45 * 0.45 - 0.09);
    in.compressor_omega = omega;
    EndwallBladeRow row;
    row.station_in = 1;
    row.station_out = 2;
    row.hub = {omega, 30, kPi * 0.6 / 30.0, 0.3, 0.0, false, 0.0};
    row.shroud = {omega, 30, kPi * 0.9 / 30.0, 0.9, 0.0, true, leak};
    in.rows.push_back(row);
    const EndwallResult r = analyze_endwalls(in);
    const EndwallRowWallResult& s = r.rows[0].shroud;
    const EndwallWallStation &a = in.stations[1].shroud, &b = in.stations[2].shroud;
    const double Vt_a = a.Ctheta - a.r * omega, Vt_b = b.Ctheta - b.r * omega;  // relative to the band (p. 193)
    CHECK(rel(s.leakage.U_leak_in, 0.5 * a.r * omega - a.r * omega) < 1e-14);  // 8-61: r omega_c / 2 - r omega
    CHECK(rel(s.leakage.U_leak_out, 0.5 * b.r * omega - b.r * omega) < 1e-14);
    auto mom_t = [](const EndwallState& x, const EndwallWallStation& w, double Vt) {
        return 2.0 * kPi * w.r * w.rho * w.Vm * (Vt * x.Q - x.q12);  // int rho V_m V_theta dy (relative to the band)
    };
    CHECK(std::abs(mom_t(s.leakage.plus, a, Vt_a) - mom_t(s.leakage.in, a, Vt_a) - leak * s.leakage.U_leak_in) < 1e-10);
    CHECK(std::abs(mom_t(s.leakage.out, b, Vt_b) - mom_t(s.leakage.minus, b, Vt_b)) < 1e-10);
    CHECK(rel(s.leakage.plus.q12 - s.leakage.in.q12, s.leakage.delta_plus * (Vt_a - s.leakage.U_leak_in)) < 1e-12);
    // a clearance on a shrouded wall is refused (its leakage is the seal's)
    EndwallInput bad = in;
    bad.rows[0].shroud.clearance = 1e-3;
    CHECK_THROWS_AS(analyze_endwalls(bad), EndwallError);
    CHECK_THROWS_AS(analyze_endwalls(bad), std::invalid_argument);
}
