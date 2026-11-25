// The end-wall boundary layers of Aungier (2003) Ch. 8 (Endwall.hpp): the profile relations against
// the integrals of their profiles, the streamwise thicknesses against their definitions, the
// closures against BoundaryLayer's, the straight annulus against the 1/7-power flat plate and the
// Sec. 5.6 march, the order of the integration, zero skew, the jump conditions, K_B against its
// definition, the seal leakage balances, the blade force defect, Figs. 8-4 and 8-6, the refusals.
#include "SolverRef.hpp"

#include "turbolenza/BoundaryLayer.hpp"
#include "turbolenza/Endwall.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace BladeGen::Turbolenza;

namespace {

constexpr double kPi = std::numbers::pi;

// Euler's Beta function: int_0^1 x^(a-1) (1-x)^(b-1) dx
double beta_fn(double a, double b) { return std::tgamma(a) * std::tgamma(b) / std::tgamma(a + b); }

double rel(double a, double b) { return std::abs(a - b) / std::max(std::abs(b), 1e-300); }

EndwallWallStation wall(double m, double r, double Vm, double Ct, double rho, double P, double mu, double omega_w = 0.0,
                        double phi = 0.0)
{
    EndwallWallStation w;
    w.z = m * std::cos(phi);  // a straight (conical) wall: dz/dm = cos(phi), dr/dm = sin(phi)
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

// a straight annulus r_h, r_s; uniform edge (Vm, Ct, rho), stationary walls
EndwallInput annulus(int N, double L, double rh, double rs, double Vm, double Ct, double rho, double mu)
{
    EndwallInput in;
    for (int k = 0; k < N; ++k) {
        const double m = L * k / (N - 1);
        EndwallStation s;
        s.hub = wall(m, rh, Vm, Ct, rho, 1e5, mu);
        s.shroud = wall(m, rs, Vm, Ct, rho, 1e5, mu);
        s.quasi_normal_length = rs - rh;
        in.stations.push_back(s);
    }
    in.mass_flow = rho * Vm * kPi * (rs * rs - rh * rh);
    return in;
}

EndwallState state_of(const EndwallLayer& L) { return {L.Q, L.theta11, L.q12}; }

// the extended profile 8-28: f = eta^0.05 + c (1 - eta)^2 eta^0.1, c = 0.1705 (1 - 20 m)
struct Extended
{
    double delta2, theta22, theta12;
};
Extended extended_integrals(double m, double n)
{
    const double c = 0.1705 * (1.0 - 20.0 * m);
    const double int_f = 1.0 / 1.05 + c * beta_fn(1.1, 3.0);
    const double int_f2 = 1.0 / 1.1 + 2.0 * c * beta_fn(1.15, 3.0) + c * c * beta_fn(1.2, 5.0);
    const double int_uf = 1.0 / (n + 1.05) + c * beta_fn(n + 1.1, 3.0);  // u = eta^n
    return {1.0 - int_f, int_f - int_f2, 1.0 / (n + 1.0) - int_uf};
}

}  // namespace

TEST_CASE("end wall: the power-law relations 8-20 .. 8-26 are the integrals of 8-18, 8-19", "[turbolenza][endwall]")
{
    double worst = 0.0;
    for (double n : {0.1, 0.2, 0.4, 0.7})
        for (double m : {0.05, 0.1, 1.0 / 7.0, 0.3, 0.7}) {
            const double delta = 2e-3;
            const EndwallProfile p = endwall_profile(delta, n, m);
            // int_0^1 of the power laws (constant density)
            const double d1 = delta * (1.0 - 1.0 / (n + 1.0));
            const double t11 = delta * (1.0 / (n + 1.0) - 1.0 / (2.0 * n + 1.0));
            const double t12 = delta * (1.0 / (n + 1.0) - 1.0 / (n + m + 1.0));
            const double d2 = delta * (1.0 - 1.0 / (m + 1.0));
            const double t22 = delta * (1.0 / (m + 1.0) - 1.0 / (2.0 * m + 1.0));
            worst = std::max({worst, rel(p.delta1, d1), rel(p.theta11, t11), rel(p.theta12, t12), rel(p.delta2, d2),
                              rel(p.theta22, t22)});
            worst = std::max(worst, rel(t11 / (delta - d1 - 2.0 * t11), n));                         // [8-20]
            worst = std::max(worst, rel(endwall_m_of_theta12(t12, delta, n), m));                   // [8-21]
            worst = std::max(worst, rel(2.0 * p.H1 * t11 / (p.H1 - 1.0), delta - d1));               // [8-23]
            worst = std::max(worst, rel(d2 / t22, endwall_H2(m)));                                  // [8-25]
        }
    CHECK(sref::measure("endwall power-law relations vs profile integrals", worst) <= 1e-13);  // measured 3.1e-15
}

TEST_CASE("end wall: the extended profile 8-28 and its relations 8-29 .. 8-32", "[turbolenza][endwall]")
{
    const double n = 0.2, delta = 1.0;
    double worst_exact = 0.0;
    for (double m : {0.05, 0.02, 0.0, -0.1, -0.3, -0.6, -0.85}) {
        const EndwallProfile p = endwall_profile(delta, n, m);
        const Extended x = extended_integrals(m, n);
        worst_exact = std::max(worst_exact, std::abs(p.delta2 - x.delta2));    // [8-29] exact
        worst_exact = std::max(worst_exact, std::abs(p.theta12 - x.theta12));  // [8-31] exact
        worst_exact = std::max(worst_exact, std::abs(endwall_m_of_theta12(x.theta12, delta, n) - m));
    }
    CHECK(sref::measure("endwall 8-29, 8-31 vs the integrals of 8-28", worst_exact) <= 1e-13);  // measured 3.3e-16
    // 8-30 as printed is not the integral: the exact theta22/delta is a quadratic in m (fit it at three m)
    const double y0 = extended_integrals(0.0, n).theta22, ym = extended_integrals(-0.5, n).theta22,
                 yp = extended_integrals(0.05, n).theta22;
    // solve y = c0 + c1 m + c2 m^2 through (0, y0), (-0.5, ym), (0.05, yp)
    const double c0 = y0;
    const double c2 = ((ym - c0) / -0.5 - (yp - c0) / 0.05) / (-0.5 - 0.05);
    const double c1 = (yp - c0) / 0.05 - c2 * 0.05;
    sref::measure("endwall exact theta22/delta: constant", c0);
    sref::measure("endwall exact theta22/delta: m coefficient", c1);
    sref::measure("endwall exact theta22/delta: m^2 coefficient", c2);
    CHECK(std::abs(c0 + 4.38e-4) < 2e-6);
    CHECK(std::abs(c1 - 0.9502) < 2e-4);
    CHECK(std::abs(c2 + 1.5126) < 2e-4);
    // the printed fit through m = 0 and m = 0.05: equal to the integral at m = 0.05 (to its 4 digits)
    CHECK(std::abs(endwall_profile(1.0, n, 0.0499999999).theta22 - yp) < 2e-6);
    const double d03 = endwall_profile(1.0, n, -0.3).theta22 / extended_integrals(-0.3, n).theta22 - 1.0;
    const double d085 = endwall_profile(1.0, n, -0.85).theta22 / extended_integrals(-0.85, n).theta22 - 1.0;
    sref::measure("endwall 8-30 printed / integral - 1 at m = -0.3", d03);
    sref::measure("endwall 8-30 printed / integral - 1 at m = -0.85", d085);
    CHECK(std::abs(d03 - 0.035) < 0.005);   // measured 3.5e-2
    CHECK(std::abs(d085 - 0.065) < 0.005);  // measured 6.5e-2
    // continuity at m = 0.05 (8-21 / 8-31 and 8-25 / 8-32)
    CHECK(std::abs(endwall_H2(0.05) - 1.1) < 1e-15);
    // 8-32 at m = 0.05 is 1.1000011: its rounded 1.684 leaves a step of 1.1e-6 against 8-25 (measured)
    CHECK(std::abs(sref::measure("endwall 8-32 at m = 0.05 minus 1.1", endwall_H2(0.05 - 1e-12) - 1.1)) < 2e-6);
    const double t_sw = 0.05 / ((n + 1.0) * (n + 1.05));
    CHECK(std::abs(endwall_m_of_theta12(t_sw * (1.0 + 1e-12), 1.0, n) - 0.05) < 1e-10);
    CHECK(std::abs(endwall_m_of_theta12(t_sw * (1.0 - 1e-12), 1.0, n) - 0.05) < 1e-10);
    CHECK(std::isinf(endwall_m_of_theta12(1.0 / (n + 1.0), 1.0, n)));  // the pole of 8-21
    // Fig. 8-4 (p. 185), read from the page: H2 = 0.40 at the curve's end m = -0.85, H2 = 0.80 at m = -0.14,
    // H2 = 0.60 at m = -0.38, the power law through (2.2, 0.6); reading uncertainty +-0.02 in H2
    CHECK(std::abs(endwall_H2(-0.85) - 0.40) < 0.02);
    CHECK(std::abs(endwall_H2(-0.14) - 0.80) < 0.02);
    CHECK(std::abs(endwall_H2(-0.38) - 0.60) < 0.02);
    CHECK(std::abs(endwall_H2(0.6) - 2.2) < 1e-15);
}

TEST_CASE("end wall: the streamwise thicknesses 8-40, 8-41 are the integrals 8-38, 8-39", "[turbolenza][endwall]")
{
    double worst = 0.0;
    for (double beta : {0.0, 0.3, -0.6, 1.2})
        for (double m : {0.1, 0.2, 0.45}) {
            const double n = 0.2, delta = 1.0;
            const EndwallProfile p = endwall_profile(delta, n, m);
            const EndwallStreamwise s = endwall_streamwise(p.delta1, p.delta2, p.theta11, p.theta12, p.theta22, beta);
            // V_s / V_e = cos^2 u + sin^2 w [8-37], u = eta^n, w = eta^m
            const double a = std::cos(beta) * std::cos(beta), b = std::sin(beta) * std::sin(beta);
            const double int_s = a / (n + 1.0) + b / (m + 1.0);
            const double int_s2 = a * a / (2.0 * n + 1.0) + b * b / (2.0 * m + 1.0) + 2.0 * a * b / (n + m + 1.0);
            worst = std::max({worst, std::abs(s.delta_star - (1.0 - int_s)), std::abs(s.theta - (int_s - int_s2))});
        }
    CHECK(sref::measure("endwall 8-40, 8-41 vs the streamwise integrals", worst) <= 1e-14);  // measured 1.7e-16
    // zero skew (m = n): theta = theta11, delta* = delta1* at any angle
    const EndwallProfile c = endwall_profile(1.0, 0.25, 0.25);
    const EndwallStreamwise z = endwall_streamwise(c.delta1, c.delta2, c.theta11, c.theta12, c.theta22, 0.9);
    CHECK(rel(z.theta, c.theta11) < 1e-14);
    CHECK(rel(z.delta_star, c.delta1) < 1e-14);
    // the >= 0 limits of p. 186
    const EndwallStreamwise neg = endwall_streamwise(0.1, -0.05, 0.05, -0.02, -0.01, 0.7);
    const EndwallStreamwise cut = endwall_streamwise(0.1, 0.0, 0.05, 0.0, 0.0, 0.7);
    CHECK(neg.theta == cut.theta);
    CHECK(neg.delta_star == cut.delta_star);
}

TEST_CASE("end wall: c_f and E are the Sec. 5.6 correlations", "[turbolenza][endwall]")
{
    double worst = 0.0;
    for (double H : {1.1, 1.3, 1.6, 2.4})
        for (double re : {250.0, 3e3, 1e5}) worst = std::max(worst, rel(ludwieg_tillmann_cf(H, re), green_cf(H, H, re)));
    CHECK(sref::measure("endwall 8-36 vs BoundaryLayer 5-134/5-135 at H = H_k", worst) <= 1e-14);  // measured 2.1e-16
    CHECK(std::abs(ludwieg_tillmann_cf(1.4, 1e4) - 0.246 * std::exp(-1.561 * 1.4) * std::pow(1e4, -0.268)) < 1e-18);
    CHECK(head_entrainment(1.4) == 0.025 * (1.4 - 1.0));  // 8-42 / 8-43 have 5-133's form
}

TEST_CASE("end wall: a straight annulus without swirl is the turbulent flat plate", "[turbolenza][endwall]")
{
    // Re_L = 1.2 * 50 * 2 / 1.8e-5 = 6.7e6 (the BoundaryLayer flat-plate case)
    const double rho = 1.2, mu = 1.8e-5, V = 50.0, L = 2.0;
    EndwallInput in = annulus(401, L, 0.5, 0.6, V, 0.0, rho, mu);
    in.inlet.blockage = 1e-9;  // p. 194: a very small B starts the layer at the transition point (8-73)
    const EndwallResult r = analyze_endwalls(in);
    const EndwallLayer& h0 = r.stations.front().hub;
    CHECK(h0.re_theta11_limited);
    CHECK(rel(h0.Re_theta11, 250.0) < 1e-14);
    const EndwallLayer& h = r.stations.back().hub;
    const EndwallLayer& s = r.stations.back().shroud;
    CHECK(rel(h.theta11, s.theta11) < 1e-13);  // r cancels on a straight annulus
    CHECK(h.q12 == 0.0);
    CHECK(rel(h.H, h.H1) < 1e-14);   // no swirl: the streamwise layer is the meridional one
    CHECK(rel(h.theta, h.theta11) < 1e-14);
    CHECK(!h.H1_limited);
    // the 1/7-power law theta = 0.036 x Re_x^-0.2 (from x = 0)
    const double re = rho * V * L / mu;
    const double theta_power = 0.036 * L * std::pow(re, -0.2);
    const double e_power = sref::measure("endwall flat plate theta11(L) / 1/7-law - 1", h.theta11 / theta_power - 1.0);
    CHECK(std::abs(e_power) <= 0.1);  // measured -5.3e-2: the power-law closure against the power law
    sref::measure("endwall flat plate H1(L)", h.H1);
    CHECK(h.H1 > 1.2);
    CHECK(h.H1 < 1.45);
    sref::measure("endwall flat plate c_f(L) / 0.0576 Re^-0.2 - 1", h.cf / (0.0576 * std::pow(re, -0.2)) - 1.0);
    CHECK(std::abs(h.cf / (0.0576 * std::pow(re, -0.2)) - 1.0) < 0.15);
    // the Sec. 5.6 march (Head / Green, laminar start) on the same edge
    EdgeConditions e;
    for (const EndwallStation& st : in.stations) e.x.push_back(st.hub.m);
    e.W.assign(e.x.size(), V);
    e.rho.assign(e.x.size(), rho);
    e.T.assign(e.x.size(), 300.0);
    e.mu.assign(e.x.size(), mu);
    e.cp = 1e12;
    e.total_temperature = 300.0;
    const BoundaryLayerProfile bl = march_surface(e, {});
    const double e_bl = sref::measure("endwall flat plate theta11(L) / Sec. 5.6 theta(L) - 1", h.theta11 / bl.theta.back() - 1.0);
    CHECK(std::abs(e_bl) <= 0.01);  // measured -8.4e-4
    sref::measure("endwall flat plate delta1*(L) / Sec. 5.6 delta*(L) - 1", h.delta1 / bl.delta_star.back() - 1.0);
}

TEST_CASE("end wall: the integration is second order", "[turbolenza][endwall]")
{
    // a diverging annulus, decelerating free-vortex swirl, stationary walls, no blades
    auto run = [](int N) {
        EndwallInput in;
        const double L = 0.5, slope = 0.05;
        for (int k = 0; k < N; ++k) {
            const double m = L * k / (N - 1);
            EndwallStation s;
            const double rh = 0.30 + slope * m, rs = 0.50 + slope * m;
            const double Vm = 150.0 - 40.0 * m, rho = 1.2 + 0.1 * m;
            s.hub = wall(m, rh, Vm, 24.0 / rh, rho, 1e5, 1.8e-5, 0.0, std::asin(slope));
            s.shroud = wall(m, rs, Vm, 24.0 / rs, rho, 1e5, 1.8e-5, 0.0, std::asin(slope));
            in.stations.push_back(s);
        }
        in.mass_flow = 1.2 * 150.0 * kPi * (0.25 - 0.09);
        return analyze_endwalls(in);
    };
    const EndwallResult ref = run(641);
    double e[3];
    const int Ns[3] = {41, 81, 161};
    for (int i = 0; i < 3; ++i) {
        const EndwallResult r = run(Ns[i]);
        e[i] = std::abs(r.stations.back().hub.theta11 - ref.stations.back().hub.theta11) / ref.stations.back().hub.theta11;
        CHECK(!r.stations.back().hub.H1_limited);
    }
    const double p1 = sref::measure("endwall order of accuracy (41 -> 81)", std::log2(e[0] / e[1]));
    const double p2 = sref::measure("endwall order of accuracy (81 -> 161)", std::log2(e[1] / e[2]));
    CHECK(std::abs(p1 - 2.0) < 0.2);
    CHECK(std::abs(p2 - 2.0) < 0.2);
    // K_B against its definition [8-76, 8-77]; station 0 carries the inlet B = 0.02 [8-72] (Fig. 8-9 starts there)
    double worst = 0.0;
    for (std::size_t k = 0; k < ref.stations.size(); k += 40) {
        const EndwallStationResult& s = ref.stations[k];
        worst = std::max(worst, std::abs(s.blockage.K_B - (1.0 - s.blockage.B_hub - s.blockage.B_shroud)));
    }
    CHECK(worst <= 1e-15);
    const EndwallStationResult& s0 = ref.stations.front();
    CHECK(!s0.hub.re_theta11_limited);
    CHECK(std::abs(s0.blockage.B - 0.02) < 1e-15);
    CHECK(std::abs(s0.blockage.K_B - 0.98) < 1e-15);
    // each wall's part is 2 pi r rho_e V_me delta1* / mdot, recomputed from the reported delta1*
    const double mdot = 1.2 * 150.0 * kPi * (0.25 - 0.09);
    const EndwallStationResult& sl = ref.stations.back();
    const double Bh = 2.0 * kPi * 0.325 * 1.25 * 130.0 * sl.hub.delta1 / mdot;
    const double Bs = 2.0 * kPi * 0.525 * 1.25 * 130.0 * sl.shroud.delta1 / mdot;
    CHECK(rel(sl.blockage.B, Bh + Bs) < 1e-14);
    CHECK(std::abs(s0.blockage.B_hub - 0.01) < 1e-15);  // B shared equally [8-72]
    sref::measure("endwall diverging annulus K_B at exit", sl.blockage.K_B);
}

TEST_CASE("end wall: zero skew is the two-dimensional layer along the edge streamline", "[turbolenza][endwall]")
{
    // beta_e = 45 deg, uniform, stationary walls; inlet m = n (collateral)
    const double rho = 1.2, mu = 1.8e-5, Vm = 100.0, Ct = 100.0, L = 0.6;
    EndwallInput sw = annulus(121, L, 0.4, 0.5, Vm, Ct, rho, mu);
    sw.inlet.m = 0.5 * (sw.inlet.H1 - 1.0);
    const EndwallResult a = analyze_endwalls(sw);
    // the same layer without swirl at V_e, along s = m / cos(beta) (8-72 then needs mdot V_e / V_m for the same delta1*)
    const double Ve = std::hypot(Vm, Ct), cb = Vm / Ve;
    EndwallInput two = annulus(121, L / cb, 0.4, 0.5, Ve, 0.0, rho, mu);
    two.mass_flow = sw.mass_flow * Ve / Vm;
    two.inlet = sw.inlet;
    const EndwallResult b = analyze_endwalls(two);
    double worst_collateral = 0.0, worst_2d = 0.0;
    for (std::size_t k = 0; k < a.stations.size(); ++k) {
        const EndwallLayer &x = a.stations[k].hub, &y = b.stations[k].hub;
        worst_collateral = std::max({worst_collateral, rel(x.theta12, x.theta11), rel(x.m, x.n)});
        worst_2d = std::max({worst_2d, rel(x.theta11, y.theta11), rel(x.Q, y.Q), rel(x.delta_star, y.delta1),
                             rel(x.theta, y.theta11)});
    }
    CHECK(sref::measure("endwall zero skew: theta12 / theta11 - 1, m / n - 1", worst_collateral) <= 1e-12);  // measured 1.8e-15
    CHECK(sref::measure("endwall zero skew vs the 2-D layer along s", worst_2d) <= 1e-12);  // measured 1.7e-15
    // the collateral condition [8-45] holds for a stationary "row" frame
    CHECK(std::abs(a.stations.back().hub.q12 - Ct * a.stations.back().hub.theta11) <= 1e-12 * Ct * a.stations.back().hub.theta11);
}

TEST_CASE("end wall: the jump conditions 8-67 .. 8-71", "[turbolenza][endwall]")
{
    // the book's algebra: a stationary-wall power-law layer C_theta = C_theta_e eta^m seen from a wall rotating at omega
    {
        const double n = 0.2, m = 0.3, Cte = 80.0, wr = 120.0;  // omega r
        const double Wte = Cte - wr;                              // [8-67]
        const double t12 = 1.0 / (n + 1.0) - 1.0 / (n + m + 1.0), d2 = 1.0 - 1.0 / (m + 1.0);
        const double t22 = 1.0 / (m + 1.0) - 1.0 / (2.0 * m + 1.0);
        // the primed thicknesses from the definitions 8-13 .. 8-15 with W_theta = C_theta - omega r
        const double t12p = Cte * t12 / Wte;
        const double d2p = (Wte - (Cte / (m + 1.0) - wr)) / Wte;
        const double t22p = ((Cte / (m + 1.0) - wr) * Wte - (Cte * Cte / (2.0 * m + 1.0) - 2.0 * wr * Cte / (m + 1.0) + wr * wr)) /
                            (Wte * Wte);
        CHECK(std::abs(Cte * t12 - Wte * t12p) < 1e-13);                                               // [8-68]
        CHECK(std::abs(Cte * d2 - Wte * d2p) < 1e-13);                                                 // [8-69]
        CHECK(std::abs(Cte * Cte * t22 - (Wte * Wte * t22p + wr * Wte * d2p)) < 1e-11);                // [8-70]
        CHECK(std::abs(Cte * Cte * (t22 + d2) - (Wte * Wte * (t22p + d2p) + 2.0 * wr * Wte * d2p)) < 1e-11);  // [8-71]
    }
    // the march across a stationary -> rotating hub at station 5 (omega = 400 rad/s, r = 0.3: omega r = 120)
    const double rho = 1.2, mu = 1.8e-5;
    EndwallInput in = annulus(21, 0.3, 0.3, 0.45, 100.0, 80.0, rho, mu);
    for (std::size_t k = 6; k < in.stations.size(); ++k) in.stations[k].hub.wall_omega = 400.0;
    const EndwallResult r = analyze_endwalls(in);
    const EndwallLayer& c = r.stations[5].hub;  // stationary frame
    const EndwallLayer w = endwall_closures(state_of(c), in.stations[5].hub, 400.0, c.m, in.settings);  // rotating frame
    // the march carries delta - delta1*, theta11 and V_theta_e theta12 (8-68) unchanged; the meridional layer is frame-free
    CHECK(w.Q == c.Q);
    CHECK(w.theta11 == c.theta11);
    CHECK(w.delta1 == c.delta1);
    CHECK(w.q12 == c.q12);
    CHECK(rel(c.Vtheta_e * c.theta12, w.Vtheta_e * w.theta12) < 1e-14);  // [8-68]
    CHECK(w.Vtheta_e == Catch::Approx(80.0 - 0.3 * 400.0).margin(1e-12));  // [8-4]; not bitwise: fused multiply-add
    CHECK(w.theta12 < 0.0);  // W_theta and C_theta of opposite sign: theta12' < 0 (p. 183)
    CHECK(w.m < 0.05);       // the extended profile 8-28
    // the profile family does not carry 8-69 and 8-71 exactly (measured, not imposed)
    const double e69 = (w.Vtheta_e * w.delta2) / (c.Vtheta_e * c.delta2) - 1.0;
    const double lhs71 = c.Vtheta_e * c.Vtheta_e * c.theta22 * (1.0 + c.H2);
    const double rhs71 = w.Vtheta_e * w.Vtheta_e * w.theta22 * (1.0 + w.H2) + 2.0 * 120.0 * w.Vtheta_e * w.delta2;
    sref::measure("endwall jump: 8-69 through the profile family, ratio - 1", e69);
    sref::measure("endwall jump: 8-71 through the profile family, ratio - 1", rhs71 / lhs71 - 1.0);
    // downstream of the jump the station's layer is in the rotating frame
    CHECK(r.stations[6].hub.wall_omega == 400.0);
    CHECK(r.stations[5].hub.wall_omega == 0.0);
    CHECK(std::isfinite(r.stations.back().blockage.K_B));
    sref::measure("endwall jump: hub theta11 at exit / at the jump", r.stations.back().hub.theta11 / c.theta11);
}

namespace {

// a shrouded stator (hub) between stations 2 and 4, turning C_theta 150 -> 60, mdot_leak at the hub
EndwallInput shrouded_stator(double mdot_leak)
{
    EndwallInput in;
    const int N = 7;
    for (int k = 0; k < N; ++k) {
        const double m = 0.03 * k;
        double Vm = 120.0, Ct = 150.0, P = 1.00e5;
        if (k >= 2) {
            const double t = std::min(1.0, (k - 2) / 2.0);
            Vm = 120.0 - 10.0 * t;
            Ct = 150.0 - 90.0 * t;
            P = 1.00e5 + 6.0e3 * t;
        }
        EndwallStation s;
        s.hub = wall(m, 0.30, Vm, Ct, 1.2, P, 1.8e-5);
        s.shroud = wall(m, 0.45, Vm, Ct * 0.30 / 0.45, 1.2, P, 1.8e-5);
        s.quasi_normal_length = 0.15;
        in.stations.push_back(s);
    }
    in.mass_flow = 1.2 * 120.0 * kPi * (0.45 * 0.45 - 0.09);
    in.compressor_omega = 1000.0;
    EndwallBladeRow row;
    row.station_in = 2;
    row.station_out = 4;
    row.hub = {0.0, 40, kPi * 0.6 / 40.0, 0.5, 0.0, true, mdot_leak};
    row.shroud = {0.0, 40, kPi * 0.9 / 40.0, 0.5, 5e-4, false, 0.0};
    in.rows.push_back(row);
    return in;
}

// rho_e V_me times: (delta - delta1*) (mass), V_me (delta - delta1* - theta11) (meridional momentum),
// V_theta_e (delta - delta1*) - V_theta_e theta12 (tangential momentum relative to the wall)
struct Fluxes
{
    double mass, mom_m, mom_t;
};
Fluxes fluxes(const EndwallState& s, const EndwallWallStation& w, double Vt)
{
    const double k = 2.0 * kPi * w.r * w.rho * w.Vm;
    return {k * s.Q, k * w.Vm * (s.Q - s.theta11), k * (Vt * s.Q - s.q12)};
}

}  // namespace

TEST_CASE("end wall: the seal leakage corrections conserve mass and momentum (Exercises 8.1 - 8.4)", "[turbolenza][endwall]")
{
    for (double leak : {0.05, -0.05}) {
        const EndwallInput in = shrouded_stator(leak);
        const EndwallResult r = analyze_endwalls(in);
        const EndwallRowWallResult& h = r.rows[0].hub;
        REQUIRE(h.shrouded);
        const EndwallWallStation &a = in.stations[2].hub, &b = in.stations[4].hub;
        // U_leak = r omega_c / 2 relative to a stationary shroud [8-61]
        CHECK(std::abs(h.leakage.U_leak_in - 0.5 * 0.3 * 1000.0) < 1e-12);
        CHECK(std::abs(h.leakage.U_leak_out - 0.5 * 0.3 * 1000.0) < 1e-12);
        const Fluxes fi = fluxes(h.leakage.in, a, a.Ctheta), fp = fluxes(h.leakage.plus, a, a.Ctheta);
        const Fluxes fm = fluxes(h.leakage.minus, b, b.Ctheta), fo = fluxes(h.leakage.out, b, b.Ctheta);
        // inlet: + mdot_leak of mass, no meridional momentum, tangential momentum mdot U_leak when it enters
        CHECK(std::abs(fp.mass - fi.mass - leak) < 1e-14);
        CHECK(std::abs(fp.mom_m - fi.mom_m) < 1e-12);
        CHECK(std::abs(fp.mom_t - fi.mom_t - (leak > 0.0 ? leak * h.leakage.U_leak_in : 0.0)) < 1e-11);
        // discharge: - mdot_leak of mass; the fluid leaves at rest relative to the wall, or enters with U_leak
        CHECK(std::abs(fo.mass - fm.mass + leak) < 1e-14);
        CHECK(std::abs(fo.mom_m - fm.mom_m) < 1e-12);
        CHECK(std::abs(fo.mom_t - fm.mom_t - (leak < 0.0 ? -leak * h.leakage.U_leak_out : 0.0)) < 1e-11);
        // the station records keep the layer outside the row
        CHECK(r.stations[2].hub.Q == h.leakage.in.Q);
        CHECK(r.stations[4].hub.Q == h.leakage.out.Q);
        CHECK(r.stations[4].hub.q12 == h.leakage.out.q12);
        CHECK(!r.rows[0].shroud.shrouded);
    }
    // a shrouded row's wall rotates with the blades (p. 193)
    EndwallInput bad = shrouded_stator(0.05);
    bad.stations[3].hub.wall_omega = 100.0;
    CHECK_THROWS_AS(analyze_endwalls(bad), EndwallError);
    // leakage needs a shroud
    bad = shrouded_stator(0.05);
    bad.rows[0].shroud.seal_leakage = 0.01;
    CHECK_THROWS_AS(analyze_endwalls(bad), EndwallError);
}

namespace {

// a rotor (omega 800 rad/s) between stations 2 and 4: hub drum rotating with it, casing stationary with a clearance
EndwallInput rotor(EndwallForceDefect model, bool correction, double inlet_blockage = 0.02)
{
    EndwallInput in;
    in.inlet.blockage = inlet_blockage;
    const int N = 8;
    const double omega = 800.0;
    for (int k = 0; k < N; ++k) {
        const double m = 0.025 * k;
        const double t = k <= 2 ? 0.0 : std::min(1.0, (k - 2) / 2.0);
        const double Vm = 150.0 - 15.0 * t, P = 1.0e5 + 9.0e3 * t, rho = 1.2 + 0.08 * t;
        EndwallStation s;
        s.hub = wall(m, 0.30, Vm, 110.0 * t, rho, P, 1.8e-5, (k >= 3 && k <= 4) ? omega : 0.0);
        s.shroud = wall(m, 0.45, Vm, 0.30 / 0.45 * 110.0 * t + 20.0 * t, rho, P, 1.8e-5);
        s.quasi_normal_length = 0.15;
        in.stations.push_back(s);
    }
    in.mass_flow = 1.2 * 150.0 * kPi * (0.45 * 0.45 - 0.09);
    in.compressor_omega = omega;
    EndwallBladeRow row;
    row.station_in = 2;
    row.station_out = 4;
    row.hub = {omega, 30, kPi * 0.6 / 30.0, 0.3, 0.0, false, 0.0};
    row.shroud = {omega, 30, kPi * 0.9 / 30.0, 0.9, 1e-3, false, 0.0};
    in.rows.push_back(row);
    in.settings.force_defect = model;
    in.settings.force_defect_correction = correction;
    return in;
}

}  // namespace

TEST_CASE("end wall: the blade force defect of Sec. 8.5", "[turbolenza][endwall]")
{
    // the closures: Fig. 8-6's curve x^3 / (1 + x^3), x = 8 delta*/g, read at delta*/g = 0.10, 0.20, 0.29 as
    // 0.35, 0.80, 0.925 (reading uncertainty +-0.02)
    CHECK(std::abs(force_defect_shape(0.8) - 0.35) < 0.02);
    CHECK(std::abs(force_defect_shape(1.6) - 0.80) < 0.02);
    CHECK(std::abs(force_defect_shape(2.32) - 0.925) < 0.02);
    CHECK(force_defect_shape(1.0) == 0.5);
    CHECK(std::abs(force_defect_factor(0.02, 10.0) - 0.5) < 1e-9);                      // F -> 1/2 [8-55]
    CHECK(std::abs(base_force_defect(0.02, 1e-3, 10.0) - (0.12 * 0.02 + 5e-4)) < 1e-12);  // nu0 -> 0.12 g + delta_c/2
    CHECK(rel(staggered_spacing(0.05, std::acos(0.6)), 0.03) < 1e-15);

    // a thin (B 0.02 at the inlet) and a thick (0.15) layer: 8 theta11 / g from about 0.1 to about 1
    for (double Bin : {0.02, 0.15}) {
        const std::string tag = Bin < 0.1 ? " (B 0.02)" : " (B 0.15)";
        const EndwallInput in = rotor(EndwallForceDefect::RowChapter8, true, Bin);
        const EndwallResult r = analyze_endwalls(in);
        for (bool hub : {true, false}) {
            const EndwallRowWallResult& w = hub ? r.rows[0].hub : r.rows[0].shroud;
            const EndwallRowWall& rw = hub ? in.rows[0].hub : in.rows[0].shroud;
            const EndwallWallStation& b = hub ? in.stations[4].hub : in.stations[4].shroud;
            const EndwallLayer& L_in = hub ? r.stations[2].hub : r.stations[2].shroud;
            const EndwallLayer& L_out = hub ? r.stations[4].hub : r.stations[4].shroud;
            // g [8-2] of the given s bar and gamma bar
            CHECK(rel(w.g, rw.pitch * std::cos(rw.stagger)) < 1e-15);
            CHECK(rel(w.theta11_mean, 0.5 * (L_in.theta11 + L_out.theta11)) < 1e-15);
            const double x = 8.0 * w.theta11_mean / w.g, x3 = x * x * x;
            CHECK(rel(w.nu0, (0.12 * w.g + 0.5 * rw.clearance) * x3 / (1.0 + x3)) < 1e-14);  // [8-44]
            CHECK(rel(w.F, 0.5 * x3 / (1.0 + x3)) < 1e-14);                                   // [8-55]
            // the collateral defect [8-45] at the discharge, W relative to the rotor
            CHECK(rel(w.collateral_defect, L_out.q12 - (b.Ctheta - b.r * 800.0) * L_out.theta11) < 1e-14);
            // converged: nu = nu0 + dnu [8-56] and dnu is the fixed point of 8-57
            CHECK(std::abs(w.nu - (w.nu0 + w.dnu)) <= 1e-11 * w.g);
            CHECK(std::abs(w.dnu - (w.dnu - w.collateral_defect / w.xi) * w.F) <= 1e-9 * std::abs(w.dnu) + 1e-12 * w.g);
            sref::measure(std::string(hub ? "endwall rotor hub nu / delta1*_out" : "endwall rotor casing nu / delta1*_out") + tag,
                          w.nu / L_out.delta1);
            sref::measure(std::string(hub ? "endwall rotor hub row iterations" : "endwall rotor casing row iterations") + tag, w.iterations);
        }
        // 8-54 with the boundary layer obtained with dnu = 0 (nu = nu0): 8-57's fixed point up to f_me of 8-8 vs 8-49
        const EndwallResult r0 = analyze_endwalls(rotor(EndwallForceDefect::RowChapter8, false, Bin));
        for (bool hub : {true, false}) {
            const EndwallRowWallResult& w = hub ? r.rows[0].hub : r.rows[0].shroud;
            const EndwallRowWallResult& w0 = hub ? r0.rows[0].hub : r0.rows[0].shroud;
            CHECK(w0.dnu == 0.0);
            CHECK(std::abs(w0.nu - w0.nu0) <= 1e-11 * w0.g);
            const double dnu54 = -w0.collateral_defect * w0.F / w0.xi;
            sref::measure(std::string(hub ? "endwall rotor hub dnu(8-57) / dnu(8-54) - 1" : "endwall rotor casing dnu(8-57) / dnu(8-54) - 1") + tag,
                          w.dnu / dnu54 - 1.0);
            CHECK(std::abs(w.dnu / dnu54 - 1.0) < 0.5);
            CHECK(std::abs(w.collateral_defect) < std::abs(w0.collateral_defect));  // toward collateral flow
            sref::measure(std::string(hub ? "endwall rotor hub collateral defect corrected / uncorrected"
                              : "endwall rotor casing collateral defect corrected / uncorrected") + tag,
                          w.collateral_defect / w0.collateral_defect);
        }
        // no force defect: nu = 0
        const EndwallResult rn = analyze_endwalls(rotor(EndwallForceDefect::None, true, Bin));
        CHECK(rn.rows[0].shroud.nu == 0.0);
        CHECK(rn.rows[0].shroud.nu0 == 0.0);
        sref::measure(std::string("endwall rotor K_B at exit") + tag, r.stations.back().blockage.K_B);
        sref::measure(std::string("endwall rotor K_B at exit without nu") + tag, rn.stations.back().blockage.K_B);
        CHECK(r.stations.back().blockage.K_B != rn.stations.back().blockage.K_B);
    }
}

TEST_CASE("end wall: refusals", "[turbolenza][endwall]")
{
    const EndwallInput good = annulus(5, 0.1, 0.3, 0.4, 100.0, 50.0, 1.2, 1.8e-5);
    CHECK_NOTHROW(analyze_endwalls(good));
    EndwallInput e = good;
    e.stations.resize(1);
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);  // one station
    e = good;
    e.stations[2].hub.m = e.stations[1].hub.m;
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);  // m not increasing
    e = good;
    e.stations[3].shroud.Vm = 0.0;
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);  // V_m = 0
    e = good;
    e.stations[3].hub.mu = std::nan("");
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);
    e = good;
    e.mass_flow = 0.0;
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);
    e = good;
    e.inlet.blockage = 1.0;
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);
    e = good;
    e.inlet.H1 = 1.0;
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);
    e = good;
    EndwallBladeRow row;
    row.station_in = 1;
    row.station_out = 3;
    row.hub = row.shroud = {0.0, 20, 0.1, 0.4, 0.0, false, 0.0};
    e.rows = {row, row};
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);  // overlapping rows
    e.rows = {row};
    e.rows[0].hub.pitch = 0.0;
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);
    e.rows = {row};
    e.rows[0].station_out = 5;
    CHECK_THROWS_AS(analyze_endwalls(e), EndwallError);  // beyond the last station
    e.rows = {row};
    CHECK_NOTHROW(analyze_endwalls(e));
    CHECK_THROWS_AS(endwall_closures({0.0, 1e-3, 0.0}, good.stations[0].hub, 0.0, 0.1, {}), EndwallError);
    CHECK_THROWS_AS(endwall_inlet_state(1e-3, 1.0, 0.1, 0.0), EndwallError);
    CHECK_THROWS_AS(endwall_profile(1e-3, 0.0, 0.1), EndwallError);
}

TEST_CASE("end wall: the step is pure and the march is the sequence of steps", "[turbolenza][endwall]")
{
    const EndwallInput in = rotor(EndwallForceDefect::RowChapter8, true);
    const EndwallResult r = analyze_endwalls(in);
    for (bool hub : {true, false}) {
        auto W = [&](std::size_t k) -> const EndwallWallStation& { return hub ? in.stations[k].hub : in.stations[k].shroud; };
        auto L = [&](std::size_t k) -> const EndwallLayer& { return hub ? r.stations[k].hub : r.stations[k].shroud; };
        EndwallLayer prev = endwall_inlet(W(0), in.inlet, in.mass_flow, in.settings);
        CHECK(prev.theta11 == L(0).theta11);
        for (std::size_t k = 1; k <= 2; ++k) {  // the blade-free steps before the row 2 -> 4
            const EndwallStepResult a = endwall_step(prev, W(k - 1), W(k), nullptr, in.compressor_omega, in.settings);
            const EndwallStepResult b = endwall_step(prev, W(k - 1), W(k), nullptr, in.compressor_omega, in.settings);
            CHECK(a.layer.theta11 == b.layer.theta11);  // pure: repeatable bit for bit
            CHECK(a.layer.q12 == b.layer.q12);
            CHECK(a.layer.theta11 == L(k).theta11);     // the march's station
            CHECK(a.layer.Q == L(k).Q);
            CHECK(!a.row);
            prev = a.layer;
        }
    }
    // a one-step row (the Ch. 7 case): the row's march over 2 -> 4 without station 3 equals endwall_step with the row
    EndwallInput one = in;
    one.stations.erase(one.stations.begin() + 3);
    one.rows[0].station_out = 3;
    const EndwallResult r1 = analyze_endwalls(one);
    const EndwallStepResult s = endwall_step(r1.stations[2].shroud, one.stations[2].shroud, one.stations[3].shroud,
                                             &one.rows[0].shroud, one.compressor_omega, one.settings);
    REQUIRE(s.row);
    CHECK(s.layer.theta11 == r1.stations[3].shroud.theta11);
    CHECK(s.layer.q12 == r1.stations[3].shroud.q12);
    CHECK(s.row_result.nu == r1.rows[0].shroud.nu);
    CHECK(s.layer.nu == s.row_result.nu);
    const EndwallStepResult s2 = endwall_step(r1.stations[2].shroud, one.stations[2].shroud, one.stations[3].shroud,
                                              &one.rows[0].shroud, one.compressor_omega, one.settings);
    CHECK(s2.layer.theta11 == s.layer.theta11);
    // one interval per row against two: the same row, measured
    sref::measure("endwall rotor casing theta11 at TE, one step / two steps - 1",
                  r1.stations[3].shroud.theta11 / r.stations[4].shroud.theta11 - 1.0);
    // endwall_blockage is the station's blockage
    const EndwallBlockage b = endwall_blockage(r1.stations[3].hub, one.stations[3].hub, r1.stations[3].shroud,
                                               one.stations[3].shroud, one.mass_flow);
    CHECK(b.K_B == r1.stations[3].blockage.K_B);
    CHECK_THROWS_AS(endwall_blockage(r1.stations[3].hub, one.stations[3].hub, r1.stations[3].shroud,
                                     one.stations[3].shroud, 0.0),
                    EndwallError);
    // a step backwards in m is refused
    CHECK_THROWS_AS(endwall_step(r1.stations[3].hub, one.stations[3].hub, one.stations[2].hub, nullptr, 0.0, one.settings),
                    EndwallError);
}
