// The blade-surface boundary layer of Aungier (2003) Sec. 5.6 (BoundaryLayer.hpp): Sutherland's law,
// the Gruschwitz closures against the profile's own integrals, the flat plate against Pohlhausen's
// closed form, the plane stagnation point against the method's own equilibrium, transition, the
// turbulent flat plate against the 1/7-power law, the H_k limit, the wake, the loss, the refusals.
#include "SolverRef.hpp"

#include "turbolenza/BoundaryLayer.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

using namespace BladeGen::Turbolenza;

namespace {

// 8-point Gauss-Legendre on [0, 1]: exact for polynomials up to degree 15
double integrate01(const auto& f)
{
    constexpr double x[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    constexpr double w[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    double s = 0.0;
    for (int k = 0; k < 4; ++k) s += w[k] * (f(0.5 - 0.5 * x[k]) + f(0.5 + 0.5 * x[k]));
    return 0.5 * s;
}

// the Gruschwitz / Pohlhausen profile W/W_e = C1 eta + C2 eta^2 + C3 eta^3 + C4 eta^4 (Eq. 5-114)
double profile(double L, double e)
{
    return (2.0 + L / 6.0) * e + (-L / 2.0) * e * e + (L / 2.0 - 2.0) * e * e * e + (1.0 - L / 6.0) * e * e * e * e;
}

// an incompressible edge: W_e^2 / (2 cp T) ~ 1e-10 (cp huge), T't = T (b0 = 1)
EdgeConditions incompressible(const std::vector<double>& x, const std::vector<double>& W, double rho, double mu)
{
    EdgeConditions e;
    e.x = x;
    e.W = W;
    e.rho.assign(x.size(), rho);
    e.T.assign(x.size(), 300.0);
    e.mu.assign(x.size(), mu);
    e.cp = 1e12;
    e.total_temperature = 300.0;
    return e;
}

std::vector<double> uniform(double a, double b, int n)
{
    std::vector<double> x(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) x[static_cast<std::size_t>(k)] = a + (b - a) * k / (n - 1);
    return x;
}

}  // namespace

TEST_CASE("boundary layer: Sutherland's law", "[turbolenza][bl]")
{
    const SutherlandLaw air;
    CHECK(sutherland_viscosity(273.15, air) == 1.716e-5);  // T = T_ref: mu_ref exactly
    // 300 K: 1.716e-5 (300/273.15)^1.5 383.55/410.4 = 1.8459e-5 Pa s (air tables: 1.846e-5)
    CHECK(std::abs(sutherland_viscosity(300.0, air) / 1.8459e-5 - 1.0) < 2e-5);
    CHECK(sutherland_viscosity(400.0, air) > sutherland_viscosity(300.0, air));  // a gas: mu rises with T
    CHECK_THROWS_AS(sutherland_viscosity(0.0, air), BoundaryLayerError);
    CHECK_THROWS_AS(sutherland_viscosity(std::nan(""), air), BoundaryLayerError);
    CHECK_THROWS_AS(sutherland_viscosity(300.0, SutherlandLaw{-1.0, 273.15, 110.4}), BoundaryLayerError);
}

TEST_CASE("boundary layer: the Gruschwitz closures are the profile's integrals", "[turbolenza][bl]")
{
    double worst_exact = 0.0, worst_F = 0.0, worst_cubic = 0.0;
    for (double L : {-12.0, -6.0, -1.5, 0.0, 2.0, 7.052, 12.0}) {
        // the edge conditions of the quartic: f(1) = 1, f'(1) = 0, f''(1) = 0, f''(0) = -Lambda
        CHECK(std::abs(profile(L, 1.0) - 1.0) < 1e-15);
        // theta / delta' = int f (1 - f) [5-115], delta_E / delta' = int f (1 - f^2) [5-116],
        // int (1 - f) = 3/10 - Lambda/120 [5-117 at M = 0], int f = 7/10 + Lambda/120 [5-138]
        const double th = integrate01([&](double e) { const double f = profile(L, e); return f * (1.0 - f); });
        const double en = integrate01([&](double e) { const double f = profile(L, e); return f * (1.0 - f * f); });
        const double vel = integrate01([&](double e) { return 1.0 - profile(L, e); });
        const double mass = integrate01([&](double e) { return profile(L, e); });
        const LaminarThicknesses t0 = laminar_thicknesses(L, 1.0, 0.0);  // M = 0: delta_h = 0, F unused
        worst_exact = std::max({worst_exact, std::abs(th - gruschwitz_theta_ratio(L)),
                                std::abs(vel - t0.H * gruschwitz_theta_ratio(L)),
                                std::abs(mass - t0.delta_minus_delta_star * gruschwitz_theta_ratio(L))});
        // delta_E: [5-116] read through delta_h = m2 delta_E: H(m2) - H(0) = m2 (E + F) / g
        const double m2 = 1e-3;
        const LaminarThicknesses t1 = laminar_thicknesses(L, 1.0, m2);
        const double e_plus_f = (t1.H - t0.H) * gruschwitz_theta_ratio(L) / m2;
        // F [5-118]: with Pr = 1 and an adiabatic wall the Crocco-Busemann relation T/T_e = 1 + m2 (1 - f^2)
        // gives int (1 - W/W_e) dy = delta' [int (1 - f) + m2 int (1 - f)(1 - f^2)], so
        // F = int (1 - f)(1 - f^2) deta, a cubic in Lambda (our derivation, not Aungier's):
        //     20983/90090 - 8989/1081080 L + 211/3243240 L^2 + 1/617760 L^3
        //   = 0.2329115 - 0.8314833 (L/100) + 0.6505840 (L/100)^2 + 1.6187516 (L/100)^3.
        // Aungier prints 0.232912, 0.831483, 0.650584 (these agree) and 17.8063 for the cubic term, in both
        // books (2003 Eq. 5-118, 2000 Eq. 15-17).  The code keeps the printed value; the test pins the
        // difference to that one coefficient.
        const double F = integrate01([&](double e) { const double f = profile(L, e); return (1.0 - f) * (1.0 - f * f); });
        const double l = L / 100.0;
        worst_F = std::max(worst_F, std::abs(e_plus_f - (en + F)));
        worst_cubic = std::max(worst_cubic, std::abs(e_plus_f - (en + F) - (17.8063 - 1.0 / 617760.0 * 1e6) * l * l * l));
    }
    CHECK(sref::measure("bl Gruschwitz polynomials vs the profile integrals", worst_exact) <= 1e-14);  // measured 4.4e-16
    // printed F against the integral on [-12, 12]: the cubic coefficient only (17.8063 against 1.61875)
    sref::measure("bl Gruschwitz F as printed vs int (1-f)(1-f^2), |L| <= 12", worst_F);  // measured 2.8e-2 at L = -12
    CHECK(sref::measure("bl Gruschwitz F as printed minus its cubic-term difference", worst_cubic) <= 2e-6);
    // the wall slope C1 = f'(0) = 2 + Lambda/6: the Centrifugal (2000) form of c_f
    CHECK(laminar_wall_factor(3.0, LaminarCfForm::Centrifugal2000) == 2.5);
    CHECK(laminar_wall_factor(3.0, LaminarCfForm::Axial2003AsPrinted) == 1.5);
}

TEST_CASE("boundary layer: Lambda of K is the unique root on [-12, 12]", "[turbolenza][bl]")
{
    // K' = g (37/315 - 3 L/945 - 5 L^2/9072) vanishes at L = 12: K(12) is the largest K
    CHECK(gruschwitz_k(12.0) > gruschwitz_k(12.0 - 1e-3));
    CHECK((gruschwitz_k(12.0) - gruschwitz_k(12.0 - 1e-4)) / 1e-4 < 1e-6);
    double prev = gruschwitz_k(-12.0), worst = 0.0;
    for (int k = 1; k <= 240; ++k) {
        const double L = -12.0 + 0.1 * k;
        const double K = gruschwitz_k(L);
        CHECK(K > prev);
        prev = K;
        if (L < 11.9) worst = std::max(worst, std::abs(gruschwitz_lambda(K) - L));
    }
    CHECK(sref::measure("bl Lambda(K(Lambda)) - Lambda on [-12, 11.9]", worst) <= 1e-9);
    CHECK(gruschwitz_lambda(0.0) == 0.0);  // K(0) = 0 exactly: the bisection's first midpoint
    CHECK(gruschwitz_lambda(1.0) == 12.0);     // beyond K(12): Eq. 5-125
    CHECK(gruschwitz_lambda(-1.0) == -12.0);   // below K(-12): separation
    CHECK_THROWS_AS(gruschwitz_lambda(std::nan("")), BoundaryLayerError);
}

TEST_CASE("boundary layer: the laminar flat plate is Pohlhausen's", "[turbolenza][bl]")
{
    // W_e, rho, mu constant; Re_L = 1.2 * 10 * 0.1 / 1.8e-5 = 6.7e4: laminar throughout (Re_theta ~ 177)
    const double rho = 1.2, mu = 1.8e-5, W = 10.0, L = 0.1;
    const std::vector<double> x = uniform(0.0, L, 4001);
    const EdgeConditions e = incompressible(x, std::vector<double>(x.size(), W), rho, mu);
    const double nu = mu / rho;
    for (LaminarCfForm form : {LaminarCfForm::Centrifugal2000, LaminarCfForm::Axial2003AsPrinted}) {
        BoundaryLayerSettings s;
        s.laminar_cf_form = form;
        const BoundaryLayerProfile p = march_surface(e, s);
        CHECK(p.transition_index == -1);
        CHECK(p.theta.front() == 0.0);
        // d theta/dx = C1 g nu / (W theta) => theta = sqrt(2 C1 g nu x / W): C1 = 2 gives Pohlhausen's
        // 0.6854 sqrt(nu x / W) (sqrt(4 * 37/315)), C1 = 1 gives 0.4847
        const double c = std::sqrt(2.0 * laminar_wall_factor(0.0, form) * 37.0 / 315.0);
        const double exact = c * std::sqrt(nu * L / W);
        const double err = std::abs(p.theta.back() / exact - 1.0);
        const bool centrifugal = form == LaminarCfForm::Centrifugal2000;
        if (centrifugal) CHECK(std::abs(c - 0.6854) < 1e-4);
        CHECK(sref::measure(centrifugal ? "bl flat plate theta (2 + L/6) vs 0.6854 sqrt(nu x/W)"
                                        : "bl flat plate theta (1 + L/6) vs 0.4847 sqrt(nu x/W)",
                            err) <= 3e-3);  // measured 3.3e-4 (4000 steps of the first-order difference)
        // H = (3/10) / (37/315) = 2.5541 and c_f sqrt(Re_x) = 2 C1 g / c
        CHECK(std::abs(p.H.back() - 94.5 / 37.0) < 1e-8);
        CHECK(std::abs(p.Lambda.back()) < 1e-300);
        const double cf_exact = 2.0 * laminar_wall_factor(0.0, form) * 37.0 / 315.0 / c / std::sqrt(W * L / nu);
        CHECK(std::abs(p.cf.back() / cf_exact - 1.0) <= 3e-3);
        // mass deficit rho W delta*
        CHECK(std::abs(p.mass_deficit.back() - rho * W * p.delta_star.back()) <= 1e-15 * p.mass_deficit.back());
    }
}

TEST_CASE("boundary layer: the plane stagnation point reaches the method's equilibrium", "[turbolenza][bl]")
{
    // W_e = a x: theta constant with (2 + H) K = C1 g, i.e. Lambda (2 g + 3/10 - Lambda/120) = 2 + Lambda/6
    // (incompressible, C1 = 2 + Lambda/6): Lambda0 = 7.052, K0 = 0.0770, theta = sqrt(K0 nu / a)
    auto residual = [](double L) {
        return L * (2.0 * gruschwitz_theta_ratio(L) + 0.3 - L / 120.0) - (2.0 + L / 6.0);
    };
    double lo = 0.0, hi = 12.0;
    for (int k = 0; k < 200; ++k) {
        const double mid = 0.5 * (lo + hi);
        (residual(mid) < 0.0 ? lo : hi) = mid;
    }
    const double L0 = 0.5 * (lo + hi), K0 = gruschwitz_k(L0);
    CHECK(std::abs(L0 - 7.052) < 1e-3);
    CHECK(std::abs(K0 - 0.0770) < 1e-4);
    const double a = 100.0, rho = 1.2, mu = 1.8e-5;
    const std::vector<double> x = uniform(0.0, 0.02, 2001);
    std::vector<double> W(x.size());
    for (std::size_t k = 0; k < x.size(); ++k) W[k] = a * x[k];
    const BoundaryLayerProfile p = march_surface(incompressible(x, W, rho, mu), {});
    CHECK(p.transition_index == -1);
    const double theta0 = std::sqrt(K0 * mu / rho / a);
    CHECK(sref::measure("bl stagnation Lambda at x = 0.02 vs Lambda0", std::abs(p.Lambda.back() - L0)) <= 0.02);  // measured 2.4e-3
    CHECK(sref::measure("bl stagnation theta at x = 0.02 vs sqrt(K0 nu/a)", std::abs(p.theta.back() / theta0 - 1.0)) <=
          1e-3);  // measured 1.1e-4
}

TEST_CASE("boundary layer: laminar separation at Lambda = -12 in a retarded flow", "[turbolenza][bl]")
{
    // Howarth's linearly retarded flow W = W0 (1 - x/L): Lambda reaches -12 (tau_w = 0 in the 2 + L/6 form);
    // the separation x/L of Pohlhausen's method is about 0.156 (Schlichting, from memory: not verified here)
    const double W0 = 10.0, Lr = 1.0, rho = 1.2, mu = 1.8e-5;
    const std::vector<double> x = uniform(0.0, 0.3, 3001);
    std::vector<double> W(x.size());
    for (std::size_t k = 0; k < x.size(); ++k) W[k] = W0 * (1.0 - x[k] / Lr);
    BoundaryLayerSettings s;
    s.transition_re_theta = 1e9;  // keep it laminar
    const BoundaryLayerProfile p = march_surface(incompressible(x, W, rho, mu), s);
    REQUIRE(p.laminar_separation_index > 0);
    const std::size_t k = static_cast<std::size_t>(p.laminar_separation_index);
    CHECK(p.Lambda[k] == -12.0);
    CHECK(p.regime[k] == BlRegime::LaminarSeparated);
    CHECK(p.Lambda[k - 1] > -12.0);
    CHECK(p.cf[k] == 0.0);  // C1 = 2 - 12/6 = 0
    sref::measure("bl Howarth laminar separation x/L (Pohlhausen ~0.156)", x[k] / Lr);
    CHECK(std::abs(x[k] / Lr - 0.156) < 0.02);
}

TEST_CASE("boundary layer: transition at Re_theta 250 and the turbulent flat plate", "[turbolenza][bl]")
{
    // Re_L = 1.2 * 50 * 2 / 1.8e-5 = 6.7e6
    const double rho = 1.2, mu = 1.8e-5, W = 50.0, L = 2.0;
    const std::vector<double> x = uniform(0.0, L, 4001);
    const BoundaryLayerProfile p = march_surface(incompressible(x, std::vector<double>(x.size(), W), rho, mu), {});
    REQUIRE(p.transition_index > 0);
    const std::size_t t = static_cast<std::size_t>(p.transition_index);
    CHECK(p.Re_theta[t - 1] <= 250.0);
    CHECK(p.Re_theta[t] > 250.0);
    CHECK(p.regime[t] == BlRegime::Transition);
    CHECK(p.regime[t - 1] == BlRegime::Laminar);
    CHECK(p.regime[t + 1] == BlRegime::Turbulent);
    // the transition station carries the laminar theta and delta - delta* = delta' (7/10 + Lambda/120) [5-136 .. 5-138]
    const LaminarThicknesses lam = laminar_thicknesses(p.Lambda[t], p.theta[t], W * W / (2.0 * 1e12 * 300.0));
    CHECK(std::abs((p.delta[t] - p.delta_star[t]) / lam.delta_minus_delta_star - 1.0) < 1e-14);
    // Re_theta 250 on the flat plate: x_t = (250 / 0.6854)^2 nu / W
    const double x_t = std::pow(250.0 / 0.6854, 2.0) * mu / rho / W;
    sref::measure("bl flat plate transition x [m] (Pohlhausen 0.0399)", x[t]);
    CHECK(std::abs(x[t] / x_t - 1.0) < 0.05);  // measured 2.7e-2 (dx = 5e-4 m: the first-order start, about 2 stations)
    // turbulent: the 1/7-power law theta = 0.036 x Re_x^-0.2 (from x = 0), c_f = 0.0576 Re_x^-0.2;
    // Head's method from the transition point is another model: compare at x = L
    const double re = rho * W * L / mu;
    const double theta_power = 0.036 * L * std::pow(re, -0.2);
    const double err_theta = std::abs(p.theta.back() / theta_power - 1.0);
    sref::measure("bl turbulent flat plate theta(L) / 1/7-law - 1", p.theta.back() / theta_power - 1.0);
    CHECK(err_theta <= 0.1);  // measured -4.9e-2: Head / Green against the power law, two models
    sref::measure("bl turbulent flat plate H(L)", p.H.back());
    CHECK(p.H.back() > 1.25);
    CHECK(p.H.back() < 1.5);
    sref::measure("bl turbulent flat plate c_f(L) / 0.0576 Re^-0.2 - 1", p.cf.back() / (0.0576 * std::pow(re, -0.2)) - 1.0);
    CHECK(p.separation_index == -1);
}

TEST_CASE("boundary layer: the H_k limit, the wake and the loss", "[turbolenza][bl]")
{
    // a turbulent layer under a strong deceleration: H_k reaches 2.4, is held there in the closures
    const double rho = 1.2, mu = 1.8e-5;
    const std::vector<double> x = uniform(0.0, 1.0, 2001);
    std::vector<double> W(x.size());
    for (std::size_t k = 0; k < x.size(); ++k) W[k] = x[k] < 0.3 ? 60.0 : 60.0 * (1.0 - 0.6 * (x[k] - 0.3));
    EdgeConditions e = incompressible(x, W, rho, mu);
    e.total_temperature = 330.0;  // T't / T_e = 1.1: Green's H = (H_k + 1) 1.1 - 1
    e.cp = 1004.5;
    const BoundaryLayerProfile p = march_surface(e, {});
    REQUIRE(p.separation_index > 0);
    const std::size_t s = static_cast<std::size_t>(p.separation_index);
    CHECK(p.regime[s] == BlRegime::Separated);
    CHECK(p.Hk[s] >= 2.4);
    CHECK(p.Hk[s - 1] < 2.4);
    CHECK(std::abs(p.H[s] - ((2.4 + 1.0) * 1.1 - 1.0)) < 1e-14);  // the limited H_k in [5-131]
    sref::measure("bl decelerated turbulent separation x [m]", x[s]);
    // the closures
    CHECK(std::abs(head_hk(head_h1(1.7)) - 1.7) < 1e-14);
    CHECK(std::isinf(head_hk(3.3)));
    CHECK(head_entrainment(2.4) == 0.025 * 1.4);
    CHECK(green_h(1.4, 1.0) == 1.4);

    // the wake at constant W: tau_w = 0 and dW = 0 keep rho W^2 theta, so theta is constant;
    // (delta - delta*) grows by dx E
    EdgeConditions w = incompressible(uniform(0.0, 0.1, 11), std::vector<double>(11, 40.0), rho, mu);
    const BoundaryLayerProfile wake = march_wake(w, 1e-3, 5e-3, 2e-3, {});
    for (std::size_t k = 1; k < 11; ++k) {
        CHECK(std::abs(wake.theta[k] / 1e-3 - 1.0) < 1e-14);
        const double q = wake.delta[k] - wake.delta_star[k], q0 = wake.delta[k - 1] - wake.delta_star[k - 1];
        CHECK(std::abs(q - q0 - 0.01 * head_entrainment(std::min(head_hk(q / 1e-3), 2.4))) < 1e-15);
        CHECK(wake.cf[k] == 0.0);
    }
    CHECK(wake.delta_star[0] == 2e-3);  // station 0 holds the given sums

    // Eqs. 5-139 .. 5-141 by hand: S cos(beta_out) = 0.1 cos(20 deg)
    const double bi = 45.0 * std::numbers::pi / 180.0, bo = 20.0 * std::numbers::pi / 180.0;
    const ProfileLoss l = profile_loss(0.002, 0.005, 0.1, bi, bo);
    const double Th = 0.002 / (0.1 * std::cos(bo)), D = 0.005 / (0.1 * std::cos(bo));
    const double r = std::cos(bi) / std::cos(bo);
    CHECK(std::abs(l.omega - r * r * (2.0 * Th + D * D) / ((1.0 - D) * (1.0 - D))) < 1e-15);
    CHECK_THROWS_AS(profile_loss(0.002, 0.1, 0.1, bi, bo), BoundaryLayerError);  // Delta* >= 1
}

TEST_CASE("boundary layer: refusals", "[turbolenza][bl]")
{
    const EdgeConditions good = incompressible({0.0, 0.01, 0.02}, {10.0, 10.0, 10.0}, 1.2, 1.8e-5);
    CHECK_NOTHROW(march_surface(good, {}));
    EdgeConditions e = good;
    e.x = {0.0, 0.01, 0.01};
    CHECK_THROWS_AS(march_surface(e, {}), BoundaryLayerError);  // x not increasing
    e = good;
    e.W[1] = 0.0;
    CHECK_THROWS_AS(march_surface(e, {}), BoundaryLayerError);  // W = 0 after the first station
    e = good;
    e.rho.pop_back();
    CHECK_THROWS_AS(march_surface(e, {}), BoundaryLayerError);
    e = good;
    e.x = {0.0};
    e.W = {10.0};
    e.rho = {1.2};
    e.T = {300.0};
    e.mu = {1.8e-5};
    CHECK_THROWS_AS(march_surface(e, {}), BoundaryLayerError);  // one station
    BoundaryLayerSettings s;
    s.separation_hk = 1.0;
    CHECK_THROWS_AS(march_surface(good, s), BoundaryLayerError);
    CHECK_THROWS_AS(march_wake(good, 0.0, 1e-3, 1e-3, {}), BoundaryLayerError);
}
