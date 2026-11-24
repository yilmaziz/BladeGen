// The blade-surface boundary layer of Aungier (2003) Sec. 5.6 = Aungier (2000) Ch. 15 (see
// BoundaryLayer.hpp for every equation and every choice that is ours).
#include "turbolenza/BoundaryLayer.hpp"

#include "turbolenza/Detail.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <string>

namespace BladeGen::Turbolenza {

using detail::nan;
using detail::py_repr;

double sutherland_viscosity(double T, const SutherlandLaw& law)
{
    if (!(std::isfinite(law.mu_ref) && law.mu_ref > 0.0 && std::isfinite(law.t_ref) && law.t_ref > 0.0 &&
          std::isfinite(law.s) && law.s > 0.0))
        throw BoundaryLayerError(std::format("the Sutherland constants mu_ref = {}, T_ref = {}, S = {} must be finite and "
                                             "positive",
                                             py_repr(law.mu_ref), py_repr(law.t_ref), py_repr(law.s)));
    if (!(std::isfinite(T) && T > 0.0))
        throw BoundaryLayerError(std::format("the temperature {} K: Sutherland's law needs T > 0", py_repr(T)));
    const double r = T / law.t_ref;
    return law.mu_ref * r * std::sqrt(r) * (law.t_ref + law.s) / (T + law.s);
}

const char* laminar_cf_form_name(LaminarCfForm f) noexcept
{
    return f == LaminarCfForm::Centrifugal2000 ? "centrifugal_2000" : "axial_2003_as_printed";
}

const char* bl_regime_name(BlRegime r) noexcept
{
    switch (r) {
    case BlRegime::Laminar: return "laminar";
    case BlRegime::LaminarSeparated: return "laminar_separated";
    case BlRegime::Transition: return "transition";
    case BlRegime::Turbulent: return "turbulent";
    case BlRegime::Separated: return "separated";
    }
    return "?";
}

// --- closures ----------------------------------------------------------------------------------

double gruschwitz_theta_ratio(double L) noexcept { return 37.0 / 315.0 - L / 945.0 - L * L / 9072.0; }  // [5-115]

double gruschwitz_k(double L) noexcept  // [5-124]
{
    const double g = gruschwitz_theta_ratio(L);
    return g * g * L;
}

// K is strictly increasing on [-12, 12): K' = g (37/315 - 3 L/945 - 5 L^2/9072), g > 0 there, and
// the bracket's root is L = 12 exactly.  So the root on [-12, 12] is unique: bisection to adjacent doubles.
double gruschwitz_lambda(double K)
{
    if (!std::isfinite(K)) throw BoundaryLayerError(std::format("K = {}: the Gruschwitz K must be finite", py_repr(K)));
    const double lo_l = -kGruschwitzLambdaLimit, hi_l = kGruschwitzLambdaLimit;
    if (K <= gruschwitz_k(lo_l)) return lo_l;  // Eq. 5-125: Lambda = -12, separation
    if (K >= gruschwitz_k(hi_l)) return hi_l;
    double lo = lo_l, hi = hi_l;
    while (true) {
        const double mid = 0.5 * (lo + hi);
        if (!(lo < mid && mid < hi)) break;
        const double k_mid = gruschwitz_k(mid);
        if (k_mid == K) return mid;
        if (k_mid < K) lo = mid;
        else hi = mid;
    }
    return std::abs(gruschwitz_k(lo) - K) <= std::abs(gruschwitz_k(hi) - K) ? lo : hi;
}

double laminar_wall_factor(double L, LaminarCfForm form) noexcept
{
    return (form == LaminarCfForm::Centrifugal2000 ? 2.0 : 1.0) + L / 6.0;  // [15-20] / [5-123]
}

LaminarThicknesses laminar_thicknesses(double L, double theta, double mach2) noexcept
{
    const double g = gruschwitz_theta_ratio(L);
    const double energy = (798048.0 - 4656.0 * L - 758.0 * L * L - 7.0 * L * L * L) / 4324320.0;  // [5-116]
    const double l = L / 100.0;
    const double F = 0.232912 - 0.831483 * l + 0.650584 * l * l + 17.8063 * l * l * l;           // [5-118]
    const double velocity = 3.0 / 10.0 - L / 120.0 + F * mach2;                                    // [5-117]
    const double enthalpy = mach2 * energy;                                                         // [5-119]
    LaminarThicknesses t;
    t.H = (enthalpy + velocity) / g;  // [5-120] / [5-115]
    t.delta_prime = theta / g;
    t.delta_star = t.H * theta;
    t.delta_minus_delta_star = t.delta_prime * (7.0 / 10.0 + L / 120.0);  // [5-138]
    return t;
}

double head_hk(double h1) noexcept
{
    if (!(h1 > 3.3)) return std::numeric_limits<double>::infinity();
    return 1.0 + std::pow(0.9 / (h1 - 3.3), 0.75);  // [5-132]
}

double head_h1(double hk) noexcept { return 3.3 + 0.9 / std::pow(hk - 1.0, 4.0 / 3.0); }

double head_entrainment(double hk) noexcept { return 0.025 * (hk - 1.0); }  // [5-133]

double green_h(double hk, double total_to_static) noexcept { return (hk + 1.0) * total_to_static - 1.0; }  // [5-131]

double green_cf(double hk, double h, double re_theta) noexcept
{
    const double inc = 0.246 * std::exp(-1.561 * hk) * std::pow(re_theta, -0.268);  // [5-134]
    return inc * hk * (hk + 1.0) / (2.0 * hk + h * (hk - 1.0));                     // [5-135]
}

// --- the march ---------------------------------------------------------------------------------

namespace {

void check_edge(const EdgeConditions& e, const char* what, bool first_speed_may_be_zero)
{
    const std::size_t n = e.x.size();
    if (n < 2) throw BoundaryLayerError(std::format("the {} has {} station(s): a march needs two", what, n));
    if (e.W.size() != n || e.rho.size() != n || e.T.size() != n || e.mu.size() != n)
        throw BoundaryLayerError(std::format("the {}'s edge conditions are not all of {} stations", what, n));
    if (!(std::isfinite(e.cp) && e.cp > 0.0 && std::isfinite(e.total_temperature) && e.total_temperature > 0.0))
        throw BoundaryLayerError(std::format("the {} needs cp > 0 and T't > 0 (got {}, {})", what, py_repr(e.cp),
                                             py_repr(e.total_temperature)));
    for (std::size_t k = 0; k < n; ++k) {
        const bool w_ok = std::isfinite(e.W[k]) && (e.W[k] > 0.0 || (k == 0 && first_speed_may_be_zero && e.W[k] == 0.0));
        if (!std::isfinite(e.x[k]) || !w_ok || !(std::isfinite(e.rho[k]) && e.rho[k] > 0.0) ||
            !(std::isfinite(e.T[k]) && e.T[k] > 0.0) || !(std::isfinite(e.mu[k]) && e.mu[k] > 0.0))
            throw BoundaryLayerError(std::format("the {}'s edge state at station {} is not admissible (x {}, W {}, rho {}, "
                                                 "T {}, mu {}): finite, W > 0, rho, T, mu > 0",
                                                 what, k, py_repr(e.x[k]), py_repr(e.W[k]), py_repr(e.rho[k]),
                                                 py_repr(e.T[k]), py_repr(e.mu[k])));
        if (k > 0 && !(e.x[k] > e.x[k - 1]))
            throw BoundaryLayerError(std::format("the {}'s arc length does not increase at station {} ({} after {})", what,
                                                 k, py_repr(e.x[k]), py_repr(e.x[k - 1])));
    }
}

void check_settings(const BoundaryLayerSettings& s)
{
    if (!(std::isfinite(s.transition_re_theta) && s.transition_re_theta > 0.0))
        throw BoundaryLayerError(std::format("transition Re_theta {}: a positive number", py_repr(s.transition_re_theta)));
    if (!(std::isfinite(s.separation_hk) && s.separation_hk > 1.0))
        throw BoundaryLayerError(std::format("separation H_k {}: a number > 1", py_repr(s.separation_hk)));
}

// The root of G on (0, inf) by bracketing from `guess` (G(0+) < 0 <= G(large)) and bisection to
// adjacent doubles; the endpoint with the smaller |G|.  `none` names the failure.
template <class Fn>
double positive_root(Fn&& G, double guess, const std::string& none)
{
    double lo = guess, hi = guess;
    int steps = 0;
    if (G(guess) < 0.0) {
        while (G(hi) < 0.0) {
            lo = hi;
            hi *= 2.0;
            if (++steps > 400 || !std::isfinite(hi)) throw BoundaryLayerError(none + ": the residual stays negative");
        }
    } else {
        while (G(lo) >= 0.0) {
            hi = lo;
            lo *= 0.5;
            if (++steps > 400 || !(lo > 0.0)) throw BoundaryLayerError(none + ": the residual stays positive");
        }
    }
    while (true) {
        const double mid = 0.5 * (lo + hi);
        if (!(lo < mid && mid < hi)) break;
        if (G(mid) < 0.0) lo = mid;
        else hi = mid;
    }
    return std::abs(G(lo)) <= std::abs(G(hi)) ? lo : hi;
}

BoundaryLayerProfile empty_profile(const EdgeConditions& e)
{
    BoundaryLayerProfile p;
    const std::size_t n = e.x.size();
    p.x = e.x;
    p.W = e.W;
    p.rho = e.rho;
    p.T = e.T;
    p.mu = e.mu;
    for (auto* v : {&p.theta, &p.delta_star, &p.delta, &p.H, &p.Hk, &p.cf, &p.Lambda, &p.Re_theta, &p.mass_deficit})
        v->assign(n, nan());
    p.regime.assign(n, BlRegime::Laminar);
    return p;
}

// The turbulent closures of station k at (theta, Q = delta - delta*), written into p.
struct TurbulentState
{
    double Hk_raw, Hk, H, delta_star, E, Re, cf;
};
TurbulentState turbulent_state(double theta, double Q, double rho, double W, double mu, double b0, double cap, bool wall)
{
    TurbulentState s{};
    s.Hk_raw = head_hk(Q / theta);
    s.Hk = std::min(s.Hk_raw, cap);  // 2003 p. 111: H_k limited to continue through separation
    s.H = green_h(s.Hk, b0);
    s.delta_star = s.H * theta;
    s.E = head_entrainment(s.Hk);
    s.Re = rho * W * theta / mu;
    s.cf = wall ? green_cf(s.Hk, s.H, s.Re) : 0.0;
    return s;
}

void store_turbulent(BoundaryLayerProfile& p, std::size_t k, double theta, double Q, const TurbulentState& s, double cap,
                     BlRegime regime)
{
    p.theta[k] = theta;
    p.delta_star[k] = s.delta_star;
    p.delta[k] = Q + s.delta_star;
    p.H[k] = s.H;
    p.Hk[k] = s.Hk_raw;
    p.cf[k] = s.cf;
    p.Re_theta[k] = s.Re;
    p.mass_deficit[k] = p.rho[k] * p.W[k] * s.delta_star;
    const bool separated = s.Hk_raw >= cap;
    p.regime[k] = (regime == BlRegime::Turbulent && separated) ? BlRegime::Separated : regime;
    if (separated && p.separation_index < 0) p.separation_index = static_cast<int>(k);
}

// Head / Green station k from the upstream station (2000 pp. 278-279: estimate theta and (delta -
// delta*), H_k, E, H and delta*, c_f, then Eqs. 15-1 and 15-24 for theta and (delta - delta*), to
// convergence).  Ours: the two station equations are solved by nested bisection, not substitution
// (which diverges where dx / theta is large: right after transition on a coarse grid).  For a theta,
// the entrainment residual Q - [N + dx rho W E(H_k(Q / theta))] / (rho W) increases with Q (E falls
// as H1 rises) and is negative at Q = 0 and >= 0 at Q = N / (rho W) + dx E_max: its root Q(theta) is
// unique.  The momentum residual theta - [M + dx tau_w - delta* rho W dW] / (rho W^2) at (theta,
// Q(theta)) is then bracketed from the upstream momentum flux held.
void turbulent_station(BoundaryLayerProfile& p, std::size_t k, double b0, double cap, bool wall, const char* what)
{
    const double rho = p.rho[k], W = p.W[k], mu = p.mu[k];
    const double dx = p.x[k] - p.x[k - 1], dW = W - p.W[k - 1];
    const double Qprev = p.delta[k - 1] - p.delta_star[k - 1];
    const double M = p.rho[k - 1] * p.W[k - 1] * p.W[k - 1] * p.theta[k - 1];  // rho W^2 theta upstream
    const double N = p.rho[k - 1] * p.W[k - 1] * Qprev;                        // rho W (delta - delta*) upstream
    const double q_max = N / (rho * W) + dx * head_entrainment(cap);          // E <= E(cap)
    const std::string where = std::format("the turbulent {} station {} (x = {:.6g} m)", what, k, p.x[k]);
    auto q_of = [&](double theta) {  // the root of the entrainment equation [5-127], Eq. 5-94
        double lo = 0.0, hi = q_max;
        while (true) {
            const double mid = 0.5 * (lo + hi);
            if (!(lo < mid && mid < hi)) break;
            const double E = turbulent_state(theta, mid, rho, W, mu, b0, cap, wall).E;
            if (mid - (N + dx * rho * W * E) / (rho * W) < 0.0) lo = mid;
            else hi = mid;
        }
        return hi;
    };
    auto G = [&](double theta) {  // the momentum equation [5-106], Eq. 5-94
        const TurbulentState s = turbulent_state(theta, q_of(theta), rho, W, mu, b0, cap, wall);
        const double tau = 0.5 * rho * W * W * s.cf;
        return theta - (M + dx * tau - s.delta_star * rho * W * dW) / (rho * W * W);
    };
    const double theta = positive_root(G, M / (rho * W * W), where + " has no momentum thickness");
    const double Q = q_of(theta);
    if (!(std::isfinite(theta) && theta > 0.0 && std::isfinite(Q) && Q > 0.0))
        throw BoundaryLayerError(std::format("{} has no positive solution (theta {}, delta - delta* {})", where,
                                             py_repr(theta), py_repr(Q)));
    store_turbulent(p, k, theta, Q, turbulent_state(theta, Q, rho, W, mu, b0, cap, wall), cap, BlRegime::Turbulent);
}

}  // namespace

BoundaryLayerProfile march_surface(const EdgeConditions& e, const BoundaryLayerSettings& settings)
{
    check_edge(e, "surface", true);
    check_settings(settings);
    BoundaryLayerProfile p = empty_profile(e);
    const std::size_t n = e.x.size();
    const double cap = settings.separation_hk;
    const LaminarCfForm form = settings.laminar_cf_form;

    // the first station: theta = 0 (2000 p. 276), K = 0, Lambda = 0
    {
        const LaminarThicknesses t = laminar_thicknesses(0.0, 0.0, e.W[0] * e.W[0] / (2.0 * e.cp * e.T[0]));
        p.theta[0] = 0.0;
        p.delta_star[0] = 0.0;
        p.delta[0] = 0.0;
        p.H[0] = t.H;
        p.Lambda[0] = 0.0;
        p.Re_theta[0] = 0.0;
        p.mass_deficit[0] = 0.0;
    }
    bool turbulent = false;
    for (std::size_t k = 1; k < n; ++k) {
        const double b0 = e.total_temperature / e.T[k];  // [5-121]
        if (turbulent) {
            turbulent_station(p, k, b0, cap, true, "surface");
            continue;
        }
        // Gruschwitz station k (2000 p. 276, steps 1-8)
        const double rho = e.rho[k], W = e.W[k], mu = e.mu[k];
        const double dx = e.x[k] - e.x[k - 1], dW = W - e.W[k - 1], dWdx = dW / dx;
        const double mach2 = W * W / (2.0 * e.cp * e.T[k]);
        const double M = e.rho[k - 1] * e.W[k - 1] * e.W[k - 1] * p.theta[k - 1];
        struct Eval
        {
            double lambda;
            LaminarThicknesses t;
            double tau;
            double F;
        };
        auto eval = [&](double theta) {
            Eval v{};
            v.lambda = gruschwitz_lambda(b0 * rho * theta * theta / mu * dWdx);  // [5-122], [5-124]
            v.t = laminar_thicknesses(v.lambda, theta, mach2);                    // [5-115 .. 5-120]
            v.tau = mu * W * laminar_wall_factor(v.lambda, form) / v.t.delta_prime;  // [15-20] / [5-123]
            v.F = (M + dx * v.tau - v.t.delta_star * rho * W * dW) / (rho * W * W);  // [5-106], Eq. 5-94
            return v;
        };
        auto G = [&](double theta) { return theta - eval(theta).F; };
        // the "safe" start (Exercise 5.6): from the first station the one-step flat-plate value of the
        // method, further on the upstream momentum flux held
        const double guess = p.theta[k - 1] == 0.0
                                 ? std::sqrt(laminar_wall_factor(0.0, form) * gruschwitz_theta_ratio(0.0) * mu * dx / (rho * W))
                                 : M / (rho * W * W);
        const double theta = positive_root(
            G, guess,
            std::format("the laminar station {} (x = {:.6g} m) has no momentum thickness (dW/W = {:.3g} over dx = {:.3g} m)",
                        k, e.x[k], dW / W, dx));
        const Eval v = eval(theta);
        const double re = rho * W * theta / mu;  // [5-126]
        p.theta[k] = theta;
        p.Lambda[k] = v.lambda;
        p.Re_theta[k] = re;
        if (re > settings.transition_re_theta) {  // transition [5-126]; theta and (delta - delta*) carried [5-136, 5-137]
            turbulent = true;
            p.transition_index = static_cast<int>(k);
            const double Q = v.t.delta_minus_delta_star;  // [5-138]
            store_turbulent(p, k, theta, Q, turbulent_state(theta, Q, rho, W, mu, b0, cap, true), cap,
                            BlRegime::Transition);
            continue;
        }
        p.delta_star[k] = v.t.delta_star;
        p.delta[k] = v.t.delta_minus_delta_star + v.t.delta_star;
        p.H[k] = v.t.H;
        p.cf[k] = 2.0 * v.tau / (rho * W * W);  // [5-109]
        p.mass_deficit[k] = rho * W * v.t.delta_star;
        const bool separated = v.lambda <= -kGruschwitzLambdaLimit;
        p.regime[k] = separated ? BlRegime::LaminarSeparated : BlRegime::Laminar;
        if (separated && p.laminar_separation_index < 0) p.laminar_separation_index = static_cast<int>(k);
    }
    return p;
}

BoundaryLayerProfile march_wake(const EdgeConditions& e, double theta0, double q0, double delta_star0,
                                const BoundaryLayerSettings& settings)
{
    check_edge(e, "wake", false);
    check_settings(settings);
    if (!(std::isfinite(theta0) && theta0 > 0.0 && std::isfinite(q0) && q0 > 0.0 && std::isfinite(delta_star0) &&
          delta_star0 >= 0.0))
        throw BoundaryLayerError(std::format("the wake starts from theta {}, delta - delta* {}, delta* {}: theta and "
                                             "delta - delta* must be positive, delta* >= 0",
                                             py_repr(theta0), py_repr(q0), py_repr(delta_star0)));
    BoundaryLayerProfile p = empty_profile(e);
    const double cap = settings.separation_hk;
    // station 0: the given sums (the closure's delta* would differ from the sum of the surfaces')
    TurbulentState s0 = turbulent_state(theta0, q0, e.rho[0], e.W[0], e.mu[0], e.total_temperature / e.T[0], cap, false);
    s0.delta_star = delta_star0;
    s0.H = delta_star0 / theta0;
    store_turbulent(p, 0, theta0, q0, s0, cap, BlRegime::Turbulent);
    for (std::size_t k = 1; k < e.x.size(); ++k)
        turbulent_station(p, k, e.total_temperature / e.T[k], cap, false, "wake");
    return p;
}

ProfileLoss profile_loss(double theta_sum, double delta_star_sum, double pitch, double beta_in, double beta_out)
{
    if (!(std::isfinite(theta_sum) && theta_sum >= 0.0 && std::isfinite(delta_star_sum) && delta_star_sum >= 0.0 &&
          std::isfinite(pitch) && pitch > 0.0 && std::isfinite(beta_in) && std::isfinite(beta_out) &&
          std::cos(beta_out) > 0.0 && std::cos(beta_in) > 0.0))
        throw BoundaryLayerError("the loss coefficient needs finite thicknesses >= 0, a positive pitch and flow angles "
                                 "within (-90, 90) deg");
    ProfileLoss l;
    const double width = pitch * std::cos(beta_out);
    l.Theta = theta_sum / width;       // [5-140]
    l.Delta_star = delta_star_sum / width;  // [5-141]
    if (!(l.Delta_star < 1.0))
        throw BoundaryLayerError(std::format("Delta* = {}: the displacement thickness fills the passage, Eq. 5-139 has no "
                                             "value",
                                             py_repr(l.Delta_star)));
    const double r = std::cos(beta_in) / std::cos(beta_out);
    l.omega = r * r * (2.0 * l.Theta + l.Delta_star * l.Delta_star) / ((1.0 - l.Delta_star) * (1.0 - l.Delta_star));  // [5-139]
    return l;
}

}  // namespace BladeGen::Turbolenza
