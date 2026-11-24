// The end-wall boundary layers of Aungier (2003) Ch. 8 (see Endwall.hpp for every equation, its
// printed page, and every choice that is ours).
#include "turbolenza/Endwall.hpp"

#include "turbolenza/BoundaryLayer.hpp"
#include "turbolenza/Detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <numbers>
#include <string>

namespace BladeGen::Turbolenza {

using detail::py_repr;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kTwoPi = 2.0 * std::numbers::pi;
constexpr double kExtendedM = 0.05;  // the switch between 8-19 and 8-28 (p. 184)

bool positive(double x) noexcept { return std::isfinite(x) && x > 0.0; }

void check_edge(const EndwallWallStation& w, double omega_w, const std::string& where)
{
    if (!(positive(w.r) && std::isfinite(w.phi) && std::isfinite(w.m) && positive(w.Vm) && std::isfinite(w.Ctheta) &&
          positive(w.rho) && std::isfinite(w.P) && positive(w.mu) && std::isfinite(omega_w)))
        throw EndwallError(std::format("{}: the wall / edge state is not admissible (r {}, phi {}, m {}, Vm {}, Ctheta {}, "
                                       "rho {}, P {}, mu {}, omega_w {}): finite, r, Vm, rho, mu > 0",
                                       where, py_repr(w.r), py_repr(w.phi), py_repr(w.m), py_repr(w.Vm),
                                       py_repr(w.Ctheta), py_repr(w.rho), py_repr(w.P), py_repr(w.mu),
                                       py_repr(omega_w)));
}

void check_settings(const EndwallSettings& s)
{
    if (!(std::isfinite(s.shape_factor_limit) && s.shape_factor_limit > 1.0))
        throw EndwallError(std::format("the shape factor limit {} must be > 1", py_repr(s.shape_factor_limit)));
    if (!positive(s.re_theta_min))
        throw EndwallError(std::format("the Re_theta limit {} must be > 0", py_repr(s.re_theta_min)));
    if (!(std::isfinite(s.theta12_precaution) && s.theta12_precaution > 0.0 && s.theta12_precaution < 1.0))
        throw EndwallError(std::format("the theta12 precaution {} must lie in (0, 1)", py_repr(s.theta12_precaution)));
    if (!(positive(s.tolerance) && s.max_iterations > 0))
        throw EndwallError("the tolerance must be > 0 and max_iterations > 0");
}

void check_row(const EndwallRowWall& w, const std::string& where)
{
    if (!(std::isfinite(w.omega) && positive(w.pitch) && std::isfinite(w.stagger) && std::cos(w.stagger) > 0.0 &&
          std::isfinite(w.clearance) && w.clearance >= 0.0 && w.blade_count >= 0))
        throw EndwallError(std::format("{}: omega {} finite, pitch {} > 0, |stagger| {} < pi/2, clearance {} >= 0", where,
                                       py_repr(w.omega), py_repr(w.pitch), py_repr(w.stagger), py_repr(w.clearance)));
    if (!std::isfinite(w.seal_leakage) || (!w.shrouded && w.seal_leakage != 0.0))
        throw EndwallError(std::format("{}: the seal leakage {} needs a shrouded wall", where, py_repr(w.seal_leakage)));
}

// the bracket of 8-6 times rho_e sin(phi): rho_e sin(phi) V_theta_e [V_theta_e (delta2* + theta22) + 2 omega_w r delta2*]
double centrifugal_term(const EndwallLayer& L, const EndwallWallStation& w) noexcept
{
    return w.rho * std::sin(w.phi) * L.Vtheta_e *
           (L.Vtheta_e * (L.delta2 + L.theta22) + 2.0 * L.wall_omega * w.r * L.delta2);
}

// the interval's edge blade forces [8-8, 8-9] (ours: interval means and difference quotients)
struct IntervalForces
{
    double dm, r_mean, f_m, f_theta;
};

IntervalForces interval_forces(const EndwallWallStation& a, const EndwallWallStation& b) noexcept
{
    IntervalForces f{};
    f.dm = b.m - a.m;
    f.r_mean = 0.5 * (a.r + b.r);
    const double rho = 0.5 * (a.rho + b.rho), vm = 0.5 * (a.Vm + b.Vm);
    const double sin_phi = 0.5 * (std::sin(a.phi) + std::sin(b.phi));
    const double c2 = 0.5 * (a.Ctheta * a.Ctheta + b.Ctheta * b.Ctheta);  // (V_theta_e + omega_w r)^2 = C_theta_e^2
    f.f_m = rho * vm * (b.Vm - a.Vm) / f.dm + (b.P - a.P) / f.dm - sin_phi / f.r_mean * rho * c2;  // [8-8]
    f.f_theta = rho * vm / f.r_mean * (b.r * b.Ctheta - a.r * a.Ctheta) / f.dm;                       // [8-9]
    return f;
}

EndwallState state_of(const EndwallLayer& L) noexcept { return {L.Q, L.theta11, L.q12}; }
EndwallState scaled(const EndwallState& s, double k) noexcept { return {s.Q * k, s.theta11 * k, s.q12 * k}; }

// 8-73: rho_e V_e theta11 / mu >= 250 (ours: the state scaled, H1 and m kept)
void impose_re_theta11(EndwallLayer& L, const EndwallWallStation& w, double held_m, const EndwallSettings& set)
{
    if (!(L.Re_theta11 < set.re_theta_min)) return;
    const EndwallState s = scaled(state_of(L), set.re_theta_min / L.Re_theta11);
    const int it = L.iterations;
    const double nu = L.nu, dnu = L.dnu;
    L = endwall_closures(s, w, L.wall_omega, held_m, set);
    L.iterations = it;
    L.nu = nu;
    L.dnu = dnu;
    L.re_theta11_limited = true;
}

// One interval k-1 -> k in the frame omega_w = b.wall_omega with the force defect nu (steps 2-6 of
// p. 195; the integration rule is ours, Endwall.hpp).  `up` is the layer at k-1 (any frame).
EndwallLayer interval(const EndwallLayer& up, const EndwallWallStation& a, const EndwallWallStation& b, double nu,
                      const EndwallSettings& set, const std::string& where)
{
    if (!(b.m > a.m))
        throw EndwallError(std::format("{}: the arc length does not increase ({} after {})", where, py_repr(b.m),
                                       py_repr(a.m)));
    const double omega = b.wall_omega;
    const EndwallLayer Ca = endwall_closures(state_of(up), a, omega, up.m, set);  // k-1 in this frame (Sec. 8.7)
    const IntervalForces f = interval_forces(a, b);
    const double dm = f.dm;

    // the fluxes under d/dm at k-1 [8-5, 8-6, 8-7]
    const double Fa1 = a.r * a.rho * a.Vm * Ca.Q;
    const double Fa2 = a.r * a.rho * a.Vm * a.Vm * Ca.theta11;
    const double Fa3 = a.r * a.r * a.rho * a.Vm * Ca.q12;
    // the other terms at k-1
    const double Ea = a.r * a.rho * Ca.V_e * Ca.E;                 // r rho_e V_e E                    [8-5]
    const double Ma = a.r * Ca.tau_mw + centrifugal_term(Ca, a);  // r tau_mw + rho_e V_theta_e sin(phi) [..] [8-6]
    const double Ta = a.r * a.r * Ca.tau_thetaw;                   // r^2 tau_theta_w                  [8-7]
    const double Da = Ca.delta1 * a.r * a.rho * a.Vm;              // delta1* r rho_e V_me
    const double dVm = b.Vm - a.Vm;                                // (dV_me/dm) dm                    [8-6]
    const double drC = b.r * b.Ctheta - a.r * a.Ctheta;            // (d(r C_theta_e)/dm) dm           [8-7 = 8-9]
    const double c1 = b.r * b.rho * b.Vm, c2 = c1 * b.Vm, c3 = b.r * c1;
    const double nu_m = dm * f.r_mean * f.f_m * nu;                 // r f_me nu_m dm                  [8-6]
    const double nu_t = dm * f.r_mean * f.r_mean * f.f_theta * nu;  // r^2 f_theta_e nu_theta dm       [8-7]

    // the initial estimate (step 1): 8-74, 8-75, H1 held, 8-23
    EndwallState X{0.0, Fa2 / c2, Fa3 / c3};
    X.Q = 2.0 * Ca.H1 * X.theta11 / (Ca.H1 - 1.0);
    for (int it = 1; it <= set.max_iterations; ++it) {
        const EndwallLayer Cb = endwall_closures(X, b, omega, Ca.m, set);  // steps 2, 4, 5
        const double Eb = b.r * b.rho * Cb.V_e * Cb.E;
        const double Mb = b.r * Cb.tau_mw + centrifugal_term(Cb, b);
        const double Tb = b.r * b.r * Cb.tau_thetaw;
        const double D = 0.5 * (Da + Cb.delta1 * b.r * b.rho * b.Vm);
        EndwallState Y;  // step 3: 8-5, 8-6, 8-7 between k-1 and k
        Y.Q = (Fa1 + 0.5 * dm * (Ea + Eb)) / c1;
        Y.theta11 = (Fa2 - D * dVm + 0.5 * dm * (Ma + Mb) + nu_m) / c2;
        Y.q12 = (Fa3 - D * drC + 0.5 * dm * (Ta + Tb) + nu_t) / c3;
        if (!(positive(Y.Q) && positive(Y.theta11) && std::isfinite(Y.q12)))
            throw EndwallError(std::format("{}: no admissible layer (delta - delta1* {}, theta11 {}, V_theta theta12 {} "
                                           "at pass {})",
                                           where, py_repr(Y.Q), py_repr(Y.theta11), py_repr(Y.q12), it));
        const double tol = set.tolerance;
        const bool done = std::abs(Y.Q - X.Q) <= tol * Y.Q && std::abs(Y.theta11 - X.theta11) <= tol * Y.theta11 &&
                          std::abs(Y.q12 - X.q12) <= tol * (std::abs(Y.q12) + Cb.V_e * Y.theta11);
        X = Y;
        if (done) {  // step 6
            EndwallLayer L = endwall_closures(X, b, omega, Ca.m, set);
            L.iterations = it;
            L.nu = nu;
            impose_re_theta11(L, b, Ca.m, set);
            return L;
        }
    }
    throw EndwallError(std::format("{}: the station did not converge in {} passes", where, set.max_iterations));
}

// A blade row over the walls w[0] (leading edge) .. w.back() (trailing edge) from the layer `in` at w[0]:
// the inlet leakage (step 1), the march with one nu iterated by 8-44, 8-55, 8-57 (steps 2-6), the
// discharge leakage (step 7).  Returns the layers at w[1] .. w.back().
struct RowMarch
{
    std::vector<EndwallLayer> layers;
    EndwallRowWallResult res;
};

RowMarch row_march(const EndwallLayer& in, const std::vector<const EndwallWallStation*>& w, const EndwallRowWall& rw,
                   double omega_c, const EndwallSettings& set, const std::string& where)
{
    check_row(rw, where);
    const std::size_t n_int = w.size() - 1;
    const EndwallWallStation &a = *w.front(), &b = *w.back();
    if (rw.shrouded)
        for (std::size_t j = 1; j < w.size(); ++j)
            if (w[j]->wall_omega != rw.omega)
                throw EndwallError(std::format("{}: a shrouded wall rotates with the blades (p. 193), but a segment of "
                                               "the row rotates at {} (blades {})",
                                               where, py_repr(w[j]->wall_omega), py_repr(rw.omega)));
    RowMarch out;
    EndwallRowWallResult& res = out.res;
    res.shrouded = rw.shrouded;
    res.g = staggered_spacing(rw.pitch, rw.stagger);  // [8-2]

    // step 1 at the row inlet: the seal leakage [8-58 .. 8-62]
    EndwallLayer start = in;
    if (rw.shrouded) {
        const double Vt = a.Ctheta - a.r * rw.omega;  // the shroud frame is the row's (p. 193)
        res.leakage.in = state_of(in);
        res.leakage.delta_plus = leakage_thickness(rw.seal_leakage, a.r, a.rho, a.Vm);  // [8-59]
        res.leakage.U_leak_in = leakage_velocity(a.r, omega_c, Vt, a.Ctheta);          // [8-61]
        res.leakage.plus = seal_leakage_inlet(res.leakage.in, res.leakage.delta_plus, Vt, res.leakage.U_leak_in);
        start.Q = res.leakage.plus.Q;
        start.theta11 = res.leakage.plus.theta11;
        start.q12 = res.leakage.plus.q12;
    }

    // xi [8-53] (ours: the row's sum of r f_theta_e dm, Endwall.hpp)
    double rf_dm = 0.0;
    for (std::size_t j = 0; j < n_int; ++j) {
        const IntervalForces f = interval_forces(*w[j], *w[j + 1]);
        rf_dm += f.r_mean * f.f_theta * f.dm;
    }
    const double r_mean = 0.5 * (a.r + b.r);
    const double Wt_in = a.Ctheta - a.r * rw.omega, Wt_out = b.Ctheta - b.r * rw.omega;  // relative to the row
    const double W_mean = 0.5 * (Wt_in + Wt_out), Vm_mean = 0.5 * (a.Vm + b.Vm);
    res.xi = (r_mean / b.r + W_mean * Wt_out / (Vm_mean * b.Vm)) * rf_dm / (b.r * b.rho * b.Vm);

    const bool defect = set.force_defect == EndwallForceDefect::RowChapter8;
    // nu starts from nu0 of the row-inlet theta11, dnu = 0 (ours)
    double dnu = 0.0, nu = defect ? base_force_defect(res.g, rw.clearance, start.theta11) : 0.0;
    out.layers.resize(n_int);
    bool converged = false;
    for (int it = 1; it <= set.max_iterations && !converged; ++it) {
        const EndwallLayer* up = &start;
        for (std::size_t j = 0; j < n_int; ++j) {
            out.layers[j] = interval(*up, *w[j], *w[j + 1], nu, set, where);
            out.layers[j].dnu = dnu;
            up = &out.layers[j];
        }
        const EndwallLayer& minus = out.layers.back();
        res.iterations = it;
        res.nu = nu;
        res.collateral_defect = minus.q12 - Wt_out * minus.theta11;  // [8-45]: V_theta_e theta12 - W_theta_e theta11
        res.theta11_mean = 0.5 * (start.theta11 + minus.theta11);
        res.F = force_defect_factor(res.g, res.theta11_mean);  // [8-55]
        if (!defect) {
            converged = true;
            break;
        }
        res.nu0 = base_force_defect(res.g, rw.clearance, res.theta11_mean);  // [8-44]
        res.dnu = 0.0;
        if (set.force_defect_correction && res.xi != 0.0)
            res.dnu = (dnu - res.collateral_defect / res.xi) * res.F;  // [8-57]
        const double nu_new = res.nu0 + res.dnu;                     // [8-56]
        converged = std::abs(nu_new - nu) <= set.tolerance * res.g;  // step 6 on nu (theta11, H1 per station)
        dnu = res.dnu;
        nu = nu_new;
    }
    if (!converged)
        throw EndwallError(std::format("{}: the blade force defect did not converge in {} iterations", where,
                                       set.max_iterations));

    // step 7 at the row discharge [8-63 .. 8-66]
    if (rw.shrouded) {
        EndwallLayer& last = out.layers.back();
        const double Vt = b.Ctheta - b.r * rw.omega;
        res.leakage.minus = state_of(last);
        res.leakage.delta_minus = leakage_thickness(rw.seal_leakage, b.r, b.rho, b.Vm);  // [8-63]
        res.leakage.U_leak_out = leakage_velocity(b.r, omega_c, Vt, b.Ctheta);          // [8-61]
        res.leakage.out = seal_leakage_discharge(res.leakage.minus, res.leakage.delta_minus, Vt, res.leakage.U_leak_out);
        const int it = last.iterations;
        const double lnu = last.nu, ldnu = last.dnu;
        last = endwall_closures(res.leakage.out, b, b.wall_omega, last.m, set);
        last.iterations = it;
        last.nu = lnu;
        last.dnu = ldnu;
    }
    return out;
}

const EndwallWallStation& wall_of(const EndwallStation& s, bool hub) noexcept { return hub ? s.hub : s.shroud; }
const EndwallRowWall& wall_of(const EndwallBladeRow& r, bool hub) noexcept { return hub ? r.hub : r.shroud; }

}  // namespace

const char* endwall_entrainment_name(EndwallEntrainment e) noexcept
{
    return e == EndwallEntrainment::FreeStreamH ? "free_stream_H" : "meridional_H1";
}

const char* endwall_force_defect_name(EndwallForceDefect f) noexcept
{
    return f == EndwallForceDefect::None ? "none" : "row_chapter8";
}

// --- profiles (Sec. 8.3) ----------------------------------------------------------------------

EndwallProfile endwall_profile(double delta, double n, double m)
{
    if (!(positive(delta) && positive(n) && std::isfinite(m)))
        throw EndwallError(std::format("the profile delta {}, n {}, m {}: delta, n > 0, m finite", py_repr(delta),
                                       py_repr(n), py_repr(m)));
    EndwallProfile p;
    p.delta = delta;
    p.n = n;
    p.m = m;
    p.delta1 = delta * n / (n + 1.0);  // [8-24]
    p.H1 = 2.0 * n + 1.0;              // [8-22]
    p.theta11 = p.delta1 / p.H1;
    if (m >= kExtendedM) {
        p.H2 = 2.0 * m + 1.0;                                 // [8-25]
        p.delta2 = delta * m / (m + 1.0);                     // [8-26]
        p.theta22 = p.delta2 / p.H2;
        p.theta12 = delta * m / ((n + 1.0) * (n + m + 1.0));  // 8-21 solved for theta12
    } else {
        p.delta2 = delta * 20.0 * m / 21.0;                        // [8-29]
        p.theta22 = delta * (0.95 * m - 1.684 * m * m);            // [8-30] as printed
        p.H2 = endwall_H2(m);                                      // [8-32]
        p.theta12 = delta * (0.05 / ((n + 1.0) * (n + 1.05)) +
                             (m - 0.05) * 6.82 / ((n + 1.1) * (n + 2.1) * (n + 3.1)));  // 8-31 solved for theta12
    }
    return p;
}

double endwall_m_of_theta12(double theta12, double delta, double n) noexcept
{
    const double t_switch = kExtendedM * delta / ((n + 1.0) * (n + 1.05));  // theta12 of m = 0.05 (both forms)
    if (theta12 >= t_switch) {
        const double den = delta - theta12 * (n + 1.0);
        if (!(den > 0.0)) return kInf;
        return theta12 * (n + 1.0) * (n + 1.0) / den;  // [8-21]
    }
    return 0.05 + (theta12 / delta - 0.05 / ((n + 1.0) * (n + 1.05))) * (n + 1.1) * (n + 2.1) * (n + 3.1) / 6.82;  // [8-31]
}

double endwall_H2(double m) noexcept
{
    return m >= kExtendedM ? 2.0 * m + 1.0 : 20.0 / (21.0 * (0.95 - 1.684 * m));  // [8-25], [8-32]
}

// --- Sec. 8.4 -----------------------------------------------------------------------------------

EndwallStreamwise endwall_streamwise(double delta1, double delta2, double theta11, double theta12, double theta22,
                                     double beta) noexcept
{
    // p. 186: every thickness on the right of 8-40 and 8-41 limited to >= 0
    delta1 = std::max(delta1, 0.0);
    delta2 = std::max(delta2, 0.0);
    theta11 = std::max(theta11, 0.0);
    theta12 = std::max(theta12, 0.0);
    theta22 = std::max(theta22, 0.0);
    const double c = std::cos(beta), s = std::sin(beta);
    const double c2 = c * c, s2 = s * s;
    EndwallStreamwise r;
    r.delta_star = delta1 * c2 + delta2 * s2;  // [8-40]
    r.theta = (theta11 + delta1) * c2 * c2 + (theta22 + delta2) * s2 * s2 + 2.0 * (theta12 + delta1) * s2 * c2 -
              r.delta_star;  // [8-41]
    return r;
}

double ludwieg_tillmann_cf(double H, double re_theta) noexcept
{
    return 0.246 * std::exp(-1.561 * H) * std::pow(re_theta, -0.268);  // [8-36]
}

// --- the state ----------------------------------------------------------------------------------

EndwallState endwall_inlet_state(double delta1, double H1, double m, double Vtheta_e)
{
    if (!(positive(delta1) && std::isfinite(H1) && H1 > 1.0 && std::isfinite(m) && std::isfinite(Vtheta_e)))
        throw EndwallError(std::format("the inlet layer delta1* {}, H1 {}, m {}, V_theta_e {}: delta1* > 0, H1 > 1, "
                                       "finite",
                                       py_repr(delta1), py_repr(H1), py_repr(m), py_repr(Vtheta_e)));
    const double n = 0.5 * (H1 - 1.0);                  // [8-22]
    const double delta = delta1 * (n + 1.0) / n;        // [8-24]
    const EndwallProfile p = endwall_profile(delta, n, m);
    return {delta - delta1, delta1 / H1, Vtheta_e * p.theta12};
}

EndwallLayer endwall_closures(const EndwallState& s, const EndwallWallStation& w, double omega_w, double held_m,
                              const EndwallSettings& set)
{
    if (!(positive(s.Q) && positive(s.theta11) && std::isfinite(s.q12) && std::isfinite(held_m)))
        throw EndwallError(std::format("the end-wall state delta - delta1* {}, theta11 {}, V_theta theta12 {} (held m {}): "
                                       "positive thicknesses, finite",
                                       py_repr(s.Q), py_repr(s.theta11), py_repr(s.q12), py_repr(held_m)));
    check_edge(w, omega_w, "the end-wall closures");
    const double lim = set.shape_factor_limit;
    EndwallLayer L;
    L.Q = s.Q;
    L.theta11 = s.theta11;
    L.q12 = s.q12;
    L.wall_omega = omega_w;

    // step 4: n and delta [8-20, 8-22 .. 8-24]; H1 <= 2.4 (p. 184)
    const double h1_raw = s.Q > 2.0 * s.theta11 ? s.Q / (s.Q - 2.0 * s.theta11) : kInf;  // 8-23 solved for H1
    L.H1_limited = !(h1_raw <= lim);
    L.H1 = std::min(h1_raw, lim);
    L.n = 0.5 * (L.H1 - 1.0);  // [8-22] (= theta11 / (delta - delta1* - 2 theta11), 8-20, unlimited)
    L.delta1 = L.H1 * s.theta11;
    L.delta = s.Q + L.delta1;  // = (n + 1)(delta - delta1*) [8-24] unlimited

    // the edge relative to the wall [8-4]
    L.Vtheta_e = w.Ctheta - w.r * omega_w;
    L.V_e = std::hypot(w.Vm, L.Vtheta_e);
    L.beta_e = std::atan2(L.Vtheta_e, w.Vm);

    // theta12 and m: the step-4 precaution, 8-21 or 8-31, H2 <= 2.4
    double theta12_used = 0.0, m = held_m;
    if (L.Vtheta_e != 0.0) {
        L.theta12 = s.q12 / L.Vtheta_e;
        const double cap = set.theta12_precaution * L.delta / (L.n + 1.0);
        L.theta12_limited = L.theta12 > cap;
        theta12_used = std::min(L.theta12, cap);
        m = endwall_m_of_theta12(theta12_used, L.delta, L.n);
    }
    const double m_max = 0.5 * (lim - 1.0);
    L.H2_limited = m > m_max;
    L.m = std::min(m, m_max);
    const EndwallProfile p = endwall_profile(L.delta, L.n, L.m);  // step 5 [8-25, 8-26] or [8-29, 8-30, 8-32]
    L.delta2 = p.delta2;
    L.theta22 = p.theta22;
    L.H2 = p.H2;
    if (L.Vtheta_e == 0.0) {
        L.theta12 = p.theta12;
        theta12_used = p.theta12;
    }

    // Sec. 8.4
    const EndwallStreamwise sw = endwall_streamwise(L.delta1, L.delta2, L.theta11, theta12_used, L.theta22, L.beta_e);
    L.delta_star = sw.delta_star;
    L.theta = sw.theta;
    const double h_raw = sw.theta > 0.0 ? sw.delta_star / sw.theta : kInf;
    L.H_limited = !(h_raw <= lim);
    L.H = std::min(h_raw, lim);
    L.Re_theta = w.rho * L.V_e * std::max(sw.theta, 0.0) / w.mu;
    L.Re_theta11 = w.rho * L.V_e * s.theta11 / w.mu;
    L.re_theta_limited = L.Re_theta < set.re_theta_min;
    L.cf = ludwieg_tillmann_cf(L.H, std::max(L.Re_theta, set.re_theta_min));  // [8-36]
    // [8-43] / [8-42]: the form of BoundaryLayer's Head entrainment [5-133], 0.025 (x - 1)
    L.E = head_entrainment(set.entrainment == EndwallEntrainment::FreeStreamH ? L.H : L.H1);
    L.tau_mw = 0.5 * L.cf * w.rho * L.V_e * w.Vm;            // [8-34]
    L.tau_thetaw = 0.5 * L.cf * w.rho * L.V_e * L.Vtheta_e;  // [8-35]
    return L;
}

// --- Sec. 8.5 -----------------------------------------------------------------------------------

double staggered_spacing(double pitch, double stagger) noexcept { return pitch * std::cos(stagger); }  // [8-2]

double force_defect_shape(double x) noexcept
{
    const double x3 = x * x * x;
    return x3 / (1.0 + x3);
}

double base_force_defect(double g, double clearance, double theta11_mean) noexcept
{
    return (0.12 * g + 0.5 * clearance) * force_defect_shape(8.0 * theta11_mean / g);  // [8-44]
}

double force_defect_factor(double g, double theta11_mean) noexcept
{
    return 0.5 * force_defect_shape(8.0 * theta11_mean / g);  // [8-55]
}

// --- Sec. 8.6 -----------------------------------------------------------------------------------

double leakage_thickness(double mdot_leak, double r, double rho, double Vm) noexcept
{
    return mdot_leak / (kTwoPi * r * rho * Vm);  // [8-59], [8-63]
}

double leakage_velocity(double r, double omega_c, double Vtheta, double Ctheta) noexcept
{
    return 0.5 * r * omega_c + Vtheta - Ctheta;  // [8-61]
}

EndwallState seal_leakage_inlet(const EndwallState& in, double dp, double Vtheta_in, double U_leak) noexcept
{
    return {in.Q + dp,                                                       // [8-58]
            in.theta11 + dp,                                                 // [8-60]
            in.q12 + dp * Vtheta_in - 0.5 * (dp + std::abs(dp)) * U_leak};   // [8-62]
}

EndwallState seal_leakage_discharge(const EndwallState& minus, double dm, double Vtheta_out, double U_leak) noexcept
{
    return {minus.Q - dm,                                                      // [8-64]
            minus.theta11 - dm,                                                // [8-65]
            minus.q12 - dm * Vtheta_out + 0.5 * (dm - std::abs(dm)) * U_leak};  // [8-66]
}

// --- the analysis (Sec. 8.8) --------------------------------------------------------------------

EndwallLayer endwall_inlet(const EndwallWallStation& w, const EndwallInlet& inlet, double mass_flow,
                           const EndwallSettings& set)
{
    check_settings(set);
    check_edge(w, w.wall_omega, "the inlet station");
    if (!positive(mass_flow)) throw EndwallError(std::format("the mass flow {} must be > 0", py_repr(mass_flow)));
    if (!(positive(inlet.blockage) && inlet.blockage < 1.0))
        throw EndwallError(std::format("the inlet blockage {} must lie in (0, 1)", py_repr(inlet.blockage)));
    if (!(std::isfinite(inlet.H1) && inlet.H1 > 1.0 && inlet.H1 <= set.shape_factor_limit))
        throw EndwallError(std::format("the inlet H1 {} must lie in (1, {}]", py_repr(inlet.H1),
                                       py_repr(set.shape_factor_limit)));
    if (!(std::isfinite(inlet.m) && 2.0 * inlet.m + 1.0 <= set.shape_factor_limit))
        throw EndwallError(std::format("the inlet m {} must be finite with 2m + 1 <= {}", py_repr(inlet.m),
                                       py_repr(set.shape_factor_limit)));
    const double delta1 = inlet.blockage * mass_flow / (2.0 * kTwoPi * w.r * w.rho * w.Vm);  // [8-72], one wall
    const EndwallState s = endwall_inlet_state(delta1, inlet.H1, inlet.m, w.Ctheta - w.r * w.wall_omega);
    EndwallLayer L = endwall_closures(s, w, w.wall_omega, inlet.m, set);
    impose_re_theta11(L, w, inlet.m, set);  // [8-73]
    return L;
}

EndwallStepResult endwall_step(const EndwallLayer& prev, const EndwallWallStation& prev_wall,
                               const EndwallWallStation& wall, const EndwallRowWall* row, double compressor_omega,
                               const EndwallSettings& set)
{
    check_settings(set);
    check_edge(prev_wall, prev_wall.wall_omega, "the step's upstream station");
    check_edge(wall, wall.wall_omega, "the step's station");
    if (!std::isfinite(compressor_omega)) throw EndwallError("the compressor rotation speed must be finite");
    EndwallStepResult r;
    if (row == nullptr) {
        r.layer = interval(prev, prev_wall, wall, 0.0, set, "the end-wall step");
        return r;
    }
    const RowMarch rm = row_march(prev, {&prev_wall, &wall}, *row, compressor_omega, set, "the end-wall row step");
    r.layer = rm.layers.back();
    r.row = true;
    r.row_result = rm.res;
    return r;
}

EndwallBlockage endwall_blockage(const EndwallLayer& hub, const EndwallWallStation& h, const EndwallLayer& shroud,
                                 const EndwallWallStation& s, double mass_flow)
{
    if (!positive(mass_flow)) throw EndwallError(std::format("the mass flow {} must be > 0", py_repr(mass_flow)));
    EndwallBlockage b;
    b.B_hub = kTwoPi * h.r * h.rho * h.Vm * hub.delta1 / mass_flow;  // [8-76]
    b.B_shroud = kTwoPi * s.r * s.rho * s.Vm * shroud.delta1 / mass_flow;
    b.B = b.B_hub + b.B_shroud;
    b.K_B = 1.0 - b.B;  // [8-77]
    return b;
}

EndwallResult analyze_endwalls(const EndwallInput& in)
{
    check_settings(in.settings);
    const int N = static_cast<int>(in.stations.size());
    if (N < 2) throw EndwallError(std::format("{} station(s): the end-wall march needs two", N));
    if (!positive(in.mass_flow)) throw EndwallError(std::format("the mass flow {} must be > 0", py_repr(in.mass_flow)));
    if (!std::isfinite(in.compressor_omega)) throw EndwallError("the compressor rotation speed must be finite");
    for (const EndwallStation& st : in.stations)
        if (!(std::isnan(st.quasi_normal_length) || positive(st.quasi_normal_length)))
            throw EndwallError(std::format("the quasi-normal length {} must be > 0 (or NaN)",
                                           py_repr(st.quasi_normal_length)));
    int last_out = 0;
    for (std::size_t r = 0; r < in.rows.size(); ++r) {
        const EndwallBladeRow& R = in.rows[r];
        if (!(R.station_in >= last_out && R.station_out > R.station_in && R.station_out < N))
            throw EndwallError(std::format("blade row {}: stations {} -> {} (rows in order, not overlapping, in < out < {})",
                                           r, R.station_in, R.station_out, N));
        if (R.hub.omega != R.shroud.omega)
            throw EndwallError(std::format("blade row {}: one rotation speed for both walls ({} and {})", r,
                                           py_repr(R.hub.omega), py_repr(R.shroud.omega)));
        last_out = R.station_out;
    }

    EndwallResult res;
    res.stations.resize(static_cast<std::size_t>(N));
    res.rows.resize(in.rows.size());
    for (bool hub : {true, false}) {
        const char* name = hub ? "hub" : "shroud";
        auto W = [&](int k) -> const EndwallWallStation& { return wall_of(in.stations[static_cast<std::size_t>(k)], hub); };
        auto L = [&](int k) -> EndwallLayer& {
            EndwallStationResult& s = res.stations[static_cast<std::size_t>(k)];
            return hub ? s.hub : s.shroud;
        };
        for (int k = 0; k < N; ++k) check_edge(W(k), W(k).wall_omega, std::format("the {} at station {}", name, k));
        L(0) = endwall_inlet(W(0), in.inlet, in.mass_flow, in.settings);
        std::size_t next_row = 0;
        int k = 0;
        while (k < N - 1) {
            if (next_row < in.rows.size() && in.rows[next_row].station_in == k) {
                const EndwallBladeRow& R = in.rows[next_row];
                std::vector<const EndwallWallStation*> w;
                for (int j = R.station_in; j <= R.station_out; ++j) w.push_back(&W(j));
                const RowMarch rm = row_march(L(k), w, wall_of(R, hub), in.compressor_omega, in.settings,
                                              std::format("the {} in blade row {}", name, next_row));
                for (std::size_t j = 0; j < rm.layers.size(); ++j) L(k + 1 + static_cast<int>(j)) = rm.layers[j];
                (hub ? res.rows[next_row].hub : res.rows[next_row].shroud) = rm.res;
                k = R.station_out;
                ++next_row;
            } else {
                L(k + 1) = interval(L(k), W(k), W(k + 1), 0.0, in.settings, std::format("the {} at station {}", name, k + 1));
                ++k;
            }
        }
    }
    for (int k = 0; k < N; ++k) {
        EndwallStationResult& s = res.stations[static_cast<std::size_t>(k)];
        const EndwallStation& st = in.stations[static_cast<std::size_t>(k)];
        s.blockage = endwall_blockage(s.hub, st.hub, s.shroud, st.shroud, in.mass_flow);
        s.K_B_geometric = 1.0 - (s.hub.delta1 + s.shroud.delta1) / st.quasi_normal_length;  // ours
    }
    return res;
}

}  // namespace BladeGen::Turbolenza
