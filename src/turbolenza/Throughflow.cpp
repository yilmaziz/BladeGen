/// @file turbolenza/Throughflow.cpp
/// @brief Aungier (2003) Chapter 7: the meridional through-flow analysis (see the header for the
///        conventions).  Every equation is cited where it is used; "ours" marks a choice the
///        chapter does not make.
///
/// The procedure (Fig. 7-3 full normal equilibrium, Fig. 7-4 simplified models):
///  1. initial stream surfaces: Eqs. 7-54 to 7-56 for uniform flow with epsilon = 0 (p. 165);
///  2. on every quasi-normal, in order: link to the upstream one (Sec. 7.3), solve the flow field
///     (Sec. 7.2), the end-wall model; repeat until the flow field converges (the linking loop);
///     the simplified models reposition their stream surfaces inside this loop, undamped (p. 168);
///  3. full model only: reposition the stream surfaces (Sec. 7.4) with Novak's damping (Eq. 7-57,
///     B* adjusted by Eq. 7-59 when the position error grows) and go back to 2 until they converge.
/// The flow field of one quasi-normal (Sec. 7.2): Eq. 7-34 with f1 from Eq. 7-16 or 7-19, f4 from
/// Eqs. 7-17/7-18 or 7-20/7-21 and Eq. 7-33, dW_m/dm in f1 from Eq. 7-29 with Eq. 7-30; its
/// solution Eqs. 7-35, 7-36 with W_m(0) found from Eq. 7-5 (choke by Eq. 7-39); the static state
/// from Eqs. 7-37, 7-38 at constant entropy.
#include "turbolenza/Throughflow.hpp"

#include "turbolenza/BoundaryLayer.hpp"
#include "turbolenza/ThroughflowNumerics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <utility>

namespace BladeGen::Turbolenza::Throughflow {

namespace {

constexpr double kPi = std::numbers::pi;
using Vec = std::vector<double>;

std::size_t us(int k) { return static_cast<std::size_t>(k); }
int ii(std::size_t k) { return static_cast<int>(k); }
bool finite(double v) { return std::isfinite(v); }

[[noreturn]] void fail(const std::string& what) { throw ThroughflowError("throughflow: " + what); }

/// Linear interpolation in a strictly increasing table covering x (no extrapolation: refused).
double interp(const Vec& xs, const Vec& fs, double x)
{
    if (x < xs.front() - 1e-12 || x > xs.back() + 1e-12) fail("interpolation outside the table");
    if (x <= xs.front()) return fs.front();
    if (x >= xs.back()) return fs.back();
    const auto it = std::upper_bound(xs.begin(), xs.end(), x);
    const std::size_t k = static_cast<std::size_t>(it - xs.begin()) - 1;
    const double t = (x - xs[k]) / (xs[k + 1] - xs[k]);
    return fs[k] + t * (fs[k + 1] - fs[k]);
}

void require_table(const Vec& span, const char* who)
{
    if (span.size() < 2) fail(std::string(who) + ": at least two span points are needed");
    if (span.front() != 0.0 || span.back() != 1.0) fail(std::string(who) + ": the span points must run from 0 to 1");
    for (std::size_t k = 1; k < span.size(); ++k)
        if (!(span[k] > span[k - 1])) fail(std::string(who) + ": the span points must increase strictly");
}

/// Geometry of one straight quasi-normal (Eq. 7-1).
struct QnGeometry
{
    double zh = 0.0, rh = 0.0, zs = 0.0, rs = 0.0;
    double ys = 0.0, lambda = 0.0, cosl = 1.0, sinl = 0.0;
};

QnGeometry qn_geometry(const QuasiNormal& q)
{
    QnGeometry g;
    g.zh = q.z_hub;
    g.rh = q.r_hub;
    g.zs = q.z_shroud;
    g.rs = q.r_shroud;
    g.ys = std::hypot(g.zs - g.zh, g.rs - g.rh);
    g.lambda = std::atan2(g.zh - g.zs, g.rs - g.rh);  // Eq. 7-1: tan(lambda) = (z_h - z_s) / (r_s - r_h)
    g.cosl = (g.rs - g.rh) / g.ys;
    g.sinl = (g.zh - g.zs) / g.ys;
    return g;
}

/// Eqs. 7-55, 7-56 with y_h = 0.
void place(const QnGeometry& g, double y, double& z, double& r)
{
    z = g.zh + (g.zs - g.zh) * y / g.ys;
    r = g.rh + (g.rs - g.rh) * y / g.ys;
}

/// The row's blade height: the mean of its LE and TE line lengths (ours).
double blade_height(const BladeRowSpec& b)
{
    return 0.5 * (std::hypot(b.le_shroud_z - b.le_hub_z, b.le_shroud_r - b.le_hub_r) +
                  std::hypot(b.te_shroud_z - b.te_hub_z, b.te_shroud_r - b.te_hub_r));
}

/// Initial positions (Fig. 7-3 / 7-4 "Initialize stream surfaces", p. 165): Eq. 7-54 for uniform
/// flow and epsilon = 0, int_0^y r dy = psi int_0^ys r dy with r linear in y (closed form, ours).
Vec initial_positions(const QnGeometry& g, int nss)
{
    Vec y(us(nss));
    const double a = (g.rs - g.rh) / (2.0 * g.ys), b = g.rh;
    const double total = 0.5 * (g.rh + g.rs) * g.ys;
    for (int j = 0; j < nss; ++j) {
        const double psi = static_cast<double>(j) / static_cast<double>(nss - 1);
        const double c = psi * total;  // a y^2 + b y = c
        y[us(j)] = 2.0 * c / (b + std::sqrt(b * b + 4.0 * a * c));
    }
    y.front() = 0.0;
    y.back() = g.ys;
    return y;
}

/// Sheet values (sheet k between surfaces k and k + 1, located at its y-midpoint) carried to the
/// surfaces at second order: linear interpolation between the two neighbouring sheet midpoints
/// inside, linear extrapolation from the two nearest midpoints at the walls (ours).
Vec sheets_to_surfaces(const Vec& y, const Vec& v)
{
    const std::size_t n = y.size();
    Vec out(n);
    auto mid = [&](std::size_t k) { return 0.5 * (y[k] + y[k + 1]); };
    auto line = [&](std::size_t k, double x) {  // through sheets k and k + 1
        const double a = mid(k), b = mid(k + 1);
        return v[k] + (v[k + 1] - v[k]) * (x - a) / (b - a);
    };
    for (std::size_t j = 0; j < n; ++j) {
        if (j == 0) out[j] = line(0, y[0]);
        else if (j + 1 == n) out[j] = line(n - 3, y[j]);
        else out[j] = line(j - 1, y[j]);
    }
    return out;
}

/// Eq. 7-58 with the square root (the book prints m = int [1 + (dr/dz)^2] dz, a misprint: arc
/// length is dm = [1 + (dr/dz)^2]^(1/2) dz, Eq. 3-5); dr/dz by Eqs. 7-67 to 7-69 and the integral
/// by Eqs. 7-72, 7-73 (ours: the same three-point formulas, p. 165 "evaluated numerically").
Vec meridional_coordinate(const Vec& z, const Vec& r)
{
    for (std::size_t i = 1; i < z.size(); ++i)
        if (!(z[i] > z[i - 1])) fail("Eq. 7-58 needs z to increase along every stream surface");
    const Vec drdz = derivative(z, r, EndRule::ThreePoint);
    Vec g(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) g[i] = std::sqrt(1.0 + drdz[i] * drdz[i]);
    return cumulative_integral(z, g);
}

/// The flow specification on one quasi-normal for one flow iteration.
struct Spec
{
    double omega = 0.0;           ///< frame
    bool known_beta = false;      ///< Eqs. 7-19 to 7-21 (else 7-16 to 7-18)
    Vec I, s, rCt, tanb, KW, loss, row_loss;
    Vec KBlade;                   ///< Ch. 12 blade-metal blockage inside the integral (1 in this version)
};

/// The per-node state at a given W_m.
struct NodeState
{
    double Wm = 0.0, Wt = 0.0, Ct = 0.0, Hrel = 0.0, h = 0.0;
    GasState g;
};

struct Geo
{
    Vec y, z, r, m;
    Vec phi, eps;               ///< geometric: the mass equations (Eqs. 7-5, 7-32, 7-54) and the outputs
    Vec phi_eq, kappa_eq, eps_eq, dphidy;  ///< what the normal-momentum equation uses (the model level)
    double cosl = 1.0;          ///< cos(lambda) = dr/dy
};

/// The state at W_m on node j; valid = false when h <= 0 (not thrown: the W_m(0) search probes).
bool node_state(const ThroughflowGas& gas, const Spec& sp, const Geo& ge, std::size_t j, double Wm, NodeState& st)
{
    st.Wm = Wm;
    const double r = ge.r[j];
    st.Wt = sp.known_beta ? Wm * sp.tanb[j] : sp.rCt[j] / r - sp.omega * r;  // W_theta = C_theta - omega r (Eq. 7-44)
    st.Ct = st.Wt + sp.omega * r;
    st.Hrel = sp.I[j] + 0.5 * (r * sp.omega) * (r * sp.omega);   // Eq. 7-37
    st.h = st.Hrel - 0.5 * (Wm * Wm + st.Wt * st.Wt);            // Eq. 7-38
    if (!(st.h > 0.0) || !finite(st.h)) return false;
    st.g = gas_state(gas, st.h, sp.s[j]);
    return true;
}

struct MassEval
{
    bool valid = false;
    double mdot = 0.0, mean_mm = 0.0;
};

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Gas, small models
// ─────────────────────────────────────────────────────────────────────────────

void ThroughflowGas::validate() const
{
    if (!(gamma > 1.0) || !finite(gamma)) fail("gamma must be > 1");
    if (!(gas_constant > 0.0) || !finite(gas_constant)) fail("R must be > 0");
    if (!(viscosity > 0.0) || !finite(viscosity)) fail("the viscosity must be > 0");
    if (!(mu_ref > 0.0) || !(t_mu > 0.0) || !(s_mu > 0.0) || !finite(mu_ref) || !finite(t_mu) || !finite(s_mu))
        fail("Sutherland's constants must be > 0");
    if (!(t_ref > 0.0) || !(p_ref > 0.0) || !finite(t_ref) || !finite(p_ref)) fail("the entropy reference must be > 0");
}

double ThroughflowGas::entropy(double T, double P) const
{
    if (!(T > 0.0) || !(P > 0.0) || !finite(T) || !finite(P)) fail("entropy of a non-positive state");
    return cp() * std::log(T / t_ref) - gas_constant * std::log(P / p_ref);  // Eq. 2-25
}

double ThroughflowGas::pressure(double h, double s) const
{
    if (!(h > 0.0) || !finite(h) || !finite(s)) fail("pressure at a non-positive enthalpy");
    const double T = h / cp();
    return p_ref * std::exp((cp() * std::log(T / t_ref) - s) / gas_constant);  // Eq. 2-25 solved for P
}

double ThroughflowGas::mu(double T) const
{
    if (!(T > 0.0) || !finite(T)) fail("viscosity at a non-positive temperature");
    if (!sutherland) return viscosity;
    SutherlandLaw law;
    law.mu_ref = mu_ref;
    law.t_ref = t_mu;
    law.s = s_mu;
    return sutherland_viscosity(T, law);
}

GasState gas_state(const ThroughflowGas& gas, double h, double s)
{
    GasState g;
    g.T = h / gas.cp();
    g.P = gas.pressure(h, s);
    g.rho = g.P / (gas.gas_constant * g.T);
    g.a = std::sqrt(gas.gamma * gas.gas_constant * g.T);  // Eq. 2-28
    return g;
}

SpecifiedBladeRow::SpecifiedBladeRow(std::vector<double> span, std::vector<double> beta2, std::vector<double> loss,
                                     std::vector<double> wake_blockage)
    : span_(std::move(span)), beta2_(std::move(beta2)), loss_(std::move(loss)), kw_(std::move(wake_blockage))
{
    require_table(span_, "SpecifiedBladeRow");
    if (beta2_.size() != span_.size() || loss_.size() != span_.size()) fail("SpecifiedBladeRow: profile lengths differ");
    if (kw_.empty()) kw_.assign(span_.size(), 1.0);
    if (kw_.size() != span_.size()) fail("SpecifiedBladeRow: K_W length differs");
    for (std::size_t k = 0; k < span_.size(); ++k) {
        if (!finite(beta2_[k]) || !(std::abs(beta2_[k]) < 0.5 * kPi)) fail("SpecifiedBladeRow: |beta2| must be < pi/2");
        if (!finite(loss_[k])) fail("SpecifiedBladeRow: non-finite loss");
        if (!(kw_[k] > 0.0 && kw_[k] <= 1.0)) fail("SpecifiedBladeRow: K_W must be in (0, 1]");
    }
}

BladeRowDischarge SpecifiedBladeRow::evaluate(const BladeRowInlet& in) const
{
    BladeRowDischarge d;
    d.beta2 = interp(span_, beta2_, in.span);
    d.loss_coefficient = interp(span_, loss_, in.span);
    d.wake_blockage = interp(span_, kw_, in.span);
    return d;
}

ConstantEndwall::ConstantEndwall(double blockage_factor) : kb_(blockage_factor)
{
    if (!(kb_ > 0.0 && kb_ <= 1.0)) fail("ConstantEndwall: K_B must be in (0, 1]");
}

EndwallResult ConstantEndwall::evaluate(const EndwallStation& /*st*/)
{
    EndwallResult r;
    r.blockage_factor = kb_;
    return r;
}

std::vector<double> distribute_row_loss(double dPt, RowLossResult::Side side, const std::vector<double>& psi)
{
    if (!finite(dPt)) fail("distribute_row_loss: non-finite Delta P_t");
    std::vector<double> out(psi.size());
    for (std::size_t j = 0; j < psi.size(); ++j) {
        if (!(psi[j] >= 0.0 && psi[j] <= 1.0)) fail("distribute_row_loss: mass fractions must be in [0, 1]");
        const double t = side == RowLossResult::Side::Shroud ? psi[j] : 1.0 - psi[j];
        out[j] = 2.0 * dPt * t;  // linear, 0 at the opposite wall, mean over psi = dPt (p. 147; ours)
    }
    return out;
}

InletProfile InletProfile::uniform(double Pt, double Tt, Swirl swirl, double value)
{
    InletProfile p;
    p.span = {0.0, 1.0};
    p.Pt = {Pt, Pt};
    p.Tt = {Tt, Tt};
    p.swirl = swirl;
    p.swirl_value = {value, value};
    return p;
}

std::string to_string(ThroughflowModel m)
{
    switch (m) {
    case ThroughflowModel::FullNormalEquilibrium: return "full normal equilibrium";
    case ThroughflowModel::Approximate: return "approximate (linear phi, kappa_m)";
    case ThroughflowModel::SimpleNonIsentropic: return "simple non-isentropic normal equilibrium";
    case ThroughflowModel::SimpleIsentropic: return "simple (isentropic) normal equilibrium";
    }
    return "?";
}

std::string to_string(ThroughflowResult::Status s)
{
    switch (s) {
    case ThroughflowResult::Status::Converged: return "converged";
    case ThroughflowResult::Status::Choked: return "choked";
    case ThroughflowResult::Status::NotConverged: return "not converged";
    }
    return "?";
}

double novak_damping(double mach_m, double y_s, double dm, double bstar, double mach_cap)
{
    if (!(y_s > 0.0) || !(dm > 0.0) || !(bstar > 0.0) || !finite(mach_m) || !(mach_cap > 0.0 && mach_cap < 1.0))
        fail("novak_damping: y_s, dm, B* must be > 0 and the Mach cap in (0, 1)");
    const double Mm = std::min(std::abs(mach_m), mach_cap);               // p. 165: M_m <= 0.95
    return 1.0 / (1.0 + (1.0 - Mm * Mm) * y_s * y_s / (bstar * dm * dm));  // Eq. 7-57
}

double BstarControl::update(double error)
{
    if (error > (1.0 + threshold) * min_error) {
        bstar = 0.5 * (bstar + 1.0);  // Eq. 7-59
        min_error = error;
    } else {
        min_error = std::min(min_error, error);
    }
    return bstar;
}

double line_crossing_fraction(double az, double ar, double bz, double br, double z1, double r1, double z2, double r2)
{
    // a + t (b - a) = p1 + u (p2 - p1): Cramer's rule for t
    const double dz = bz - az, dr = br - ar, ez = z2 - z1, er = r2 - r1;
    const double det = dz * (-er) - dr * (-ez);
    if (!(std::abs(det) > 1e-300)) fail("a stream surface is parallel to a blade edge line");
    const double fz = z1 - az, fr = r1 - ar;
    return (fz * (-er) - fr * (-ez)) / det;
}

void ThroughflowCase::validate() const
{
    gas.validate();
    if (!(mass_flow > 0.0) || !finite(mass_flow)) fail("the mass flow must be > 0");
    if (!finite(compressor_omega)) fail("non-finite compressor speed");
    const int nqn = ii(quasi_normals.size());
    if (nqn < 3) fail("at least three quasi-normals are needed");
    for (int i = 0; i < nqn; ++i) {
        const auto& q = quasi_normals[us(i)];
        for (double v : {q.z_hub, q.r_hub, q.z_shroud, q.r_shroud, q.hub_wall_omega, q.shroud_wall_omega})
            if (!finite(v)) fail("quasi-normal " + std::to_string(i) + " has a non-finite value");
        for (const auto& v : {q.phi_hub, q.kappa_hub, q.phi_shroud, q.kappa_shroud})
            if (v && !finite(*v)) fail("quasi-normal " + std::to_string(i) + " has a non-finite wall slope or curvature");
        if (!(q.r_hub > 0.0) || !(q.r_shroud > 0.0)) fail("quasi-normal " + std::to_string(i) + ": r must be > 0");
        if (q.inside_blade_row)
            fail("quasi-normal " + std::to_string(i) + " lies inside a blade row: Chapter 7 places them outside (p. 155)");
        if (!(std::hypot(q.z_shroud - q.z_hub, q.r_shroud - q.r_hub) > 0.0))
            fail("quasi-normal " + std::to_string(i) + " has zero length");
        if (i > 0) {
            const auto& p = quasi_normals[us(i - 1)];
            if (!(q.z_hub > p.z_hub) || !(q.z_shroud > p.z_shroud))
                fail("the quasi-normal end points must move downstream (z increasing) on both walls");
        }
    }
    const auto& o = options;
    if (o.n_stream_surfaces < 3) fail("at least three stream surfaces are needed");
    if (!(o.bstar > 0.0) || !finite(o.bstar)) fail("B* must be > 0");
    if (!(o.bstar_growth_threshold >= 0.0) || !finite(o.bstar_growth_threshold)) fail("the B* growth threshold must be >= 0");
    if (!(o.max_mach_damping > 0.0 && o.max_mach_damping < 1.0)) fail("the damping Mach limit must be in (0, 1)");
    if (!(o.min_one_minus_mm2 > 0.0 && o.min_one_minus_mm2 <= 1.0)) fail("the Eq. 7-30 limit must be in (0, 1]");
    if (!(o.flow_tolerance > 0.0) || !(o.position_tolerance > 0.0) || !(o.mass_tolerance > 0.0)) fail("tolerances must be > 0");
    if (o.max_flow_iterations < 1 || o.max_outer_iterations < 1) fail("iteration limits must be >= 1");
    if (o.choke_repositionings < 0) fail("choke_repositionings must be >= 0");
    if (!(o.curvature_damping >= 0.0) || !finite(o.curvature_damping)) fail("the curvature damping D must be >= 0");
    // inlet
    require_table(inlet.span, "inlet profile");
    const std::size_t ns = inlet.span.size();
    if (inlet.Pt.size() != ns || inlet.Tt.size() != ns || inlet.swirl_value.size() != ns) fail("inlet profile lengths differ");
    for (std::size_t k = 0; k < ns; ++k) {
        if (!(inlet.Pt[k] > 0.0) || !(inlet.Tt[k] > 0.0) || !finite(inlet.Pt[k]) || !finite(inlet.Tt[k]))
            fail("inlet Pt and Tt must be > 0");
        if (!finite(inlet.swirl_value[k])) fail("non-finite inlet swirl");
        if (inlet.swirl == InletProfile::Swirl::Alpha && !(std::abs(inlet.swirl_value[k]) < 0.5 * kPi))
            fail("inlet |alpha| must be < pi/2");
    }
    // rows
    std::vector<int> used;
    for (std::size_t k = 0; k < rows.size(); ++k) {
        const auto& b = rows[k];
        const std::string who = "row " + std::to_string(k) + " (" + b.name + ")";
        if (b.qn_inlet < 0 || b.qn_inlet > nqn - 2) fail(who + ": qn_inlet out of range");
        if (b.qn_outlet != -1 && b.qn_outlet != b.qn_inlet + 1)
            fail(who + ": quasi-normals inside a blade row (Ch. 12) are not in this version");
        if (!(b.hub_clearance >= 0.0) || !(b.tip_clearance >= 0.0) || !finite(b.hub_clearance) || !finite(b.tip_clearance))
            fail(who + ": clearances must be >= 0");
        if (std::find(used.begin(), used.end(), b.qn_inlet) != used.end()) fail(who + ": two rows between the same quasi-normals");
        if (!used.empty() && !(b.qn_inlet > used.back()))
            fail(who + ": the rows must be listed in the flow direction (N_row of Eq. 6-87 is their order)");
        used.push_back(b.qn_inlet);
        if (!b.model) fail(who + ": no blade-row model");
        if (b.n_blades < 1) fail(who + ": n_blades must be >= 1");
        if (b.s_row != 1 && b.s_row != -1) fail(who + ": s_row must be +1 or -1");
        if (!finite(b.omega)) fail(who + ": non-finite omega");
        if (b.sections.size() < 2) fail(who + ": at least two sections are needed");
        Vec sp;
        for (const auto& s : b.sections) {
            sp.push_back(s.span);
            if (!(s.chord > 0.0) || !finite(s.chord)) fail(who + ": the chord must be > 0");
        }
        require_table(sp, "blade row sections");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// The solver
// ─────────────────────────────────────────────────────────────────────────────

namespace {

/// Stored per quasi-normal across passes.
struct QnData
{
    QnGeometry g;
    Vec y;                  ///< stream-surface positions (current)
    bool solved = false;    ///< has a flow solution (this or an earlier pass)
    Vec Wm;                 ///< last solution
    Vec H, rCt, s, T, P, rho, a, Ct;  ///< absolute linking data of the last solution
    WallState wall[2];      ///< hub, shroud of the last solution
    double KB = 1.0;
    ThroughflowStation out;
};

class Solver
{
public:
    explicit Solver(const ThroughflowCase& c)
        : c_(c), o_(c.options), nqn_(ii(c.quasi_normals.size())), nss_(c.options.n_stream_surfaces)
    {
        endwall_ = c.endwall ? c.endwall : std::make_shared<ConstantEndwall>(1.0);
        q_.resize(us(nqn_));
        for (int i = 0; i < nqn_; ++i) {
            q_[us(i)].g = qn_geometry(c.quasi_normals[us(i)]);
            q_[us(i)].y = initial_positions(q_[us(i)].g, nss_);
        }
        leak_.assign(c.rows.size(), 0.0);
        row_at_.assign(us(nqn_), -1);  // row between i - 1 and i
        row_from_.assign(us(nqn_), -1);  // row between i and i + 1
        for (std::size_t k = 0; k < c.rows.size(); ++k) {
            row_at_[us(c.rows[k].qn_inlet + 1)] = ii(k);
            row_from_[us(c.rows[k].qn_inlet)] = ii(k);
        }
        // the wall meridional coordinates (fixed: the wall stream surfaces stay on the contours)
        Vec zh(us(nqn_)), rh(us(nqn_)), zs(us(nqn_)), rs(us(nqn_));
        for (int i = 0; i < nqn_; ++i) {
            zh[us(i)] = q_[us(i)].g.zh;
            rh[us(i)] = q_[us(i)].g.rh;
            zs[us(i)] = q_[us(i)].g.zs;
            rs[us(i)] = q_[us(i)].g.rs;
        }
        m_hub_ = meridional_coordinate(zh, rh);
        m_shroud_ = meridional_coordinate(zs, rs);
        // the end-wall slope and curvature at every quasi-normal (the Approximate model): given, or
        // from the wall points by Eqs. 7-2 and 7-9 with Eq. 7-67 inside and Eq. 7-70 at the ends (ours)
        const Vec* mw[2] = {&m_hub_, &m_shroud_};
        const Vec* rw[2] = {&rh, &rs};
        for (int side = 0; side < 2; ++side) {
            const Vec drdm = derivative(*mw[side], *rw[side], EndRule::TwoPoint);
            Vec phi(us(nqn_));
            for (int i = 0; i < nqn_; ++i) {
                if (!(std::abs(drdm[us(i)]) < 1.0)) fail("Eq. 7-2: |dr/dm| >= 1 on a wall");
                phi[us(i)] = std::asin(drdm[us(i)]);
            }
            const Vec dphidm = derivative(*mw[side], phi, EndRule::TwoPoint);
            wall_phi_[side].resize(us(nqn_));
            wall_kappa_[side].resize(us(nqn_));
            for (int i = 0; i < nqn_; ++i) {
                const auto& q = c.quasi_normals[us(i)];
                const auto& gphi = side == 0 ? q.phi_hub : q.phi_shroud;
                const auto& gkap = side == 0 ? q.kappa_hub : q.kappa_shroud;
                wall_phi_[side][us(i)] = gphi ? *gphi : phi[us(i)];
                wall_kappa_[side][us(i)] = gkap ? *gkap : -dphidm[us(i)];
            }
        }
    }

    ThroughflowResult run();

private:
    bool full() const { return c_.model == ThroughflowModel::FullNormalEquilibrium; }
    bool simple() const
    {
        return c_.model == ThroughflowModel::SimpleNonIsentropic || c_.model == ThroughflowModel::SimpleIsentropic;
    }

    /// m, phi, kappa_m of every stream surface from the current positions (full model, p. 165:
    /// Eq. 7-58, then Eqs. 7-2 and 7-9 by Eq. 7-67 inside and Eq. 7-70 at the ends, p. 171).
    void surface_geometry(std::vector<Vec>& m, std::vector<Vec>& phi, std::vector<Vec>& kappa) const;
    /// The geometry of quasi-normal i at positions y for the model level.
    Geo geometry(int i, const Vec& y) const;
    /// One quasi-normal: the linking loop of Figs. 7-3 / 7-4.  Returns false when choked.
    bool solve_qn(int i, int pass, int& iterations, bool& converged);
    /// The flow spec of quasi-normal i (linking, Sec. 7.3) at the current estimate.
    Spec link(int i, const Geo& ge, const Vec& Wm, const Vec& rho_est);
    /// Eqs. 7-54 and the inversion to the mass fractions (Sec. 7.4).
    Vec reposition(const Geo& ge, const Vec& q) const;
    /// Mass-flow integrand K_W r rho W_m cos(eps) at the nodes (Eq. 7-5 / 7-54).
    Vec mass_integrand(const Geo& ge, const Spec& sp, const std::vector<NodeState>& st) const;
    MassEval mass_at(const Geo& ge, const Spec& sp, const LinearOdeSolution& ode, double w0, double KB,
                     std::vector<NodeState>& st) const;
    double solve_w0(const Geo& ge, const Spec& sp, const LinearOdeSolution& ode, double KB, double guess,
                    bool& choked, std::vector<NodeState>& st) const;
    void fill_station(int i, const Geo& ge, const Spec& sp, const std::vector<NodeState>& st, double KB,
                      int iterations, bool converged, bool choked, bool mach_limit);
    EndwallStation endwall_station(int i, int pass, const Geo& ge, const Spec& sp, const std::vector<NodeState>& st) const;
    /// The state of stream surface j as a wall state (m: its wall arc length).
    WallState wall_state(const Geo& ge, const Spec& sp, std::size_t j, const NodeState& x, double m) const;

    const ThroughflowCase& c_;
    const ThroughflowOptions& o_;
    int nqn_, nss_;
    std::shared_ptr<EndwallModel> endwall_;
    std::vector<QnData> q_;
    std::vector<int> row_at_, row_from_;
    Vec leak_;  ///< the latest seal leakage flow of every row (RowLossModel)
    Vec m_hub_, m_shroud_;
    Vec wall_phi_[2], wall_kappa_[2];  ///< hub, shroud: given or from the wall points (QuasiNormal)
    // full model: stream-surface geometry of the current pass, [j][i]
    std::vector<Vec> sm_, sphi_, skappa_, kappa_prev_;
};

void Solver::surface_geometry(std::vector<Vec>& m, std::vector<Vec>& phi, std::vector<Vec>& kappa) const
{
    m.assign(us(nss_), Vec(us(nqn_)));
    phi = m;
    kappa = m;
    for (int j = 0; j < nss_; ++j) {
        Vec z(us(nqn_)), r(us(nqn_));
        for (int i = 0; i < nqn_; ++i) place(q_[us(i)].g, q_[us(i)].y[us(j)], z[us(i)], r[us(i)]);
        m[us(j)] = meridional_coordinate(z, r);                      // Eq. 7-58
        const Vec drdm = derivative(m[us(j)], r, EndRule::TwoPoint);  // Eq. 7-2 by Eqs. 7-67, 7-70
        for (int i = 0; i < nqn_; ++i) {
            const double sp = drdm[us(i)];
            if (!(std::abs(sp) < 1.0)) fail("Eq. 7-2: |dr/dm| >= 1 on stream surface " + std::to_string(j));
            phi[us(j)][us(i)] = std::asin(sp);
        }
        const Vec dphidm = derivative(m[us(j)], phi[us(j)], EndRule::TwoPoint);  // Eq. 7-9 by Eqs. 7-67, 7-70
        for (int i = 0; i < nqn_; ++i) kappa[us(j)][us(i)] = -dphidm[us(i)];
    }
}

Geo Solver::geometry(int i, const Vec& y) const
{
    const auto& g = q_[us(i)].g;
    Geo ge;
    ge.y = y;
    ge.cosl = g.cosl;
    const std::size_t n = y.size();
    ge.z.resize(n);
    ge.r.resize(n);
    ge.phi.assign(n, 0.0);
    ge.eps.assign(n, 0.0);
    ge.phi_eq.assign(n, 0.0);
    ge.kappa_eq.assign(n, 0.0);
    ge.eps_eq.assign(n, 0.0);
    ge.dphidy.assign(n, 0.0);
    ge.m.assign(n, 0.0);
    for (std::size_t j = 1; j < n; ++j)
        if (!(y[j] > y[j - 1])) fail("stream surfaces cross on quasi-normal " + std::to_string(i));
    for (std::size_t j = 0; j < n; ++j) place(g, y[j], ge.z[j], ge.r[j]);
    switch (c_.model) {
    case ThroughflowModel::FullNormalEquilibrium:
        for (std::size_t j = 0; j < n; ++j) {
            ge.phi[j] = sphi_[j][us(i)];
            ge.kappa_eq[j] = skappa_[j][us(i)];
            ge.m[j] = sm_[j][us(i)];
        }
        ge.dphidy = derivative(ge.y, ge.phi, EndRule::ThreePoint);  // p. 172: Eqs. 7-67 to 7-69 in y
        break;
    case ThroughflowModel::Approximate: {  // p. 168: phi and kappa_m linear in y from hub to shroud
        const double ph = wall_phi_[0][us(i)], ps = wall_phi_[1][us(i)];
        const double kh = wall_kappa_[0][us(i)], ks = wall_kappa_[1][us(i)];
        for (std::size_t j = 0; j < n; ++j) {
            const double t = y[j] / g.ys;
            ge.phi[j] = ph + t * (ps - ph);
            ge.kappa_eq[j] = kh + t * (ks - kh);
            ge.dphidy[j] = (ps - ph) / g.ys;
        }
        break;
    }
    case ThroughflowModel::SimpleNonIsentropic:
    case ThroughflowModel::SimpleIsentropic: {
        // epsilon = kappa_m = 0 in the normal-momentum equation only (p. 168); the mass equations keep
        // the geometric phi, from Eq. 7-2 by Eq. 7-70 on the marching pair (i - 1, i) (p. 172: the data
        // downstream are not known when quasi-normal i is solved).  Quasi-normal 0 has no upstream
        // pair: there phi is linear in y between the end-wall angles (ours, as the Approximate model)
        if (i == 0) {
            for (std::size_t j = 0; j < n; ++j)
                ge.phi[j] = wall_phi_[0][0] + y[j] / g.ys * (wall_phi_[1][0] - wall_phi_[0][0]);
            break;
        }
        const auto& go = q_[us(i - 1)].g;
        const Vec& yo = q_[us(i - 1)].y;
        for (std::size_t j = 0; j < n; ++j) {
            double zo = 0.0, ro = 0.0;
            place(go, yo[j], zo, ro);
            ge.phi[j] = std::atan2(ge.r[j] - ro, ge.z[j] - zo);  // sin(phi) = dr/dm, dm the chord of the pair (Eq. 7-70)
        }
        break;
    }
    }
    for (std::size_t j = 0; j < n; ++j) {
        ge.eps[j] = ge.phi[j] - g.lambda;  // Eq. 7-3
        if (!simple()) {
            ge.phi_eq[j] = ge.phi[j];
            ge.eps_eq[j] = ge.eps[j];
        }
    }
    return ge;
}

Spec Solver::link(int i, const Geo& ge, const Vec& Wm, const Vec& rho_est)
{
    const auto& gas = c_.gas;
    const std::size_t n = us(nss_);
    Spec sp;
    sp.KW.assign(n, 1.0);
    sp.KBlade.assign(n, 1.0);
    sp.loss.assign(n, 0.0);
    sp.row_loss.assign(n, 0.0);
    sp.I.resize(n);
    sp.s.resize(n);
    if (i == 0) {  // inlet boundary conditions (p. 162), interpolated at y / y_s (ours)
        const auto& in = c_.inlet;
        sp.omega = 0.0;
        sp.known_beta = in.swirl == InletProfile::Swirl::Alpha;
        sp.rCt.assign(n, 0.0);
        sp.tanb.assign(n, 0.0);
        for (std::size_t j = 0; j < n; ++j) {
            const double t = std::clamp(ge.y[j] / q_[0].g.ys, 0.0, 1.0);
            const double Pt = interp(in.span, in.Pt, t), Tt = interp(in.span, in.Tt, t);
            const double v = interp(in.span, in.swirl_value, t);
            sp.I[j] = gas.cp() * Tt;  // I = H in the absolute frame
            sp.s[j] = gas.entropy(Tt, Pt);
            if (sp.known_beta) sp.tanb[j] = std::tan(v);  // p. 162: C, beta, H for W, beta', I with omega = 0
            else sp.rCt[j] = ge.r[j] * v;
        }
        return sp;
    }
    const QnData& up = q_[us(i - 1)];
    const int k = row_at_[us(i)];
    if (k < 0) {  // no blade row: Eqs. 7-40 to 7-42, (rC_theta)_2 = (rC_theta)_1, H_2 = H_1, s_2 = s_1
        sp.omega = 0.0;
        sp.known_beta = false;
        sp.rCt = up.rCt;
        sp.tanb.assign(n, 0.0);
        sp.I = up.H;
        sp.s = up.s;
        return sp;
    }
    // a blade row between i - 1 and i: Eqs. 7-43 to 7-49 with the row model (Sec. 7.3)
    const BladeRowSpec& row = c_.rows[us(k)];
    const double w = row.omega;
    sp.omega = w;
    sp.known_beta = true;
    sp.rCt.assign(n, 0.0);
    sp.tanb.resize(n);
    Vec secspan;
    for (const auto& s : row.sections) secspan.push_back(s.span);
    auto section_at = [&](double span) {
        BladeSectionGeometry gsec;
        auto pick = [&](auto field) {
            Vec v;
            for (const auto& s : row.sections) v.push_back(s.*field);
            return interp(secspan, v, span);
        };
        gsec.span = span;
        gsec.chord = pick(&BladeSectionGeometry::chord);
        gsec.stagger = pick(&BladeSectionGeometry::stagger);
        gsec.camber = pick(&BladeSectionGeometry::camber);
        gsec.kappa1 = pick(&BladeSectionGeometry::kappa1);
        gsec.kappa2 = pick(&BladeSectionGeometry::kappa2);
        gsec.thickness_ratio = pick(&BladeSectionGeometry::thickness_ratio);
        gsec.max_camber_position = pick(&BladeSectionGeometry::max_camber_position);
        gsec.max_thickness_position = pick(&BladeSectionGeometry::max_thickness_position);
        return gsec;
    };
    Vec upz(n), upr(n), Hrel2(n), s1v(n), qref(n);
    for (std::size_t j = 0; j < n; ++j) place(up.g, up.y[j], upz[j], upr[j]);
    for (std::size_t j = 0; j < n; ++j) {
        BladeRowInlet in;
        in.row = k;
        in.stream_surface = ii(j);
        in.n_stream_surfaces = nss_;
        in.mass_fraction = static_cast<double>(j) / static_cast<double>(nss_ - 1);
        in.omega = w;
        in.n_blades = row.n_blades;
        in.z1 = upz[j];
        in.r1 = upr[j];
        in.z2 = ge.z[j];
        in.r2 = ge.r[j];
        if (j == 0) in.span = 0.0;  // the wall stream surfaces are the walls (ours)
        else if (j + 1 == n) in.span = 1.0;
        else {
            const double sle = line_crossing_fraction(row.le_hub_z, row.le_hub_r, row.le_shroud_z, row.le_shroud_r, in.z1,
                                                      in.r1, in.z2, in.r2);
            const double ste = line_crossing_fraction(row.te_hub_z, row.te_hub_r, row.te_shroud_z, row.te_shroud_r, in.z1,
                                                      in.r1, in.z2, in.r2);
            in.span = 0.5 * (sle + ste);  // ours (Exercise 7.1)
            if (!(in.span > 0.0 && in.span < 1.0))
                fail(row.name + ": stream surface " + std::to_string(j) + " leaves the blade span");
        }
        in.geometry = section_at(in.span);
        in.pitch = 2.0 * kPi * 0.5 * (in.r1 + in.r2) / static_cast<double>(row.n_blades);
        in.solidity = in.geometry.chord / in.pitch;
        // inlet relative conditions from the upstream absolute flow (Eqs. 7-43 to 7-46)
        const double H1 = up.H[j], rCt1 = up.rCt[j], s1 = up.s[j], Wm1 = up.Wm[j];
        const double Ct1 = up.Ct[j];
        in.I1 = H1 - w * rCt1;                       // Eq. 7-43
        in.Ctheta1 = Ct1;
        in.Wtheta1 = Ct1 - in.r1 * w;                // Eq. 7-44
        in.Wm1 = Wm1;
        in.W1 = std::hypot(Wm1, in.Wtheta1);         // Eq. 7-45
        const double h1 = gas.cp() * up.T[j];
        const double Hrel1 = h1 + 0.5 * in.W1 * in.W1;  // Eq. 7-46
        in.T1 = up.T[j];
        in.P1 = up.P[j];
        in.rho1 = up.rho[j];
        in.a1 = up.a[j];
        in.s1 = s1;
        in.Tt1_rel = Hrel1 / gas.cp();
        in.Pt1_rel = gas.pressure(Hrel1, s1);
        in.beta1 = std::atan2(in.Wtheta1, Wm1);
        in.alpha1 = std::atan2(Ct1, Wm1);
        in.mach1_rel = in.W1 / in.a1;
        in.mach1_abs = std::hypot(Wm1, Ct1) / in.a1;
        in.reynolds = in.rho1 * in.W1 * in.geometry.chord / gas.mu(in.T1);
        in.h1 = h1;
        in.Hrel1 = Hrel1;
        in.mu1 = gas.mu(in.T1);
        in.s_row = row.s_row;
        in.is_igv = row.is_igv;
        in.profile_family = row.profile_family;
        in.blade_height = blade_height(row);
        in.row_number = k + 1;
        in.mass_flow = c_.mass_flow;
        in.Wm2 = Wm[j];
        in.rho2 = rho_est[j];
        const BladeRowDischarge d = row.model->evaluate(in);
        if (!finite(d.beta2) || !(std::abs(d.beta2) < 0.5 * kPi)) fail(row.name + ": the model returned |beta2| >= pi/2");
        if (!finite(d.loss_coefficient)) fail(row.name + ": the model returned a non-finite loss");
        if (!(d.wake_blockage > 0.0 && d.wake_blockage <= 1.0)) fail(row.name + ": the model returned K_W outside (0, 1]");
        // discharge: Eqs. 7-47, 7-48
        sp.I[j] = in.I1;                                                   // Eq. 7-47
        Hrel2[j] = sp.I[j] + 0.5 * (ge.r[j] * w) * (ge.r[j] * w);          // Eq. 7-48
        s1v[j] = s1;
        qref[j] = in.Pt1_rel - in.P1;  // the reference pressure of Eq. 7-49 on this surface
        sp.tanb[j] = std::tan(d.beta2);
        sp.KW[j] = d.wake_blockage;
        sp.loss[j] = d.loss_coefficient;
    }
    // Eq. 7-49 with the per-surface losses: P't2 = P't2,id - omegabar (P't1 - P1), then s2 from (H'2, P't2)
    auto discharge_entropy = [&]() {
        for (std::size_t j = 0; j < n; ++j) {
            const double Pt2id = gas.pressure(Hrel2[j], s1v[j]);  // ideal: inlet entropy (p. 163)
            const double Pt2 = Pt2id - sp.loss[j] * qref[j];      // Eq. 7-49
            if (!(Pt2 > 0.0)) fail(row.name + ": the loss gives a non-positive discharge total pressure");
            sp.s[j] = gas.entropy(Hrel2[j] / gas.cp(), Pt2);
        }
    };
    discharge_entropy();
    sp.row_loss.assign(n, 0.0);
    if (row.row_loss) {
        // the row-level losses (Ch. 6 Secs. 6.9, 6.10): wall states at both quasi-normals, the discharge at
        // the current estimate with the per-surface losses only (ours)
        RowLossInput rl;
        rl.row = k;
        rl.omega = w;
        rl.n_blades = row.n_blades;
        rl.s_row = row.s_row;
        rl.shrouded = row.shrouded;
        rl.hub_clearance = row.hub_clearance;
        rl.tip_clearance = row.tip_clearance;
        rl.mass_flow = c_.mass_flow;
        rl.blade_height = blade_height(row);
        rl.hub_section = row.sections.front();
        rl.shroud_section = row.sections.back();
        rl.row_number = k + 1;
        rl.hub1 = up.wall[0];
        rl.shroud1 = up.wall[1];
        for (WallState* w1 : {&rl.hub1, &rl.shroud1}) w1->Wtheta = w1->Ctheta - w1->r * w;  // the row's frame (Eq. 7-44)
        NodeState x;
        for (int side = 0; side < 2; ++side) {
            const std::size_t j = side == 0 ? 0 : n - 1;
            if (!node_state(gas, sp, ge, j, Wm[j], x)) fail(row.name + ": invalid discharge estimate at a wall");
            WallState& ws = side == 0 ? rl.hub2 : rl.shroud2;
            ws = wall_state(ge, sp, j, x, side == 0 ? m_hub_[us(i)] : m_shroud_[us(i)]);
        }
        const RowLossResult lr = row.row_loss->evaluate(rl);
        if (!finite(lr.clearance_dPt) || !finite(lr.seal_dPt) || !finite(lr.leakage_mass_flow))
            fail(row.name + ": the row-level loss model returned a non-finite value");
        leak_[us(k)] = lr.leakage_mass_flow;
        Vec psi(n);
        for (std::size_t j = 0; j < n; ++j) psi[j] = static_cast<double>(j) / static_cast<double>(n - 1);
        const Vec a = distribute_row_loss(lr.clearance_dPt, lr.clearance_side, psi);
        const Vec b = distribute_row_loss(lr.seal_dPt, lr.seal_side, psi);
        for (std::size_t j = 0; j < n; ++j) {
            sp.row_loss[j] = (a[j] + b[j]) / qref[j];  // Delta P_t,j / (P't1 - P1)_j (ours)
            sp.loss[j] += sp.row_loss[j];
        }
        discharge_entropy();
    }
    return sp;
}

Vec Solver::mass_integrand(const Geo& ge, const Spec& sp, const std::vector<NodeState>& st) const
{
    Vec q(ge.y.size());
    for (std::size_t j = 0; j < q.size(); ++j) q[j] = sp.KW[j] * sp.KBlade[j] * ge.r[j] * st[j].g.rho * st[j].Wm * std::cos(ge.eps[j]);
    return q;
}

MassEval Solver::mass_at(const Geo& ge, const Spec& sp, const LinearOdeSolution& ode, double w0, double KB,
                         std::vector<NodeState>& st) const
{
    MassEval e;
    const std::size_t n = ge.y.size();
    const Vec W = ode.evaluate(w0);
    Vec mm(n);
    for (std::size_t j = 0; j < n; ++j) {
        if (!node_state(c_.gas, sp, ge, j, W[j], st[j])) return e;
        mm[j] = W[j] / st[j].g.a;  // Eq. 7-27
    }
    e.valid = true;
    e.mdot = 2.0 * kPi * KB * integral(ge.y, mass_integrand(ge, sp, st));  // Eq. 7-5
    e.mean_mm = integral(ge.y, mm) / (ge.y.back() - ge.y.front());       // Eq. 7-39
    return e;
}

double Solver::solve_w0(const Geo& ge, const Spec& sp, const LinearOdeSolution& ode, double KB, double guess,
                        bool& choked, std::vector<NodeState>& st) const
{
    // W_m(0), the constant of integration of Eq. 7-35, from Eq. 7-5 (p. 161), limited by the choke
    // criterion Eq. 7-39 (p. 161: "easily employed to limit the value of W_m(0)") or Eq. 12-37.  Ours: the
    // search starts where the smallest W_m is 0 (the stream-tube form Eq. 7-33 allows W_m = 0, a
    // reverse flow is refused), brackets the target upward, and bisects (Illinois).
    choked = false;
    const double target = c_.mass_flow;
    const std::size_t n = ge.y.size();
    double lo = -std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < n; ++j) lo = std::max(lo, -ode.G[j] / ode.F[j]);
    MassEval elo = mass_at(ge, sp, ode, lo, KB, st);
    if (!elo.valid) fail("no valid state on a quasi-normal even at the smallest meridional velocity");
    if (elo.mdot > target * (1.0 + o_.mass_tolerance)) fail("the mass flow needs reverse flow on a quasi-normal");
    // below the choke limit: Eq. 7-39 (mean M_m < 1), or Eq. 12-37 (d mdot / d W_m(0) > 0, by a
    // forward difference over 1e-7 max(1, |W_m(0)|), ours)
    std::vector<NodeState> probe(st.size());
    auto ok = [&](double x, const MassEval& e) {
        if (!e.valid) return false;
        if (o_.choke_test == ChokeTest::MeanMeridionalMach) return e.mean_mm < 1.0;
        const MassEval e2 = mass_at(ge, sp, ode, x + 1e-7 * std::max(1.0, std::abs(x)), KB, probe);
        return e2.valid && e2.mdot > e.mdot;
    };
    double step = std::max(std::abs(guess - lo), 1.0);
    double a = lo, hi = lo;
    MassEval ea = elo, ehi = elo;
    bool bracket = false;
    for (int it = 0; it < 200; ++it) {
        hi = a + step;
        ehi = mass_at(ge, sp, ode, hi, KB, st);
        if (!ok(hi, ehi)) break;
        if (ehi.mdot >= target) {
            bracket = true;
            break;
        }
        a = hi;
        ea = ehi;
        step *= 2.0;
    }
    if (!bracket) {
        if (ok(hi, ehi)) fail("the meridional velocity search did not bracket the mass flow");
        // the choke limit between a (ok) and hi (not ok): bisection on Eq. 7-39 / validity
        double b = hi;
        for (int it = 0; it < 200 && (b - a) > 1e-14 * std::max(1.0, std::abs(b)); ++it) {
            const double mid = 0.5 * (a + b);
            const MassEval em = mass_at(ge, sp, ode, mid, KB, st);
            if (ok(mid, em)) {
                a = mid;
                ea = em;
            } else b = mid;
        }
        if (ea.mdot < target * (1.0 - o_.mass_tolerance)) {
            choked = true;
            (void)mass_at(ge, sp, ode, a, KB, st);
            return a;
        }
        hi = a;
        ehi = ea;
        a = lo;
        ea = elo;
        if (ehi.mdot < target) {  // within the tolerance below the target at the limit
            (void)mass_at(ge, sp, ode, hi, KB, st);
            return hi;
        }
    }
    // Illinois regula falsi on [a, hi]: mdot(a) < target <= mdot(hi)
    double fa = ea.mdot - target, fb = ehi.mdot - target, b = hi;
    int side = 0;
    double x = b;
    for (int it = 0; it < 300; ++it) {
        x = (a * fb - b * fa) / (fb - fa);
        if (!(x > a && x < b)) x = 0.5 * (a + b);
        const MassEval ex = mass_at(ge, sp, ode, x, KB, st);
        if (!ex.valid) fail("invalid state inside the W_m(0) bracket");
        const double fx = ex.mdot - target;
        if (std::abs(fx) <= o_.mass_tolerance * target || (b - a) <= 1e-15 * std::max(1.0, std::abs(b))) break;
        if ((fx < 0.0) == (fa < 0.0)) {
            a = x;
            fa = fx;
            if (side == -1) fb *= 0.5;
            side = -1;
        } else {
            b = x;
            fb = fx;
            if (side == +1) fa *= 0.5;
            side = +1;
        }
    }
    (void)mass_at(ge, sp, ode, x, KB, st);
    return x;
}

Vec Solver::reposition(const Geo& ge, const Vec& q) const
{
    // Eq. 7-54 (K_B omitted, p. 164): mdot(y) = 2 pi int_0^y K_W r rho W_m cos(eps) dy; the stream
    // surfaces go to y with mdot(y) = psi_j mdot(y_s), "by interpolation from this function": ours,
    // the three-point series of Eqs. 7-72 / 7-73 integrated to y inside the segment (so the
    // interpolant is consistent with the integral itself), solved by bisection.
    const Vec& y = ge.y;
    const std::size_t n = y.size();
    const Vec M = cumulative_integral(y, q);
    for (std::size_t k = 0; k + 1 < n; ++k)
        if (!(M[k + 1] > M[k])) fail("the mass-flow function Eq. 7-54 does not increase across a stream tube");
    Vec ynew(n);
    ynew.front() = y.front();
    ynew.back() = y.back();
    for (std::size_t j = 1; j + 1 < n; ++j) {
        const double t = static_cast<double>(j) / static_cast<double>(n - 1) * M.back();
        std::size_t k = static_cast<std::size_t>(std::upper_bound(M.begin(), M.end(), t) - M.begin());
        k = std::min(std::max<std::size_t>(k, 1), n - 1) - 1;
        double a = y[k], b = y[k + 1];
        const double need = t - M[k];
        for (int it = 0; it < 200 && (b - a) > 1e-15 * y.back(); ++it) {
            const double mid = 0.5 * (a + b);
            if (partial_segment_integral(y, q, ii(k), mid) < need) a = mid;
            else b = mid;
        }
        ynew[j] = 0.5 * (a + b);
    }
    return ynew;
}

WallState Solver::wall_state(const Geo& ge, const Spec& sp, std::size_t j, const NodeState& x, double m) const
{
    WallState w;
    w.z = ge.z[j];
    w.r = ge.r[j];
    w.m = m;
    w.phi = ge.phi[j];
    w.Cm = x.Wm;
    w.Ctheta = x.Ct;
    w.Wtheta = x.Wt;
    w.rho = x.g.rho;
    w.P = x.g.P;
    w.T = x.g.T;
    w.a = x.g.a;
    w.mu = c_.gas.mu(x.g.T);
    const double H = sp.I[j] + sp.omega * ge.r[j] * x.Ct;  // Eq. 7-51
    w.Tt = H / c_.gas.cp();
    w.Pt = c_.gas.pressure(H, sp.s[j]);
    return w;
}

EndwallStation Solver::endwall_station(int i, int pass, const Geo& ge, const Spec& sp,
                                       const std::vector<NodeState>& st) const
{
    EndwallStation e;
    e.qn = i;
    e.n_qn = nqn_;
    e.pass = pass;
    e.mass_flow = c_.mass_flow;
    e.compressor_omega = c_.compressor_omega;
    e.y_s = q_[us(i)].g.ys;
    e.lambda = q_[us(i)].g.lambda;
    e.frame_omega = sp.omega;
    e.row_upstream = row_at_[us(i)];
    e.row_downstream = row_from_[us(i)];
    const std::size_t last = ge.y.size() - 1;
    const auto& qn = c_.quasi_normals[us(i)];
    auto edge = [&](std::size_t j, double m, double wall_omega) {
        EndwallEdge w;
        static_cast<WallState&>(w) = wall_state(ge, sp, j, st[j], m);
        w.wall_omega = wall_omega;
        w.Vtheta = st[j].Ct - ge.r[j] * wall_omega;  // Eq. 8-4
        return w;
    };
    e.hub = edge(0, m_hub_[us(i)], qn.hub_wall_omega);
    e.shroud = edge(last, m_shroud_[us(i)], qn.shroud_wall_omega);
    // the row data: for the row upstream the boundary-layer step across the row ends here
    auto row_data = [&](int k, bool discharge) {
        EndwallRow w;
        if (k < 0) return w;
        const BladeRowSpec& b = c_.rows[us(k)];
        w.row = k;
        w.row_number = k + 1;
        w.omega = b.omega;
        w.n_blades = b.n_blades;
        w.shrouded = b.shrouded;
        w.leakage_mass_flow = leak_[us(k)];
        const double Z = static_cast<double>(b.n_blades);
        for (int side = 0; side < 2; ++side) {
            EndwallRowWall& ww = side == 0 ? w.hub : w.shroud;
            const std::size_t j = side == 0 ? 0 : last;
            const Vec& mw = side == 0 ? m_hub_ : m_shroud_;
            ww.section = side == 0 ? b.sections.front() : b.sections.back();
            ww.stagger_mean = ww.section.stagger;
            ww.clearance = side == 0 ? b.hub_clearance : b.tip_clearance;
            if (discharge) {
                const WallState& w1 = q_[us(i - 1)].wall[side];
                ww.r_inlet = w1.r;
                ww.r_discharge = ge.r[j];
                ww.Wtheta_inlet = w1.Ctheta - w1.r * b.omega;  // Eq. 8-67 in the row frame
                ww.Wtheta_discharge = st[j].Ct - ge.r[j] * b.omega;
                ww.delta_m = mw[us(i)] - mw[us(i - 1)];
                ww.delta_rCtheta = ge.r[j] * st[j].Ct - w1.r * w1.Ctheta;
                ww.pitch_mean = 2.0 * kPi * 0.5 * (ww.r_inlet + ww.r_discharge) / Z;
            } else {
                ww.r_inlet = ge.r[j];
                ww.Wtheta_inlet = st[j].Ct - ge.r[j] * b.omega;
                ww.pitch_mean = 2.0 * kPi * ge.r[j] / Z;  // only the inlet is known here
            }
        }
        return w;
    };
    e.upstream = row_data(e.row_upstream, true);
    e.downstream = row_data(e.row_downstream, false);
    return e;
}

void Solver::fill_station(int i, const Geo& ge, const Spec& sp, const std::vector<NodeState>& st, double KB,
                          int iterations, bool converged, bool choked, bool mach_limit)
{
    QnData& d = q_[us(i)];
    const auto& gas = c_.gas;
    const std::size_t n = ge.y.size();
    for (std::size_t j = 0; j < n; ++j)
        if (!(st[j].Wm > 0.0))
            fail("quasi-normal " + std::to_string(i) + ": W_m <= 0 on stream surface " + std::to_string(j));
    ThroughflowStation& s = d.out;
    s = ThroughflowStation{};
    s.index = i;
    s.lambda = d.g.lambda;
    s.y_s = d.g.ys;
    s.frame_omega = sp.omega;
    s.row_upstream = row_at_[us(i)];
    s.K_B_endwall = KB;
    s.choked = choked;
    s.mach_limit_active = mach_limit;
    s.converged = converged;
    s.iterations = iterations;
    // the sheets' y-widths carried to the surfaces (CONVENTIONS)
    Vec width(n - 1);
    for (std::size_t k = 0; k + 1 < n; ++k) width[k] = ge.y[k + 1] - ge.y[k];
    const Vec D = sheets_to_surfaces(ge.y, width);
    const double dpsi = 1.0 / static_cast<double>(n - 1);
    d.Wm.resize(n);
    d.H.resize(n);
    d.rCt.resize(n);
    d.s.resize(n);
    d.T.resize(n);
    d.P.resize(n);
    d.rho.resize(n);
    d.a.resize(n);
    d.Ct.resize(n);
    Vec mm(n);
    for (std::size_t j = 0; j < n; ++j) {
        const NodeState& x = st[j];
        ThroughflowNode nd;
        nd.z = ge.z[j];
        nd.r = ge.r[j];
        nd.y = ge.y[j];
        nd.m = ge.m[j];
        nd.mass_fraction = static_cast<double>(j) / static_cast<double>(n - 1);
        nd.phi = ge.phi[j];
        nd.kappa_m = ge.kappa_eq[j];
        nd.epsilon = ge.eps[j];
        nd.Wm = x.Wm;
        nd.Ctheta = x.Ct;  // Eq. 7-50 on a rotor discharge
        nd.Wtheta = x.Wt;
        nd.C = std::hypot(x.Wm, x.Ct);  // Eqs. 7-52, 7-53
        nd.W = std::hypot(x.Wm, x.Wt);
        nd.alpha = std::atan2(x.Ct, x.Wm);
        nd.beta = std::atan2(x.Wt, x.Wm);
        nd.P = x.g.P;
        nd.T = x.g.T;
        nd.rho = x.g.rho;
        nd.a = x.g.a;
        nd.h = x.h;
        nd.mu = gas.mu(x.g.T);
        nd.mach = nd.C / x.g.a;
        nd.mach_rel = nd.W / x.g.a;
        nd.mach_m = x.Wm / x.g.a;
        nd.mach_theta = x.Ct / x.g.a;
        nd.s = sp.s[j];
        nd.I = sp.I[j];
        nd.Hrel = x.Hrel;
        nd.H = sp.I[j] + sp.omega * ge.r[j] * x.Ct;  // Eq. 7-51
        nd.rCtheta = ge.r[j] * x.Ct;
        nd.Tt = nd.H / gas.cp();
        nd.Pt = gas.pressure(nd.H, sp.s[j]);
        nd.Tt_rel = x.Hrel / gas.cp();
        nd.Pt_rel = gas.pressure(x.Hrel, sp.s[j]);
        nd.K_W = sp.KW[j];
        nd.K_B_blade = sp.KBlade[j];
        nd.b = std::cos(ge.eps[j]) * D[j];  // the area term of Eq. 7-32 for a sheet of dpsi at the surface
        nd.b_mass = dpsi * c_.mass_flow /
                    (2.0 * kPi * ge.r[j] * KB * sp.KW[j] * sp.KBlade[j] * x.g.rho * x.Wm);  // Eq. 12-39
        nd.loss_coefficient = sp.loss[j];
        nd.row_loss_coefficient = sp.row_loss[j];
        s.nodes.push_back(nd);
        mm[j] = nd.mach_m;
        d.Wm[j] = x.Wm;
        d.H[j] = nd.H;
        d.rCt[j] = nd.rCtheta;
        d.s[j] = nd.s;
        d.T[j] = nd.T;
        d.P[j] = nd.P;
        d.rho[j] = nd.rho;
        d.a[j] = nd.a;
        d.Ct[j] = nd.Ctheta;
    }
    // the stream sheets (j, j + 1) at their mean surfaces (CONVENTIONS)
    const double dmdot_sheet = c_.mass_flow / static_cast<double>(n - 1);
    for (std::size_t j = 0; j + 1 < n; ++j) {
        ThroughflowSheetNode sh;
        sh.z = 0.5 * (ge.z[j] + ge.z[j + 1]);
        sh.r = 0.5 * (ge.r[j] + ge.r[j + 1]);
        sh.b_geometric = std::cos(0.5 * (ge.eps[j] + ge.eps[j + 1])) * (ge.y[j + 1] - ge.y[j]);
        const double rrw = 0.5 * (ge.r[j] * st[j].g.rho * st[j].Wm + ge.r[j + 1] * st[j + 1].g.rho * st[j + 1].Wm);
        sh.b = dmdot_sheet / (2.0 * kPi * rrw);
        s.sheets.push_back(sh);
    }
    d.wall[0] = wall_state(ge, sp, 0, st[0], m_hub_[us(i)]);
    d.wall[1] = wall_state(ge, sp, n - 1, st[n - 1], m_shroud_[us(i)]);
    s.mass_flow = 2.0 * kPi * KB * integral(ge.y, mass_integrand(ge, sp, st));  // Eq. 7-5
    s.mass_error = (s.mass_flow - c_.mass_flow) / c_.mass_flow;
    s.mean_mach_m = integral(ge.y, mm) / d.g.ys;  // Eq. 7-39
    d.KB = KB;
    d.solved = true;
}

bool Solver::solve_qn(int i, int pass, int& iterations, bool& converged)
{
    QnData& d = q_[us(i)];
    const auto& gas = c_.gas;
    const std::size_t n = us(nss_);
    Vec y = d.y;
    Geo ge = geometry(i, y);
    // the first estimate: this quasi-normal's last solution, else the upstream one (W_m2 = W_m1,
    // p. 163), else uniform W_m = mdot / (rho_t A) with A of Eq. 7-60 (ours, as Eq. 7-61)
    Vec Wm(n), rho(n);
    double KB = 1.0;
    if (d.solved) {
        Wm = d.Wm;
        rho = d.rho;
        KB = d.KB;
    } else if (i > 0) {
        Wm = q_[us(i - 1)].Wm;
        rho = q_[us(i - 1)].rho;
        KB = q_[us(i - 1)].KB;
    } else {
        const double Pt = c_.inlet.Pt.front(), Tt = c_.inlet.Tt.front();
        const double rho_t = Pt / (gas.gas_constant * Tt);
        const double A = kPi * (d.g.rs + d.g.rh) * d.g.ys;
        // ours: an estimate beyond half the limiting speed sqrt(2 H) is replaced by that half (a first
        // iterate only; the mass-flow iteration and the choke test decide)
        Wm.assign(n, std::min(c_.mass_flow / (rho_t * A), 0.5 * std::sqrt(2.0 * gas.cp() * Tt)));
        rho.assign(n, rho_t);
    }
    std::vector<NodeState> st(n);
    bool choked = false;
    converged = false;
    Vec f1(n), f4(n);
    for (iterations = 1; iterations <= o_.max_flow_iterations; ++iterations) {
        // Link To Upstream Quasi-Normal (Sec. 7.3), at the current discharge estimate
        const Spec sp = link(i, ge, Wm, rho);
        // the state at the current estimate (the coefficients of Eq. 7-15)
        for (std::size_t j = 0; j < n; ++j)
            if (!node_state(gas, sp, ge, j, Wm[j], st[j]))
                fail("quasi-normal " + std::to_string(i) + ": the estimate gives h <= 0");
        const Vec dIdy = derivative(ge.y, sp.I, EndRule::ThreePoint);  // p. 172: Eqs. 7-67 to 7-69 in y
        const Vec dsdy = derivative(ge.y, sp.s, EndRule::ThreePoint);
        Vec rtb(n);
        for (std::size_t j = 0; j < n; ++j) rtb[j] = sp.known_beta ? ge.r[j] * sp.tanb[j] : sp.rCt[j];
        const Vec dXdy = derivative(ge.y, rtb, EndRule::ThreePoint);  // d(r tan beta')/dy or d(r C_theta)/dy
        // Eqs. 7-31 to 7-33 on the stream tubes between the surfaces (the sheets): sheet k carries the
        // share dpsi = 1 / (n - 1) of the K_B-free mass function of Eq. 7-54, dmdot = dpsi mdot / K_B
        // (ours: Eq. 7-32 omits K_B, Eq. 7-5 has it), and dA = 2 pi r K_W cos(eps) dy (Eq. 7-32) with the
        // mean of the sheet's two surfaces; its rho dA / dmdot (= 1 / W_m of the sheet, Eq. 7-31) is
        // carried to the surfaces at second order (ours)
        Vec g_tube(n, 0.0);
        if (o_.stream_tube_form) {
            const double dm = c_.mass_flow / (static_cast<double>(n - 1) * KB);
            Vec gk(n - 1);
            for (std::size_t k = 0; k + 1 < n; ++k) {
                auto area = [&](std::size_t j) {  // rho 2 pi r K_W K_B_blade cos(eps), per unit dy
                    return st[j].g.rho * 2.0 * kPi * ge.r[j] * sp.KW[j] * sp.KBlade[j] * std::cos(ge.eps[j]);
                };
                const double dA_rho = 0.5 * (area(k) + area(k + 1)) * (ge.y[k + 1] - ge.y[k]);
                gk[k] = dA_rho / dm;
            }
            g_tube = sheets_to_surfaces(ge.y, gk);
        }
        bool mach_limit = false;
        for (std::size_t j = 0; j < n; ++j) {
            const NodeState& x = st[j];
            // the normal-momentum equation's epsilon and kappa_m (0 in the simple models, p. 168)
            const double r = ge.r[j], ce = std::cos(ge.eps_eq[j]), se = std::sin(ge.eps_eq[j]);
            const double kap = ge.kappa_eq[j];
            // (1 / W_m) dW_m/dm from Eq. 7-29 with the limit Eq. 7-30 (Mtheta = C_theta / a, Eq. 7-28)
            const double Mm = x.Wm / x.g.a, Mt = x.Ct / x.g.a;
            const double lim = std::max(1.0 - Mm * Mm, o_.min_one_minus_mm2);
            if (1.0 - Mm * Mm < o_.min_one_minus_mm2 && se != 0.0) mach_limit = true;
            const double dlnW =
                (-(1.0 + Mt * Mt) * std::sin(ge.phi_eq[j]) / r - ge.dphidy[j] / ce - kap * std::tan(ge.eps_eq[j])) / lim;
            const double entropy_term = c_.model == ThroughflowModel::SimpleIsentropic ? 0.0 : x.g.T * dsdy[j];
            double f2 = 0.0, f3 = 0.0;
            if (sp.known_beta) {
                const double tb = sp.tanb[j], cb2 = 1.0 / (1.0 + tb * tb);  // cos^2 beta'
                f1[j] = cb2 * (-kap * ce - tb / r * dXdy[j] + se * dlnW);  // Eq. 7-19
                f2 = -2.0 * sp.omega * cb2 * tb * ge.cosl;  // Eq. 7-20: -2 omega cos(beta') sin(beta') cos(lambda)
                f3 = cb2 * (dIdy[j] - entropy_term);       // Eq. 7-21
            } else {
                f1[j] = -kap * ce + se * dlnW;                     // Eq. 7-16
                f2 = 0.0;                                          // Eq. 7-17
                f3 = dIdy[j] - entropy_term - x.Wt / r * dXdy[j];  // Eq. 7-18
            }
            // Eq. 7-33 (Aungier's choice, p. 160), or Eq. 7-15 as it stands: f3 / W_m
            f4[j] = o_.stream_tube_form ? f2 + f3 * g_tube[j] : f2 + f3 / x.Wm;
        }
        const LinearOdeSolution ode = solve_linear_ode(ge.y, f1, f4);  // Eqs. 7-35, 7-36
        bool ch = false;
        const double w0 = solve_w0(ge, sp, ode, KB, Wm.front(), ch, st);
        choked = ch;
        const Vec Wn = ode.evaluate(w0);
        double change = 0.0, scale = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            change = std::max(change, std::abs(Wn[j] - Wm[j]));
            scale = std::max(scale, std::abs(Wn[j]));
        }
        change /= scale;
        Wm = Wn;
        for (std::size_t j = 0; j < n; ++j) rho[j] = st[j].g.rho;
        // simplified models: reposition inside the quasi-normal loop, undamped (p. 168)
        double dpos = 0.0;
        if (!full()) {
            const Vec yn = reposition(ge, mass_integrand(ge, sp, st));
            for (std::size_t j = 0; j < n; ++j) dpos = std::max(dpos, std::abs(yn[j] - y[j]) / d.g.ys);
            y = yn;
        }
        // Compute End-Wall Boundary Layer Data (Chapter 8)
        const EndwallResult ew = endwall_->evaluate(endwall_station(i, pass, ge, sp, st));
        if (!(ew.blockage_factor > 0.0 && ew.blockage_factor <= 1.0)) fail("the end-wall model returned K_B outside (0, 1]");
        const double dkb = std::abs(ew.blockage_factor - KB);
        KB = ew.blockage_factor;
        const bool done = change <= o_.flow_tolerance && dpos <= o_.position_tolerance && dkb <= o_.flow_tolerance;
        if (done || choked) {
            // the reported solution is the one of this iteration, on the positions it was solved on
            // (the simplified models' last move is below the position tolerance)
            converged = done;
            d.y = ge.y;
            fill_station(i, ge, sp, st, KB, iterations, converged, choked, mach_limit);
            return !choked;
        }
        if (iterations == o_.max_flow_iterations) {
            d.y = ge.y;
            fill_station(i, ge, sp, st, KB, iterations, false, choked, mach_limit);
            return true;
        }
        if (!full()) ge = geometry(i, y);
    }
    return true;  // not reached: the last iteration returns
}

ThroughflowResult Solver::run()
{
    ThroughflowResult res;
    res.model = c_.model;
    BstarControl bc;
    bc.bstar = o_.bstar;
    bc.threshold = o_.bstar_growth_threshold;
    double bstar = o_.bstar;
    int repositionings = 0;
    const int passes = full() ? o_.max_outer_iterations : 1;
    for (int pass = 0; pass < passes; ++pass) {
        if (full()) {
            surface_geometry(sm_, sphi_, skappa_);
            // Ch. 12 Eq. 12-38 (D = 0, Chapter 7's default, leaves kappa_m as computed)
            if (o_.curvature_damping > 0.0 && !kappa_prev_.empty())
                for (int j = 0; j < nss_; ++j)
                    for (int i = 0; i < nqn_; ++i)
                        skappa_[us(j)][us(i)] = (skappa_[us(j)][us(i)] + o_.curvature_damping * kappa_prev_[us(j)][us(i)]) /
                                                (1.0 + o_.curvature_damping);
            kappa_prev_ = skappa_;
        }
        endwall_->begin_pass(pass, nqn_);
        ThroughflowIteration rec;
        rec.pass = pass;
        rec.bstar = bstar;
        int solved_upto = nqn_;
        bool choked = false, failed = false;
        for (int i = 0; i < nqn_; ++i) {
            int its = 0;
            bool conv = false;
            const bool ok = solve_qn(i, pass, its, conv);
            rec.flow_iterations += its;
            rec.max_mass_error = std::max(rec.max_mass_error, std::abs(q_[us(i)].out.mass_error));
            if (!ok) {
                choked = true;
                res.choked_station = i;
                solved_upto = i + 1;
                break;
            }
            if (!conv) {
                failed = true;
                solved_upto = i + 1;
                break;
            }
        }
        rec.choked = choked;
        res.stations.clear();
        for (int i = 0; i < solved_upto; ++i) res.stations.push_back(q_[us(i)].out);
        if (failed) {
            res.history.push_back(rec);
            res.status = ThroughflowResult::Status::NotConverged;
            res.message = "the flow field of quasi-normal " + std::to_string(solved_upto - 1) + " did not converge";
            break;
        }
        if (!full()) {  // Fig. 7-4: one marching pass
            res.history.push_back(rec);
            res.status = choked ? ThroughflowResult::Status::Choked : ThroughflowResult::Status::Converged;
            if (choked) res.message = "choked at quasi-normal " + std::to_string(res.choked_station);
            break;
        }
        if (choked && repositionings >= o_.choke_repositionings) {  // Fig. 7-3: "Have stream surfaces been repositioned?"
            res.history.push_back(rec);
            res.status = ThroughflowResult::Status::Choked;
            res.message = "choked at quasi-normal " + std::to_string(res.choked_station) + " after " +
                          std::to_string(repositionings) + " repositionings";
            break;
        }
        // Reposition Stream Surfaces (Sec. 7.4) on the quasi-normals that have a flow solution
        std::vector<Vec> ycalc(us(nqn_));
        double err = 0.0;
        for (int i = 0; i < nqn_; ++i) {
            QnData& d = q_[us(i)];
            if (!d.solved) continue;
            const Geo ge = geometry(i, d.y);
            Vec qv(us(nss_));
            for (int j = 0; j < nss_; ++j) {
                const auto& nd = d.out.nodes[us(j)];
                qv[us(j)] = nd.K_W * nd.K_B_blade * ge.r[us(j)] * nd.rho * nd.Wm * std::cos(ge.eps[us(j)]);
            }
            ycalc[us(i)] = reposition(ge, qv);
            for (int j = 0; j < nss_; ++j) err = std::max(err, std::abs(ycalc[us(i)][us(j)] - d.y[us(j)]) / d.g.ys);
        }
        rec.position_error = err;
        if (!choked && err <= o_.position_tolerance) {  // "Stream surfaces converged?"
            res.history.push_back(rec);
            res.status = ThroughflowResult::Status::Converged;
            break;
        }
        if (o_.adjust_bstar) bstar = bc.update(err);  // Eq. 7-59 (BstarControl)
        rec.bstar = bstar;
        // Novak's damping, Eq. 7-57, with M_m <= 0.95 (p. 165); dm the smallest meridional spacing to
        // the adjacent quasi-normals along the stream surface (Eq. 7-58 coordinates of this pass)
        double fmin = 1.0;
        for (int i = 0; i < nqn_; ++i) {
            QnData& d = q_[us(i)];
            if (!d.solved) continue;
            for (int j = 1; j + 1 < nss_; ++j) {
                const Vec& mj = sm_[us(j)];
                double dm = std::numeric_limits<double>::infinity();
                if (i > 0) dm = std::min(dm, mj[us(i)] - mj[us(i - 1)]);
                if (i + 1 < nqn_) dm = std::min(dm, mj[us(i + 1)] - mj[us(i)]);
                const double F = novak_damping(d.out.nodes[us(j)].mach_m, d.g.ys, dm, bstar, o_.max_mach_damping);
                fmin = std::min(fmin, F);
                d.y[us(j)] += F * (ycalc[us(i)][us(j)] - d.y[us(j)]);
            }
        }
        rec.min_damping = fmin;
        res.history.push_back(rec);
        ++repositionings;
        if (pass + 1 == passes) {
            res.status = ThroughflowResult::Status::NotConverged;
            res.message = "the stream surfaces did not converge in " + std::to_string(passes) + " passes";
        }
    }
    // stream surfaces through the solved quasi-normals (Eq. 7-58 on the reported positions)
    const int ns = ii(res.stations.size());
    for (int j = 0; j < nss_; ++j) {
        ThroughflowStreamSurface s;
        s.mass_fraction = static_cast<double>(j) / static_cast<double>(nss_ - 1);
        for (int i = 0; i < ns; ++i) {
            const auto& nd = res.stations[us(i)].nodes[us(j)];
            s.z.push_back(nd.z);
            s.r.push_back(nd.r);
            s.b.push_back(nd.b);
            s.b_mass.push_back(nd.b_mass);
        }
        s.m = ns >= 3 ? meridional_coordinate(s.z, s.r) : Vec(s.z.size(), 0.0);
        res.stream_surfaces.push_back(std::move(s));
    }
    for (int i = 0; i < ns; ++i)
        for (int j = 0; j < nss_; ++j) res.stations[us(i)].nodes[us(j)].m = res.stream_surfaces[us(j)].m[us(i)];
    for (int j = 0; j + 1 < nss_; ++j) {
        ThroughflowStreamSheet sh;
        sh.mass_fraction = (static_cast<double>(j) + 0.5) / static_cast<double>(nss_ - 1);
        for (int i = 0; i < ns; ++i) {
            const auto& x = res.stations[us(i)].sheets[us(j)];
            sh.z.push_back(x.z);
            sh.r.push_back(x.r);
            sh.b.push_back(x.b);
            sh.b_geometric.push_back(x.b_geometric);
        }
        sh.m = ns >= 3 ? meridional_coordinate(sh.z, sh.r) : Vec(sh.z.size(), 0.0);
        res.stream_sheets.push_back(std::move(sh));
    }
    return res;
}

}  // namespace

ThroughflowResult solve_throughflow(const ThroughflowCase& c)
{
    c.validate();
    Solver s(c);
    return s.run();
}

}  // namespace BladeGen::Turbolenza::Throughflow
