// Throughflow.hpp, equation by equation: tests that fail when one term of Aungier's (2003) Ch. 7 is
// wrong (the geometric epsilon in the mass equations of every model, Eqs. 7-19/7-16, 7-20, 7-29,
// 7-30, 7-43 to 7-49, 7-57, 7-59, 7-70, 7-39, K_B / K_W, Eq. 8-4, the Fig. 7-3 choke path).  Every
// measured number is printed ("[throughflow-eq] ...").
#include "turbolenza/BoundaryLayer.hpp"
#include "turbolenza/Throughflow.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <numbers>
#include <vector>

namespace TF = BladeGen::Turbolenza::Throughflow;

namespace {

constexpr double kPt = 101325.0, kTt = 288.15, kPi = std::numbers::pi;

const TF::ThroughflowModel kModels[] = {TF::ThroughflowModel::FullNormalEquilibrium, TF::ThroughflowModel::Approximate,
                                        TF::ThroughflowModel::SimpleNonIsentropic, TF::ThroughflowModel::SimpleIsentropic};

/// Quasi-normals with hub points (z_k, rh(z_k)), z_k = k dz, and shroud points (z_k - lean, rs(z_k - lean)).
TF::ThroughflowCase channel(TF::ThroughflowModel model, int nqn, double dz, double lean, double (*rh)(double),
                            double (*rs)(double))
{
    TF::ThroughflowCase c;
    c.model = model;
    for (int k = 0; k < nqn; ++k) {
        TF::QuasiNormal q;
        q.z_hub = k * dz;
        q.r_hub = rh(q.z_hub);
        q.z_shroud = q.z_hub - lean;
        q.r_shroud = rs(q.z_shroud);
        c.quasi_normals.push_back(q);
    }
    c.inlet = TF::InletProfile::uniform(kPt, kTt);
    c.mass_flow = 50.0;
    return c;
}

double flat_h(double) { return 0.3; }
double flat_s(double) { return 0.5; }
double cone_h(double z) { return 0.3 + z; }
double cone_s(double z) { return 0.5 + z; }
double bump_h(double z) { return 0.3 + 0.15 * std::sin(kPi * z); }
double bump_s(double z) { return 0.5 + 0.15 * std::sin(kPi * z); }
double diffuser_s(double z) { return 0.5 + z * std::tan(20.0 * kPi / 180.0); }

/// The geometric phi of node (i, j) recomputed by the test from the reported positions, by the rule
/// of the model level (straight walls of angles phi_h, phi_s): the simple models Eq. 7-2 by Eq. 7-70
/// on the marching pair (i - 1, i), wall-interpolated on quasi-normal 0; the approximate model linear
/// in y between the walls; the full model Eq. 7-2 with dr/dm by Eq. 7-67 inside and Eq. 7-70 at the
/// ends along the stream surface (m by Eq. 7-58 from the result).
double own_phi(const TF::ThroughflowResult& r, TF::ThroughflowModel model, std::size_t i, std::size_t j, double phi_h,
               double phi_s)
{
    const auto& st = r.stations[i];
    const auto& nd = st.nodes[j];
    const double t = nd.y / st.y_s;
    const auto& sf = r.stream_surfaces[j];
    const std::size_t n = sf.r.size();
    switch (model) {
    case TF::ThroughflowModel::Approximate: return phi_h + t * (phi_s - phi_h);
    case TF::ThroughflowModel::SimpleNonIsentropic:
    case TF::ThroughflowModel::SimpleIsentropic:
        if (i == 0) return phi_h + t * (phi_s - phi_h);
        return std::atan2(sf.r[i] - sf.r[i - 1], sf.z[i] - sf.z[i - 1]);
    case TF::ThroughflowModel::FullNormalEquilibrium: break;
    }
    double drdm = 0.0;
    if (i == 0) drdm = (sf.r[1] - sf.r[0]) / (sf.m[1] - sf.m[0]);
    else if (i + 1 == n) drdm = (sf.r[i] - sf.r[i - 1]) / (sf.m[i] - sf.m[i - 1]);
    else {
        const double m1 = sf.m[i - 1], m2 = sf.m[i], m3 = sf.m[i + 1], f1 = sf.r[i - 1], f2 = sf.r[i], f3 = sf.r[i + 1];
        drdm = ((m3 - m2) * (f2 - f1) / (m2 - m1) + (m2 - m1) * (f3 - f2) / (m3 - m2)) / (m3 - m1);
    }
    return std::asin(drdm);
}

/// mdot through quasi-normal i from the node data by the trapezoid rule, with the test's own lambda
/// (Eq. 7-1) and epsilon = phi - lambda (Eq. 7-3), phi by own_phi: 2 pi K_B int K_W r rho W_m cos(eps) dy
/// (Eq. 7-5).
double independent_mass(const TF::ThroughflowResult& r, std::size_t i, const TF::QuasiNormal& q, TF::ThroughflowModel model,
                        double phi_h, double phi_s)
{
    const auto& st = r.stations[i];
    const double lambda = std::atan2(q.z_hub - q.z_shroud, q.r_shroud - q.r_hub);
    double m = 0.0;
    for (std::size_t j = 0; j + 1 < st.nodes.size(); ++j) {
        const auto& a = st.nodes[j];
        const auto& b = st.nodes[j + 1];
        const double fa = a.K_W * a.r * a.rho * a.Wm * std::cos(own_phi(r, model, i, j, phi_h, phi_s) - lambda);
        const double fb = b.K_W * b.r * b.rho * b.Wm * std::cos(own_phi(r, model, i, j + 1, phi_h, phi_s) - lambda);
        m += 0.5 * (fa + fb) * (b.y - a.y);
    }
    return 2.0 * kPi * st.K_B_endwall * m;
}

TF::InletProfile swirl(double rh, double rs, int n, double (*ct)(double))
{
    TF::InletProfile p;
    for (int k = 0; k < n; ++k) {
        const double t = k / (n - 1.0);
        p.span.push_back(t);
        p.Pt.push_back(kPt);
        p.Tt.push_back(kTt);
        p.swirl_value.push_back(ct(rh + t * (rs - rh)));
    }
    return p;
}

TF::BladeRowSpec row_between(int qn_inlet, double omega, const std::shared_ptr<const TF::BladeRowModel>& model,
                             double le_h_z, double le_h_r, double le_s_z, double le_s_r, double te_h_z, double te_h_r,
                             double te_s_z, double te_s_r)
{
    TF::BladeRowSpec row;
    row.qn_inlet = qn_inlet;
    row.omega = omega;
    row.n_blades = 20;
    row.le_hub_z = le_h_z;
    row.le_hub_r = le_h_r;
    row.le_shroud_z = le_s_z;
    row.le_shroud_r = le_s_r;
    row.te_hub_z = te_h_z;
    row.te_hub_r = te_h_r;
    row.te_shroud_z = te_s_z;
    row.te_shroud_r = te_s_r;
    TF::BladeSectionGeometry g;
    g.chord = 0.05;
    row.sections = {g, g};
    row.sections[1].span = 1.0;
    row.model = model;
    return row;
}

std::shared_ptr<const TF::BladeRowModel> specified(double beta2, double loss, double kw = 1.0)
{
    return std::make_shared<TF::SpecifiedBladeRow>(std::vector<double>{0.0, 1.0}, std::vector<double>{beta2, beta2},
                                                   std::vector<double>{loss, loss}, std::vector<double>{kw, kw});
}

/// A row-level loss model with a fixed clearance Delta P_t at the shroud; records its last input.
class FixedRowLoss final : public TF::RowLossModel
{
public:
    explicit FixedRowLoss(double dpt) : dpt_(dpt) {}
    TF::RowLossResult evaluate(const TF::RowLossInput& in) const override
    {
        *last_ = in;
        TF::RowLossResult r;
        r.clearance_dPt = dpt_;
        r.clearance_side = TF::RowLossResult::Side::Shroud;
        r.leakage_mass_flow = 0.25;
        return r;
    }
    double dpt_;
    std::shared_ptr<TF::RowLossInput> last_ = std::make_shared<TF::RowLossInput>();
};

class RecordingEndwall final : public TF::EndwallModel
{
public:
    TF::EndwallResult evaluate(const TF::EndwallStation& st) override
    {
        if (st.qn == watch) seen = st;
        return {};
    }
    int watch = 0;
    TF::EndwallStation seen;
};

double max_mass_error(const TF::ThroughflowResult& r)
{
    double e = 0.0;
    for (const auto& s : r.stations) e = std::max(e, std::abs(s.mass_error));
    return e;
}

}  // namespace

TEST_CASE("throughflow-eq: every model keeps the geometric epsilon in Eq. 7-5 (cone, leaned quasi-normals)",
          "[turbolenza][throughflow]")
{
    struct Geometry
    {
        const char* name;
        double (*rh)(double);
        double (*rs)(double);
        double lean;
    };
    const Geometry geoms[] = {{"45-degree cone, radial", cone_h, cone_s, 0.0},
                              {"straight annulus, leaned 20 deg", flat_h, flat_s, 0.2 * std::tan(20.0 * kPi / 180.0)},
                              {"45-degree cone, leaned 18.4 deg", cone_h, cone_s, 0.05}};
    for (const auto& gm : geoms) {
        std::vector<double> mean_wm;
        for (auto model : kModels) {
            auto c = channel(model, 7, 0.1, gm.lean, gm.rh, gm.rs);
            c.options.n_stream_surfaces = 41;
            const auto r = TF::solve_throughflow(c);
            REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
            // the test's own phi (own_phi), not the code's
            const double phi_true = gm.rh == flat_h ? 0.0 : kPi / 4.0;
            double e = 0.0, eown = 0.0;
            for (std::size_t i = 0; i < r.stations.size(); ++i) {
                e = std::max(e, std::abs(independent_mass(r, i, c.quasi_normals[i], model, phi_true, phi_true) / c.mass_flow - 1.0));
                for (std::size_t j = 0; j < r.stations[i].nodes.size(); ++j)
                    eown = std::max(eown, std::abs(r.stations[i].nodes[j].phi - own_phi(r, model, i, j, phi_true, phi_true)));
            }
            CHECK(eown <= 1e-9);  // the code's phi is the model's rule (measured: see the report)
            const auto& st = r.stations[3];
            double wm = 0.0;
            for (const auto& nd : st.nodes) wm += nd.Wm / static_cast<double>(st.nodes.size());
            mean_wm.push_back(wm);
            // the geometric phi: exact on the straight walls of every quasi-normal (the first included),
            // and the interior surfaces within the measured deviation of the channel's angle
            double ephi = 0.0, ewall = 0.0;
            for (const auto& s : r.stations) {
                for (const auto& nd : s.nodes) ephi = std::max(ephi, std::abs(nd.phi - phi_true));
                ewall = std::max({ewall, std::abs(s.nodes.front().phi - phi_true), std::abs(s.nodes.back().phi - phi_true)});
            }
            CHECK(ewall <= 1e-12);
            std::printf("[throughflow-eq] %s, %s: max |mdot(trapezoid, own eps) / mdot - 1| = %.3e, mean W_m(QN 3) = %.4f, "
                        "max |phi - phi_wall| %.2e, max |phi - own phi| %.2e\n", gm.name, TF::to_string(model).c_str(), e, wm,
                        ephi, eown);
            CHECK(ephi <= 2e-2);  // measured: see the report
            CHECK(e <= 2e-3);  // the trapezoid rule on 41 surfaces (measured: see the report)
            CHECK(max_mass_error(r) <= 1e-11);
            if (gm.rh == flat_h) {  // uniform flow: rho W_m pi (rs^2 - rh^2) = mdot exactly
                const auto& nd = st.nodes[20];
                CHECK(std::abs(nd.rho * nd.Wm * kPi * (0.25 - 0.09) / c.mass_flow - 1.0) <= 1e-9);
            }
        }
        const auto [lo, hi] = std::minmax_element(mean_wm.begin(), mean_wm.end());
        std::printf("[throughflow-eq] %s: the four models' mean W_m spread %.3e\n", gm.name, (*hi - *lo) / *hi);
        CHECK((*hi - *lo) / *hi <= 0.02);
    }
}

TEST_CASE("throughflow-eq: Eq. 7-19 against Eq. 7-16 on a curved channel (inlet by alpha and by C_theta)",
          "[turbolenza][throughflow]")
{
    // the same physical inlet, once as an angle (beta' known, Eqs. 7-19 to 7-21 with omega = 0, p. 162)
    // and once as the swirl velocity that solution has (C_theta known, Eqs. 7-16 to 7-18): the two
    // forms of the normal-momentum equation must give the same quasi-normal 0, curvature included
    auto a = channel(TF::ThroughflowModel::FullNormalEquilibrium, 21, 0.05, 0.0, bump_h, bump_s);
    a.options.n_stream_surfaces = 21;
    a.inlet = TF::InletProfile::uniform(kPt, kTt, TF::InletProfile::Swirl::Alpha, 0.4);
    const auto ra = TF::solve_throughflow(a);
    REQUIRE(ra.status == TF::ThroughflowResult::Status::Converged);
    const auto& s0 = ra.stations[0];
    double ealpha = 0.0, kmax = 0.0;
    for (const auto& nd : s0.nodes) {
        ealpha = std::max(ealpha, std::abs(nd.alpha - 0.4));  // tan(alpha) = C_theta / C_m
        kmax = std::max(kmax, std::abs(nd.kappa_m));
    }
    CHECK(ealpha <= 1e-12);
    auto b = a;
    b.inlet.span.clear();
    b.inlet.Pt.clear();
    b.inlet.Tt.clear();
    b.inlet.swirl_value.clear();
    b.inlet.swirl = TF::InletProfile::Swirl::Ctheta;
    for (const auto& nd : s0.nodes) {
        b.inlet.span.push_back(nd.y / s0.y_s);
        b.inlet.Pt.push_back(kPt);
        b.inlet.Tt.push_back(kTt);
        b.inlet.swirl_value.push_back(nd.Ctheta);
    }
    b.inlet.span.front() = 0.0;
    b.inlet.span.back() = 1.0;
    const auto rb = TF::solve_throughflow(b);
    REQUIRE(rb.status == TF::ThroughflowResult::Status::Converged);
    double d = 0.0;
    for (std::size_t j = 0; j < s0.nodes.size(); ++j)
        d = std::max(d, std::abs(rb.stations[0].nodes[j].Wm / s0.nodes[j].Wm - 1.0));
    std::printf("[throughflow-eq] alpha inlet vs C_theta inlet on the bump: max |alpha - 0.4| %.2e, max |kappa_m| %.3f 1/m, "
                "max |dW_m| / W_m on QN 0 = %.3e\n", ealpha, kmax, d);
    CHECK(kmax >= 0.05);  // the test needs curvature (measured 0.12 1/m)
    CHECK(d <= 1e-4);     // the discrete d(r tan beta')/dy against d(r C_theta)/dy (measured 2.0e-6)
}

TEST_CASE("throughflow-eq: radial and leaned quasi-normals give the same flow on a curved channel (Eqs. 7-3, 7-29)",
          "[turbolenza][throughflow]")
{
    // the epsilon terms of Eqs. 7-16, 7-19 and 7-29 vanish on radial quasi-normals of a nearly axial
    // channel only where phi = 0; with the shroud ends leaned upstream the same hub points must see
    // the same W_m (up to the discretization).  A free vortex at the true radii of quasi-normal 0
    // makes the two inlets the same flow
    for (auto model : {TF::ThroughflowModel::FullNormalEquilibrium, TF::ThroughflowModel::Approximate}) {
        auto rad = channel(model, 21, 0.05, 0.0, bump_h, bump_s);
        auto lean = channel(model, 21, 0.05, 0.04, bump_h, bump_s);
        for (auto* c : {&rad, &lean}) {
            c->options.n_stream_surfaces = 21;
            c->inlet = swirl(c->quasi_normals[0].r_hub, c->quasi_normals[0].r_shroud, 401, [](double r) { return 45.0 / r; });
        }
        const auto rr = TF::solve_throughflow(rad);
        const auto rl = TF::solve_throughflow(lean);
        REQUIRE(rr.status == TF::ThroughflowResult::Status::Converged);
        REQUIRE(rl.status == TF::ThroughflowResult::Status::Converged);
        double d = 0.0, emax = 0.0;
        for (std::size_t i = 1; i + 1 < rr.stations.size(); ++i) {
            d = std::max(d, std::abs(rl.stations[i].nodes.front().Wm / rr.stations[i].nodes.front().Wm - 1.0));
            emax = std::max(emax, std::abs(rl.stations[i].nodes.front().epsilon));
        }
        std::printf("[throughflow-eq] bump, %s: radial vs leaned quasi-normals, max |dW_m| / W_m on the hub %.3e "
                    "(max |eps| hub %.3f rad)\n", TF::to_string(model).c_str(), d, emax);
        CHECK(emax >= 0.2);
        // measured 2.0e-3 (full: the discretization of two different grids), 8.7e-3 (approximate: its
        // linear phi and kappa_m depend on the quasi-normal's orientation)
        CHECK(d <= (model == TF::ThroughflowModel::FullNormalEquilibrium ? 5e-3 : 2e-2));
    }
}

TEST_CASE("throughflow-eq: the (1 + M_theta^2) sin(phi) / r term of Eq. 7-29: radial and leaned quasi-normals on a swirling cone",
          "[turbolenza][throughflow]")
{
    // a 45-degree cone with a strong free vortex (M_theta ~ 0.6 at the hub against M_m ~ 0.2): the
    // dW_m/dm term of Eq. 7-29 enters f1 with sin(epsilon), so radial (eps = 45 deg) and leaned
    // (eps = 26.6 deg) quasi-normals weight it differently; the hub W_m at the same hub points must agree
    for (auto model : {TF::ThroughflowModel::FullNormalEquilibrium, TF::ThroughflowModel::Approximate}) {
        auto rad = channel(model, 7, 0.1, 0.0, cone_h, cone_s);
        auto lean = channel(model, 7, 0.1, 0.05, cone_h, cone_s);
        for (auto* c : {&rad, &lean}) {
            c->options.n_stream_surfaces = 21;
            c->inlet = swirl(c->quasi_normals[0].r_hub, c->quasi_normals[0].r_shroud, 401, [](double r) { return 60.0 / r; });
        }
        const auto rr = TF::solve_throughflow(rad);
        const auto rl = TF::solve_throughflow(lean);
        REQUIRE(rr.status == TF::ThroughflowResult::Status::Converged);
        REQUIRE(rl.status == TF::ThroughflowResult::Status::Converged);
        double d = 0.0, mt = 0.0;
        for (std::size_t i = 1; i < rr.stations.size(); ++i) {
            d = std::max(d, std::abs(rl.stations[i].nodes.front().Wm / rr.stations[i].nodes.front().Wm - 1.0));
            mt = std::max(mt, rr.stations[i].nodes.front().mach_theta);
        }
        std::printf("[throughflow-eq] swirling cone, %s: radial vs leaned, max |dW_m| / W_m on the hub %.3e (M_theta hub %.3f)\n",
                    TF::to_string(model).c_str(), d, mt);
        CHECK(mt >= 0.4);
        // measured 4.1e-4 (full), 6.3e-3 (approximate: its linear phi and kappa_m depend on the
        // quasi-normal's orientation)
        CHECK(d <= (model == TF::ThroughflowModel::FullNormalEquilibrium ? 2e-3 : 1.5e-2));
    }
}

TEST_CASE("throughflow-eq: Eq. 7-70 gives the stream-surface slope at the first and last quasi-normal",
          "[turbolenza][throughflow]")
{
    auto c = channel(TF::ThroughflowModel::FullNormalEquilibrium, 21, 0.05, 0.0, bump_h, bump_s);
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    // the last pass's geometry is the one of the reported positions to the position tolerance
    double e = 0.0;
    const std::size_t L = r.stations.size() - 1;
    for (std::size_t j = 0; j < r.stream_surfaces.size(); ++j) {
        const auto& s = r.stream_surfaces[j];
        const double p0 = std::asin((s.r[1] - s.r[0]) / (s.m[1] - s.m[0]));       // Eq. 7-2 by Eq. 7-70
        const double pL = std::asin((s.r[L] - s.r[L - 1]) / (s.m[L] - s.m[L - 1]));
        e = std::max({e, std::abs(r.stations[0].nodes[j].phi - p0), std::abs(r.stations[L].nodes[j].phi - pL)});
    }
    std::printf("[throughflow-eq] Eq. 7-70 at the end quasi-normals: max |phi - two-point slope| = %.3e rad\n", e);
    CHECK(e <= 1e-5);  // the positions move by < 1e-7 y_s in the last pass
}

TEST_CASE("throughflow-eq: Novak's damping Eq. 7-57 by hand, and in the solver", "[turbolenza][throughflow]")
{
    // 1 / F = 1 + (1 - M^2) y_s^2 / (B* dm^2): M = 0.5, y_s = 0.2, dm = 0.1, B* = 8: 1 / F = 1.375
    CHECK(std::abs(TF::novak_damping(0.5, 0.2, 0.1, 8.0) - 1.0 / 1.375) <= 1e-15);
    // M_m capped at 0.95 (p. 165): M = 0.99 gives 1 / F = 1 + 0.0975 * 4 / 8
    CHECK(std::abs(TF::novak_damping(0.99, 0.2, 0.1, 8.0) - 1.0 / 1.04875) <= 1e-15);
    CHECK(TF::novak_damping(0.99, 0.2, 0.1, 8.0) == TF::novak_damping(0.95, 0.2, 0.1, 8.0));
    CHECK(std::abs(TF::novak_damping(0.3, 0.2, 0.05, 2.0) - 1.0 / (1.0 + 0.91 * 16.0 / 2.0)) <= 1e-15);
    CHECK_THROWS_AS(TF::novak_damping(0.5, 0.2, 0.0, 8.0), TF::ThroughflowError);
    // the solver: a forced vortex in a straight annulus (the stream surfaces move from the uniform
    // start); the smallest F of the last damped pass from the nodes of the converged solution
    auto c = channel(TF::ThroughflowModel::FullNormalEquilibrium, 7, 0.1, 0.0, flat_h, flat_s);
    c.inlet = swirl(0.3, 0.5, 2, [](double r) { return 200.0 * r; });
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    REQUIRE(r.history.size() >= 2);
    double fmin = 1.0;
    for (const auto& st : r.stations)
        for (std::size_t j = 1; j + 1 < st.nodes.size(); ++j)
            fmin = std::min(fmin, TF::novak_damping(st.nodes[j].mach_m, st.y_s, 0.1, 8.0));
    const double rec = r.history[r.history.size() - 2].min_damping;
    std::printf("[throughflow-eq] Novak F: solver %.10f, by hand from the converged nodes %.10f\n", rec, fmin);
    CHECK(std::abs(rec - fmin) <= 1e-6);
}

TEST_CASE("throughflow-eq: Eq. 7-59 against the running minimum of the position error", "[turbolenza][throughflow]")
{
    TF::BstarControl bc;
    bc.bstar = 8.0;
    bc.threshold = 0.05;
    CHECK(bc.update(1.0) == 8.0);
    CHECK(bc.update(0.9) == 8.0);
    CHECK(bc.update(0.93) == 8.0);   // +3.3 % over the minimum
    CHECK(bc.update(0.96) == 4.5);   // +6.7 % over the minimum 0.9: (8 + 1) / 2, although only +3.2 % per pass
    CHECK(bc.update(0.97) == 4.5);   // the minimum restarts at 0.96
    CHECK(bc.update(1.02) == 2.75);  // (4.5 + 1) / 2
    CHECK(bc.update(0.5) == 2.75);
    // in the solver: with B* = 1e4 the first change is (1e4 + 1) / 2 exactly
    auto c = channel(TF::ThroughflowModel::FullNormalEquilibrium, 21, 0.05, 0.0, bump_h, bump_s);
    c.inlet = swirl(0.3, 0.5, 2, [](double r) { return 150.0 * r; });
    c.options.bstar = 1e4;
    const auto r = TF::solve_throughflow(c);
    double first = 0.0;
    for (const auto& h : r.history)
        if (h.bstar != 1e4) {
            first = h.bstar;
            break;
        }
    std::printf("[throughflow-eq] B* = 1e4: first adjusted B* = %.4f, %zu passes, %s\n", first, r.history.size(),
                TF::to_string(r.status).c_str());
    CHECK(first == 5000.5);
}

TEST_CASE("throughflow-eq: a rotor with pre-swirl on a cone, Eqs. 7-43 to 7-49 surface by surface, with a row-level loss",
          "[turbolenza][throughflow]")
{
    // 10-degree cone, free-vortex pre-swirl r C_theta = 20, rotor omega = -400 rad/s, beta'_2 = 1.0,
    // omegabar = 0.04 and a clearance loss of 500 Pa at the shroud
    static const double t10 = std::tan(10.0 * kPi / 180.0);
    auto c = channel(TF::ThroughflowModel::SimpleNonIsentropic, 7, 0.1, 0.0, [](double z) { return 0.3 + t10 * z; },
                     [](double z) { return 0.5 + t10 * z; });
    c.options.n_stream_surfaces = 21;
    c.inlet = swirl(0.3, 0.5, 401, [](double r) { return 20.0 / r; });
    auto rl = std::make_shared<FixedRowLoss>(500.0);
    auto row = row_between(2, -400.0, specified(1.0, 0.04), 0.23, 0.3 + t10 * 0.23, 0.23, 0.5 + t10 * 0.23, 0.27,
                           0.3 + t10 * 0.27, 0.27, 0.5 + t10 * 0.27);
    row.row_loss = rl;
    c.rows = {row};
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    const auto& gas = c.gas;
    const double w = -400.0;
    const auto& q1 = r.stations[2];
    const auto& q2 = r.stations[3];
    double eI = 0.0, es = 0.0, erow = 0.0, dr = 0.0;
    for (std::size_t j = 0; j < q1.nodes.size(); ++j) {
        const auto& n1 = q1.nodes[j];
        const auto& n2 = q2.nodes[j];
        const double I = n1.H - w * n1.rCtheta;                        // Eq. 7-43
        const double Wt1 = n1.Ctheta - n1.r * w;                       // Eq. 7-44
        const double Hrel1 = gas.cp() * n1.T + 0.5 * (n1.Wm * n1.Wm + Wt1 * Wt1);  // Eqs. 7-45, 7-46
        const double Pt1 = gas.pressure(Hrel1, n1.s);
        const double Hrel2 = I + 0.5 * (n2.r * w) * (n2.r * w);        // Eqs. 7-47, 7-48 with r_2
        const double Pt2id = gas.pressure(Hrel2, n1.s);
        const double psi = static_cast<double>(j) / static_cast<double>(q1.nodes.size() - 1);
        const double row_loss = 2.0 * 500.0 * psi / (Pt1 - n1.P);      // the linear spread (distribute_row_loss)
        const double Pt2 = Pt2id - (0.04 + row_loss) * (Pt1 - n1.P);   // Eq. 7-49
        const double s2 = gas.entropy(Hrel2 / gas.cp(), Pt2);
        eI = std::max(eI, std::abs(n2.I - I) / std::abs(I));
        es = std::max(es, std::abs(n2.s - s2));
        erow = std::max(erow, std::abs(n2.row_loss_coefficient - row_loss));
        dr = std::max(dr, std::abs(n2.r - n1.r));
        // the identity behind it: s2 - s1 = -R ln(P't2 / P't2,id)
        CHECK(std::abs((n2.s - n1.s) + gas.gas_constant * std::log(Pt2 / Pt2id)) <= 1e-6);
    }
    std::printf("[throughflow-eq] rotor with pre-swirl: max |I2 - (H1 - omega rC_theta1)| / I %.2e, max |s2 - Eq. 7-49| %.2e J/(kg K), "
                "max |row loss - spread| %.2e, max |r2 - r1| %.4f m\n", eI, es, erow, dr);
    CHECK(eI <= 1e-12);
    CHECK(es <= 1e-6);
    CHECK(erow <= 1e-12);
    CHECK(dr >= 0.01);  // r1 != r2: Eq. 7-48 must take r_2
    // the row-level model got the inlet wall states in the row's frame
    const auto& in = *rl->last_;
    CHECK(std::abs(in.hub1.Wtheta - (in.hub1.Ctheta - in.hub1.r * w)) <= 1e-12);
    CHECK(std::abs(in.shroud2.Wtheta - (in.shroud2.Ctheta - in.shroud2.r * w)) <= 1e-9);
    CHECK(in.row_number == 1);
    CHECK(in.mass_flow == 50.0);
}

TEST_CASE("throughflow-eq: Eq. 7-20's cos(lambda): the constant-beta' rotor closed form on leaned quasi-normals",
          "[turbolenza][throughflow]")
{
    // simple model (epsilon = kappa_m = 0 in the momentum equation): dW_m/dy = cos(lambda) dW_m/dr, so the
    // closed form in r of the radial case holds on leaned quasi-normals too
    const double omega = -300.0, beta2 = 0.7, lean = 0.2 * std::tan(20.0 * kPi / 180.0);
    auto c = channel(TF::ThroughflowModel::SimpleNonIsentropic, 7, 0.1, lean, flat_h, flat_s);
    c.options.n_stream_surfaces = 41;
    c.rows = {row_between(2, omega, specified(beta2, 0.0), 0.22, 0.3, 0.22 - lean, 0.5, 0.28, 0.3, 0.28 - lean, 0.5)};
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    const auto& st = r.stations[3];
    const double t = std::tan(beta2), k = t * t / (1.0 + t * t), A = 2.0 * omega * t / ((1.0 + t * t) * (1.0 + k));
    const double C = (st.nodes.front().Wm + A * 0.3) * std::pow(0.3, k);
    double e = 0.0;
    for (const auto& nd : st.nodes) e = std::max(e, std::abs(nd.Wm - (C * std::pow(nd.r, -k) - A * nd.r)) / nd.Wm);
    std::printf("[throughflow-eq] rotor on 20-deg leaned quasi-normals: max |W_m - closed form| / W_m = %.3e\n", e);
    CHECK(e <= 1e-5);  // measured: see the report (the radial case: 2.8e-7 at 41 surfaces)
}

TEST_CASE("throughflow-eq: near and beyond choke on a 45-degree cone, all four models (Eq. 7-39 on a converged profile)",
          "[turbolenza][throughflow]")
{
    // Radial quasi-normals through 45-degree flow.  Two 1-D chokes of quasi-normal 0:
    //  * through the radial line itself, uniform flow: rho* a* pi (rs^2 - rh^2) cos(45 deg);
    //  * through the cross-section normal to the walls at the hub point of quasi-normal 0 (the
    //    narrowest section the flow at that quasi-normal has passed): rho* a* 2 pi int r dn, n from
    //    (0, 0.3) along (-sin 45, cos 45) to the shroud wall, 0.875 of the first.
    // The simple models (epsilon = kappa_m = 0 in the momentum equation: W_m uniform along the radial
    // line) must choke at the first; the full and approximate ones (W_m falls along the radial line
    // with the sin(eps) dW_m/dm term of Eqs. 7-16, 7-29) near the second.  Every choke must be the one
    // of a converged profile, at mean M_m = 1 (Eq. 7-39), with valid states (never h <= 0), and the
    // same for every mdot beyond it.
    const TF::ThroughflowGas gas;
    const double g = gas.gamma, Tst = kTt * 2.0 / (g + 1.0);
    const double rho_st = kPt / (gas.gas_constant * kTt) * std::pow(2.0 / (g + 1.0), 1.0 / (g - 1.0));
    const double flux = rho_st * std::sqrt(g * gas.gas_constant * Tst);
    const double mdot_radial = flux * kPi * (0.25 - 0.09) * std::cos(kPi / 4.0);
    const double tn = 0.2 / (std::cos(kPi / 4.0) + std::sin(kPi / 4.0));  // the normal's length to the shroud wall
    const double mdot_normal = flux * 2.0 * kPi * (0.3 * tn + std::cos(kPi / 4.0) * tn * tn / 2.0);
    std::printf("[throughflow-eq] 45-degree cone: 1-D choke through the radial QN 0 %.4f kg/s, through the normal "
                "section at its hub point %.4f kg/s (ratio %.4f)\n", mdot_radial, mdot_normal, mdot_normal / mdot_radial);
    const bool simple_model[] = {false, false, true, true};
    for (std::size_t k = 0; k < 4; ++k) {
        const auto model = kModels[k];
        // below every choke: converged, the target mass
        for (double mdot : {70.0, 72.0, 74.0}) {
            auto c = channel(model, 7, 0.1, 0.0, cone_h, cone_s);
            c.mass_flow = mdot;
            const auto r = TF::solve_throughflow(c);
            INFO(r.message);
            REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
            CHECK(max_mass_error(r) <= 1e-11);
            std::printf("[throughflow-eq] cone mdot %.0f, %s: converged, QN 0 mean M_m %.4f, %d flow iterations on QN 0\n",
                        mdot, TF::to_string(model).c_str(), r.stations[0].mean_mach_m, r.stations[0].iterations);
        }
        // the simple models converge up to their (radial) choke
        if (simple_model[k])
            for (double mdot : {78.0, 82.0, 85.0}) {
                auto c = channel(model, 7, 0.1, 0.0, cone_h, cone_s);
                c.mass_flow = mdot;
                const auto r = TF::solve_throughflow(c);
                REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
                CHECK(max_mass_error(r) <= 1e-11);
                std::printf("[throughflow-eq] cone mdot %.0f, %s: converged, QN 0 mean M_m %.4f\n", mdot,
                            TF::to_string(model).c_str(), r.stations[0].mean_mach_m);
            }
        // beyond the choke
        std::vector<double> mc;
        const std::vector<double> beyond = simple_model[k] ? std::vector<double>{90.0, 100.0, 150.0, 200.0}
                                                           : std::vector<double>{78.0, 82.0, 85.0, 90.0, 100.0, 150.0, 200.0};
        for (double mdot : beyond) {
            auto c = channel(model, 7, 0.1, 0.0, cone_h, cone_s);
            c.mass_flow = mdot;
            const auto r = TF::solve_throughflow(c);
            INFO(r.message);
            REQUIRE(r.status == TF::ThroughflowResult::Status::Choked);
            const auto& st = r.stations.back();
            CHECK(r.choked_station == 0);
            CHECK(st.converged);                            // a choke only on a converged profile
            CHECK(std::abs(st.mean_mach_m - 1.0) <= 1e-6);  // Eq. 7-39 holds at the choke
            for (const auto& nd : st.nodes) CHECK(nd.h > 0.0);  // never an invalid state
            double mmax = 0.0;
            for (const auto& nd : st.nodes) mmax = std::max(mmax, nd.mach_m);
            mc.push_back(st.mass_flow);
            std::printf("[throughflow-eq] cone mdot %.0f, %s: choked, choke mdot %.4f (%.4f of the radial, %.4f of the normal "
                        "1-D), max M_m %.4f, %d iterations, %zu passes\n", mdot, TF::to_string(model).c_str(), st.mass_flow,
                        st.mass_flow / mdot_radial, st.mass_flow / mdot_normal, mmax, st.iterations, r.history.size());
        }
        const auto [lo, hi] = std::minmax_element(mc.begin(), mc.end());
        CHECK((*hi - *lo) / *hi <= 1e-6);  // the same choke for every mdot input
        if (simple_model[k]) CHECK(std::abs(*lo / mdot_radial - 1.0) <= 1e-4);  // measured 1e-6
        else CHECK(std::abs(*lo / mdot_normal - 1.0) <= 0.02);                  // measured +1.0 % (full), -0.2 % (approx.)
    }
}

TEST_CASE("throughflow-eq: the limit of Eq. 7-30 near choke (45-degree cone at 74 kg/s)", "[turbolenza][throughflow]")
{
    // just below the full and approximate models' choke (75.8, 74.9) the hub nodes of quasi-normal 0
    // reach M_m > 0.95: the limit max(1 - M_m^2, 0.1) of Eq. 7-30 acts; without it (0.1 -> 1e-6) the flow
    // differs or does not converge
    for (auto model : {TF::ThroughflowModel::FullNormalEquilibrium, TF::ThroughflowModel::Approximate}) {
        auto c = channel(model, 7, 0.1, 0.0, cone_h, cone_s);
        c.mass_flow = 74.0;
        const auto r1 = TF::solve_throughflow(c);
        REQUIRE(r1.status == TF::ThroughflowResult::Status::Converged);
        bool active = false;
        double mmax = 0.0;
        for (const auto& st : r1.stations) {
            active = active || st.mach_limit_active;
            for (const auto& nd : st.nodes) mmax = std::max(mmax, nd.mach_m);
        }
        auto c2 = c;
        c2.options.min_one_minus_mm2 = 1e-6;
        double d = -1.0;
        std::string status = "refused";
        try {
            const auto r2 = TF::solve_throughflow(c2);
            status = TF::to_string(r2.status);
            if (r2.status == TF::ThroughflowResult::Status::Converged) {
                d = 0.0;
                for (std::size_t i = 0; i < r1.stations.size(); ++i)
                    for (std::size_t j = 0; j < r1.stations[i].nodes.size(); ++j)
                        d = std::max(d, std::abs(r2.stations[i].nodes[j].Wm / r1.stations[i].nodes[j].Wm - 1.0));
            }
        } catch (const TF::ThroughflowError&) {
        }
        std::printf("[throughflow-eq] cone mdot 74, %s: max M_m %.4f, Eq. 7-30 active %d; without the limit: %s, "
                    "max |dW_m| / W_m %.3e\n", TF::to_string(model).c_str(), mmax, active ? 1 : 0, status.c_str(), d);
        CHECK(active);
        CHECK((d < 0.0 || d >= 1e-5));
    }
}

TEST_CASE("throughflow-eq: the full model ends on a choke only after a repositioning (Fig. 7-3)", "[turbolenza][throughflow]")
{
    for (int need : {0, 1, 2}) {
        auto c = channel(TF::ThroughflowModel::FullNormalEquilibrium, 7, 0.1, 0.0, flat_h, flat_s);
        c.mass_flow = 500.0;
        c.options.choke_repositionings = need;
        const auto r = TF::solve_throughflow(c);
        std::printf("[throughflow-eq] choke with %d repositionings required: %s after %zu passes (%s)\n", need,
                    TF::to_string(r.status).c_str(), r.history.size(), r.message.c_str());
        CHECK(r.status == TF::ThroughflowResult::Status::Choked);
        CHECK(r.history.size() == static_cast<std::size_t>(need + 1));
    }
}

TEST_CASE("throughflow-eq: K_W in Eqs. 7-5 and 7-32, K_B in the dmdot of Eq. 7-33 (a stator, an entropy gradient)",
          "[turbolenza][throughflow]")
{
    // a loss-free stator with K_W = 0.8, K_B = 0.85 and an inlet entropy gradient (f3 != 0 on the
    // discharge): (a) the discharge mass from the nodes with K_W and K_B; (b) the stream-tube form
    // (rho dA / dmdot, Eq. 7-33) agrees with f3 / W_m (Eq. 7-15)
    std::vector<TF::ThroughflowResult> res;
    for (bool tube : {true, false}) {
        auto c = channel(TF::ThroughflowModel::SimpleNonIsentropic, 7, 0.1, 0.0, flat_h, flat_s);
        c.options.n_stream_surfaces = 41;
        c.options.stream_tube_form = tube;
        c.inlet.Pt = {0.98 * kPt, 1.02 * kPt};
        c.endwall = std::make_shared<TF::ConstantEndwall>(0.85);
        c.rows = {row_between(2, 0.0, specified(0.3, 0.0, 0.8), 0.23, 0.3, 0.23, 0.5, 0.27, 0.3, 0.27, 0.5)};
        res.push_back(TF::solve_throughflow(c));
        REQUIRE(res.back().status == TF::ThroughflowResult::Status::Converged);
        const double m3 = independent_mass(res.back(), 3, c.quasi_normals[3], c.model, 0.0, 0.0);
        std::printf("[throughflow-eq] stator K_W 0.8, K_B 0.85, %s: discharge mdot(nodes) / mdot - 1 = %.3e (K_W %.2f)\n",
                    tube ? "Eq. 7-33" : "f3 / W_m", m3 / c.mass_flow - 1.0, res.back().stations[3].nodes[5].K_W);
        CHECK(std::abs(m3 / c.mass_flow - 1.0) <= 1e-4);
    }
    double d = 0.0;
    for (std::size_t j = 0; j < res[0].stations[3].nodes.size(); ++j)
        d = std::max(d, std::abs(res[0].stations[3].nodes[j].Wm / res[1].stations[3].nodes[j].Wm - 1.0));
    std::printf("[throughflow-eq] stator: Eq. 7-33 against f3 / W_m on the discharge: max |dW_m| / W_m = %.3e\n", d);
    CHECK(d <= 1e-4);
}

TEST_CASE("throughflow-eq: Eq. 7-39 averages W_m / a: a swirling flow chokes at a mean M_m of 1", "[turbolenza][throughflow]")
{
    auto c = channel(TF::ThroughflowModel::SimpleNonIsentropic, 5, 0.1, 0.0, flat_h, flat_s);
    c.inlet = swirl(0.3, 0.5, 401, [](double r) { return 60.0 / r; });
    c.mass_flow = 500.0;
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Choked);
    const auto& st = r.stations.back();
    double mean = 0.0;  // trapezoid of W_m / a from the nodes
    for (std::size_t j = 0; j + 1 < st.nodes.size(); ++j)
        mean += 0.5 * (st.nodes[j].mach_m + st.nodes[j + 1].mach_m) * (st.nodes[j + 1].y - st.nodes[j].y);
    mean /= st.y_s;
    std::printf("[throughflow-eq] swirling choke: mdot %.4f, mean M_m %.6f (station), %.6f (trapezoid), max M %.4f\n",
                st.mass_flow, st.mean_mach_m, mean, st.nodes.back().mach);
    CHECK(std::abs(st.mean_mach_m - 1.0) <= 1e-6);
    CHECK(std::abs(mean - 1.0) <= 2e-3);
}

TEST_CASE("throughflow-eq: Eq. 8-4 and the row data handed to the end-wall model", "[turbolenza][throughflow]")
{
    auto c = channel(TF::ThroughflowModel::SimpleNonIsentropic, 7, 0.1, 0.0, flat_h, flat_s);
    c.inlet = swirl(0.3, 0.5, 401, [](double r) { return 10.0 / r; });
    c.quasi_normals[3].hub_wall_omega = -300.0;
    c.compressor_omega = -300.0;
    c.rows = {row_between(2, -300.0, specified(1.0, 0.0), 0.23, 0.3, 0.23, 0.5, 0.27, 0.3, 0.27, 0.5)};
    auto ew = std::make_shared<RecordingEndwall>();
    ew->watch = 3;
    c.endwall = ew;
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    const auto& e = ew->seen;
    CHECK(e.hub.wall_omega == -300.0);
    CHECK(e.shroud.wall_omega == 0.0);
    CHECK(std::abs(e.hub.Vtheta - (e.hub.Ctheta - e.hub.r * -300.0)) <= 1e-12);  // Eq. 8-4
    CHECK(e.shroud.Vtheta == e.shroud.Ctheta);
    CHECK(e.compressor_omega == -300.0);
    CHECK(e.upstream.row == 0);
    CHECK(e.upstream.row_number == 1);
    CHECK(std::abs(e.upstream.hub.Wtheta_inlet - (r.stations[2].nodes.front().Ctheta - 0.3 * -300.0)) <= 1e-9);
    CHECK(std::abs(e.upstream.hub.delta_rCtheta - (r.stations[3].nodes.front().rCtheta - r.stations[2].nodes.front().rCtheta)) <=
          1e-9);
    CHECK(std::abs(e.hub.mu - BladeGen::Turbolenza::sutherland_viscosity(e.hub.T)) <= 1e-18);
    std::printf("[throughflow-eq] end-wall station 3: hub V_theta %.4f, C_theta %.4f, mu %.4e at T %.2f K\n", e.hub.Vtheta,
                e.hub.Ctheta, e.hub.mu, e.hub.T);
}

TEST_CASE("throughflow-eq: the Approximate model's wall slopes from the wall points, and its dphi/dy (Eq. 7-29)",
          "[turbolenza][throughflow]")
{
    // a conical diffuser: flat hub, a straight 20-degree shroud, radial quasi-normals, wall slopes NOT
    // given: computed from the wall points (exact for straight walls); the approximate model's linear
    // phi then follows the full model's surfaces closely
    auto approx = channel(TF::ThroughflowModel::Approximate, 7, 0.1, 0.0, flat_h, diffuser_s);
    auto full = channel(TF::ThroughflowModel::FullNormalEquilibrium, 7, 0.1, 0.0, flat_h, diffuser_s);
    const auto ra = TF::solve_throughflow(approx);
    const auto rf = TF::solve_throughflow(full);
    REQUIRE(ra.status == TF::ThroughflowResult::Status::Converged);
    REQUIRE(rf.status == TF::ThroughflowResult::Status::Converged);
    double ephi = 0.0, d = 0.0;
    for (std::size_t i = 0; i < ra.stations.size(); ++i) {
        ephi = std::max({ephi, std::abs(ra.stations[i].nodes.front().phi),
                         std::abs(ra.stations[i].nodes.back().phi - 20.0 * kPi / 180.0)});
        for (std::size_t j = 0; j < ra.stations[i].nodes.size(); ++j)
            d = std::max(d, std::abs(ra.stations[i].nodes[j].Wm / rf.stations[i].nodes[j].Wm - 1.0));
    }
    std::printf("[throughflow-eq] diffuser: computed wall phi error %.2e rad; approximate vs full max |dW_m| / W_m %.3e\n",
                ephi, d);
    CHECK(ephi <= 1e-12);
    CHECK(d <= 1e-2);  // measured: see the report
}

TEST_CASE("throughflow-eq: the gas (Eq. 2-25 both ways, Sutherland's law)", "[turbolenza][throughflow]")
{
    TF::ThroughflowGas gas;
    for (double T : {200.0, 288.15, 400.0})
        for (double P : {5e4, 101325.0, 3e5}) {
            const double s = gas.entropy(T, P);
            CHECK(std::abs(gas.pressure(gas.cp() * T, s) / P - 1.0) <= 1e-13);
        }
    CHECK(gas.mu(273.15) == 1.716e-5);
    CHECK(std::abs(gas.mu(300.0) - BladeGen::Turbolenza::sutherland_viscosity(300.0)) <= 1e-20);
    gas.sutherland = false;
    CHECK(gas.mu(300.0) == gas.viscosity);
    CHECK_THROWS_AS(gas.mu(0.0), TF::ThroughflowError);
}
