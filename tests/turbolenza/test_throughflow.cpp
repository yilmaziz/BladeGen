// Throughflow.hpp: Aungier (2003) Ch. 7, the meridional through-flow analysis, against closed forms
// (radial equilibrium of a free vortex, a forced vortex and a constant relative angle in a rotating
// frame; uniform flow; a 45-degree cone for Eq. 7-58), its own defining equations (mass on every
// quasi-normal, Eq. 7-5; the stream-tube thickness), the model levels of Sec. 7.6, the interfaces and
// the refusals.  Every measured number is printed ("[throughflow] ...").
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

constexpr double kPt = 101325.0, kTt = 288.15;

/// A straight cylindrical annulus r_h .. r_s with radial quasi-normals at z = k dz.
TF::ThroughflowCase annulus(TF::ThroughflowModel model, double rh = 0.3, double rs = 0.5, int nqn = 7, double dz = 0.1)
{
    TF::ThroughflowCase c;
    c.model = model;
    for (int k = 0; k < nqn; ++k) {
        TF::QuasiNormal q;
        q.z_hub = q.z_shroud = k * dz;
        q.r_hub = rh;
        q.r_shroud = rs;
        c.quasi_normals.push_back(q);
    }
    c.inlet = TF::InletProfile::uniform(kPt, kTt);
    c.mass_flow = 50.0;
    return c;
}

/// Inlet C_theta(r) tabulated at n span points of the radial inlet quasi-normal.
TF::InletProfile swirl_profile(double rh, double rs, int n, double (*ct)(double))
{
    TF::InletProfile p;
    for (int k = 0; k < n; ++k) {
        const double t = k / (n - 1.0), r = rh + t * (rs - rh);
        p.span.push_back(t);
        p.Pt.push_back(kPt);
        p.Tt.push_back(kTt);
        p.swirl_value.push_back(ct(r));
    }
    return p;
}

double max_mass_error(const TF::ThroughflowResult& r)
{
    double e = 0.0;
    for (const auto& s : r.stations) e = std::max(e, std::abs(s.mass_error));
    return e;
}

const TF::ThroughflowModel kModels[] = {TF::ThroughflowModel::FullNormalEquilibrium, TF::ThroughflowModel::Approximate,
                                        TF::ThroughflowModel::SimpleNonIsentropic, TF::ThroughflowModel::SimpleIsentropic};

/// An end-wall model that records its calls and returns a fixed K_B.
class RecordingEndwall final : public TF::EndwallModel
{
public:
    explicit RecordingEndwall(double kb) : kb_(kb) {}
    void begin_pass(int pass, int n_qn) override
    {
        passes.push_back(pass);
        nqn = n_qn;
        last_qn = -1;
    }
    TF::EndwallResult evaluate(const TF::EndwallStation& st) override
    {
        if (st.qn < last_qn) order_ok = false;  // never goes back within a pass
        last_qn = st.qn;
        ++calls;
        last = st;
        TF::EndwallResult r;
        r.blockage_factor = kb_;
        return r;
    }
    double kb_;
    std::vector<int> passes;
    int nqn = 0, last_qn = -1, calls = 0;
    bool order_ok = true;
    TF::EndwallStation last;
};

}  // namespace

TEST_CASE("throughflow: uniform axial flow in a straight annulus (all four models)", "[turbolenza][throughflow]")
{
    for (auto model : kModels) {
        auto c = annulus(model);
        const auto r = TF::solve_throughflow(c);
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        // closed form: W_m uniform, mdot = rho W_m pi (rs^2 - rh^2) with the static state of W_m
        const auto& st = r.stations.back();
        double wmin = 1e300, wmax = 0.0;
        for (const auto& s : r.stations)
            for (const auto& nd : s.nodes) {
                wmin = std::min(wmin, nd.Wm);
                wmax = std::max(wmax, nd.Wm);
            }
        const auto& nd = st.nodes[3];
        const double area = std::numbers::pi * (0.5 * 0.5 - 0.3 * 0.3);
        const double mdot = nd.rho * nd.Wm * area;
        std::printf("[throughflow] uniform, %s: W_m = %.10f, spread %.2e, mdot(closed form) - mdot = %.2e, mass error %.2e\n",
                    TF::to_string(model).c_str(), nd.Wm, (wmax - wmin) / wmax, mdot - c.mass_flow, max_mass_error(r));
        CHECK((wmax - wmin) / wmax <= 1e-11);   // measured 1.7e-16 (full), 1.3e-12 (simplified: repositioning round-off)
        CHECK(std::abs(mdot - c.mass_flow) <= 1e-9 * c.mass_flow);
        CHECK(max_mass_error(r) <= 1e-11);
        // the sheets' b (actual sheet mass) equals their geometric thickness here (r linear in y: the
        // three-point mass integral is exact); the surface tubes' Eq. 12-39 form differs from the
        // Eq. 7-32 area term by the tube centroid offset (half tubes at the walls: first order)
        for (const auto& x : st.sheets) CHECK(std::abs(x.b - x.b_geometric) <= 1e-9 * x.b);
        double ew = 0.0, ei = 0.0;
        for (std::size_t j = 0; j < st.nodes.size(); ++j) {
            const double e = std::abs(st.nodes[j].b_mass / st.nodes[j].b - 1.0);
            double& slot = (j == 0 || j + 1 == st.nodes.size()) ? ew : ei;
            slot = std::max(slot, e);
        }
        std::printf("[throughflow] uniform, %s: |b_mass / b - 1| walls %.3e, interior %.3e\n", TF::to_string(model).c_str(),
                    ew, ei);
        CHECK(ew <= 1e-2);   // measured 4.6e-3 (n_ss = 11)
        CHECK(ei <= 3e-3);   // measured 1.4e-3
        // Eq. 7-58 on a cylinder: m = z
        for (std::size_t i = 0; i < r.stream_surfaces[5].m.size(); ++i)
            CHECK(std::abs(r.stream_surfaces[5].m[i] - 0.1 * static_cast<double>(i)) <= 1e-12);
    }
}

TEST_CASE("throughflow: free vortex r C_theta = const gives uniform W_m (simple radial equilibrium)",
          "[turbolenza][throughflow]")
{
    // 7-18 with I and s uniform and d(r C_theta)/dy = 0: f3 = 0, so W_m is uniform.  The inlet
    // profile is linear between its 401 points, so r C_theta is constant only to ~1e-7 (ours)
    for (auto model : kModels) {
        auto c = annulus(model);
        c.inlet = swirl_profile(0.3, 0.5, 401, [](double r) { return 32.0 / r; });
        const auto r = TF::solve_throughflow(c);
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        double spread = 0.0;
        for (const auto& s : r.stations) {
            double a = 1e300, b = 0.0;
            for (const auto& nd : s.nodes) {
                a = std::min(a, nd.Wm);
                b = std::max(b, nd.Wm);
                CHECK(std::abs(nd.rCtheta - 32.0) <= 1e-4);  // measured 1.4e-5 (the profile interpolation)
            }
            spread = std::max(spread, (b - a) / b);
        }
        std::printf("[throughflow] free vortex, %s: max (W_max - W_min) / W_max = %.3e, mass error %.2e\n",
                    TF::to_string(model).c_str(), spread, max_mass_error(r));
        CHECK(spread <= 1e-6);  // measured 1.3e-8 (the profile interpolation)
        CHECK(max_mass_error(r) <= 1e-11);
    }
}

TEST_CASE("throughflow: forced vortex C_theta = Omega r against W_m^2 = W_mh^2 - 2 Omega^2 (r^2 - r_h^2)",
          "[turbolenza][throughflow]")
{
    // uniform H and s: W_m dW_m/dr = -(C_theta / r) d(r C_theta)/dr = -2 Omega^2 r (Eq. 7-18 with
    // y = r - r_h); exact for a compressible gas too, since rho does not enter
    for (bool tube : {true, false}) {
    double prev = 0.0;
    for (int nss : {6, 11, 21, 41}) {
        auto c = annulus(TF::ThroughflowModel::SimpleNonIsentropic);
        c.options.n_stream_surfaces = nss;
        c.options.stream_tube_form = tube;
        c.inlet = swirl_profile(0.3, 0.5, 2, [](double r) { return 150.0 * r; });  // linear: exact
        const auto r = TF::solve_throughflow(c);
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        const auto& st = r.stations[3];
        const double wh = st.nodes.front().Wm;
        double e = 0.0;
        for (const auto& nd : st.nodes) {
            const double ref = std::sqrt(wh * wh - 2.0 * 150.0 * 150.0 * (nd.r * nd.r - 0.09));
            e = std::max(e, std::abs(nd.Wm - ref) / ref);
        }
        std::printf("[throughflow] forced vortex, %s, n_ss = %2d: max |W_m - closed form| / W_m = %.3e (W_m hub %.4f)\n",
                    tube ? "Eq. 7-33" : "f3 / W_m", nss, e, wh);
        // measured ratios: Eq. 7-33 2.33, 2.98, 3.40 (toward 4: second order); f3 / W_m 5.1, 6.0, 6.8
        if (prev > 0.0) CHECK(prev / e >= (tube ? 2.2 : 4.5));
        prev = e;
        CHECK(max_mass_error(r) <= 1e-11);
    }
    CHECK(prev <= (tube ? 2e-4 : 1e-5));  // n_ss = 41: measured 1.1e-4 and 4.2e-6
    }
}

TEST_CASE("throughflow: the order of both f4 forms on a 45-degree cone with swirl (self-convergence)",
          "[turbolenza][throughflow]")
{
    // W_m at psi = 0, 0.5, 1 on quasi-normal 3 against 81 surfaces, simple model, forced-vortex inlet
    for (bool tube : {true, false}) {
        auto run = [&](int nss) {
            auto c = annulus(TF::ThroughflowModel::SimpleNonIsentropic);
            for (auto& q : c.quasi_normals) {
                q.r_hub = 0.3 + q.z_hub;
                q.r_shroud = 0.5 + q.z_shroud;
            }
            c.options.n_stream_surfaces = nss;
            c.options.stream_tube_form = tube;
            c.inlet = swirl_profile(0.3, 0.5, 2, [](double r) { return 150.0 * r; });
            const auto r = TF::solve_throughflow(c);
            REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
            const auto& n = r.stations[3].nodes;
            return std::vector<double>{n.front().Wm, n[n.size() / 2].Wm, n.back().Wm};
        };
        const auto ref = run(81);
        double prev = 0.0;
        for (int nss : {11, 21, 41}) {
            const auto v = run(nss);
            double e = 0.0;
            for (std::size_t k = 0; k < 3; ++k) e = std::max(e, std::abs(v[k] / ref[k] - 1.0));
            std::printf("[throughflow] cone order, %s, n_ss = %2d: max |W_m / W_m(81) - 1| at psi 0, 0.5, 1 = %.3e\n",
                        tube ? "Eq. 7-33" : "f3 / W_m", nss, e);
            if (prev > 0.0) CHECK(prev / e >= 2.5);  // measured: see the report
            prev = e;
        }
    }
}

TEST_CASE("throughflow: rotating frame, constant beta' (Eqs. 7-19 to 7-21) against the closed form",
          "[turbolenza][throughflow]")
{
    // a loss-free rotor with constant discharge beta' = b in a straight annulus, axial inflow: I and s
    // uniform, so (1 + t^2) dW_m/dr = -t^2 W_m / r - 2 omega t, t = tan b, with the solution
    // W_m = C r^-k - 2 omega t r / ((1 + t^2)(1 + k)), k = t^2 / (1 + t^2)
    const double omega = -300.0, beta2 = 0.7;  // omega < 0: the BladeGen compressor sense; beta2 > 0
    // the closed form has no curvature: the simple model is checked; the full model moves the stream
    // surfaces radially behind the rotor, so kappa_m, epsilon != 0 there (its departure is printed)
    for (auto model : {TF::ThroughflowModel::SimpleNonIsentropic, TF::ThroughflowModel::FullNormalEquilibrium}) {
    const bool simple = model == TF::ThroughflowModel::SimpleNonIsentropic;
    double prev = 0.0;
    for (int nss : {6, 11, 21, 41}) {
        auto c = annulus(model);
        c.options.n_stream_surfaces = nss;
        TF::BladeRowSpec row;
        row.name = "rotor";
        row.qn_inlet = 2;
        row.omega = omega;
        row.n_blades = 20;
        row.le_hub_z = row.le_shroud_z = 0.23;
        row.te_hub_z = row.te_shroud_z = 0.27;
        row.le_hub_r = row.te_hub_r = 0.3;
        row.le_shroud_r = row.te_shroud_r = 0.5;
        TF::BladeSectionGeometry g;
        g.chord = 0.05;
        g.span = 0.0;
        row.sections.push_back(g);
        g.span = 1.0;
        row.sections.push_back(g);
        row.model = std::make_shared<TF::SpecifiedBladeRow>(std::vector<double>{0.0, 1.0}, std::vector<double>{beta2, beta2},
                                                            std::vector<double>{0.0, 0.0});
        c.rows.push_back(row);
        const auto r = TF::solve_throughflow(c);
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        const auto& st = r.stations[3];
        CHECK(st.frame_omega == omega);
        const double t = std::tan(beta2), k = t * t / (1.0 + t * t), A = 2.0 * omega * t / ((1.0 + t * t) * (1.0 + k));
        const double C = (st.nodes.front().Wm + A * 0.3) * std::pow(0.3, k);
        double e = 0.0, ds = 0.0, dI = 0.0;
        for (const auto& nd : st.nodes) {
            const double ref = C * std::pow(nd.r, -k) - A * nd.r;
            e = std::max(e, std::abs(nd.Wm - ref) / ref);
            CHECK(std::abs(nd.beta - beta2) <= 1e-12);
            ds = std::max(ds, std::abs(nd.s - st.nodes.front().s));
            dI = std::max(dI, std::abs(nd.I - st.nodes.front().I));
        }
        // the absolute discharge (Eqs. 7-50, 7-51): C_theta = W_m tan beta' + r omega, H = I + omega r C_theta
        const auto& nd = st.nodes[static_cast<std::size_t>(nss / 2)];
        CHECK(std::abs(nd.Ctheta - (nd.Wm * t + nd.r * omega)) <= 1e-9);
        CHECK(std::abs(nd.H - (nd.I + omega * nd.r * nd.Ctheta)) <= 1e-6);
        std::printf("[throughflow] rotor beta' = const, %s, n_ss = %2d: max |W_m - closed form| / W_m = %.3e, |ds| %.1e, "
                    "|dI| %.1e, outer passes %zu\n", TF::to_string(model).c_str(), nss, e, ds, dI, r.history.size());
        if (prev > 0.0 && simple) CHECK(prev / e >= 1.5);  // measured: see the report
        prev = e;
        CHECK(ds <= 1e-9);
        CHECK(max_mass_error(r) <= 1e-11);
        // rothalpy is conserved through the row (Eq. 7-47): I2 = H1 (axial inflow)
        CHECK(std::abs(st.nodes.front().I - r.stations[2].nodes.front().H) <= 1e-6);
    }
    if (simple) CHECK(prev <= 1e-2);
    }
}

TEST_CASE("throughflow: the four model levels on a straight annulus, with and without an entropy gradient",
          "[turbolenza][throughflow]")
{
    // no curvature, epsilon = 0: the models differ only in the entropy term (Sec. 7.6)
    TF::InletProfile pr = swirl_profile(0.3, 0.5, 2, [](double r) { return 200.0 * r; });
    auto run = [&](TF::ThroughflowModel m, bool entropy) {
        auto c = annulus(m);
        c.inlet = pr;
        if (entropy) c.inlet.Pt = {0.95 * kPt, 1.05 * kPt};  // Tt uniform: s falls toward the shroud
        return TF::solve_throughflow(c);
    };
    for (bool entropy : {false, true}) {
        const auto ref = run(TF::ThroughflowModel::FullNormalEquilibrium, entropy);
        for (auto model : kModels) {
            const auto r = run(model, entropy);
            REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
            double d = 0.0;
            for (std::size_t i = 0; i < r.stations.size(); ++i)
                for (std::size_t j = 0; j < r.stations[i].nodes.size(); ++j)
                    d = std::max(d, std::abs(r.stations[i].nodes[j].Wm - ref.stations[i].nodes[j].Wm) /
                                        ref.stations[i].nodes[j].Wm);
            std::printf("[throughflow] straight annulus, entropy gradient %s, %s vs full: max |dW_m| / W_m = %.3e\n",
                        entropy ? "yes" : "no", TF::to_string(model).c_str(), d);
            // measured 1.8e-7 and 5.9e-8: the full model's positions are converged to 1e-7 of y_s only
            if (!entropy || model != TF::ThroughflowModel::SimpleIsentropic) CHECK(d <= 1e-6);
            else CHECK(d >= 1e-2);  // dropping T ds/dy matters: measured in the report
        }
    }
}

TEST_CASE("throughflow: Eq. 7-58 with its square root on a 45-degree conical annulus", "[turbolenza][throughflow]")
{
    // hub r = 0.3 + z, shroud r = 0.5 + z, radial quasi-normals: on the wall stream surfaces
    // m = sqrt(2) (z - z_0) exactly; the printed form without the root would give 2 (z - z_0)
    auto c = annulus(TF::ThroughflowModel::FullNormalEquilibrium);
    for (std::size_t k = 0; k < c.quasi_normals.size(); ++k) {
        auto& q = c.quasi_normals[k];
        q.r_hub = 0.3 + q.z_hub;
        q.r_shroud = 0.5 + q.z_shroud;
        q.phi_hub = q.phi_shroud = std::numbers::pi / 4.0;
    }
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    double e = 0.0;
    for (int j : {0, 10})
        for (std::size_t i = 0; i < r.stations.size(); ++i)
            e = std::max(e, std::abs(r.stream_surfaces[static_cast<std::size_t>(j)].m[i] - std::sqrt(2.0) * 0.1 * static_cast<double>(i)));
    double ephi = 0.0;
    for (const auto& s : r.stations)
        for (int j : {0, 10}) ephi = std::max(ephi, std::abs(s.nodes[static_cast<std::size_t>(j)].phi - std::numbers::pi / 4.0));
    std::printf("[throughflow] 45-degree cone: max |m - sqrt(2) dz| = %.3e on the walls, max |phi - 45 deg| = %.3e, %zu passes\n",
                e, ephi, r.history.size());
    CHECK(e <= 1e-12);
    CHECK(ephi <= 1e-9);
    CHECK(max_mass_error(r) <= 1e-11);
}

TEST_CASE("throughflow: Novak damping, the position error history and the B* adjustment (Eqs. 7-57, 7-59)",
          "[turbolenza][throughflow]")
{
    // a cone with closely spaced quasi-normals (y_s / dm = 4): curvature couples the quasi-normals
    auto make = [](double bstar) {
        auto c = annulus(TF::ThroughflowModel::FullNormalEquilibrium, 0.3, 0.5, 21, 0.05);
        for (auto& q : c.quasi_normals) {
            const double z = q.z_hub, dr = 0.15 * std::sin(std::numbers::pi * z);  // a bump in both walls
            q.r_hub = 0.3 + dr;
            q.r_shroud = 0.5 + dr;
        }
        c.inlet = swirl_profile(0.3, 0.5, 2, [](double r) { return 150.0 * r; });
        c.options.bstar = bstar;
        return c;
    };
    const auto r8 = TF::solve_throughflow(make(8.0));
    REQUIRE(r8.status == TF::ThroughflowResult::Status::Converged);
    std::printf("[throughflow] Novak B* = 8: %zu passes; position errors:", r8.history.size());
    for (std::size_t k = 0; k < r8.history.size(); k += std::max<std::size_t>(1, r8.history.size() / 12))
        std::printf(" %.2e", r8.history[k].position_error);
    std::printf(" ... %.2e; min F %.3f; final B* %.3f\n", r8.history.back().position_error, r8.history.front().min_damping,
                r8.history.back().bstar);
    // the error falls by orders of magnitude to the tolerance
    CHECK(r8.history.back().position_error <= 1e-7);
    CHECK(r8.history.front().position_error >= 1e3 * r8.history.back().position_error);
    // with little damping the error grows and Eq. 7-59 lowers B*; the run still converges
    const auto rbig = TF::solve_throughflow(make(1e4));
    std::printf("[throughflow] Novak B* = 1e4: status %s, %zu passes, final B* %.3f, final error %.2e\n",
                TF::to_string(rbig.status).c_str(), rbig.history.size(), rbig.history.back().bstar,
                rbig.history.back().position_error);
    CHECK(rbig.history.back().bstar < 1e4);
    CHECK(rbig.status == TF::ThroughflowResult::Status::Converged);
    // without Eq. 7-59 the same case does not converge
    auto cfix = make(1e4);
    cfix.options.adjust_bstar = false;
    cfix.options.max_outer_iterations = 60;
    bool diverged = false;
    try {
        const auto rf = TF::solve_throughflow(cfix);
        diverged = rf.status != TF::ThroughflowResult::Status::Converged;
        std::printf("[throughflow] Novak B* = 1e4 without Eq. 7-59: status %s after %zu passes, error %.2e\n",
                    TF::to_string(rf.status).c_str(), rf.history.size(), rf.history.back().position_error);
    } catch (const TF::ThroughflowError& e) {
        diverged = true;
        std::printf("[throughflow] Novak B* = 1e4 without Eq. 7-59: refused (%s)\n", e.what());
    }
    CHECK(diverged);
}

TEST_CASE("throughflow: the end-wall model interface (call order, K_B in Eq. 7-5 only)", "[turbolenza][throughflow]")
{
    auto c = annulus(TF::ThroughflowModel::FullNormalEquilibrium);
    auto ew = std::make_shared<RecordingEndwall>(0.95);
    c.endwall = ew;
    const auto r = TF::solve_throughflow(c);
    REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
    CHECK(ew->order_ok);
    CHECK(ew->passes.size() == r.history.size());
    CHECK(ew->nqn == 7);
    CHECK(ew->calls >= 7 * static_cast<int>(r.history.size()));
    CHECK(ew->last.qn == 6);
    CHECK(ew->last.mass_flow == 50.0);
    CHECK(ew->last.hub.r == 0.3);
    CHECK(ew->last.shroud.r == 0.5);
    CHECK(ew->last.hub.mu == c.gas.mu(ew->last.hub.T));  // Sutherland at the edge static temperature
    // K_B multiplies the whole integral: rho W_m (the uniform flow) carries mdot / K_B
    const auto& nd = r.stations[3].nodes[4];
    const double area = std::numbers::pi * (0.25 - 0.09);
    CHECK(r.stations[3].K_B_endwall == 0.95);
    CHECK(std::abs(nd.rho * nd.Wm * area * 0.95 - 50.0) <= 1e-9 * 50.0);
    // the sheets carry the actual mass: b = K_B times their geometric b
    for (const auto& x : r.stations[3].sheets) CHECK(std::abs(x.b - 0.95 * x.b_geometric) <= 1e-9 * x.b);
    CHECK(max_mass_error(r) <= 1e-11);
}

TEST_CASE("throughflow: row-level losses are distributed linearly and enter Eq. 7-49", "[turbolenza][throughflow]")
{
    const std::vector<double> psi = {0.0, 0.25, 0.5, 0.75, 1.0};
    const auto a = TF::distribute_row_loss(100.0, TF::RowLossResult::Side::Shroud, psi);
    const auto b = TF::distribute_row_loss(100.0, TF::RowLossResult::Side::Hub, psi);
    CHECK(a.front() == 0.0);
    CHECK(a.back() == 200.0);
    CHECK(b.front() == 200.0);
    CHECK(b.back() == 0.0);
    double mean = 0.0;  // trapezoid = exact for linear
    for (std::size_t k = 0; k + 1 < psi.size(); ++k) mean += 0.5 * (a[k] + a[k + 1]) * 0.25;
    CHECK(std::abs(mean - 100.0) <= 1e-12);
    CHECK_THROWS_AS(TF::distribute_row_loss(1.0, TF::RowLossResult::Side::Hub, {0.0, 1.5}), TF::ThroughflowError);
}

TEST_CASE("throughflow: refusals", "[turbolenza][throughflow]")
{
    auto ok = annulus(TF::ThroughflowModel::SimpleNonIsentropic);
    CHECK_NOTHROW(ok.validate());
    auto c = ok;
    c.quasi_normals.resize(2);
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);
    c = ok;
    c.mass_flow = 0.0;
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);
    c = ok;
    c.quasi_normals[3].inside_blade_row = true;
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);
    c = ok;
    c.quasi_normals[3].z_hub = c.quasi_normals[2].z_hub;
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);
    c = ok;
    c.options.n_stream_surfaces = 2;
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);
    c = ok;
    c.inlet.Pt = {kPt, -1.0};
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);
    // a row: no model, inside quasi-normals, overlapping rows
    TF::BladeRowSpec row;
    row.qn_inlet = 2;
    row.n_blades = 10;
    TF::BladeSectionGeometry g;
    g.chord = 0.05;
    row.sections = {g, g};
    row.sections[1].span = 1.0;
    c = ok;
    c.rows = {row};
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);  // no model
    row.model = std::make_shared<TF::SpecifiedBladeRow>(std::vector<double>{0.0, 1.0}, std::vector<double>{0.1, 0.1},
                                                        std::vector<double>{0.0, 0.0});
    row.qn_outlet = 4;
    c.rows = {row};
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);  // quasi-normals inside the row (Ch. 12)
    row.qn_outlet = -1;
    c.rows = {row, row};
    CHECK_THROWS_AS(TF::solve_throughflow(c), TF::ThroughflowError);
    // the specified row's profiles
    CHECK_THROWS_AS(TF::SpecifiedBladeRow({0.0, 0.9}, {0.1, 0.1}, {0.0, 0.0}), TF::ThroughflowError);
    CHECK_THROWS_AS(TF::SpecifiedBladeRow({0.0, 1.0}, {0.1, 2.0}, {0.0, 0.0}), TF::ThroughflowError);
    CHECK_THROWS_AS(TF::SpecifiedBladeRow({0.0, 1.0}, {0.1, 0.1}, {0.0, 0.0}, {1.0, 1.5}), TF::ThroughflowError);
    CHECK_THROWS_AS(TF::ConstantEndwall(0.0), TF::ThroughflowError);
    // a mass flow beyond choke: reported, not thrown (Fig. 7-4 ends on it)
    c = ok;
    c.mass_flow = 500.0;
    const auto r = TF::solve_throughflow(c);
    CHECK(r.status == TF::ThroughflowResult::Status::Choked);
    CHECK(r.choked_station == 0);
    std::printf("[throughflow] mdot = 500 kg/s in the annulus: %s at quasi-normal %d, computed mdot %.4f, mean M_m %.4f\n",
                TF::to_string(r.status).c_str(), r.choked_station, r.stations.back().mass_flow, r.stations.back().mean_mach_m);
    CHECK(r.stations.back().mean_mach_m >= 1.0 - 1e-9);
    // the same with Ch. 12's test, Eq. 12-37: the maximum of mdot(W_m(0)) is the choke
    c.options.choke_test = TF::ChokeTest::MassFlowDerivative;
    const auto r2 = TF::solve_throughflow(c);
    CHECK(r2.status == TF::ThroughflowResult::Status::Choked);
    std::printf("[throughflow] the same with Eq. 12-37: choke mdot %.4f, mean M_m %.4f\n", r2.stations.back().mass_flow,
                r2.stations.back().mean_mach_m);
    CHECK(r2.stations.back().mass_flow >= r.stations.back().mass_flow - 1e-6);  // 7-39 stops at or before the maximum
}
