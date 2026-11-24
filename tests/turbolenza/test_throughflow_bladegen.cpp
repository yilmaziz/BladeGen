// ThroughflowBladeGen.hpp: the through-flow on BladeGen's default channel and blade (Aungier 2003
// Ch. 7): the automatic quasi-normals, the exact wall slope and curvature, the four model levels on the
// curved channel (where full normal equilibrium and the linear approximation differ), and a rotating
// SpecifiedBladeRow case converging.  Every measured number is printed ("[throughflow-bg] ...").
#include "turbolenza/ThroughflowBladeGen.hpp"

#include "blade/FlowPath.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <numbers>
#include <vector>

namespace TF = BladeGen::Turbolenza::Throughflow;
namespace D = BladeGen::Design;

namespace {

TF::ThroughflowCase channel_case(const TF::BladeGenThroughflowGeometry& g, TF::ThroughflowModel m, double mdot)
{
    TF::ThroughflowCase c;
    c.model = m;
    c.quasi_normals = g.quasi_normals;
    c.inlet = TF::InletProfile::uniform(101325.0, 288.15);
    c.mass_flow = mdot;
    c.options.max_outer_iterations = 2000;
    return c;
}

double max_mass_error(const TF::ThroughflowResult& r)
{
    double e = 0.0;
    for (const auto& s : r.stations) e = std::max(e, std::abs(s.mass_error));
    return e;
}

}  // namespace

TEST_CASE("throughflow bladegen: the wall slope and curvature are those of the wall spline", "[turbolenza][throughflow]")
{
    // a straight wall: phi = atan(dr/dz), kappa = 0; a circle-like check through a B-spline line
    const auto line = BladeGen::Blade::MeridionalCurve::line({0.0, 1.0}, {2.0, 3.0}).curve();
    double phi = 0.0, kappa = 1.0;
    TF::wall_slope_curvature(line, 0.5 * (line.t0() + line.t1()), phi, kappa);
    CHECK(std::abs(phi - std::numbers::pi / 4.0) <= 1e-12);
    CHECK(std::abs(kappa) <= 1e-9);
    // the default hub (Bell & Mehta, rising toward the TE): phi > 0 inside, 0 on the extensions
    const BladeGen::Blade::FlowPath fp(D::BladeDesign::defaults().meridional);
    TF::wall_slope_curvature(fp.hub(), 1.5, phi, kappa);
    CHECK(phi > 0.0);
    TF::wall_slope_curvature(fp.hub(), 0.5, phi, kappa);
    CHECK(std::abs(phi) <= 1e-12);
    // kappa = -dphi/dm by a central difference of phi along the wall
    const double u = 1.3, du = 1e-5;
    double p1 = 0.0, p2 = 0.0, k0 = 0.0, kd = 0.0;
    TF::wall_slope_curvature(fp.hub(), u - du, p1, kd);
    TF::wall_slope_curvature(fp.hub(), u + du, p2, kd);
    TF::wall_slope_curvature(fp.hub(), u, phi, k0);
    const auto a = fp.hub().eval(u - du), b = fp.hub().eval(u + du);
    const double dm = (b - a).norm();
    std::printf("[throughflow-bg] hub kappa at u = 1.3: %.6e 1/mm, -dphi/dm by differences %.6e\n", k0, -(p2 - p1) / dm);
    CHECK(std::abs(k0 + (p2 - p1) / dm) <= 1e-6 * std::max(1.0, std::abs(k0)) + 1e-9);
}

TEST_CASE("throughflow bladegen: the automatic quasi-normals of the default design", "[turbolenza][throughflow]")
{
    D::Evaluator ev;
    const auto d = D::BladeDesign::defaults();
    TF::BladeGenThroughflowOptions opt;
    opt.rpm = 15000.0;
    opt.n_inlet = 3;
    opt.n_exit = 3;
    const auto g = TF::throughflow_geometry(d, ev, opt);
    REQUIRE(g.quasi_normals.size() == 10u);  // far inlet, 3, before LE, after TE, 3, far exit
    CHECK(g.qn_before_le == 4);
    CHECK(g.qn_after_te == 5);
    for (const auto& q : g.quasi_normals) CHECK(q.z_hub == q.z_shroud);  // radial
    CHECK(g.quasi_normals[4].z_hub < -0.015);  // upstream of both LE points (hub LE at -15 mm)
    CHECK(g.quasi_normals[5].z_hub > 0.015);   // downstream of both TE points
    CHECK(g.row.omega == -2.0 * std::numbers::pi * 15000.0 / 60.0);  // the sign of the header: omega = -2 pi rpm / 60
    CHECK(g.row.s_row == 1);
    CHECK(g.row.sections.size() == d.layers.size());
    for (const auto& s : g.row.sections) {
        CHECK(s.kappa1 > s.kappa2);  // a compressor rotor in book signs
        CHECK(s.chord > 0.0);
        CHECK(std::abs(s.camber - (s.kappa1 - s.kappa2)) <= 1e-15);
    }
    std::printf("[throughflow-bg] default design: %zu quasi-normals, z = %.4f .. %.4f m; mid section kappa1 %.2f kappa2 %.2f "
                "stagger %.2f deg, chord %.4f m, t/c %.4f\n", g.quasi_normals.size(), g.quasi_normals.front().z_hub,
                g.quasi_normals.back().z_hub, g.row.sections[2].kappa1 * 180 / std::numbers::pi,
                g.row.sections[2].kappa2 * 180 / std::numbers::pi, g.row.sections[2].stagger * 180 / std::numbers::pi,
                g.row.sections[2].chord, g.row.sections[2].thickness_ratio);
    // the equal-parameter layout
    opt.layout = TF::BladeGenThroughflowOptions::Layout::EqualParameter;
    const auto ge = TF::throughflow_geometry(d, ev, opt);
    CHECK(ge.quasi_normals.size() == 10u);
    CHECK(ge.quasi_normals[4].z_hub < ge.quasi_normals[4].z_shroud);  // leans like the LE line
}

TEST_CASE("throughflow bladegen: the model levels on the curved default channel (no blade row)", "[turbolenza][throughflow]")
{
    D::Evaluator ev;
    TF::BladeGenThroughflowOptions opt;
    const auto g = TF::throughflow_geometry(D::BladeDesign::defaults(), ev, opt);
    std::vector<TF::ThroughflowResult> res;
    for (auto m : {TF::ThroughflowModel::FullNormalEquilibrium, TF::ThroughflowModel::Approximate,
                   TF::ThroughflowModel::SimpleNonIsentropic, TF::ThroughflowModel::SimpleIsentropic}) {
        const auto r = TF::solve_throughflow(channel_case(g, m, 3.0));
        std::printf("[throughflow-bg] channel %s: %s, %zu passes, errors", TF::to_string(m).c_str(), TF::to_string(r.status).c_str(),
                    r.history.size());
        for (std::size_t k = 0; k < r.history.size(); k += std::max<std::size_t>(1, r.history.size() / 10))
            std::printf(" %.2e(F %.3f, B* %.2f)", r.history[k].position_error, r.history[k].min_damping, r.history[k].bstar);
        std::printf("\n");
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        CHECK(max_mass_error(r) <= 1e-11);
        res.push_back(r);
    }
    // at the quasi-normal before the LE (the channel's curvature is largest around it)
    const std::size_t i = static_cast<std::size_t>(g.qn_before_le);
    auto diff = [&](std::size_t a, std::size_t b) {
        double d = 0.0;
        for (std::size_t k = 0; k < res[0].stations.size(); ++k)
            for (std::size_t j = 0; j < res[a].stations[k].nodes.size(); ++j)
                d = std::max(d, std::abs(res[a].stations[k].nodes[j].Wm - res[b].stations[k].nodes[j].Wm) /
                                    res[b].stations[k].nodes[j].Wm);
        return d;
    };
    const auto& full = res[0].stations[i];
    std::printf("[throughflow-bg] curved channel, quasi-normal %zu: W_m hub / shroud: full %.3f / %.3f, approximate %.3f / %.3f, "
                "simple %.3f / %.3f; kappa_m hub %.3f shroud %.3f 1/m (full), %.3f %.3f (approx.)\n", i,
                full.nodes.front().Wm, full.nodes.back().Wm, res[1].stations[i].nodes.front().Wm,
                res[1].stations[i].nodes.back().Wm, res[2].stations[i].nodes.front().Wm, res[2].stations[i].nodes.back().Wm,
                full.nodes.front().kappa_m, full.nodes.back().kappa_m, res[1].stations[i].nodes.front().kappa_m,
                res[1].stations[i].nodes.back().kappa_m);
    const double d01 = diff(1, 0), d02 = diff(2, 0), d23 = diff(3, 2);
    std::printf("[throughflow-bg] curved channel: max |dW_m| / W_m: approximate vs full %.3e, simple vs full %.3e, "
                "isentropic vs non-isentropic %.3e; full model %zu passes\n", d01, d02, d23, res[0].history.size());
    CHECK(d01 >= 1e-3);   // they differ (measured, see the report)
    // every station lies on the cylindrical extensions (Bell & Mehta walls: phi = kappa = 0 there), so the
    // approximate model sees no curvature and equals the simple one on this channel
    CHECK(d02 >= d01);
    CHECK(d23 <= 1e-12);  // uniform inlet: no entropy gradient, the two simple models agree
    // the full model: the surface tubes' Eq. 12-39 form against the Eq. 7-32 area term (first order at
    // the walls, where the half tube's centroid is off the wall), and the sheets' b against K_B b_geometric
    double ew = 0.0, ei = 0.0, es = 0.0;
    for (const auto& s : res[0].stations) {
        for (std::size_t j = 0; j < s.nodes.size(); ++j) {
            const double e = std::abs(s.nodes[j].b_mass / s.nodes[j].b - 1.0);
            double& slot = (j == 0 || j + 1 == s.nodes.size()) ? ew : ei;
            slot = std::max(slot, e);
        }
        for (const auto& x : s.sheets) es = std::max(es, std::abs(x.b / (s.K_B_endwall * x.b_geometric) - 1.0));
    }
    std::printf("[throughflow-bg] full model: |b_mass / b - 1| walls %.3e, interior %.3e; sheets |b / b_geometric - 1| %.3e\n",
                ew, ei, es);
    CHECK(ew <= 0.08);   // measured 3.9e-2
    CHECK(ei <= 0.03);   // measured 1.1e-2
    CHECK(es <= 0.01);   // measured 7.2e-4
}

TEST_CASE("throughflow bladegen: a rotating SpecifiedBladeRow on the default design converges", "[turbolenza][throughflow]")
{
    D::Evaluator ev;
    TF::BladeGenThroughflowOptions opt;
    opt.rpm = 15000.0;
    auto g = TF::throughflow_geometry(D::BladeDesign::defaults(), ev, opt);
    // the discharge angle: the layers' exit metal angle (zero deviation) and 5 % loss (ours)
    std::vector<double> span, beta2, loss;
    for (const auto& s : g.row.sections) {
        span.push_back(s.span);
        beta2.push_back(s.kappa2);
        loss.push_back(0.05);
    }
    g.row.model = std::make_shared<TF::SpecifiedBladeRow>(span, beta2, loss);
    for (auto m : {TF::ThroughflowModel::FullNormalEquilibrium, TF::ThroughflowModel::Approximate,
                   TF::ThroughflowModel::SimpleNonIsentropic}) {
        auto c = channel_case(g, m, 3.0);
        c.rows = {g.row};
        const auto r = TF::solve_throughflow(c);
        INFO(r.message);
        REQUIRE(r.status == TF::ThroughflowResult::Status::Converged);
        CHECK(max_mass_error(r) <= 1e-11);
        const auto& in = r.stations[static_cast<std::size_t>(g.qn_before_le)];
        const auto& out = r.stations[static_cast<std::size_t>(g.qn_after_te)];
        const auto& mid_in = in.nodes[5];
        const auto& mid_out = out.nodes[5];
        // the rotor does work (Euler): H2 - H1 = omega (r C_theta)_2 - omega (r C_theta)_1 > 0
        double dH = 0.0, prs = 0.0;
        for (std::size_t j = 0; j < in.nodes.size(); ++j) {
            CHECK(std::abs((out.nodes[j].H - in.nodes[j].H) - g.row.omega * (out.nodes[j].rCtheta - in.nodes[j].rCtheta)) <=
                  1e-6 * std::abs(out.nodes[j].H - in.nodes[j].H) + 1e-6);
            dH = std::max(dH, out.nodes[j].H - in.nodes[j].H);
            prs = std::max(prs, out.nodes[j].Pt / in.nodes[j].Pt);
            CHECK(out.nodes[j].H > in.nodes[j].H);
            CHECK(out.nodes[j].s > in.nodes[j].s);  // 5 % loss
        }
        const double inlet_beta = mid_in.Wm > 0.0 ? std::atan2(mid_in.Ctheta - mid_in.r * g.row.omega, mid_in.Wm) : 0.0;
        std::printf("[throughflow-bg] rotor %s: %zu passes, mid-span beta'1 %.2f deg (kappa1 %.2f), beta'2 %.2f deg, "
                    "W_m1 %.2f W_m2 %.2f m/s, max dH %.1f J/kg, max Pt2/Pt1 %.4f\n", TF::to_string(m).c_str(),
                    r.history.size(), inlet_beta * 180 / std::numbers::pi, g.row.sections[2].kappa1 * 180 / std::numbers::pi,
                    mid_out.beta * 180 / std::numbers::pi, mid_in.Wm, mid_out.Wm, dH, prs);
        CHECK(inlet_beta > 0.0);  // the header's sign: rpm > 0 gives beta'_1 > 0
        if (m == TF::ThroughflowModel::FullNormalEquilibrium) {
            std::printf("[throughflow-bg] rotor full model, position error by pass:");
            for (const auto& h : r.history) std::printf(" %.1e", h.position_error);
            std::printf("\n");
            // b(m) of the mid stream surface for the blade-to-blade solver
            const auto& ss = r.stream_surfaces[5];
            std::printf("[throughflow-bg] mid stream surface: m / mm, r / mm, b / mm:");
            for (std::size_t k = 0; k < ss.m.size(); ++k) std::printf(" (%.2f %.2f %.3f)", 1e3 * ss.m[k], 1e3 * ss.r[k], 1e3 * ss.b[k]);
            std::printf("\n");
        }
    }
}
