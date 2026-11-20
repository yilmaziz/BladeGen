// Tests of fit/InverseFit: the published NACA 65-(12)10 table (as fitted by blade2d).
#include "fit/InverseFit.hpp"
#include "profile/Quality.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using namespace BladeGen;
using namespace BladeGen::Fit;
using Catch::Approx;
using Math::Vec2;

namespace {
void read_naca(const std::string& file, std::vector<Vec2>& ss, std::vector<Vec2>& ps)
{
    std::ifstream f(std::string(BLADEGEN_SOURCE_DIR) + "/tests/data/" + file);
    REQUIRE(f.good());
    std::string line;
    bool header = false;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (!header) { header = true; continue; }
        std::stringstream s(line);
        std::string c[4];
        for (auto& x : c) std::getline(s, x, ',');
        auto add = [](std::vector<Vec2>& v, const std::string& a, const std::string& b) {
            if (a.empty() || b.empty() || a.find("nan") != std::string::npos) return;
            const Vec2 q(std::stod(a), std::stod(b));
            if (v.empty() || (q - v.back()).norm() > 0.0) v.push_back(q);
        };
        add(ss, c[0], c[1]);
        add(ps, c[2], c[3]);
    }
}

/// blade2d naca65/fit_naca65_1210.py search box
std::vector<FitParam> naca_spec(const Profile::ProfileParams& s)
{
    std::vector<FitParam> v = {
        {"le.radius", 0, 0.0015, 0.025}, {"le.stretch_ab", 0, 0.3, 4.0}, {"le.asym", 0, -0.008, 0.008},
        {"le.beta_deg", 0, 10.0, 60.0}, {"le.d_beta_ss_deg", 0, 5.0, 75.0}, {"le.d_beta_ps_deg", 0, 5.0, 80.0},
        {"te.radius", 0, 0.0001, 0.012}, {"te.stretch_ab", 0, 0.2, 4.0}, {"te.asym", 0, -0.005, 0.005},
        {"te.beta_deg", 0, 5.0, 40.0}, {"te.d_beta_ss_deg", 0, 0.0, 25.0}, {"te.d_beta_ps_deg", 0, 0.0, 25.0},
        {"ss.handle_le_frac", 0, 0.02, 0.6}, {"ss.ctrl2_mprime", 0, 0.1, 0.85}, {"ss.ctrl2_theta", 0, 0.04, 0.32},
        {"ss.handle_te_frac", 0, 0.04, 0.7}, {"thickness.t_max", 0, 0.05, 0.16}, {"thickness.mprime_max", 0, 0.2, 0.6},
        {"thickness.t_le", 0, 0.002, 0.08}, {"thickness.t_te", 0, 0.0005, 0.02}, {"thickness.shape_le", 0, 0.2, 3.0},
        {"thickness.shape_te", 0, 0.3, 4.0}, {"ps.handle_le_frac", 0, 0.002, 0.5}, {"ps.handle_te_frac", 0, 0.005, 0.5},
    };
    for (auto& p : v) p.x0 = std::clamp(Profile::get_param(s, p.name), p.lo, p.hi);
    return v;
}
}  // namespace

TEST_CASE("side distances vanish on the profile's own points", "[fit]")
{
    const Profile::Profile2D p(Profile::ProfileParams::naca65_1210());
    const auto ss = p.sample_side(Profile::Side::Suction, 20, 60).xy, ps = p.sample_side(Profile::Side::Pressure, 20, 60).xy;
    std::vector<double> dss, dps;
    side_distances(p, ss, ps, dss, dps, 600, 2000);
    for (double d : dss) REQUIRE(d < 2e-6);
    for (double d : dps) REQUIRE(d < 2e-6);
}

TEST_CASE("normalize_contour finds LE / TE, sides and the frame", "[fit]")
{
    // a known profile, rotated by 30 deg, scaled by 2 and shifted
    const Profile::Profile2D p(Profile::ProfileParams::naca65_1210());
    const auto loop = p.sample_loop(30, 100).xy;
    const double g = Math::deg2rad(30.0);
    std::vector<Vec2> in;
    for (const auto& q : loop) in.emplace_back(1.0 + 2.0 * (std::cos(g) * q.x() - std::sin(g) * q.y()),
                                               -0.5 + 2.0 * (std::sin(g) * q.x() + std::cos(g) * q.y()));
    const auto n = normalize_contour(in);
    // the longest chord of the contour is not exactly the model chord P_LE -> P_TE
    REQUIRE(n.stagger == Approx(g).margin(6e-3));
    REQUIRE(n.chord == Approx(2.0).epsilon(2e-3));
    REQUIRE_FALSE(n.mirrored);
    REQUIRE(n.ss[n.ss.size() / 2].y() > n.ps[n.ps.size() / 2].y());
    REQUIRE((n.to_input(n.ss.front()) - n.le).norm() < 1e-12);
    // the mirror image is detected
    for (auto& q : in) q.y() = -q.y();
    REQUIRE(normalize_contour(in).mirrored);
}

TEST_CASE("inverse fit of the published NACA 65-(12)10 coordinates", "[fit][slow]")
{
    std::vector<Vec2> ss, ps;
    read_naca("naca65_1210_XY_dense241.csv", ss, ps);
    REQUIRE(ss.size() > 200);
    for (auto mode : {Profile::ConnectorMode::G2Legacy, Profile::ConnectorMode::Fair}) {
        auto start = mode == Profile::ConnectorMode::G2Legacy ? Profile::ProfileParams::naca65_1210_blade2d()
                                                              : Profile::ProfileParams::naca65_1210();
        // start away from the answer: perturb the shape
        start.le.radius *= 1.3;
        start.thickness.t_max *= 0.9;
        start.ss.ctrl2_theta *= 1.1;
        FitOptions opt;
        opt.params = naca_spec(start);
        opt.stages = 2;
        const auto r = fit_profile(ss, ps, start, opt);
        INFO(Profile::to_string(mode) << ": RMS " << 100 * r.rms << " % chord, max " << 100 * r.max << " %, "
                                      << r.evaluations << " evaluations, " << r.seconds << " s");
        // blade2d: RMS 0.0095 %, max 0.04 % chord on this target
        REQUIRE(r.rms < 2.0e-4);
        REQUIRE(r.max < 1.0e-3);
    }
}
