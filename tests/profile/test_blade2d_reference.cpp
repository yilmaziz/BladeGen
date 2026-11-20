// Port check: ConnectorMode::G2Legacy reproduces the blade2d Python engine.
// Reference: tests/data/blade2d_naca65_reference.csv from
// tools/reference/make_blade2d_reference.py (blade2d profile_c2_connect.py).
// blade2d takes arc length and normals from 2000 samples; BladeGen integrates
// them exactly, so the PS control points differ at the 1e-6 level.
#include "profile/Profile2D.hpp"

#include <catch2/catch_test_macros.hpp>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace BladeGen::Profile;
using BladeGen::Math::Vec2;

namespace {
std::map<std::string, std::vector<Vec2>> load()
{
    std::ifstream f(std::string(BLADEGEN_SOURCE_DIR) + "/tests/data/blade2d_naca65_reference.csv");
    REQUIRE(f.good());
    std::map<std::string, std::vector<Vec2>> m;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("kind", 0) == 0) continue;
        std::stringstream ss(line);
        std::string kind, i, x, y;
        std::getline(ss, kind, ',');
        std::getline(ss, i, ',');
        std::getline(ss, x, ',');
        std::getline(ss, y, ',');
        m[kind].emplace_back(std::stod(x), std::stod(y));
    }
    return m;
}
}  // namespace

TEST_CASE("G2Legacy reproduces blade2d's NACA 65-(12)10 profile", "[profile][reference]")
{
    const auto ref = load();
    const Profile2D p(ProfileParams::naca65_1210_blade2d());
    // suction side: identical construction, identical control points
    const auto& ssc = p.suction_side().ctrl();
    REQUIRE(ref.at("ss_ctrl").size() == ssc.size());
    for (std::size_t i = 0; i < ssc.size(); ++i) REQUIRE((ref.at("ss_ctrl")[i] - ssc[i]).norm() < 1e-12);
    // pressure side control points: blade2d's sampled arc length vs exact
    const auto& psc = p.pressure_side().ctrl();
    REQUIRE(ref.at("ps_ctrl").size() == psc.size());
    double dmax = 0.0;
    for (std::size_t i = 0; i < psc.size(); ++i) dmax = std::max(dmax, (ref.at("ps_ctrl")[i] - psc[i]).norm());
    INFO("max PS control point difference " << dmax);
    REQUIRE(dmax < 2e-6);
    // curves at blade2d's sample parameters
    const int N = 2000;
    for (std::size_t k = 0; k < ref.at("ss").size(); ++k) {
        const double u = static_cast<double>(20 * k) / (N - 1);
        REQUIRE((p.suction_side().eval(u) - ref.at("ss")[k]).norm() < 1e-12);
        REQUIRE((p.pressure_side().eval(u) - ref.at("ps")[k]).norm() < 2e-6);
    }
    // LE / TE connectors (201 uniform samples in blade2d)
    double dle = 0.0, dte = 0.0;
    for (std::size_t k = 0; k < ref.at("le").size(); ++k) {
        const double u = static_cast<double>(k) / 200.0;
        dle = std::max(dle, (p.le_connector().eval(u) - ref.at("le")[k]).norm());
        dte = std::max(dte, (p.te_connector().eval(u) - ref.at("te")[k]).norm());
    }
    INFO("LE " << dle << " TE " << dte);
    REQUIRE(dle < 2e-6);
    REQUIRE(dte < 2e-6);
}
