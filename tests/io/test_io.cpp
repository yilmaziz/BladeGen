// Tests of io/Rotor67 (against the published tables and the Python construction),
// io/BladeJson (round trips, blade2d import) and the Rotor 67 section fit.
#include "fit/SectionFit.hpp"
#include "io/BladeJson.hpp"
#include "io/Rotor67.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>

using namespace BladeGen;
using Catch::Approx;
using Math::Vec2;
using Math::Vec3;

namespace {
std::string src(const std::string& p) { return std::string(BLADEGEN_SOURCE_DIR) + "/" + p; }
const IO::Rotor67Data& rotor67()
{
    static const auto d = IO::read_rotor67(src("data/rotor67/Rotor_14_Sections_and_Flowpath_EN.txt"));
    return d;
}
}  // namespace

TEST_CASE("Rotor 67 tables: sections, flowpath and notes", "[rotor67]")
{
    const auto& d = rotor67();
    REQUIRE(d.sections.size() == 14);
    REQUIRE(d.n_blades == 22);
    REQUIRE(d.tip_clearance == Approx(0.001016));
    REQUIRE(d.hub.size() == 80);
    REQUIRE(d.rotating_hub[0] == Approx(-0.01374));
    for (const auto& s : d.sections) {
        REQUIRE(s.z.size() == 35);
        // the edge-circle centres lie on the stream line at their meridional distance
        for (const auto* c : {&s.le, &s.te}) {
            const Vec2 p = s.streamline.at_m(c->mc);
            // (the published centres are within ~10 um of the stream line: table precision)
            REQUIRE(std::hypot(p.x() - c->zc, p.y() - c->rc) < 2e-5);
        }
        // surfaces meet the tabulated end slopes: tan(beta) = r dtheta/dm
        const double sl = s.surf1.derivs(0.0, 1)[1][0];
        REQUIRE(std::atan(s.r.front() * sl) == Approx(s.beta1.front()).margin(1e-12));
        // tangency points lie on the circles
        const double mt = s.m_tan[0][0];
        REQUIRE(std::hypot(mt - s.le.mc, s.le.rc * (s.theta(1, mt) - s.le.thc)) == Approx(s.le.rho).epsilon(1e-4));  // radius tabulated to 5 digits
    }
    REQUIRE_THROWS(IO::read_rotor67(src("no/such/file.txt")));
}

TEST_CASE("Rotor 67 published sections agree with the Python construction (rotor67_sections.py)", "[rotor67]")
{
    std::ifstream f(src("tests/data/rotor67_sections_xyz_subset.csv"));
    REQUIRE(f.good());
    std::map<int, std::vector<Vec3>> py;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream s(line);
        std::string c[5];
        for (auto& x : c) std::getline(s, x, ',');
        py[std::stoi(c[0])].emplace_back(std::stod(c[2]), std::stod(c[3]), std::stod(c[4]));
    }
    REQUIRE(py.size() == 3);
    for (const auto& [idx, pts] : py) {
        const auto ours = rotor67().sections[static_cast<std::size_t>(idx - 1)].profile_xyz(161, 33);
        REQUIRE(ours.size() == pts.size());
        double dmax = 0.0;
        for (std::size_t i = 0; i < pts.size(); ++i) dmax = std::max(dmax, (ours[i] - pts[i]).norm());
        INFO("section " << idx << ": max difference " << 1e6 * dmax << " um");
        REQUIRE(dmax < 2e-6);  // the CSV is rounded to 1e-8 m; spline variants differ by < 1 um
    }
}

TEST_CASE("blade JSON round trip and blade2d profile import", "[json]")
{
    const auto def = Blade::BladeDefinition::demo_rotor();
    const std::string path = std::string(BLADEGEN_TEST_OUTPUT) + "/demo_rotor.json";
    std::filesystem::create_directories(BLADEGEN_TEST_OUTPUT);
    IO::write_blade_json(def, path);
    const auto back = IO::read_blade_json(path);
    REQUIRE(back.stations.size() == def.stations.size());
    const Blade::Blade a(def), b(back);
    const auto sa = a.section(0.37, 10, 30), sb = b.section(0.37, 10, 30);
    for (std::size_t i = 0; i < sa.xyz.size(); ++i) REQUIRE((sa.xyz[i] - sb.xyz[i]).norm() < 1e-14);
    // blade2d's flat blade_params.json -> the same profile as the built-in preset
    const auto p = IO::read_profile_json(src("tests/data/naca65_1210_blade_params.json"));
    const auto q = Profile::ProfileParams::naca65_1210_blade2d();
    REQUIRE(p.connector == Profile::ConnectorMode::G2Legacy);
    for (const auto& i : Profile::param_table()) REQUIRE(Profile::get_param(p, i.name) == Approx(Profile::get_param(q, i.name)));
    REQUIRE_THROWS_AS(IO::profile_from_json(nlohmann::json{{"le.radius", -1.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS(IO::profile_from_json(nlohmann::json{{"no.such", 1.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS(IO::blade_from_json(nlohmann::json{{"name", "x"}}), std::invalid_argument);
}

TEST_CASE("Rotor 67 section 7 is reproduced by the inverse fit to a few micrometres", "[rotor67][slow]")
{
    const auto& s = rotor67().sections[6];
    Fit::FitOptions opt;
    opt.fairness_weight = 0.3;
    const auto r = Fit::fit_section(s.profile_conformal(161, 33), s.streamline, 0.5, Profile::ProfileParams::cda_start(), opt);
    INFO("RMS " << 1e6 * r.rms_3d << " um, max " << 1e6 * r.max_3d << " um, " << r.fit.seconds << " s");
    REQUIRE(r.rms_3d < 5e-6);
    REQUIRE(r.max_3d < 2e-5);
    // the placement reproduces the frame: LE point of the fitted section = published LE
    const auto map = Blade::make_section_map(s.streamline, r.station.placement);
    REQUIRE((map.conformal(Vec2(0, 0)) - r.frame.le).norm() < 1e-12);
    REQUIRE((map.conformal(Vec2(1, 0)) - r.frame.te).norm() < 1e-12);
}
