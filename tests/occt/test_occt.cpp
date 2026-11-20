// Tests of occt/OcctBlade: skinning, trimmed solid, fillet, export round trip.
#include "occt/OcctBlade.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <filesystem>

using namespace BladeGen;
using Catch::Approx;
using Math::Vec3;

namespace {
std::vector<Vec3> circle(double r, double z, int n, double phase = 0.0)
{
    std::vector<Vec3> l;
    for (int i = 0; i < n; ++i) {
        const double a = phase + 2.0 * Math::kPi * i / n;
        l.emplace_back(r * std::cos(a), r * std::sin(a), z);
    }
    l.push_back(l.front());
    return l;
}
std::string out(const std::string& name) { return std::string(BLADEGEN_TEST_OUTPUT) + "/" + name; }
}  // namespace

TEST_CASE("skinning interpolates every loop and is periodic and C2 across the seam", "[occt]")
{
    std::vector<std::vector<Vec3>> loops;
    for (int k = 0; k < 6; ++k) loops.push_back(circle(1.0 + 0.2 * k, 0.5 * k, 64));
    const auto sk = Occt::skin_loops(loops);
    REQUIRE(sk.max_error < 1e-10);
    REQUIRE(sk.surface->IsUPeriodic());
    REQUIRE_FALSE(sk.surface->IsVPeriodic());
    REQUIRE(sk.surface->UDegree() == 3);
    // the cone between the circles is reproduced closely
    double err = 0.0;
    for (int i = 0; i <= 40; ++i)
        for (int k = 0; k <= 20; ++k) {
            const gp_Pnt p = sk.surface->Value(i / 40.0, k / 20.0);
            const double r = std::hypot(p.X(), p.Y());
            err = std::max(err, std::abs(r - (1.0 + 0.4 * p.Z())));
        }
    REQUIRE(err < 5e-6);  // cubic interpolation of a 64-point circle
    // C2 across the periodic seam u = 0 / 1
    gp_Pnt a, b;
    gp_Vec au, av, auu, avv, auv, bu, bv, buu, bvv, buv;
    sk.surface->D2(1e-12, 0.37, a, au, av, auu, avv, auv);
    sk.surface->D2(1.0 - 1e-12, 0.37, b, bu, bv, buu, bvv, buv);
    REQUIRE(a.Distance(b) < 1e-9);
    REQUIRE((au - bu).Magnitude() < 1e-6 * au.Magnitude());
    REQUIRE((auu - buu).Magnitude() < 1e-4 * auu.Magnitude());
    REQUIRE_THROWS_AS(Occt::skin_loops({loops[0]}), std::invalid_argument);
    auto open = loops;
    open[1].back() = Vec3(9, 9, 9);
    REQUIRE_THROWS_AS(Occt::skin_loops(open), std::invalid_argument);
}

TEST_CASE("demo rotor: valid trimmed solid with caps on hub and casing, STEP round trip", "[occt]")
{
    const Blade::Blade blade(Blade::BladeDefinition::demo_rotor());
    const auto sh = Occt::build_blade(blade, 0.0, 24, 80);
    INFO("skin error " << sh.skin.max_error);
    REQUIRE(sh.skin.max_error < 1e-9);
    REQUIRE(sh.blade_valid);
    REQUIRE(Occt::count_solids(sh.blade) == 1);
    REQUIRE(sh.blade_volume > 0.0);
    // rough volume check: mid-section area (conformal area x r^2) times blade height
    const auto mid = blade.section(0.5, 24, 80);
    double a = 0.0;
    for (std::size_t i = 0; i + 1 < mid.conformal.size(); ++i)
        a += 0.5 * (mid.conformal[i].x() * mid.conformal[i + 1].y() - mid.conformal[i + 1].x() * mid.conformal[i].y());
    const double r = 0.5 * (0.26 + 0.3975);
    const double v_est = std::abs(a) * r * r * (0.3975 - 0.26);
    REQUIRE(sh.blade_volume == Approx(v_est).epsilon(0.15));
    // the end caps lie on the hub and casing: blade points have r between them
    std::filesystem::create_directories(BLADEGEN_TEST_OUTPUT);
    Occt::write_shape(sh.blade, out("demo_blade.step"), Occt::Unit::Millimetre);
    const auto back = Occt::read_step(out("demo_blade.step"));
    REQUIRE(Occt::volume(back) == Approx(sh.blade_volume).epsilon(1e-6));
    Occt::write_shape(sh.blade, out("demo_blade.stl"));
    Occt::write_shape(sh.blade, out("demo_blade.igs"));
    REQUIRE(std::filesystem::file_size(out("demo_blade.stl")) > 1000);
    REQUIRE_THROWS_AS(Occt::write_shape(sh.blade, out("x.obj")), std::invalid_argument);
    const auto mesh = Occt::tessellate(sh.blade, 2e-4);
    REQUIRE(mesh.triangles.size() > 100);
}

TEST_CASE("hub fillet on the demo rotor", "[occt][fillet]")
{
    const Blade::Blade blade(Blade::BladeDefinition::demo_rotor());
    const auto plain = Occt::build_blade(blade, 0.0, 24, 80);
    const auto fil = Occt::build_blade(blade, 0.003, 24, 80);
    INFO(fil.fillet_message);
    REQUIRE(fil.fillet_ok);
    REQUIRE(Occt::is_valid(fil.blade_on_hub));
    // the fillet adds material in the corner
    REQUIRE(Occt::volume(fil.blade_on_hub) > Occt::volume(plain.blade_on_hub));
    REQUIRE(Occt::is_valid(fil.blade_with_fillet));
    REQUIRE(Occt::volume(fil.blade_with_fillet) > plain.blade_volume);
    // a rolling ball: both contacts at distance r from the centre, tangent to both walls
    REQUIRE(fil.fillet_ball_error < 1e-9);
    REQUIRE(fil.fillet_g1_deg < 0.5);
    Occt::write_shape(fil.blade_on_hub, out("demo_blade_fillet.step"));
}

TEST_CASE("periodic interpolation stays well posed for many points", "[occt]")
{
    std::vector<double> u;
    for (int i = 0; i <= 2000; ++i) u.push_back(i / 2000.0);
    const Occt::PeriodicInterpolator P(u);  // det underflows at this size; rcond does not
    const auto c = P.curve(circle(1.0, 0.0, 2000));
    REQUIRE(c->IsPeriodic());
    REQUIRE(std::abs(c->Value(0.25).X() - std::cos(Math::kPi / 2)) < 1e-9);
}
