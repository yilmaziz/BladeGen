// Tests of profile/Quality: good profiles pass, bad ones are flagged but still built.
#include "profile/Quality.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>

using namespace BladeGen::Profile;
using BladeGen::Math::Vec2;

namespace {
bool has(const QualityReport& r, const std::string& code)
{
    return std::any_of(r.issues.begin(), r.issues.end(), [&](const Issue& i) { return i.code == code; });
}
}  // namespace

TEST_CASE("the default NACA 65 profile (fairness-weighted fit) has no geometry errors", "[quality]")
{
    const auto r = check_quality(Profile2D(ProfileParams::naca65_1210()));
    INFO(r.summary());
    REQUIRE(r.ok());
}

TEST_CASE("a fair profile has no geometry errors", "[quality]")
{
    auto prm = ProfileParams::naca65_1210();
    prm.te = {0.003, 1.2, 0.0, 10.0, 6.4, 4.0, 1.0};
    prm.ps.handle_le_frac = 0.06;
    const Profile2D p(prm);
    const auto r = check_quality(p);
    INFO(r.summary());
    REQUIRE(r.ok());
    REQUIRE(r.self_intersections == 0);
    REQUIRE(r.min_thickness > 0.0);
    REQUIRE(r.ss.kappa.size() == r.ss.xy.size());
}

TEST_CASE("the published blade2d NACA 65 fit has an S-shaped TE connector (found by the curvature check)", "[quality]")
{
    // blade2d's fit put the TE radius on its lower bound (0.000255 chord) with a zero
    // PS wedge; the connector then dips to a concave curvature before the TE point.
    for (auto mode : {ConnectorMode::Fair, ConnectorMode::C2Parametric, ConnectorMode::G2Legacy}) {
        auto prm = ProfileParams::naca65_1210_blade2d();
        prm.connector = mode;
        const auto r = check_quality(Profile2D(prm));
        INFO(r.summary());
        REQUIRE(has(r, "edge_inflection"));
    }
}

TEST_CASE("a suction side pulled below its chord is flagged but still built", "[quality]")
{
    auto prm = ProfileParams::naca65_1210();
    prm.ss.ctrl2_theta = -0.05;  // concave SS, PS crosses it
    const Profile2D p(prm);      // no throw: the geometry is generated
    const auto r = check_quality(p);
    INFO(r.summary());
    REQUIRE_FALSE(r.ok());
    REQUIRE((has(r, "ss_inflection") || has(r, "self_intersection") || has(r, "negative_thickness")));
}

TEST_CASE("an over-tensioned LE connector loops and is reported as an error", "[quality]")
{
    auto prm = ProfileParams::naca65_1210();
    prm.le.tension = 6.0;
    const Profile2D p(prm);
    const auto r = check_quality(p);
    INFO(r.summary());
    REQUIRE((has(r, "edge_inflection") || has(r, "self_intersection")));
    REQUIRE(r.count(Severity::Error) > 0);
}

TEST_CASE("self-intersection counter on simple polygons", "[quality]")
{
    const std::vector<Vec2> square = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}};
    REQUIRE(count_self_intersections(square) == 0);
    const std::vector<Vec2> bow = {{0, 0}, {1, 1}, {1, 0}, {0, 1}, {0, 0}};
    REQUIRE(count_self_intersections(bow) == 1);
}
