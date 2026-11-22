// Tests of design/BladeDesign: five layers on the meridional channel, conformal placement,
// metal angle, 3-D mesh on hub and casing, JSON (v2 and the v1.0.1 format), interpolation.
#include "design/BladeDesign.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Design;
using Catch::Approx;
using Math::Vec2;
using Math::Vec3;

namespace {
double dist_to(const std::vector<Vec2>& poly, const Vec2& q)
{
    double d = 1e300;
    for (std::size_t i = 1; i < poly.size(); ++i) {
        const Vec2 a = poly[i - 1], ab = poly[i] - a;
        const double t = std::clamp((q - a).dot(ab) / ab.squaredNorm(), 0.0, 1.0);
        d = std::min(d, (q - (a + t * ab)).norm());
    }
    return d;
}
Vec2 xr_of(const Vec3& p) { return {p.x(), std::hypot(p.y(), p.z())}; }
}  // namespace

TEST_CASE("Natural spline weights: interpolation, partition of unity, linear exactness", "[design]")
{
    const std::vector<double> s = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (std::size_t k = 0; k < s.size(); ++k) {
        const auto w = spline_weights(s, s[k]);
        for (std::size_t i = 0; i < s.size(); ++i) CHECK(w[i] == Approx(i == k ? 1.0 : 0.0).margin(1e-14));
    }
    for (double q : {0.1, 0.4, 0.93}) {
        const auto w = spline_weights(s, q);
        double sum = 0.0, lin = 0.0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            sum += w[i];
            lin += w[i] * (2.0 + 3.0 * s[i]);
        }
        CHECK(sum == Approx(1.0).margin(1e-14));
        CHECK(lin == Approx(2.0 + 3.0 * q).margin(1e-13));
    }
    CHECK(spline_weights({0.3}, 0.9)[0] == 1.0);
    CHECK(spline_weights({0.0, 1.0}, 0.25)[1] == Approx(0.25));
    CHECK_THROWS_AS(spline_weights({}, 0.5), std::invalid_argument);
    CHECK_THROWS_AS(spline_weights({0.0, 0.0, 1.0}, 0.5), std::invalid_argument);
}

TEST_CASE("Two-sided tanh rows: ends, monotone, symmetric, clustered", "[design]")
{
    const auto u = tanh_two_sided(21, 2.0);
    CHECK(u.front() == 0.0);
    CHECK(u.back() == 1.0);
    for (std::size_t i = 1; i < u.size(); ++i) CHECK(u[i] > u[i - 1]);
    for (std::size_t i = 0; i < u.size(); ++i) CHECK(u[i] + u[u.size() - 1 - i] == Approx(1.0).margin(1e-14));
    CHECK(u[1] < 0.05);  // first cell smaller than uniform (0.05)
    const auto v = tanh_two_sided(5, 0.0);
    CHECK(v[1] == Approx(0.25));
    CHECK_THROWS_AS(tanh_two_sided(1, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(tanh_two_sided(5, -1.0), std::invalid_argument);
}

TEST_CASE("Default design: five layers placed between the LE and TE lines, mesh on the walls", "[design]")
{
    const auto d = BladeDesign::defaults();
    REQUIRE(d.layers.size() == 5);
    CHECK(d.layers[1].span == 0.25);
    CHECK(d.layers[3].span == 0.75);
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    const Blade::FlowPath fp(d.meridional);
    const auto hub = fp.surface(0.0).polyline(2000), cas = fp.surface(1.0).polyline(2000);
    for (std::size_t k = 0; k < 5; ++k) {
        const auto& L = r.layers[k];
        CHECK((L.surface.le_point() - fp.le_point(d.layers[k].span)).norm() < 1e-12);
        CHECK(L.mfrac.front() == Approx(0.0).margin(1e-12));
        CHECK(L.mfrac.back() == Approx(1.0).margin(1e-9));
        // the camber starts at the LE point and ends at the TE point of the layer
        CHECK((L.camber_xr.front() - fp.le_point(d.layers[k].span)).norm() < 1e-9);
        CHECK((L.camber_xr.back() - fp.te_point(d.layers[k].span)).norm() < 1e-7);
        CHECK(L.beta_le > L.beta_te);  // compressor: the camber turns the flow toward the axis
        CHECK(L.section->quality.ok());
        CHECK(L.solidity_3d > 0.0);
    }
    // metal angle = stagger + camber slope; hub more cambered than the tip
    CHECK(r.layers[0].beta_le - r.layers[0].beta_te > r.layers[4].beta_le - r.layers[4].beta_te);
    // mesh rows: row 0 on the hub, the last row on the casing
    REQUIRE(r.blade.nj == d.mesh.n_rows);
    REQUIRE(r.blade.ni == static_cast<int>(r.layers[0].section->loop.xy.size()));
    for (int i = 0; i < r.blade.ni; i += 7) {
        CHECK(dist_to(hub, xr_of(r.blade.at(i, 0))) < 1e-5);  // reference polyline sag ~4e-6
        CHECK(dist_to(cas, xr_of(r.blade.at(i, r.blade.nj - 1))) < 1e-5);
    }
    for (const auto& p : r.hub_wall.p) CHECK(dist_to(hub, xr_of(p)) < 1e-5);
    CHECK(r.hub_wall.p.size() == static_cast<std::size_t>(d.mesh.n_wall_m * d.mesh.n_wall_theta));
    // metal-angle field between the layer values
    CHECK(r.beta.value.size() == static_cast<std::size_t>(r.beta.ni * r.beta.nj));
    CHECK(r.beta.hi > r.beta.lo);
    CHECK(r.beta.value.front() == Approx(r.layers[0].beta_le).margin(1e-9));
    const auto loops = section_loops(r);
    CHECK(loops.size() == static_cast<std::size_t>(r.blade.nj));
}

TEST_CASE("Metal angle on the stream surface equals beta from the conformal plane (3-D check)", "[design]")
{
    auto d = BladeDesign::defaults();
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    const auto& L = r.layers[2];
    for (std::size_t i = 10; i + 10 < L.mprime.size(); i += 15) {
        // 3-D camber tangent between neighbouring camber points vs meridional and circumferential directions
        const auto P = [&](std::size_t k) { return L.surface.cartesian(L.mprime[k], Math::deg2rad(L.theta_deg[k])); };
        const Vec3 t = P(i + 1) - P(i - 1);
        const double h = 1e-7, mp = L.mprime[i], th = Math::deg2rad(L.theta_deg[i]);
        const Vec3 em = (L.surface.cartesian(mp + h, th) - L.surface.cartesian(mp - h, th)).normalized();
        const Vec3 et = (L.surface.cartesian(mp, th + h) - L.surface.cartesian(mp, th - h)).normalized();
        CHECK(Math::rad2deg(std::atan2(t.dot(et), t.dot(em))) == Approx(L.beta_deg[i]).margin(0.3));
    }
}

TEST_CASE("Cylindrical walls: the layers are wrapped on cylinders", "[design]")
{
    auto d = BladeDesign::defaults();
    d.meridional.hub = {-15.0, 50.0, 15.0, 50.0};
    d.meridional.hub_law = Blade::WallLaw::Linear;
    d.meridional.casing = {-15.0, 80.0, 15.0, 80.0};  // both walls cylinders (the default casing is not)
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    for (std::size_t k = 0; k < 5; ++k) {
        const auto& L = r.layers[k];
        REQUIRE(L.surface.cylindrical());
        const double rad = 50.0 + 30.0 * d.layers[k].span;
        for (std::size_t i = 0; i < L.xyz.size(); i += 9) {
            CHECK(std::hypot(L.xyz[i].y(), L.xyz[i].z()) == Approx(rad).epsilon(1e-12));
            CHECK(L.xyz[i].x() == Approx(-15.0 + rad * L.conformal[i].x()).margin(1e-10));
        }
        // axial chord 30 mm: c' cos(gamma) = 30 / r
        CHECK(L.chord_mprime * std::cos(L.gamma) == Approx(30.0 / rad).epsilon(1e-12));
    }
}

TEST_CASE("Design JSON: round trip, v1.0.1 import, invalid input", "[design]")
{
    auto d = BladeDesign::defaults();
    d.layers[1].method = Method::Custom;
    d.layers[1].custom.le.radius = 0.011;
    d.meridional.casing_law = Blade::WallLaw::Linear;
    d.n_blades = 31;
    const auto j = to_json(d);
    const auto e = design_from_json(j);
    CHECK(to_json(e) == j);
    CHECK(e.layers[1].custom.le.radius == 0.011);

    nlohmann::json v1;
    v1["format"] = "bladegen-design-1";
    v1["stations"] = nlohmann::json::array();
    const double cl[] = {1.5, 1.2, 0.8};
    for (int k = 0; k < 3; ++k)
        v1["stations"].push_back({{"name", k == 0 ? "Hub" : k == 1 ? "Mid" : "Tip"}, {"span", 0.5 * k},
                                  {"method", "NACA 65"}, {"cl0", cl[k]}, {"t_max_over_c", 0.1},
                                  {"variant", "emery_base"}, {"stagger_deg", 30.0 + 10 * k}, {"solidity", 1.2}});
    const auto f = design_from_json(v1);
    REQUIRE(f.layers.size() == 5);
    CHECK(f.layers[1].naca.cl0 == Approx(1.35));
    CHECK(f.layers[3].stagger_deg == Approx(45.0));
    CHECK(f.layers[4].naca.variant == Naca::Variant::EmeryBase);

    auto bad = j;
    bad["format"] = "nonsense";
    CHECK_THROWS_AS(design_from_json(bad), std::invalid_argument);
    bad = j;
    bad["layers"][2]["span"] = 0.1;  // not increasing
    CHECK_THROWS_AS(design_from_json(bad), std::invalid_argument);
    bad = j;
    bad["layers"][0]["method"] = "DCA";
    CHECK_THROWS_AS(design_from_json(bad), std::invalid_argument);
    CHECK_THROWS(read_design("does/not/exist.json"));
}

TEST_CASE("Evaluator: invalid designs are reported, the cache reuses unchanged layers", "[design]")
{
    Evaluator ev;
    auto d = BladeDesign::defaults();
    d.meridional.casing.r_le = 10.0;
    auto r = ev.evaluate(d);
    CHECK_FALSE(r.ok());
    CHECK_FALSE(r.errors.empty());

    d = BladeDesign::defaults();
    d.layers[2].method = Method::Custom;
    d.layers[2].custom.le.radius = 0.0;
    r = ev.evaluate(d);
    CHECK_FALSE(r.ok());

    d = BladeDesign::defaults();
    const auto a = ev.evaluate(d);
    REQUIRE(a.ok());
    d.layers[4].stagger_deg += 2.0;  // placement only: every 2-D section from the cache
    const auto b = ev.evaluate(d);
    REQUIRE(b.ok());
    for (std::size_t k = 0; k < 5; ++k) CHECK(a.layers[k].section.get() == b.layers[k].section.get());
    CHECK(b.layers[4].beta_le == Approx(a.layers[4].beta_le + 2.0).margin(1e-9));
}

TEST_CASE("NACA 65 + Custom layer: fitted custom parameters close to the section", "[design][slow]")
{
    Evaluator ev;
    auto d = BladeDesign::defaults();
    d.layers[2].method = Method::Naca65Custom;
    fit_layer(d.layers[2], ev);
    REQUIRE(d.layers[2].fitted_to);
    CHECK(d.layers[2].fit_rms < 1e-3);
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    REQUIRE(r.layers[2].section->deviation_rms);
    CHECK(*r.layers[2].section->deviation_rms < 1.5e-3);
}

TEST_CASE("Rotor 67 layers: exact at the 14 fitted sections, blended between, invalid span rejected", "[design][rotor67]")
{
    const auto& t = rotor67_sections();
    REQUIRE(t.size() == 14);
    for (std::size_t k = 0; k < t.size(); ++k) {
        const auto l = rotor67_layer(t[k].span);
        CHECK(l.stagger_deg == Approx(t[k].stagger_deg).margin(1e-12));
        CHECK(Profile::to_vector(l.profile) == Profile::to_vector(t[k].profile));
        if (k > 0) CHECK(t[k].stagger_deg > t[k - 1].stagger_deg);  // stagger rises toward the tip
    }
    const double s = 0.5 * (t[3].span + t[4].span);
    const auto m = rotor67_layer(s);
    CHECK(m.stagger_deg == Approx(0.5 * (t[3].stagger_deg + t[4].stagger_deg)));
    CHECK(m.profile.le.radius == Approx(0.5 * (t[3].profile.le.radius + t[4].profile.le.radius)));
    CHECK_THROWS_AS(rotor67_layer(-0.1), std::invalid_argument);
    CHECK_THROWS_AS(rotor67_layer(1.1), std::invalid_argument);
    CHECK_THROWS_AS(rotor67_layer(std::nan("")), std::invalid_argument);
}

TEST_CASE("Default design: Rotor 67 custom layers, scaled by the meridional path", "[design][rotor67]")
{
    // the metal angles of the fitted sections agree with the published surface angles (TP-2879 hub:
    // surfaces 44 / 32 deg at the LE, -32 / -14 deg at the TE; tip: 68 / 63 and 49 / 59 deg)
    auto d = BladeDesign::defaults();
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    CHECK(r.layers[0].beta_le == Approx(38.0).margin(6.0));
    CHECK(r.layers[0].beta_te == Approx(-23.0).margin(6.0));
    CHECK(r.layers[4].beta_le == Approx(65.6).margin(6.0));
    CHECK(r.layers[4].beta_te == Approx(53.8).margin(6.0));
    for (std::size_t k = 0; k < 5; ++k) {
        CHECK(d.layers[k].method == Method::Custom);
        CHECK(d.layers[k].stagger_deg == Approx(rotor67_layer(d.layers[k].span).stagger_deg));
        CHECK(r.layers[k].section->quality.count(Profile::Severity::Error) == 0);
        // control polygons and the thickness law of a custom layer
        const auto& sec = *r.layers[k].section;
        REQUIRE(sec.control);
        REQUIRE(sec.control->ss.size() == 5);
        CHECK(sec.control->ss[2].x() == Approx(d.layers[k].custom.ss.ctrl2_mprime));
        CHECK(sec.control->ss[2].y() == Approx(d.layers[k].custom.ss.ctrl2_theta));
        CHECK(sec.control->ps.size() == 4 + d.layers[k].custom.ps.mprime_ctrl.size());
        REQUIRE(sec.thickness_law);
        const auto mx = *std::max_element(sec.thickness_law->law.begin(), sec.thickness_law->law.end());
        CHECK(mx == Approx(d.layers[k].custom.thickness.t_max).epsilon(1e-2));  // 120 law stations
    }
    // the chord follows the channel: twice the axial length, twice the chord (tip layer on a cylinder)
    auto c = d;
    c.meridional.hub.z_le = c.meridional.casing.z_le = -15.0;
    c.meridional.hub.z_te = c.meridional.casing.z_te = 15.0;
    c.meridional.casing.r_le = c.meridional.casing.r_te = 80.0;
    const auto rc = ev.evaluate(c);
    REQUIRE(rc.ok());
    REQUIRE(rc.layers[4].surface.cylindrical());
    auto e = c;
    e.meridional.hub.z_le = e.meridional.casing.z_le = -30.0;
    e.meridional.hub.z_te = e.meridional.casing.z_te = 30.0;
    const auto r2 = ev.evaluate(e);
    REQUIRE(r2.ok());
    CHECK(r2.layers[4].chord_mprime == Approx(2.0 * rc.layers[4].chord_mprime).epsilon(1e-12));  // unrolled chord
    CHECK(r2.layers[4].beta_le == Approx(rc.layers[4].beta_le).margin(1e-9));  // angles do not scale
    const auto n = BladeDesign::naca65_defaults();
    for (const auto& l : n.layers) CHECK(l.method == Method::Naca65);
    CHECK(ev.evaluate(n).ok());
}

