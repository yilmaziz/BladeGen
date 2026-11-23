// Tests of the stacking of the sections (design/BladeDesign.hpp, "Stacking"; Vos et al. 2025, Sec. 4):
// the prescribed edges, and a reference point (centre of gravity, LE or TE) stacked with variant 1
// (same Z and theta: a radial stacking line) or variant 2 (closest point to the previous layer's
// reference point: along the stream-surface normals).  Each check is against the defining equation:
// the reference point on its target, the closest-point condition (C(m') - Q) . dC/dm' = 0 and a
// brute-force distance scan, and the cylinder, where both variants must give the same blade.
#include "design/BladeDesign.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen;
using namespace BladeGen::Design;
using Catch::Approx;
using Math::Vec2;
using Math::Vec3;

namespace {

BladeDesign stacked(Stacking::Variant v, Stacking::Reference ref, double z_ref, double theta_ref_deg)
{
    auto d = BladeDesign::defaults();
    d.stacking = {Stacking::Mode::ReferencePoint, v, ref, z_ref, theta_ref_deg};
    return d;
}

void cylindrical(BladeDesign& d)
{
    d.meridional.hub = {-15.0, 50.0, 15.0, 50.0};
    d.meridional.casing = {-15.0, 80.0, 15.0, 80.0};
}

/// The mesh rows at the hub and tip are the placed first and last layers.
void check_end_rows(const DesignResult& r)
{
    const auto& L0 = r.layers.front();
    const auto& L1 = r.layers.back();
    for (int i = 0; i < r.blade.ni; i += 7) {
        CHECK((r.blade.at(i, 0) - L0.xyz[static_cast<std::size_t>(i)]).norm() < 1e-9);
        CHECK((r.blade.at(i, r.blade.nj - 1) - L1.xyz[static_cast<std::size_t>(i)]).norm() < 1e-9);
    }
    CHECK((r.stacking_xr.front() - L0.stack_xr).norm() < 1e-9);
    CHECK((r.stacking_xr.back() - L1.stack_xr).norm() < 1e-9);
}

}  // namespace

TEST_CASE("StreamSurface::dxr: the derivative of xr, of length r", "[design][stacking]")
{
    const Blade::FlowPath fp(BladeDesign::defaults().meridional);
    for (const double span : {0.0, 0.4, 1.0}) {
        const auto S = fp.surface(span);
        REQUIRE_FALSE(S.cylindrical());
        const double a = S.mprime_inlet(), b = S.mprime_outlet();
        for (int i = 1; i < 40; ++i) {
            const double mp = a + (b - a) * (i + 0.37) / 41.0, h = 1e-6;
            const Vec2 fd = (S.xr(mp + h) - S.xr(mp - h)) / (2.0 * h);
            const Vec2 d = S.dxr(mp);
            CHECK((d - fd).norm() < 1e-6 * d.norm());
            CHECK(d.norm() == Approx(S.xr(mp).y()).epsilon(1e-5));  // d(x, r)/dm' = r t, |t| = 1
        }
        // beyond the ends: the end slopes (xr continues straight)
        CHECK((S.dxr(a - 1.0) - S.dxr(a)).norm() < 1e-12);
        CHECK((S.dxr(b + 1.0) - S.dxr(b)).norm() < 1e-12);
    }
    auto d = BladeDesign::defaults();
    cylindrical(d);
    const auto C = Blade::FlowPath(d.meridional).surface(0.5);
    REQUIRE(C.cylindrical());
    CHECK((C.dxr(0.1) - Vec2(65.0, 0.0)).norm() < 1e-12);
}

TEST_CASE("closest_mprime: the foot point against a brute-force scan; an end of the channel is an error", "[design][stacking]")
{
    const Blade::FlowPath fp(BladeDesign::defaults().meridional);
    const auto H = fp.surface(0.0);
    for (const double span : {0.25, 0.7, 1.0})
        for (const double x : {-8.0, 0.0, 6.0}) {
            const auto S = fp.surface(span);
            const Vec2 Q = H.xr(mprime_at_x(H, x));
            const double mp = closest_mprime(S, Q);
            const Vec2 e = S.xr(mp) - Q, t = S.dxr(mp);
            CHECK(std::abs(e.dot(t)) < 1e-10 * e.norm() * t.norm());  // Q lies on the normal through the foot point
            // no point of the streamline is closer
            const double a = S.mprime_inlet(), b = S.mprime_outlet();
            double best = 1e300, at = 0.0;
            const int n = 200000;
            for (int i = 0; i <= n; ++i) {
                const double m = a + (b - a) * i / n;
                if (const double v = (S.xr(m) - Q).norm(); v < best) {
                    best = v;
                    at = m;
                }
            }
            CHECK(e.norm() <= best + 1e-12);
            CHECK(std::abs(mp - at) < 2.0 * (b - a) / n);
        }
    CHECK_THROWS_AS(closest_mprime(fp.surface(0.5), Vec2(-1000.0, 70.0)), std::invalid_argument);
    CHECK_THROWS_AS(closest_mprime(fp.surface(0.5), Vec2(1000.0, 70.0)), std::invalid_argument);
}

TEST_CASE("Stacking on the prescribed edges (the default): the placement of v1.0.3", "[design][stacking]")
{
    const auto d = BladeDesign::defaults();
    REQUIRE(d.stacking.mode == Stacking::Mode::Edges);
    const Blade::FlowPath fp(d.meridional);
    const auto st = resolve_stacking(d, fp);
    CHECK(st.mode == Stacking::Mode::Edges);
    CHECK(st.x.empty());
    const auto S = fp.surface(0.3);
    for (const bool mirror : {false, true}) {
        // chord scale and loop are not used on the edges
        const auto a = nominal_placement(st, S, 33.0, -1.5, 1.7, {}, mirror);
        const auto b = placement(33.0, -1.5, S.mprime_te(), mirror);
        CHECK(a.gamma == b.gamma);
        CHECK(a.c == b.c);
        CHECK(a.m0 == b.m0);
        CHECK(a.t0 == b.t0);
        CHECK(a.sign == b.sign);
    }
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    for (const auto& L : r.layers) {  // the stacking point is the LE, on the LE line
        CHECK((L.stack_point - Vec2(L.mprime_le, L.theta_le)).norm() < 1e-15);
        CHECK(L.mprime_le == 0.0);
        CHECK(L.chord_mprime == Approx(L.surface.mprime_te() / std::cos(L.gamma)).epsilon(1e-14));
    }
    REQUIRE(r.stacking_xr.size() == r.row_span.size());
    for (std::size_t j = 0; j < r.row_span.size(); ++j) CHECK((r.stacking_xr[j] - r.blade_le[j]).norm() < 1e-14);
}

TEST_CASE("Variant 1, centre of gravity: every section's centre of gravity at (Z_ref, theta_ref), a radial line", "[design][stacking]")
{
    const double z_ref = 0.5, th_ref = 3.0;
    const auto d = stacked(Stacking::Variant::SameZTheta, Stacking::Reference::CG, z_ref, th_ref);
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    const Blade::FlowPath fp(d.meridional);
    double shift = 0.0;  // how far the chord-frame centroid alone would miss the target
    for (std::size_t k = 0; k < r.layers.size(); ++k) {
        const auto& L = r.layers[k];
        REQUIRE_FALSE(L.surface.cylindrical());
        // the centre of gravity, computed independently from the placed loop, is on the target
        CHECK(L.cg_xr.x() == Approx(z_ref).margin(1e-10));
        CHECK(L.cg_theta == Approx(Math::deg2rad(th_ref)).margin(1e-12));
        CHECK((L.stack_point - Vec2(L.cg_mprime, L.cg_theta)).norm() < 1e-15);
        CHECK(std::atan2(L.cg_xyz.y(), L.cg_xyz.z()) == Approx(Math::deg2rad(th_ref)).margin(1e-12));
        // the chord of the prescribed edges (chord scale 1) and the layer's stagger
        CHECK(L.chord_mprime == Approx(L.surface.mprime_te() / std::cos(L.gamma)).epsilon(1e-14));
        CHECK(Math::rad2deg(L.gamma) == Approx(d.layers[k].stagger_deg).margin(1e-12));
        // the conformal map does not carry the centroid over: placing the chord-frame centroid on the target misses
        const Vec2 qc = polygon_centroid(L.section->loop.xy).centroid;
        auto p = placement(d.layers[k].stagger_deg, 0.0, L.surface.mprime_te(), d.mirror);
        const Vec2 naive = Vec2(mprime_at_x(L.surface, z_ref), Math::deg2rad(th_ref)) - p(qc);
        shift = std::max(shift, std::abs(naive.x() - L.mprime_le));
    }
    INFO("largest m' correction of the centre-of-gravity iteration " << shift);
    CHECK(shift > 1e-5);
    // every mesh row's centre of gravity is on the radial line too
    REQUIRE(r.stacking_xr.size() == r.row_span.size());
    for (std::size_t j = 0; j < r.row_span.size(); ++j) {
        CHECK(r.stacking_xr[j].x() == Approx(z_ref).margin(1e-10));
        CHECK(std::atan2(r.stacking_xyz[j].y(), r.stacking_xyz[j].z()) == Approx(Math::deg2rad(th_ref)).margin(1e-12));
    }
    check_end_rows(r);
}

TEST_CASE("Variant 1 with the LE and with the TE as reference point", "[design][stacking]")
{
    Evaluator ev;
    {
        const auto r = ev.evaluate(stacked(Stacking::Variant::SameZTheta, Stacking::Reference::LE, -12.0, -4.0));
        REQUIRE(r.ok());
        for (const auto& L : r.layers) {
            CHECK(L.camber_xr.front().x() == Approx(-12.0).margin(1e-10));
            CHECK(L.theta_le == Approx(Math::deg2rad(-4.0)).margin(1e-14));
            CHECK((L.stack_point - Vec2(L.mprime_le, L.theta_le)).norm() < 1e-15);
        }
        // the casing's LE is upstream of its LE line (a flaring channel: hub LE Z -15, casing LE Z -11)
        CHECK(r.layers.back().mprime_le < 0.0);
        for (const auto& p : r.stacking_xr) CHECK(p.x() == Approx(-12.0).margin(1e-10));
        check_end_rows(r);
    }
    {
        const auto r = ev.evaluate(stacked(Stacking::Variant::SameZTheta, Stacking::Reference::TE, 12.0, 6.0));
        REQUIRE(r.ok());
        for (const auto& L : r.layers) {
            CHECK(L.camber_xr.back().x() == Approx(12.0).margin(1e-9));
            CHECK(L.theta_deg.back() == Approx(6.0).margin(1e-9));
        }
        for (const auto& p : r.stacking_xr) CHECK(p.x() == Approx(12.0).margin(1e-10));
        check_end_rows(r);
    }
}

TEST_CASE("Variant 2, centre of gravity: the closest-point chain, along the stream-surface normals", "[design][stacking]")
{
    const double th_ref = -2.0;
    const auto d = stacked(Stacking::Variant::ClosestPoint, Stacking::Reference::CG, 0.0, th_ref);
    Evaluator ev;
    const auto r = ev.evaluate(d);
    REQUIRE(r.ok());
    CHECK(r.layers.front().cg_xr.x() == Approx(0.0).margin(1e-10));
    for (std::size_t k = 0; k < r.layers.size(); ++k) {
        const auto& L = r.layers[k];
        CHECK(L.cg_theta == Approx(Math::deg2rad(th_ref)).margin(1e-12));
        CHECK((L.stack_point - Vec2(L.cg_mprime, L.cg_theta)).norm() < 1e-15);
        if (k == 0) continue;
        // (C_k(m') - Q_{k-1}) . dC_k/dm' = 0: the previous centre of gravity is on this layer's normal
        const Vec2 Q = r.layers[k - 1].cg_xr;
        const Vec2 e = L.cg_xr - Q, t = L.surface.dxr(L.cg_mprime);
        CHECK(std::abs(e.dot(t)) < 1e-9 * e.norm() * t.norm());
        // and no point of the layer's streamline is closer to it
        const auto& S = L.surface;
        double best = 1e300;
        for (int i = 0; i <= 50000; ++i) best = std::min(best, (S.xr(S.mprime_inlet() + (S.mprime_outlet() - S.mprime_inlet()) * i / 50000.0) - Q).norm());
        CHECK(e.norm() <= best + 1e-12);
    }
    // in the flaring channel the stacking line leans: variant 1 would keep Z
    INFO("tip centre of gravity Z " << r.layers.back().cg_xr.x());
    CHECK(std::abs(r.layers.back().cg_xr.x()) > 0.1);
    // the rows take Z of their centre of gravity from the spline through the layers
    const auto st = resolve_stacking(d, Blade::FlowPath(d.meridional));
    for (std::size_t j = 0; j < r.row_span.size(); ++j) CHECK(r.stacking_xr[j].x() == Approx(st.x_at(r.row_span[j])).margin(1e-10));
    for (std::size_t k = 0; k < r.layers.size(); ++k) CHECK(st.x[k] == Approx(r.layers[k].cg_xr.x()).margin(1e-10));
    check_end_rows(r);
}

TEST_CASE("Cylindrical walls: variants 1 and 2 give the same blade; the centre of gravity is the chord-frame centroid", "[design][stacking]")
{
    auto d1 = stacked(Stacking::Variant::SameZTheta, Stacking::Reference::CG, 1.0, 2.0);
    cylindrical(d1);
    auto d2 = d1;
    d2.stacking.variant = Stacking::Variant::ClosestPoint;
    Evaluator ev;
    const auto a = ev.evaluate(d1), b = ev.evaluate(d2);
    REQUIRE(a.ok());
    REQUIRE(b.ok());
    REQUIRE(a.blade.p.size() == b.blade.p.size());
    double gap = 0.0;
    for (std::size_t i = 0; i < a.blade.p.size(); ++i) gap = std::max(gap, (a.blade.p[i] - b.blade.p[i]).norm());
    CHECK(gap < 1e-9);
    for (std::size_t k = 0; k < a.layers.size(); ++k) {
        const auto& L = a.layers[k];
        REQUIRE(L.surface.cylindrical());
        // the placement that puts the chord-frame centroid on the target is the converged one
        auto p = placement(d1.layers[k].stagger_deg, 0.0, L.surface.mprime_te(), d1.mirror);
        const Vec2 qc = polygon_centroid(L.section->loop.xy).centroid;
        const Vec2 sh = Vec2(mprime_at_x(L.surface, 1.0), Math::deg2rad(2.0)) - p(qc);
        p.m0 += sh.x();
        p.t0 += sh.y();
        for (std::size_t i = 0; i < L.conformal.size(); i += 5) CHECK((L.conformal[i] - p(L.section->loop.xy[i])).norm() < 1e-12);
    }
}

TEST_CASE("Chord scale: c' scales with it, the reference point stays; not used on the edges", "[design][stacking]")
{
    auto d = stacked(Stacking::Variant::SameZTheta, Stacking::Reference::CG, 0.0, 0.0);
    Evaluator ev;
    const auto a = ev.evaluate(d);
    d.layers[2].chord_scale = 1.3;
    const auto b = ev.evaluate(d);
    REQUIRE(a.ok());
    REQUIRE(b.ok());
    const auto& M = b.layers[2];
    CHECK(M.chord_mprime == Approx(1.3 * a.layers[2].chord_mprime).epsilon(1e-14));
    CHECK(M.cg_xr.x() == Approx(0.0).margin(1e-10));
    CHECK(M.cg_theta == Approx(0.0).margin(1e-12));
    CHECK(M.section_area > 1.5 * a.layers[2].section_area);  // ~1.3^2
    for (const std::size_t k : {0u, 4u}) CHECK((b.layers[k].xyz.front() - a.layers[k].xyz.front()).norm() < 1e-12);
    // on the prescribed edges the chord follows the LE / TE lines
    d.stacking.mode = Stacking::Mode::Edges;
    const auto e = ev.evaluate(d);
    REQUIRE(e.ok());
    CHECK(e.layers[2].chord_mprime == Approx(e.layers[2].surface.mprime_te() / std::cos(e.layers[2].gamma)).epsilon(1e-14));
}

TEST_CASE("Re-stagger of a reference-point stacked blade: rigid about the axis through its hub section", "[design][stacking]")
{
    auto d = stacked(Stacking::Variant::ClosestPoint, Stacking::Reference::CG, 0.0, 1.0);
    Evaluator ev;
    const auto a = ev.evaluate(d);
    REQUIRE(a.ok());
    d.restagger = {5.0, Restagger::Axis::Hub, 0.3, 0.0, 0.0};
    const auto b = ev.evaluate(d);
    REQUIRE(b.ok());
    const Blade::FlowPath fp(d.meridional);
    const Vec2 axis = restagger_axis(d, fp);
    // the axis passes through the stacked hub section's chord point 0.3
    const auto& H = a.layers[0];
    const Vec2 P0(H.mprime_le + 0.3 * H.chord_mprime * std::cos(H.gamma), H.theta_le + 0.3 * H.chord_mprime * std::sin(H.gamma));
    CHECK(mprime_at_x(H.surface, axis.x()) == Approx(P0.x()).margin(1e-12));
    CHECK(axis.y() == Approx(P0.y()).margin(1e-14));
    const double dd = Math::deg2rad(5.0);
    for (std::size_t k = 0; k < a.layers.size(); ++k) {
        const Vec2 P(mprime_at_x(a.layers[k].surface, axis.x()), axis.y());
        for (std::size_t i = 0; i < a.layers[k].conformal.size(); i += 4) {
            const Vec2 v = a.layers[k].conformal[i] - P;
            const Vec2 turned = P + Vec2(std::cos(dd) * v.x() - std::sin(dd) * v.y(), std::sin(dd) * v.x() + std::cos(dd) * v.y());
            CHECK((b.layers[k].conformal[i] - turned).norm() < 1e-11);
        }
        CHECK(b.layers[k].beta_le == Approx(a.layers[k].beta_le + 5.0).margin(1e-11));
    }
    check_end_rows(b);
}

TEST_CASE("Stacking: invalid input, design errors, JSON round trip", "[design][stacking]")
{
    auto d = stacked(Stacking::Variant::SameZTheta, Stacking::Reference::CG, 0.0, 0.0);
    d.layers[1].chord_scale = 0.1;
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    d.layers[1].chord_scale = std::nan("");
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    d.layers[1].chord_scale = 1.0;
    d.stacking.theta_ref_deg = 200.0;
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    d.stacking.theta_ref_deg = 0.0;
    d.stacking.z_ref = std::nan("");
    CHECK_THROWS_AS(d.validate(), std::invalid_argument);
    // a reference point outside the channel is a design error (reported, not thrown)
    Evaluator ev;
    for (const auto v : {Stacking::Variant::SameZTheta, Stacking::Variant::ClosestPoint}) {
        d.stacking.variant = v;
        d.stacking.z_ref = 1e4;
        const auto r = ev.evaluate(d);
        CHECK_FALSE(r.ok());
        REQUIRE_FALSE(r.errors.empty());
        CHECK(r.errors.front().find("stacking") != std::string::npos);
    }
    // the centre of gravity needs a profile
    const Blade::FlowPath fp(d.meridional);
    d.stacking.z_ref = 0.0;
    const auto st = resolve_stacking(d, fp);
    CHECK_THROWS_AS(nominal_placement(st, fp.surface(0.0), 30.0, 0.0, 1.0, {}, false), std::invalid_argument);
    CHECK_THROWS_AS(nominal_placement(st, fp.surface(0.0), 30.0, 0.0, 0.0, {}, false), std::invalid_argument);
    // JSON
    d.stacking = {Stacking::Mode::ReferencePoint, Stacking::Variant::ClosestPoint, Stacking::Reference::TE, 2.5, -7.0};
    d.layers[3].chord_scale = 1.15;
    const auto e = design_from_json(to_json(d));
    CHECK(e.stacking.mode == Stacking::Mode::ReferencePoint);
    CHECK(e.stacking.variant == Stacking::Variant::ClosestPoint);
    CHECK(e.stacking.reference == Stacking::Reference::TE);
    CHECK(e.stacking.z_ref == 2.5);
    CHECK(e.stacking.theta_ref_deg == -7.0);
    CHECK(e.layers[3].chord_scale == 1.15);
    auto j = to_json(d);
    j.erase("stacking");  // designs saved before v1.0.4: the prescribed edges
    CHECK(design_from_json(j).stacking.mode == Stacking::Mode::Edges);
    for (const char* key : {"mode", "variant", "reference"}) {
        j = to_json(d);
        j["stacking"][key] = "radial";
        CHECK_THROWS_AS(design_from_json(j), std::invalid_argument);
    }
    for (auto m : {Stacking::Mode::Edges, Stacking::Mode::ReferencePoint}) CHECK(stacking_mode_from_string(to_string(m)) == m);
    for (auto v : {Stacking::Variant::SameZTheta, Stacking::Variant::ClosestPoint}) CHECK(stacking_variant_from_string(to_string(v)) == v);
    for (auto r : {Stacking::Reference::CG, Stacking::Reference::LE, Stacking::Reference::TE})
        CHECK(stacking_reference_from_string(to_string(r)) == r);
}
