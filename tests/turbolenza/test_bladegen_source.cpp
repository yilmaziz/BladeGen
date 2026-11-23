// BladeGenSource.hpp: a BladeGen design layer -> TURBOlenza's BladeSection (the quintic fit of
// BladeGen's own side curves, the blade angles, the edge metadata).
#include "RefData.hpp"

#include "turbolenza/BladeGenSource.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>

namespace T = BladeGen::Turbolenza;
namespace D = BladeGen::Design;
namespace P = BladeGen::Profile;

namespace {

struct Deviation
{
    double distance = 0.0;      ///< max |normal distance| of the check points from the fitted curve
    double at_x = 0.0;          ///< x/c of the worst point
    double curvature = 0.0;     ///< max |kappa_fit - kappa| / max(|kappa|, 1)
    double curvature_at_x = 0.0;
    double smooth = 0.0;        ///< the same where BladeGen's curvature varies by < 1 % over the local window
    double smooth_at_x = 0.0;
    double overshoot = 0.0;     ///< fit curvature outside BladeGen's range over the window, / max(|kappa|, 1)
    double overshoot_at_x = 0.0;
    int inflections_fit = 0, inflections_ref = 0;
};

int sign_changes(const std::vector<double>& k)
{
    int n = 0;
    double last = 0.0;
    for (double v : k) {
        if (std::abs(v) < 1e-9) continue;  // a zero passes no inflection by itself
        if (last != 0.0 && (v > 0.0) != (last > 0.0)) ++n;
        last = v;
    }
    return n;
}

/// The check points (a sample independent of the fit's) against the fitted curve: each point
/// projected onto the curve (Newton from its arc-length fraction).  The window of a check point:
/// its knot span of the fit and one span width on either side.
Deviation deviation(const T::ParametricBladeCurve& c, const T::SideSample& s)
{
    Deviation d;
    const std::size_t n = s.xy.size();
    std::vector<double> s_len(n, 0.0), u(n), kfit(n);
    for (std::size_t k = 1; k < n; ++k)
        s_len[k] = s_len[k - 1] + std::hypot(s.xy[k].m - s.xy[k - 1].m, s.xy[k].y - s.xy[k - 1].y);
    for (std::size_t k = 0; k < n; ++k) {
        const auto pr = c.project(s.xy[k], s_len[k] / s_len.back());
        u[k] = pr[0];
        kfit[k] = c.curvature(pr[0]);
        if (std::abs(pr[1]) > d.distance) {
            d.distance = std::abs(pr[1]);
            d.at_x = s.xy[k].m;
        }
    }
    std::vector<double> knots = {0.0};
    const auto ik = c.interior_knots();
    knots.insert(knots.end(), ik.begin(), ik.end());
    knots.push_back(1.0);
    for (std::size_t k = 0; k < n; ++k) {
        const auto it = std::upper_bound(knots.begin(), knots.end(), u[k]);
        const auto pos = std::max<std::ptrdiff_t>(it - knots.begin() - 1, 0);
        const std::size_t a = std::min(knots.size() - 2, static_cast<std::size_t>(pos));
        const double w = knots[a + 1] - knots[a];
        double lo = s.kappa[k], hi = s.kappa[k];
        // u of the check points is monotone (projections of an ordered sample)
        for (std::size_t q = k; q-- > 0 && u[q] >= knots[a] - w;) {
            lo = std::min(lo, s.kappa[q]);
            hi = std::max(hi, s.kappa[q]);
        }
        for (std::size_t q = k + 1; q < n && u[q] <= knots[a + 1] + w; ++q) {
            lo = std::min(lo, s.kappa[q]);
            hi = std::max(hi, s.kappa[q]);
        }
        const double scale = std::max(std::abs(s.kappa[k]), 1.0);
        const double e = std::abs(kfit[k] - s.kappa[k]) / scale;
        if (e > d.curvature) {
            d.curvature = e;
            d.curvature_at_x = s.xy[k].m;
        }
        if ((hi - lo) / scale < 0.01 && e > d.smooth) {
            d.smooth = e;
            d.smooth_at_x = s.xy[k].m;
        }
        const double o = std::max({0.0, lo - kfit[k], kfit[k] - hi}) / scale;
        if (o > d.overshoot) {
            d.overshoot = o;
            d.overshoot_at_x = s.xy[k].m;
        }
    }
    d.inflections_fit = sign_changes(kfit);
    d.inflections_ref = sign_changes(s.kappa);
    return d;
}

D::Layer naca_layer(double cl0, double tc, BladeGen::Naca::Variant v)
{
    D::Layer l;
    l.method = D::Method::Naca65;
    l.naca.cl0 = cl0;
    l.naca.t_max_over_c = tc;
    l.naca.variant = v;
    return l;
}

D::Layer custom_layer(const BladeGen::Profile::ProfileParams& p)
{
    D::Layer l;
    l.method = D::Method::Custom;
    l.custom = p;
    return l;
}

/// Curvature bounds of one layer (measured first, see the printed lines; bound ~2.5x measured).
struct CurvatureBounds
{
    double smooth = 5e-2;       ///< curvature error where BladeGen's curvature is smooth over the window
    double overshoot = 5e-2;    ///< the wiggle measure
    int extra_inflections = 0;
};

void check_layer(const D::Layer& layer, const char* name, const CurvatureBounds& bound = {})
{
    D::Evaluator ev;
    T::SectionFitReport rep;
    const T::BladeSection b = T::section_from_layer(layer, ev, &rep);
    // the end points are exact (fixed control points; the clamped ends are those points)
    for (const auto* c : {&b.suction, &b.pressure}) {
        CHECK(c->degree() == 5);
        CHECK(c->evaluate(0.0).m == 0.0);
        CHECK(c->evaluate(0.0).y == 0.0);
        CHECK(c->evaluate(1.0).m == 1.0);
        CHECK(c->evaluate(1.0).y == 0.0);
    }
    CHECK(b.suction.reference_length() == Catch::Approx(rep.suction.arc_length).epsilon(1e-15));
    CHECK(b.suction_parameters.size() == 201);
    CHECK(b.suction_parameters.front() == 0.0);
    CHECK(b.suction_parameters.back() == 1.0);
    // the suction side is the +y side
    CHECK(b.suction.evaluate(0.4).y > b.pressure.evaluate(0.4).y);
    // an independent check sample (other counts, so other points) of BladeGen's own curve
    const auto ss = T::sample_layer_side(layer, ev, P::Side::Suction, 777, 3001);
    const auto ps = T::sample_layer_side(layer, ev, P::Side::Pressure, 777, 3001);
    for (int side = 0; side < 2; ++side) {
        const auto& r = side == 0 ? rep.suction : rep.pressure;
        const Deviation d = deviation(side == 0 ? b.suction : b.pressure, side == 0 ? ss : ps);
        std::printf("[bladegen_source] %s %s: %d ctrl (%d rounds, %d samples), fit residual %.2e, fit curvature error %.2e; "
                    "check: distance %.2e c at x/c %.5f; curvature error %.2e at x/c %.5f, where smooth %.2e at x/c %.5f, "
                    "overshoot %.2e at x/c %.5f; inflections %d / %d\n",
                    name, side == 0 ? "SS" : "PS", r.control_points, r.refinements, r.samples, r.max_residual,
                    r.max_relative_curvature_error, d.distance, d.at_x, d.curvature, d.curvature_at_x, d.smooth, d.smooth_at_x,
                    d.overshoot, d.overshoot_at_x, d.inflections_fit, d.inflections_ref);
        // the requirement: within 1e-7 chord of BladeGen's curve
        CHECK(d.distance < 1e-7);
        // the curvature follows BladeGen's.  BladeGen's pieces join G2, not G3: at the LE / TE
        // fairs, the Kovach blend and the custom connectors dkappa/ds jumps (by ~1e4 .. 1e5 / c^2),
        // which a C4 spline rounds over its narrowest knot span, with a ringing of a few per cent
        // of max(|kappa|, 1) next to the join (measured: the printed lines).  Away from the joins
        // the error is below the fit's target 1e-3; the inflections are BladeGen's own.
        CHECK(d.smooth < bound.smooth);
        CHECK(d.overshoot < bound.overshoot);
        CHECK(d.inflections_fit >= d.inflections_ref);
        CHECK(d.inflections_fit <= d.inflections_ref + bound.extra_inflections);
    }
}

}  // namespace

TEST_CASE("the quintic fit reproduces a circular arc and its curvature", "[turbolenza][bladegen]")
{
    // closed form: the arc of radius R through (0, 0) and (1, 0), bulging toward +y, LE -> TE
    // (clockwise: signed curvature -1/R)
    const double R = 0.8;
    const double yc = -std::sqrt(R * R - 0.25);
    const double a0 = std::atan2(-yc, -0.5), a1 = std::atan2(-yc, 0.5);  // angles of the end points from the centre
    T::SideSample s;
    const int n = 4000;
    for (int k = 0; k < n; ++k) {
        const double a = a0 + (a1 - a0) * k / (n - 1);
        s.xy.push_back({0.5 + R * std::cos(a), yc + R * std::sin(a)});
        s.kappa.push_back(-1.0 / R);
    }
    s.xy.front() = {0.0, 0.0};
    s.xy.back() = {1.0, 0.0};
    T::SideFitReport rep;
    const auto c = T::fit_side(s, {}, &rep);
    CHECK(rep.max_residual < 1e-8);
    // the chord-length sum of 4000 points is short of the arc by n kappa^2 h^3 / 24 ~ 5e-9 (relative)
    CHECK(rep.arc_length == Catch::Approx(R * (a0 - a1)).epsilon(1e-8));
    for (double u : {0.0, 0.013, 0.25, 0.5, 0.77, 0.999, 1.0}) {
        const auto x = c.evaluate(u);
        CHECK(std::abs(std::hypot(x.m - 0.5, x.y - yc) - R) < 1e-8);
        CHECK(c.curvature(u) == Catch::Approx(-1.0 / R).epsilon(1e-4));
        // u is the arc-length fraction: the point lies at the angle a0 + u (a1 - a0)
        const double a = std::atan2(x.y - yc, x.m - 0.5);
        CHECK(std::abs(a - (a0 + u * (a1 - a0))) < 1e-7);
    }
}

TEST_CASE("fit_side and surface_parameters refuse degenerate input", "[turbolenza][bladegen]")
{
    T::SideSample s;
    for (int k = 0; k < 50; ++k) {
        s.xy.push_back({k / 49.0, 0.0});
        s.kappa.push_back(0.0);
    }
    CHECK_THROWS_AS(T::fit_side(s), std::invalid_argument);  // too few samples
    s = {};
    for (int k = 0; k < 500; ++k) {
        const double x = k / 499.0;
        s.xy.push_back({x, 0.1 * x * (1.0 - x)});
        s.kappa.push_back(0.0);
    }
    s.xy.back() = {1.0, 1e-6};  // does not end at the TE point
    CHECK_THROWS_AS(T::fit_side(s), std::invalid_argument);
    s.xy.back() = {1.0, 0.0};
    s.kappa.pop_back();  // curvature missing
    CHECK_THROWS_AS(T::fit_side(s), std::invalid_argument);
    s.kappa.push_back(0.0);
    T::SideFitOptions bad;
    bad.tolerance = 0.0;
    CHECK_THROWS_AS(T::fit_side(s, bad), std::invalid_argument);
    // a straight chord is no blade side (it is not longer than the chord)
    T::SideSample flat;
    for (int k = 0; k < 500; ++k) {
        flat.xy.push_back({k / 499.0, 0.0});
        flat.kappa.push_back(0.0);
    }
    CHECK_THROWS_AS(T::fit_side(flat), std::invalid_argument);
    CHECK_THROWS_AS(T::surface_parameters(1), std::invalid_argument);
    const auto u = T::surface_parameters(5);
    REQUIRE(u.size() == 5);
    CHECK(u[2] == Catch::Approx(0.5).margin(1e-15));
    CHECK(u[1] == Catch::Approx(0.5 * (1.0 - std::cos(std::numbers::pi / 4))).margin(1e-15));
}

TEST_CASE("a NACA 65-(12A10)10 layer: fit within 1e-7 c, curvature followed, equivalent-arc angles", "[turbolenza][bladegen]")
{
    const auto layer = naca_layer(1.2, 0.10, BladeGen::Naca::Variant::KovachSandercock);
    check_layer(layer, "NACA 65-(12A10)10 kovach_sandercock");
    D::Evaluator ev;
    T::SectionFitReport rep;
    const auto b = T::section_from_layer(layer, ev, &rep);
    // Aungier Eq. 4-6: tan(camber / 4) = 0.1103 cl0, chi_le = chi_te = camber / 2
    const double camber = 4.0 * std::atan(0.1103 * 1.2);
    CHECK(b.chi_le == Catch::Approx(camber / 2.0).epsilon(1e-15));
    CHECK(b.chi_te == Catch::Approx(camber / 2.0).epsilon(1e-15));
    CHECK(b.designation == "NACA 65-(12A10)10 (kovach_sandercock)");
    // nose circle R_LE = 0.00687 (t / 0.10)^2 through (0, 0); Kovach TE circle 0.008 c through (1, 0)
    CHECK(b.le_radius == Catch::Approx(0.00687).epsilon(1e-12));
    CHECK(std::hypot(b.le_center.m, b.le_center.y) == Catch::Approx(b.le_radius).epsilon(1e-9));
    CHECK(b.te_radius == Catch::Approx(0.008).epsilon(1e-12));
    CHECK(std::hypot(b.te_center.m - 1.0, b.te_center.y) == Catch::Approx(b.te_radius).epsilon(1e-9));
}

TEST_CASE("other NACA 65 layers: thin, highly cambered, the sharp Emery TE", "[turbolenza][bladegen]")
{
    // The thin section's Kovach blend (station 0.60 +- 0.005: a 0.5 deg corner rounded over 0.01 c)
    // and Emery's sharp TE closing piece have the steepest curvature changes: wider bounds
    // (measured: overshoot 0.15 at the blend, 0.10 at the sharp TE with one extra inflection there).
    check_layer(naca_layer(0.4, 0.05, BladeGen::Naca::Variant::KovachSandercock), "NACA 65-(4A10)05 kovach_sandercock",
                {0.15, 0.4, 0});
    check_layer(naca_layer(2.0, 0.12, BladeGen::Naca::Variant::EmeryBase), "NACA 65-(20A10)12 emery_base", {0.15, 0.4, 1});
}

TEST_CASE("custom layers: fit within 1e-7 c, curvature followed, the camber-line angles", "[turbolenza][bladegen]")
{
    const auto p = BladeGen::Profile::ProfileParams::naca65_1210();
    const auto layer = custom_layer(p);
    check_layer(layer, "Custom naca65_1210");
    D::Evaluator ev;
    const auto b = T::section_from_layer(layer, ev);
    CHECK(b.le_radius == p.le.radius);
    CHECK(b.designation == "Custom");
    // chi: BladeGen's own camber-line end angles, the LE / TE metal angles its designer reports
    const auto sec = ev.section(layer, D::MeshOptions{});
    const double le_meas = sec->camber.angle_deg.front(), te_meas = -sec->camber.angle_deg.back();
    CHECK(b.chi_le == Catch::Approx(le_meas * std::numbers::pi / 180.0).epsilon(1e-15));
    CHECK(b.chi_te == Catch::Approx(te_meas * std::numbers::pi / 180.0).epsilon(1e-15));
    std::printf("[bladegen_source] Custom naca65_1210: chi_le %.4f deg (edge axis %.4f, wedge bisector %.4f), "
                "chi_te %.4f deg (%.4f, %.4f)\n",
                le_meas, p.le.beta_deg, p.le.beta_deg + 0.5 * (p.le.d_beta_ss_deg - p.le.d_beta_ps_deg), te_meas, p.te.beta_deg,
                p.te.beta_deg + 0.5 * (p.te.d_beta_ss_deg - p.te.d_beta_ps_deg));
    // this profile is BladeGen's fit of NACA 65-(12A10)10: a cambered section whose camber-line
    // angles lie between the A10 line's equivalent arc (15.1 deg) and its steeper tangents
    CHECK(b.chi_le > 0.15);
    CHECK(b.chi_te > 0.15);
    CHECK(b.chi_le + b.chi_te < 1.0);
    // the same through BladeGen's placed layer: beta_LE - stagger
    D::BladeDesign design = D::BladeDesign::defaults();
    design.layers[2].method = D::Method::Custom;
    design.layers[2].custom = p;
    const auto res = ev.evaluate(design);
    REQUIRE(res.ok());
    CHECK(res.layers[2].beta_le - design.layers[2].stagger_deg == Catch::Approx(le_meas).margin(1e-9));
    CHECK(design.layers[2].stagger_deg - res.layers[2].beta_te == Catch::Approx(te_meas).margin(1e-9));
    // the NASA Rotor 67 mid-span profile (the default design's Mid layer)
    check_layer(custom_layer(D::rotor67_layer(0.5).profile), "Custom Rotor 67 mid");
}

TEST_CASE("BladeGen's NACA 65-(12A10)10 against v1's frozen C2 curves", "[turbolenza][bladegen][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(fam, "grid_family");
    const auto v1 = ref::surface_curves(fam.at("cascade").at("active_curves"));
    D::Evaluator ev;
    const auto b = T::section_from_layer(naca_layer(1.2, 0.10, BladeGen::Naca::Variant::KovachSandercock), ev);
    // the same construction steps (Aungier Table 4-1, A10 mean line, Kovach taper), but
    // BladeGen's own LE / TE fairs: measured, reported, bounded loosely
    for (int side = 0; side < 2; ++side) {
        const auto& ours = side == 0 ? b.suction : b.pressure;
        const auto& theirs = side == 0 ? v1.suction_unstaggered : v1.pressure_unstaggered;
        double worst = 0.0, at = 0.0, body = 0.0, body_at = 0.0;
        for (int k = 0; k <= 4000; ++k) {
            const double u = 0.5 * (1.0 - std::cos(std::numbers::pi * k / 4000.0));
            const auto x = ours.evaluate(u);
            const double d = std::abs(theirs.project(x, u)[1]);
            if (d > worst) {
                worst = d;
                at = x.m;
            }
            if (x.m > 0.02 && x.m < 0.98 && d > body) {
                body = d;
                body_at = x.m;
            }
        }
        std::printf("[bladegen_source] BladeGen vs v1 C2 %s: max distance %.3e c at x/c %.5f; 0.02 < x/c < 0.98: %.3e c at x/c %.5f\n",
                    side == 0 ? "suction" : "pressure", worst, at, body, body_at);
        CHECK(worst < 1e-3);
    }
}
