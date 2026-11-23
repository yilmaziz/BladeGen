// Tests of io/SectionExport (bladegen export-sections): format and version fields, a
// bit-identical round trip of every double, the (m', theta) / (z, r, theta) pairs against
// BladeGen's own transform, and the NACA 65 cascade's chord, stagger and pitch read back
// from the points.
#include "io/BladeJson.hpp"
#include "io/SectionExport.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

using namespace BladeGen;
using Catch::Approx;
using Math::Vec2;
using Math::Vec3;
using nlohmann::ordered_json;

namespace {

std::string src(const std::string& p) { return std::string(BLADEGEN_SOURCE_DIR) + "/" + p; }

const Blade::BladeDefinition& rotor67()
{
    static const auto d = IO::read_blade_json(src("data/rotor67/rotor67_bladegen.json"));
    return d;
}

const ordered_json& rotor67_hub()
{
    static const auto j = IO::export_design_sections(rotor67(), {}, "rotor67_bladegen.json", "test");
    return j;
}

IO::CascadeExportOptions naca_1210()
{
    IO::CascadeExportOptions o;
    o.naca = {1.2, 0.10, Naca::Variant::KovachSandercock};
    o.chord = 0.5;
    o.stagger_deg = 33.0;
    o.solidity = 1.0;
    return o;
}

/// Same structure, same keys in the same order, every double bit-identical.
void require_identical(const ordered_json& a, const ordered_json& b, std::size_t& n_doubles)
{
    if (a.is_number_float()) {
        REQUIRE(b.is_number_float());
        const auto x = std::bit_cast<std::uint64_t>(a.get<double>()), y = std::bit_cast<std::uint64_t>(b.get<double>());
        if (x != y) FAIL("double changed in the round trip: " << a.get<double>() << " -> " << b.get<double>());
        ++n_doubles;
        return;
    }
    if (a.is_number_integer()) {  // signed and unsigned JSON integers compare by value
        REQUIRE(b.is_number_integer());
        REQUIRE(a == b);
        return;
    }
    REQUIRE(a.type() == b.type());
    if (a.is_object()) {
        REQUIRE(a.size() == b.size());
        auto ib = b.begin();
        for (auto ia = a.begin(); ia != a.end(); ++ia, ++ib) {
            REQUIRE(ia.key() == ib.key());
            require_identical(ia.value(), ib.value(), n_doubles);
        }
    } else if (a.is_array()) {
        REQUIRE(a.size() == b.size());
        for (std::size_t k = 0; k < a.size(); ++k) require_identical(a[k], b[k], n_doubles);
    } else {
        REQUIRE(a == b);
    }
}

std::size_t round_trip(const ordered_json& j)
{
    const auto back = ordered_json::parse(IO::dump_sections(j));
    std::size_t n = 0;
    require_identical(j, back, n);
    return n;
}

Vec2 v2(const ordered_json& a) { return {a[0].get<double>(), a[1].get<double>()}; }

}  // namespace

TEST_CASE("export-sections: format, version and the fields of a blade design", "[export]")
{
    const auto j = ordered_json::parse(IO::dump_sections(rotor67_hub()));
    REQUIRE(j["format"] == "bladegen.sections");
    REQUIRE(j["version"].is_number_integer());
    REQUIRE(j["version"] == 1);
    REQUIRE(j["kind"] == "blade_design");
    REQUIRE(j["generator"]["name"] == "BladeGen");
    REQUIRE(j["generator"]["version"] == "test");
    REQUIRE(j["source"]["design_file"] == "rotor67_bladegen.json");
    REQUIRE(j["source"]["design_name"] == "NASA_Rotor_67");
    REQUIRE(j["units"]["length"] == "m");
    REQUIRE(j["units"]["angle"] == "rad");
    REQUIRE(j["n_blades"] == 22);
    REQUIRE(j["pitch_angle"].get<double>() == 2.0 * Math::kPi / 22.0);
    for (const char* w : {"hub", "shroud"}) REQUIRE(j["channel"][w]["zr"].size() == 2001);
    // hub below the shroud everywhere (z, r)
    for (std::size_t i = 0; i < 2001; i += 10) REQUIRE(j["channel"]["hub"]["zr"][i][1] < j["channel"]["shroud"]["zr"][i][1]);

    REQUIRE(j["sections"].size() == 1);
    const auto& s = j["sections"][0];
    REQUIRE(s["span"].get<double>() == 0.0);
    // Profile2D::sample_loop(40, 160): 2 * 40 + 160 - 2 points from the LE to the TE point
    for (const char* side : {"suction_side", "pressure_side"}) {
        const auto& sd = s[side];
        REQUIRE(sd["n_points"] == 238);
        REQUIRE(sd["mprime_theta"].size() == 238);
        REQUIRE(sd["zrtheta"].size() == 238);
        REQUIRE(sd["piece"].size() == 238);
        REQUIRE(sd["mprime_theta"].front() == s["leading_edge"]["mprime_theta"]);
        REQUIRE(sd["mprime_theta"].back() == s["trailing_edge"]["mprime_theta"]);
        REQUIRE(sd["zrtheta"].front() == s["leading_edge"]["zrtheta"]);
        REQUIRE(sd["zrtheta"].back() == s["trailing_edge"]["zrtheta"]);
        // LE connector, body, TE connector in that order
        const int body = std::string(side) == "suction_side" ? 1 : 3;
        REQUIRE(sd["piece"].front() == 0);
        REQUIRE(sd["piece"][100] == body);
        REQUIRE(sd["piece"].back() == 2);
        // m' increases from the LE to the TE over the body (no loop back)
        REQUIRE(sd["mprime_theta"].back()[0] > sd["mprime_theta"].front()[0]);
    }
    // the LE / TE points are the placement's m'_LE and m'_TE
    REQUIRE(s["leading_edge"]["mprime_theta"][0].get<double>() == Approx(s["mprime_le"].get<double>()).margin(1e-12));
    REQUIRE(s["trailing_edge"]["mprime_theta"][0].get<double>() == Approx(s["mprime_te"].get<double>()).margin(1e-12));
    REQUIRE(s["leading_edge"]["mprime_theta"][1].get<double>() == Approx(s["theta_le"].get<double>()).margin(1e-12));
    // the suction side faces +theta (no mirror): at mid-body the SS lies above the PS
    REQUIRE(s["mirror"] == false);
    // stagger from the chord points: tan(gamma) = d theta / d m'
    const Vec2 le = v2(s["leading_edge"]["mprime_theta"]), te = v2(s["trailing_edge"]["mprime_theta"]);
    REQUIRE(std::atan2(te.y() - le.y(), te.x() - le.x()) == Approx(s["stagger"].get<double>()).margin(1e-12));
    REQUIRE((te - le).norm() == Approx(s["chord_mprime"].get<double>()).margin(1e-12));

    // several spans, and the output carries no directory
    IO::DesignExportOptions o;
    o.spans = {0.0, 0.5, 1.0};
    const auto k = IO::export_design_sections(rotor67(), o, "rotor67_bladegen.json", "test");
    REQUIRE(k["sections"].size() == 3);
    REQUIRE(k["sections"][2]["span"].get<double>() == 1.0);
    const std::string text = IO::dump_sections(k);
    REQUIRE(text.find(BLADEGEN_SOURCE_DIR) == std::string::npos);
    REQUIRE(text.find(":\\") == std::string::npos);
}

TEST_CASE("export-sections: every double reads back bit-identical", "[export]")
{
    // awkward values: shortest forms, 17-digit forms, subnormal, extremes, signed zero
    const ordered_json special = ordered_json::array({0.1, 1.0 / 3.0, 2.0 / 3.0, Math::kPi, 1.0, -0.0, 0.0, 1e23,
                                                      5e-324, 2.2250738585072014e-308, 1.7976931348623157e308,
                                                      -123456789.12345679, 9007199254740993.0, 0.30000000000000004});
    REQUIRE(round_trip(special) == special.size());
    // a whole number stays a floating-point JSON number
    const auto one = ordered_json::parse(IO::dump_sections(ordered_json::array({1.0, -2.0, 1e22})));
    for (const auto& v : one) REQUIRE(v.is_number_float());
    // -0.0 keeps its sign
    REQUIRE(std::signbit(ordered_json::parse(IO::dump_sections(ordered_json::array({-0.0})))[0].get<double>()));

    // the whole blade-design file and the whole cascade file, through a file on disk
    REQUIRE(round_trip(rotor67_hub()) > 3000);
    const auto c = IO::export_naca65_cascade(naca_1210(), "test");
    REQUIRE(round_trip(c) > 1000);
    const std::string path = std::string(BLADEGEN_TEST_OUTPUT) + "/naca65_cascade_roundtrip.json";
    IO::write_sections(c, path);
    std::ifstream f(path, std::ios::binary);
    std::stringstream buf;
    buf << f.rdbuf();
    std::size_t n = 0;
    require_identical(c, ordered_json::parse(buf.str()), n);
    REQUIRE(n > 1000);
    // 17 significant digits: 0.1 is written as 0.10000000000000001
    REQUIRE(IO::dump_sections(ordered_json::array({0.1})).find("0.10000000000000001") != std::string::npos);

    // non-finite numbers are refused
    REQUIRE_THROWS_AS(IO::dump_sections(ordered_json::array({std::numeric_limits<double>::quiet_NaN()})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(IO::dump_sections(ordered_json::array({std::numeric_limits<double>::infinity()})),
                      std::invalid_argument);
}

namespace {

/// Every (m', theta) point of a written section against (z, r, theta): BladeGen's own
/// transform m' -> m -> (z, r) on `sl`, the curve the streamline was written from.  The 3-D
/// points are within `max_move` of BladeGen's 3-D section (0: the same points).
void check_pairs(const ordered_json& s, const Blade::MeridionalCurve& sl, const Blade::Section3D& sec, double max_move)
{
    const std::size_t ite = (sec.xyz.size() - 1) / 2;
    for (const char* side : {"suction_side", "pressure_side"}) {
        const bool ss = std::string(side) == "suction_side";
        const auto& mt = s[side]["mprime_theta"];
        const auto& zrt = s[side]["zrtheta"];
        REQUIRE(mt.size() == ite + 1);
        for (std::size_t i = 0; i < mt.size(); ++i) {
            const double mp = mt[i][0], th = mt[i][1];
            const double z = zrt[i][0], r = zrt[i][1], th3 = zrt[i][2];
            REQUIRE(th3 == th);  // the same theta
            const Vec2 p = sl.at_m(sl.m_of_mprime(mp));
            REQUIRE(z == Approx(p.x()).margin(1e-14));
            REQUIRE(r == Approx(p.y()).margin(1e-14));
            REQUIRE(sl.mprime_of_m(sl.m_of_mprime(mp)) == Approx(mp).margin(1e-11));
            const std::size_t k = ss ? i : (i == 0 ? 0 : 2 * ite - i);
            const Vec3 x = Blade::cartesian(z, r, th);
            if (max_move == 0.0) {
                REQUIRE((x - sec.xyz[k]).norm() < 1e-14);
                REQUIRE((Vec2(mp, th) - sec.conformal[k]).norm() < 1e-15);
            } else {
                // moved onto the written curve, to the foot point: the offset is normal to it
                const Vec2 q = sec.streamline.at_m(sec.streamline.m_of_mprime(sec.conformal[k].x()));
                const Vec2 off = q - p;
                REQUIRE(off.norm() <= max_move);
                REQUIRE(std::abs(off.dot(sl.tangent_at_m(sl.m_of_mprime(mp)))) < 1e-12);
                REQUIRE(th == sec.conformal[k].y());
            }
        }
    }
    // LE / TE and the placement values describe the written points
    const Vec2 le = v2(s["leading_edge"]["mprime_theta"]), te = v2(s["trailing_edge"]["mprime_theta"]);
    REQUIRE(le.x() == Approx(s["mprime_le"].get<double>()).margin(1e-12));
    REQUIRE(te.x() == Approx(s["mprime_te"].get<double>()).margin(1e-12));
    REQUIRE(s["m_le_fraction"].get<double>() * sl.length() == Approx(sl.m_of_mprime(le.x())).margin(1e-12));
    REQUIRE(s["m_te_fraction"].get<double>() * sl.length() == Approx(sl.m_of_mprime(te.x())).margin(1e-12));
}

/// The written streamline: inlet to outlet of `sl`, equal steps in m no larger than
/// `max_step`, and m' = int dm / r by the trapezoidal rule over its own table, whole range.
void check_streamline(const ordered_json& st, const Blade::MeridionalCurve& sl, double max_step)
{
    const std::size_t n = st["zr"].size();
    REQUIRE(st["n_points"] == n);
    REQUIRE(st["m"].size() == n);
    REQUIRE(st["mprime"].size() == n);
    REQUIRE(st["m"][0].get<double>() == 0.0);
    REQUIRE(st["mprime"][0].get<double>() == Approx(0.0).margin(1e-15));  // m' = 0 at the inlet
    REQUIRE(st["m"].back().get<double>() == Approx(sl.length()).margin(1e-15));
    double integral = 0.0, worst = 0.0;
    for (std::size_t i = 1; i < n; ++i) {
        const double dm = st["m"][i].get<double>() - st["m"][i - 1].get<double>();
        REQUIRE(dm <= max_step);
        integral += 0.5 * dm * (1.0 / st["zr"][i][1].get<double>() + 1.0 / st["zr"][i - 1][1].get<double>());
        worst = std::max(worst, std::abs(integral - st["mprime"][i].get<double>()));
        const double chord = (v2(st["zr"][i]) - v2(st["zr"][i - 1])).norm();
        REQUIRE(dm == Approx(chord).epsilon(1e-6));
    }
    REQUIRE(worst < 1e-7 * st["mprime"].back().get<double>());
}

}  // namespace

TEST_CASE("export-sections: Rotor 67 hub on the channel's hub curve, inlet to outlet", "[export]")
{
    const auto j = ordered_json::parse(IO::dump_sections(rotor67_hub()));
    const auto& s = j["sections"][0];
    const auto& hub = rotor67().channel.hub();
    const Blade::Blade blade(rotor67());
    const auto sec = blade.section(0.0, 40, 160);

    // the streamline is the hub curve itself (one curve, so to round-off), over its whole
    // length, and no coarser than the former blade-length table (201 points, 0.47 mm)
    const auto& st = s["streamline"];
    REQUIRE(st["zr"].size() == 2001);
    for (std::size_t i = 0; i < st["zr"].size(); ++i) {
        const double m = st["m"][i];
        REQUIRE((v2(st["zr"][i]) - hub.at_m(m)).norm() == 0.0);
        REQUIRE(st["mprime"][i].get<double>() == hub.mprime_of_m(m));
        REQUIRE(v2(st["zr"][i]) == v2(j["channel"]["hub"]["zr"][i]));  // the channel table: the same points
    }
    REQUIRE(st["zr"].front()[0].get<double>() < -0.18);  // inlet, upstream of the blade (LE at z ~ 0)
    REQUIRE(st["zr"].back()[0].get<double>() > 0.35);    // outlet
    check_streamline(st, hub, 0.47e-3);

    // the section points: moved from the published section-1 surface onto the hub curve
    // (TP-2879 tables: the two curves are ~10 um apart); (m', theta) <-> (z, r, theta) exact on it
    check_pairs(s, hub, sec, 1.5e-5);
    // m' from the inlet: the LE lies about 0.18 m downstream of it
    REQUIRE(s["mprime_le"].get<double>() > 1.5);
    REQUIRE(s["m_le_fraction"].get<double>() > 0.3);
}

TEST_CASE("export-sections: a channel design's sections are BladeGen's own, inlet to outlet", "[export]")
{
    auto def = Blade::BladeDefinition::demo_rotor();
    def.channel = Blade::MeridionalChannel::straight(-0.1, 0.25, 0.25, 0.27, 0.40, 0.395);  // with inlet / outlet
    IO::DesignExportOptions o;
    o.spans = {0.0, 0.5, 1.0};
    const auto j = ordered_json::parse(IO::dump_sections(IO::export_design_sections(def, o, "demo", "test")));
    const Blade::Blade blade(def);
    for (std::size_t k = 0; k < 3; ++k) {
        const double span = o.spans[k];
        const auto& s = j["sections"][k];
        const auto sl = blade.streamline_at(span);
        const auto sec = blade.section(span, 40, 160);
        check_pairs(s, sl, sec, 0.0);  // BladeGen's own points, unchanged
        check_streamline(s["streamline"], sl, 0.47e-3);
        // the streamline written is the one BladeGen placed the section on (one curve, to round-off)
        for (std::size_t i = 0; i < s["streamline"]["zr"].size(); i += 50)
            REQUIRE((v2(s["streamline"]["zr"][i]) - sl.at_m(s["streamline"]["m"][i].get<double>())).norm() == 0.0);
        REQUIRE(s["streamline"]["zr"].front()[0].get<double>() == Approx(-0.1).margin(1e-12));
        REQUIRE(s["streamline"]["zr"].back()[0].get<double>() == Approx(0.25).margin(1e-12));
        REQUIRE(s["mprime_le"].get<double>() == sec.mprime_le);
        REQUIRE(s["chord_mprime"].get<double>() == sec.chord_mprime);
    }
    // span 0 of a channel design: BladeGen's channel streamline at span 0 against the hub curve
    const auto sl0 = blade.streamline_at(0.0);
    double worst = 0.0;
    for (std::size_t i = 0; i < j["sections"][0]["streamline"]["zr"].size(); ++i) {
        const double m = j["sections"][0]["streamline"]["m"][i];
        const double mh = m * def.channel.hub().length() / sl0.length();
        worst = std::max(worst, (sl0.at_m(m) - def.channel.hub().at_m(mh)).norm());
    }
    REQUIRE(worst < 1e-9);
}

TEST_CASE("export-sections: NACA 65 cascade chord, stagger and pitch from the points", "[export]")
{
    const auto o = naca_1210();
    const auto j = ordered_json::parse(IO::dump_sections(IO::export_naca65_cascade(o, "test")));
    REQUIRE(j["format"] == "bladegen.sections");
    REQUIRE(j["version"] == 1);
    REQUIRE(j["kind"] == "naca65_cascade");
    REQUIRE(j["source"]["naca"]["designation"] == "NACA 65-(12A10)10");
    REQUIRE(j["source"]["naca"]["variant"] == "kovach_sandercock");

    const auto& ss = j["suction_side"]["xy"];
    const auto& ps = j["pressure_side"]["xy"];
    // 40 LE-fair points, 101 stations, 40 TE-piece points, shared ends: 2 * 40 + 101 - 2
    REQUIRE(ss.size() == 179);
    REQUIRE(ps.size() == 179);
    REQUIRE(ss.front() == ps.front());
    REQUIRE(ss.back() == ps.back());
    REQUIRE(ss.front() == j["leading_edge"]["xy"]);
    REQUIRE(ss.back() == j["trailing_edge"]["xy"]);

    const Vec2 le = v2(ss.front()), te = v2(ss.back());
    REQUIRE(le.norm() < 1e-15);
    const double chord = (te - le).norm(), stagger = std::atan2(te.y() - le.y(), te.x() - le.x());
    REQUIRE(chord == Approx(0.5).margin(1e-15));
    REQUIRE(stagger == Approx(Math::deg2rad(33.0)).margin(1e-14));
    REQUIRE(j["cascade"]["stagger"].get<double>() == Math::deg2rad(33.0));
    // pitch = chord / solidity, from the points' chord; periodic translation along +y
    REQUIRE(j["cascade"]["pitch"].get<double>() == Approx(chord / j["cascade"]["solidity"].get<double>()).margin(1e-15));
    REQUIRE(j["cascade"]["periodic_translation"][0].get<double>() == 0.0);
    REQUIRE(j["cascade"]["periodic_translation"][1].get<double>() == 0.5);

    // rotate back by -stagger and scale by 1 / chord: BladeGen's own unit-chord section
    const Naca::Naca65Section sec(o.naca);
    const double cg = std::cos(stagger), sg = std::sin(stagger);
    for (const auto sd : {Profile::Side::Suction, Profile::Side::Pressure}) {
        const auto ref = sec.sample_side(sd, 40, 101);
        const auto& side = j[sd == Profile::Side::Suction ? "suction_side" : "pressure_side"];
        for (std::size_t i = 0; i < ref.xy.size(); ++i) {
            const Vec2 p = v2(side["xy"][i]) / chord;
            const Vec2 q(p.x() * cg + p.y() * sg, -p.x() * sg + p.y() * cg);
            REQUIRE((q - ref.xy[i]).norm() < 1e-14);
            REQUIRE(v2(side["xi_eta"][i]) == ref.xy[i]);
        }
    }
    // the suction side is the +y side of the staggered chord: positive cross product at mid-chord
    const Vec2 d = (te - le) / chord;
    REQUIRE(Math::cross2(d, v2(ss[90]) - le) > 0.0);
    REQUIRE(Math::cross2(d, v2(ps[90]) - le) < Math::cross2(d, v2(ss[90]) - le));

    // the sharp-TE variant and a negative stagger
    auto e = o;
    e.naca.variant = Naca::Variant::EmeryBase;
    e.stagger_deg = -20.0;
    const auto k = IO::export_naca65_cascade(e, "test");
    REQUIRE(k["section"]["sharp_te"] == true);
    const Vec2 te2 = v2(k["trailing_edge"]["xy"]);
    REQUIRE(std::atan2(te2.y(), te2.x()) == Approx(Math::deg2rad(-20.0)).margin(1e-14));
}

TEST_CASE("export-sections: invalid input throws", "[export]")
{
    const auto& d = rotor67();
    for (double span : {-0.01, 1.01, std::numeric_limits<double>::quiet_NaN()}) {
        IO::DesignExportOptions o;
        o.spans = {span};
        REQUIRE_THROWS_AS(IO::export_design_sections(d, o, "x", "test"), std::invalid_argument);
    }
    IO::DesignExportOptions o;
    o.spans.clear();
    REQUIRE_THROWS_AS(IO::export_design_sections(d, o, "x", "test"), std::invalid_argument);
    o = {};
    o.n_edge = 1;
    REQUIRE_THROWS_AS(IO::export_design_sections(d, o, "x", "test"), std::invalid_argument);
    o = {};
    o.n_channel = 1;
    REQUIRE_THROWS_AS(IO::export_design_sections(d, o, "x", "test"), std::invalid_argument);

    auto c = naca_1210();
    for (double chord : {0.0, -0.5, std::numeric_limits<double>::infinity()}) {
        c.chord = chord;
        REQUIRE_THROWS_AS(IO::export_naca65_cascade(c, "test"), std::invalid_argument);
    }
    c = naca_1210();
    c.solidity = 0.0;
    REQUIRE_THROWS_AS(IO::export_naca65_cascade(c, "test"), std::invalid_argument);
    c = naca_1210();
    for (double g : {90.0, -90.0, std::numeric_limits<double>::quiet_NaN()}) {
        c.stagger_deg = g;
        REQUIRE_THROWS_AS(IO::export_naca65_cascade(c, "test"), std::invalid_argument);
    }
    c = naca_1210();
    c.n_stations = 1;
    REQUIRE_THROWS_AS(IO::export_naca65_cascade(c, "test"), std::invalid_argument);
    c = naca_1210();
    c.naca.t_max_over_c = -0.1;
    REQUIRE_THROWS(IO::export_naca65_cascade(c, "test"));
}
