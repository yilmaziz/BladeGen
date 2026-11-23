#include "io/SectionExport.hpp"

#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace BladeGen::IO {

using nlohmann::ordered_json;
using Math::Vec2;

namespace {

ordered_json pair(double a, double b) { return ordered_json::array({a, b}); }
ordered_json triple(double a, double b, double c) { return ordered_json::array({a, b, c}); }

ordered_json units()
{
    return {{"length", "m"}, {"angle", "rad"}, {"mprime", "1 (m' = integral of dm / r, dimensionless)"}};
}

ordered_json generator(const std::string& version) { return {{"name", "BladeGen"}, {"version", version}}; }

/// Index of the TE point in Profile2D::sample_loop(n_edge, n_body): half LE connector
/// (n_edge points), SS body (n_body - 1 more), half TE connector (n_edge - 1 more).
std::size_t loop_te_index(int n_edge, int n_body)
{
    return static_cast<std::size_t>(2 * n_edge + n_body - 3);
}

void check_counts(int n_edge, int n_body, const char* what)
{
    if (n_edge < 2 || n_body < 2) throw std::invalid_argument(std::string(what) + ": need >= 2 points per piece");
}

/// The stream surface written for a section: the whole channel, inlet to outlet.
///  * Designs without their own section streamlines: Blade::streamline_at(span), the
///    channel streamline on which BladeGen places the section itself.
///  * Designs with their own (blade-length) section streamlines, e.g. Rotor 67, whose
///    published section surfaces stop at the blade: the channel's hub curve at span 0, its
///    shroud curve at span 1, the channel streamline in between.  BladeGen's section points
///    are then moved onto that curve (foot_m).
///  * Own streamlines and no channel: the section's own streamline (blade extent only).
struct ExportSurface
{
    Blade::MeridionalCurve curve;
    bool moved = false;  ///< section points moved from their own streamline onto `curve`
};

ExportSurface export_surface(const Blade::BladeDefinition& def, const Blade::Blade& blade, double span)
{
    const bool own = def.stations.front().streamline.has_value();
    if (!own || def.channel.hub().empty()) return {blade.streamline_at(span), false};
    if (span == 0.0) return {def.channel.hub(), true};
    if (span == 1.0) return {def.channel.shroud(), true};
    return {def.channel.streamline(span), true};
}

/// Meridional distance m of the foot point of q on the curve: Gauss-Newton on
/// (P(m) - q) . T(m) = 0, started from the nearest tabulated point.
double foot_m(const Blade::MeridionalCurve& c, const std::vector<double>& m_tab, const std::vector<Vec2>& p_tab,
              const Vec2& q)
{
    std::size_t best = 0;
    for (std::size_t i = 1; i < p_tab.size(); ++i)
        if ((p_tab[i] - q).squaredNorm() < (p_tab[best] - q).squaredNorm()) best = i;
    double m = m_tab[best];
    for (int it = 0; it < 60; ++it) {
        const double dm = (q - c.at_m(m)).dot(c.tangent_at_m(m));
        m += dm;
        if (std::abs(dm) <= 1e-15 * c.length()) break;
    }
    return m;
}

}  // namespace

ordered_json export_design_sections(const Blade::BladeDefinition& def, const DesignExportOptions& opt,
                                    const std::string& source_name, const std::string& generator_version)
{
    check_counts(opt.n_edge, opt.n_body, "export_design_sections");
    if (opt.n_channel < 2) throw std::invalid_argument("export_design_sections: n_channel >= 2");
    if (opt.spans.empty()) throw std::invalid_argument("export_design_sections: no span requested");
    for (double s : opt.spans)
        if (!(s >= 0.0 && s <= 1.0)) throw std::invalid_argument("export_design_sections: span must lie in [0, 1]");
    const Blade::Blade blade(def);  // validates the definition
    const bool own = def.stations.front().streamline.has_value();
    const bool full = !own || !def.channel.hub().empty();
    const bool moved = own && full;

    ordered_json j;
    j["format"] = kSectionsFormat;
    j["version"] = kSectionsVersion;
    j["kind"] = "blade_design";
    j["generator"] = generator(generator_version);
    j["source"] = {{"design_file", source_name}, {"design_name", def.name}};
    j["units"] = units();
    j["convention"] = {
        {"machine_axis", "z, the axis of rotation (BladeGen's internal x axis); r is the radius from it"},
        {"theta", "circumferential angle [rad]: the Cartesian point is (z, r sin(theta), r cos(theta)) in BladeGen's "
                  "internal (x axial, y, z) frame, theta = 0 on its +z axis; in the designer's machine labels "
                  "(Z axial, X = r sin(theta), Y = r cos(theta)) theta turns from +Y toward +X, i.e. clockwise seen "
                  "from downstream (from +Z looking toward -Z)"},
        {"mprime", std::string(full ? "m' = integral of dm / r along the section's streamline, which runs through the "
                                      "whole channel: m = 0 and m' = 0 at the channel inlet (the streamline's first "
                                      "point), not at the LE; "
                                    : "m' = integral of dm / r along the section's streamline, which covers the blade "
                                      "only (the design has no channel): m = 0 and m' = 0 at its first point; ") +
                       "m is the meridional arc length in the (z, r) plane; the section points' mprime_theta, "
                       "mprime_le and mprime_te use the same origin; beyond the streamline's ends (m < 0 or m > its "
                       "length) it continues straight along its end tangent, m' = m'_end + ln(r / r_end) / "
                       "(dr/dm)_end, or (m - m_end) / r_end when dr/dm = 0"},
        {"blade_to_blade",
         std::string("(m', theta): the conformal plane of the stream surface; the point (m', theta) lies at (z, r) = "
                     "streamline(m(m')) with the same theta. ") +
             (moved ? "The design has its own section streamlines over the blade only; the streamline written here is "
                      "the channel's (its hub curve at span 0, its shroud curve at span 1, the channel streamline in "
                      "between), and every BladeGen section point is moved onto it, to its foot point in the (z, r) "
                      "plane with theta kept, i.e. by the distance between the two curves"
                    : "The section points are BladeGen's own, on the streamline written here")},
        {"placement", std::string("BladeGen places the profile by m' = mprime_le + chord_mprime (xi cos(stagger) - s eta "
                                  "sin(stagger)), theta = theta_le + chord_mprime (xi sin(stagger) + s eta "
                                  "cos(stagger)), (xi, eta) the unstaggered unit-chord profile with LE (0, 0) and TE "
                                  "(1, 0), s = +1 (s = -1 when mirror); mprime_le, mprime_te, chord_mprime and stagger "
                                  "are those of the written LE and TE points, m_le_fraction and m_te_fraction their m "
                                  "over the streamline length") +
                          (moved ? " (after the move onto the channel streamline the formula holds only to the "
                                   "distance between the curves)"
                                 : "")},
        {"stagger", "chord angle in the (m', theta) plane from the +m' axis toward +theta [rad]"},
        {"suction_side", "on the +theta side of the chord (the -theta side when mirror is true)"},
        {"point_order", "each side runs from the LE point to the TE point; both sides start at the same LE point "
                        "and end at the same TE point"},
        {"piece", "per point: 0 LE connector, 1 suction-side body, 2 TE connector, 3 pressure-side body"},
        {"pitch_angle", "2 pi / n_blades; the neighbouring blades lie at theta +- pitch_angle"},
    };
    j["n_blades"] = def.n_blades;
    j["pitch_angle"] = 2.0 * Math::kPi / def.n_blades;

    ordered_json channel = nullptr;
    if (!def.channel.hub().empty()) {
        channel = ordered_json::object();
        for (const auto* w : {&def.channel.hub(), &def.channel.shroud()}) {
            ordered_json zr = ordered_json::array();
            for (const auto& p : w->sample(opt.n_channel)) zr.push_back(pair(p.x(), p.y()));
            channel[w == &def.channel.hub() ? "hub" : "shroud"] = {{"n_points", opt.n_channel}, {"zr", zr}};
        }
    }
    j["channel"] = channel;

    const std::size_t ite = loop_te_index(opt.n_edge, opt.n_body);
    ordered_json sections = ordered_json::array();
    for (double span : opt.spans) {
        const auto s = blade.section(span, opt.n_edge, opt.n_body);
        if (s.conformal.size() != 2 * ite + 1) throw std::logic_error("export_design_sections: unexpected loop layout");
        const auto surf = export_surface(def, blade, span);
        const auto& sl = surf.curve;

        ordered_json st;
        st["n_points"] = opt.n_channel;
        ordered_json zr = ordered_json::array(), m = ordered_json::array(), mp = ordered_json::array();
        std::vector<double> m_tab;
        std::vector<Vec2> p_tab;
        for (int i = 0; i < opt.n_channel; ++i) {
            const double mi = sl.length() * i / (opt.n_channel - 1);
            const Vec2 p = sl.at_m(mi);
            m_tab.push_back(mi);
            p_tab.push_back(p);
            zr.push_back(pair(p.x(), p.y()));
            m.push_back(mi);
            mp.push_back(sl.mprime_of_m(mi));
        }
        st["zr"] = zr;
        st["m"] = m;
        st["mprime"] = mp;

        // (m', theta) of every loop point on the written streamline
        std::vector<Vec2> conf(s.conformal.size());
        for (std::size_t k = 0; k < conf.size(); ++k) {
            if (!surf.moved) {
                conf[k] = s.conformal[k];
            } else {  // BladeGen's point on its own streamline, then its foot point on the written one
                const Vec2 q = s.streamline.at_m(s.streamline.m_of_mprime(s.conformal[k].x()));
                conf[k] = {sl.mprime_of_m(foot_m(sl, m_tab, p_tab, q)), s.conformal[k].y()};
            }
        }
        // the loop point k in both coordinate systems (one evaluation of m -> (z, r))
        auto point = [&](std::size_t k, ordered_json& a, ordered_json& b) {
            const Vec2& q = conf[k];
            const Vec2 p = sl.at_m(sl.m_of_mprime(q.x()));
            a = pair(q.x(), q.y());
            b = triple(p.x(), p.y(), q.y());
        };
        auto side = [&](const std::vector<std::size_t>& idx) {
            ordered_json a = ordered_json::array(), b = ordered_json::array(), pc = ordered_json::array();
            for (std::size_t k : idx) {
                ordered_json u, v;
                point(k, u, v);
                a.push_back(u);
                b.push_back(v);
                pc.push_back(s.piece[k]);
            }
            return ordered_json{{"n_points", idx.size()}, {"mprime_theta", a}, {"zrtheta", b}, {"piece", pc}};
        };
        std::vector<std::size_t> ss, ps{0};
        for (std::size_t k = 0; k <= ite; ++k) ss.push_back(k);
        for (std::size_t k = 2 * ite - 1; k >= ite; --k) ps.push_back(k);  // the loop's PS runs TE -> LE

        ordered_json le_a, le_b, te_a, te_b;
        point(0, le_a, le_b);
        point(ite, te_a, te_b);
        // placement values of the written LE / TE points (BladeGen's own when the points are not moved)
        double mp_le = s.mprime_le, mp_te = sl.mprime_of_m(s.placement.m_te * sl.length());
        double chord = s.chord_mprime, stagger = Math::deg2rad(s.placement.stagger_deg);
        double f_le = s.placement.m_le, f_te = s.placement.m_te;
        if (surf.moved) {
            const Vec2 d = conf[ite] - conf[0];
            mp_le = conf[0].x();
            mp_te = conf[ite].x();
            chord = d.norm();
            stagger = std::atan2(d.y(), d.x());
            f_le = sl.m_of_mprime(mp_le) / sl.length();
            f_te = sl.m_of_mprime(mp_te) / sl.length();
        }
        sections.push_back({
            {"span", span},
            {"stagger", stagger},
            {"theta_le", Math::deg2rad(s.placement.theta_le_deg)},
            {"mirror", s.placement.mirror},
            {"m_le_fraction", f_le},
            {"m_te_fraction", f_te},
            {"mprime_le", mp_le},
            {"mprime_te", mp_te},
            {"chord_mprime", chord},
            {"leading_edge", {{"mprime_theta", le_a}, {"zrtheta", le_b}}},
            {"trailing_edge", {{"mprime_theta", te_a}, {"zrtheta", te_b}}},
            {"suction_side", side(ss)},
            {"pressure_side", side(ps)},
            {"streamline", st},
        });
    }
    j["sections"] = sections;
    return j;
}

ordered_json export_naca65_cascade(const CascadeExportOptions& opt, const std::string& generator_version)
{
    check_counts(opt.n_edge, opt.n_stations, "export_naca65_cascade");
    if (!(opt.chord > 0.0) || !std::isfinite(opt.chord)) throw std::invalid_argument("export_naca65_cascade: chord must be > 0");
    if (!(opt.solidity > 0.0) || !std::isfinite(opt.solidity))
        throw std::invalid_argument("export_naca65_cascade: solidity must be > 0");
    if (!(std::abs(opt.stagger_deg) < 90.0)) throw std::invalid_argument("export_naca65_cascade: |stagger| must be below 90 deg");
    const Naca::Naca65Section sec(opt.naca);  // validates the NACA parameters
    const auto& md = sec.metadata();
    const double c = opt.chord, g = Math::deg2rad(opt.stagger_deg), cg = std::cos(g), sg = std::sin(g);
    const double pitch = c / opt.solidity;
    auto place = [&](const Vec2& q) { return pair(c * (q.x() * cg - q.y() * sg), c * (q.x() * sg + q.y() * cg)); };

    ordered_json j;
    j["format"] = kSectionsFormat;
    j["version"] = kSectionsVersion;
    j["kind"] = "naca65_cascade";
    j["generator"] = generator(generator_version);
    j["source"] = {{"naca", {{"designation", md.designation},
                             {"variant", Naca::to_string(opt.naca.variant)},
                             {"cl0", opt.naca.cl0},
                             {"t_max_over_c", opt.naca.t_max_over_c}}}};
    j["units"] = {{"length", "m"}, {"angle", "rad"}};
    j["convention"] = {
        {"frame", "planar cascade: x axial (along the machine axis), y tangential (the pitch direction)"},
        {"profile", "unstaggered unit-chord section (xi, eta): LE point (0, 0), TE point (1, 0), camber and the "
                    "suction side toward +eta"},
        {"placement", "x = chord (xi cos(stagger) - eta sin(stagger)), y = chord (xi sin(stagger) + eta cos(stagger)): "
                      "scaled by the chord and rotated about the LE point, which stays at (0, 0)"},
        {"stagger", "angle of the chord line from the +x (axial) axis toward +y, positive counter-clockwise in the "
                    "(x, y) plane; the TE point is chord (cos(stagger), sin(stagger))"},
        {"pitch", "chord / solidity; the neighbouring blades are the same points translated by "
                  "+- periodic_translation = (0, +- pitch)"},
        {"point_order", "each side runs from the LE point to the TE point; both sides start at the same LE point "
                        "and end at the same TE point"},
        {"piece", "per point: 0 LE fair, 1 suction-side body, 2 TE piece, 3 pressure-side body"},
        {"continuity", "curvature-continuous (G2) where the LE fair and the TE piece join the body, C4 inside the "
                       "body; kovach_sandercock: round TE (G2 fair connector), emery_base: sharp TE"},
    };
    j["cascade"] = {
        {"chord", c},
        {"stagger", g},
        {"stagger_deg", opt.stagger_deg},
        {"solidity", opt.solidity},
        {"pitch", pitch},
        {"periodic_translation", pair(0.0, pitch)},
    };
    j["section"] = {
        {"le_radius", c * md.le_radius},
        {"te_radius", c * md.te_radius},
        {"max_thickness", c * md.max_thickness},
        {"sharp_te", sec.sharp_te()},
        {"n_edge", opt.n_edge},
        {"n_stations", opt.n_stations},
    };
    ordered_json sides = ordered_json::object();
    Vec2 le{0.0, 0.0}, te{1.0, 0.0};
    for (const auto sd : {Profile::Side::Suction, Profile::Side::Pressure}) {
        const auto cs = sec.sample_side(sd, opt.n_edge, opt.n_stations);
        ordered_json xy = ordered_json::array(), xe = ordered_json::array(), pc = ordered_json::array();
        for (std::size_t k = 0; k < cs.xy.size(); ++k) {
            xy.push_back(place(cs.xy[k]));
            xe.push_back(pair(cs.xy[k].x(), cs.xy[k].y()));
            pc.push_back(cs.piece[k]);
        }
        if (sd == Profile::Side::Suction) {
            le = cs.xy.front();
            te = cs.xy.back();
        }
        sides[sd == Profile::Side::Suction ? "suction_side" : "pressure_side"] =
            ordered_json{{"n_points", cs.xy.size()}, {"xy", xy}, {"xi_eta", xe}, {"piece", pc}};
    }
    j["leading_edge"] = {{"xy", place(le)}};
    j["trailing_edge"] = {{"xy", place(te)}};
    j["suction_side"] = sides["suction_side"];
    j["pressure_side"] = sides["pressure_side"];
    return j;
}

// ─────────────────────────────────────────────────────────────────────────────

namespace {

void write_string(std::ostream& o, const std::string& s) { o << ordered_json(s).dump(); }

void write_double(std::ostream& o, double v)
{
    if (!std::isfinite(v)) throw std::invalid_argument("dump_sections: non-finite number");
    std::ostringstream t;
    t.imbue(std::locale::classic());
    t.precision(17);
    t << v;
    std::string s = t.str();
    if (s.find_first_of(".eE") == std::string::npos) s += ".0";  // stays a floating-point number in JSON
    o << s;
}

bool numeric_array(const ordered_json& a)
{
    if (!a.is_array()) return false;
    for (const auto& e : a)
        if (!e.is_number()) return false;
    return true;
}

void write(std::ostream& o, const ordered_json& j, int indent)
{
    const std::string pad(static_cast<std::size_t>(indent + 2), ' '), end(static_cast<std::size_t>(indent), ' ');
    switch (j.type()) {
    case ordered_json::value_t::object: {
        if (j.empty()) { o << "{}"; return; }
        o << "{\n";
        bool first = true;
        for (auto it = j.begin(); it != j.end(); ++it) {
            o << (first ? "" : ",\n") << pad;
            first = false;
            write_string(o, it.key());
            o << ": ";
            write(o, it.value(), indent + 2);
        }
        o << "\n" << end << "}";
        return;
    }
    case ordered_json::value_t::array: {
        if (j.empty()) { o << "[]"; return; }
        if (numeric_array(j)) {  // a point or a table: one line
            o << "[";
            for (std::size_t k = 0; k < j.size(); ++k) {
                if (k) o << ", ";
                write(o, j[k], indent);
            }
            o << "]";
            return;
        }
        o << "[\n";
        for (std::size_t k = 0; k < j.size(); ++k) {
            o << (k ? ",\n" : "") << pad;
            write(o, j[k], indent + 2);
        }
        o << "\n" << end << "]";
        return;
    }
    case ordered_json::value_t::number_float: write_double(o, j.get<double>()); return;
    case ordered_json::value_t::string: write_string(o, j.get<std::string>()); return;
    default: o << j.dump(); return;  // integers, booleans, null
    }
}

}  // namespace

std::string dump_sections(const ordered_json& j)
{
    std::ostringstream o;
    o.imbue(std::locale::classic());
    write(o, j, 0);
    o << "\n";
    return o.str();
}

void write_sections(const ordered_json& j, const std::string& path)
{
    const std::string text = dump_sections(j);
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << text;
    if (!f) throw std::runtime_error("error writing " + path);
}

}  // namespace BladeGen::IO
