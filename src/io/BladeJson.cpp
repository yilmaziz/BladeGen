#include "io/BladeJson.hpp"

#include <fstream>
#include <stdexcept>

namespace BladeGen::IO {

using nlohmann::json;
using Profile::ProfileParams;

json to_json(const ProfileParams& p)
{
    json j = json::object();
    for (const auto& i : Profile::param_table()) j[i.name] = Profile::get_param(p, i.name);
    j["connector"] = Profile::to_string(p.connector);
    j["ps.mprime_ctrl"] = p.ps.mprime_ctrl;
    j["thickness.enforce_c1"] = p.thickness.enforce_c1;
    return j;
}

ProfileParams profile_from_json(const json& j)
{
    if (!j.is_object()) throw std::invalid_argument("profile JSON must be an object");
    ProfileParams p;
    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string& k = it.key();
        if (k.empty() || k[0] == '_') continue;  // metadata
        if (k == "connector") p.connector = Profile::connector_mode_from_string(it->get<std::string>());
        else if (k == "ps.mprime_ctrl") p.ps.mprime_ctrl = it->get<std::vector<double>>();
        else if (k == "thickness.enforce_c1") p.thickness.enforce_c1 = it->get<bool>();
        else Profile::set_param(p, k, it->get<double>());
    }
    Profile::validate(p);
    return p;
}

ProfileParams profile_from_blade2d(const json& j)
{
    static const std::vector<std::pair<const char*, const char*>> map = {
        {"le_radius", "le.radius"}, {"le_stretch_ab", "le.stretch_ab"}, {"le_asym", "le.asym"},
        {"le_beta_deg", "le.beta_deg"}, {"le_d_beta_ss_deg", "le.d_beta_ss_deg"}, {"le_d_beta_ps_deg", "le.d_beta_ps_deg"},
        {"te_radius", "te.radius"}, {"te_stretch_ab", "te.stretch_ab"}, {"te_asym", "te.asym"},
        {"te_beta_deg", "te.beta_deg"}, {"te_d_beta_ss_deg", "te.d_beta_ss_deg"}, {"te_d_beta_ps_deg", "te.d_beta_ps_deg"},
        {"ss_handle_le_frac", "ss.handle_le_frac"}, {"ss_ctrl2_mprime", "ss.ctrl2_mprime"},
        {"ss_ctrl2_theta", "ss.ctrl2_theta"}, {"ss_handle_te_frac", "ss.handle_te_frac"},
        {"th_t_max", "thickness.t_max"}, {"th_mprime_max", "thickness.mprime_max"}, {"th_t_le", "thickness.t_le"},
        {"th_t_te", "thickness.t_te"}, {"th_shape_le", "thickness.shape_le"}, {"th_shape_te", "thickness.shape_te"},
        {"ps_handle_le_frac", "ps.handle_le_frac"}, {"ps_handle_te_frac", "ps.handle_te_frac"},
    };
    ProfileParams p = ProfileParams::naca65_1210_blade2d();  // G2Legacy, blade2d tensions
    for (const auto& [from, to] : map) {
        if (!j.contains(from)) throw std::invalid_argument(std::string("blade2d JSON: missing '") + from + "'");
        Profile::set_param(p, to, j.at(from).get<double>());
    }
    std::vector<double> st;
    for (int i = 0; j.contains("ps_mprime_ctrl_" + std::to_string(i)); ++i)
        st.push_back(j.at("ps_mprime_ctrl_" + std::to_string(i)).get<double>());
    if (!st.empty()) p.ps.mprime_ctrl = st;
    if (j.contains("th_enforce_c1")) p.thickness.enforce_c1 = j.at("th_enforce_c1").get<bool>();
    Profile::validate(p);
    return p;
}

json to_json(const Blade::MeridionalCurve& c)
{
    json pts = json::array();
    for (const auto& q : c.curve().ctrl()) pts.push_back({q.x(), q.y()});
    return {{"degree", c.curve().degree()}, {"knots", c.curve().knots()}, {"control_points", pts}};
}

Blade::MeridionalCurve curve_from_json(const json& j)
{
    auto points = [](const json& a) {
        std::vector<Math::Vec2> v;
        for (const auto& q : a) {
            if (!q.is_array() || q.size() != 2) throw std::invalid_argument("curve JSON: points must be [x, r] pairs");
            v.emplace_back(q[0].get<double>(), q[1].get<double>());
        }
        return v;
    };
    if (j.contains("control_points"))
        return Blade::MeridionalCurve(Math::BSpline2(j.at("degree").get<int>(), j.at("knots").get<std::vector<double>>(),
                                                     points(j.at("control_points"))));
    if (j.contains("points")) return Blade::MeridionalCurve::through_points(points(j.at("points")));
    throw std::invalid_argument("curve JSON needs 'control_points' (+ degree, knots) or 'points'");
}

json to_json(const Blade::BladeDefinition& d)
{
    json st = json::array();
    for (const auto& s : d.stations) {
        json o = {{"span", s.span},
                  {"profile", to_json(s.profile)},
                  {"placement",
                   {{"m_le", s.placement.m_le}, {"m_te", s.placement.m_te}, {"stagger_deg", s.placement.stagger_deg},
                    {"theta_le_deg", s.placement.theta_le_deg}, {"mirror", s.placement.mirror}}}};
        if (s.streamline) o["streamline"] = to_json(*s.streamline);
        st.push_back(o);
    }
    json j = {{"name", d.name},
              {"n_blades", d.n_blades},
              {"tip_clearance", d.tip_clearance},
              {"interpolation", d.interpolation == Blade::SpanInterpolation::Cubic ? "cubic" : "linear"},
              {"span_extension", d.span_extension},
              {"n_sections", d.n_sections},
              {"stations", st}};
    if (!d.channel.hub().empty()) j["channel"] = {{"hub", to_json(d.channel.hub())}, {"shroud", to_json(d.channel.shroud())}};
    return j;
}

Blade::BladeDefinition blade_from_json(const json& j)
{
    Blade::BladeDefinition d;
    d.name = j.value("name", std::string("blade"));
    d.n_blades = j.value("n_blades", 20);
    d.tip_clearance = j.value("tip_clearance", 0.0);
    const std::string interp = j.value("interpolation", std::string("cubic"));
    if (interp != "cubic" && interp != "linear") throw std::invalid_argument("interpolation must be 'cubic' or 'linear'");
    d.interpolation = interp == "cubic" ? Blade::SpanInterpolation::Cubic : Blade::SpanInterpolation::Linear;
    d.span_extension = j.value("span_extension", 0.04);
    d.n_sections = j.value("n_sections", 15);
    if (j.contains("channel"))
        d.channel = Blade::MeridionalChannel(curve_from_json(j.at("channel").at("hub")), curve_from_json(j.at("channel").at("shroud")));
    if (!j.contains("stations")) throw std::invalid_argument("blade JSON: 'stations' missing");
    for (const auto& s : j.at("stations")) {
        Blade::BladeStation b;
        b.span = s.at("span").get<double>();
        b.profile = profile_from_json(s.at("profile"));
        const auto& pl = s.at("placement");
        b.placement.m_le = pl.at("m_le").get<double>();
        b.placement.m_te = pl.at("m_te").get<double>();
        b.placement.stagger_deg = pl.at("stagger_deg").get<double>();
        b.placement.theta_le_deg = pl.value("theta_le_deg", 0.0);
        b.placement.mirror = pl.value("mirror", false);
        if (s.contains("streamline")) b.streamline = curve_from_json(s.at("streamline"));
        d.stations.push_back(std::move(b));
    }
    d.validate();
    return d;
}

namespace {
json read_file(const std::string& path)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    try {
        return json::parse(f);
    } catch (const json::parse_error& e) {
        throw std::invalid_argument(path + ": " + e.what());
    }
}
void write_file(const json& j, const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << j.dump(2) << "\n";
}
}  // namespace

void write_blade_json(const Blade::BladeDefinition& d, const std::string& path) { write_file(to_json(d), path); }
Blade::BladeDefinition read_blade_json(const std::string& path) { return blade_from_json(read_file(path)); }
void write_profile_json(const ProfileParams& p, const std::string& path) { write_file(to_json(p), path); }

ProfileParams read_profile_json(const std::string& path)
{
    const auto j = read_file(path);
    return j.contains("le_radius") ? profile_from_blade2d(j) : profile_from_json(j);
}

}  // namespace BladeGen::IO
