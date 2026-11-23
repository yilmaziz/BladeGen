#include "design/BladeDesign.hpp"

#include "fit/InverseFit.hpp"
#include "fit/NacaFit.hpp"
#include "io/BladeJson.hpp"
#include "math/Interpolation.hpp"
#include "math/Numerics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <future>
#include <sstream>
#include <stdexcept>

namespace BladeGen::Design {

namespace {

using json = nlohmann::json;
constexpr const char* kFormat = "bladegen-design-2";
constexpr const char* kFormatV1 = "bladegen-design-1";
constexpr std::size_t kCacheSize = 48;

void require(bool ok, const std::string& msg)
{
    if (!ok) throw std::invalid_argument("BladeDesign: " + msg);
}

json naca_json(const Naca::Naca65Params& p)
{
    return {{"cl0", p.cl0}, {"t_max_over_c", p.t_max_over_c}, {"variant", Naca::to_string(p.variant)}};
}

Naca::Naca65Params naca_from(const json& j)
{
    Naca::Naca65Params p;
    p.cl0 = j.value("cl0", p.cl0);
    p.t_max_over_c = j.value("t_max_over_c", p.t_max_over_c);
    p.variant = Naca::variant_from_string(j.value("variant", Naca::to_string(p.variant)));
    return p;
}

json wall_json(const Blade::WallPoints& w, Blade::WallLaw law)
{
    return {{"z_le", w.z_le}, {"r_le", w.r_le}, {"z_te", w.z_te}, {"r_te", w.r_te}, {"law", Blade::to_string(law)}};
}

void wall_from(const json& j, Blade::WallPoints& w, Blade::WallLaw& law)
{
    w.z_le = j.value("z_le", w.z_le);
    w.r_le = j.value("r_le", w.r_le);
    w.z_te = j.value("z_te", w.z_te);
    w.r_te = j.value("r_te", w.r_te);
    law = Blade::wall_law_from_string(j.value("law", Blade::to_string(law)));
}

std::string cache_key(const Layer& l, const MeshOptions& m)
{
    std::string k = to_string(l.method) + "|" + naca_json(l.naca).dump() + "|" + std::to_string(m.n_edge) + "," +
                    std::to_string(m.n_body);
    if (l.uses_custom()) k += "|" + IO::to_json(l.custom).dump();
    if (l.fitted_to) k += "|fit" + naca_json(*l.fitted_to).dump();
    return k;
}

std::string naca_key(const Naca::Naca65Params& p) { return naca_json(p).dump(); }

/// Linear interpolation of y(x) at xq (x increasing; clamped).
double lerp_table(const std::vector<double>& x, const std::vector<double>& y, double xq)
{
    if (x.size() == 1) return y.front();
    if (xq <= x.front()) return y.front();
    if (xq >= x.back()) return y.back();
    const auto it = std::upper_bound(x.begin(), x.end(), xq);
    const auto k = static_cast<std::size_t>(it - x.begin());
    const double w = (xq - x[k - 1]) / (x[k] - x[k - 1]);
    return y[k - 1] + w * (y[k] - y[k - 1]);
}

}  // namespace

std::string to_string(Method m)
{
    switch (m) {
    case Method::Naca65: return "NACA 65";
    case Method::Custom: return "Custom";
    case Method::Naca65Custom: return "NACA 65 + Custom";
    }
    return "?";
}

Method method_from_string(const std::string& s)
{
    if (s == "NACA 65" || s == "naca65") return Method::Naca65;
    if (s == "Custom" || s == "custom") return Method::Custom;
    if (s == "NACA 65 + Custom" || s == "naca65_custom") return Method::Naca65Custom;
    throw std::invalid_argument("unknown construction method '" + s + "'");
}

namespace {
#include "design/Rotor67Stations.inc"
}  // namespace

const std::vector<Rotor67Layer>& rotor67_sections()
{
    static const std::vector<Rotor67Layer> table = [] {
        std::vector<Rotor67Layer> t;
        for (const char* txt : kRotor67Stations) {
            const auto j = json::parse(txt);
            Rotor67Layer l;
            l.span = j.at("span").get<double>();
            l.stagger_deg = j.at("stagger_deg").get<double>();
            l.theta_le_deg = j.at("theta_le_deg").get<double>();
            l.profile = IO::profile_from_json(j.at("profile"));
            t.push_back(std::move(l));
        }
        return t;
    }();
    return table;
}

Rotor67Layer rotor67_layer(double span)
{
    if (!std::isfinite(span) || span < 0.0 || span > 1.0) throw std::invalid_argument("rotor67_layer: span must be in [0, 1]");
    const auto& t = rotor67_sections();
    std::size_t k = 1;
    while (k + 1 < t.size() && t[k].span < span) ++k;
    const auto& a = t[k - 1];
    const auto& b = t[k];
    const double w = std::clamp((span - a.span) / (b.span - a.span), 0.0, 1.0);
    Rotor67Layer l;
    l.span = span;
    l.stagger_deg = a.stagger_deg + w * (b.stagger_deg - a.stagger_deg);
    l.theta_le_deg = a.theta_le_deg + w * (b.theta_le_deg - a.theta_le_deg);
    l.profile = w == 0.0 ? a.profile : w == 1.0 ? b.profile : Profile::blend(a.profile, b.profile, w);
    return l;
}

BladeDesign BladeDesign::defaults()
{
    BladeDesign d = naca65_defaults();
    d.name = "rotor67_layers";
    for (auto& l : d.layers) {
        const auto r = rotor67_layer(l.span);
        l.method = Method::Custom;
        l.custom = r.profile;
        l.stagger_deg = r.stagger_deg;
        l.theta_le_deg = r.theta_le_deg;
    }
    return d;
}

BladeDesign BladeDesign::naca65_defaults()
{
    BladeDesign d;
    d.name = "axial_rotor";
    const char* names[] = {"Hub", "25 %", "Mid", "75 %", "Tip"};
    const double cl0[] = {1.5, 1.35, 1.2, 1.0, 0.8};
    const double tc[] = {0.10, 0.09, 0.08, 0.07, 0.06};
    const double stagger[] = {25.0, 31.0, 37.0, 42.5, 48.0};
    const double sol[] = {1.6, 1.45, 1.3, 1.15, 1.0};
    for (int k = 0; k < 5; ++k) {
        Layer l;
        l.name = names[k];
        l.span = 0.25 * k;
        l.naca.cl0 = cl0[k];
        l.naca.t_max_over_c = tc[k];
        l.stagger_deg = stagger[k];
        l.solidity = sol[k];
        d.layers.push_back(l);
    }
    return d;
}

std::string to_string(Restagger::Axis a)
{
    switch (a) {
    case Restagger::Axis::Hub: return "hub";
    case Restagger::Axis::Mid: return "mid";
    case Restagger::Axis::Tip: return "tip";
    case Restagger::Axis::Point: return "point";
    }
    return "hub";
}

Restagger::Axis restagger_axis_from_string(const std::string& s)
{
    if (s == "hub") return Restagger::Axis::Hub;
    if (s == "mid") return Restagger::Axis::Mid;
    if (s == "tip") return Restagger::Axis::Tip;
    if (s == "point") return Restagger::Axis::Point;
    throw std::invalid_argument("unknown re-stagger axis '" + s + "' (hub, mid, tip, point)");
}

std::string to_string(Stacking::Mode m) { return m == Stacking::Mode::Edges ? "edges" : "reference_point"; }
std::string to_string(Stacking::Variant v) { return v == Stacking::Variant::SameZTheta ? "same_z_theta" : "closest_point"; }
std::string to_string(Stacking::Reference r)
{
    switch (r) {
    case Stacking::Reference::CG: return "cg";
    case Stacking::Reference::LE: return "le";
    case Stacking::Reference::TE: return "te";
    }
    return "cg";
}

Stacking::Mode stacking_mode_from_string(const std::string& s)
{
    if (s == "edges") return Stacking::Mode::Edges;
    if (s == "reference_point") return Stacking::Mode::ReferencePoint;
    throw std::invalid_argument("unknown stacking mode '" + s + "' (edges, reference_point)");
}

Stacking::Variant stacking_variant_from_string(const std::string& s)
{
    if (s == "same_z_theta") return Stacking::Variant::SameZTheta;
    if (s == "closest_point") return Stacking::Variant::ClosestPoint;
    throw std::invalid_argument("unknown stacking variant '" + s + "' (same_z_theta, closest_point)");
}

Stacking::Reference stacking_reference_from_string(const std::string& s)
{
    if (s == "cg") return Stacking::Reference::CG;
    if (s == "le") return Stacking::Reference::LE;
    if (s == "te") return Stacking::Reference::TE;
    throw std::invalid_argument("unknown stacking reference '" + s + "' (cg, le, te)");
}

void BladeDesign::validate() const
{
    meridional.validate();
    require(!layers.empty(), "need at least one layer");
    for (std::size_t i = 0; i < layers.size(); ++i) {
        const auto& l = layers[i];
        require(std::isfinite(l.span) && l.span >= 0.0 && l.span <= 1.0, "layer spans must be in [0, 1]");
        if (i > 0) require(l.span > layers[i - 1].span, "layer spans must increase");
        require(std::abs(l.stagger_deg) < 80.0, "|stagger| must be below 80 deg");
        require(std::abs(l.stagger_deg + restagger.angle_deg) < 80.0, "|stagger + re-stagger| must be below 80 deg");
        require(l.solidity > 0.0, "solidity must be > 0");
        require(std::isfinite(l.chord_scale) && l.chord_scale >= 0.2 && l.chord_scale <= 5.0, "chord scale in [0.2, 5]");
        l.naca.validate();
        if (l.uses_custom()) Profile::validate(l.custom);
    }
    require(n_blades >= 1 && n_blades <= 500, "n_blades in [1, 500]");
    require(std::isfinite(restagger.angle_deg) && std::abs(restagger.angle_deg) <= 30.0, "re-stagger in [-30, 30] deg");
    require(std::isfinite(restagger.chord_fraction) && restagger.chord_fraction >= 0.0 && restagger.chord_fraction <= 1.0,
            "re-stagger pivot chord fraction in [0, 1]");
    require(std::isfinite(restagger.z) && std::isfinite(restagger.theta_deg), "re-stagger point must be finite");
    require(std::isfinite(stacking.z_ref) && std::isfinite(stacking.theta_ref_deg) && std::abs(stacking.theta_ref_deg) <= 180.0,
            "stacking reference point: finite Z, theta in [-180, 180] deg");
    require(mesh.n_edge >= 4 && mesh.n_body >= 8, "mesh: n_edge >= 4, n_body >= 8");
    require(mesh.n_rows >= 2 && mesh.n_rows <= 400, "mesh: n_rows in [2, 400]");
    require(mesh.n_wall_m >= 4 && mesh.n_wall_theta >= 2 && mesh.n_field_chord >= 3, "mesh: wall / field resolution");
    require(mesh.row_clustering >= 0.0 && mesh.row_clustering <= 10.0, "mesh: row clustering in [0, 10]");
}

// ── JSON ─────────────────────────────────────────────────────────────────────

json to_json(const BladeDesign& d)
{
    json j;
    j["format"] = kFormat;
    j["name"] = d.name;
    j["n_blades"] = d.n_blades;
    j["mirror"] = d.mirror;
    j["restagger"] = {{"angle_deg", d.restagger.angle_deg}, {"axis", to_string(d.restagger.axis)},
                      {"chord_fraction", d.restagger.chord_fraction}, {"z", d.restagger.z}, {"theta_deg", d.restagger.theta_deg}};
    j["stacking"] = {{"mode", to_string(d.stacking.mode)}, {"variant", to_string(d.stacking.variant)},
                     {"reference", to_string(d.stacking.reference)}, {"z_ref", d.stacking.z_ref},
                     {"theta_ref_deg", d.stacking.theta_ref_deg}};
    const auto& m = d.meridional;
    j["meridional"] = {{"hub", wall_json(m.hub, m.hub_law)},
                       {"casing", wall_json(m.casing, m.casing_law)},
                       {"inlet_extension", m.inlet_extension},
                       {"outlet_extension", m.outlet_extension}};
    j["mesh"] = {{"n_edge", d.mesh.n_edge}, {"n_body", d.mesh.n_body}, {"n_rows", d.mesh.n_rows},
                 {"row_clustering", d.mesh.row_clustering}, {"n_wall_m", d.mesh.n_wall_m},
                 {"n_wall_theta", d.mesh.n_wall_theta}, {"n_field_chord", d.mesh.n_field_chord}};
    j["layers"] = json::array();
    for (const auto& l : d.layers) {
        json q = {{"name", l.name}, {"span", l.span}, {"method", to_string(l.method)}, {"naca", naca_json(l.naca)},
                  {"custom", IO::to_json(l.custom)}, {"stagger_deg", l.stagger_deg},
                  {"theta_le_deg", l.theta_le_deg}, {"solidity", l.solidity}, {"chord_scale", l.chord_scale},
                  {"fit_rms", l.fit_rms}};
        q["fitted_to"] = l.fitted_to ? naca_json(*l.fitted_to) : json(nullptr);
        j["layers"].push_back(q);
    }
    return j;
}

namespace {

Layer layer_from(const json& q)
{
    Layer l;
    l.name = q.value("name", l.name);
    l.span = q.value("span", l.span);
    l.method = method_from_string(q.value("method", to_string(l.method)));
    if (q.contains("naca")) l.naca = naca_from(q["naca"]);
    if (q.contains("custom") && q["custom"].is_object()) l.custom = IO::profile_from_json(q["custom"]);
    l.stagger_deg = q.value("stagger_deg", l.stagger_deg);
    l.theta_le_deg = q.value("theta_le_deg", l.theta_le_deg);
    l.solidity = q.value("solidity", l.solidity);
    l.chord_scale = q.value("chord_scale", l.chord_scale);
    if (q.contains("fit_rms") && q["fit_rms"].is_number()) l.fit_rms = q["fit_rms"].get<double>();
    if (q.contains("fitted_to") && q["fitted_to"].is_object()) l.fitted_to = naca_from(q["fitted_to"]);
    return l;
}

/// v1.0.1 station (Python section designer).
Layer layer_from_v1(const json& q)
{
    Layer l;
    l.name = q.value("name", l.name);
    l.span = q.value("span", l.span);
    l.method = method_from_string(q.value("method", to_string(l.method)));
    l.naca.cl0 = q.value("cl0", l.naca.cl0);
    l.naca.t_max_over_c = q.value("t_max_over_c", l.naca.t_max_over_c);
    l.naca.variant = Naca::variant_from_string(q.value("variant", Naca::to_string(l.naca.variant)));
    if (q.contains("custom") && q["custom"].is_object()) l.custom = IO::profile_from_json(q["custom"]);
    l.stagger_deg = q.value("stagger_deg", l.stagger_deg);
    l.solidity = q.value("solidity", l.solidity);
    if (q.contains("fit_rms") && q["fit_rms"].is_number()) l.fit_rms = q["fit_rms"].get<double>();
    if (q.contains("fitted_to") && q["fitted_to"].is_array() && q["fitted_to"].size() == 3) {
        Naca::Naca65Params f;
        f.cl0 = q["fitted_to"][0].get<double>();
        f.t_max_over_c = q["fitted_to"][1].get<double>();
        f.variant = Naca::variant_from_string(q["fitted_to"][2].get<std::string>());
        l.fitted_to = f;
    }
    return l;
}

/// A layer between two others: NACA parameters and custom parameters blended.
Layer blend_layers(const Layer& a, const Layer& b, double w, const std::string& name, double span)
{
    Layer l = w < 0.5 ? a : b;
    l.name = name;
    l.span = span;
    l.naca.cl0 = a.naca.cl0 + w * (b.naca.cl0 - a.naca.cl0);
    l.naca.t_max_over_c = a.naca.t_max_over_c + w * (b.naca.t_max_over_c - a.naca.t_max_over_c);
    l.custom = Profile::blend(a.custom, b.custom, w);
    l.stagger_deg = a.stagger_deg + w * (b.stagger_deg - a.stagger_deg);
    l.theta_le_deg = a.theta_le_deg + w * (b.theta_le_deg - a.theta_le_deg);
    l.solidity = a.solidity + w * (b.solidity - a.solidity);
    if (a.uses_custom() && b.uses_custom()) l.method = Method::Custom;
    else if (!a.uses_custom() || !b.uses_custom()) l.method = Method::Naca65;
    l.fitted_to.reset();
    l.fit_rms = 0.0;
    return l;
}

}  // namespace

BladeDesign design_from_json(const json& j)
{
    const std::string fmt = j.value("format", std::string(kFormat));
    BladeDesign d = BladeDesign::defaults();
    if (fmt == kFormatV1) {
        const auto& st = j.at("stations");
        if (st.size() != 3) throw std::invalid_argument("design v1: expected 3 stations (hub, mid, tip)");
        const Layer hub = layer_from_v1(st[0]), mid = layer_from_v1(st[1]), tip = layer_from_v1(st[2]);
        d.layers = {hub, blend_layers(hub, mid, 0.5, "25 %", 0.25), mid, blend_layers(mid, tip, 0.5, "75 %", 0.75), tip};
        d.layers[0].span = 0.0;
        d.layers[2].span = 0.5;
        d.layers[4].span = 1.0;
        d.validate();
        return d;
    }
    if (fmt != kFormat) throw std::invalid_argument("unknown design format '" + fmt + "'");
    d.name = j.value("name", d.name);
    d.n_blades = j.value("n_blades", d.n_blades);
    d.mirror = j.value("mirror", d.mirror);
    if (j.contains("restagger")) {
        const auto& q = j.at("restagger");
        auto& R = d.restagger;
        R.angle_deg = q.value("angle_deg", R.angle_deg);
        R.axis = restagger_axis_from_string(q.value("axis", to_string(R.axis)));
        R.chord_fraction = q.value("chord_fraction", R.chord_fraction);
        R.z = q.value("z", R.z);
        R.theta_deg = q.value("theta_deg", R.theta_deg);
    }
    if (j.contains("stacking")) {  // designs saved before v1.0.4: the prescribed edges
        const auto& q = j.at("stacking");
        auto& St = d.stacking;
        St.mode = stacking_mode_from_string(q.value("mode", to_string(St.mode)));
        St.variant = stacking_variant_from_string(q.value("variant", to_string(St.variant)));
        St.reference = stacking_reference_from_string(q.value("reference", to_string(St.reference)));
        St.z_ref = q.value("z_ref", St.z_ref);
        St.theta_ref_deg = q.value("theta_ref_deg", St.theta_ref_deg);
    }
    if (j.contains("meridional")) {
        const auto& m = j["meridional"];
        if (m.contains("hub")) wall_from(m["hub"], d.meridional.hub, d.meridional.hub_law);
        if (m.contains("casing")) wall_from(m["casing"], d.meridional.casing, d.meridional.casing_law);
        d.meridional.inlet_extension = m.value("inlet_extension", d.meridional.inlet_extension);
        d.meridional.outlet_extension = m.value("outlet_extension", d.meridional.outlet_extension);
    }
    if (j.contains("mesh")) {
        const auto& m = j["mesh"];
        d.mesh.n_edge = m.value("n_edge", d.mesh.n_edge);
        d.mesh.n_body = m.value("n_body", d.mesh.n_body);
        d.mesh.n_rows = m.value("n_rows", d.mesh.n_rows);
        d.mesh.row_clustering = m.value("row_clustering", d.mesh.row_clustering);
        d.mesh.n_wall_m = m.value("n_wall_m", d.mesh.n_wall_m);
        d.mesh.n_wall_theta = m.value("n_wall_theta", d.mesh.n_wall_theta);
        d.mesh.n_field_chord = m.value("n_field_chord", d.mesh.n_field_chord);
    }
    if (j.contains("layers")) {
        d.layers.clear();
        for (const auto& q : j["layers"]) d.layers.push_back(layer_from(q));
    }
    d.validate();
    return d;
}

void write_design(const BladeDesign& d, const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << to_json(d).dump(2) << "\n";
}

BladeDesign read_design(const std::string& path)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    return design_from_json(json::parse(f));
}

// ── interpolation helpers ────────────────────────────────────────────────────

std::vector<double> spline_weights(const std::vector<double>& spans, double s)
{
    const std::size_t n = spans.size();
    if (n == 0) throw std::invalid_argument("spline_weights: no spans");
    std::vector<double> w(n, 0.0);
    if (n == 1) {
        w[0] = 1.0;
        return w;
    }
    for (std::size_t i = 1; i < n; ++i)
        if (!(spans[i] > spans[i - 1])) throw std::invalid_argument("spline_weights: spans must increase");
    if (n == 2) {
        const double t = (s - spans[0]) / (spans[1] - spans[0]);
        w[0] = 1.0 - t;
        w[1] = t;
        return w;
    }
    const double sq = std::clamp(s, spans.front(), spans.back());
    for (std::size_t k = 0; k < n; ++k) {
        std::vector<Math::Vec<1>> e(n, Math::Vec<1>(0.0));
        e[k] = Math::Vec<1>(1.0);
        const auto c = Math::interpolate_cubic<1>(spans, e, Math::EndCondition<1>::natural(), Math::EndCondition<1>::natural());
        w[k] = c.eval(sq)[0];
    }
    return w;
}

std::vector<double> tanh_two_sided(int n, double d)
{
    if (n < 2) throw std::invalid_argument("tanh_two_sided: n >= 2");
    if (!(d >= 0.0) || !std::isfinite(d)) throw std::invalid_argument("tanh_two_sided: d must be finite and >= 0");
    std::vector<double> u(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double xi = static_cast<double>(i) / (n - 1);
        u[static_cast<std::size_t>(i)] = d < 1e-8 ? xi : 0.5 * (1.0 + std::tanh(d * (xi - 0.5)) / std::tanh(0.5 * d));
    }
    u.front() = 0.0;
    u.back() = 1.0;
    return u;
}

// ── 2-D sections ─────────────────────────────────────────────────────────────

Section2D build_section2d(const Layer& layer, const Naca::Naca65Section* naca, const MeshOptions& mesh)
{
    Section2D s;
    if (layer.uses_custom()) {
        const Profile::Profile2D prof(layer.custom);
        s.quality = Profile::check_quality(prof);
        s.loop = prof.sample_loop(mesh.n_edge, mesh.n_body);
        s.designation = layer.method == Method::Custom ? "Custom" : "NACA 65 + Custom (fitted)";
        s.le_radius = layer.custom.le.radius;
        s.control = Section2D::Control{prof.le_connector().ctrl(), prof.suction_side().ctrl(), prof.te_connector().ctrl(),
                                       prof.pressure_side().ctrl()};
        s.thickness_law = prof.thickness(120);
        s.xi_max = prof.xi_max();
        s.te_radius = layer.custom.te.radius;
        if (layer.method == Method::Naca65Custom && naca) {
            const auto nq = Naca::check_quality(*naca);
            std::vector<Vec2> ss, ps;
            for (std::size_t i = 0; i < nq.ss.xy.size(); i += 3) ss.push_back(nq.ss.xy[i]);
            for (std::size_t i = 0; i < nq.ps.xy.size(); i += 3) ps.push_back(nq.ps.xy[i]);
            std::vector<double> dss, dps;
            Fit::side_distances(prof, ss, ps, dss, dps);
            double acc = 0.0;
            for (double v : dss) acc += v * v;
            for (double v : dps) acc += v * v;
            s.deviation_rms = std::sqrt(acc / static_cast<double>(dss.size() + dps.size()));
        }
    } else {
        if (!naca) throw std::logic_error("build_section2d: NACA 65 layer without its section");
        s.quality = Naca::check_quality(*naca);
        s.loop = naca->sample_loop(mesh.n_edge, mesh.n_body);
        s.designation = naca->metadata().designation + " (" + Naca::to_string(layer.naca.variant) + ")";
        s.le_radius = naca->metadata().le_radius;
        s.te_radius = naca->metadata().te_radius;
        s.equivalent_arc_deg = std::array<double, 2>{Math::rad2deg(naca->metadata().chi_le), Math::rad2deg(naca->metadata().chi_te)};
    }
    s.camber = Profile::camber_line(s.quality.ss.xy, s.quality.ps.xy, 101);
    s.area = std::abs(Profile::polygon_area(s.loop.xy));
    return s;
}

std::shared_ptr<const Naca::Naca65Section> Evaluator::naca(const Naca::Naca65Params& p)
{
    const auto key = naca_key(p);
    {
        std::lock_guard<std::mutex> g(mu_);
        if (auto it = naca_.find(key); it != naca_.end()) return it->second;
    }
    auto sec = std::make_shared<const Naca::Naca65Section>(p);  // built outside the lock (parallel layers)
    std::lock_guard<std::mutex> g(mu_);
    if (naca_.size() > kCacheSize) naca_.clear();
    return naca_.emplace(key, std::move(sec)).first->second;
}

std::shared_ptr<const Section2D> Evaluator::section(const Layer& layer, const MeshOptions& mesh)
{
    const auto key = cache_key(layer, mesh);
    {
        std::lock_guard<std::mutex> g(mu_);
        if (auto it = sections_.find(key); it != sections_.end()) return it->second;
    }
    std::shared_ptr<const Naca::Naca65Section> ns;
    if (!layer.uses_custom() || layer.method == Method::Naca65Custom) ns = naca(layer.naca);
    auto sec = std::make_shared<const Section2D>(build_section2d(layer, ns.get(), mesh));
    std::lock_guard<std::mutex> g(mu_);
    if (!sections_.count(key)) {
        sections_[key] = sec;
        order_.push_back(key);
        while (order_.size() > kCacheSize) {
            sections_.erase(order_.front());
            order_.erase(order_.begin());
        }
    }
    return sec;
}

void Evaluator::clear()
{
    std::lock_guard<std::mutex> g(mu_);
    sections_.clear();
    naca_.clear();
    order_.clear();
}

// ── centres of gravity ───────────────────────────────────────────────────────

PolygonCentroid polygon_centroid(const std::vector<Vec2>& p)
{
    if (p.size() < 3) throw std::invalid_argument("polygon_centroid: need at least 3 points");
    double a2 = 0.0;
    Vec2 c = Vec2::Zero();
    const Vec2 o = p.front();  // relative to the first point: no cancellation far from the origin
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Vec2 u = p[i] - o, v = p[(i + 1) % p.size()] - o;
        const double cr = u.x() * v.y() - v.x() * u.y();
        a2 += cr;
        c += cr * (u + v);
    }
    if (!(std::abs(a2) > 0.0)) throw std::invalid_argument("polygon_centroid: zero area");
    PolygonCentroid r;
    r.area = 0.5 * std::abs(a2);
    r.centroid = o + c / (3.0 * a2);
    return r;
}

MassProperties mass_properties(const Grid& g)
{
    if (g.ni < 7 || g.nj < 2 || g.p.size() != static_cast<std::size_t>(g.ni * g.nj))
        throw std::invalid_argument("mass_properties: need closed row loops (ni >= 7, nj >= 2)");
    double V = 0.0;
    Vec3 M = Vec3::Zero();
    const Vec3 o = g.at(0, 0);  // tetrahedra from a point near the blade: no cancellation
    auto tri = [&](const Vec3& a, const Vec3& b, const Vec3& c) {
        const Vec3 pa = a - o, pb = b - o, pc = c - o;
        const double v = pa.dot(pb.cross(pc)) / 6.0;
        V += v;
        M += v * (pa + pb + pc) / 4.0;
    };
    const int ni = g.ni, nj = g.nj, N = ni - 1;
    for (int j = 0; j + 1 < nj; ++j)
        for (int i = 0; i + 1 < ni; ++i) {
            tri(g.at(i, j), g.at(i + 1, j), g.at(i + 1, j + 1));
            tri(g.at(i, j), g.at(i + 1, j + 1), g.at(i, j + 1));
        }
    // caps: strips between the points i and N - i (suction / pressure side), oriented against the sides
    for (int i = 0; i < N / 2; ++i) {
        const int b = i + 1, c = N - i - 1, e = N - i;
        if (c > b) {
            tri(g.at(i, 0), g.at(c, 0), g.at(b, 0));
            tri(g.at(i, nj - 1), g.at(b, nj - 1), g.at(c, nj - 1));
        }
        tri(g.at(i, 0), g.at(e, 0), g.at(c, 0));
        tri(g.at(i, nj - 1), g.at(c, nj - 1), g.at(e, nj - 1));
    }
    if (!(std::abs(V) > 0.0)) throw std::invalid_argument("mass_properties: zero volume");
    MassProperties r;
    r.volume = std::abs(V);
    r.centroid = o + M / V;
    return r;
}

// ── placement ────────────────────────────────────────────────────────────────

Vec2 Placement::operator()(const Vec2& q) const
{
    const double cg = std::cos(gamma), sn = std::sin(gamma);
    return {m0 + c * (q.x() * cg - sign * q.y() * sn), t0 + c * (q.x() * sn + sign * q.y() * cg)};
}

Placement restaggered(const Placement& p, double d_deg, const Vec2& pivot)
{
    if (d_deg == 0.0) return p;
    const double d = Math::deg2rad(d_deg);
    Placement q = p;
    q.gamma = p.gamma + d;
    const Vec2 v = Vec2(p.m0, p.t0) - pivot;
    q.m0 = pivot.x() + std::cos(d) * v.x() - std::sin(d) * v.y();
    q.t0 = pivot.y() + std::sin(d) * v.x() + std::cos(d) * v.y();
    return q;
}

Placement placement(double stagger_deg, double theta_le_deg, double mprime_te, bool mirror, double restagger_deg, const Vec2& pivot)
{
    Placement p;
    const double g0 = Math::deg2rad(stagger_deg);
    p.gamma = g0;
    p.c = mprime_te / std::cos(g0);
    p.sign = mirror ? -1.0 : 1.0;
    p.m0 = 0.0;
    p.t0 = Math::deg2rad(theta_le_deg);
    return restaggered(p, restagger_deg, pivot);
}

double mprime_at_x(const Blade::StreamSurface& S, double x)
{
    const double a = S.mprime_inlet(), b = S.mprime_outlet();
    const auto f = [&](double mp) { return S.xr(mp).x() - x; };
    if ((f(a) > 0.0) == (f(b) > 0.0)) {
        std::ostringstream m;
        m << "Z = " << x << " mm is outside the stream surface at span " << S.span() << " (Z " << S.xr(a).x() << " .. "
          << S.xr(b).x() << ")";
        throw std::invalid_argument(m.str());
    }
    return Math::brent_root(f, a, b, 1e-13);
}

double closest_mprime(const Blade::StreamSurface& S, const Vec2& q)
{
    const double a = S.mprime_inlet(), b = S.mprime_outlet();
    constexpr int n = 400;
    const auto at = [&](int i) { return a + (b - a) * i / n; };
    const auto f = [&](double mp) {
        const Vec2 e = S.xr(mp) - q;
        return e.dot(e);
    };
    int best = 0;
    double fb = f(a);
    for (int i = 1; i <= n; ++i)
        if (const double v = f(at(i)); v < fb) {
            fb = v;
            best = i;
        }
    if (best == 0 || best == n) {
        std::ostringstream m;
        m << "the point of the stream surface at span " << S.span() << " closest to (Z " << q.x() << ", R " << q.y()
          << ") is an end of the channel";
        throw std::invalid_argument(m.str());
    }
    // d/dm' |xr - q|^2 / 2 = (xr - q) . dxr: negative before the minimum, positive after it
    const auto g = [&](double mp) { return (S.xr(mp) - q).dot(S.dxr(mp)); };
    const double lo = at(best - 1), hi = at(best + 1);
    if ((g(lo) > 0.0) == (g(hi) > 0.0)) throw std::invalid_argument("closest_mprime: the closest point is not bracketed");
    return Math::brent_root(g, lo, hi, 1e-14);
}

double StackingLine::x_at(double span) const
{
    if (x.empty() || x.size() != spans.size()) throw std::logic_error("StackingLine: no reference points");
    if (variant == Stacking::Variant::SameZTheta) return x.front();
    const auto w = spline_weights(spans, span);
    double v = 0.0;
    for (std::size_t k = 0; k < x.size(); ++k) v += w[k] * x[k];
    return v;
}

StackingLine resolve_stacking(const BladeDesign& d, const Blade::FlowPath& fp)
{
    StackingLine st;
    const auto& in = d.stacking;
    st.mode = in.mode;
    st.variant = in.variant;
    st.reference = in.reference;
    st.theta = Math::deg2rad(in.theta_ref_deg);
    for (const auto& l : d.layers) st.spans.push_back(l.span);
    if (st.mode == Stacking::Mode::Edges) return st;
    const auto S0 = fp.surface(st.spans.front());
    Vec2 Q = S0.xr(mprime_at_x(S0, in.z_ref));
    st.x.push_back(in.z_ref);
    for (std::size_t k = 1; k < st.spans.size(); ++k) {
        const auto S = fp.surface(st.spans[k]);
        if (st.variant == Stacking::Variant::SameZTheta) {
            (void)mprime_at_x(S, in.z_ref);  // the radial line must cross every layer
            st.x.push_back(in.z_ref);
        } else {
            Q = S.xr(closest_mprime(S, Q));
            st.x.push_back(Q.x());
        }
    }
    return st;
}

SectionCG section_cg(const Blade::StreamSurface& S, const std::vector<Vec2>& conformal)
{
    std::vector<Vec2> dev;
    dev.reserve(conformal.size());
    for (const auto& q : conformal) dev.emplace_back(S.m_of_mprime(q.x()), S.xr(q.x()).y() * q.y());
    const auto pc = polygon_centroid(dev);
    // m' of the centroid's m (m(m') increases; beyond the ends it continues straight)
    const auto f = [&](double mp) { return S.m_of_mprime(mp) - pc.centroid.x(); };
    double lo = S.mprime_inlet(), hi = S.mprime_outlet();
    const double w = hi - lo;
    for (int i = 0; i < 20 && f(lo) > 0.0; ++i) lo -= w;
    for (int i = 0; i < 20 && f(hi) < 0.0; ++i) hi += w;
    SectionCG r;
    r.area = pc.area;
    r.mprime = Math::brent_root(f, lo, hi, 1e-14);
    r.theta = pc.centroid.y() / S.xr(r.mprime).y();
    return r;
}

Placement nominal_placement(const StackingLine& st, const Blade::StreamSurface& S, double stagger_deg, double theta_le_deg,
                            double chord_scale, const std::vector<Vec2>& loop, bool mirror)
{
    if (st.mode == Stacking::Mode::Edges) return placement(stagger_deg, theta_le_deg, S.mprime_te(), mirror);
    if (!(chord_scale > 0.0)) throw std::invalid_argument("nominal_placement: chord scale must be > 0");
    // c' = chord_scale m'_TE / cos(gamma); the reference point q_ref onto the target T
    Placement p = placement(stagger_deg, 0.0, chord_scale * S.mprime_te(), mirror);
    const Vec2 T(mprime_at_x(S, st.x_at(S.span())), st.theta);
    const Vec2 q = st.reference == Stacking::Reference::LE   ? Vec2(0.0, 0.0)
                   : st.reference == Stacking::Reference::TE ? Vec2(1.0, 0.0)
                                                             : polygon_centroid(loop).centroid;
    const Vec2 at = p(q);
    p.m0 = T.x() - at.x();
    p.t0 = T.y() - at.y();
    if (st.reference != Stacking::Reference::CG) return p;
    // the physical centre of gravity onto T: shift in (m', theta) until it is there
    std::vector<Vec2> conf(loop.size());
    for (int it = 0; it < 60; ++it) {
        for (std::size_t i = 0; i < loop.size(); ++i) conf[i] = p(loop[i]);
        const auto c = section_cg(S, conf);
        const Vec2 dlt = T - Vec2(c.mprime, c.theta);
        p.m0 += dlt.x();
        p.t0 += dlt.y();
        if (std::max(std::abs(dlt.x()), std::abs(dlt.y())) < 1e-13) return p;
    }
    throw std::invalid_argument("the centre-of-gravity stacking did not converge");
}

LayerResult place_layer(const Layer& layer, std::shared_ptr<const Section2D> sec, const Blade::FlowPath& fp,
                        int n_blades, bool mirror, double restagger_deg, const Vec2& axis, const StackingLine* stacking)
{
    LayerResult r;
    r.section = sec;
    r.surface = fp.surface(layer.span);
    const auto& S = r.surface;
    const StackingLine edges{};
    const StackingLine& st = stacking ? *stacking : edges;
    const Vec2 pivot = restagger_deg == 0.0 ? Vec2::Zero() : Vec2(mprime_at_x(S, axis.x()), axis.y());
    const auto conf = restaggered(
        nominal_placement(st, S, layer.stagger_deg, layer.theta_le_deg, layer.chord_scale, sec->loop.xy, mirror), restagger_deg, pivot);
    r.gamma = conf.gamma;
    r.chord_mprime = conf.c;
    r.mprime_le = conf.m0;
    r.theta_le = conf.t0;
    const Vec2 te = conf({1.0, 0.0});
    r.conformal.reserve(sec->loop.xy.size());
    r.xyz.reserve(sec->loop.xy.size());
    for (const auto& q : sec->loop.xy) {
        r.conformal.push_back(conf(q));
        r.xyz.push_back(S.cartesian(r.conformal.back().x(), r.conformal.back().y()));
    }
    const auto& cam = sec->camber;
    const std::size_t n = cam.point.size();
    const double m_le = S.m_of_mprime(conf.m0);
    for (std::size_t i = 0; i < n; ++i) {
        const Vec2 mt = conf(cam.point[i]);
        r.mprime.push_back(mt.x() - conf.m0);
        r.m.push_back(S.m_of_mprime(mt.x()) - m_le);
        r.mfrac.push_back((mt.x() - conf.m0) / (te.x() - conf.m0));
        r.theta_deg.push_back(Math::rad2deg(mt.y()));
        // the in-plane map is a similarity: a direction at angle a becomes gamma + sign * a
        r.beta_deg.push_back(Math::rad2deg(conf.gamma) + conf.sign * cam.angle_deg[i]);
        r.thickness.push_back(cam.thickness[i]);
        r.xc.push_back(cam.point[i].x());
        r.camber_xr.push_back(S.xr(mt.x()));
    }
    r.beta_le = r.beta_deg.front();
    r.beta_te = r.beta_deg.back();
    r.theta_te_deg = r.theta_deg.back();
    r.r_le = S.xr(conf.m0).y();
    r.r_te = S.xr(te.x()).y();
    const Vec3 pa = S.cartesian(conf.m0, conf.t0), pb = S.cartesian(te.x(), te.y());
    r.chord = (pb - pa).norm();
    r.pitch_le = 2.0 * Math::kPi * r.r_le / n_blades;
    r.solidity_3d = r.chord / (Math::kPi * (r.r_le + r.r_te) / n_blades);
    // centre of gravity: area centroid in the (m, r theta) plane, mapped back to the surface
    const auto cg = section_cg(S, r.conformal);
    r.section_area = cg.area;
    r.cg_mprime = cg.mprime;
    r.cg_theta = cg.theta;
    r.cg_xr = S.xr(r.cg_mprime);
    r.cg_xyz = S.cartesian(r.cg_mprime, r.cg_theta);
    // stacking point
    if (st.mode == Stacking::Mode::Edges || st.reference == Stacking::Reference::LE) r.stack_point = {conf.m0, conf.t0};
    else if (st.reference == Stacking::Reference::TE) r.stack_point = te;
    else r.stack_point = {r.cg_mprime, r.cg_theta};
    r.stack_xr = S.xr(r.stack_point.x());
    return r;
}

namespace {

/// The row at span s between the layers: stagger, lean, chord scale and (when `loops` is given)
/// the unstaggered profile, by natural spline weights over the layer spans (exact at the layers).
struct RowSpec
{
    double stagger_deg = 0.0, theta_le_deg = 0.0, chord_scale = 0.0;
    std::vector<Vec2> loop;
};

RowSpec row_spec(const BladeDesign& d, const std::vector<const std::vector<Vec2>*>& loops, double s)
{
    std::vector<double> spans;
    for (const auto& l : d.layers) spans.push_back(l.span);
    const auto w = spline_weights(spans, s);
    RowSpec r;
    for (std::size_t k = 0; k < d.layers.size(); ++k) {
        r.stagger_deg += w[k] * d.layers[k].stagger_deg;
        r.theta_le_deg += w[k] * d.layers[k].theta_le_deg;
        r.chord_scale += w[k] * d.layers[k].chord_scale;
    }
    if (loops.size() == d.layers.size()) {
        r.loop.assign(loops.front()->size(), Vec2::Zero());
        for (std::size_t k = 0; k < loops.size(); ++k)
            for (std::size_t i = 0; i < r.loop.size(); ++i) r.loop[i] += w[k] * (*loops[k])[i];
    }
    return r;
}

Vec2 restagger_axis_impl(const BladeDesign& d, const Blade::FlowPath& fp, const StackingLine& st,
                         const std::vector<const std::vector<Vec2>*>& loops)
{
    const auto& R = d.restagger;
    if (R.axis == Restagger::Axis::Point) return {R.z, Math::deg2rad(R.theta_deg)};
    const double s0 = d.layers.front().span, s1 = d.layers.back().span;
    const double s = R.axis == Restagger::Axis::Hub ? s0 : R.axis == Restagger::Axis::Tip ? s1 : 0.5 * (s0 + s1);
    // the nominal (un-re-staggered) reference section, interpolated in span like the mesh rows
    const auto row = row_spec(d, loops, s);
    const auto S = fp.surface(s);
    const Vec2 P = nominal_placement(st, S, row.stagger_deg, row.theta_le_deg, row.chord_scale, row.loop, d.mirror)({R.chord_fraction, 0.0});
    return {S.xr(P.x()).x(), P.y()};
}

}  // namespace

Vec2 restagger_axis(const BladeDesign& d, const Blade::FlowPath& fp)
{
    const auto st = resolve_stacking(d, fp);
    std::vector<std::shared_ptr<const Section2D>> secs;
    std::vector<const std::vector<Vec2>*> loops;
    if (st.mode == Stacking::Mode::ReferencePoint && st.reference == Stacking::Reference::CG) {
        Evaluator ev;
        for (const auto& l : d.layers) {
            secs.push_back(ev.section(l, d.mesh));
            loops.push_back(&secs.back()->loop.xy);
        }
    }
    return restagger_axis_impl(d, fp, st, loops);
}

// ── design ───────────────────────────────────────────────────────────────────

bool DesignResult::ok() const
{
    if (!errors.empty()) return false;
    return std::all_of(layers.begin(), layers.end(), [](const LayerResult& l) { return l.error.empty(); });
}

DesignResult Evaluator::evaluate(const BladeDesign& d)
{
    const auto t0 = std::chrono::steady_clock::now();
    DesignResult out;
    std::optional<Blade::FlowPath> fp;
    try {
        d.validate();
        fp.emplace(d.meridional);
    } catch (const std::exception& e) {
        out.errors.push_back(e.what());
        return out;
    }
    // 2-D sections: cached ones directly, the others in parallel
    const std::size_t nl = d.layers.size();
    std::vector<std::shared_ptr<const Section2D>> secs(nl);
    std::vector<std::string> errs(nl);
    std::vector<std::future<void>> jobs;
    for (std::size_t k = 0; k < nl; ++k)
        jobs.push_back(std::async(std::launch::async, [&, k] {
            try {
                secs[k] = section(d.layers[k], d.mesh);
            } catch (const std::exception& e) {
                errs[k] = e.what();
            }
        }));
    for (auto& j : jobs) j.get();
    // stacking
    StackingLine st;
    try {
        st = resolve_stacking(d, *fp);
    } catch (const std::exception& e) {
        out.errors.push_back(std::string("stacking: ") + e.what());
        return out;
    }
    std::vector<const std::vector<Vec2>*> loops;
    if (std::all_of(secs.begin(), secs.end(), [](const auto& p) { return p != nullptr; }))
        for (const auto& p : secs) loops.push_back(&p->loop.xy);
    // re-stagger axis
    Vec2 axis = Vec2::Zero();
    const double dg = d.restagger.angle_deg;
    if (dg != 0.0) {
        try {
            axis = restagger_axis_impl(d, *fp, st, loops);
            // the axis must cross every stream surface between hub and casing
            for (const double s : {0.0, 1.0}) (void)mprime_at_x(fp->surface(s), axis.x());
        } catch (const std::exception& e) {
            out.errors.push_back(std::string("re-stagger: ") + e.what());
            return out;
        }
    }
    // placement on the stream surfaces
    out.layers.resize(nl);
    for (std::size_t k = 0; k < nl; ++k) {
        if (!errs[k].empty()) {
            out.layers[k].error = errs[k];
            continue;
        }
        try {
            out.layers[k] = place_layer(d.layers[k], secs[k], *fp, d.n_blades, d.mirror, dg, axis, &st);
        } catch (const std::exception& e) {
            out.layers[k].error = e.what();
        }
    }
    // meridional view
    const auto hub = fp->surface(0.0), cas = fp->surface(1.0);
    out.hub = hub.polyline();
    out.casing = cas.polyline();
    for (const auto& l : out.layers)
        out.streamlines.push_back(l.error.empty() ? l.surface.polyline() : std::vector<Vec2>{});
    out.le_line = {fp->le_point(0.0), fp->le_point(1.0)};
    out.te_line = {fp->te_point(0.0), fp->te_point(1.0)};

    if (!out.ok()) {
        out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return out;
    }
    // 3-D mesh: rows between the layers
    std::vector<double> spans;
    for (const auto& l : d.layers) spans.push_back(l.span);
    const auto rows = tanh_two_sided(d.mesh.n_rows, d.mesh.row_clustering);
    const double s_lo = spans.front(), s_hi = spans.back();
    const std::size_t np = out.layers.front().section->loop.xy.size();
    for (const auto& l : out.layers)
        if (l.section->loop.xy.size() != np) {
            out.errors.push_back("layers have different loop layouts");
            return out;
        }
    out.blade.ni = static_cast<int>(np);
    out.blade.nj = d.mesh.n_rows;
    out.blade.p.reserve(np * rows.size());
    const int nf = d.mesh.n_field_chord;
    out.beta.ni = out.thickness.ni = nf;
    out.beta.nj = out.thickness.nj = d.mesh.n_rows;
    // each layer's beta / thickness resampled on the meridional chord fraction grid
    std::vector<std::vector<double>> beta_k(nl), thk_k(nl);
    for (std::size_t k = 0; k < nl; ++k) {
        const auto& L = out.layers[k];
        for (int i = 0; i < nf; ++i) {
            const double f = static_cast<double>(i) / (nf - 1);
            beta_k[k].push_back(lerp_table(L.mfrac, L.beta_deg, f));
            thk_k[k].push_back(lerp_table(L.mfrac, L.thickness, f));
        }
    }
    double theta_hub = 0.0, theta_cas = 0.0;
    std::vector<Vec2> row_conf(np);
    for (std::size_t j = 0; j < rows.size(); ++j) {
        const double s = s_lo + rows[j] * (s_hi - s_lo);
        out.row_span.push_back(s);
        const auto w = spline_weights(spans, s);
        const auto row = row_spec(d, loops, s);
        const auto S = fp->surface(s);
        Placement conf;
        try {
            const Vec2 pivot = dg == 0.0 ? Vec2::Zero() : Vec2(mprime_at_x(S, axis.x()), axis.y());
            conf = restaggered(nominal_placement(st, S, row.stagger_deg, row.theta_le_deg, row.chord_scale, row.loop, d.mirror), dg, pivot);
        } catch (const std::exception& e) {
            out.errors.push_back(std::string("stacking: mesh row at span ") + std::to_string(s) + ": " + e.what());
            out.blade = {};
            return out;
        }
        const double mp_le = conf.m0, mp_te = conf({1.0, 0.0}).x();
        double th_mean = 0.0;
        for (std::size_t i = 0; i < np; ++i) {
            row_conf[i] = conf(row.loop[i]);
            out.blade.p.push_back(S.cartesian(row_conf[i].x(), row_conf[i].y()));
            th_mean += row_conf[i].y();
        }
        th_mean /= static_cast<double>(np);
        // the row's LE / TE and stacking point
        out.blade_le.push_back(S.xr(mp_le));
        out.blade_te.push_back(S.xr(mp_te));
        Vec2 sp(conf.m0, conf.t0);
        if (st.mode == Stacking::Mode::ReferencePoint && st.reference == Stacking::Reference::TE) sp = conf({1.0, 0.0});
        if (st.mode == Stacking::Mode::ReferencePoint && st.reference == Stacking::Reference::CG) {
            const auto c = section_cg(S, row_conf);
            sp = {c.mprime, c.theta};
        }
        out.stacking_xr.push_back(S.xr(sp.x()));
        out.stacking_xyz.push_back(S.cartesian(sp.x(), sp.y()));
        if (j == 0) theta_hub = th_mean;
        if (j + 1 == rows.size()) theta_cas = th_mean;
        for (int i = 0; i < nf; ++i) {
            const double f = static_cast<double>(i) / (nf - 1);
            double b = 0.0, t = 0.0;
            for (std::size_t k = 0; k < nl; ++k) {
                b += w[k] * beta_k[k][static_cast<std::size_t>(i)];
                t += w[k] * thk_k[k][static_cast<std::size_t>(i)];
            }
            const Vec2 xr = S.xr(mp_le + f * (mp_te - mp_le));
            out.beta.xr.push_back(xr);
            out.beta.value.push_back(b);
            out.thickness.xr.push_back(xr);
            out.thickness.value.push_back(t);
        }
    }
    for (auto* f : {&out.beta, &out.thickness}) {
        const auto [lo, hi] = std::minmax_element(f->value.begin(), f->value.end());
        f->lo = *lo;
        f->hi = *hi;
    }
    // hub and casing walls over one pitch around the blade
    const double pitch = 2.0 * Math::kPi / d.n_blades;
    auto wall = [&](const Blade::StreamSurface& S, double th_c, Grid& G) {
        G.ni = d.mesh.n_wall_m;
        G.nj = d.mesh.n_wall_theta;
        G.p.clear();
        for (int j = 0; j < G.nj; ++j) {
            const double th = th_c + pitch * (static_cast<double>(j) / (G.nj - 1) - 0.5);
            for (int i = 0; i < G.ni; ++i) {
                const double mp = S.mprime_inlet() + (S.mprime_outlet() - S.mprime_inlet()) * i / (G.ni - 1);
                G.p.push_back(S.cartesian(mp, th));
            }
        }
    };
    wall(hub, theta_hub, out.hub_wall);
    wall(cas, theta_cas, out.casing_wall);
    // the blade's volume and centre of gravity
    const auto mp = mass_properties(out.blade);
    out.blade_volume = mp.volume;
    out.blade_cg = mp.centroid;
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

void fit_layer(Layer& layer, Evaluator& ev)
{
    const auto sec = ev.naca(layer.naca);
    const auto r = Fit::fit_naca65(*sec);
    layer.custom = r.params;
    layer.fitted_to = layer.naca;
    layer.fit_rms = r.rms;
    layer.method = Method::Naca65Custom;
}

std::vector<std::vector<Vec3>> section_loops(const DesignResult& r)
{
    if (r.blade.p.empty()) throw std::invalid_argument("section_loops: no blade mesh");
    std::vector<std::vector<Vec3>> loops;
    for (int j = 0; j < r.blade.nj; ++j) {
        std::vector<Vec3> l;
        for (int i = 0; i < r.blade.ni; ++i) l.push_back(r.blade.at(i, j));
        loops.push_back(std::move(l));
    }
    return loops;
}

}  // namespace BladeGen::Design
