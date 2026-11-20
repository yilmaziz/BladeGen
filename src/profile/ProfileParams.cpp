#include "profile/ProfileParams.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <unordered_map>

namespace BladeGen::Profile {

std::string to_string(ConnectorMode m)
{
    switch (m) {
        case ConnectorMode::Fair: return "Fair";
        case ConnectorMode::C2Parametric: return "C2Parametric";
        case ConnectorMode::G2Legacy: return "G2Legacy";
    }
    return "Fair";
}

ConnectorMode connector_mode_from_string(const std::string& s)
{
    if (s == "Fair" || s == "fair") return ConnectorMode::Fair;
    if (s == "C2Parametric" || s == "C2" || s == "c2") return ConnectorMode::C2Parametric;
    if (s == "G2Legacy" || s == "g2" || s == "G2" || s == "c2_quintic" || s == "clamped_quintic")
        return ConnectorMode::G2Legacy;
    throw std::invalid_argument("unknown connector mode '" + s + "' (Fair, C2Parametric or G2Legacy)");
}

ProfileParams ProfileParams::naca65_1210()
{
    // BladeGen's own inverse fit of NACA 65-(12)10 (blade2d naca65_1210_XY_dense241.csv)
    // with Fair connectors and fairness weight 1: RMS 0.0104 % chord, no geometry
    // errors (the blade2d fit below has an S-shaped TE connector).
    ProfileParams p;
    p.le = {0.012627, 1.257999, 0.004407, 13.058742, 24.067893, 11.360914, 1.0};
    p.te = {0.0002, 4.0, -0.000125, 26.350313, 11.924809, 0.0, 1.0};
    p.ss = {0.094282, 0.52934, 0.178684, 0.026252};
    p.thickness = {0.110394, 0.40528, 0.036028, 0.001928, 0.683756, 1.147683, false};
    p.ps = {0.010024, 0.048277, {0.2, 0.45, 0.8}};
    return p;
}

ProfileParams ProfileParams::naca65_1210_blade2d()
{
    // blade2d naca65/naca65_1210_blade_params.json (inverse fit to NACA 65-(12)10)
    // with the blade2d connector: the geometry of the Python engine.
    ProfileParams p;
    p.le = {0.012073, 1.198788, 0.003422, 15.561637, 22.417183, 15.638344, 0.28};
    p.te = {0.000255, 4.0, -0.000526, 24.659155, 6.386629, 0.0, 0.20};
    p.ss = {0.092724, 0.52168, 0.177553, 0.041191};
    p.thickness = {0.11027, 0.405561, 0.035922, 0.002126, 0.685326, 1.145254, false};
    p.ps = {0.006057, 0.05123, {0.2, 0.45, 0.8}};
    p.connector = ConnectorMode::G2Legacy;  // le_alpha_* = 0.28, te_alpha_* = 0.20
    return p;
}

ProfileParams ProfileParams::cda_start()
{
    // blade2d CDA/fit_cda_c2.py SPEC, initial values x0.
    ProfileParams p;
    p.le = {0.010, 2.0, 0.0, 17.0, 16.0, 16.0, 1.0};
    p.te = {0.004, 1.2, 0.0, 9.0, 1.0, 1.0, 1.0};
    p.ss = {0.16, 0.40, 0.15, 0.29};
    p.thickness = {0.080, 0.30, 0.030, 0.006, 1.0, 1.5, false};
    p.ps = {0.14, 0.18, {0.20, 0.45, 0.80}};
    return p;
}

namespace {

using Getter = std::function<double&(ProfileParams&)>;

struct Entry
{
    ParamInfo info;
    Getter ref;
};

const std::vector<Entry>& entries()
{
    static const std::vector<Entry> e = {
        // family, label, suggested bounds (blade2d fit_cda_c2.py SPEC, widened where useful)
        {{"le.radius", "Leading edge", "LE radius r", 0.001, 0.06}, [](ProfileParams& p) -> double& { return p.le.radius; }},
        {{"le.stretch_ab", "Leading edge", "LE stretch a/b", 0.5, 6.0}, [](ProfileParams& p) -> double& { return p.le.stretch_ab; }},
        {{"le.asym", "Leading edge", "LE asymmetry", -0.02, 0.02}, [](ProfileParams& p) -> double& { return p.le.asym; }},
        {{"le.beta_deg", "Leading edge", "LE axis angle beta [deg]", -10.0, 45.0}, [](ProfileParams& p) -> double& { return p.le.beta_deg; }},
        {{"le.d_beta_ss_deg", "Leading edge", "LE wedge SS [deg]", 0.0, 40.0}, [](ProfileParams& p) -> double& { return p.le.d_beta_ss_deg; }},
        {{"le.d_beta_ps_deg", "Leading edge", "LE wedge PS [deg]", 0.0, 40.0}, [](ProfileParams& p) -> double& { return p.le.d_beta_ps_deg; }},
        {{"le.tension", "Leading edge", "LE connector tension", 0.2, 3.0}, [](ProfileParams& p) -> double& { return p.le.tension; }},
        {{"te.radius", "Trailing edge", "TE radius r", 0.0002, 0.03}, [](ProfileParams& p) -> double& { return p.te.radius; }},
        {{"te.stretch_ab", "Trailing edge", "TE stretch a/b", 0.2, 6.0}, [](ProfileParams& p) -> double& { return p.te.stretch_ab; }},
        {{"te.asym", "Trailing edge", "TE asymmetry", -0.01, 0.01}, [](ProfileParams& p) -> double& { return p.te.asym; }},
        {{"te.beta_deg", "Trailing edge", "TE axis angle beta [deg]", -30.0, 40.0}, [](ProfileParams& p) -> double& { return p.te.beta_deg; }},
        {{"te.d_beta_ss_deg", "Trailing edge", "TE wedge SS [deg]", 0.0, 30.0}, [](ProfileParams& p) -> double& { return p.te.d_beta_ss_deg; }},
        {{"te.d_beta_ps_deg", "Trailing edge", "TE wedge PS [deg]", 0.0, 30.0}, [](ProfileParams& p) -> double& { return p.te.d_beta_ps_deg; }},
        {{"te.tension", "Trailing edge", "TE connector tension", 0.2, 3.0}, [](ProfileParams& p) -> double& { return p.te.tension; }},
        {{"ss.handle_le_frac", "Suction side", "SS LE handle", 0.005, 0.6}, [](ProfileParams& p) -> double& { return p.ss.handle_le_frac; }},
        {{"ss.ctrl2_mprime", "Suction side", "SS CP2 m'", 0.1, 0.8}, [](ProfileParams& p) -> double& { return p.ss.ctrl2_mprime; }},
        {{"ss.ctrl2_theta", "Suction side", "SS CP2 theta", -0.1, 0.4}, [](ProfileParams& p) -> double& { return p.ss.ctrl2_theta; }},
        {{"ss.handle_te_frac", "Suction side", "SS TE handle", 0.005, 0.6}, [](ProfileParams& p) -> double& { return p.ss.handle_te_frac; }},
        {{"thickness.t_max", "Thickness distribution", "max thickness t_max", 0.005, 0.2}, [](ProfileParams& p) -> double& { return p.thickness.t_max; }},
        {{"thickness.mprime_max", "Thickness distribution", "t_max position m'", 0.1, 0.8}, [](ProfileParams& p) -> double& { return p.thickness.mprime_max; }},
        {{"thickness.t_le", "Thickness distribution", "thickness at LE side t_le", 0.0, 0.1}, [](ProfileParams& p) -> double& { return p.thickness.t_le; }},
        {{"thickness.t_te", "Thickness distribution", "thickness at TE side t_te", 0.0, 0.04}, [](ProfileParams& p) -> double& { return p.thickness.t_te; }},
        {{"thickness.shape_le", "Thickness distribution", "forward exponent", 0.2, 4.0}, [](ProfileParams& p) -> double& { return p.thickness.shape_le; }},
        {{"thickness.shape_te", "Thickness distribution", "aft exponent", 0.2, 5.0}, [](ProfileParams& p) -> double& { return p.thickness.shape_te; }},
        {{"ps.handle_le_frac", "Pressure side", "PS LE handle", 0.001, 0.6}, [](ProfileParams& p) -> double& { return p.ps.handle_le_frac; }},
        {{"ps.handle_te_frac", "Pressure side", "PS TE handle", 0.001, 0.6}, [](ProfileParams& p) -> double& { return p.ps.handle_te_frac; }},
    };
    return e;
}

const Entry& find(const std::string& name)
{
    static const auto index = [] {
        std::unordered_map<std::string, std::size_t> m;
        const auto& e = entries();
        for (std::size_t i = 0; i < e.size(); ++i) m.emplace(e[i].info.name, i);
        return m;
    }();
    const auto it = index.find(name);
    if (it == index.end()) throw std::invalid_argument("unknown profile parameter '" + name + "'");
    return entries()[it->second];
}

}  // namespace

const std::vector<ParamInfo>& param_table()
{
    static const std::vector<ParamInfo> t = [] {
        std::vector<ParamInfo> v;
        for (const auto& e : entries()) v.push_back(e.info);
        return v;
    }();
    return t;
}

double get_param(const ProfileParams& p, const std::string& name)
{
    auto copy = p;
    return find(name).ref(copy);
}

void set_param(ProfileParams& p, const std::string& name, double value)
{
    if (!std::isfinite(value)) throw std::invalid_argument("set_param: non-finite value for " + name);
    find(name).ref(p) = value;
}

std::vector<double> to_vector(const ProfileParams& p)
{
    auto copy = p;
    std::vector<double> v;
    for (const auto& e : entries()) v.push_back(e.ref(copy));
    return v;
}

void from_vector(ProfileParams& p, const std::vector<double>& v)
{
    const auto& e = entries();
    if (v.size() != e.size()) throw std::invalid_argument("from_vector: wrong parameter count");
    for (std::size_t i = 0; i < e.size(); ++i) e[i].ref(p) = v[i];
}

ProfileParams blend(const ProfileParams& a, const ProfileParams& b, double w)
{
    ProfileParams out = a;
    const auto va = to_vector(a), vb = to_vector(b);
    std::vector<double> v(va.size());
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = va[i] + w * (vb[i] - va[i]);
    from_vector(out, v);
    out.thickness.enforce_c1 = (w < 0.5) ? a.thickness.enforce_c1 : b.thickness.enforce_c1;
    if (a.ps.mprime_ctrl.size() == b.ps.mprime_ctrl.size())
        for (std::size_t i = 0; i < out.ps.mprime_ctrl.size(); ++i)
            out.ps.mprime_ctrl[i] = a.ps.mprime_ctrl[i] + w * (b.ps.mprime_ctrl[i] - a.ps.mprime_ctrl[i]);
    return out;
}

}  // namespace BladeGen::Profile
