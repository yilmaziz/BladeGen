// tests/turbolenza/RefData.hpp — reading the TURBOlenza v1 reference data.
//
// The JSON files in tests/turbolenza/data are written by tools/turbolenza_reference/dump_v1.py
// from the v1 Python solver (they are not committed: run the script to make them).  A test that
// needs one calls ref::load(name) and SKIPs when it is missing.  Arrays are {"shape", "data"},
// flattened in C order (row major), which is the layout of turbolenza::Field.
#pragma once

#include "turbolenza/Curve.hpp"
#include "turbolenza/Field.hpp"
#include "turbolenza/GridData.hpp"

#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace ref {

using nlohmann::json;
namespace T = BladeGen::Turbolenza;

inline std::string data_dir() { return std::string(BLADEGEN_SOURCE_DIR) + "/tests/turbolenza/data/"; }

/// The reference file `name` (without .json), or nullopt when it has not been generated.
inline std::optional<json> load(const std::string& name)
{
    std::ifstream in(data_dir() + name + ".json");
    if (!in) return std::nullopt;
    return json::parse(in);
}

/// Loads the file or SKIPs the test (the reference data is generated, not committed).
#define TURBOLENZA_REF_OR_SKIP(var, name)                                                              \
    auto var##_opt = ref::load(name);                                                                  \
    if (!var##_opt) SKIP("reference data " name ".json missing: run tools/turbolenza_reference/dump_v1.py"); \
    const ref::json& var = *var##_opt

inline double num(const json& j)
{
    if (j.is_string()) {
        const auto s = j.get<std::string>();
        if (s == "nan") return std::numeric_limits<double>::quiet_NaN();
        if (s == "inf") return std::numeric_limits<double>::infinity();
        if (s == "-inf") return -std::numeric_limits<double>::infinity();
    }
    if (j.is_null()) return std::numeric_limits<double>::quiet_NaN();
    return j.get<double>();
}

inline std::vector<int> shape(const json& a) { return a.at("shape").get<std::vector<int>>(); }

inline std::vector<double> vec(const json& a)
{
    std::vector<double> out;
    const json& d = a.at("data");
    out.reserve(d.size());
    for (const auto& x : d) out.push_back(num(x));
    return out;
}

inline T::Field field(const json& a)
{
    const auto s = shape(a);
    REQUIRE(s.size() == 2);
    T::Field f(s[0], s[1]);
    f.v = vec(a);
    return f;
}

/// (n, 2) array -> points.
inline std::vector<T::Point> points(const json& a)
{
    const auto s = shape(a);
    REQUIRE((s.size() == 2 && s[1] == 2));
    const auto d = vec(a);
    std::vector<T::Point> p(static_cast<std::size_t>(s[0]));
    for (std::size_t k = 0; k < p.size(); ++k) p[k] = {d[2 * k], d[2 * k + 1]};
    return p;
}

inline T::WallSide wall_side(const json& w)
{
    T::WallSide s;
    for (const auto& name : w.at("surface").at("strings"))
        s.surface.push_back(name.get<std::string>() == "suction" ? T::Surface::Suction : T::Surface::Pressure);
    s.parameter = vec(w.at("parameter"));
    s.u = vec(w.at("u"));
    s.contour = vec(w.at("contour"));
    for (double h : vec(w.at("hidden"))) s.hidden.push_back(static_cast<int>(h));
    return s;
}

/// A v1 Grid (grid_json of dump_v1.py) as turbolenza::Grid.  blade_curves stays empty: a loaded
/// grid carries the wall-curve derivatives (wall_d) but not the curves.
inline T::Grid grid(const json& g)
{
    T::Grid out;
    const json& b = g.at("boundaries");
    auto& B = out.boundaries;
    B.chord = num(b.at("chord"));
    B.pitch = num(b.at("pitch"));
    const auto t = vec(b.at("periodic_translation"));
    B.periodic_translation = {t[0], t[1]};
    B.inlet_angle = num(b.at("inlet_angle"));
    B.exit_angle = num(b.at("exit_angle"));
    B.m_in = num(b.at("m_in"));
    B.m_le = num(b.at("m_le"));
    B.m_te = num(b.at("m_te"));
    B.m_out = num(b.at("m_out"));
    B.lower_blade_boundary = points(b.at("lower_blade_boundary"));
    B.upper_blade_boundary = points(b.at("upper_blade_boundary"));
    B.upstream_m = vec(b.at("upstream_m"));
    B.upstream_beta = vec(b.at("upstream_beta"));
    B.upstream_lower_boundary = points(b.at("upstream_lower_boundary"));
    B.upstream_upper_boundary = points(b.at("upstream_upper_boundary"));
    B.downstream_m = vec(b.at("downstream_m"));
    B.downstream_beta = vec(b.at("downstream_beta"));
    B.downstream_lower_boundary = points(b.at("downstream_lower_boundary"));
    B.downstream_upper_boundary = points(b.at("downstream_upper_boundary"));
    B.kappa_le = B.upstream_beta.back();
    B.kappa_te = B.downstream_beta.front();
    out.m = vec(g.at("m"));
    out.eta = vec(g.at("eta"));
    const auto ns = shape(g.at("nodes"));
    REQUIRE(ns.size() == 3);
    const auto nd = vec(g.at("nodes"));
    out.node_m = T::Field(ns[0], ns[1]);
    out.node_y = T::Field(ns[0], ns[1]);
    for (std::size_t k = 0; k < out.node_m.size(); ++k) {
        out.node_m.v[k] = nd[2 * k];
        out.node_y.v[k] = nd[2 * k + 1];
    }
    out.first_blade_i = g.at("first_blade_i").get<int>();
    out.last_blade_i = g.at("last_blade_i").get<int>();
    const json& wi = g.at("wall_identity");
    for (double s : vec(wi.at("stations"))) out.wall_identity.stations.push_back(static_cast<int>(s));
    out.wall_identity.lower = wall_side(wi.at("lower"));
    out.wall_identity.upper = wall_side(wi.at("upper"));
    if (g.contains("wall_curve_derivative")) {
        for (int k = 0; k < 3; ++k) {
            const json& d = g.at("wall_curve_derivative").at(std::to_string(k));
            out.wall_d[k].lower = points(d.at("lower"));
            out.wall_d[k].upper = points(d.at("upper"));
        }
    }
    return out;
}

/// A v1 ParametricBladeCurve (knots, control_points, reference_length, degree).
inline T::ParametricBladeCurve curve(const json& c)
{
    return {vec(c.at("knots")), points(c.at("control_points")), num(c.at("reference_length")),
            c.at("degree").get<int>()};
}

/// v1 SurfaceCurves.
inline T::SurfaceCurves surface_curves(const json& s)
{
    return {curve(s.at("suction")), curve(s.at("pressure")), curve(s.at("suction_unstaggered")),
            curve(s.at("pressure_unstaggered")), vec(s.at("suction_parameters")), vec(s.at("pressure_parameters"))};
}

/// max |a - b| / max(|b|) over the entries (NaN in both at the same place counts as equal;
/// NaN in one only gives +inf).
inline double rel_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    REQUIRE(a.size() == b.size());
    double num_ = 0.0, den = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        if (std::isnan(a[k]) || std::isnan(b[k])) {
            if (std::isnan(a[k]) != std::isnan(b[k])) return std::numeric_limits<double>::infinity();
            continue;
        }
        num_ = std::max(num_, std::abs(a[k] - b[k]));
        den = std::max(den, std::abs(b[k]));
    }
    return den > 0.0 ? num_ / den : num_;
}
inline double rel_diff(const T::Field& a, const T::Field& b)
{
    REQUIRE(a.same_shape(b));
    return rel_diff(a.v, b.v);
}

}  // namespace ref
