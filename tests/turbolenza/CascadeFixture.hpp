// tests/turbolenza/CascadeFixture.hpp — the v1 reference cascade (NACA 65-(12A10)10, C2, solidity 1,
// stagger 33 deg) rebuilt by the port: the BladeSection from v1's blade-frame curves, then
// make_cascade.  Shared by test_cascade.cpp and test_grid*.cpp.
#pragma once

#include "RefData.hpp"

#include "turbolenza/Cascade.hpp"

#include <cmath>
#include <optional>
#include <string>

namespace fixture {

namespace T = BladeGen::Turbolenza;

inline constexpr double kPi = 3.141592653589793;

/// A quintic Bezier curve (no interior knot) through the given control points.
inline T::ParametricBladeCurve bezier(const std::vector<T::Point>& P)
{
    const int p = static_cast<int>(P.size()) - 1;
    std::vector<double> t(static_cast<std::size_t>(2 * p + 2), 0.0);
    for (int k = p + 1; k < 2 * p + 2; ++k) t[static_cast<std::size_t>(k)] = 1.0;
    return {t, P, 1.0, p};
}


/// A thin cambered synthetic blade (two quintic Bezier curves): no v1 data needed.
inline T::BladeSection synthetic_section()
{
    T::BladeSection s;
    s.designation = "synthetic quintic blade";
    s.suction = bezier({{0.0, 0.0}, {0.0, 0.02}, {0.2, 0.09}, {0.6, 0.10}, {0.9, 0.03}, {1.0, 0.0}});
    s.pressure = bezier({{0.0, 0.0}, {0.0, -0.01}, {0.2, 0.0}, {0.6, 0.04}, {0.9, 0.01}, {1.0, 0.0}});
    for (int k = 0; k <= 80; ++k) {
        const double u = 0.5 * (1.0 - std::cos(kPi * k / 80.0));
        s.suction_parameters.push_back(u);
        s.pressure_parameters.push_back(u);
    }
    s.suction_parameters.front() = s.pressure_parameters.front() = 0.0;
    s.suction_parameters.back() = s.pressure_parameters.back() = 1.0;
    s.chi_le = 0.3;
    s.chi_te = 0.25;
    return s;
}


/// The BladeSection of the v1 cascade (blade frame), from grid_family.json's cascade.
inline T::BladeSection v1_section(const ref::json& cascade)
{
    const ref::json& c = cascade.at("active_curves");
    const ref::json& a = cascade.at("active");
    const double stagger = ref::num(cascade.at("stagger"));
    T::BladeSection s;
    s.designation = cascade.at("metadata").at("designation").get<std::string>();
    s.suction = ref::curve(c.at("suction_unstaggered"));
    s.pressure = ref::curve(c.at("pressure_unstaggered"));
    s.suction_parameters = ref::vec(c.at("suction_parameters"));
    s.pressure_parameters = ref::vec(c.at("pressure_parameters"));
    s.chi_le = ref::num(cascade.at("chi_le"));
    s.chi_te = ref::num(cascade.at("chi_te"));
    s.le_radius = ref::num(a.at("leading_edge_radius"));
    s.te_radius = ref::num(a.at("trailing_edge_radius"));
    // metadata only: the cascade-frame centres turned back to the blade frame
    const auto le = ref::vec(a.at("leading_edge_center")), te = ref::vec(a.at("trailing_edge_center"));
    s.le_center = T::rotate_by_stagger({le[0], le[1]}, -stagger);
    s.te_center = T::rotate_by_stagger({te[0], te[1]}, -stagger);
    return s;
}

/// The reference data and the port's cascade of it (solidity 1, stagger 33 deg), built once.
struct Reference
{
    ref::json family;
    T::Cascade cascade;
};

inline const Reference* reference()
{
    static std::optional<Reference> cached;
    static bool tried = false;
    if (!tried) {
        tried = true;
        auto fam = ref::load("grid_family");
        if (fam) {
            T::Cascade c = T::make_cascade(v1_section(fam->at("cascade")), 1.0, 33.0);
            cached = Reference{std::move(*fam), std::move(c)};
        }
    }
    return cached ? &*cached : nullptr;
}

#define TURBOLENZA_CASCADE_OR_SKIP(var)                                                                     \
    const fixture::Reference* var = fixture::reference();                                                  \
    if (!var) SKIP("reference data grid_family.json missing: run tools/turbolenza_reference/dump_v1.py")

}  // namespace fixture
