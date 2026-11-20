/// @file io/BladeJson.hpp
/// @brief JSON files for profiles and blades (nlohmann/json).
///
/// Profile: a flat object of the scalar parameters (Profile::param_table names,
/// e.g. "le.radius") plus "connector", "ps.mprime_ctrl", "thickness.enforce_c1".
///
/// Blade:
/// @code
/// { "name": "...", "n_blades": 22, "tip_clearance": 0.001, "interpolation": "cubic",
///   "span_extension": 0.04, "n_sections": 15,
///   "channel": { "hub": CURVE, "shroud": CURVE },
///   "stations": [ { "span": 0.0, "profile": PROFILE,
///                   "placement": { "m_le": 0.3, "m_te": 0.7, "stagger_deg": 35,
///                                  "theta_le_deg": 0, "mirror": false },
///                   "streamline": CURVE  (optional) } ] }
/// @endcode
/// CURVE is either {"degree", "knots", "control_points": [[x, r], ...]} (exact
/// B-spline) or {"points": [[x, r], ...]} (C2 cubic through the points).  Lengths in
/// metres, angles in degrees.
///
/// blade2d compatibility: profile_from_blade2d() reads the flat
/// "<name>_blade_params.json" written by blade2d's inverse fitters.
#pragma once

#include "blade/Blade.hpp"
#include "profile/ProfileParams.hpp"

#include <nlohmann/json.hpp>

#include <string>

namespace BladeGen::IO {

[[nodiscard]] nlohmann::json to_json(const Profile::ProfileParams& p);
[[nodiscard]] Profile::ProfileParams profile_from_json(const nlohmann::json& j);
/// blade2d blade_params.json (le_radius, th_t_max, ps_mprime_ctrl_0, ...) -> G2Legacy profile.
[[nodiscard]] Profile::ProfileParams profile_from_blade2d(const nlohmann::json& j);

[[nodiscard]] nlohmann::json to_json(const Blade::MeridionalCurve& c);
[[nodiscard]] Blade::MeridionalCurve curve_from_json(const nlohmann::json& j);

[[nodiscard]] nlohmann::json to_json(const Blade::BladeDefinition& d);
[[nodiscard]] Blade::BladeDefinition blade_from_json(const nlohmann::json& j);

void write_blade_json(const Blade::BladeDefinition& d, const std::string& path);
[[nodiscard]] Blade::BladeDefinition read_blade_json(const std::string& path);
void write_profile_json(const Profile::ProfileParams& p, const std::string& path);
/// Reads a BladeGen profile JSON, or a blade2d blade_params.json (detected by its keys).
[[nodiscard]] Profile::ProfileParams read_profile_json(const std::string& path);

}  // namespace BladeGen::IO
