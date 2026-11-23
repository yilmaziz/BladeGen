/// @file io/SectionExport.hpp
/// @brief Versioned JSON of section points for downstream tools (grid generators):
///        `bladegen export-sections`.
///
/// Format "bladegen.sections", version 1.  Two kinds:
///
///  * "blade_design": sections of a blade design (BladeDefinition) at requested span
///    fractions.  Per section, the suction- and pressure-side points from the LE point to
///    the TE point in the blade-to-blade plane (m', theta) and the same points in 3-D
///    (z, r, theta) on the section's stream surface, from one evaluation: the points are
///    Blade::section()'s loop (Profile2D::sample_loop, the layout of the 3-D skin), with
///    (z, r) = streamline(m(m')) exactly as SectionMap::cartesian computes it.  The section's
///    streamline is written with its m and m' tables, and the hub / shroud channel in (z, r).
///  * "naca65_cascade": a planar NACA 65 cascade section in cascade coordinates (axial x,
///    tangential y), scaled by the chord, staggered about the LE point and with the pitch
///    chord / solidity.
///
/// Coordinates: z (and x of the cascade) is the machine axis, BladeGen's internal x; r the
/// radius; theta the circumferential angle of BladeGen, with the Cartesian point
/// (x_int, y_int, z_int) = (z, r sin theta, r cos theta), i.e. theta = 0 on +z_int.  In the
/// designer's machine labels (Z axial, X = r sin theta, Y = r cos theta) theta turns from +Y
/// toward +X.  Lengths in metres, angles in radians; m' = \int dm / r is dimensionless.
///
/// Numbers: every double is written with 17 significant digits (std::setprecision(17), the
/// "C" locale), which reads back bit-identically with any correctly rounded strtod.
#pragma once

#include "blade/Blade.hpp"
#include "naca/Naca65.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace BladeGen::IO {

inline constexpr const char* kSectionsFormat = "bladegen.sections";
inline constexpr int kSectionsVersion = 1;

/// Sampling of a blade-design export.
struct DesignExportOptions
{
    std::vector<double> spans{0.0};  ///< span fractions in [0, 1]; 0 = hub section
    int n_edge = 40;                 ///< points per half edge connector (Profile2D::sample_loop)
    int n_body = 160;                ///< points on each side's body
    int n_channel = 201;             ///< points of the hub / shroud and streamline tables (equal steps in m)
};

/// A planar NACA 65 cascade.
struct CascadeExportOptions
{
    Naca::Naca65Params naca{};
    double chord = 1.0;        ///< [m], > 0
    double stagger_deg = 0.0;  ///< chord angle from the +x (axial) axis toward +y, |stagger| < 90
    double solidity = 1.0;     ///< chord / pitch, > 0
    int n_edge = 40;           ///< points of each LE and TE piece per side (Naca65Section::sample_side)
    int n_stations = 101;      ///< points on each side's body (cosine-clustered in sqrt(x))
};

/// The export of a blade design.  `source_name` is written as the source (give a file
/// name, not a path); `generator_version` is the BladeGen version.
[[nodiscard]] nlohmann::ordered_json export_design_sections(const Blade::BladeDefinition& def,
                                                            const DesignExportOptions& opt,
                                                            const std::string& source_name,
                                                            const std::string& generator_version);

/// The export of a NACA 65 cascade section.
[[nodiscard]] nlohmann::ordered_json export_naca65_cascade(const CascadeExportOptions& opt,
                                                           const std::string& generator_version);

/// Serialise with 17 significant digits per double (objects indented, number arrays on
/// one line).  Throws std::invalid_argument on a non-finite number.
[[nodiscard]] std::string dump_sections(const nlohmann::ordered_json& j);
void write_sections(const nlohmann::ordered_json& j, const std::string& path);

}  // namespace BladeGen::IO
