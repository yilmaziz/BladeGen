/// @file blade/Blade.hpp
/// @brief 3-D blades: profiles placed on stream surfaces and stacked hub to tip.
///
/// One section = a Profile2D (unstaggered, chord 1) placed in the conformal plane
/// (m', theta) of its streamline:
///
///     m'    = m'_LE + c' ( xi cos(gamma) - s eta sin(gamma) )
///     theta = theta_LE + c' ( xi sin(gamma) + s eta cos(gamma) )
///
/// with gamma the stagger angle (from the m' axis toward +theta), s = +1 when the
/// suction side faces +theta (s = -1 mirrors the profile), and the conformal chord
/// c' fixed by the meridional positions of the LE and TE:  c' cos(gamma) =
/// m'(TE) - m'(LE).  The point (m', theta) then maps to Cartesian
///
///     x = x(m'),  y = r(m') sin(theta),  z = r(m') cos(theta)     (axis = x)
///
/// which is the convention of the blade2d / Rotor 67 scripts (theta = 0 on +z).
/// Because the map is conformal, the profile angles in (m', theta) are the true
/// flow angles of the stream surface.
#pragma once

#include "blade/Meridional.hpp"
#include "profile/Profile2D.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace BladeGen::Blade {

using Math::Vec3;

/// Where and how a profile sits on its stream surface.
struct SectionPlacement
{
    double m_le = 0.25;        ///< LE meridional position, fraction of the streamline length
    double m_te = 0.75;        ///< TE meridional position, fraction of the streamline length
    double stagger_deg = 30.0; ///< chord angle from the m' axis toward +theta [deg]
    double theta_le_deg = 0.0; ///< circumferential position of the LE point [deg] (lean)
    bool mirror = false;       ///< true: suction side toward -theta
};

/// A design station: span fraction, profile and placement (and optionally its own
/// streamline; otherwise the channel streamline at that span is used).
struct BladeStation
{
    double span = 0.0;
    Profile::ProfileParams profile{};
    SectionPlacement placement{};
    std::optional<MeridionalCurve> streamline{};
};

enum class SpanInterpolation { Linear, Cubic };

struct BladeDefinition
{
    std::string name = "blade";
    MeridionalChannel channel{};
    std::vector<BladeStation> stations;   ///< ordered hub -> tip (span increasing)
    int n_blades = 20;
    SpanInterpolation interpolation = SpanInterpolation::Cubic;
    /// Sections are built over span in [-ext, 1 + ext] so the blade surface crosses
    /// hub and casing (the solid is then trimmed on them).
    double span_extension = 0.04;
    int n_sections = 15;                  ///< sections between span 0 and 1 for skinning
    double tip_clearance = 0.0;           ///< radial-normal gap below the casing [m]

    void validate() const;
    /// A demonstration compressor rotor (NACA 65-like profiles, twisted, 0.25 - 0.40 m).
    [[nodiscard]] static BladeDefinition demo_rotor();
};

/// One built section: its profile, placement and points.
struct Section3D
{
    double span = 0.0;
    Profile::ProfileParams params{};
    SectionPlacement placement{};
    MeridionalCurve streamline{};
    double mprime_le = 0.0;    ///< absolute m' of the LE point on the streamline
    double chord_mprime = 0.0; ///< conformal chord c'
    std::vector<Vec2> unstaggered;   ///< loop in (xi, eta)
    std::vector<Vec2> conformal;     ///< loop in (m', theta [rad])
    std::vector<Vec3> xyz;           ///< loop in Cartesian [m]
    std::vector<int> piece;          ///< Profile::Piece of every loop point
};

/// Map unstaggered (xi, eta) to (m', theta) and Cartesian for one placement.
struct SectionMap
{
    double mprime_le = 0.0, chord = 0.0, gamma = 0.0, theta_le = 0.0, sign = 1.0;
    const MeridionalCurve* streamline = nullptr;

    [[nodiscard]] Vec2 conformal(const Vec2& xieta) const;
    [[nodiscard]] Vec3 cartesian(const Vec2& mprime_theta) const;
};

[[nodiscard]] SectionMap make_section_map(const MeridionalCurve& sl, const SectionPlacement& pl);

/// Build one section on a streamline (n_edge / n_body: Profile2D::sample_loop layout).
[[nodiscard]] Section3D build_section(const Profile::ProfileParams& prof, const SectionPlacement& pl,
                                      const MeridionalCurve& sl, double span, int n_edge = 40, int n_body = 160);

/// Cartesian point of (x, r, theta) with the BladeGen convention (axis = x).
[[nodiscard]] inline Vec3 cartesian(double x, double r, double theta)
{
    return {x, r * std::sin(theta), r * std::cos(theta)};
}

class Blade
{
public:
    explicit Blade(BladeDefinition def);

    [[nodiscard]] const BladeDefinition& definition() const noexcept { return def_; }
    /// Profile and placement interpolated in span (clamped to the end stations).
    [[nodiscard]] Profile::ProfileParams profile_at(double span) const;
    [[nodiscard]] SectionPlacement placement_at(double span) const;
    [[nodiscard]] MeridionalCurve streamline_at(double span) const;
    [[nodiscard]] Section3D section(double span, int n_edge = 40, int n_body = 160) const;
    /// Sections for skinning: span from -ext to 1 + ext.
    [[nodiscard]] std::vector<Section3D> sections(int n_edge = 40, int n_body = 160) const;

private:
    [[nodiscard]] double interp(const std::vector<double>& y, double span) const;

    BladeDefinition def_;
    std::vector<double> spans_;
};

}  // namespace BladeGen::Blade
