/// @file io/Rotor67.hpp
/// @brief NASA Rotor 67: the 14 published blade sections and the flowpath.
///
/// Source: A. J. Strazisar et al., "Laser anemometer measurements in a transonic
/// axial-flow fan rotor", NASA TP-2879 (1989), blade-section and flowpath tables
/// (transcribed in data/rotor67/Rotor_14_Sections_and_Flowpath_EN.txt; cm, rad, deg).
///
/// Section construction (port of blade2d/newblade rotor67/rotor67_sections.py):
///   m          arc length along the (Z, R) stream line through the K points
///              (C2 cubic, chord-length parameters), m = 0 at K = 1;
///   surfaces   theta_j(m): C2 cubic through the K points, end slopes clamped to
///              tan(BETA_j) / R (BETA_j from the meridional direction);
///   edges      circles in the local plane (m - m_c, s = R_c (theta - theta_c));
///              tangency = foot of the perpendicular from the centre to each surface;
///   profile    surface 1 LE -> TE, TE arc, surface 2 TE -> LE, LE arc.
/// Lengths are returned in metres.
#pragma once

#include "blade/Meridional.hpp"
#include "math/BSpline.hpp"

#include <string>
#include <vector>

namespace BladeGen::IO {

using Math::Vec2;
using Math::Vec3;

struct EdgeCircle
{
    double rho = 0.0;  ///< radius [m]
    double rc = 0.0;   ///< centre radius R [m]
    double zc = 0.0;   ///< centre axial Z [m]
    double mc = 0.0;   ///< centre meridional distance [m]
    double thc = 0.0;  ///< centre angle [rad]
};

struct Rotor67Section
{
    int index = 0;                    ///< 1 (hub) ... 14 (tip)
    std::vector<double> z, r;         ///< stream-line points [m]
    std::vector<double> th1, th2;     ///< surface angles [rad]
    std::vector<double> beta1, beta2; ///< surface angles from the meridional direction [rad]
    EdgeCircle le, te;
    Blade::MeridionalCurve streamline;
    std::vector<double> m;            ///< meridional distance of the K points [m]
    Math::BSpline<1> surf1, surf2;    ///< theta_j(m)
    double m_tan[2][2] = {{0, 0}, {0, 0}};  ///< tangency m of surface j (0,1) at LE (0) / TE (1)

    /// Closed profile in (m, theta) (n_surface points per surface, n_arc per arc).
    [[nodiscard]] std::vector<Vec2> profile_m_theta(int n_surface = 161, int n_arc = 33) const;
    /// The same loop in (m', theta) of the section's stream line.
    [[nodiscard]] std::vector<Vec2> profile_conformal(int n_surface = 161, int n_arc = 33) const;
    /// Cartesian (x = Z, y = R sin theta, z = R cos theta) [m].
    [[nodiscard]] std::vector<Vec3> profile_xyz(int n_surface = 161, int n_arc = 33) const;
    [[nodiscard]] double theta(int surface, double m) const;
};

struct Rotor67Data
{
    std::vector<Rotor67Section> sections;  ///< hub -> tip
    std::vector<Vec2> hub, shroud;          ///< flowpath (Z, R) [m]
    int n_blades = 22;
    double tip_clearance = 0.0;             ///< [m]
    double rotating_hub[2] = {0, 0};        ///< Z range of the rotating hub [m]

    [[nodiscard]] Blade::MeridionalCurve hub_curve() const;
    [[nodiscard]] Blade::MeridionalCurve shroud_curve() const;
};

/// Read and construct everything from the published table file.
[[nodiscard]] Rotor67Data read_rotor67(const std::string& path);

}  // namespace BladeGen::IO
