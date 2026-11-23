/// @file turbolenza/Cascade.hpp
/// @brief Cascade generation: one blade section set at a stagger and repeated at a pitch,
///        with its geometric throat.  The blade itself comes from BladeGen.
///
/// Port of the cascade part of v1 geometry.py (make_cascade, place_cascade,
/// find_geometric_throat / closest_points, validate_cascade's frame checks; Aungier 2003
/// Secs. 4.9, 4.10).  v1's own NACA 65 construction is NOT ported: TURBOlenza v2 takes the blade
/// section from BladeGen (BladeGenSource.hpp) or, for verification, v1's frozen C2 curves
/// (tests/turbolenza/data).
///
/// Frames and signs (v1, unchanged).  Blade frame: x along the chord, y normal to it, LE (0, 0),
/// TE (1, 0), camber toward +y, chord 1; the suction surface is the +y side.  Cascade frame: x
/// axial (the grid's m), y tangential.  The blade turns counterclockwise by the stagger about its
/// LE, X = x cos(stagger) - y sin(stagger), Y = x sin(stagger) + y cos(stagger), and its neighbour
/// is the same blade moved by periodic_translation = (0, pitch), pitch = 1 / solidity.  Angles in
/// radians from the axial direction, positive toward +y.  The blade angles kappa_le = stagger +
/// chi_le and kappa_te = stagger - chi_te are where the grid's side lines meet the blade (v1: the
/// equivalent circular arc of NACA 65, Aungier Eq. 4-6; a BladeGen custom profile: its camber-line
/// angles at the LE and TE).
#pragma once

#include "turbolenza/Curve.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// A cascade or blade that breaks a geometric invariant (v1 GeometryError).
class GeometryError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// One blade section in the blade frame (chord 1, LE (0, 0), TE (1, 0), suction side +y): the
/// input of the cascade generation.
struct BladeSection
{
    std::string designation;        ///< e.g. "NACA 65-(12A10)10 (kovach_sandercock)", for reports
    ParametricBladeCurve suction;   ///< LE -> TE, u in [0, 1]; X(0) = (0, 0), X(1) = (1, 0)
    ParametricBladeCurve pressure;  ///< LE -> TE
    std::vector<double> suction_parameters;   ///< u of the sampled surface points (v1: where the C1 vertices are)
    std::vector<double> pressure_parameters;
    double chi_le = 0.0;            ///< LE blade angle from the chord line [rad] (> 0 toward +y)
    double chi_te = 0.0;            ///< TE blade angle from the chord line [rad] (> 0 toward -y)
    double le_radius = 0.0;         ///< chord fractions; metadata for reports
    double te_radius = 0.0;
    Point le_center{0.0, 0.0}, te_center{1.0, 0.0};
};

/// v1 CascadeCoordinates: every length of a cascade in one unit (here always chords: the solver
/// takes the cascade chord-normalized and the true chord as c_ref).
struct CascadeCoordinates
{
    double chord = 1.0;
    double pitch = 1.0;
    std::vector<Point> suction_surface_unstaggered, pressure_surface_unstaggered;
    std::vector<Point> suction_surface, pressure_surface;   ///< cascade frame, LE -> TE
    Point leading_edge_point{0, 0}, trailing_edge_point{0, 0};
    Point leading_edge_center{0, 0};
    double leading_edge_radius = 0.0;
    Point trailing_edge_center{0, 0};
    double trailing_edge_radius = 0.0;
    Point periodic_translation{0, 1};   ///< neighbour blade = reference blade + this
    double throat_opening = 0.0;
    Point throat_point_reference{0, 0};  ///< on the suction surface of the reference blade
    Point throat_point_neighbor{0, 0};   ///< on the pressure surface of the neighbour
};

/// v1 Cascade (C2 only, chord-normalized: v1's active == normalized).
struct Cascade
{
    std::string designation;
    double solidity = 1.0;
    double stagger = 0.0;               ///< radians
    double chi_le = 0.0, chi_te = 0.0;
    double camber_angle = 0.0;          ///< chi_le + chi_te
    double kappa_le = 0.0;              ///< stagger + chi_le
    double kappa_te = 0.0;              ///< stagger - chi_te
    double throat_over_pitch = 0.0;     ///< geometric (Aungier Sec. 4.9): the preferred value
    double throat_over_chord = 0.0;
    CascadeCoordinates active;
    std::shared_ptr<const SurfaceCurves> active_curves;  ///< the authoritative curves (cascade + blade frame)
};

/// The checked cascade of `section` at `solidity` (chord / pitch) and `stagger_deg` (from the
/// axial direction; the one degree-to-radian conversion).  Throws std::invalid_argument for an
/// invalid input and GeometryError when neighbouring blades overlap or an invariant breaks.
[[nodiscard]] Cascade make_cascade(const BladeSection& section, double solidity, double stagger_deg);

/// Shortest distance between two polylines treated as smooth curves (cubic splines through the
/// points), with its two end points (v1 closest_points).
struct ClosestPoints
{
    double distance = 0.0;
    Point a, b;
};
[[nodiscard]] ClosestPoints closest_points(const std::vector<Point>& curve_a, const std::vector<Point>& curve_b);

}  // namespace BladeGen::Turbolenza
