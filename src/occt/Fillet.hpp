/// @file occt/Fillet.hpp
/// @brief Constant-radius rolling-ball fillet between a blade and the hub.
///
/// OCCT's BRepFilletAPI fails on this junction (walking failure) — the same
/// finding as the blade2d rotor67_passage_solid.py script — so the fillet is
/// constructed explicitly, as in that script:
///
///   spine   at every station u around the blade root the ball centre
///           c = S(u, v*) + r n(u, v*)  (n: unit normal of the blade surface on the
///           fluid side) is placed at distance r from the hub, i.e. on the hub
///           contour offset outward by r in the meridional plane; v* by Brent.
///   contact the foot of c on the hub (meridional normal) and S(u, v*) on the blade;
///   surface the circular arcs from hub contact to blade contact, skinned into a
///           periodic C2 B-spline surface (u around the blade, v along the arc),
///           tangent to both walls (G1) by construction;
///   B-rep   the hub face gets the hub contact loop as a hole, the blade face is
///           rebuilt from the blade contact loop upward, and the fillet face shares
///           both loops with them (exact shared edges, no tangent booleans).
#pragma once

#include "blade/Meridional.hpp"

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include <string>
#include <vector>

namespace BladeGen::Occt {

struct FilletOptions
{
    double radius = 0.0;    ///< ball radius [m]
    int stations = 1440;    ///< ball positions around the root (G1 error ~0.2 deg at 1440)
    int arc_rows = 15;      ///< points across each arc (>= 3)
    double tolerance = 1e-6;
};

struct FilletResult
{
    bool ok = false;
    std::string message;
    TopoDS_Shape hub_blade;      ///< hub sector + blade + fillet, one solid
    TopoDS_Shape blade_fillet;   ///< blade + fillet standing on its hub footprint
    TopoDS_Face fillet_face;
    double max_ball_error = 0.0; ///< max | |c - contact| - r | over both contacts
    double max_g1_deg = 0.0;     ///< max angle between fillet and wall normals at the contacts
    double min_arc_deg = 0.0;    ///< smallest arc opening (pi - corner angle)
    double max_arc_deg = 0.0;
};

/// blade: the blade solid (lateral face + caps), blade_face: its skin face (extended
/// below the hub), hub_solid: the hub sector, hub: the hub contour (x, r).
[[nodiscard]] FilletResult rolling_ball_fillet(const TopoDS_Shape& blade, const TopoDS_Face& blade_face,
                                               const TopoDS_Shape& hub_solid, const Blade::MeridionalCurve& hub,
                                               const FilletOptions& opt);

}  // namespace BladeGen::Occt
