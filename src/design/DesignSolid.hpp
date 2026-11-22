/// @file design/DesignSolid.hpp
/// @brief OCCT solid of a designer blade: the mesh rows skinned by the C2 NURBS surface of
///        occt/OcctBlade (extended past hub and casing, then trimmed by the revolved walls).
#pragma once

#include "design/BladeDesign.hpp"
#include "occt/OcctBlade.hpp"

#include <gp_Trsf.hxx>

namespace BladeGen::Design {

/// Builds the solid from a design (any construction method, NACA 65 sections included).
/// `scale` converts the design length unit to metres (1e-3 for mm).  The design is
/// re-sampled with `rows` uniform spanwise rows for the skin.
/// `blade_only`: Occt::SolidOptions::blade_only (the viewer).
[[nodiscard]] Occt::BladeShapes build_design_solid(const BladeDesign& d, Evaluator& ev, double scale = 1e-3, int rows = 11,
                                                   double fillet_radius = 0.0, bool blade_only = false);

/// The OCCT 3-D view's geometry: one blade solid or all `n_blades` copies around the annulus,
/// and the hub and casing as surfaces of revolution of the whole meridional channel, over one
/// pitch centred on the blade (single blade) or over 360 deg (annulus).  Metres.
struct SceneOptions
{
    bool annulus = false;  ///< all n_blades blades (else one blade)
    int rows = 11;         ///< spanwise rows of the skin (build_design_solid)
};

struct Scene
{
    Occt::BladeShapes shapes;  ///< the single blade (shapes.blade) and its skin (blade_only: no properties)
    TopoDS_Shape blades;       ///< compound of 1 or n_blades located copies of shapes.blade (shared geometry)
    TopoDS_Shape hub;          ///< hub face of revolution
    TopoDS_Shape casing;       ///< casing face of revolution
    int n_blades = 0;
    int copies = 0;            ///< blades in `blades`
    double pitch = 0.0;        ///< 2 pi / n_blades [rad]
    double theta_blade = 0.0;  ///< theta of the blade's inside point (mid row, mid chord) [rad]
    double theta_start = 0.0;  ///< the walls start at this theta ...
    double sweep = 0.0;        ///< ... and sweep this angle toward +theta [rad]
};

/// Rotation of a copy by `angle` about the machine axis toward +theta (y = r sin theta, z = r cos theta).
[[nodiscard]] gp_Trsf theta_rotation(double angle);

/// Builds the scene; throws std::invalid_argument for n_blades < 1 or a design with errors.
[[nodiscard]] Scene build_scene(const BladeDesign& d, Evaluator& ev, const SceneOptions& opt = {});

/// Wavefront OBJ of the blade, hub and casing grids (quads, design units).
void write_obj(const DesignResult& r, const std::string& path);

/// CSV of every layer: unstaggered (x/c, y/c), conformal (m', theta), Cartesian (x, y, z).
void write_layer_csv(const BladeDesign& d, const DesignResult& r, const std::string& path);

}  // namespace BladeGen::Design
