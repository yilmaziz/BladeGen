/// @file design/DesignSolid.hpp
/// @brief OCCT solid of a designer blade: the mesh rows skinned by the C2 NURBS surface of
///        occt/OcctBlade (extended past hub and casing, then trimmed by the revolved walls).
#pragma once

#include "design/BladeDesign.hpp"
#include "occt/OcctBlade.hpp"

namespace BladeGen::Design {

/// Builds the solid from a design (any construction method, NACA 65 sections included).
/// `scale` converts the design length unit to metres (1e-3 for mm).  The design is
/// re-sampled with `rows` uniform spanwise rows for the skin.
[[nodiscard]] Occt::BladeShapes build_design_solid(const BladeDesign& d, Evaluator& ev, double scale = 1e-3, int rows = 11,
                                                   double fillet_radius = 0.0);

/// Wavefront OBJ of the blade, hub and casing grids (quads, design units).
void write_obj(const DesignResult& r, const std::string& path);

/// CSV of every layer: unstaggered (x/c, y/c), conformal (m', theta), Cartesian (x, y, z).
void write_layer_csv(const BladeDesign& d, const DesignResult& r, const std::string& path);

}  // namespace BladeGen::Design
