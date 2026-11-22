/// @file design/DesignSolid.hpp
/// @brief OCCT solid of a designer blade: the mesh rows skinned by the C2 NURBS surface of
///        occt/OcctBlade (extended past hub and casing, then trimmed by the revolved walls).
#pragma once

#include "design/BladeDesign.hpp"
#include "occt/OcctBlade.hpp"

#include <gp_Trsf.hxx>

#include <array>
#include <vector>

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

/// Triangles with smooth vertex normals (metres).
struct TriangleMesh
{
    std::vector<Vec3> p;                      ///< vertices
    std::vector<Vec3> n;                      ///< unit vertex normals
    std::vector<std::array<int, 3>> tri;      ///< vertex indices, counter-clockwise seen along n
};

/// The fast display of the OCCT 3-D view, made from the design's own mesh rows (no OCCT
/// modelling, no Boolean, no BRepMesh).  The blade is the row grid in triangles (it passes
/// through every mesh point; flat between them) with flat caps on the hub and tip rows; hub and
/// casing are grids on their surfaces of revolution over one pitch or 360 deg, as in Scene.
/// Normals face out of the blade and into the passage.  Metres, like Scene.
struct SceneMesh
{
    TriangleMesh blade;                          ///< one blade (the viewer turns copies by k * pitch)
    TriangleMesh hub, casing;
    std::vector<std::vector<Vec3>> blade_lines;  ///< wireframe: the rows and every few loop columns
    std::vector<std::vector<Vec3>> hub_lines, casing_lines;
    std::vector<std::vector<Vec3>> blade_edges;  ///< outline: hub and tip rows, LE and TE lines
    std::vector<std::vector<Vec3>> hub_edges, casing_edges;
    int n_blades = 0;
    int copies = 0;
    double pitch = 0.0, theta_blade = 0.0, theta_start = 0.0, sweep = 0.0;
};

/// Builds the fast display from an evaluated design (r = Evaluator::evaluate(d), r.ok()).
/// Throws std::invalid_argument for a result without a blade mesh or n_blades < 1.
[[nodiscard]] SceneMesh build_scene_mesh(const BladeDesign& d, const DesignResult& r, bool annulus);

/// theta of a design's blade: the inside point (mid row, mid chord) [rad].
[[nodiscard]] double blade_theta(const DesignResult& r);

/// Rotation of a copy by `angle` about the machine axis toward +theta (y = r sin theta, z = r cos theta).
[[nodiscard]] gp_Trsf theta_rotation(double angle);

/// Builds the scene; throws std::invalid_argument for n_blades < 1 or a design with errors.
[[nodiscard]] Scene build_scene(const BladeDesign& d, Evaluator& ev, const SceneOptions& opt = {});

/// Wavefront OBJ of the blade, hub and casing grids (quads, design units).
void write_obj(const DesignResult& r, const std::string& path);

/// CSV of every layer: unstaggered (x/c, y/c), conformal (m', theta), Cartesian (x, y, z).
void write_layer_csv(const BladeDesign& d, const DesignResult& r, const std::string& path);

}  // namespace BladeGen::Design
