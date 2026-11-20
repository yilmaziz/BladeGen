/// @file occt/OcctBlade.hpp
/// @brief OpenCASCADE geometry of a blade: C2 NURBS skin, trimmed solid, hub
///        fillet, export and tessellation.
///
/// Skinning (L. Piegl & W. Tiller, The NURBS Book, §10.3 "skinned surfaces"):
///   1. every section loop is interpolated by a periodic cubic B-spline with the
///      SAME parameters (the average of the sections' chord-length parameters,
///      eq. 9.18), so all section curves share one knot vector;
///   2. the poles are interpolated spanwise, pole index by pole index, by C2
///      cubic B-splines with common (averaged chord-length) span parameters.
/// The surface is periodic in u (around the profile), C2 in u and v, and passes
/// through every section exactly.
///
/// Solid: the blade face (extended past hub and casing) splits a sector of the
/// annulus between hub and casing (or the tip-gap line), so the blade solid's end
/// caps lie exactly on the hub and casing surfaces of revolution.
#pragma once

#include "blade/Blade.hpp"
#include "blade/Meridional.hpp"

#include <Eigen/Dense>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include <array>
#include <string>
#include <vector>

namespace BladeGen::Occt {

using Math::Vec2;
using Math::Vec3;

struct SkinResult
{
    Handle(Geom_BSplineSurface) surface;
    std::vector<double> u;      ///< common chordwise parameters of the loop points
    std::vector<double> v;      ///< span parameter of every section
    double max_error = 0.0;     ///< max |S(u_j, v_k) - loop_k[j]| (interpolation check)
};

/// C2 periodic cubic interpolation at fixed parameters u_0 < ... < u_N (u_N closes
/// the loop).  OCCT's GeomAPI_Interpolate in periodic mode puts a double knot at the
/// seam, so its curves are only C1 there; here every knot is simple (C2 everywhere).
/// The interpolation matrix is factored once and reused for every point set.
class PeriodicInterpolator
{
public:
    explicit PeriodicInterpolator(const std::vector<double>& u);
    /// Curve through pts[0..N-1] at u_0..u_{N-1} (a closing pts[N] is ignored).
    [[nodiscard]] Handle(Geom_BSplineCurve) curve(const std::vector<Vec3>& pts) const;

private:
    int n_;
    TColStd_Array1OfReal knots_;
    TColStd_Array1OfInteger mults_;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu_;
};

/// Skin closed loops (last point == first, the same point count in every loop).
/// `u` (optional, N + 1 increasing values from 0 to 1) fixes the chordwise
/// parameters; by default they are the averaged normalised chord lengths.
[[nodiscard]] SkinResult skin_loops(const std::vector<std::vector<Vec3>>& loops,
                                    const std::vector<double>& u = {});

/// Loops extended past the first / last one: quadratic (C2) extrapolation of the
/// skin along its span parameter by `before` / `after` times the mean section spacing.
[[nodiscard]] std::vector<std::vector<Vec3>> extend_loops(const std::vector<std::vector<Vec3>>& loops,
                                                          double before, double after);

struct SolidOptions
{
    Blade::MeridionalCurve hub;      ///< hub contour (x, r)
    Blade::MeridionalCurve casing;   ///< casing contour (x, r)
    double tip_clearance = 0.0;      ///< blade tip = casing offset inward by this along its normal
    double axial_margin = 0.25;      ///< annulus sector extends this many axial chords past LE / TE
    double sector_margin_deg = 4.0;  ///< sector extends this past the blade's theta range
    double hub_bore_radius = 0.0;    ///< hub metal down to this radius (0: 85 % of the hub radius)
    double fillet_radius = 0.0;      ///< blade / hub fillet radius (0: none)
};

/// Everything the OCCT layer builds for one blade.
struct BladeShapes
{
    SkinResult skin;
    TopoDS_Face blade_face;   ///< the (extended) skin face
    TopoDS_Shape blade;       ///< blade solid between hub and casing / tip line
    TopoDS_Shape hub;         ///< hub sector solid (bore to hub surface)
    TopoDS_Shape blade_on_hub;///< hub sector fused with the blade, filleted when requested
    TopoDS_Shape blade_with_fillet; ///< blade + fillet standing on its hub footprint (fillet only)
    TopoDS_Face fillet_face;  ///< the fillet surface (fillet only)
    bool fillet_ok = false;   ///< the fillet was built and the result is valid
    std::string fillet_message;
    double fillet_g1_deg = 0.0;     ///< max normal angle fillet / walls at the contact loops
    double fillet_ball_error = 0.0; ///< max deviation of the contact distances from r
    double blade_volume = 0.0;
    double blade_area = 0.0;
    bool blade_valid = false;
};

/// Build the skin, the trimmed blade solid and (if requested) the filleted hub+blade.
[[nodiscard]] BladeShapes build_shapes(const std::vector<std::vector<Vec3>>& loops, const SolidOptions& opt,
                                       const Vec3& inside_point);

/// Same from a Blade (its sections, channel hub and shroud, tip clearance).
[[nodiscard]] BladeShapes build_blade(const Blade::Blade& blade, double fillet_radius = 0.0, int n_edge = 40,
                                      int n_body = 160);

/// A point strictly inside a section's profile (mid-chord, between SS and PS), 3-D.
[[nodiscard]] Vec3 inside_point(const Blade::Section3D& s);

/// Add a constant-radius fillet between a blade and the hub (BRepFilletAPI).
/// Returns the hub+blade solid; `ok` false (and the unfilleted fuse) on failure.
[[nodiscard]] TopoDS_Shape fillet_blade_hub(const TopoDS_Shape& blade, const TopoDS_Shape& hub, double radius,
                                            bool& ok, std::string& message);

// ── properties, checks, export ───────────────────────────────────────────────
[[nodiscard]] double volume(const TopoDS_Shape& s);
[[nodiscard]] double area(const TopoDS_Shape& s);
[[nodiscard]] bool is_valid(const TopoDS_Shape& s);
[[nodiscard]] int count_faces(const TopoDS_Shape& s);
[[nodiscard]] int count_solids(const TopoDS_Shape& s);

enum class Unit { Metre, Millimetre };
/// STEP AP214 (by extension .step/.stp), IGES (.igs/.iges), STL (.stl, binary),
/// BREP (.brep).  Geometry is in metres; `unit` is the unit written to the file.
void write_shape(const TopoDS_Shape& s, const std::string& path, Unit unit = Unit::Millimetre);
[[nodiscard]] TopoDS_Shape read_step(const std::string& path);

/// Triangulation of a shape for viewing: vertices (m), triangles, face id per triangle.
struct Mesh
{
    std::vector<Vec3> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::vector<int> face_of_triangle;
};
[[nodiscard]] Mesh tessellate(const TopoDS_Shape& s, double linear_deflection, double angular_deflection = 0.3);

/// Minimum distance between two shapes (BRepExtrema).
[[nodiscard]] double distance(const TopoDS_Shape& a, const TopoDS_Shape& b);
/// Minimum distance from points to a shape.
[[nodiscard]] std::vector<double> distances(const std::vector<Vec3>& pts, const TopoDS_Shape& s);

}  // namespace BladeGen::Occt
