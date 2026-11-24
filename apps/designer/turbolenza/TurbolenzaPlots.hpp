// TurbolenzaPlots: what the TURBOlenza window draws, as PlotWidget series and fills.  Drawing
// only: every number comes from BladeGen (the profile), the TURBOlenza core (cascade, grid,
// solution) or the run records; nothing here computes flow.
//
//   profile     the blade section (suction side blue, pressure side red, camber line)
//   cascade     three blades of the cascade with the geometric throat and the inlet direction
//   grid        the (m, eta) grid of one passage: station and eta lines, the first and last blade
//               stations, the two blades
//   residuals   density, velocity, Kutta changes and the linear residual per pass (log y)
//   contours    a node field as filled contour bands with iso-lines: marching squares on the
//               structured grid cells (linear interpolation along each cell edge; a saddle cell is
//               resolved by its centre value), periodic copies one pitch up and down, the blades on top, a colour bar
//               with the turbo map
//   loading     a blade-surface quantity against x/c, suction side blue, pressure side red
//   boundary    the boundary layer of both surfaces (H, theta/c, delta*/c, c_f, Re_theta) with
//   layer       Aungier's Fig. 5-15 overlaid; the coupling history (change of m, beta_out per k)
//   mesh study  beta_out against dm/c (or the node count), one labelled marker per level; the
//               blade loading of every level on one axes (one colour per level)
#pragma once

#include "PlotWidget.hpp"

#include "design/BladeDesign.hpp"
#include "turbolenza/Case.hpp"
#include "turbolenza/MeshStudy.hpp"

#include <QString>

#include <vector>

namespace Designer::Tl {

namespace T = BladeGen::Turbolenza;

inline const QColor kSuction(31, 119, 180);   // v1 plots.py SUCTION_COLOR
inline const QColor kPressure(214, 39, 40);   // v1 PRESSURE_COLOR

void drawProfile(PlotWidget* w, const BladeGen::Design::Section2D& s, const QString& title);
void drawCascade(PlotWidget* w, const T::Cascade& c, double inlet_angle_deg);
void drawGrid(PlotWidget* w, const T::Grid& g);
/// Not in v1: the stream sheet of a blade on its stream surface, r(m) / r_LE and b(m) / b(m_lo) against
/// m / c on [m_lo, m_hi] (chords), the blade's LE (m = 0) and TE (m = m_te) marked.
void drawStreamSheet(PlotWidget* w, const T::StreamSheet& sheet, double m_lo, double m_hi, double m_te);
/// Not in v1: the inlet velocity triangle (TURBOlenza's frame: m to the right, +theta up), C = W + U.
void drawVelocityTriangle(PlotWidget* w, const T::VelocityTriangle& t);

/// The residual history (the same records the run printed).  tol: density, velocity, Kutta
/// tolerances (a NaN skips the line).
void drawResiduals(PlotWidget* w, const std::vector<T::PassRecord>& records, double tol_rho, double tol_w, double tol_k,
                   const QString& status);

enum class ContourField { Pressure, Mach };
/// bands: number of colour bands; isolines: thin lines at the band edges.
void drawContours(PlotWidget* w, const T::SolutionData& s, ContourField f, int bands, bool isolines, const QString& title);

enum class LoadingQuantity { Mach, Pressure, DunavantS };
/// inviscid: a viscous run's k = 0 loading, drawn dashed beside the reported (coupled) one.
void drawLoading(PlotWidget* w, const T::SolutionData& s, LoadingQuantity q, const QString& title,
                 const T::ViscousData* inviscid = nullptr);

/// The boundary layer (Aungier 2003 Sec. 5.6) of both surfaces against x/c, transition (o) and
/// separation (x) marked; overlay (H only): Aungier (2003) Fig. 5-15, digitized
/// (tests/turbolenza/validation/aungier2003_fig5_15), as markers.
enum class BlQuantity { H, ThetaOverC, DeltaStarOverC, Cf, ReTheta };
void drawBoundaryLayer(PlotWidget* w, const T::ViscousData& v, double chord, BlQuantity q, bool overlay, const QString& title);
/// The coupling history: max relative change of m per BL iteration (log y) with the tolerance.
void drawBlChange(PlotWidget* w, const T::ViscousData& v, double tolerance);
/// The coupling history: beta_out of each BL iteration's inviscid solve.
void drawBlBetaOut(PlotWidget* w, const T::ViscousData& v);

/// Mesh study: beta_out [deg] of every level with a solution against dm/c (or the total nodes when
/// x_nodes), labelled with the level; unconverged levels as red crosses, named in the title.
void drawStudyBetaOut(PlotWidget* w, const std::vector<T::MeshStudyEntry>& levels, bool x_nodes);
/// Mesh study: the blade loading of every level with surfaces, one colour per level (legend),
/// suction side solid, pressure side dashed.
void drawStudyLoading(PlotWidget* w, const std::vector<T::MeshStudyEntry>& levels, LoadingQuantity q);

/// The closed outline of blade k (k = 0 the reference blade) from the grid boundaries (v1
/// plots.plot_mach_contour: the suction samples reversed, then the pressure samples).
[[nodiscard]] QPolygonF bladeOutline(const T::PassageBoundaries& b, int k);

/// Marching squares on one field (exposed for the self-test): the filled band polygons and the
/// iso-line segments (pairs of points) of the levels lo + (hi - lo) b / bands, b = 1 .. bands - 1.
struct ContourGeometry
{
    QVector<Fill> fills;
    QVector<QPointF> segments;  ///< a, b, a, b, ...
};
[[nodiscard]] ContourGeometry contourGeometry(const T::Field& x, const T::Field& y, const T::Field& v, double lo, double hi,
                                              int bands, double dy);
[[nodiscard]] QColor bandColor(int b, int bands);

}  // namespace Designer::Tl
