/// @file turbolenza/GridData.hpp
/// @brief The data of one blade-passage grid: the contract between grid generation
///        (Grid.hpp, port of v1 grid.py) and every solver module (metrics .. postprocess).
///
/// Port of the v1 dataclasses PassageBoundaries, WallSide, WallIdentity and Grid
/// (turbolenza/src/turbolenza/b2b/grid.py; Aungier 2003 Secs. 5.2, 5.3, Figs. 5-3, 5-5).
/// Plain data: the generator fills it and checks its invariants; the solver modules only read
/// it.  Where v1 evaluates a curve on demand (Grid.wall_curve_derivative), the grid stores the
/// values once, at build time (`wall_d`), so a module never needs the curves themselves.
///
/// Conventions (v1 grid.py, unchanged):
///  * planar cascade: m is the axial coordinate, y = r theta the tangential one; lengths in the
///    active unit of the cascade (chord 1 in the solver: the cascade is chord-normalized and the
///    true chord enters as c_ref in the coefficients);
///  * eta = (theta - theta0(m)) / (theta1(m) - theta0(m)) (Aungier Eq. 5-2): eta = 0 is the lower
///    wall, the suction surface of the reference blade, eta = 1 the upper wall, the pressure
///    surface of its +pitch neighbour (Fig. 5-3);
///  * angles in radians from the meridional (axial) direction, positive toward +y;
///  * i along m (n_m stations, constant dm), j along eta (n_eta nodes, constant deta);
///  * stations 0 .. first_blade_i - 1 are the inlet periodic stations, first_blade_i ..
///    last_blade_i the blade stations, the rest the exit periodic stations.
#pragma once

#include "turbolenza/Field.hpp"

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

struct SurfaceCurves;  // Curve.hpp: the authoritative C2 blade surfaces
class StreamSheet;    // StreamSheet.hpp: r(m), b(m) of an axisymmetric stream sheet (not in v1)

enum class Region : int { Inlet = 0, Passage = 1, Exit = 2 };   // v1 grid.INLET, PASSAGE, EXIT
enum class Surface : int { Suction = 0, Pressure = 1 };

/// Where the first and the last blade station lie (Aungier 2003 Sec. 5.3, pp. 89-90, Fig. 5-5;
/// 2000 Sec. 13.4, pp. 247-249, Fig. 13-3).  Not in v1: an EXPERIMENT option of v2, so the user
/// can test Aungier's statement that the "touching" layout can cause local numerical instability.
///  * Inside (default, v1, Aungier's recommendation): both stations inside the passage, one node
///    on each surface; the side lines meet the blade at the chord-end LE / TE points.
///  * Touching (Fig. 5-5 "outside", Fig. 13-3 "External Leading-Edge Grid Point"): the first
///    blade station passes through the nose tip (the blade's least m), the last through the tail
///    end (its greatest m); both wall nodes of such a station are that one point (the upper node
///    one pitch up: zero blade thickness), and the side lines end there.  Every station strictly
///    between is a normal blade station.  The grid alone carries the layout to the solver, whose
///    modules are v1's, unchanged: see Grid::wall_d and WallIdentity at the touching stations.
enum class EdgeNodes : int { Inside = 0, Touching = 1 };
[[nodiscard]] inline const char* edge_nodes_name(EdgeNodes e) noexcept
{
    return e == EdgeNodes::Touching ? "touching" : "inside";
}

/// v1 PassageBoundaries: the side boundaries of one passage, far inlet to far exit.
struct PassageBoundaries
{
    double chord = 1.0;                 ///< active.chord, the unit of the extension lengths
    double pitch = 1.0;                 ///< active.pitch
    Point periodic_translation{0.0, 1.0};  ///< (0, pitch)
    double inlet_angle = 0.0;           ///< far-inlet flow angle = upstream_beta[0]
    double exit_angle = 0.0;            ///< far-exit flow angle = downstream_beta.back()
    double m_in = 0.0, m_le = 0.0, m_te = 0.0, m_out = 0.0;
    double kappa_le = 0.0, kappa_te = 0.0;  ///< blade angles at which the side lines meet the blade
    std::vector<Point> lower_blade_boundary;  ///< suction surface samples, LE -> TE
    std::vector<Point> upper_blade_boundary;  ///< pressure surface samples + periodic_translation
    std::shared_ptr<const SurfaceCurves> blade_curves;  ///< the C2 curves (v2 is C2 only: never null)
    std::vector<double> upstream_m, upstream_beta;      ///< m_in .. m_le; inlet_angle .. kappa_le
    std::vector<Point> upstream_lower_boundary, upstream_upper_boundary;
    std::vector<double> downstream_m, downstream_beta;  ///< m_te .. m_out; kappa_te .. exit_angle
    std::vector<Point> downstream_lower_boundary, downstream_upper_boundary;
    /// Not in v1: the axisymmetric stream sheet the passage lies on (y = r_ref theta; the side lines
    /// integrate r d theta / dm = tan(beta), StreamSheet.hpp).  Null: the planar cascade (v1).
    std::shared_ptr<const StreamSheet> sheet;
};

/// v1 WallSide: one wall's nodes at the blade stations, in station order.
struct WallSide
{
    std::vector<Surface> surface;   ///< the named surface each node lies on
    std::vector<double> parameter;  ///< 0 at the LE point, 1 at the TE point (C2: u itself)
    std::vector<double> u;          ///< parameter of the node on its surface's curve
    std::vector<double> contour;    ///< tau along the wall's contour from the nose tip
    std::vector<int> hidden;        ///< further crossings of the station line behind the node
};

/// v1 WallIdentity.
struct WallIdentity
{
    std::vector<int> stations;      ///< first_blade_i .. last_blade_i
    WallSide lower, upper;
};

/// v1 Grid (frozen dataclass + its properties).
struct Grid
{
    PassageBoundaries boundaries;
    std::vector<double> m;          ///< n_m uniform stations, m[0] = m_in, m.back() = m_out
    std::vector<double> eta;        ///< n_eta uniform, eta[0] = 0, eta.back() = 1
    Field node_m, node_y;           ///< the nodes (v1 nodes[..., 0] and nodes[..., 1])
    int first_blade_i = 0, last_blade_i = 0;
    WallIdentity wall_identity;
    /// d^k X / du^k (k = 0, 1, 2) of the authoritative curve at each wall node of the blade
    /// stations (v1 Grid.wall_curve_derivative(k)); [k].lower[b] for station first_blade_i + b.
    /// k = 0: the wall nodes themselves (the upper wall on the +pitch neighbour).
    struct WallDerivative
    {
        std::vector<Point> lower, upper;
    };
    ///
    /// EdgeNodes::Touching, at the two touching stations only (entries 0 and last_blade_i -
    /// first_blade_i): k = 0 is the touching point itself (upper: + periodic_translation), but
    /// k = 1 and k = 2 are NOT the curve's derivatives (the curve's tangent is vertical at the nose
    /// tip and a round tail end: m is no coordinate along the wall there).  They are the mesh's
    /// assigned edge direction, the side line at its end parametrized by arc length s: X_s =
    /// (cos kappa, sin kappa), X_ss = (d beta/ds) (-sin kappa, cos kappa), d beta/ds = cos(kappa)
    /// d beta/dm of the side line's linear-angle law, d beta/dm = (far angle - kappa) / (m_far -
    /// m_junction) (kappa_le, inlet_angle, m_in at the LE station; kappa_te, exit_angle, m_out at
    /// the TE station).  The unchanged metrics then give tan(beta0) = tan(beta1) = y_s / m_s =
    /// tan(kappa) and d tan(beta)/dm = (y_ss m_s - y_s m_ss) / m_s^3 = sec^2(kappa) d beta/dm there
    /// (Aungier 2003 p. 86: the side boundary turns to the blade angle at the blade).
    WallDerivative wall_d[3];
    /// The layout of the first and last blade stations (not in v1; see EdgeNodes).  Touching:
    /// stations first_blade_i and last_blade_i are the touching stations, the rest as Inside.
    EdgeNodes edge_nodes = EdgeNodes::Inside;
    /// Touching only: the curve the nose tip ([0]) and the tail end ([1]) actually lie on, and their
    /// u there.  The wall identity LABELS each touching node with its wall's surface on the next
    /// station (lower suction, upper pressure), so the unchanged Kutta condition (Eq. 5-34 at the
    /// last blade station: j = 0 against j = n_eta - 1) and post-processing read one suction and
    /// one pressure node; its u is this u, on this curve.
    std::array<Surface, 2> touching_surface{Surface::Suction, Surface::Suction};
    std::array<double, 2> touching_u{0.0, 0.0};

    [[nodiscard]] int n_m() const noexcept { return static_cast<int>(m.size()); }
    [[nodiscard]] int n_eta() const noexcept { return static_cast<int>(eta.size()); }
    [[nodiscard]] double dm() const { return (m.back() - m.front()) / static_cast<double>(m.size() - 1); }
    [[nodiscard]] double deta() const { return (eta.back() - eta.front()) / static_cast<double>(eta.size() - 1); }
    [[nodiscard]] int last_inlet_i() const noexcept { return first_blade_i - 1; }
    [[nodiscard]] int first_exit_i() const noexcept { return last_blade_i + 1; }
    [[nodiscard]] Region region(int i) const noexcept
    {
        return i < first_blade_i ? Region::Inlet : (i > last_blade_i ? Region::Exit : Region::Passage);
    }
    [[nodiscard]] bool has_wall_curves() const noexcept { return true; }  // v2: C2 only
    [[nodiscard]] double m_le() const noexcept { return boundaries.m_le; }
    [[nodiscard]] double m_te() const noexcept { return boundaries.m_te; }
    [[nodiscard]] double m_in() const noexcept { return boundaries.m_in; }
    [[nodiscard]] double m_out() const noexcept { return boundaries.m_out; }
    [[nodiscard]] double chord() const noexcept { return boundaries.chord; }
    [[nodiscard]] double pitch() const noexcept { return boundaries.pitch; }
    [[nodiscard]] double inlet_angle() const noexcept { return boundaries.inlet_angle; }
    [[nodiscard]] double exit_angle() const noexcept { return boundaries.exit_angle; }
    [[nodiscard]] double m_first_blade() const { return m[static_cast<std::size_t>(first_blade_i)]; }
    [[nodiscard]] double m_last_blade() const { return m[static_cast<std::size_t>(last_blade_i)]; }
    [[nodiscard]] double actual_le_offset() const { return m_first_blade() - m_le(); }
    [[nodiscard]] double actual_te_offset() const { return m_last_blade() - m_te(); }
    [[nodiscard]] bool touching() const noexcept { return edge_nodes == EdgeNodes::Touching; }
    /// Station i is a touching station (Touching layout, i = first_blade_i or last_blade_i).
    [[nodiscard]] bool touching_station(int i) const noexcept
    {
        return touching() && (i == first_blade_i || i == last_blade_i);
    }
};

}  // namespace BladeGen::Turbolenza
