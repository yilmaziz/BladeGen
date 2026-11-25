/// @file turbolenza/Grid.hpp
/// @brief Grid generation for one blade passage: port of v1 grid.py (C2 path; sections 1-9 of
///        its docstring), the computational (m, eta) grid with its inlet and exit extensions
///        (Aungier 2003 Secs. 5.2, 5.3, Figs. 5-3, 5-5).  The data is in GridData.hpp.
///
/// Kept from v1 unchanged: constant dm and deta (Aungier p. 87); the side lines whose angle beta
/// (not tan beta) runs linearly in m from the blade angle at the LE / TE to the flow angle at the
/// far station (Aungier Sec. 5.3, p. 86; 2000 Sec. 13.1, p. 237), so the grid lines at the far
/// stations lie along the flow; the exit extension rebuilt for every discharge angle of the Kutta
/// loop; the blade-passage topology read from the station cuts; the five geometry-driven standard
/// levels.  Aungier's rule for the edges (2003 Sec. 5.3, pp. 89-90, Fig. 5-5; 2000 Sec. 13.4,
/// pp. 247-249, Figs. 13-3, 13-4): the first and last grid points on the blade surfaces lie inside
/// the blade passage, one node on each surface, never at the LE or TE itself.  v2 enforces it on
/// every grid it hands to the solver (the standard levels by construction; an explicit
/// n_passage x n_eta grid is refused when its first or last blade station does not cut exactly
/// one suction and one pressure surface: v1's require_ss_ps).
///
/// EXPERIMENT (not in v1): GridSettings::edge_nodes = EdgeNodes::Touching builds Aungier's other
/// layout (2003 Fig. 5-5 "outside"; 2000 Fig. 13-3), which he reports can cause local numerical
/// instability: the first blade station through the nose tip m_tip, the last through the tail end
/// m_tail, both wall nodes of such a station on that one point (zero blade thickness), the side
/// lines ending there (anchored at the nose tip / tail end instead of the chord-end LE / TE
/// points; the idea of v1's scripts/sideline_anchoring_study.py "tip" anchor).  The lattice: dm =
/// (m_tail - m_tip) / n_cells (standard levels: the same target dm per level, candidate band,
/// CELL_ASPECT eta rule and whole-cell extensions; explicit: n_cells = n_passage), stations at
/// m_tip + k dm.  The stations next to the touching ones (first_blade_i + 1, last_blade_i - 1)
/// must each cut one suction and one pressure surface, or the grid is refused.  The Grid
/// invariants that forbid a touching station (two distinct wall nodes; the TOUCH rule that makes
/// a line through a single blade point periodic) are relaxed for those two stations only.
/// Only the grid knows the layout; the solver modules are v1's, unchanged.  So the grid hands them
/// a touching station in their own terms (GridData.hpp): Grid::wall_d[1], [2] there are the
/// mesh's assigned edge direction (the side line's end at kappa_le / kappa_te), not the curve's
/// vertical tangent, and the wall identity labels the lower node suction, the upper pressure.
#pragma once

#include "turbolenza/Cascade.hpp"
#include "turbolenza/GridData.hpp"

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 config.GridSettings: explicit (n_passage, n_eta; level empty) or a standard level
/// (level one of kGridLevels; n_passage, n_eta empty).  Extensions in true blade chords.
struct GridSettings
{
    std::optional<int> n_passage;
    std::optional<int> n_eta;
    double upstream_extension_chords = 1.0;
    double downstream_extension_chords = 1.0;
    std::optional<std::string> level;
    EdgeNodes edge_nodes = EdgeNodes::Inside;  ///< not in v1: the EXPERIMENT layout (see the file comment)
};

/// v1 GRID_LEVELS.
inline const std::array<const char*, 5> kGridLevels = {"very_coarse", "coarse", "medium", "fine", "very_fine"};

/// v1 GridError: the grid cannot be built, or a Grid would break the blade-passage topology.
class GridError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 GridTopologyError: an explicit grid whose first or last blade station does not cut exactly
/// one suction and one pressure surface (Aungier's edge rule above).
class GridTopologyError : public GridError
{
public:
    GridTopologyError(const std::string& what, std::string edge_, std::optional<int> n_passage_,
                      std::optional<int> nearest_lower_, std::optional<int> nearest_higher_, std::optional<int> finest_,
                      std::optional<int> ceiling_)
        : GridError(what), edge(std::move(edge_)), n_passage(n_passage_), nearest_lower(nearest_lower_),
          nearest_higher(nearest_higher_), finest(finest_), ceiling(ceiling_)
    {
    }
    std::string edge;  ///< "LE", "TE" or "LE and TE"
    std::optional<int> n_passage, nearest_lower, nearest_higher, finest, ceiling;
};

/// v1 rejection reasons of a standard-level candidate (REJECT_*).
inline constexpr const char* kRejectLE = "invalid LE cut";
inline constexpr const char* kRejectTE = "invalid TE cut";
inline constexpr const char* kRejectBranch = "wrong branch identity";
inline constexpr const char* kRejectExtension = "extension outside tolerance";
inline constexpr const char* kRejectEta = "eta-resolution issue";
inline constexpr const char* kRejectOther = "other topology failure";

/// v1 CandidateRecord: one candidate spacing tried for a standard level.
struct CandidateRecord
{
    int n_cells = 0;
    double dm_over_chord = 0.0;
    int rank = 0;          ///< 1: nearest the target
    bool accepted = false;
    std::string reason;    ///< one of kReject*, empty when accepted
    std::string detail;
};

/// v1 GridGenerationError: no candidate of a standard level passes (the reasons are counted).
class GridGenerationError : public GridError
{
public:
    using GridError::GridError;
    GridGenerationError(const std::string& what, std::optional<std::string> level_,
                        std::vector<CandidateRecord> candidates_, std::vector<std::pair<std::string, int>> reasons_)
        : GridError(what), level(std::move(level_)), candidates(std::move(candidates_)), reasons(std::move(reasons_))
    {
    }
    std::optional<std::string> level;                  ///< the level asked for; empty: the geometry itself
    std::vector<CandidateRecord> candidates;           ///< every candidate tried
    std::vector<std::pair<std::string, int>> reasons;  ///< rejection reason -> count, first-seen order
};

/// v1 StandardGridInfo (the numbers of a standard grid; lengths in chords).  v1's geometry_id
/// (a sha256 fingerprint of the cascade, for messages) is not ported.
struct StandardGridInfo
{
    std::string level;
    double target_dm_over_chord = 0.0;
    int n_blade_cells = 0;
    double dm = 0.0, dm_over_chord = 0.0, deta = 0.0;
    int n_m = 0, n_eta = 0, total_nodes = 0;
    int first_blade_i = 0, last_blade_i = 0, n_blade_stations = 0;
    double first_ss_x_c = 0.0, first_ps_x_c = 0.0, last_ss_x_c = 0.0, last_ps_x_c = 0.0;
    double le_gap_over_chord = 0.0, te_gap_over_chord = 0.0;
    double m_first_minus_m_le_over_chord = 0.0, m_te_minus_m_last_over_chord = 0.0;
    double upstream_extension_chords = 0.0, downstream_extension_chords = 0.0;
    double cell_aspect = 0.0;
    std::vector<std::string> checks;
    // added with the port (the v1 fields of the same names)
    double m_first = 0.0, m_last = 0.0;
    /// (m, y) on the reference blade (the pressure-side wall node moved back one pitch)
    Point first_ss_wall_point, first_ps_wall_point, last_ss_wall_point, last_ps_wall_point;
    std::string admissibility;
};

struct StandardGrid
{
    Grid grid;
    StandardGridInfo info;
    // added with the port: v1 CandidateStatistics
    std::vector<CandidateRecord> candidates;  ///< every candidate tried, in search order
    double target_dm_over_chord = 0.0;
    std::array<double, 2> band{0.0, 0.0};     ///< dm / chord range searched
};

// ── v1 grid.py constants ─────────────────────────────────────────────────────
inline constexpr double kTouch = 1e-12;                    ///< TOUCH, chords
inline constexpr double kMediumDmOverChord = 0.020;        ///< MEDIUM_DM_OVER_CHORD
inline constexpr double kCandidateBand = 0.10;             ///< CANDIDATE_BAND
inline constexpr double kTransitionClearance = 0.10;       ///< TRANSITION_CLEARANCE, of dm
inline constexpr double kCellAspect = 1.8;                 ///< CELL_ASPECT = pitch deta / dm
inline constexpr int kMinEtaCells = 6;                     ///< MIN_ETA_CELLS
inline constexpr double kNominalExtensionChords = 1.0;     ///< NOMINAL_EXTENSION_CHORDS
inline constexpr double kExtensionToleranceChords = 0.05;  ///< EXTENSION_TOLERANCE_CHORDS
inline constexpr double kUniformSpacingTolerance = 1e-13;  ///< UNIFORM_SPACING_TOLERANCE, chords

/// v1 LEVEL_TARGET_DM_OVER_CHORD: 0.020 sqrt(2)^(2 - k) for level k of kGridLevels;
/// std::invalid_argument for an unknown level.
[[nodiscard]] double level_target_dm_over_chord(const std::string& level);

/// v1 PassageTransitions: the blade ends along m (grid.py section 9), in the active unit.
struct PassageTransitions
{
    double m_tip = 0.0, m_le = 0.0, m_start = 0.0, m_end = 0.0, m_te = 0.0, m_tail = 0.0, chord = 1.0;
    struct Interval
    {
        double lo = 0.0, hi = 0.0;
        int suction = 0, pressure = 0;
    };
    std::vector<Interval> coverage;  ///< every m interval between branch ends, with its crossing counts

    [[nodiscard]] double length() const noexcept { return m_end - m_start; }
    [[nodiscard]] double le_window() const noexcept { return m_start - m_tip; }
    [[nodiscard]] double te_window() const noexcept { return m_tail - m_end; }
    [[nodiscard]] double le_face() const noexcept { return 0.5 * (m_tip + m_start); }
    [[nodiscard]] double te_face() const noexcept { return 0.5 * (m_end + m_tail); }
    [[nodiscard]] double lattice_length() const noexcept { return te_face() - le_face(); }
    [[nodiscard]] bool valid_cut(double m) const noexcept { return m_start < m && m < m_end; }
};

/// v1 passage_transitions.  GridGenerationError when no SS + PS-valid interval holds the middle of
/// the chord or the passage is closed.
[[nodiscard]] PassageTransitions passage_transitions(const Cascade& cascade);

/// v1 EdgeStationTopology: how the first ("LE") or the last ("TE") blade station cuts the blade.
struct EdgeStationTopology
{
    std::string edge;
    int station = 0;
    double m = 0.0;
    int crossings = 0, suction = 0, pressure = 0;
    double branch_point = 0.0;  ///< m_le or m_te
    double extreme = 0.0;       ///< m of the nose tip or the tail end
    bool on_branch_point = false;
    [[nodiscard]] bool valid() const noexcept
    {
        return crossings == 2 && suction == 1 && pressure == 1 && !on_branch_point;
    }
};

/// v1 wall_station_topology: (LE, TE) of a built grid.
[[nodiscard]] std::array<EdgeStationTopology, 2> wall_station_topology(const Grid& grid);

/// v1 make_passage_boundaries.  exit_angle empty: kappa_te.
[[nodiscard]] PassageBoundaries make_passage_boundaries(const Cascade& cascade, const GridSettings& settings,
                                                        double inlet_angle, std::optional<double> exit_angle);

/// v1 make_grid: explicit or standard-level grid (settings.level).  require_ss_ps: refuse an
/// explicit grid whose edge stations break Aungier's rule (GridTopologyError).  The solver always
/// builds with require_ss_ps = true.
[[nodiscard]] Grid make_grid(const Cascade& cascade, const GridSettings& settings, double inlet_angle,
                             std::optional<double> exit_angle, bool require_ss_ps = false);

/// v1 build_grid: the standard grid of one level (with its info).
[[nodiscard]] StandardGrid build_grid(const Cascade& cascade, const std::string& level, double inlet_angle,
                                      std::optional<double> exit_angle);

/// v1 build_grid with the nominal (upstream, downstream) extension lengths in chords.
[[nodiscard]] StandardGrid build_grid(const Cascade& cascade, const std::string& level, double inlet_angle,
                                      std::optional<double> exit_angle, std::array<double, 2> extension_chords);
/// The same in the given edge-node layout (EdgeNodes::Inside: exactly the overload above).
[[nodiscard]] StandardGrid build_grid(const Cascade& cascade, const std::string& level, double inlet_angle,
                                      std::optional<double> exit_angle, std::array<double, 2> extension_chords,
                                      EdgeNodes edge_nodes);

/// v1 rebuild_exit: the same grid with its exit extension turned to exit_angle; only the exit
/// stations' nodes move.
[[nodiscard]] Grid rebuild_exit(const Grid& grid, double exit_angle);
/// v1 rebuild_exit of PassageBoundaries: exit_angle, downstream_beta and the downstream lines.
[[nodiscard]] PassageBoundaries rebuild_exit(const PassageBoundaries& boundaries, double exit_angle);

/// v1 validate_uniform_grid_topology: GridTopologyError unless the first and last blade stations
/// each cut one suction and one pressure surface.  A Touching grid (not in v1): unless its first
/// and last blade stations lie on the nose tip and the tail end and the stations next to them
/// (first_blade_i + 1, last_blade_i - 1) each cut one suction and one pressure surface.
void validate_uniform_grid_topology(const Grid& grid);
/// The same with v1's keywords: the n_passage that built the grid (for the nearest admissible
/// values in the message and the error) and the stagger [rad] that names the cascade.
void validate_uniform_grid_topology(const Grid& grid, std::optional<int> n_passage, std::optional<double> stagger);

/// v1 find_valid_n_passage: the admissible n_passage among n_values, ascending (no grid is built);
/// in the layout settings.edge_nodes.
[[nodiscard]] std::vector<int> find_valid_n_passage(const Cascade& cascade, const GridSettings& settings,
                                                    double inlet_angle, const std::vector<int>& n_values);
/// v1 n_passage_ceiling: the first n_passage from which every finer one fails; empty when the blade
/// has neither an LE nor a TE window.  Touching: the first n_passage whose dm = (m_tail - m_tip) /
/// n_passage does not exceed the wider window (station first_blade_i + 1 or last_blade_i - 1 then
/// lies in it).
[[nodiscard]] std::optional<int> n_passage_ceiling(const Cascade& cascade, const GridSettings& settings,
                                                   double inlet_angle);

/// v1 validate_standard_grid: the invariants of a standard grid; the names of the checks passed, or
/// GridError naming the first that fails.  v1's last check, "metrics with positive S", calls
/// metrics.py and is not run by the grid module (the metrics come after the grid in v2).
[[nodiscard]] std::vector<std::string> validate_standard_grid(const Grid& grid);

/// Not in v1: the passage width S at the far inlet station, r_in (theta1 - theta0) = sigma(m_in)
/// pitch on a stream sheet (Eq. 5-9); the pitch of the planar cascade (v1's S1).
[[nodiscard]] double inlet_passage_width(const Grid& grid);
/// Not in v1: the radius of the far inlet station [chords]; empty for the planar cascade.
[[nodiscard]] std::optional<double> inlet_radius(const Grid& grid);

}  // namespace BladeGen::Turbolenza
