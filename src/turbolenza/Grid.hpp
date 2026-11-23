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
#pragma once

#include "turbolenza/Cascade.hpp"
#include "turbolenza/GridData.hpp"

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
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

/// v1 GridGenerationError: no candidate of a standard level passes (the reasons are counted).
class GridGenerationError : public GridError
{
public:
    using GridError::GridError;
};

/// v1 StandardGridInfo (the numbers of a standard grid; lengths in chords).
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
};

struct StandardGrid
{
    Grid grid;
    StandardGridInfo info;
};

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

/// v1 rebuild_exit: the same grid with its exit extension turned to exit_angle; only the exit
/// stations' nodes move.
[[nodiscard]] Grid rebuild_exit(const Grid& grid, double exit_angle);

/// v1 validate_uniform_grid_topology: GridTopologyError unless the first and last blade stations
/// each cut one suction and one pressure surface.
void validate_uniform_grid_topology(const Grid& grid);

}  // namespace BladeGen::Turbolenza
