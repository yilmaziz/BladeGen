/// @file turbolenza/MeshStudy.hpp
/// @brief A mesh-independence study on the run layer: the same case solved on several standard
///        grid levels (kGridLevels), one after the other, and its files.  No numerics here: each
///        level is build_case + run_case of a copy of the base settings with only the grid level
///        changed.
#pragma once

#include "turbolenza/Case.hpp"

#include <atomic>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// One level of the study.
struct MeshStudyEntry
{
    std::string level;
    int n_m = 0, n_eta = 0, total_nodes = 0;     ///< 0: the case was not built
    double dm_over_chord = std::numeric_limits<double>::quiet_NaN();  ///< the standard grid's dm / chord; NaN: not built
    bool converged = false;
    int passes = 0;                              ///< RunResult::iterations
    double beta_out_deg = std::numeric_limits<double>::quiet_NaN();  ///< the reported solution's exit angle; NaN: no solution
    double mach_max = std::numeric_limits<double>::quiet_NaN();  ///< NaN: no pass
    double wall_time_s = 0.0;                    ///< build + solve
    std::string error;                           ///< empty: converged; else the build or solve error
    std::optional<Case> c;                       ///< empty: build_case failed
    std::optional<RunResult> result;             ///< empty: build_case failed
};

struct MeshStudyResult
{
    std::vector<MeshStudyEntry> levels;          ///< in the order run
    bool cancelled = false;                      ///< the remaining levels were skipped
};

/// Called from the thread that runs the study.
struct MeshStudyCallbacks
{
    /// the case of a level is built (its grid size known), its solve starts
    std::function<void(const std::string& level, int index, int n_m, int n_eta)> level_started;
    /// every pass of the level's solve (run_case's observer)
    std::function<void(const std::string& level, const PassRecord&)> pass;
    /// the level is done (converged, failed or cancelled)
    std::function<void(const MeshStudyEntry&)> level_finished;
};

/// For each level in order: a copy of `base` with grid_level = level (n_passage, n_eta cleared;
/// every other setting unchanged), build_case + run_case.  A level that fails (CaseError at build,
/// a solver error, not converged) is recorded with its error and the study goes on.  `cancel`
/// (may be null) stops the current solve after its pass and skips the remaining levels.
[[nodiscard]] MeshStudyResult run_mesh_study(const CaseSettings& base, const BladeSection& section, const std::vector<std::string>& levels,
                                             const MeshStudyCallbacks& callbacks = {}, const std::atomic<bool>* cancel = nullptr);

/// mesh_study.csv (level, n_m, n_eta, nodes, dm_over_chord, converged, passes, beta_out_deg,
/// mach_max, wall_s, error) in `directory` and export_run's files of every built level in
/// `directory`/<level>/.  An existing file is overwritten only when `overwrite`.  Returns the
/// files written, relative to `directory`.
std::vector<std::string> write_mesh_study(const MeshStudyResult& result, const std::string& directory, bool overwrite = false);

}  // namespace BladeGen::Turbolenza
