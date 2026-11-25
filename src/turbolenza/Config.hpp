/// @file turbolenza/Config.hpp
/// @brief Solver inputs and numerical settings as plain structs: nothing here computes anything.
///
/// Port of v1 config.py (TURBOlenza v1 commit 67599a6).  The grid settings are Grid.hpp's
/// GridSettings; the gas is Eos.hpp's PerfectGas.  Production paths only: the v1 study options
/// of IterationSettings (kutta_eval, kutta_eval_offset_chords, kutta_eval_interpolation,
/// transition_form, kutta_method) are not ported; v2 always runs their defaults ("station",
/// "expanded", "current").
///
/// Frames and units (v1): the inlet total pressure and temperature are ABSOLUTE-frame totals;
/// the inlet Mach number and the flow angles are relative to the blade row.  An angle a person
/// types is in degrees and its name ends in _deg; inside, radians.
///
/// Tolerances and pass limits have no default in v1 ("no number chosen for one cascade carries
/// over unnoticed"): here they start as NaN / 0, which the solver refuses or never meets.
#pragma once

#include "turbolenza/Eos.hpp"
#include "turbolenza/Grid.hpp"

#include <limits>
#include <optional>

namespace BladeGen::Turbolenza {

/// v1 KUTTA_RELAXATION: lambda_K of the project baseline (the undamped Aungier correction).
inline constexpr double kKuttaRelaxation = 1.0;

/// v1 FlowConditions: the operating point.
struct FlowConditions
{
    std::optional<double> inlet_mach;           ///< far-inlet relative Mach number; empty: mass-flow override only
    double inlet_angle_deg = 0.0;               ///< relative flow angle at the far inlet
    double exit_angle_guess_deg = 0.0;          ///< first discharge angle (the fixed one at a fixed angle)
    double inlet_absolute_total_pressure = 101325.0;
    double inlet_absolute_total_temperature = 288.15;
    double angular_speed = 0.0;                 ///< omega; only 0 is solved
    std::optional<double> stream_surface_radius;            ///< rotating rows only
    std::optional<double> inlet_absolute_tangential_velocity;  ///< rotating rows only
};

/// v1 IterationSettings (production fields).
struct IterationSettings
{
    double density_tolerance = std::numeric_limits<double>::quiet_NaN();
    int max_density_iterations = 0;
    double kutta_tolerance = std::numeric_limits<double>::quiet_NaN();
    int max_kutta_iterations = 0;
    double density_relaxation = kDensityRelaxation;     ///< lambda_rho, 0 < lambda_rho <= 1
    std::optional<double> velocity_tolerance;           ///< empty: density_tolerance
    double kutta_relaxation = kKuttaRelaxation;         ///< lambda_K, 0 < lambda_K <= 1
    bool stop_at_sonic = true;
};

/// v1 SolverConfig: everything the solver needs besides the geometry.
struct SolverConfig
{
    GridSettings grid;
    FlowConditions flow;
    PerfectGas gas;
    IterationSettings iteration;
};

}  // namespace BladeGen::Turbolenza
