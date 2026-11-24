/// @file turbolenza/Postprocess.hpp
/// @brief Results derived from a finished Solution: blade surfaces, pressure coefficients,
///        discharge angle, summary.  Nothing here solves or feeds back into a solve.
///
/// Port of v1 postprocess.py (TURBOlenza v1 commit 67599a6):
///  * blade surfaces over the blade stations only; which wall is which surface from the wall
///    identity (a wall whose surface changes along the blade is a PostprocessError); Mach and
///    static pressure of the reported evaluation at the wall nodes;
///  * x/c = (P - P_LE) . (P_TE - P_LE) / |P_TE - P_LE|^2 (the blade-frame chordwise coordinate; the
///    upper wall moved back by periodic_translation); m/c = (m - m_le) / chord beside it;
///  * pressure coefficients on the specified inlet state: Dunavant et al. 1955 S = (Pt1 - p) / q1
///    (default) and Cp = (p - p1) / q1, q1 = rho1 W1^2 / 2;
///  * discharge angle = degrees(grid.exit_angle); the far-exit velocity angle
///    atan2(int W_theta deta, int W_m deta) as a check (trapezoidal rule).
/// read_surface_points (experimental CSV) is not ported (the window reads data files).
#pragma once

#include "turbolenza/Solver.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 PostprocessError.
class PostprocessError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

enum class PressureConvention { Dunavant, InletStatic };  ///< v1 DUNAVANT, INLET_STATIC
[[nodiscard]] const char* pressure_label(PressureConvention c) noexcept;

/// v1 BladeSurfaces.
struct BladeSurfaces
{
    std::vector<int> stations;
    std::vector<double> m, m_c;
    int ss_j = 0, ps_j = 0;
    std::vector<double> x_c_ss, x_c_ps, mach_ss, mach_ps, pressure_ss, pressure_ps;
};

/// v1 SurfaceCoefficient.
struct SurfaceCoefficient
{
    PressureConvention convention = PressureConvention::Dunavant;
    std::string label;
    std::vector<double> x_c_ss, values_ss, x_c_ps, values_ps;
};

/// v1 SolutionSummary.
struct SolutionSummary
{
    std::string limitation;
    bool converged = false;
    int iterations = 0;
    double mass_flow_per_depth = 0.0;
    std::string mass_flow_source;
    std::optional<double> inlet_mach;
    double inlet_angle_deg = 0.0, discharge_angle_deg = 0.0, exit_velocity_angle_deg = 0.0;
    double exit_angle_mismatch_deg = 0.0, turning_deg = 0.0;
    double surface_mach_max = 0.0, surface_mach_min = 0.0;
    double ss_peak_mach = 0.0, ss_peak_x_c = 0.0, ps_peak_mach = 0.0, ps_peak_x_c = 0.0;
    double kutta_residual = 0.0, normalized_kutta_residual = 0.0;
    double mach_max = 0.0;
    std::array<int, 2> mach_max_node{0, 0};
    double density_min = 0.0, density_max = 0.0, pressure_min = 0.0, pressure_max = 0.0;
    double temperature_min = 0.0, temperature_max = 0.0;
    double mass_flow_error_max = 0.0, mass_flow_error_blade_median = 0.0, linear_residual = 0.0;
    double density_change = 0.0, velocity_change = 0.0;
};

[[nodiscard]] BladeSurfaces blade_surfaces(const Solution& solution);
/// x/c of points of the reference blade (cascade frame) between the LE and TE points.
[[nodiscard]] std::vector<double> chordwise_coordinate(const std::vector<Point>& points, const Point& leading_edge,
                                                       const Point& trailing_edge);
[[nodiscard]] SurfaceCoefficient surface_pressure_coefficient(const Solution& solution,
                                                              PressureConvention convention = PressureConvention::Dunavant);
/// The coefficient of `convention` for static pressures p on the inlet state (PostprocessError
/// without an inlet state or with q1 <= 0).
[[nodiscard]] std::vector<double> pressure_coefficient_values(const std::vector<double>& p,
                                                              const std::optional<InletState>& inlet,
                                                              PressureConvention convention);
[[nodiscard]] double discharge_angle_deg(const Solution& solution);
[[nodiscard]] double exit_velocity_angle_deg(const Solution& solution);
[[nodiscard]] SolutionSummary summarize(const Solution& solution);
[[nodiscard]] std::string format_summary(const SolutionSummary& summary);

}  // namespace BladeGen::Turbolenza
