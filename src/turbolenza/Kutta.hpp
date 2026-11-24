/// @file turbolenza/Kutta.hpp
/// @brief The trailing-edge Kutta condition and the discharge-angle correction.
///
/// Port of v1 kutta.py, production rule (TURBOlenza v1 commit 67599a6; Aungier 2003 Eq. 5-34;
/// Aungier 2000 Sec. 13.3, Eqs. 13-57, 13-59 in this project's angle convention):
///  * where: i = last_blade_i, the SS and PS node read from the wall identity (KuttaError unless
///    one suction and one pressure node);
///  * residual: R_K = W_SS - W_PS, eps_K = |R_K| / (0.5 (|W_SS| + |W_PS|)) (0 when both speeds 0);
///  * correction: Delta(tan beta_out) = 2 (W0 - W1) dm / (S (cos beta0 + cos beta1) int_0^1 W_m deta),
///    W0, W1 the speeds on the eta = 0 / 1 walls at the TE station, S, beta0, beta1 there (cos from
///    the project's tan(beta) = dy/dm, 1 / sqrt(1 + tan^2)), the integral at the far exit station by
///    the trapezoidal rule; dimensionless (normalized S and dm);
///  * next angle: atan(tan(beta_out) + lambda_K Delta), 0 < lambda_K <= 1; a zero step returns the
///    angle exactly.
/// The v1 study options (kutta_eval "fixed_m", kutta_method "aungier_13_57*") are not ported.
#pragma once

#include "turbolenza/Flow.hpp"
#include "turbolenza/GridData.hpp"
#include "turbolenza/Metrics.hpp"

#include <stdexcept>

namespace BladeGen::Turbolenza {

/// v1 KuttaError.
class KuttaError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 TrailingEdgeStation.
struct TrailingEdgeStation
{
    int i = 0;
    double m = 0.0, m_offset = 0.0;
    int ss_j = 0, ps_j = 0;
};

/// v1 KuttaResidual.
struct KuttaResidual
{
    TrailingEdgeStation station;
    double W_ss = 0.0, W_ps = 0.0;
    double residual = 0.0;          ///< R_K = W_ss - W_ps
    double normalized = 0.0;        ///< eps_K
    double wall_difference = 0.0;   ///< W0 - W1
};

/// v1 ExitAngleCorrection (undamped).
struct ExitAngleCorrection
{
    double wall_difference = 0.0, delta_c_theta = 0.0, exit_meridional_integral = 0.0;
    double passage_width = 0.0, dm = 0.0, cos_beta0 = 0.0, cos_beta1 = 0.0;
    double delta_tan = 0.0;
};

[[nodiscard]] TrailingEdgeStation trailing_edge_station(const Grid& grid);
[[nodiscard]] KuttaResidual kutta_residual(const FlowField& flow, const Grid& grid);
[[nodiscard]] ExitAngleCorrection exit_angle_correction(const KuttaResidual& residual, const FlowField& flow,
                                                        const GridMetrics& metrics);
[[nodiscard]] double next_exit_angle(double exit_angle, const ExitAngleCorrection& correction, double relaxation);

}  // namespace BladeGen::Turbolenza
