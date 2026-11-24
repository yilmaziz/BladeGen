/// @file turbolenza/Metrics.hpp
/// @brief The geometric fields of Aungier's transformed blade-to-blade equations on a grid.
///
/// Port of v1 metrics.py (TURBOlenza v1 commit 67599a6; Aungier 2003 Secs. 5.2, 5.3; Eqs. 5-1,
/// 5-2, 5-3, 5-9).  Pure geometry, computed once per grid (again after rebuild_exit).
///
///  * planar cascade: y = r theta, so r d theta / dm = dy/dm, no radius anywhere;
///  * eta = (y - y0(m)) / S(m), S = y1 - y0 the tangential passage width [Eq. 5-9]; det J = S;
///  * tan(beta) = (1 - eta) tan(beta0) + eta tan(beta1) [Eq. 5-3], tan(beta0/1) = dy0/dm, dy1/dm;
///  * d/dm of the boundaries in three runs (inlet side lines 0 .. f - 1, blade walls f .. l, exit
///    side lines l + 1 .. n_m - 1), second order: central inside, (-3, 4, -1) / (2 dm) at a run's
///    first station and (3, -4, 1) / (2 dm) at its last; no stencil crosses a transition; a run of
///    fewer than three stations is a MetricsError;
///  * C2 blade walls (v2 is C2 only): the authoritative curve's own derivatives at the wall nodes
///    (Grid::wall_d), tan(beta) = y_u / m_u and d tan(beta)/dm = (y_uu m_u - y_u m_uu) / m_u^3,
///    refused where |m_u| <= sqrt(eps) |X_u|;
///  * cos(beta) <= sqrt(eps) anywhere is a MetricsError (1 / cos^2 singular); nothing is clipped.
///
/// Not in v1: a grid on an axisymmetric stream sheet (PassageBoundaries::sheet; StreamSheet.hpp).
/// Its nodes are (m, y = r_ref theta), so with sigma(m) = r(m) / r_ref at each station:
///  * S = r (theta1 - theta0) = sigma (y1 - y0)                      [Eq. 5-9]
///  * tan(beta0/1) = r d theta0/1 / dm = sigma dy0/1 / dm            [Eq. 5-3]
///    (the blade walls: sigma y_u / m_u; d tan / dm = sigma' y_u / m_u + sigma (y_uu m_u - y_u m_uu) / m_u^3,
///    sigma' = sin(phi) / r_ref; the side lines: the run differences of tan as v1)
///  * dS/dm = sin(phi) (theta1 - theta0) + tan(beta1) - tan(beta0)
///  * the station values r, b / b(m_in) and sin(phi) = dr/dm [Eq. 5-15] for the coefficients, the
///    flow recovery, the gas state and the Kutta correction.
/// The checks of the grid contract (S = pitch on the extensions, the eta map) stay on y1 - y0.
/// Without a sheet nothing of this runs: the fields are v1's, bit for bit.
#pragma once

#include "turbolenza/Field.hpp"
#include "turbolenza/GridData.hpp"

#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 EPSILON, ROUND_OFF (relative tolerance of the contract checks), SINGULAR_COS.
inline constexpr double kEpsilon = std::numeric_limits<double>::epsilon();
inline constexpr double kRoundOff = 64.0 * kEpsilon;
[[nodiscard]] double singular_cos() noexcept;

/// v1 MetricsError.
class MetricsError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 GridMetrics.  Station fields have n_m entries, node fields are (n_m, n_eta).
struct GridMetrics
{
    std::shared_ptr<const Grid> grid;   ///< the grid they belong to (by reference, as v1)
    std::vector<double> S;              ///< y1 - y0 [Eq. 5-9]
    std::vector<double> dS_dm;          ///< tan_beta1 - tan_beta0
    std::vector<double> tan_beta0, tan_beta1;
    Field tan_beta, beta, sin_beta, cos_beta;
    std::vector<double> dtan_beta0_dm, dtan_beta1_dm;
    Field dtan_beta_dm;
    /// Not in v1: the stream sheet at the stations (StreamSheet.hpp).  The planar cascade: 1, 1, 1
    /// and 0, the neutral factors (so every product with them is v1's, bit for bit).
    std::vector<double> radius_ratio;   ///< sigma = r / r_ref
    std::vector<double> radius;         ///< r [chords] (planar: 1, no radius)
    std::vector<double> thickness;      ///< b / b(m_in): Aungier's b, 1 at the far inlet
    std::vector<double> sin_phi;        ///< dr/dm [Eq. 5-15]
    bool on_sheet = false;              ///< the grid lies on a stream sheet
    /// The station values; a GridMetrics filled by hand (v1 reference data) without them: planar.
    [[nodiscard]] double sigma_at(int i) const noexcept { return at(radius_ratio, i, 1.0); }
    [[nodiscard]] double radius_at(int i) const noexcept { return at(radius, i, 1.0); }
    [[nodiscard]] double thickness_at(int i) const noexcept { return at(thickness, i, 1.0); }
    [[nodiscard]] double sin_phi_at(int i) const noexcept { return at(sin_phi, i, 0.0); }

    [[nodiscard]] const std::vector<double>& jacobian_det() const noexcept { return S; }
    [[nodiscard]] const std::vector<double>& dtan_beta_deta() const noexcept { return dS_dm; }
    [[nodiscard]] static double at(const std::vector<double>& v, int i, double neutral) noexcept
    {
        return v.empty() ? neutral : v[static_cast<std::size_t>(i)];
    }
};

/// v1 build_metrics: checks the grid contract (shapes, finite, uniform m and eta, nodes on their
/// stations and on the linear eta map, S = pitch on the extensions, three stations per run, the
/// C2 wall derivatives at the wall nodes) and computes the fields.  MetricsError otherwise.
[[nodiscard]] GridMetrics build_metrics(std::shared_ptr<const Grid> grid);

/// v1 _runs: the three m regions (start, stop) with stop - start >= 3 each, else MetricsError.
struct MRun
{
    const char* name;
    int start, stop;
};
[[nodiscard]] std::vector<MRun> m_runs(int n_m, int first_blade_i, int last_blade_i);

/// v1 _d_dm: d/dm of station values (or of node values, column by column: `stride` = n_eta)
/// inside each run, second order.  `values` has n_m * stride entries, row major.
[[nodiscard]] std::vector<double> d_dm_runs(const std::vector<double>& values, int stride, double dm,
                                            const std::vector<MRun>& runs);

}  // namespace BladeGen::Turbolenza
