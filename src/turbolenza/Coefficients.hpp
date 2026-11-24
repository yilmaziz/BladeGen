/// @file turbolenza/Coefficients.hpp
/// @brief The coefficient fields A, B, C, D, E of Aungier's stream-function equation.
///
/// Port of v1 coefficients.py (TURBOlenza v1 commit 67599a6; Aungier 2003 Eqs. 5-20 to 5-26).
/// Expanded form (Eq. 5-21): A psi_eta,eta - 2 B psi_eta,m + C psi_m,m + D psi_eta + E psi_m = 0
/// (planar, non-rotating: 2 S omega sin(phi) = 0), with mu = mdot' / rho, mdot' = mdot / b the
/// mass flow through one passage per unit depth:
///     A = mu / (S cos^2 beta)  [5-22]    B = mu tan(beta)  [5-23]    C = mu S  [5-24]
///     D = dA/d eta - dB/dm     [5-25]    E = dC/dm - dB/d eta  [5-26]
/// D and E are second-order differences of the complete fields: in m one stencil over the whole
/// range (central, one-sided (-3, 4, -1) / (2 dm) at i = 0 and (3, -4, 1) / (2 dm) at n_m - 1),
/// crossing the transitions; in eta central inside, one-sided into the fluid at the walls of the
/// blade stations, and the periodic central difference (f[1] - f[n_eta - 2]) / (2 deta) at both
/// images of the periodic stations.
/// Lengths: the grid is chord-normalized; c_ref (metres) is the one physical scale, S = c_ref S*,
/// dm = c_ref dm*; c_ref empty = normalized mode (L_ref = 1 chord).  A grid whose chord is not 1
/// is refused (the chord would be applied twice).  The per-node maps run on the solve's threads.
#pragma once

#include "turbolenza/Field.hpp"
#include "turbolenza/Metrics.hpp"

#include <optional>
#include <stdexcept>

namespace BladeGen::Turbolenza {

/// v1 CoefficientError.
class CoefficientError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 CoefficientFields.
struct CoefficientFields
{
    Field A, B, C, D, E;
    double mass_flow_per_depth = 0.0;   ///< mdot' (kg/(m s), or reference units without c_ref)
    std::optional<double> c_ref;        ///< physical chord [m]; empty: normalized mode

    [[nodiscard]] double length_scale() const noexcept { return c_ref ? *c_ref : 1.0; }
    [[nodiscard]] bool physical_scale_defined() const noexcept { return c_ref.has_value(); }
};

/// v1 build_coefficients.  rho: finite, > 0, (n_m, n_eta), periodic (to ROUND_OFF max rho) on
/// the periodic stations; mdot' finite > 0; c_ref empty or finite > 0.  CoefficientError otherwise.
[[nodiscard]] CoefficientFields build_coefficients(const GridMetrics& metrics, const Field& rho,
                                                   double mass_flow_per_depth,
                                                   std::optional<double> c_ref = std::nullopt);

}  // namespace BladeGen::Turbolenza
