/// @file turbolenza/Viscous.hpp
/// @brief The boundary layer of Aungier (2003) Sec. 5.6 on a solved blade-to-blade flow, and its
///        two-way coupling to the potential flow by transpiration.
///
/// Not in v1.  The layer itself is BoundaryLayer.hpp (Aungier's method).  This file reads its
/// edge conditions from a Solution, turns its displacement into the equivalent inviscid flow and
/// iterates the two.  Aungier computes the layer one way only (Fig. 5-15: "using results from the
/// inviscid potential flow analysis results as the boundary layer edge conditions", 2003 p. 112);
/// everything about the coupling is ours.
///
/// Edge conditions (ours: Aungier names the edge state, not where it is read):
///  * the stations of a surface are its wall nodes on the blade stations first_blade_i ..
///    last_blade_i; the first one is Aungier's "first station (blade leading edge)" (2000 p. 276),
///    where theta = 0 (2003 Fig. 5-5: his first blade grid point lies inside the passage, as ours);
///  * x: the arc length of the wall node's own C2 surface curve from its LE point (u = 0), by
///    Gauss-Legendre quadrature of |X'(u)| (8 points on 8 parts of every knot span), times the chord;
///  * W_e = |W| at the wall node (the reported evaluation), rho_e and T_e the static state of that
///    speed (eos), mu = Sutherland(T_e) (BoundaryLayer.hpp), T't = H'/cp, cp of the gas;
///  * x/c the blade-frame chordwise coordinate of postprocess (Postprocess.hpp); f the blade-station
///    fraction (m_i - m_first) / (m_last - m_first), the abscissa of Aungier's Figs. 5-6 and 5-15 as we
///    read them (PORT.md "The Fig. 5-15 study"; Aungier does not state it).
///
/// The equivalent inviscid flow (Lighthill, M. J., 1958, "On displacement thickness", J. Fluid Mech.
/// 4, 383-392): the grid stays; the mass deficit m = rho_e W_e delta* of each wall node enters the
/// wall relations, psi = -m_lower/mdot' on the lower wall (j = 0) and 1 + m_upper/mdot' on the upper
/// (j = n_eta - 1), and the exit periodic stations carry the wake's: jump = 1 + m_wake/mdot'
/// (BoundaryConditions.hpp).  The wake (WakeModel):
///  * Constant: m_wake = m_lower + m_upper at the last blade station, at every exit station (the
///    jump is then continuous with the last blade station's psi_upper - psi_lower by construction);
///  * Marched: Eqs. 5-106 and 5-127 with tau_w = 0 (BoundaryLayer.hpp march_wake) along the lower
///    periodic line j = 0 from the last blade station: x the length of the j = 0 node polyline from
///    that station's lower wall node; station 0 holds the sums of both surfaces' theta, delta - delta*
///    and delta* and the mean of the two wall states (equal speeds by the Kutta condition); further
///    on, the edge state of the j = 0 node; the far exit station takes the jump the downstream
///    psi_m = 0 relation implies, 4/3 m(N-1) - 1/3 m(N-2) (its psi is eliminated by that relation).
///
/// The coupling loop (ours): BL iteration k = 0, 1, ...: the inviscid solve (density + Kutta, or
/// the fixed angle) with the transpiration m_k (k = 0: none, the inviscid solution); the layers
/// from its wall values; the mismatch e_k = max |m_BL - m_k| / max |m_BL| over both surfaces' wall
/// nodes; converged when e_k <= tolerance (one way: after k = 0, Aungier's analysis); else
/// m_{k+1} = m_k + omega (m_BL - m_k) for every wall node and exit station (the under-relaxation of
/// delta* applied to the mass deficit the wall relations carry), and the next inviscid solve starts
/// from solve k's density and grid (its exit angle).  k = max_iterations without convergence: not
/// converged (reported).  A failure of a called module in an iteration k >= 1 ends the run with
/// BoundaryLayerCouplingError (no damping beyond omega is added).
#pragma once

#include "turbolenza/BoundaryLayer.hpp"
#include "turbolenza/Case.hpp"
#include "turbolenza/Solver.hpp"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

struct ViscousSettings
{
    BlCoupling coupling = BlCoupling::TwoWay;
    WakeModel wake = WakeModel::Constant;
    BoundaryLayerSettings layer;
    SutherlandLaw viscosity;
    double relaxation = 0.2;   ///< omega_BL: 0.5 diverges on the Dunavant case (PORT.md "The Fig. 5-15 study")
    double tolerance = 1e-4;
    int max_iterations = 60;
};

/// The layer of one wall.
struct WallLayer
{
    int j = 0;                          ///< 0 (lower) or n_eta - 1 (upper)
    Surface surface = Surface::Suction;
    std::vector<int> stations;          ///< first_blade_i .. last_blade_i
    std::vector<double> x_c;
    std::vector<double> f;              ///< the blade-station fraction of the stations (blade_station_fraction)
    EdgeConditions edge;
    BoundaryLayerProfile profile;
};

/// The wake along the exit periodic line (marched model only; empty profile for the constant one).
struct WakeLayer
{
    bool marched = false;
    std::vector<int> stations;          ///< last_blade_i .. n_m - 1
    std::vector<double> m_c;            ///< (m - m_te) / chord
    EdgeConditions edge;
    BoundaryLayerProfile profile;
};

/// The mass deficits per station (n_m each): lower and upper at the blade stations, wake at the
/// exit stations, 0 elsewhere [kg/(m s)].
struct MassDeficit
{
    std::vector<double> lower, upper, wake;
};

struct ViscousAnalysis
{
    WallLayer lower, upper;
    WakeLayer wake;
    MassDeficit deficit;                ///< what this solution's layers ask the wall relations for
    ProfileLoss loss;
    [[nodiscard]] const WallLayer& suction() const noexcept { return lower.surface == Surface::Suction ? lower : upper; }
    [[nodiscard]] const WallLayer& pressure() const noexcept { return lower.surface == Surface::Suction ? upper : lower; }
};

/// The blade-station fraction f = (m_i - m_first) / (m_last - m_first) of the blade stations
/// i = first_blade_i .. last_blade_i: 0 at the first, 1 at the last.  BoundaryLayerError unless the grid
/// has two or more blade stations with increasing m.
[[nodiscard]] std::vector<double> blade_station_fraction(const Grid& grid);

/// The arc length (chord units) of each wall node of wall j along its surface curve from the LE.
[[nodiscard]] std::vector<double> wall_arc_length(const Grid& grid, int j);

/// The layers of both walls (and the marched wake) of a solution.  gas and the relative total
/// temperature T't are the solve's.  BoundaryLayerError when the layer cannot be computed.
[[nodiscard]] ViscousAnalysis analyze_boundary_layers(const Solution& solution, const PerfectGas& gas,
                                                      double total_temperature, const ViscousSettings& settings);

/// The wall psi and periodic jumps of a mass deficit (BoundaryConditions.hpp).
[[nodiscard]] Transpiration transpiration_of(const Grid& grid, const MassDeficit& deficit, double mass_flow_per_depth);

/// e = max |a - b| / max |a| over the blade stations of both walls (0 when a is 0 everywhere).
[[nodiscard]] double deficit_change(const Grid& grid, const MassDeficit& a, const MassDeficit& b);

/// The outcome of the coupled iteration (no exception escapes: errors are in `error`).
struct ViscousSolve
{
    std::shared_ptr<const Solution> solution;   ///< the reported solve (the last one), or the last completed pass
    std::shared_ptr<const Solution> inviscid;   ///< k = 0
    std::optional<ViscousAnalysis> analysis;    ///< the last layers computed (of `analysis_solution`)
    /// The solution `analysis` was computed from: `solution` itself, or, when a later solve failed, the
    /// previous iteration's solution (the last one whose layers exist).
    std::shared_ptr<const Solution> analysis_solution;
    int analysis_iteration = -1;                ///< the BL iteration k of `analysis` (-1: none)
    std::vector<BlIteration> history;
    MassDeficit used;                           ///< the deficit `solution` was solved with
    bool converged = false;
    std::string error;                          ///< "<ErrorType>: message"; empty when converged
    std::optional<int> failed_pass;
};

/// The coupled iteration from a built grid (coupled: the density + Kutta loop, else the fixed angle).
[[nodiscard]] ViscousSolve solve_viscous(std::shared_ptr<const Grid> grid, const SolverConfig& config, double c_ref,
                                         std::optional<double> mass_flow_per_depth, const SolveOptions& options,
                                         const ViscousSettings& settings, bool coupled);

}  // namespace BladeGen::Turbolenza
