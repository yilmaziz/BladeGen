/// @file turbolenza/Eos.hpp
/// @brief The thermodynamic layer of the density loop: relative speed W -> static state, and the
///        relaxed density of the next pass.
///
/// Port of v1 eos.py (TURBOlenza v1 commit 67599a6; Aungier 2003 Ch. 5, the energy relation
/// h = H' - W^2 / 2 with H' = I + (r omega)^2 / 2, I = cp Tt_in - omega r C_theta_in).
///
///  * the gas: thermally and calorically perfect, p = rho R T, cp = gamma R / (gamma - 1),
///    h = cp T, a = sqrt(gamma R T) (v1 CaloricallyPerfectGas);
///  * isentropic to the inlet total state: p = Pt_in (T / Tt_in)^(gamma / (gamma - 1));
///  * static_state_from_velocity: h = H' - W^2 / 2, T = h / cp, p, rho_EOS = p / (R T), a; T <= 0
///    is a ThermodynamicStateError (nothing clipped);
///  * initial_density: Aungier's first guess, the state at W = 0 (Pt / (R Tt) when stationary);
///  * relax_density: rho_next = (1 - lambda) rho + lambda rho_target, 0 < lambda <= 1
///    (lambda = 1 Aungier's direct replacement, 0.5 the v1 project baseline).
/// The per-node maps run on the threads of the current solve (Detail.hpp).
#pragma once

#include "turbolenza/Field.hpp"

#include <optional>
#include <stdexcept>

namespace BladeGen::Turbolenza {

inline constexpr double kDensityRelaxation = 0.5;          ///< v1 DENSITY_RELAXATION (project baseline)
inline constexpr double kAungierDirectReplacement = 1.0;   ///< v1 AUNGIER_DIRECT_REPLACEMENT

/// v1 ThermodynamicStateError.
class ThermodynamicStateError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 CaloricallyPerfectGas (alias PerfectGas).  The default is air as in CaseSettings.
struct PerfectGas
{
    double gamma = 1.4;
    double gas_constant = 287.05;

    PerfectGas() = default;
    /// Throws ThermodynamicStateError unless gamma > 1 and R > 0, both finite (v1 __post_init__).
    PerfectGas(double gamma_, double gas_constant_);

    [[nodiscard]] double cp() const noexcept { return gamma * gas_constant / (gamma - 1.0); }
    [[nodiscard]] double cv() const noexcept { return gas_constant / (gamma - 1.0); }
};

/// v1 TotalState: the inlet absolute total state, the entropy reference.
struct TotalState
{
    double pressure = 0.0, temperature = 0.0;
    TotalState() = default;
    /// Throws ThermodynamicStateError unless both are finite and positive.
    TotalState(double pressure_, double temperature_);
};

/// v1 BladeRowThermoReference (built by blade_row_reference).
struct BladeRowReference
{
    TotalState inlet_total_state;
    double rothalpy = 0.0;
    double omega = 0.0;
    std::optional<double> radius;

    [[nodiscard]] double blade_speed() const noexcept { return omega == 0.0 ? 0.0 : *radius * omega; }
    [[nodiscard]] double relative_total_enthalpy() const noexcept
    {
        const double U = blade_speed();
        return rothalpy + 0.5 * (U * U);
    }
};

/// v1 blade_row_reference: I = cp Tt_in - omega r C_theta_in; a stationary row needs neither
/// radius nor swirl, a rotating one both (ThermodynamicStateError otherwise).
[[nodiscard]] BladeRowReference blade_row_reference(const TotalState& total, const PerfectGas& gas, double omega = 0.0,
                                                    std::optional<double> radius = std::nullopt,
                                                    std::optional<double> inlet_c_theta = std::nullopt);

/// v1 StaticState of a field of speeds.
struct StaticState
{
    Field temperature, pressure, density, enthalpy, sound_speed;
};

/// The static state at one speed (v1 static_state_from_velocity on a number).
struct ScalarState
{
    double temperature = 0.0, pressure = 0.0, density = 0.0, enthalpy = 0.0, sound_speed = 0.0;
};

/// v1 static_state_from_velocity on a field: W finite and >= 0, T > 0 everywhere, else
/// ThermodynamicStateError with v1's message.
[[nodiscard]] StaticState static_state_from_velocity(const Field& W, const PerfectGas& gas,
                                                     const BladeRowReference& reference);
[[nodiscard]] ScalarState static_state_at(double W, const PerfectGas& gas, const BladeRowReference& reference);

/// v1 initial_density: the relative total density (the state at W = 0), uniform.
[[nodiscard]] double initial_density(const PerfectGas& gas, const BladeRowReference& reference);
[[nodiscard]] Field initial_density(const PerfectGas& gas, const BladeRowReference& reference, int n_m, int n_eta);

/// v1 relax_density: (1 - lambda) old + lambda target; 0 < lambda <= 1, both finite, positive,
/// of one shape (ThermodynamicStateError otherwise).
[[nodiscard]] Field relax_density(const Field& density_old, const Field& density_target,
                                  double relaxation = kDensityRelaxation);

}  // namespace BladeGen::Turbolenza
