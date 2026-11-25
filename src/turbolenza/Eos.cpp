// Port of v1 eos.py (see Eos.hpp).  Operation order as v1's NumPy expressions.
#include "turbolenza/Eos.hpp"

#include "turbolenza/Detail.hpp"

#include <cmath>
#include <format>
#include <string>
#include <utility>

namespace BladeGen::Turbolenza {

using detail::py_repr;

namespace {

double real(double value, const char* name)
{
    if (!std::isfinite(value)) throw ThermodynamicStateError(std::format("{} is {}: it must be finite", name, py_repr(value)));
    return value;
}

void finite_field(const Field& f, const char* name)
{
    int bad = 0;
    for (double v : f.v)
        if (!std::isfinite(v)) ++bad;
    if (bad) throw ThermodynamicStateError(std::format("{} is not finite at {} node(s)", name, bad));
}

}  // namespace

PerfectGas::PerfectGas(double gamma_, double gas_constant_)
{
    const double g = real(gamma_, "gamma"), R = real(gas_constant_, "the gas constant R");
    if (!(g > 1.0)) throw ThermodynamicStateError(std::format("gamma = {}: a perfect gas needs gamma > 1", py_repr(g)));
    if (!(R > 0.0)) throw ThermodynamicStateError(std::format("R = {}: the gas constant must be positive", py_repr(R)));
    gamma = g;
    gas_constant = R;
}

TotalState::TotalState(double pressure_, double temperature_)
{
    const double p = real(pressure_, "the total pressure"), T = real(temperature_, "the total temperature");
    if (!(p > 0.0))
        throw ThermodynamicStateError(std::format("the total pressure is {}: it must be positive", py_repr(p)));
    if (!(T > 0.0))
        throw ThermodynamicStateError(std::format("the total temperature is {}: it must be positive", py_repr(T)));
    pressure = p;
    temperature = T;
}

// v1 blade_row_reference: I = H_in - omega r C_theta_in, H_in = cp Tt_in.
BladeRowReference blade_row_reference(const TotalState& total, const PerfectGas& gas, double omega,
                                      std::optional<double> radius, std::optional<double> inlet_c_theta)
{
    omega = real(omega, "omega");
    const double H_in = gas.cp() * total.temperature;
    if (omega == 0.0) {
        if (radius) {
            const double r = real(*radius, "the radius");
            if (!(r > 0.0))
                throw ThermodynamicStateError(
                    std::format("the stream-surface radius is {}: it must be positive", py_repr(r)));
        }
        return BladeRowReference{total, H_in, 0.0, radius, {}};
    }
    if (!radius || !inlet_c_theta)
        throw ThermodynamicStateError("a rotating row (omega != 0) needs radius and inlet_c_theta (the inlet absolute "
                                      "tangential velocity) for its rothalpy");
    const double r = real(*radius, "the radius");
    const double c_theta = real(*inlet_c_theta, "the inlet absolute tangential velocity");
    if (!(r > 0.0))
        throw ThermodynamicStateError(std::format("the stream-surface radius is {}: it must be positive", py_repr(r)));
    return BladeRowReference{total, H_in - omega * r * c_theta, omega, r, {}};
}

ScalarState static_state_at(double W, const PerfectGas& gas, const BladeRowReference& reference)
{
    if (!std::isfinite(W)) throw ThermodynamicStateError("the speed W is not finite at 1 node(s)");
    if (W < 0.0) throw ThermodynamicStateError("the speed W is negative at 1 node(s): W is a magnitude");
    const double H_rel = reference.relative_total_enthalpy();
    const double cp = gas.cp();
    const double h = H_rel - 0.5 * (W * W);
    const double T = h / cp;
    if (!(T > 0.0))
        throw ThermodynamicStateError(std::format(
            "T <= 0 at 1 node(s) (largest W {:.6g}): the kinetic energy W^2 / 2 reaches the relative total enthalpy "
            "H' = {:.6g}, i.e. W >= sqrt(2 H') = {:.6g}. Nothing is clipped",
            W, H_rel, std::sqrt(std::max(2.0 * H_rel, 0.0))));
    const TotalState& total = reference.inlet_total_state;
    const double exponent = gas.gamma / (gas.gamma - 1.0);
    ScalarState s;
    s.enthalpy = h;
    s.temperature = T;
    s.pressure = total.pressure * std::pow(T / total.temperature, exponent);
    s.density = s.pressure / (gas.gas_constant * T);
    s.sound_speed = std::sqrt(gas.gamma * gas.gas_constant * T);
    return s;
}

// v1 static_state_from_velocity: h = H' - W^2 / 2, T = h / cp, p = Pt (T / Tt)^(g / (g - 1)),
// rho = p / (R T), a = sqrt(gamma R T).
StaticState static_state_from_velocity(const Field& W, const PerfectGas& gas, const BladeRowReference& reference)
{
    finite_field(W, "the speed W");
    int negative = 0;
    for (double w : W.v)
        if (w < 0.0) ++negative;
    if (negative)
        throw ThermodynamicStateError(
            std::format("the speed W is negative at {} node(s): W is a magnitude", negative));
    const double H_rel = reference.relative_total_enthalpy();
    const double cp = gas.cp();
    StaticState s{Field(W.n_m, W.n_eta), Field(W.n_m, W.n_eta), Field(W.n_m, W.n_eta), Field(W.n_m, W.n_eta),
                  Field(W.n_m, W.n_eta)};
    const int n = static_cast<int>(W.size());
    const bool stations = !reference.station_blade_speed.empty();
    if (stations && static_cast<int>(reference.station_blade_speed.size()) != W.n_m)
        throw ThermodynamicStateError(std::format("the reference has {} station blade speeds, the field {} stations",
                                                  reference.station_blade_speed.size(), W.n_m));
    const int n_eta = W.n_eta;
    TURBOLENZA_PARALLEL_FOR
    for (int k = 0; k < n; ++k) {
        const std::size_t q = static_cast<std::size_t>(k);
        const double w = W.v[q];
        // Eq. 5-51: h = I + (r omega)^2 / 2 - W^2 / 2, at the node's station on a rotating stream sheet
        const double H = stations ? reference.relative_total_enthalpy_at(k / n_eta) : H_rel;
        const double h = H - 0.5 * (w * w);
        s.enthalpy.v[q] = h;
        s.temperature.v[q] = h / cp;
    }
    int bad = 0;
    double W_max = 0.0;
    for (std::size_t q = 0; q < W.size(); ++q) {
        if (!(s.temperature.v[q] > 0.0)) ++bad;
        W_max = std::max(W_max, W.v[q]);
    }
    if (bad)
        throw ThermodynamicStateError(std::format(
            "T <= 0 at {} node(s) (largest W {:.6g}): the kinetic energy W^2 / 2 reaches the relative total enthalpy "
            "H' = {:.6g}, i.e. W >= sqrt(2 H') = {:.6g}. Nothing is clipped",
            bad, W_max, H_rel, std::sqrt(std::max(2.0 * H_rel, 0.0))));
    const TotalState& total = reference.inlet_total_state;
    const double exponent = gas.gamma / (gas.gamma - 1.0);
    const double gR = gas.gamma * gas.gas_constant;
    const double R = gas.gas_constant;
    TURBOLENZA_PARALLEL_FOR
    for (int k = 0; k < n; ++k) {
        const std::size_t q = static_cast<std::size_t>(k);
        const double T = s.temperature.v[q];
        const double p = total.pressure * std::pow(T / total.temperature, exponent);
        s.pressure.v[q] = p;
        s.density.v[q] = p / (R * T);
        s.sound_speed.v[q] = std::sqrt(gR * T);
    }
    return s;
}

double initial_density(const PerfectGas& gas, const BladeRowReference& reference)
{
    return static_state_at(0.0, gas, reference).density;
}

Field initial_density(const PerfectGas& gas, const BladeRowReference& reference, int n_m, int n_eta)
{
    return Field(n_m, n_eta, initial_density(gas, reference));
}

Field relax_density(const Field& density_old, const Field& density_target, double relaxation)
{
    const double lam = real(relaxation, "the density relaxation");
    if (!(0.0 < lam && lam <= 1.0))
        throw ThermodynamicStateError(
            std::format("the density relaxation is {}: it must satisfy 0 < lambda <= 1", py_repr(lam)));
    finite_field(density_old, "the old density");
    finite_field(density_target, "the target density");
    if (!density_old.same_shape(density_target))
        throw ThermodynamicStateError(std::format("the old density has shape ({}, {}), the target ({}, {})",
                                                  density_old.n_m, density_old.n_eta, density_target.n_m,
                                                  density_target.n_eta));
    for (const auto& [name, f] : {std::pair<const char*, const Field*>{"old", &density_old},
                                  std::pair<const char*, const Field*>{"target", &density_target}})
        for (double v : f->v)
            if (!(v > 0.0)) throw ThermodynamicStateError(std::format("the {} density is not positive everywhere", name));
    Field out(density_old.n_m, density_old.n_eta);
    const double keep = 1.0 - lam;
    const int n = static_cast<int>(out.size());
    TURBOLENZA_PARALLEL_FOR
    for (int k = 0; k < n; ++k) {
        const std::size_t q = static_cast<std::size_t>(k);
        out.v[q] = keep * density_old.v[q] + lam * density_target.v[q];
    }
    return out;
}

}  // namespace BladeGen::Turbolenza
