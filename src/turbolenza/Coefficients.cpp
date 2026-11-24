// Port of v1 coefficients.py (see Coefficients.hpp): build_coefficients, _mass_flow,
// _reference_chord, _density, _metric_inputs, _finite, _d_dm, _d_deta.
#include "turbolenza/Coefficients.hpp"

#include "turbolenza/Detail.hpp"

#include <cmath>
#include <format>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

using detail::py_repr;
using detail::sz;

namespace {

double mass_flow(double value)
{
    if (!(std::isfinite(value) && value > 0.0))
        throw CoefficientError(std::format("mass_flow_per_depth = {}: the mass flow through one passage per unit depth "
                                           "must be finite and positive",
                                           py_repr(value)));
    return value;
}

std::optional<double> reference_chord(std::optional<double> value)
{
    if (!value) return std::nullopt;
    if (!(std::isfinite(*value) && *value > 0.0))
        throw CoefficientError(std::format("c_ref = {}: the physical chord must be a finite positive length in metres, "
                                           "or None for normalized mode",
                                           py_repr(*value)));
    return value;
}

void check_density(const Field& rho, int n_m, int n_eta, const std::vector<char>& blade)
{
    if (rho.n_m != n_m || rho.n_eta != n_eta)
        throw CoefficientError(std::format("rho has shape ({}, {}), the grid needs the node field ({}, {}) (n_m, n_eta)",
                                           rho.n_m, rho.n_eta, n_m, n_eta));
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (!std::isfinite(rho(i, j)))
                throw CoefficientError(std::format("rho is not finite at node ({}, {}): {}", i, j, py_repr(rho(i, j))));
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (!(rho(i, j) > 0.0))
                throw CoefficientError(std::format("rho = {} at node ({}, {}): the density must be strictly positive",
                                                   py_repr(rho(i, j)), i, j));
    double rho_max = rho.v[0];
    for (double r : rho.v) rho_max = std::max(rho_max, r);
    const double tolerance = kRoundOff * rho_max;
    double worst = -1.0;
    int at = -1;
    bool bad = false;
    for (int i = 0; i < n_m; ++i) {
        if (blade[sz(i)]) continue;
        const double jump = std::abs(rho(i, 0) - rho(i, n_eta - 1));
        if (!(jump <= tolerance)) bad = true;
        if (jump > worst) {
            worst = jump;
            at = i;
        }
    }
    if (bad)
        throw CoefficientError(std::format("rho at the periodic station {} is {:.17g} at eta = 0 and {:.17g} at eta = 1: "
                                           "the two nodes are one point and its periodic image, so rho must agree there "
                                           "(to {:.3g})",
                                           at, rho(at, 0), rho(at, n_eta - 1), tolerance));
}

void finite(const char* name, const Field& f, const GridMetrics& metrics, const Field& rho)
{
    for (int i = 0; i < f.n_m; ++i)
        for (int j = 0; j < f.n_eta; ++j)
            if (!std::isfinite(f(i, j)))
                throw CoefficientError(std::format("{} is not finite at node ({}, {}) (S* = {:.6g}, cos(beta) = {:.3g}, "
                                                   "rho = {:.6g}): the coefficient is undefined there, and nothing is "
                                                   "clipped to hide it",
                                                   name, i, j, metrics.S[sz(i)], metrics.cos_beta(i, j), rho(i, j)));
}

// v1 _d_dm: one stencil over the whole m range.
double d_dm_at(const Field& f, int i, int j, double dm)
{
    const int N = f.n_m - 1;
    if (i == 0) return (-3.0 * f(0, j) + 4.0 * f(1, j) - f(2, j)) / (2.0 * dm);
    if (i == N) return (3.0 * f(N, j) - 4.0 * f(N - 1, j) + f(N - 2, j)) / (2.0 * dm);
    return (f(i + 1, j) - f(i - 1, j)) / (2.0 * dm);
}

// v1 _d_deta: central inside; blade stations one-sided at the walls; periodic stations the
// periodic central difference at both images.
double d_deta_at(const Field& f, int i, int j, double deta, bool blade)
{
    const int J = f.n_eta - 1;
    if (j > 0 && j < J) return (f(i, j + 1) - f(i, j - 1)) / (2.0 * deta);
    if (!blade) return (f(i, 1) - f(i, J - 1)) / (2.0 * deta);
    if (j == 0) return (-3.0 * f(i, 0) + 4.0 * f(i, 1) - f(i, 2)) / (2.0 * deta);
    return (3.0 * f(i, J) - 4.0 * f(i, J - 1) + f(i, J - 2)) / (2.0 * deta);
}

}  // namespace

CoefficientFields build_coefficients(const GridMetrics& metrics, const Field& rho, double mass_flow_per_depth,
                                     std::optional<double> c_ref, double omega)
{
    if (!metrics.grid) throw CoefficientError("build_coefficients takes a GridMetrics with its grid");
    const double mdot = mass_flow(mass_flow_per_depth);
    c_ref = reference_chord(c_ref);
    const double length = c_ref ? *c_ref : 1.0;
    const Grid& grid = *metrics.grid;
    if (!(std::abs(grid.chord() - 1.0) <= kRoundOff))
        throw CoefficientError(std::format("the grid's chord is {:.17g}, not 1: coefficients.py takes a chord-normalized "
                                           "grid and applies the physical chord once, as c_ref. Build the cascade with "
                                           "physical_chord=None and pass c_ref={:.17g} instead",
                                           grid.chord(), grid.chord()));
    const int n_m = grid.n_m(), n_eta = grid.n_eta();
    if (n_eta < 3)
        throw CoefficientError(std::format("the grid has {} eta nodes; a second-order eta difference needs three", n_eta));
    std::vector<char> blade(sz(n_m), 0);
    for (int i = grid.first_blade_i; i <= grid.last_blade_i; ++i) blade[sz(i)] = 1;
    check_density(rho, n_m, n_eta, blade);

    // v1 _metric_inputs
    for (int i = 0; i < n_m; ++i) {
        const double s = metrics.S[sz(i)];
        if (!(std::isfinite(s) && s > 0.0))
            throw CoefficientError(
                std::format("S* = {} at station {}: the passage spacing must be finite and positive", py_repr(s), i));
    }
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (!std::isfinite(metrics.tan_beta(i, j)))
                throw CoefficientError(
                    std::format("tan(beta) = {} at node ({}, {}): it must be finite", py_repr(metrics.tan_beta(i, j)), i, j));
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (!(std::isfinite(metrics.cos_beta(i, j)) && metrics.cos_beta(i, j) > 0.0))
                throw CoefficientError(std::format("cos(beta) = {} at node ({}, {}): it must be finite and positive",
                                                   py_repr(metrics.cos_beta(i, j)), i, j));
    std::vector<double> S(sz(n_m));
    for (int i = 0; i < n_m; ++i) S[sz(i)] = length * metrics.S[sz(i)];

    CoefficientFields out{Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta), Field(n_m, n_eta),
                          Field(n_m, n_eta), mdot, c_ref, Field(), 0.0};
    const int n = n_m * n_eta;
    const bool sheet = metrics.on_sheet;
    if (!std::isfinite(omega))
        throw CoefficientError(std::format("omega = {}: the rotation speed must be finite", py_repr(omega)));
    if (omega != 0.0 && !sheet)
        throw CoefficientError("a rotating row (omega != 0) needs a grid on a stream sheet: its radius r(m) and "
                               "sin(phi) = dr/dm carry the rotation (Eqs. 5-21, 5-51)");
    TURBOLENZA_PARALLEL_FOR
    for (int k = 0; k < n; ++k) {
        const std::size_t q = static_cast<std::size_t>(k);
        const double s = S[sz(k / n_eta)];
        // rho b on a stream sheet (b / b_in, Eqs. 5-22 .. 5-24); v1's rho otherwise
        const double r = sheet ? rho.v[q] * metrics.thickness_at(k / n_eta) : rho.v[q];
        const double c = metrics.cos_beta.v[q], t = metrics.tan_beta.v[q];
        out.A.v[q] = mdot / (s * r * (c * c));  // Eq. 5-22
        out.B.v[q] = mdot * t / r;              // Eq. 5-23
        out.C.v[q] = mdot * s / r;              // Eq. 5-24
    }
    out.omega = omega;
    if (omega != 0.0) {  // not in v1: the rotation source of Eq. 5-21, 2 S omega sin(phi) [Eqs. 5-13 .. 5-15]
        out.Q = Field(n_m, n_eta);
        for (int i = 0; i < n_m; ++i)
            for (int j = 0; j < n_eta; ++j) out.Q(i, j) = 2.0 * S[sz(i)] * omega * metrics.sin_phi_at(i);
    }
    finite("A", out.A, metrics, rho);
    finite("B", out.B, metrics, rho);
    finite("C", out.C, metrics, rho);
    const double dm = length * grid.dm(), deta = grid.deta();
    TURBOLENZA_PARALLEL_FOR
    for (int k = 0; k < n; ++k) {
        const int i = k / n_eta, j = k % n_eta;
        const bool b = blade[sz(i)] != 0;
        const std::size_t q = static_cast<std::size_t>(k);
        out.D.v[q] = d_deta_at(out.A, i, j, deta, b) - d_dm_at(out.B, i, j, dm);  // Eq. 5-25
        out.E.v[q] = d_dm_at(out.C, i, j, dm) - d_deta_at(out.B, i, j, deta, b);  // Eq. 5-26
    }
    finite("D", out.D, metrics, rho);
    finite("E", out.E, metrics, rho);
    return out;
}

}  // namespace BladeGen::Turbolenza
