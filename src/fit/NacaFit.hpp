/// @file fit/NacaFit.hpp
/// @brief "NACA 65 + Custom": BladeGen profile parameters fitted to a NACA 65 section.
///
/// The parametric NACA 65 section (naca/Naca65.hpp) is the target of the inverse fit
/// (fit/InverseFit.hpp); the result is an ordinary ProfileParams, so every custom knob
/// can then edit the NACA geometry.  The start is BladeGen's 65-(12)10 fit scaled to the
/// requested Cl0 and thickness (edge angles with the equivalent-arc angles, thickness
/// values with t, the LE radius with t^2, the TE radius with the variant).
#pragma once

#include "fit/InverseFit.hpp"
#include "naca/Naca65.hpp"

namespace BladeGen::Fit {

/// Start parameters for the fit of a NACA 65 section.
[[nodiscard]] Profile::ProfileParams naca65_start(const Naca::Naca65Params& p);

/// Fit options for a NACA 65 target: fairness weight 1, two LM stages of 40 iterations,
/// 50 / 140 samples per connector / side (RMS ~1e-4 c in 3-5 s; more iterations do not
/// improve it).
[[nodiscard]] FitOptions naca65_fit_options();

/// Fit custom profile parameters to the section (target: n_target points per side,
/// clustered at both edges).  opt.params empty: default_fit_params(start).
[[nodiscard]] FitResult fit_naca65(const Naca::Naca65Section& s, const FitOptions& opt = naca65_fit_options(),
                                   int n_target = 100);

}  // namespace BladeGen::Fit
