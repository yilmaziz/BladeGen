// tests/turbolenza/CouplingFixture.hpp — through-flow cases with the Chapter 6 and Chapter 8 adapters
// (ThroughflowCoupling.hpp) attached, shared by the coupling and the sheet tests.  Defined in
// test_throughflow_coupling.cpp.
#pragma once

#include "turbolenza/Throughflow.hpp"

namespace coupling {

namespace TF = BladeGen::Turbolenza::Throughflow;

struct StageOptions
{
    TF::ThroughflowModel model = TF::ThroughflowModel::Approximate;
    int n_stream_surfaces = 11;
    bool endwall = true;    ///< Ch. 8 (ChapterEightEndwall); false: no blockage (K_B = 1)
    bool row_loss = true;   ///< Ch. 6 Secs. 6.9, 6.10 (rotor tip gap, shrouded stator seal)
    bool mirror = false;    ///< every angle, swirl and omega sign-changed, s_row too (theta -> -theta)
    double rpm = 5000.0;
    bool reynolds_correction = false;  ///< Ch. 6's Eqs. 6-100 .. 6-102 (K2 with K_Re of Re_c = rho1 W1 c / mu(T1))
};

/// The synthetic stage: inlet guide vane, rotor and shrouded stator on a straight annulus r = 0.25 .. 0.35 m,
/// 7 radial quasi-normals, mdot = 22 kg/s, every row Chapter 6 (CorrelationBladeRow, CorrelationRowLoss).
[[nodiscard]] TF::ThroughflowCase synthetic_stage(const StageOptions& o = {});
inline constexpr int kIgv = 0, kRotor = 1, kStator = 2;

/// BladeGen's default design (one rotor) at 15000 rpm, 3 kg/s, radial quasi-normals, with Chapter 6 (its custom
/// profile taken as NACA 65 thickness on the NACA 65 mean line, ours) and Chapter 8, a 0.3 mm tip gap.
[[nodiscard]] TF::ThroughflowCase bladegen_rotor_case(TF::ThroughflowModel model, int* qn_before_le = nullptr);

/// A rotor (SpecifiedBladeRow: beta'_2 = 0.6 rad, omegabar = 0.03) on a 15-degree conical annulus (both walls
/// sloped, radial quasi-normals at z = 0, 0.05, 0.10, 0.16, 0.22 m, the row between the third and fourth), inlet
/// C_theta = 20 m/s, 22 kg/s, 5000 rpm, no end-wall model, the approximate model.
[[nodiscard]] TF::ThroughflowCase cone_rotor_case();

}  // namespace coupling
