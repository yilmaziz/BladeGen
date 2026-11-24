/// @file turbolenza/ThroughflowCoupling.hpp
/// @brief The integration adapters of Aungier's compressor analysis (Axial-Flow Compressors, ASME Press,
///        2003): Chapter 6's cascade correlations as Chapter 7's BladeRowModel and RowLossModel, and
///        Chapter 8's end-wall boundary layers as Chapter 7's EndwallModel.  The third adapter, a Chapter 7
///        stream surface as the blade-to-blade solver's stream sheet (Ch. 12 p. 287), is ThroughflowSheet.hpp.
///
/// Nothing here is a new correlation or a new equation: every number comes from Correlations.hpp (Ch. 6),
/// Endwall.hpp (Ch. 8), Throughflow.hpp (Ch. 7) or the blade-to-blade case (Case.hpp, Ch. 5).  The adapters
/// map one module's state into the other's frame and back.  Every mapping the book does not state is
/// marked "ours" here, once.
///
/// FRAMES (Throughflow.hpp CONVENTIONS, binding):
///  * Chapter 7 and Chapter 8 use Aungier's frame with a signed omega: W_theta = C_theta - omega r, angles
///    from the meridional direction signed toward +theta.  TURBOlenza's and BladeGen's theta is this theta and
///    omega = -2 pi rpm / 60.  Both modules take each other's quantities verbatim (no sign changes).
///  * Chapter 6's correlations are written in the positive-camber CASCADE frame (Fig. 6-1).  The row sign
///    s_row (BladeRowSpec::s_row) maps the book frame onto it: cascade angle = s_row * book angle, for the
///    stagger, the metal angles, the flow angles, W_theta and omega (Correlations.hpp does the multiplication;
///    CascadeGeometry::row_sign = s_row), and the returned beta'_2 = s_row * beta_2,cascade (p. 128's rule,
///    generalised to every row: the q3d contract Sec. 0.2).  A guide vane of negative camber after s_row gets
///    p. 128's own flip on top (Correlations.hpp).
///  * The blade-to-blade solver's blade has its suction side toward +y (Cascade.hpp).  A row with s_row = -1
///    is handed to it mirrored (theta -> -theta: rpm and the inlet swirl sign-changed, the caller's blade the
///    mirrored section); its angles times s_row are the book's (ours, the q3d contract Sec. 0.3).
///
/// CHAPTER 6 -> CHAPTER 7 (CorrelationBladeRow, CorrelationRowLoss), per stream surface j of a row between
/// the quasi-normals QN1 (inlet) and QN2 (discharge), both outside the blade (p. 155):
///  * inlet (Correlations InletState): r1, W_m1, W_theta1 = C_theta1 - omega r1 (Eq. 7-44), omega, M'_1 =
///    W1 / a1, rho1 and mu1 = mu(T1) by Sutherland's law at the static T1 (the through-flow's gas): Re_c =
///    rho1 W1 c / mu1 (ours, the book names Re_c only; Correlations.hpp computes it from these);
///  * discharge estimate (DischargeEstimate): r2 and the through-flow's current W_m2 (p. 163);
///  * the gas of Eqs. 6-53, 6-59, 6-60, 6-70 .. 6-75, 6-91 from the inlet state itself: R = P1 / (rho1 T1),
///    gamma = a1^2 rho1 / P1 (exact for the through-flow's perfect gas; ours, so one gas serves both);
///  * the section on the surface (BladeRowInlet::geometry, the through-flow's linear interpolation in span of
///    BladeRowSpec::sections, which a BladeGen layer set fills through ThroughflowBladeGen.hpp or a plain table
///    by hand): chord, camber theta = kappa_1 - kappa_2, t_b / c, a / c (0.5 for the circular arc and the NACA 65
///    equivalent arc, Eq. 4-35, whatever the section's measured value; ours); solidity c / s with the through-flow's
///    pitch s = 2 pi r_mean / Z, r_mean = (r1 + r2) / 2; o / s = throat / s (ours: o is taken where the
///    section gives it, s at r_mean);
///  * the stagger (StaggerSource): MetalAngles (default, ours) places the chord line so that the Ch. 4
///    camberline of camber theta has the section's metal angles exactly, gamma = kappa_1 - chi_1(theta, a / c)
///    in the cascade frame (for a circular arc gamma = (kappa_1 + kappa_2) / 2): Ch. 6 measures incidence and
///    deviation from kappa_1, kappa_2 (6-2, 6-3), so the true metal angles are kept; Geometric passes the
///    section's own stagger, and Ch. 4's kappa then differ from the section's by the mismatch (reported);
///  * the row data of Eqs. 6-99, 6-100: s_bar / h_bar with s_bar = 2 pi r_bar / Z, r_bar the mean of the four
///    LE and TE end radii of the row (ours) and h_bar the through-flow's blade height;
///  * out: beta'_2 (book frame), omega_bar = CascadePerformance::loss (profile, Mach, shock and off-design,
///    already the Eq. 7-49 form, Delta P't / (P't1 - P1) at QN1), K_W = 1 - B_wake (Eq. 7-6).
///  * the row-level losses (Secs. 6.9, 6.10, endwall_leakage): the gap-side wall states at QN1 and QN2 from
///    the through-flow (absolute C_theta, book sign), the gap radius the mean of the two wall radii (ours), the
///    wall section's chord and stagger; a tip gap where the row has a clearance (shroud or hub, one only), a
///    shroud seal on CorrelationRowLossOptions::seal_wall for a shrouded row (ours: the hub, a stator's inner
///    band) with the seal geometry of the options; mdot_leak signed + from discharge to inlet, for Ch. 8.
///
/// CHAPTER 8 -> CHAPTER 7 (ChapterEightEndwall): Fig. 7-3 calls it once per quasi-normal per flow iteration,
/// in order, and repeats a quasi-normal until its flow converges.  The adapter keeps the march of the pass:
/// at QN 0 endwall_inlet (Sec. 8.8, B of EndwallInlet split equally, Eq. 8-72), at QN i > 0 endwall_step from
/// the layer the last call for QN i - 1 left (the converged one), with the row data when a blade row lies
/// between i - 1 and i (the step spans it, p. 303), then endwall_blockage: K_B = 1 - B (Eqs. 8-76, 8-77) of
/// QN i, the only value returned (p. 195; no force defect and no loss enter Ch. 7, p. 175).  A repeated call
/// for QN i recomputes from the same layer i - 1: the same input gives the same K_B.  The edge is the end-wall
/// stream surface at the quasi-normal foot (p. 164), velocities absolute (C_m, C_theta), the wall frame
/// omega_w of the segment ending there (QuasiNormal::hub_wall_omega / shroud_wall_omega), omega_c =
/// ThroughflowCase::compressor_omega (Eq. 8-61).  Row walls: s_bar and gamma_bar of the through-flow, the
/// clearance at that wall; a shrouded row carries its shroud band on seal_wall (ours: the hub), where the
/// clearance is 0 (its leakage is the seal's, Sec. 8.6) and the seal leakage is the row-level model's.  No
/// damping of K_B (p. 195 asks for "some limits and numerical damping" without a form; none is applied, ours).
///
/// CHAPTER 7 -> BLADE-TO-BLADE: ThroughflowSheet.hpp (a header of its own: Correlations.hpp and Case.hpp,
/// through Flow.hpp, both declare BladeGen::Turbolenza::InletState, so no translation unit can include both).
#pragma once

#include "turbolenza/Correlations.hpp"
#include "turbolenza/Endwall.hpp"
#include "turbolenza/Throughflow.hpp"

#include <memory>
#include <optional>
#include <vector>

namespace BladeGen::Turbolenza::Throughflow {

// ─────────────────────────────────────────────────────────────────────────────
// Chapter 6 -> Chapter 7: the blade row on one stream surface
// ─────────────────────────────────────────────────────────────────────────────

/// Where the cascade stagger handed to Ch. 6 comes from (see the file comment).
enum class StaggerSource
{
    MetalAngles,  ///< gamma = kappa_1 - chi_1(theta, a/c) of the Ch. 4 camberline: the section's kappa kept (ours)
    Geometric,    ///< the section's chord-line stagger as it is
};

struct CorrelationRowOptions
{
    /// The profile family and camberline of the row; empty: from BladeRowSpec::profile_family ("NACA65":
    /// NACA thickness on the NACA 65 mean line; "C4", "DCA": circular arc; "A4K6": NACA thickness on the
    /// A4K6 line).  Any other family (BladeGen's "custom") must be given here: refused otherwise.
    std::optional<ProfileFamily> profile;
    std::optional<Camberline> camberline;
    StaggerSource stagger = StaggerSource::MetalAngles;
    CorrelationOptions correlations;  ///< Aungier's defaults (Correlations.hpp)
};

/// Chapter 6 as Chapter 7's BladeRowModel (see the file comment).  Built for one row; evaluate is a pure
/// function of its input.  Throws ThroughflowError for a row it cannot serve, CorrelationError (from
/// Correlations.hpp) for an input the correlations refuse.
class CorrelationBladeRow final : public BladeRowModel
{
public:
    /// The row's Z, LE and TE end points (s_bar of Eq. 6-99) and profile family.  ThroughflowError unless
    /// Z > 0, the four LE / TE radii are > 0 and the family is known (or given by the options).
    explicit CorrelationBladeRow(const BladeRowSpec& row, CorrelationRowOptions options = {});
    [[nodiscard]] BladeRowDischarge evaluate(const BladeRowInlet& in) const override;
    /// The full Ch. 6 result for the same input: evaluate returns its beta2, loss and k_w.
    [[nodiscard]] CascadePerformance performance(const BladeRowInlet& in) const;

    // the mapping (public for the tests)
    [[nodiscard]] CascadeGeometry cascade_geometry(const BladeRowInlet& in) const;
    [[nodiscard]] static InletState inlet_state(const BladeRowInlet& in);
    [[nodiscard]] static DischargeEstimate discharge_estimate(const BladeRowInlet& in);
    /// R = P1 / (rho1 T1), gamma = a1^2 rho1 / P1.  ThroughflowError for a non-physical state.
    [[nodiscard]] static Gas gas_of(const BladeRowInlet& in);
    [[nodiscard]] double mean_pitch() const noexcept { return mean_pitch_; }  ///< s_bar [m]
    [[nodiscard]] ProfileFamily profile() const noexcept { return profile_; }
    [[nodiscard]] Camberline camberline() const noexcept { return camberline_; }

private:
    CorrelationRowOptions opt_;
    ProfileFamily profile_ = ProfileFamily::Naca;
    Camberline camberline_ = Camberline::Naca65;
    double mean_pitch_ = 0.0;
};

/// The seal of a shrouded row (Sec. 6.10, Fig. 6-22).
struct SealGeometry
{
    int fins = 0;               ///< N
    double fin_thickness = 0.0; ///< t [m]
    double fin_pitch = 0.0;     ///< p [m]
    std::optional<double> radius;  ///< r_seal [m]; empty: the mean of the seal wall's radii at QN1 and QN2 (ours)
};

struct CorrelationRowLossOptions
{
    EndWall seal_wall = EndWall::Hub;  ///< the wall of a shrouded row's band and seal (ours: a stator's inner band)
    SealGeometry seal;
};

/// Chapter 6's row-level losses (Secs. 6.9, 6.10: endwall_leakage) as Chapter 7's RowLossModel (see the file
/// comment).  An unshrouded row: the tip gap at the wall with a clearance (none: no loss); a shrouded row: the
/// seal on seal_wall, with the clearance of that wall.  ThroughflowError for a row with gaps at both walls or a
/// shrouded row with a clearance at its other wall; CorrelationError for what endwall_leakage refuses.
class CorrelationRowLoss final : public RowLossModel
{
public:
    explicit CorrelationRowLoss(CorrelationRowLossOptions options = {});
    [[nodiscard]] RowLossResult evaluate(const RowLossInput& in) const override;
    /// The details of the same evaluation (no gap: a default RowLeakage).
    [[nodiscard]] RowLeakage leakage(const RowLossInput& in) const;

private:
    CorrelationRowLossOptions opt_;
};

/// The Ch. 6 wall state of a through-flow wall state (absolute C_theta, book sign).
[[nodiscard]] Turbolenza::WallState correlation_wall_state(const WallState& w);

/// A row with CorrelationBladeRow and CorrelationRowLoss attached (the row's own data in their constructors).
[[nodiscard]] BladeRowSpec with_correlations(BladeRowSpec row, const CorrelationRowOptions& row_options = {},
                                             const CorrelationRowLossOptions& loss_options = {});

// ─────────────────────────────────────────────────────────────────────────────
// Chapter 8 -> Chapter 7: the end-wall blockage
// ─────────────────────────────────────────────────────────────────────────────

struct EndwallCouplingOptions
{
    EndwallInlet inlet;                ///< B, H1, m at QN 0 (Sec. 8.8)
    EndwallSettings settings;          ///< Ch. 8's model choices (Endwall.hpp)
    EndWall seal_wall = EndWall::Hub;  ///< the band of a shrouded row (as CorrelationRowLossOptions)
};

/// Chapter 8 as Chapter 7's EndwallModel (see the file comment).
class ChapterEightEndwall final : public EndwallModel
{
public:
    explicit ChapterEightEndwall(EndwallCouplingOptions options = {});
    void begin_pass(int pass, int n_qn) override;
    /// ThroughflowError for a call out of order (QN i without QN i - 1 in this pass); EndwallError from Ch. 8.
    [[nodiscard]] EndwallResult evaluate(const EndwallStation& st) override;

    /// The march of the current (last) pass, per quasi-normal: what the last call for it computed.
    struct Record
    {
        bool valid = false;
        EndwallWallStation hub_wall, shroud_wall;
        EndwallLayer hub, shroud;
        bool row = false;  ///< the step to this quasi-normal spanned a blade row
        EndwallRowWallResult hub_row, shroud_row;
        EndwallBlockage blockage;
    };
    [[nodiscard]] const std::vector<Record>& records() const noexcept { return rec_; }
    [[nodiscard]] int pass() const noexcept { return pass_; }

    // the mapping (public for the tests)
    [[nodiscard]] static EndwallWallStation wall_station(const EndwallEdge& e);
    [[nodiscard]] Turbolenza::EndwallRowWall row_wall(const EndwallRow& row, EndWall wall) const;

private:
    EndwallCouplingOptions opt_;
    int pass_ = -1;
    std::vector<Record> rec_;
};

}  // namespace BladeGen::Turbolenza::Throughflow
