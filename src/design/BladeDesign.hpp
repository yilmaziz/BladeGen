/// @file design/BladeDesign.hpp
/// @brief The designer's blade: meridional configuration, five design layers (hub, 25 %,
///        mid, 75 %, tip), each with its own construction method, and everything the views
///        show (blade-to-blade loops, camber, metal angle, thickness, 3-D mesh, the
///        meridional metal-angle field).
///
/// Coordinates.  Every layer profile is built unstaggered, chord 1, LE point (0, 0), TE
/// point (1, 0), suction side +y (Profile2D and Naca65Section share this frame).  It is
/// placed in its layer's conformal plane (m', theta) with the stagger gamma (from the m'
/// axis toward +theta, Vos et al. 2025, Eq. 3) and the conformal chord c' fixed by the
/// meridional LE / TE positions: c' cos(gamma) = m'_TE - m'_LE (blade/Blade.hpp).  The
/// (m', theta) plane is the main design coordinate system; the views can also show
/// (m, r theta) (Vos et al. 2025, Sec. 2, option 2).  The conformal map to 3-D is the
/// cylindrical wrap on constant-radius layers and the map onto the layer's surface of
/// revolution otherwise (blade/FlowPath.hpp).
///
/// Re-stagger.  The whole blade can be re-staggered by an angle d (BladeDesign::restagger): it
/// turns by d about an axis, the radial line through a pivot at axial position x_a and angle
/// theta_a.  The pivot is a point of the hub, mid-span or tip section's chord line (chord fraction
/// x_p, 0 = LE; rotors are re-staggered about the hub LE, the default) or any (Z, theta) point.
/// Every placed section, and every mesh row, turns rigidly in its conformal plane by d about the
/// point where the axis crosses its stream surface, (m'(x_a), theta_a).  The sections stay on
/// their stream surfaces and their chord c' stays the nominal one, so the LE and TE leave the
/// meridional LE / TE lines.  The conformal map keeps angles, so every metal angle on every layer
/// changes by exactly d.  On a cylindrical layer this is the rigid rotation about the axis of
/// the unrolled (x, r theta) section.  A positive d closes the blade (higher stagger: more stall
/// margin), a negative one opens it (more flow).
///
/// Metal angle beta: the angle of the camber line from the meridional direction,
/// beta = atan(d theta / d m') = atan(r d theta / d m) (conformal: the true angle on the
/// stream surface).  Lengths in the unit of the meridional input (mm in the designer).
///
/// Mesh.  Blade rows between the layers interpolate the unstaggered loops point by point
/// (every layer is sampled with the same point layout) and the stagger / LE angle by a
/// natural cubic spline in span; spanwise rows are clustered toward hub and tip by a
/// symmetric two-sided tanh stretching (Vinokur 1983, J. Comput. Phys. 50:215; the
/// symmetric case u = 1/2 [1 + tanh(d (xi - 1/2)) / tanh(d / 2)]).
#pragma once

#include "blade/FlowPath.hpp"
#include "naca/Naca65.hpp"
#include "profile/Camber.hpp"
#include "profile/ProfileParams.hpp"
#include "profile/Quality.hpp"

#include <nlohmann/json.hpp>

#include <array>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace BladeGen::Design {

using Math::Vec2;
using Math::Vec3;

enum class Method { Naca65, Custom, Naca65Custom };

[[nodiscard]] std::string to_string(Method m);
[[nodiscard]] Method method_from_string(const std::string& s);

/// One design layer (a profile section on its stream surface).
struct Layer
{
    std::string name = "Mid";
    double span = 0.5;                       ///< 0 hub, 1 casing
    Method method = Method::Naca65;
    Naca::Naca65Params naca{};
    Profile::ProfileParams custom = Profile::ProfileParams::naca65_1210();
    std::optional<Naca::Naca65Params> fitted_to{};  ///< NACA 65 + Custom: the section the custom parameters were fitted to
    double fit_rms = 0.0;
    double stagger_deg = 30.0;               ///< chord angle from the m' axis toward +theta
    double theta_le_deg = 0.0;               ///< circumferential position of the LE point (lean)
    double solidity = 1.25;                  ///< 2-D cascade solidity (exports); the 3-D one follows from n_blades

    /// The profile used for the geometry: NACA 65 (section) or custom parameters.
    [[nodiscard]] bool uses_custom() const { return method == Method::Custom || (method == Method::Naca65Custom && fitted_to); }
};

/// Re-stagger of the whole blade (see the file comment).
struct Restagger
{
    enum class Axis { Hub, Mid, Tip, Point };
    double angle_deg = 0.0;       ///< + closes, - opens
    Axis axis = Axis::Hub;        ///< pivot: on the hub / mid-span / tip chord line, or the point below
    double chord_fraction = 0.0;  ///< pivot on the reference section's chord line (0 LE, 1 TE)
    double z = 0.0;               ///< Axis::Point: axial position of the radial axis [mm]
    double theta_deg = 0.0;       ///< Axis::Point: its angle [deg]
};
[[nodiscard]] std::string to_string(Restagger::Axis a);
[[nodiscard]] Restagger::Axis restagger_axis_from_string(const std::string& s);

struct MeshOptions
{
    int n_edge = 16;          ///< loop points per half edge connector
    int n_body = 48;          ///< loop points per side body
    int n_rows = 25;          ///< spanwise blade rows (hub .. tip)
    double row_clustering = 1.5;  ///< tanh stretching d (0 = uniform)
    int n_wall_m = 48;        ///< hub / casing grid points along the meridional direction
    int n_wall_theta = 14;    ///< hub / casing grid points across one pitch
    int n_field_chord = 41;   ///< metal-angle field: chordwise points
};

struct BladeDesign
{
    std::string name = "blade";
    Blade::MeridionalConfig meridional{};
    std::vector<Layer> layers;               ///< hub -> tip, spans increasing
    int n_blades = 18;
    bool mirror = false;                     ///< true: suction side toward -theta
    Restagger restagger{};                   ///< re-stagger of the whole blade
    MeshOptions mesh{};

    /// Default: five Custom layers (0, 25, 50, 75, 100 % span) with the NASA Rotor 67 profiles
    /// (the 14 TP-2879 sections fitted with the custom profile, blended in span), their stagger
    /// and LE angle, on the default meridional channel.  The profiles are chord-normalised, so
    /// they scale with the channel: c' cos(stagger) = m'_TE - m'_LE of each layer.  The NACA 65
    /// knobs of every layer hold the values of naca65_defaults().
    [[nodiscard]] static BladeDesign defaults();
    /// Five NACA 65 layers: camber and thickness fall and stagger rises from hub to tip.
    [[nodiscard]] static BladeDesign naca65_defaults();
    void validate() const;
};

[[nodiscard]] nlohmann::json to_json(const BladeDesign& d);
/// Reads "bladegen-design-2"; also "bladegen-design-1" (the v1.0.1 section designer: hub,
/// mid, tip; the 25 % / 75 % layers are blended from their neighbours).
[[nodiscard]] BladeDesign design_from_json(const nlohmann::json& j);
void write_design(const BladeDesign& d, const std::string& path);
[[nodiscard]] BladeDesign read_design(const std::string& path);

/// NASA Rotor 67 at span fraction s in [0, 1]: the custom profile, stagger and LE angle,
/// linearly blended between the two neighbouring fitted sections (exact at the 14 sections).
struct Rotor67Layer
{
    double span = 0.0, stagger_deg = 0.0, theta_le_deg = 0.0;
    Profile::ProfileParams profile{};
};
[[nodiscard]] Rotor67Layer rotor67_layer(double span);
[[nodiscard]] const std::vector<Rotor67Layer>& rotor67_sections();

/// The map of one section from its unstaggered frame (chord 1, LE (0, 0), TE (1, 0)) to the
/// conformal plane (m', theta), a similarity:
///   m'    = m0 + c' (x cos g - s y sin g),   theta = t0 + c' (x sin g + s y cos g)
/// with s = -1 when mirrored.  Nominal: g = gamma, c' = m'_TE / cos(gamma), (m0, t0) = (0, theta_LE).
/// Re-staggered by d about the conformal point P: g = gamma + d, the same c', and
///   (m0, t0) = P + R(d) ((0, theta_LE) - P),   R(d) the rotation by d (m' toward +theta),
/// so P stays where it was and every direction turns by d.
struct Placement
{
    double gamma = 0.0;  ///< chord angle from the m' axis [rad] (re-stagger included)
    double c = 1.0;      ///< conformal chord c'
    double m0 = 0.0;     ///< m' of the placed LE point
    double t0 = 0.0;     ///< theta of the placed LE point [rad]
    double sign = 1.0;   ///< -1: mirrored
    [[nodiscard]] Vec2 operator()(const Vec2& q) const;
};
[[nodiscard]] Placement placement(double stagger_deg, double theta_le_deg, double mprime_te, bool mirror,
                                  double restagger_deg = 0.0, const Vec2& pivot = Vec2::Zero());

/// The re-stagger axis of a design: the radial line at (x_a, theta_a [rad]).  Throws
/// std::invalid_argument when the reference section's pivot cannot be found.
[[nodiscard]] Vec2 restagger_axis(const BladeDesign& d, const Blade::FlowPath& fp);
/// m' where the stream surface reaches the axial position x (throws when it does not).
[[nodiscard]] double mprime_at_x(const Blade::StreamSurface& S, double x);

/// Natural cubic spline weights: f(s) = sum_k w_k f_k through (spans_k, f_k).
[[nodiscard]] std::vector<double> spline_weights(const std::vector<double>& spans, double s);
/// Symmetric two-sided tanh distribution on [0, 1], n points (d = 0: uniform).
[[nodiscard]] std::vector<double> tanh_two_sided(int n, double d);

// ─────────────────────────────────────────────────────────────────────────────
// Results
// ─────────────────────────────────────────────────────────────────────────────

/// 2-D analysis of one layer's profile (independent of its placement; cached).
struct Section2D
{
    Profile::QualityReport quality;   ///< sides LE -> TE, curvature, issues
    Profile::CurveSample loop;        ///< closed loop with the mesh layout (LE, SS, TE, PS)
    Profile::CamberLine camber;
    std::string designation;          ///< "NACA 65-(12A10)10 (kovach_sandercock)" / "Custom"
    double le_radius = 0.0, te_radius = 0.0, area = 0.0;
    std::optional<double> deviation_rms; ///< NACA 65 + Custom: RMS distance custom -> NACA sides
    /// NACA 65 sections: equivalent circular-arc chi_LE, chi_TE [deg] (Aungier Eq. 4-6).  The
    /// A10 mean line's own slope is logarithmic at its ends, so its LE / TE tangent is steeper.
    std::optional<std::array<double, 2>> equivalent_arc_deg;
    /// Custom layers: the B-spline control polygons (LE connector, SS, TE connector, PS; the
    /// sides LE -> TE) and the thickness law t(xi) at its SS stations.
    struct Control
    {
        std::vector<Vec2> le, ss, te, ps;
    };
    std::optional<Control> control;
    std::optional<Profile::ThicknessSample> thickness_law;
    double xi_max = 0.0;
};

/// One layer placed on its stream surface.
struct LayerResult
{
    std::string error;                ///< not empty: the layer could not be built
    std::shared_ptr<const Section2D> section;
    Blade::StreamSurface surface;
    double gamma = 0.0;               ///< stagger [rad] (re-stagger included)
    double chord_mprime = 0.0;        ///< conformal chord c'
    double mprime_le = 0.0;           ///< m' of the placed LE point (0 unless re-staggered)
    double theta_le = 0.0;            ///< theta of the placed LE point [rad]
    double r_le = 0.0, r_te = 0.0;
    double chord = 0.0;               ///< physical chord (3-D distance LE -> TE point)
    double pitch_le = 0.0;            ///< 2 pi r_LE / n_blades
    double solidity_3d = 0.0;         ///< chord / mean pitch
    std::vector<Vec2> conformal;      ///< loop in (m', theta)
    std::vector<Vec3> xyz;            ///< loop in Cartesian
    // camber on the stream surface
    std::vector<double> mprime;       ///< m' of the camber points (from the placed LE point)
    std::vector<double> m;            ///< meridional arc length of the camber points (from the placed LE point)
    std::vector<double> mfrac;        ///< meridional chord fraction (0 placed LE, 1 placed TE)
    std::vector<double> theta_deg;    ///< theta of the camber points
    std::vector<double> beta_deg;     ///< metal angle along the camber line
    std::vector<double> thickness;    ///< t / c along the camber line (chord-fraction x below)
    std::vector<double> xc;           ///< chord fraction of the camber points (unstaggered x)
    std::vector<Vec2> camber_xr;      ///< meridional (x, r) of the camber points
    double beta_le = 0.0, beta_te = 0.0, theta_te_deg = 0.0;
};

/// Structured surface grid, i fastest: p[i + ni * j].
struct Grid
{
    int ni = 0, nj = 0;
    std::vector<Vec3> p;
    [[nodiscard]] const Vec3& at(int i, int j) const { return p[static_cast<std::size_t>(i + ni * j)]; }
};

/// Scalar field on the meridional plane, i chordwise, j spanwise.
struct MeridionalField
{
    int ni = 0, nj = 0;
    std::vector<Vec2> xr;
    std::vector<double> value;
    double lo = 0.0, hi = 0.0;
};

struct DesignResult
{
    std::vector<LayerResult> layers;
    std::vector<std::string> errors;           ///< design-level problems (meridional input, ...)
    // meridional view
    std::vector<Vec2> hub, casing;             ///< wall polylines (inlet end -> outlet end)
    std::vector<std::vector<Vec2>> streamlines;///< one per layer
    std::vector<Vec2> le_line, te_line;
    // 3-D mesh (empty when a layer failed)
    Grid blade;                                ///< i around the loop, j hub -> tip
    Grid hub_wall, casing_wall;                ///< i meridional, j across one pitch
    std::vector<double> row_span;
    MeridionalField beta;                      ///< metal angle [deg] on the meridional plane
    MeridionalField thickness;                 ///< t / c on the meridional plane
    double seconds = 0.0;

    [[nodiscard]] bool ok() const;
};

/// Evaluates designs; keeps the 2-D sections of recent layer profiles (NACA 65 sections
/// take ~0.2 s to build), so a change in one layer recomputes only that layer.  The
/// uncached layers are built in parallel.  Thread-safe.
class Evaluator
{
public:
    [[nodiscard]] DesignResult evaluate(const BladeDesign& d);
    [[nodiscard]] std::shared_ptr<const Section2D> section(const Layer& layer, const MeshOptions& mesh);
    [[nodiscard]] std::shared_ptr<const Naca::Naca65Section> naca(const Naca::Naca65Params& p);
    void clear();

private:
    std::mutex mu_;
    std::map<std::string, std::shared_ptr<const Section2D>> sections_;
    std::map<std::string, std::shared_ptr<const Naca::Naca65Section>> naca_;
    std::vector<std::string> order_;  // insertion order of sections_ (bounded cache)
};

/// Builds the 2-D section of one layer (no cache).
[[nodiscard]] Section2D build_section2d(const Layer& layer, const Naca::Naca65Section* naca, const MeshOptions& mesh);

/// Places one layer (2-D section already built) on its stream surface; re-staggered by
/// `restagger_deg` about the radial axis at (x_a, theta_a) = `axis` (restagger_axis).
[[nodiscard]] LayerResult place_layer(const Layer& layer, std::shared_ptr<const Section2D> sec,
                                      const Blade::FlowPath& fp, int n_blades, bool mirror,
                                      double restagger_deg = 0.0, const Vec2& axis = Vec2::Zero());

/// Fits the custom parameters of a NACA 65 + Custom layer to its NACA section (seconds).
void fit_layer(Layer& layer, Evaluator& ev);

/// The blade mesh rows as closed loops (same layout), e.g. for the OCCT solid / STEP export.
[[nodiscard]] std::vector<std::vector<Vec3>> section_loops(const DesignResult& r);

}  // namespace BladeGen::Design
