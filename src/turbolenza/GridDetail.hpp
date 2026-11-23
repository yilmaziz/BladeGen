/// @file turbolenza/GridDetail.hpp
/// @brief Internals of the grid generation (Grid.hpp), declared for the grid sources and the
///        tests: the reference blade read as its C2 curves (v1 grid._CurveRow), the station cuts
///        (v1 _Cut), the turning points of a curve (v1 _turning_points) and NumPy helpers whose
///        exact semantics the port needs.
///
/// Turning points (v1 _turning_points, grid.py section 7).  v1 builds the u derivative of one
/// coordinate as a scipy BSpline of degree p - 1, converts it to piecewise polynomials
/// (PPoly.from_spline: on each knot span the Taylor coefficients at the span's left knot) and
/// takes PPoly.roots(extrapolate=False): every real root of each span polynomial on its closed
/// span, plus a span join where the two neighbouring pieces change sign (scipy's discontinuity
/// roots).  The port forms the same span polynomials, f(s) = sum_k X^(k+1)(t_l) / k! s^k,
/// s = u - t_l, from the curve's own derivatives at t_l, and finds their real roots on [0, h] by
/// isolation: the roots of f' (recursively, down to degree 1) split [0, h] into pieces on which
/// f is monotone, and each piece whose ends differ in sign holds exactly one root, bracketed and
/// bisected to round-off.  That finds every root of odd multiplicity (the only ones v1 keeps: a
/// root of even order is no turn) and any root at which f is exactly zero.  scipy's companion-
/// matrix roots can also report a nearly double root (|imaginary part| <= 1e-10) that the
/// isolation does not; v1 then drops it in its sign-change test, so the turning points agree.
/// The candidates then go through v1's own steps unchanged: 0 < u < 1, unique, merged within
/// _TURN_MERGE = 1e-12, kept where the derivative changes sign across +- _TURN_BRACKET = 1e-9,
/// and bisected 64 times on the curve's derivative.
#pragma once

#include "turbolenza/Curve.hpp"
#include "turbolenza/Grid.hpp"
#include "turbolenza/GridData.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza::detail {

inline constexpr int kBisections = 64;        ///< v1 _BISECTIONS
inline constexpr double kTurnBracket = 1e-9;  ///< v1 _TURN_BRACKET
inline constexpr double kTurnMerge = 1e-12;   ///< v1 _TURN_MERGE
inline constexpr int kSamples = 16;           ///< v1 _SAMPLES
inline constexpr int kZooms = 16;             ///< v1 _ZOOMS

// ── NumPy / Python semantics ─────────────────────────────────────────────────
/// np.sign: -1, 0 or 1.
[[nodiscard]] inline double np_sign(double x) noexcept { return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); }
/// np.mod for floats (the result takes the sign of b).
[[nodiscard]] double np_mod(double a, double b);
/// Python round() of a float: to nearest, ties to even.
[[nodiscard]] long long py_round(double x);
/// np.linspace(a, b, n): a + k (b - a) / (n - 1), the last one b exactly.
[[nodiscard]] std::vector<double> linspace(double a, double b, int n);
/// np.argsort(kind="stable").
[[nodiscard]] std::vector<std::size_t> stable_argsort(const std::vector<double>& v);
/// np.unique of doubles: sorted, exact duplicates removed.
[[nodiscard]] std::vector<double> unique(std::vector<double> v);
/// printf %.{precision}g (Python's format spec "g").
[[nodiscard]] std::string fmt_g(double x, int precision = 6);
[[nodiscard]] std::string fmt_e(double x, int precision);

// ── polynomials and turning points ───────────────────────────────────────────
/// Real roots in [0, h] of sum_k c[k] s^k (see the file comment): every root of odd
/// multiplicity, and every point where the polynomial evaluates to exactly zero at a
/// critical point or an end; ascending.  A polynomial that is identically zero has none.
[[nodiscard]] std::vector<double> polynomial_roots(const std::vector<double>& c, double h);
/// v1 _turning_points(curve, axis): u in (0, 1) where coordinate axis (0: m, 1: y) turns.
[[nodiscard]] std::vector<double> turning_points(const ParametricBladeCurve& curve, int axis);

// ── the station cuts (v1 _Cut) ───────────────────────────────────────────────
struct Cut
{
    std::vector<double> m;                     ///< the stations
    std::vector<int> count;                    ///< crossings per station
    std::vector<std::vector<int>> segment;     ///< per station: every branch, by ascending crossing height
    std::vector<std::vector<double>> height;   ///< the same order: crossing height, +inf past count
    std::vector<std::vector<double>> u;        ///< the same order: curve parameter, NaN past count
    std::vector<double> metal;                 ///< length of the reference blade's cut on the line
    double touch = 0.0;                        ///< TOUCH chords in the active unit

    [[nodiscard]] bool cuts_blade(std::size_t i) const { return metal[i] > touch; }
    /// numpy take_along_axis(segment, count - 1): the last column when count is 0
    [[nodiscard]] int top_segment(std::size_t i) const;
};

/// v1 _CurveRow: the reference blade of a C2 cascade as its two authoritative curves, cut into
/// branches monotone in m and y (grid.py section 7).
class CurveRow
{
public:
    explicit CurveRow(const PassageBoundaries& b);

    double touch = 0.0;
    std::shared_ptr<const SurfaceCurves> owner;           ///< keeps the curves alive
    std::array<const ParametricBladeCurve*, 2> curves{};  ///< suction, pressure (cascade frame)
    std::array<double, 2> lengths{};
    double total = 0.0;
    std::vector<int> curve_of;  ///< per branch: 0 suction, 1 pressure
    std::vector<double> u0, u1, m0, m1, lo_m, hi_m;
    double tau_tip = 0.0, tau_tail = 0.0, m_tip = 0.0, m_tail = 0.0;
    std::vector<bool> lower_contour;

    [[nodiscard]] int n_branches() const noexcept { return static_cast<int>(curve_of.size()); }
    [[nodiscard]] Surface surface_of(int branch) const noexcept
    {
        return curve_of[static_cast<std::size_t>(branch)] == 0 ? Surface::Suction : Surface::Pressure;
    }
    [[nodiscard]] std::array<double, 2> extremes() const noexcept { return {m_tip, m_tail}; }
    [[nodiscard]] double tau(int branch, double u) const;
    [[nodiscard]] Point point(int branch, double u) const;
    /// u where the branch crosses x = m (v1 solve: 64 bisections of [u0, u1]).
    [[nodiscard]] double solve(int branch, double m) const;
    /// y where the branch crosses x = m (v1 smooth).
    [[nodiscard]] double smooth(int branch, double m) const { return point(branch, solve(branch, m)).y; }
    /// GridError for an odd crossing count.
    [[nodiscard]] Cut cut(const std::vector<double>& m) const;
    /// v1 identity: the wall identity of stations first .. first + n - 1 of the cut (rows [first, last]).
    [[nodiscard]] WallIdentity identity(const Cut& cut, int first, int last) const;

    struct Pieces
    {
        std::vector<double> lo, hi;
        std::vector<int> top, bottom;
    };
    [[nodiscard]] const Pieces& pieces() const;
    /// v1 widest_cut: (span, m).
    [[nodiscard]] std::array<double, 2> widest_cut() const;
    /// v1 _Row.narrowest: (width, m).
    [[nodiscard]] std::array<double, 2> narrowest(double pitch) const;

private:
    mutable std::vector<Pieces> pieces_;  // cached (0 or 1 entry)
};

/// v1 Grid.__post_init__ on the given row: the topology, width and wall-identity checks
/// (GridError), then wall_identity and wall_d filled.  The nodes, m, eta, first/last_blade_i
/// and boundaries must be set.
void finish_grid(Grid& grid, const CurveRow& row);
/// The same with the cut of grid.m already made.
void finish_grid(Grid& grid, const CurveRow& row, const Cut& cut);

/// v1 _grid_on_stations.
[[nodiscard]] Grid grid_on_stations(const PassageBoundaries& b, const CurveRow& row, const std::vector<double>& m,
                                    int n_eta, const Cut* cut = nullptr);
/// v1 _check_passage.
void check_passage(const Cascade& cascade, const PassageBoundaries& b, const CurveRow& row);
/// v1 _stations.
[[nodiscard]] std::vector<double> stations(const PassageBoundaries& b, int n_passage);
/// v1 _mean_tan(a, b) for one a.
[[nodiscard]] double mean_tan(double a, double b);
/// v1 _side_line for one m: (beta, y).
[[nodiscard]] std::array<double, 2> side_line(double m, double fraction, const Point& junction,
                                              double junction_angle, double far_angle);

}  // namespace BladeGen::Turbolenza::detail
