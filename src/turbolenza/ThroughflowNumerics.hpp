/// @file turbolenza/ThroughflowNumerics.hpp
/// @brief The numerical approximations of the meridional through-flow analysis: three-point
///        differences and integrals on non-uniform points, and the solution of the linear
///        normal-momentum ODE.
///
/// Aungier, Axial-Flow Compressors (ASME Press, 2003), Sec. 7.8 (pp. 171-172) and Sec. 7.2
/// (pp. 160-161), equations read from the page images:
///
///  * Eq. 7-67 (interior point 2 of three successive points 1, 2, 3):
///        f'_2 = [ (m3 - m2)(f2 - f1) / (m2 - m1) + (m2 - m1)(f3 - f2) / (m3 - m2) ] / (m3 - m1);
///  * Eq. 7-68 (first point):
///        f'_1 = [ (m3 - m1)(f2 - f1) / (m2 - m1) - (m2 - m1)(f3 - f1) / (m3 - m1) ] / (m3 - m2);
///  * Eq. 7-69 (last point):
///        f'_3 = [ (m3 - m1)(f3 - f2) / (m3 - m2) - (m3 - m2)(f3 - f1) / (m3 - m1) ] / (m2 - m1);
///  * Eq. 7-70 (two-point end rule, preferred by Aungier for the stream-surface slopes and
///    curvatures "for this application, which can often induce numerical stability problems"):
///        f'_1 = (f2 - f1) / (m2 - m1)  (and the mirror image at the last point; reading: the end
///        point takes the slope of its own pair);
///  * Eq. 7-71 (interior second derivative):
///        f''_2 = 2 / (y3 - y1) [ (f3 - f2) / (y3 - y2) - (f2 - f1) / (y2 - y1) ];
///  * Eq. 7-72 (integral over a segment whose right end, point 2, is interior):
///        int_{y1}^{y2} f dy = ( f2 - f'_2 (y2 - y1) / 2 + f''_2 (y2 - y1)^2 / 6 ) (y2 - y1);
///  * Eq. 7-73 (integral over the last segment, point 3 the last point, point 2 interior):
///        int_{y2}^{y3} f dy = ( f2 + f'_2 (y3 - y2) / 2 + f''_2 (y3 - y2)^2 / 6 ) (y3 - y2);
///    with f'_2 from Eq. 7-67 and f''_2 from Eq. 7-71; integrals from 0 to y_s are the sum over the
///    segments (p. 172).  All of them are exact for quadratics (the truncated Taylor series).
///  * Eqs. 7-34 to 7-36: dW/dy = f1(y) W + f4(y) has the solution
///        W(y) = W(0) F(y) + F(y) int_0^y f4 / F dy,   F(y) = exp( int_0^y f1 dy ),
///    the integrals by Eqs. 7-72 / 7-73 (ours: the same three-point integrals Aungier uses for
///    every integral in y).
///
/// The functions refuse (std::invalid_argument) fewer than three points, abscissae that are not
/// strictly increasing, mismatched lengths and non-finite values: nothing is extrapolated.
#pragma once

#include <vector>

namespace BladeGen::Turbolenza::Throughflow {

/// How the end points of a derivative are taken.
enum class EndRule
{
    ThreePoint,  ///< Eqs. 7-68 / 7-69
    TwoPoint,    ///< Eq. 7-70
};

/// df/dx at every point: Eq. 7-67 inside, Eq. 7-68 / 7-69 or Eq. 7-70 at the ends.
[[nodiscard]] std::vector<double> derivative(const std::vector<double>& x, const std::vector<double>& f,
                                             EndRule ends = EndRule::ThreePoint);

/// d2f/dx2 at the interior points (Eq. 7-71); the returned vector has the size of x with NaN at
/// both ends (Aungier gives no end formula and the integrals never use one).
[[nodiscard]] std::vector<double> second_derivative(const std::vector<double>& x, const std::vector<double>& f);

/// The integral of f over segment k = [x_k, x_(k+1)]: Eq. 7-72 with point 2 = k + 1 when k + 1 is
/// interior, Eq. 7-73 with point 2 = k for the last segment.
[[nodiscard]] double segment_integral(const std::vector<double>& x, const std::vector<double>& f, int k);

/// The integral of f from x_k to xi, xi in [x_k, x_(k+1)]: the same truncated Taylor series as
/// segment_integral, integrated to xi (so xi = x_(k+1) gives segment_integral exactly).  Used to
/// invert the mass-flow function of Eq. 7-54 inside a segment (ours).
[[nodiscard]] double partial_segment_integral(const std::vector<double>& x, const std::vector<double>& f, int k,
                                              double xi);

/// The running integral I_k = int_{x_0}^{x_k} f dx (I_0 = 0) by summation of the segments.
[[nodiscard]] std::vector<double> cumulative_integral(const std::vector<double>& x, const std::vector<double>& f);

/// The total integral from x_0 to x_(n-1).
[[nodiscard]] double integral(const std::vector<double>& x, const std::vector<double>& f);

/// The pieces of Eq. 7-35: W(y) = W(0) F(y) + G(y), G(y) = F(y) int_0^y f4 / F dy.
struct LinearOdeSolution
{
    std::vector<double> F;  ///< Eq. 7-36
    std::vector<double> G;
    [[nodiscard]] std::vector<double> evaluate(double w0) const;
};

/// Eqs. 7-35 / 7-36 on the points y (y_0 = 0 is not required: the integrals start at y_0).
[[nodiscard]] LinearOdeSolution solve_linear_ode(const std::vector<double>& y, const std::vector<double>& f1,
                                                 const std::vector<double>& f4);

}  // namespace BladeGen::Turbolenza::Throughflow
