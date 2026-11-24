/// @file turbolenza/Detail.hpp
/// @brief Internal helpers of the TURBOlenza solver port (not part of the API the window uses).
///
///  * the reductions of v1 in v1's order: NumPy's pairwise summation (numpy add.reduce on a
///    contiguous float64 array, numpy/_core/src/umath/loops_utils.h.src pairwise_sum), np.trapezoid
///    and np.median, so that a serial sum gives the same bits as v1 where the summands agree;
///  * Python's float repr (CSV and JSON output, error messages), math.radians / math.degrees;
///  * the thread count of the per-node loops: a thread-local number set by ThreadScope for the
///    duration of one solve, read by the OpenMP `num_threads` clause (TURBOLENZA_PARALLEL_FOR).
///    Only per-node maps are parallel (PORT.md rule 4): every node is computed by the same
///    instructions whatever thread runs it, so the answer does not depend on the thread count.
#pragma once

#include "turbolenza/Field.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza::detail {

[[nodiscard]] inline std::size_t sz(int k) noexcept { return static_cast<std::size_t>(k); }
[[nodiscard]] inline int isz(std::size_t k) noexcept { return static_cast<int>(k); }
[[nodiscard]] inline double nan() noexcept { return std::numeric_limits<double>::quiet_NaN(); }

/// math.radians: x * (pi / 180) (CPython's degToRad); math.degrees: x * (180 / pi).
[[nodiscard]] inline double radians(double deg) noexcept { return deg * (std::numbers::pi / 180.0); }
[[nodiscard]] inline double degrees(double rad) noexcept { return rad * (180.0 / std::numbers::pi); }

/// NumPy's pairwise sum of n contiguous values (the order of np.sum / add.reduce on float64).
[[nodiscard]] inline double pairwise_sum(const double* a, std::size_t n) noexcept
{
    if (n < 8) {
        double res = -0.0;
        for (std::size_t i = 0; i < n; ++i) res += a[i];
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (std::size_t j = 0; j < 8; ++j) r[j] = a[j];
        std::size_t i = 8;
        for (; i < n - (n % 8); i += 8)
            for (std::size_t j = 0; j < 8; ++j) r[j] += a[i + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i];
        return res;
    }
    std::size_t n2 = n / 2;
    n2 -= n2 % 8;
    return pairwise_sum(a, n2) + pairwise_sum(a + n2, n - n2);
}
[[nodiscard]] inline double pairwise_sum(const std::vector<double>& a) noexcept { return pairwise_sum(a.data(), a.size()); }

/// np.trapezoid(y, x): add.reduce(diff(x) * (y[1:] + y[:-1]) / 2.0).
[[nodiscard]] inline double trapezoid(const double* y, const std::vector<double>& x)
{
    std::vector<double> t(x.size() > 0 ? x.size() - 1 : 0);
    for (std::size_t k = 0; k < t.size(); ++k) t[k] = (x[k + 1] - x[k]) * (y[k + 1] + y[k]) / 2.0;
    return pairwise_sum(t);
}

/// np.median of a non-empty list: the middle value, or the mean of the two middle values.
[[nodiscard]] inline double median(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return (v[n / 2 - 1] + v[n / 2]) / 2.0;
}

/// Python's repr(float): the shortest round-trip digits, positional for 1e-4 <= |x| < 1e16
/// (with ".0" when integral), else d.ddde+XX; "nan", "inf", "-inf".
[[nodiscard]] inline std::string py_repr(double x)
{
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    char buf[64];
    const auto res = std::to_chars(buf, buf + sizeof(buf), x, std::chars_format::scientific);
    std::string s(buf, res.ptr);
    std::string sign;
    if (!s.empty() && s[0] == '-') {
        sign = "-";
        s.erase(0, 1);
    }
    const auto epos = s.find('e');
    std::string mant = s.substr(0, epos);
    const int exp10 = std::stoi(s.substr(epos + 1));
    std::string digits;
    for (char c : mant)
        if (c != '.') digits += c;
    const int decpt = exp10 + 1;  // digits d1 d2 ... with the point after decpt digits
    if (decpt > -4 && decpt <= 16) {
        std::string out;
        if (decpt <= 0) {
            out = "0." + std::string(sz(-decpt), '0') + digits;
        } else if (sz(decpt) >= digits.size()) {
            out = digits + std::string(sz(decpt) - digits.size(), '0') + ".0";
        } else {
            out = digits.substr(0, sz(decpt)) + "." + digits.substr(sz(decpt));
        }
        return sign + out;
    }
    std::string m = digits.substr(0, 1);
    if (digits.size() > 1) m += "." + digits.substr(1);
    const int e = decpt - 1;
    std::string es = std::to_string(e < 0 ? -e : e);
    if (es.size() < 2) es = "0" + es;
    return sign + m + "e" + (e < 0 ? "-" : "+") + es;
}

/// The thread count of the per-node loops of the current solve (1 outside a solve).
inline thread_local int t_threads = 1;
[[nodiscard]] inline int threads() noexcept { return t_threads; }

/// Sets the per-node loop thread count of the calling thread for one scope.
class ThreadScope
{
public:
    explicit ThreadScope(int n) noexcept : saved_(t_threads) { t_threads = n < 1 ? 1 : n; }
    ~ThreadScope() { t_threads = saved_; }
    ThreadScope(const ThreadScope&) = delete;
    ThreadScope& operator=(const ThreadScope&) = delete;

private:
    int saved_;
};

}  // namespace BladeGen::Turbolenza::detail

// A per-node map over a signed int index: `TURBOLENZA_PARALLEL_FOR for (int k = 0; k < n; ++k)`.
#if defined(TURBOLENZA_HAVE_OPENMP)
#if defined(_MSC_VER) && !defined(__clang__)
#define TURBOLENZA_PARALLEL_FOR                                                                                    \
    __pragma(omp parallel for schedule(static) num_threads(::BladeGen::Turbolenza::detail::threads()) if (::BladeGen::Turbolenza::detail::threads() > 1))
#else
#define TURBOLENZA_PARALLEL_FOR_STR(x) _Pragma(#x)
#define TURBOLENZA_PARALLEL_FOR                                                                                    \
    TURBOLENZA_PARALLEL_FOR_STR(omp parallel for schedule(static) num_threads(::BladeGen::Turbolenza::detail::threads()) if (::BladeGen::Turbolenza::detail::threads() > 1))
#endif
#else
#define TURBOLENZA_PARALLEL_FOR
#endif
