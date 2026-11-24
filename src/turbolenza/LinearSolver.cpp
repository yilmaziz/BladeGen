// Port of v1 linear_solver.py (see LinearSolver.hpp): solve, residual_norm, solve_linear_system,
// _checked_matrix, _checked_rhs; the factorization is Eigen SparseLU or MKL PARDISO.
#include "turbolenza/LinearSolver.hpp"

#include "turbolenza/Detail.hpp"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseLU>

#ifdef TURBOLENZA_HAVE_OPENMP
#include <omp.h>
#endif
#ifdef TURBOLENZA_HAVE_PARDISO
#include <mkl.h>
#include <oneapi/tbb/global_control.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <string>

namespace BladeGen::Turbolenza {

using detail::sz;

const char* linear_solver_name(LinearSolverKind kind) noexcept
{
    return kind == LinearSolverKind::Pardiso ? "PARDISO (Intel MKL)" : "SparseLU (Eigen, COLAMD)";
}

namespace {

void checked_matrix(const SparseMatrix& F)
{
    if (F.rows() != F.cols())
        throw LinearSolveError(
            std::format("the matrix is ({}, {}): a linear system needs a square matrix", F.rows(), F.cols()));
    if (F.rows() == 0) throw LinearSolveError("the matrix is empty");
    long long bad = 0;
    for (int c = 0; c < F.outerSize(); ++c)
        for (SparseMatrix::InnerIterator it(F, c); it; ++it)
            if (!std::isfinite(it.value())) ++bad;
    if (bad) throw LinearSolveError(std::format("the matrix has {} non-finite entries", bad));
}

void checked_rhs(const std::vector<double>& rhs, long long n, const char* name)
{
    if (static_cast<long long>(rhs.size()) != n)
        throw LinearSolveError(std::format("the {} has shape ({},), the matrix needs ({},)", name, rhs.size(), n));
    long long bad = 0;
    for (double v : rhs)
        if (!std::isfinite(v)) ++bad;
    if (bad) throw LinearSolveError(std::format("the {} has {} non-finite values", name, bad));
}

// Eigen's dense kernels inside SparseLU parallelize with OpenMP over omp_get_max_threads(); hold
// them to one thread so that the factorization is the same whatever the solve's thread count.
class SerialEigen
{
public:
    SerialEigen() noexcept
    {
#ifdef TURBOLENZA_HAVE_OPENMP
        saved_ = omp_get_max_threads();
        omp_set_num_threads(1);
#endif
    }
    ~SerialEigen()
    {
#ifdef TURBOLENZA_HAVE_OPENMP
        omp_set_num_threads(saved_);
#endif
    }
    SerialEigen(const SerialEigen&) = delete;
    SerialEigen& operator=(const SerialEigen&) = delete;

private:
    int saved_ = 1;
};

std::vector<double> sparse_lu(const SparseMatrix& F, const std::vector<double>& Q)
{
    const SerialEigen serial;
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> lu;
    lu.analyzePattern(F);
    lu.factorize(F);
    if (lu.info() != Eigen::Success)
        throw LinearSolveError(std::format("the {} x {} matrix is singular: SparseLU reports '{}'. Nothing is "
                                           "regularized or pinned: the equations do not determine Psi",
                                           F.rows(), F.cols(), lu.lastErrorMessage()));
    const Eigen::Map<const Eigen::VectorXd> b(Q.data(), static_cast<Eigen::Index>(Q.size()));
    const Eigen::VectorXd x = lu.solve(b);
    if (lu.info() != Eigen::Success)
        throw LinearSolveError(std::format("the {} x {} matrix is singular: SparseLU reports '{}'. Nothing is "
                                           "regularized or pinned: the equations do not determine Psi",
                                           F.rows(), F.cols(), lu.lastErrorMessage()));
    return std::vector<double>(x.data(), x.data() + x.size());
}

#ifdef TURBOLENZA_HAVE_PARDISO
std::vector<double> pardiso_lu(const SparseMatrix& F, const std::vector<double>& Q, int threads)
{
    Eigen::SparseMatrix<double, Eigen::RowMajor, MKL_INT> A = F;  // CSR, sorted columns
    A.makeCompressed();
    MKL_INT n = static_cast<MKL_INT>(A.rows());
    MKL_INT mtype = 11, nrhs = 1, maxfct = 1, mnum = 1, msglvl = 0, error = 0, phase = 13;
    void* pt[64] = {};
    MKL_INT iparm[64] = {};
    pardisoinit(pt, &mtype, iparm);
    iparm[34] = 1;  // zero-based ia, ja
    std::vector<MKL_INT> perm(sz(static_cast<int>(n)), 0);
    std::vector<double> b(Q), x(Q.size(), 0.0);
    {
        const oneapi::tbb::global_control limit(oneapi::tbb::global_control::max_allowed_parallelism,
                                                static_cast<std::size_t>(threads < 1 ? 1 : threads));
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, A.valuePtr(), A.outerIndexPtr(), A.innerIndexPtr(),
                perm.data(), &nrhs, iparm, &msglvl, b.data(), x.data(), &error);
        MKL_INT release = -1, release_error = 0;
        pardiso(pt, &maxfct, &mnum, &mtype, &release, &n, nullptr, A.outerIndexPtr(), A.innerIndexPtr(), perm.data(),
                &nrhs, iparm, &msglvl, nullptr, nullptr, &release_error);
    }
    if (error != 0)
        throw LinearSolveError(std::format("the {} x {} matrix is singular or cannot be factorized: PARDISO reports "
                                           "error {}. Nothing is regularized or pinned: the equations do not determine "
                                           "Psi",
                                           n, n, static_cast<long long>(error)));
    if (iparm[13] > 0)
        throw LinearSolveError(std::format("the {} x {} matrix is numerically singular: PARDISO perturbed {} pivot(s). "
                                           "Nothing is regularized or pinned: the equations do not determine Psi",
                                           n, n, static_cast<long long>(iparm[13])));
    return x;
}
#endif

// r = F x - Q row by row in increasing column order (the order of SciPy's CSR product).
std::vector<double> residual(const SparseMatrix& F, const std::vector<double>& x, const std::vector<double>& Q)
{
    const Eigen::SparseMatrix<double, Eigen::RowMajor, int> R = F;
    std::vector<double> r(Q.size(), 0.0);
    for (int i = 0; i < R.outerSize(); ++i) {
        double s = 0.0;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor, int>::InnerIterator it(R, i); it; ++it)
            s += it.value() * x[sz(static_cast<int>(it.col()))];
        r[sz(i)] = s - Q[sz(i)];
    }
    return r;
}

}  // namespace

std::vector<double> solve(const SparseMatrix& matrix, const std::vector<double>& rhs, LinearSolverKind kind, int threads)
{
    checked_matrix(matrix);
    checked_rhs(rhs, matrix.rows(), "rhs");
    std::vector<double> x;
    if (kind == LinearSolverKind::Pardiso) {
#ifdef TURBOLENZA_HAVE_PARDISO
        x = pardiso_lu(matrix, rhs, threads);
#else
        (void)threads;
        throw LinearSolveError("PARDISO was asked for, but this build has no Intel MKL (pardiso_available() is false)");
#endif
    } else {
        x = sparse_lu(matrix, rhs);
    }
    if (static_cast<long long>(x.size()) != matrix.rows())
        throw LinearSolveError(std::format("the solve returned shape ({},), expected ({},)", x.size(), matrix.rows()));
    long long bad = 0, first = -1;
    for (std::size_t k = 0; k < x.size(); ++k)
        if (!std::isfinite(x[k])) {
            if (first < 0) first = static_cast<long long>(k);
            ++bad;
        }
    if (bad)
        throw LinearSolveError(std::format("the solution is not finite at {} of {} unknowns (first: {}, {}): overflow "
                                           "or a numerically singular matrix",
                                           bad, x.size(), first, detail::py_repr(x[static_cast<std::size_t>(first)])));
    return x;
}

double residual_norm(const SparseMatrix& matrix, const std::vector<double>& solution, const std::vector<double>& rhs)
{
    checked_matrix(matrix);
    checked_rhs(rhs, matrix.rows(), "rhs");
    checked_rhs(solution, matrix.cols(), "solution");
    double m = 0.0;
    for (double v : residual(matrix, solution, rhs)) m = std::max(m, std::abs(v));
    return m;
}

LinearSolveResult solve_linear_system(const AssemblyResult& assembly, LinearSolverKind kind, int threads)
{
    const SparseMatrix& F = assembly.matrix;
    const std::vector<double>& Q = assembly.rhs;
    LinearSolveResult out;
    out.active_solution = solve(F, Q, kind, threads);
    const std::vector<double>& x = out.active_solution;
    try {
        out.psi = reconstruct(assembly, x);
    } catch (const MatrixAssemblyError& error) {
        throw LinearSolveError(std::string("the solved unknowns do not rebuild psi: ") + error.what());
    }
    for (double v : out.psi.v)
        if (!std::isfinite(v))
            throw LinearSolveError(std::format("the reconstructed psi has shape ({}, {}) (grid ({}, {})) or is not finite",
                                               out.psi.n_m, out.psi.n_eta, assembly.n_m, assembly.n_eta));
    const std::vector<double> r = residual(F, x, Q);
    double r_inf = 0.0, x_inf = 0.0, q_inf = 0.0;
    std::vector<double> r2(r.size()), q2(Q.size());
    for (std::size_t k = 0; k < r.size(); ++k) {
        r_inf = std::max(r_inf, std::abs(r[k]));
        x_inf = std::max(x_inf, std::abs(x[k]));
        q_inf = std::max(q_inf, std::abs(Q[k]));
        r2[k] = r[k] * r[k];
        q2[k] = Q[k] * Q[k];
    }
    const double r_l2 = std::sqrt(detail::pairwise_sum(r2)), q_l2 = std::sqrt(detail::pairwise_sum(q2));
    double F_inf = 0.0;
    {
        const Eigen::SparseMatrix<double, Eigen::RowMajor, int> R = F;
        for (int i = 0; i < R.outerSize(); ++i) {
            double s = 0.0;
            for (Eigen::SparseMatrix<double, Eigen::RowMajor, int>::InnerIterator it(R, i); it; ++it)
                s += std::abs(it.value());
            F_inf = std::max(F_inf, s);
        }
    }
    const double scale = F_inf * x_inf + q_inf;
    double implied = 0.0;
    for (double v : implied_residuals(assembly, out.psi)) implied = std::max(implied, std::abs(v));
    out.residual_inf = r_inf;
    out.residual_l2 = r_l2;
    out.relative_residual = q_l2 > 0.0 ? r_l2 / q_l2 : r_l2;
    out.backward_error = scale > 0.0 ? r_inf / scale : 0.0;
    out.implied_residual = implied;
    out.rows = static_cast<int>(F.rows());
    out.cols = static_cast<int>(F.cols());
    out.nnz = static_cast<long long>(F.nonZeros());
    out.finite = true;
    out.psi_active_min = *std::min_element(x.begin(), x.end());
    out.psi_active_max = *std::max_element(x.begin(), x.end());
    return out;
}

}  // namespace BladeGen::Turbolenza
