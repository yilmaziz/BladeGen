/// @file turbolenza/LinearSolver.hpp
/// @brief Solves the assembled sparse system F Psi = Q and rebuilds the whole psi field.
///
/// Port of v1 linear_solver.py (TURBOlenza v1 commit 67599a6; the system of Aungier 2003 Sec. 5.4
/// and Aungier 2000 Ch. 13, unchanged; only the linear algebra is ours).  v1 used SciPy's spsolve
/// (SuperLU, COLAMD, partial pivoting).  v2 offers two sparse direct LU factorizations of the same
/// F and Q, with no iterative method, row scaling of the assembled system, regularization or
/// pinned node:
///  * LinearSolverKind::SparseLU: Eigen::SparseLU with COLAMD ordering (serial; Eigen's own dense
///    kernels are held to one thread, so the answer does not depend on the thread count);
///  * LinearSolverKind::Pardiso: Intel MKL PARDISO (real nonsymmetric, mtype 11, its default
///    iparm: METIS ordering, its internal scaling and matching, pivot perturbation 1e-13,
///    iterative refinement), on MKL's TBB threading layer limited to the solve's thread count
///    by oneapi::tbb::global_control.  Only when built with TURBOLENZA_HAVE_PARDISO.
///    A different factorization: it agrees with SparseLU to round-off, not bit for bit.
/// Failures (LinearSolveError): not square or empty; sizes that disagree; non-finite entries; a
/// singular matrix (SparseLU: a zero pivot; PARDISO: an error code or a perturbed pivot); a
/// non-finite solution; a reconstructed psi that is not finite.  A large residual is reported,
/// not refused.
#pragma once

#include "turbolenza/Case.hpp"
#include "turbolenza/Field.hpp"
#include "turbolenza/MatrixAssembly.hpp"

#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 LinearSolveError.
class LinearSolveError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 LinearSolveResult.
struct LinearSolveResult
{
    std::vector<double> active_solution;  ///< Psi, in matrix_assembly's active order
    Field psi;                            ///< (n_m, n_eta), the far stations rebuilt
    double residual_inf = 0.0;            ///< ||F Psi - Q||_inf
    double residual_l2 = 0.0;
    double relative_residual = 0.0;       ///< ||r||_2 / ||Q||_2 (||r||_2 if Q = 0)
    double backward_error = 0.0;          ///< ||r||_inf / (||F||_inf ||Psi||_inf + ||Q||_inf)
    double implied_residual = 0.0;        ///< max |far-station jump residual| on psi
    int rows = 0, cols = 0;
    long long nnz = 0;
    bool finite = true;
    double psi_active_min = 0.0, psi_active_max = 0.0;
};

/// The name of a solver kind for reports ("SparseLU (Eigen, COLAMD)", "PARDISO (Intel MKL)").
[[nodiscard]] const char* linear_solver_name(LinearSolverKind kind) noexcept;

/// v1 solve: Psi of matrix Psi = rhs, checked before and after.  `threads` is the PARDISO
/// thread count (SparseLU is serial).
[[nodiscard]] std::vector<double> solve(const SparseMatrix& matrix, const std::vector<double>& rhs,
                                        LinearSolverKind kind = LinearSolverKind::SparseLU, int threads = 1);

/// v1 residual_norm: ||matrix solution - rhs||_inf.
[[nodiscard]] double residual_norm(const SparseMatrix& matrix, const std::vector<double>& solution,
                                   const std::vector<double>& rhs);

/// v1 solve_linear_system: the solve, reconstruct, and the diagnostics.
[[nodiscard]] LinearSolveResult solve_linear_system(const AssemblyResult& assembly,
                                                    LinearSolverKind kind = LinearSolverKind::SparseLU,
                                                    int threads = 1);

}  // namespace BladeGen::Turbolenza
