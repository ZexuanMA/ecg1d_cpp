#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "tdvp_solver.hpp"
#include "permutation.hpp"
#include <vector>
#include <random>
#include <tuple>

namespace ecg1d {

// Compute H and S matrix elements for a basis pair (i,j) without u factors.
// Sums over permutations, applies physical constant coefficients.
std::pair<Cd, Cd> compute_HS_ij(const BasisParams& bi, const BasisParams& bj,
                                 const PermutationSet& perms,
                                 const HamiltonianTerms& terms);

// Build full H and S matrices (exploiting Hermitian symmetry)
std::pair<MatrixXcd, MatrixXcd> build_HS(const std::vector<BasisParams>& basis,
                                          const PermutationSet& perms,
                                          const HamiltonianTerms& terms);

// Result from generalized eigenvalue solve: energy + optional variance
struct EigenResult {
    double energy = std::numeric_limits<double>::infinity();
    double variance = -1.0;   // σ² = <H²>/<ψ|ψ> - E², negative = not computed
    int n_kept = 0;           // number of S eigenvectors kept after truncation
};

// Solve generalized eigenvalue problem via S^{-1/2} transformation.
// Returns lowest eigenvalue, or +inf if ill-conditioned.
// Uses eigenvalue truncation: discards S eigenvectors with w_i < w_max * rcond_trunc.
// max_cond: maximum allowed condition number after truncation (safety check)
// E_lower_bound: reject energies below this (ghost state protection)
// rcond_trunc: relative cutoff for S eigenvalue truncation (default 1e-10)
double lowest_energy(const MatrixXcd& H, const MatrixXcd& S,
                     double max_cond = 1e10,
                     double E_lower_bound = -1e10,
                     double rcond_trunc = 1e-10);

// Same as lowest_energy but also returns variance diagnostic
EigenResult lowest_energy_full(const MatrixXcd& H, const MatrixXcd& S,
                                double max_cond = 1e10,
                                double E_lower_bound = -1e10,
                                double rcond_trunc = 1e-10);

// Set u values from ground-state eigenvector of the generalized eigenvalue problem.
void set_u_from_eigenvector(std::vector<BasisParams>& basis,
                            const MatrixXcd& H, const MatrixXcd& S,
                            double max_cond = 1e10,
                            double rcond_trunc = 1e-10);

// Check if a trial basis function has excessive overlap with existing basis.
// Returns true if max normalized overlap > threshold (i.e., should reject).
bool has_excessive_overlap(const MatrixXcd& S, int k, double threshold = 0.99);

// Check if overlap matrix S is well-conditioned (w_min/w_max > 1/max_cond).
// This is the primary guard against ghost states: reject any trial that would
// make S too ill-conditioned, preventing the S^{-1/2} amplification of noise.
// Returns true if safe, false if should reject.
// Also returns w_min via optional pointer for diagnostics.
bool s_well_conditioned(const MatrixXcd& S, double max_cond = 1e10, double* w_min_out = nullptr);

// Physics-informed random basis for N=2 (CM/relative coordinate parameterization)
BasisParams random_basis_2particle(std::mt19937_64& rng, int name = 0);

// Small perturbation of existing basis function
BasisParams perturb_basis(const BasisParams& base, std::mt19937_64& rng,
                          double scale, int N);

// SVM Phase 1: Greedy basis construction
struct SvmResult {
    std::vector<BasisParams> basis;
    MatrixXcd H;
    MatrixXcd S;
};

SvmResult svm_build_basis(int N, int K_max, int n_trials,
                           const HamiltonianTerms& terms,
                           int seed = 42,
                           double E_lower_bound = 0.0);

// SVM Phase 2: Stochastic refinement
SvmResult stochastic_refine(std::vector<BasisParams> basis,
                             MatrixXcd H, MatrixXcd S,
                             const PermutationSet& perms,
                             int N, const HamiltonianTerms& terms,
                             int n_trials = 500, int max_rounds = 30,
                             int seed = 123,
                             double E_lower_bound = 0.0);

} // namespace ecg1d
