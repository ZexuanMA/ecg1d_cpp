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

// Solve generalized eigenvalue problem via S^{-1/2} transformation.
// Returns lowest eigenvalue, or +inf if ill-conditioned.
double lowest_energy(const MatrixXcd& H, const MatrixXcd& S);

// Set u values from ground-state eigenvector of the generalized eigenvalue problem.
void set_u_from_eigenvector(std::vector<BasisParams>& basis,
                            const MatrixXcd& H, const MatrixXcd& S);

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
                           int seed = 42);

// SVM Phase 2: Stochastic refinement
SvmResult stochastic_refine(std::vector<BasisParams> basis,
                             MatrixXcd H, MatrixXcd S,
                             const PermutationSet& perms,
                             int N, const HamiltonianTerms& terms,
                             int n_trials = 500, int max_rounds = 30,
                             int seed = 123);

} // namespace ecg1d
