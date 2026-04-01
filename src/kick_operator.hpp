#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"

namespace ecg1d {

// Build the kick matrix K_ij = <phi_i| exp(-i*kappa*sum_a cos(2*k_L*x_a)) |phi_j>
// Uses Jacobi-Anger expansion: e^{-i*kappa*cos(theta)} = sum_n (-i)^n J_n(kappa) e^{i*n*theta}
// For Gaussian basis, <phi_i|e^{i*n*2*k_L*x_a}|phi_j> / M_G = exp(-n^2*k_L^2*K_inv(a,a) + i*n*2*k_L*mu(a))
//
// The total kick for N particles factorizes: kick_kernel = prod_a kick_kernel_a
// because cos(2*k_L*x_a) acts on different particle coordinates.
MatrixXcd build_kick_matrix(const std::vector<BasisParams>& basis,
                            const PermutationSet& perms,
                            double kappa, double k_L,
                            int n_bessel = 20);

// Apply an instantaneous kick to the wave function:
//   |psi'> = e^{-i*kappa*V_kick} |psi>
// This updates the linear coefficients u via: u' = S^{-1} K u
// where K is the kick matrix and S is the overlap matrix.
// Updates u in-place and returns the kick fidelity (norm_after/norm_before).
// Fidelity = 1.0 means perfect (unitary), >1 means basis can't represent kicked state.
double apply_analytic_kick(std::vector<BasisParams>& basis,
                           const PermutationSet& perms,
                           double kappa, double k_L,
                           int n_bessel = 20);

// Free evolution with fixed basis: solve i S du/dt = H u exactly.
// Uses eigendecomposition: u(t) = V exp(-i D t) V^{-1} u(0)
// where (H, S) generalized eigenproblem gives eigenvalues D and eigenvectors V.
// Only u coefficients change; A, B, R stay fixed.
void free_evolve_fixed_basis(std::vector<BasisParams>& basis,
                              const MatrixXcd& H, const MatrixXcd& S,
                              double T_duration);

// Augment a ground-state-optimized basis with momentum-carrying basis functions.
// For each existing basis function, creates copies with momenta p = n * 2 * k_L
// for n = -n_mom, ..., -1, +1, ..., +n_mom.
// The momentum is encoded via complex R: R_a → R_a + i*n*k_L (with B_aa set to b_val).
// Returns the augmented basis (original + momentum copies).
std::vector<BasisParams> augment_basis_with_momentum(
    const std::vector<BasisParams>& basis,
    double k_L, int n_mom = 2, double b_val = 0.5,
    double max_cond = 1e6);

} // namespace ecg1d
