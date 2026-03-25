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
// Returns the new u coefficients (basis functions unchanged).
void apply_analytic_kick(std::vector<BasisParams>& basis,
                         const PermutationSet& perms,
                         double kappa, double k_L,
                         int n_bessel = 20);

} // namespace ecg1d
