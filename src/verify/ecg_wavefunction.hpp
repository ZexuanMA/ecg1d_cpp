#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include <Eigen/Dense>
#include <vector>

namespace ecg1d {
namespace verify {

// Single-particle position-space density on a 1D grid for the ECG many-body state.
// For N=1: n(x) = |psi(x)|^2 / <psi|psi>.
// For N=2: marginal of |psi(x_1,x_2)|^2 over x_{1-particle_index}, normalized by S.
// Works for any N; uses Gaussian marginalization formula with K_eff = 1/(K^{-1})_{aa}.
// Output integrates to 1 for a bosonic symmetric state (per-particle probability).
Eigen::VectorXd ecg_single_particle_density(
    const std::vector<BasisParams>& basis,
    const PermutationSet& perms,
    const Eigen::VectorXd& x_points,
    int particle_index = 0);

// Evaluate the (single-particle, N=1) ECG wavefunction psi(x) on a grid.
//   psi(x) = sum_i u_i * exp(-(A_i + B_i) x^2 + 2 R_i B_i x - R_i B_i R_i).
// Output is NOT normalized to unity; its integral of |psi|^2 equals <psi|psi> = S.
// Throws runtime_error if called with N != 1.
Eigen::VectorXcd ecg_wavefunction_1p(
    const std::vector<BasisParams>& basis,
    const Eigen::VectorXd& x_points);

} // namespace verify
} // namespace ecg1d
