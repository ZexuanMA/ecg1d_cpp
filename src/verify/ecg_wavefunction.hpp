#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include <Eigen/Dense>
#include <vector>

namespace ecg1d {
namespace verify {

// Single-particle position-space density on a 1D grid for the ECG symmetric
// many-body state. Returns n(x) integrating to 1 (per-particle probability).
//
// For N=1: n(x) = |psi(x)|^2 / <psi|psi>.
// For N>=2: AVERAGES the particle-a marginals over a = 0..N-1, which is
// required to get the true bosonic-symmetric density when the basis functions
// themselves are not symmetric under particle swap (A_ii can differ across
// i). If a single a were used, the result would be the "half-symmetrized"
// <psi_L|delta(x-x_a)|psi_R>/<psi_L|psi_R>, which is wrong for ECG bases
// with off-diagonal or asymmetric A.
//
// The `particle_index` parameter is kept for ABI compatibility but ignored.
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

// Evaluate the momentum probability density n(k) = |psi_tilde(k)|^2 / <psi|psi>
// for the N=1 ECG wavefunction, where
//   psi_tilde(k) = (2*pi)^(-1/2) * integral psi(x) exp(-i k x) dx.
// Uses the analytic Gaussian Fourier transform of each basis function:
//   int exp(-alpha x^2 + beta x + gamma) dx = sqrt(pi/alpha) * exp(gamma + beta^2/(4 alpha))
// with alpha = A+B, beta = 2 R B - i k, gamma = -R B R.
// Output integrates to 1 (the true momentum probability density).
// This differs from the `momentum_distribution` function in
// src/observables.hpp, which computes the Fourier transform of the POSITION
// DENSITY, |psi|^2, i.e. the one-body form factor rather than the momentum
// probability density.
// Throws runtime_error if called with N != 1.
Eigen::VectorXd ecg_momentum_density_1p(
    const std::vector<BasisParams>& basis,
    const Eigen::VectorXd& k_points);

// Single-particle momentum probability density |phi_tilde(k)|^2 for the
// NON-INTERACTING bosonic ground state, extracted from the single-particle
// position density n(x). Valid when the exact ground state is a product
// state with a REAL NON-NEGATIVE single-particle orbital phi(x), which
// holds for the N=2 trap + kicking problem here (no inter-particle term).
// In that case phi(x) = sqrt(n(x)), and
//   |phi_tilde(k)|^2 = | (2 pi)^{-1/2} * integral phi(x) exp(-i k x) dx |^2.
// Output integrates to 1 on the supplied k-grid.
// Do NOT use for interacting systems (ground state may have off-diagonal
// long-range order that makes sqrt(n(x)) not a valid single-particle orbital).
Eigen::VectorXd momentum_density_from_density(
    const Eigen::VectorXd& x_points,
    const Eigen::VectorXd& n_x,
    const Eigen::VectorXd& k_points);

} // namespace verify
} // namespace ecg1d
