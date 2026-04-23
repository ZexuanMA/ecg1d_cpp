#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include "tdvp_solver.hpp"
#include <vector>

namespace ecg1d {

// WARNING: despite the name, this function computes the Fourier transform
// of the single-particle position density |psi(x)|^2 (the one-body form
// factor / static structure factor), NOT the momentum probability density
// |psi_tilde(k)|^2. Integrates to 2*pi*|psi(0)|^2, not 1. For N=1 the true
// momentum probability density is provided by
// `ecg_momentum_density_1p` in src/verify/ecg_wavefunction.hpp. The correct
// N>=2 generalization (requires one-body reduced density matrix) is TODO.
//
// n(k) = sum_{ij} conj(u_i) u_j sum_p sign(p) * M_G * exp(-k^2/(4 K_aa)) * exp(+i k mu_a) / S
Eigen::VectorXd momentum_distribution(
    const std::vector<BasisParams>& basis,
    const PermutationSet& perms,
    const Eigen::VectorXd& k_points,
    int particle_index = 0);

// Compute energy components separately
struct EnergyComponents {
    double total;
    double kinetic;
    double harmonic;
    double interaction;  // delta + gaussian
    double kicking;
    double overlap;
};

EnergyComponents compute_energy_components(
    const std::vector<BasisParams>& basis,
    const HamiltonianTerms& terms);

} // namespace ecg1d
