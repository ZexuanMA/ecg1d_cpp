#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include "tdvp_solver.hpp"
#include <vector>

namespace ecg1d {

// Compute single-particle momentum distribution n(k) on a discrete k-grid.
// n(k) = sum_{ij} conj(u_i) u_j sum_p sign(p) * (Gaussian FT factor)
// Returns n(k) values at the given k_points.
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
