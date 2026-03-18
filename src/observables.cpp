#include "observables.hpp"
#include "hamiltonian.hpp"
#include "pair_cache.hpp"
#include "tdvp_solver.hpp"
#include <cmath>
#include <Eigen/Dense>

namespace ecg1d {

Eigen::VectorXd momentum_distribution(
    const std::vector<BasisParams>& basis,
    const PermutationSet& perms,
    const Eigen::VectorXd& k_points,
    int particle_index) {

    int K = static_cast<int>(basis.size());
    int N = basis[0].N();
    int n_k = static_cast<int>(k_points.size());

    // Compute overlap for normalization
    Cd S = overlap(basis);

    Eigen::VectorXd nk = Eigen::VectorXd::Zero(n_k);

    for (int ki = 0; ki < n_k; ki++) {
        double k = k_points(ki);
        Cd nk_val(0, 0);

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                Cd u_prod = std::conj(basis[i].u) * basis[j].u;

                for (int p = 0; p < perms.SN; p++) {
                    double sign = static_cast<double>(perms.signs[p]);
                    PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);

                    // For particle 'a' (particle_index), the momentum space
                    // density involves the Fourier transform of the Gaussian.
                    // n(k) contribution = sign * M_G * exp(-k^2 / (4 * K_aa))
                    //                     * exp(i * k * mu_a) / sqrt(K_aa)
                    // where K_aa is the (a,a) element of the K matrix and
                    // mu_a is the a-th component of mu = 0.5 * K^{-1} * b.
                    int a = particle_index;
                    Cd K_aa = c.K(a, a);
                    Cd mu_a = c.mu(a);

                    // Gaussian FT: integral of exp(-K_aa * x^2 + b_a * x - i*k*x)
                    // = sqrt(pi/K_aa) * exp((b_a - ik)^2 / (4*K_aa)) * M_G / sqrt(pi/K_aa)
                    // Simplified: M_G * exp(-k^2/(4*K_aa) + i*k*mu_a) * correction
                    // The ratio to M_G is exp(-k^2/(4*K_aa)) * exp(i*k*mu_a)
                    Cd phase = std::exp(Cd(0, 1) * k * mu_a);
                    Cd gauss_factor = std::exp(-k * k / (4.0 * K_aa));
                    Cd contrib = sign * c.M_G * gauss_factor * phase;

                    nk_val += u_prod * contrib;
                }
            }
        }

        // Normalize by overlap
        nk(ki) = (nk_val / S).real();
    }

    return nk;
}

EnergyComponents compute_energy_components(
    const std::vector<BasisParams>& basis,
    const HamiltonianTerms& terms) {

    int N = basis[0].N();
    Cd S = overlap(basis);
    EnergyComponents ec;
    ec.overlap = S.real();

    ec.kinetic = 0.0;
    ec.harmonic = 0.0;
    ec.interaction = 0.0;
    ec.kicking = 0.0;

    if (terms.kinetic) {
        Cd T = kinetic_energy_functional(basis);
        ec.kinetic = (T / S).real();
    }
    if (terms.harmonic) {
        Cd V = Harmonic_functional(basis);
        ec.harmonic = (V / S).real();
    }
    if (terms.delta && N >= 2) {
        Cd D = Delta_contact_functional(basis);
        ec.interaction += (D / S).real();
    }
    if (terms.gaussian && N >= 2) {
        Cd G = Gaussian_interaction_functional(basis);
        ec.interaction += (G / S).real();
    }
    if (terms.kicking) {
        Cd K = kicking_term_functional(basis);
        ec.kicking = (K / S).real();
    }

    ec.total = ec.kinetic + ec.harmonic + ec.interaction + ec.kicking;
    return ec;
}

} // namespace ecg1d
