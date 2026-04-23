#include "ecg_wavefunction.hpp"
#include "pair_cache.hpp"
#include "hamiltonian.hpp"
#include <cmath>
#include <stdexcept>

namespace ecg1d {
namespace verify {

Eigen::VectorXd ecg_single_particle_density(
    const std::vector<BasisParams>& basis,
    const PermutationSet& perms,
    const Eigen::VectorXd& x_points,
    int /*particle_index*/) {
    // For the ECG symmetric N-particle state the code uses the convention
    //   |psi_R> = sum_j u_j sum_p |phi_j^p>    (ket fully permutation-summed)
    //   |psi_L> = sum_i u_i |phi_i>            (bra not permuted)
    //   S_code = <psi_L | psi_R> = sum_{ij,p} u_i^* u_j M_G(i,j,p)
    // The true bosonic overlap is <psi_R|psi_R> = N! * S_code, but the code's
    // self-consistent convention (overlap/energy/gradient routines) all use
    // S_code, so we match it here for the density normalization.
    //
    // The single-particle density <psi_R|delta(x-x_a)|psi_R>/<psi_R|psi_R>
    // works out (by pushing one permutation through the delta) to:
    //   n(x) = (1/N) * sum_{a'=0..N-1} n_half(x; a')
    // where n_half(x; a) = (1/S_code) * sum_{ij,p} u_i^* u_j sign(p) *
    //                     <phi_i | delta(x - x_a) | phi_j^p>
    // is the "half-symmetrized" contribution with bra fixed and ket summed.
    // That's what the original implementation computed for a single a — but
    // for N>=2 with basis functions that are not themselves symmetric under
    // particle swap (A_11 != A_22), n_half(x; 0) != n_half(x; 1), and using
    // only a=0 gives the wrong result. Averaging over all a restores full
    // bosonic symmetry.
    //
    // The `particle_index` argument is therefore ignored (kept for ABI). The
    // output still integrates to 1 per particle.
    int K  = static_cast<int>(basis.size());
    int nx = static_cast<int>(x_points.size());
    int N  = basis[0].N();

    Cd S = overlap(basis);
    Eigen::VectorXd n = Eigen::VectorXd::Zero(nx);

    for (int a = 0; a < N; a++) {
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                Cd u_prod = std::conj(basis[i].u) * basis[j].u;
                for (int p = 0; p < perms.SN; p++) {
                    double sign = static_cast<double>(perms.signs[p]);
                    PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);

                    // Marginal Gaussian for particle a (Schur complement):
                    //   K_tilde = 1 / (K^{-1})_{aa}
                    //   mu_a    = (K^{-1} b / 2)_a
                    //   density = sign * u_i^* u_j * M_G * sqrt(K_tilde/pi)
                    //             * exp(-K_tilde * (x - mu_a)^2)
                    Cd K_aa_inv = c.K_inv(a, a);
                    Cd K_tilde  = 1.0 / K_aa_inv;
                    Cd mu_a     = c.mu(a);
                    Cd sqrt_Ktilde_over_pi = std::sqrt(K_tilde / M_PI);
                    Cd prefac = sign * u_prod * c.M_G * sqrt_Ktilde_over_pi;

                    for (int xi = 0; xi < nx; xi++) {
                        Cd dx = x_points(xi) - mu_a;
                        n(xi) += (prefac * std::exp(-K_tilde * dx * dx)).real();
                    }
                }
            }
        }
    }

    double Sr = S.real();
    if (std::abs(Sr) > 0) n /= Sr;
    if (N > 1) n /= static_cast<double>(N);       // average over particles
    return n;
}

Eigen::VectorXcd ecg_wavefunction_1p(
    const std::vector<BasisParams>& basis,
    const Eigen::VectorXd& x_points) {
    if (basis.empty() || basis[0].N() != 1) {
        throw std::runtime_error("ecg_wavefunction_1p requires N=1");
    }
    int K  = static_cast<int>(basis.size());
    int nx = static_cast<int>(x_points.size());
    Eigen::VectorXcd psi = Eigen::VectorXcd::Zero(nx);

    // phi_i(x) = u_i * exp(-(A_i + B_i) x^2 + 2 R_i B_i x - R_i B_i R_i)
    //   (N=1 so A_i, B_i, R_i are scalars.)
    for (int i = 0; i < K; i++) {
        Cd u_i = basis[i].u;
        Cd A_i = basis[i].A(0, 0);
        Cd B_i = basis[i].B(0, 0);
        Cd R_i = basis[i].R(0);
        Cd AB  = A_i + B_i;
        Cd lin = 2.0 * R_i * B_i;
        Cd cst = -R_i * B_i * R_i;
        for (int j = 0; j < nx; j++) {
            double x = x_points(j);
            Cd exponent = -AB * x * x + lin * x + cst;
            psi(j) += u_i * std::exp(exponent);
        }
    }
    return psi;
}

Eigen::VectorXd ecg_momentum_density_1p(
    const std::vector<BasisParams>& basis,
    const Eigen::VectorXd& k_points) {
    if (basis.empty() || basis[0].N() != 1) {
        throw std::runtime_error("ecg_momentum_density_1p requires N=1");
    }
    int K  = static_cast<int>(basis.size());
    int nk = static_cast<int>(k_points.size());

    double Sr = overlap(basis).real();
    Eigen::VectorXd nk_out(nk);

    Cd inv_sqrt_2pi = Cd(1.0 / std::sqrt(2.0 * M_PI), 0.0);

    for (int ki = 0; ki < nk; ki++) {
        double k = k_points(ki);
        Cd psi_tilde(0.0, 0.0);
        for (int i = 0; i < K; i++) {
            Cd u_i = basis[i].u;
            Cd A_i = basis[i].A(0, 0);
            Cd B_i = basis[i].B(0, 0);
            Cd R_i = basis[i].R(0);
            Cd alpha = A_i + B_i;                // decay coefficient
            // integrand: exp(-alpha x^2 + (2 R B - i k) x - R B R)
            //   -> sqrt(pi/alpha) * exp(-R B R + (2 R B - i k)^2 / (4 alpha))
            Cd beta  = 2.0 * R_i * B_i - Cd(0.0, 1.0) * k;
            Cd gamma = -R_i * B_i * R_i;
            Cd exponent = gamma + (beta * beta) / (4.0 * alpha);
            Cd sqrt_pi_over_alpha = std::sqrt(Cd(M_PI, 0.0) / alpha);
            psi_tilde += u_i * sqrt_pi_over_alpha * std::exp(exponent);
        }
        psi_tilde *= inv_sqrt_2pi;
        nk_out(ki) = std::norm(psi_tilde);
    }

    if (Sr > 0) nk_out /= Sr;
    return nk_out;
}

Eigen::VectorXd momentum_density_from_density(
    const Eigen::VectorXd& x_points,
    const Eigen::VectorXd& n_x,
    const Eigen::VectorXd& k_points) {
    int nx = static_cast<int>(x_points.size());
    int nk = static_cast<int>(k_points.size());
    if (n_x.size() != nx) {
        throw std::runtime_error("momentum_density_from_density: size mismatch");
    }

    // phi(x) = sqrt(max(n(x), 0)) — single-particle orbital assuming the
    // ground state is real and non-negative (non-interacting bosonic GS).
    Eigen::VectorXd phi(nx);
    for (int j = 0; j < nx; j++) {
        phi(j) = std::sqrt(std::max(n_x(j), 0.0));
    }

    // Direct DFT on the x-grid (trapezoidal rule, uniform dx assumed).
    double dx = x_points(1) - x_points(0);
    double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    Eigen::VectorXd nk_out(nk);
    for (int ki = 0; ki < nk; ki++) {
        double k = k_points(ki);
        Cd phi_tilde(0.0, 0.0);
        for (int j = 0; j < nx; j++) {
            Cd phase = std::exp(Cd(0.0, -k * x_points(j)));
            phi_tilde += phi(j) * phase;
        }
        phi_tilde *= dx * inv_sqrt_2pi;
        nk_out(ki) = std::norm(phi_tilde);
    }
    return nk_out;
}

} // namespace verify
} // namespace ecg1d
