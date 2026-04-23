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
    int particle_index) {

    int K = static_cast<int>(basis.size());
    int n_x = static_cast<int>(x_points.size());
    int a = particle_index;

    Cd S = overlap(basis);
    Eigen::VectorXd n = Eigen::VectorXd::Zero(n_x);

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            Cd u_prod = std::conj(basis[i].u) * basis[j].u;
            for (int p = 0; p < perms.SN; p++) {
                double sign = static_cast<double>(perms.signs[p]);
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);

                // Marginal Gaussian for particle a:
                //   K_tilde = 1 / (K^{-1})_{aa}
                //   mu_a    = (K^{-1} b / 2)_a   [already stored in c.mu(a)]
                //   prefactor = sign * u_i^* u_j * M_G * sqrt(K_tilde / pi)
                //   density  = prefactor * exp(-K_tilde * (x - mu_a)^2)
                Cd K_aa_inv = c.K_inv(a, a);
                Cd K_tilde  = 1.0 / K_aa_inv;
                Cd mu_a     = c.mu(a);
                Cd sqrt_Ktilde_over_pi = std::sqrt(K_tilde / M_PI);
                Cd prefac = sign * u_prod * c.M_G * sqrt_Ktilde_over_pi;

                for (int xi = 0; xi < n_x; xi++) {
                    Cd dx = x_points(xi) - mu_a;
                    Cd contrib = prefac * std::exp(-K_tilde * dx * dx);
                    n(xi) += contrib.real();
                }
            }
        }
    }

    double Sr = S.real();
    if (std::abs(Sr) > 0) n /= Sr;
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

} // namespace verify
} // namespace ecg1d
