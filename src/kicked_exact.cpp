#include "kicked_exact.hpp"
#include "physical_constants.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

namespace ecg1d {

using Cd = std::complex<double>;

KickedExactResult kicked_exact_1particle(
    double omega, double m,
    double k_L, double kappa,
    double T_period,
    int n_kicks,
    int N_grid,
    double z_max)
{
    // Grid: z_j = -z_max + j * dz, j = 0, ..., N_grid-1
    double dz = 2.0 * z_max / (N_grid - 1);
    Eigen::VectorXd z(N_grid);
    for (int j = 0; j < N_grid; j++) {
        z(j) = -z_max + j * dz;
    }

    // Build H_free as tridiagonal matrix (finite difference + harmonic trap)
    // -1/(2m) d²ψ/dz² ≈ -1/(2m) (ψ_{j+1} - 2ψ_j + ψ_{j-1}) / dz²
    // V(z) = 1/2 m ω² z²
    double coeff = 1.0 / (2.0 * m * dz * dz);
    Eigen::MatrixXd H_free = Eigen::MatrixXd::Zero(N_grid, N_grid);
    for (int j = 0; j < N_grid; j++) {
        H_free(j, j) = 2.0 * coeff + 0.5 * m * omega * omega * z(j) * z(j);
        if (j > 0)         H_free(j, j-1) = -coeff;
        if (j < N_grid-1)  H_free(j, j+1) = -coeff;
    }

    // Diagonalize H_free: H = U Λ U^T (real symmetric)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_free);
    Eigen::VectorXd eigenvalues = es.eigenvalues();
    Eigen::MatrixXd U = es.eigenvectors();

    // Precompute exp(-i λ_k T) for free propagator
    Eigen::VectorXcd phase(N_grid);
    for (int k = 0; k < N_grid; k++) {
        phase(k) = std::exp(Cd(0, -eigenvalues(k) * T_period));
    }

    // Precompute kick phase: exp(-i κ cos(2 k_L z_j))
    Eigen::VectorXcd kick_phase(N_grid);
    for (int j = 0; j < N_grid; j++) {
        kick_phase(j) = std::exp(Cd(0, -kappa * std::cos(2.0 * k_L * z(j))));
    }

    // Build kinetic energy matrix for observable (same finite difference)
    // Only need diagonal and off-diagonal for <T>
    // But easier: compute <T> = <ψ|T|ψ> using T matrix directly
    // T_jj = coeff, T_{j,j±1} = -coeff/2 ... actually T = -1/(2m) d²/dz²
    // T_jj = 2*coeff/2 = coeff (wait, let me be careful)
    // The kinetic part of H_free is: T_{jj} = 2*coeff, T_{j,j±1} = -coeff
    // But H_free = T + V, so T = H_free - V
    Eigen::MatrixXd T_mat = Eigen::MatrixXd::Zero(N_grid, N_grid);
    for (int j = 0; j < N_grid; j++) {
        T_mat(j, j) = 2.0 * coeff;
        if (j > 0)         T_mat(j, j-1) = -coeff;
        if (j < N_grid-1)  T_mat(j, j+1) = -coeff;
    }

    // Initial state: ground state of H_free (column 0 of U)
    Eigen::VectorXcd psi = U.col(0).cast<std::complex<double>>();
    // Normalize (should already be, but just in case)
    psi /= std::sqrt((psi.adjoint() * psi)(0).real() * dz);
    // Actually eigenvectors are normalized to unit L2 in discrete sense,
    // but we need continuous normalization: ∫|ψ|²dz = 1 → Σ|ψ_j|²·dz = 1
    double norm0 = (psi.adjoint() * psi)(0).real() * dz;
    psi /= std::sqrt(norm0);

    KickedExactResult result;

    // Record initial energy
    {
        Cd E_val = (psi.adjoint() * (H_free.cast<Cd>() * psi))(0) * dz;
        Cd T_val = (psi.adjoint() * (T_mat.cast<Cd>() * psi))(0) * dz;
        double norm = (psi.adjoint() * psi)(0).real() * dz;
        result.kick_energies.push_back(E_val.real());
        result.kick_kinetic.push_back(T_val.real());
        result.kick_norm.push_back(norm);
    }

    // Main evolution loop
    for (int n = 0; n < n_kicks; n++) {
        // 1. Free evolution: ψ → exp(-iH_free T) ψ
        //    = U exp(-iΛT) U† ψ
        Eigen::VectorXcd coeff_vec = U.transpose().cast<Cd>() * psi;  // U† ψ
        for (int k = 0; k < N_grid; k++) {
            coeff_vec(k) *= phase(k);  // exp(-iλ_k T)
        }
        psi = U.cast<Cd>() * coeff_vec;  // U * ...

        // 2. Apply kick: ψ(z_j) → exp(-iκ cos(2k_L z_j)) ψ(z_j)
        for (int j = 0; j < N_grid; j++) {
            psi(j) *= kick_phase(j);
        }

        // 3. Record observables
        Cd E_val = (psi.adjoint() * (H_free.cast<Cd>() * psi))(0) * dz;
        Cd T_val = (psi.adjoint() * (T_mat.cast<Cd>() * psi))(0) * dz;
        double norm = (psi.adjoint() * psi)(0).real() * dz;

        result.kick_energies.push_back(E_val.real());
        result.kick_kinetic.push_back(T_val.real());
        result.kick_norm.push_back(norm);
    }

    return result;
}

void run_kicked_exact_test(int n_kicks) {
    std::cout << "\n=== Exact Kicked Rotor (N=2, no interaction, finite difference) ===" << std::endl;

    double omega_val = omega;
    double mass_val = mass;
    double k_L_val = k_L;
    double kappa_val = kappa;
    double T_period = 1.0;

    std::cout << "Parameters: omega=" << omega_val << ", mass=" << mass_val
              << ", k_L=" << k_L_val << ", kappa=" << kappa_val
              << ", T=" << T_period << ", n_kicks=" << n_kicks << std::endl;

    int N_grid = 512;
    double z_max = 15.0;
    std::cout << "Grid: N=" << N_grid << ", z_max=" << z_max << std::endl;

    auto res = kicked_exact_1particle(omega_val, mass_val, k_L_val, kappa_val,
                                       T_period, n_kicks, N_grid, z_max);

    // N=2 non-interacting: E_2 = 2 * E_1
    std::cout << "\n" << std::setw(6) << "kick"
              << std::setw(16) << "E_1particle"
              << std::setw(16) << "E_2particle"
              << std::setw(16) << "T_1particle"
              << std::setw(12) << "norm" << std::endl;

    for (int n = 0; n <= n_kicks; n++) {
        std::cout << std::setw(6) << n
                  << std::setw(16) << std::setprecision(8) << res.kick_energies[n]
                  << std::setw(16) << std::setprecision(8) << 2.0 * res.kick_energies[n]
                  << std::setw(16) << std::setprecision(8) << res.kick_kinetic[n]
                  << std::setw(12) << std::setprecision(8) << res.kick_norm[n]
                  << std::endl;
    }

    // Summary
    std::cout << "\nSummary:" << std::endl;
    std::cout << "  Initial E (1-particle): " << std::setprecision(10) << res.kick_energies[0] << std::endl;
    std::cout << "  Final E (1-particle):   " << std::setprecision(10) << res.kick_energies[n_kicks] << std::endl;
    std::cout << "  Energy absorbed:        " << std::setprecision(10)
              << res.kick_energies[n_kicks] - res.kick_energies[0] << std::endl;
    std::cout << "  Norm drift:             " << std::scientific
              << std::abs(res.kick_norm[n_kicks] - 1.0) << std::defaultfloat << std::endl;
    std::cout << "  Initial E (2-particle): " << std::setprecision(10) << 2.0 * res.kick_energies[0] << std::endl;
    std::cout << "  Final E (2-particle):   " << std::setprecision(10) << 2.0 * res.kick_energies[n_kicks] << std::endl;
}

} // namespace ecg1d
