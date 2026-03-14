#include "tdvp_solver.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_gradient.hpp"
#include "derivatives.hpp"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

namespace ecg1d {

MatrixXcd assemble_C(const std::vector<AlphaIndex>& alpha_z_list,
                     const std::vector<BasisParams>& basis) {
    int d = static_cast<int>(alpha_z_list.size());
    MatrixXcd C_mat = MatrixXcd::Zero(d, d);
    for (int a = 0; a < d; a++) {
        const auto& alpha = alpha_z_list[a];
        for (int b = 0; b < d; b++) {
            const auto& beta = alpha_z_list[b];
            C_mat(a, b) = calculate_C(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                       beta.a1, beta.a2, beta.a3, beta.a4,
                                       basis);
        }
    }
    return C_mat;
}

Cd grad_H_for_alpha(const AlphaIndex& alpha, bool Real,
                    const std::vector<BasisParams>& basis,
                    const HamiltonianTerms& terms) {
    int N = basis[0].N();
    Cd result(0, 0);

    if (terms.kinetic)
        result += calculate_Hamiltonian_kinetic_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                        Real, basis);
    if (terms.harmonic)
        result += calculate_Hamiltonian_harmonic_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                         Real, basis);
    if (terms.delta && N >= 2)
        result += calculate_Hamiltonian_delta_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                      Real, basis);
    if (terms.gaussian && N >= 2)
        result += calculate_Hamiltonian_gaussian_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                         Real, basis);
    if (terms.kicking)
        result += calculate_Hamiltonian_kicking_partial(alpha.a1, alpha.a2, alpha.a3, alpha.a4,
                                                         Real, basis);

    return result;
}

VectorXcd assemble_grad(const std::vector<AlphaIndex>& alpha_z_list,
                        bool Real,
                        const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms) {
    int d = static_cast<int>(alpha_z_list.size());
    VectorXcd g = VectorXcd::Zero(d);
    for (int a = 0; a < d; a++) {
        g(a) = grad_H_for_alpha(alpha_z_list[a], Real, basis, terms);
    }
    return g;
}

void update_basis_function(std::vector<BasisParams>& basis,
                           const VectorXcd& dz, double dtao,
                           const std::vector<AlphaIndex>& alpha_z_list,
                           int updata_constant) {
    int length = static_cast<int>(alpha_z_list.size());
    for (int i = 0; i < length; i++) {
        const auto& idx = alpha_z_list[i];
        int dz_idx = i - updata_constant;
        if (idx.a1 == 1 || dz_idx < 0 || dz_idx >= dz.size()) continue;

        if (idx.a1 == 2) {
            basis[idx.a2].B(idx.a3, idx.a3) += dz(dz_idx) * dtao;
        } else if (idx.a1 == 3) {
            basis[idx.a2].R(idx.a3) += dz(dz_idx) * dtao;
        } else if (idx.a1 == 4) {
            Cd val = basis[idx.a2].A(idx.a3, idx.a4) + dz(dz_idx) * dtao;
            basis[idx.a2].A(idx.a3, idx.a4) = val;
            if (idx.a3 != idx.a4) {
                basis[idx.a2].A(idx.a4, idx.a3) = val;
            }
        }
    }
}

Cd compute_total_energy(const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms) {
    int N = basis[0].N();
    Cd S = overlap(basis);
    Cd E_num(0, 0);

    if (terms.kinetic)  E_num += kinetic_energy_functional(basis);
    if (terms.harmonic) E_num += Harmonic_functional(basis);
    if (terms.delta && N >= 2)    E_num += Delta_contact_functional(basis);
    if (terms.gaussian && N >= 2) E_num += Gaussian_interaction_functional(basis);
    if (terms.kicking)  E_num += kicking_term_functional(basis);

    return E_num / S;
}

TdvpStepResult tdvp_step(const std::vector<AlphaIndex>& alpha_z_list,
                          const std::vector<BasisParams>& basis,
                          double dtao,
                          const HamiltonianTerms& terms) {
    int updata_constant = 0;
    for (const auto& a : alpha_z_list) {
        if (a.a1 == 1) updata_constant++;
    }

    MatrixXcd C_1 = assemble_C(alpha_z_list, basis);
    MatrixXcd C_bar = C_1.conjugate().transpose();

    VectorXcd g_bar = assemble_grad(alpha_z_list, false, basis, terms);

    int d_update = static_cast<int>(alpha_z_list.size()) - updata_constant;
    MatrixXcd C_bar_update = C_bar.bottomRightCorner(d_update, d_update);
    VectorXcd g_bar_update = g_bar.tail(d_update);

    // Solve C_bar_update * dz = -g_bar_update using SVD with rcond threshold
    Eigen::JacobiSVD<MatrixXcd> svd(C_bar_update, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double rcond = 1e-4;
    double max_sv = svd.singularValues()(0);
    double threshold = rcond * max_sv;
    Eigen::VectorXd sv = svd.singularValues();

    VectorXcd rhs = -g_bar_update;
    VectorXcd Utrhs = svd.matrixU().adjoint() * rhs;
    VectorXcd sinv_utrhs(sv.size());
    for (int i = 0; i < sv.size(); i++) {
        sinv_utrhs(i) = (sv(i) > threshold) ? Utrhs(i) / sv(i) : Cd(0, 0);
    }
    VectorXcd dz = svd.matrixV() * sinv_utrhs;

    double norm_dz = dz.norm();
    Cd E_pre = compute_total_energy(basis, terms);

    std::vector<BasisParams> basis_return = basis;

    for (int test = 0; test < 100; test++) {
        std::vector<BasisParams> basis_trial = basis_return;
        update_basis_function(basis_trial, dz, dtao, alpha_z_list, updata_constant);
        Cd E_new = compute_total_energy(basis_trial, terms);

        if (E_new.real() < E_pre.real()) {
            return {basis_trial, E_new, E_pre, norm_dz, dtao};
        }
        dtao *= 0.5;
        if (dtao < 1e-20) {
            return {basis_return, E_pre, E_pre, norm_dz, 0.0};
        }
    }

    std::vector<BasisParams> basis_trial = basis_return;
    update_basis_function(basis_trial, dz, dtao, alpha_z_list, updata_constant);
    Cd E_new = compute_total_energy(basis_trial, terms);
    return {basis_trial, E_new, E_pre, norm_dz, dtao};
}

void evolution(const std::vector<AlphaIndex>& alpha_z_list,
               std::vector<BasisParams>& basis,
               double dtao, int max_steps, double tol,
               const HamiltonianTerms& terms) {
    double current_dtao = dtao;

    for (int step = 0; step < max_steps; step++) {
        auto result = tdvp_step(alpha_z_list, basis, current_dtao, terms);
        basis = result.basis;

        std::cout << "Step " << step
                  << ": E=" << std::setprecision(8) << result.E_new.real()
                  << ", dE=" << std::scientific << std::setprecision(2)
                  << (result.E_new.real() - result.E_pre.real())
                  << ", |dz|=" << result.norm_dz
                  << ", dtao=" << result.used_dtao
                  << std::defaultfloat << std::endl;

        if (result.used_dtao == 0.0) {
            std::cout << "Line search failed (step size too small). Stopping." << std::endl;
            break;
        }
        current_dtao = result.used_dtao;

        if (result.norm_dz < tol) {
            std::cout << "Converged: Gradient flow magnitude " << result.norm_dz
                      << " < " << tol << " after " << step << " steps." << std::endl;
            break;
        }
    }
}

} // namespace ecg1d
