#include "tdvp_solver.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_gradient.hpp"
#include "derivatives.hpp"
#include "svm.hpp"
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
                          const HamiltonianTerms& terms,
                          const SolverConfig& config,
                          const PermutationSet* perms) {
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

    // A2: Tikhonov regularization
    if (config.lambda_C > 0.0) {
        C_bar_update += config.lambda_C * MatrixXcd::Identity(d_update, d_update);
    }

    // Solve C_bar_update * dz = -g_bar_update using SVD with rcond threshold
    Eigen::JacobiSVD<MatrixXcd> svd(C_bar_update, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double max_sv = svd.singularValues()(0);
    double threshold = config.rcond * max_sv;
    Eigen::VectorXd sv = svd.singularValues();

    // Build diagnostics
    TdvpDiagnostics diag;
    diag.max_sv_C = sv(0);
    diag.min_sv_C = sv(sv.size() - 1);
    diag.cond_C = (diag.min_sv_C > 0) ? diag.max_sv_C / diag.min_sv_C : 1e30;
    diag.effective_rank = 0;
    for (int i = 0; i < sv.size(); i++) {
        if (sv(i) > threshold) diag.effective_rank++;
    }

    // Compute smallest singular value of overlap matrix S
    {
        int K = static_cast<int>(basis.size());
        int N = basis[0].N();
        MatrixXcd S_mat = MatrixXcd::Zero(K, K);
        for (int i = 0; i < K; i++) {
            for (int j = i; j < K; j++) {
                // Use overlap functional approach
                Cd s_ij(0, 0);
                PermutationSet local_perms = PermutationSet::generate(N);
                for (int p = 0; p < local_perms.SN; p++) {
                    double sign = static_cast<double>(local_perms.signs[p]);
                    PairCache c = PairCache::build(basis[i], basis[j], local_perms.matrices[p]);
                    s_ij += sign * c.M_G;
                }
                S_mat(i, j) = s_ij;
                S_mat(j, i) = std::conj(s_ij);
            }
        }
        Eigen::MatrixXd Ss = (0.5 * (S_mat + S_mat.adjoint())).real();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ss);
        if (es.info() == Eigen::Success) {
            diag.min_sv_S = es.eigenvalues()(0);
        }
    }

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

        // Always re-solve u if requested and perms provided
        if (config.resolve_u && perms != nullptr) {
            auto [H_trial, S_trial] = build_HS(basis_trial, *perms, terms);
            set_u_from_eigenvector(basis_trial, H_trial, S_trial);
            double E_eigen = lowest_energy(H_trial, S_trial);
            if (E_eigen < E_pre.real()) {
                return {basis_trial, Cd(E_eigen, 0.0), E_pre, norm_dz, dtao, diag};
            }
        } else {
            Cd E_new = compute_total_energy(basis_trial, terms);
            if (E_new.real() < E_pre.real()) {
                return {basis_trial, E_new, E_pre, norm_dz, dtao, diag};
            }
        }

        dtao *= 0.5;
        if (dtao < 1e-20) break;
    }

    // Exploratory tiny step before giving up
    double tiny_dtao = 1e-6;
    {
        std::vector<BasisParams> basis_trial = basis_return;
        update_basis_function(basis_trial, dz, tiny_dtao, alpha_z_list, updata_constant);

        if (config.resolve_u && perms != nullptr) {
            auto [H_trial, S_trial] = build_HS(basis_trial, *perms, terms);
            set_u_from_eigenvector(basis_trial, H_trial, S_trial);
            double E_eigen = lowest_energy(H_trial, S_trial);
            if (E_eigen < E_pre.real()) {
                return {basis_trial, Cd(E_eigen, 0.0), E_pre, norm_dz, tiny_dtao, diag};
            }
        } else {
            Cd E_new = compute_total_energy(basis_trial, terms);
            if (E_new.real() < E_pre.real()) {
                return {basis_trial, E_new, E_pre, norm_dz, tiny_dtao, diag};
            }
        }
    }

    return {basis_return, E_pre, E_pre, norm_dz, 0.0, diag};
}

void evolution(const std::vector<AlphaIndex>& alpha_z_list,
               std::vector<BasisParams>& basis,
               double dtao, int max_steps, double tol,
               const HamiltonianTerms& terms,
               const SolverConfig& config,
               const PermutationSet* perms) {
    double current_dtao = dtao;
    double dtao_max = (config.dtao_max > 0.0) ? config.dtao_max : 100.0 * dtao;

    // Use a mutable copy of config for adaptive lambda
    SolverConfig active_config = config;

    // Tracking for stagnation detection and diagnostics
    double best_energy = std::numeric_limits<double>::infinity();
    int steps_since_improvement = 0;
    int total_line_search_failures = 0;
    int stagnation_events = 0;
    double min_cond_seen = std::numeric_limits<double>::infinity();
    int energy_converged_count = 0;
    double initial_energy = 0.0;
    double prev_energy = std::numeric_limits<double>::infinity();

    for (int step = 0; step < max_steps; step++) {
        auto result = tdvp_step(alpha_z_list, basis, current_dtao, terms, active_config, perms);

        double E_new = result.E_new.real();
        double dE = result.E_new.real() - result.E_pre.real();
        if (step == 0) initial_energy = result.E_pre.real();

        // Track diagnostics
        if (result.diag.cond_C < min_cond_seen) min_cond_seen = result.diag.cond_C;

        // A4: Enhanced diagnostic output
        std::cout << "Step " << step
                  << ": E=" << std::setprecision(10) << E_new
                  << ", dE=" << std::scientific << std::setprecision(2) << dE
                  << ", |dz|=" << result.norm_dz
                  << ", dtao=" << result.used_dtao
                  << ", cond(C)=" << std::setprecision(2) << result.diag.cond_C
                  << ", rank=" << result.diag.effective_rank
                  << ", lam=" << std::scientific << active_config.lambda_C
                  << std::defaultfloat << std::endl;

        if (result.used_dtao == 0.0) {
            total_line_search_failures++;

            // Recovery: try with increased regularization before stopping
            if (active_config.adaptive_lambda &&
                active_config.lambda_C < active_config.lambda_max) {
                active_config.lambda_C = std::min(active_config.lambda_C * 10.0,
                                                   active_config.lambda_max);
                std::cout << "  Line search failed, increasing lambda_C to "
                          << std::scientific << active_config.lambda_C
                          << std::defaultfloat << std::endl;
                continue;
            }

            // Try a larger step size as last resort
            if (current_dtao < dtao_max * 0.01) {
                current_dtao = dtao;  // Reset to initial step size
                std::cout << "  Line search failed, resetting dtao to " << current_dtao << std::endl;
                continue;
            }

            std::cout << "Line search failed. Stopping." << std::endl;
            break;
        }

        basis = result.basis;

        // Track best energy for stagnation detection
        if (E_new < best_energy - config.energy_tol) {
            best_energy = E_new;
            steps_since_improvement = 0;
            // Restore lambda if we had increased it
            if (active_config.lambda_C > config.lambda_C) {
                active_config.lambda_C = config.lambda_C;
            }
        } else {
            steps_since_improvement++;
        }

        // Energy-based convergence: several consecutive steps with tiny |dE|
        if (std::abs(dE) < config.energy_tol) {
            energy_converged_count++;
            if (energy_converged_count >= 5) {
                std::cout << "Converged: |dE| < " << std::scientific << config.energy_tol
                          << " for 5 consecutive steps." << std::defaultfloat << std::endl;
                break;
            }
        } else {
            energy_converged_count = 0;
        }

        // Stagnation recovery
        if (config.stagnation_window > 0 &&
            steps_since_improvement >= config.stagnation_window) {
            stagnation_events++;
            std::cout << "  Stagnation detected (no improvement for "
                      << steps_since_improvement << " steps)." << std::endl;

            // Strategy 1: increase regularization temporarily
            if (active_config.adaptive_lambda &&
                active_config.lambda_C < active_config.lambda_max) {
                active_config.lambda_C = std::min(active_config.lambda_C * 10.0,
                                                   active_config.lambda_max);
                std::cout << "  Increasing lambda_C to " << std::scientific
                          << active_config.lambda_C << std::defaultfloat << std::endl;
            }

            // Strategy 2: try a larger step size
            current_dtao = std::min(current_dtao * 5.0, dtao_max);
            std::cout << "  Boosting dtao to " << current_dtao << std::endl;
            steps_since_improvement = 0;

            // Give up after 3 stagnation events
            if (stagnation_events >= 3) {
                std::cout << "Multiple stagnation events. Stopping." << std::endl;
                break;
            }
        }

        // Step size growth
        if (config.dtao_grow > 1.0) {
            current_dtao = std::min(result.used_dtao * config.dtao_grow, dtao_max);
        } else {
            current_dtao = result.used_dtao;
        }

        if (result.norm_dz < tol) {
            std::cout << "Converged: |dz| " << result.norm_dz
                      << " < " << tol << " after " << step << " steps." << std::endl;
            break;
        }

        prev_energy = E_new;
    }

    // Diagnostic summary
    Cd E_final = compute_total_energy(basis, terms);
    std::cout << "\n--- Evolution summary ---" << std::endl;
    std::cout << "  Initial E:         " << std::setprecision(12) << initial_energy << std::endl;
    std::cout << "  Final E:           " << std::setprecision(12) << E_final.real() << std::endl;
    std::cout << "  Total dE:          " << std::scientific << (E_final.real() - initial_energy) << std::endl;
    std::cout << "  Line search fails: " << total_line_search_failures << std::endl;
    std::cout << "  Stagnation events: " << stagnation_events << std::endl;
    std::cout << "  Min cond(C):       " << std::scientific << min_cond_seen << std::defaultfloat << std::endl;
}

} // namespace ecg1d
