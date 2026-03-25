#include "kicked_evolution.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_gradient.hpp"
#include "derivatives.hpp"
#include "svm.hpp"
#include "pair_cache.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

namespace ecg1d {

RealtimeStepResult realtime_tdvp_step(
    const std::vector<AlphaIndex>& alpha_z_list,
    const std::vector<BasisParams>& basis,
    double dt,
    const HamiltonianTerms& terms,
    const SolverConfig& config,
    const PermutationSet* perms) {

    int updata_constant = 0;
    for (const auto& a : alpha_z_list) {
        if (a.a1 == 1) updata_constant++;
    }

    // Assemble C matrix and gradient
    MatrixXcd C_1 = assemble_C(alpha_z_list, basis);
    MatrixXcd C_bar = C_1.conjugate().transpose();
    VectorXcd g_bar = assemble_grad(alpha_z_list, false, basis, terms);

    int d_update = static_cast<int>(alpha_z_list.size()) - updata_constant;
    MatrixXcd C_bar_update = C_bar.bottomRightCorner(d_update, d_update);
    VectorXcd g_bar_update = g_bar.tail(d_update);

    // Tikhonov regularization
    if (config.lambda_C > 0.0) {
        C_bar_update += config.lambda_C * MatrixXcd::Identity(d_update, d_update);
    }

    // SVD solve
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

    // Real-time: rhs = -i * g (multiply gradient by -i for Schrodinger evolution)
    VectorXcd rhs = Cd(0, -1) * g_bar_update;
    VectorXcd Utrhs = svd.matrixU().adjoint() * rhs;
    VectorXcd sinv_utrhs(sv.size());
    for (int i = 0; i < sv.size(); i++) {
        sinv_utrhs(i) = (sv(i) > threshold) ? Utrhs(i) / sv(i) : Cd(0, 0);
    }
    VectorXcd dz = svd.matrixV() * sinv_utrhs;

    double norm_dz = dz.norm();

    // Apply update (no line search for real-time evolution)
    std::vector<BasisParams> basis_new = basis;
    update_basis_function(basis_new, dz, dt, alpha_z_list, updata_constant);

    // Re-solve u via eigenvalue if requested
    if (config.resolve_u && perms != nullptr) {
        auto [H_new, S_new] = build_HS(basis_new, *perms, terms);
        double E_eigen = lowest_energy(H_new, S_new);
        if (std::isfinite(E_eigen)) {
            set_u_from_eigenvector(basis_new, H_new, S_new);
            return {basis_new, Cd(E_eigen, 0.0), norm_dz, dt, diag};
        }
    }

    // Renormalize u to preserve norm after parameter update.
    // When A/B/R change, S changes, and u†Su drifts. This corrects it.
    {
        int Kb = static_cast<int>(basis_new.size());
        VectorXcd u_vec(Kb);
        for (int i = 0; i < Kb; i++) u_vec(i) = basis_new[i].u;

        // Compute norm before (with old basis) and after (with new basis)
        VectorXcd u_old(Kb);
        for (int i = 0; i < Kb; i++) u_old(i) = basis[i].u;

        // We need S_old and S_new. For efficiency, compute S_new from
        // the basis_new overlap. S_old norm is approximately preserved
        // from the previous step.
        if (perms != nullptr) {
            auto [H_new, S_new] = build_HS(basis_new, *perms, terms);
            Cd norm_new = (u_vec.adjoint() * S_new * u_vec)(0);

            auto [H_old, S_old] = build_HS(basis, *perms, terms);
            Cd norm_old = (u_old.adjoint() * S_old * u_old)(0);

            if (norm_new.real() > 1e-15 && norm_old.real() > 1e-15) {
                double scale = std::sqrt(norm_old.real() / norm_new.real());
                for (int i = 0; i < Kb; i++) {
                    basis_new[i].u *= scale;
                }
            }
        }
    }

    Cd E_new = compute_total_energy(basis_new, terms);
    return {basis_new, E_new, norm_dz, dt, diag};
}

void free_evolve(std::vector<BasisParams>& basis,
                 const std::vector<AlphaIndex>& alpha_z_list,
                 double dt_step, double T_duration,
                 const HamiltonianTerms& terms,
                 const SolverConfig& config,
                 const PermutationSet* perms) {

    int n_steps = static_cast<int>(std::ceil(T_duration / dt_step));
    double dt_actual = T_duration / n_steps;

    for (int s = 0; s < n_steps; s++) {
        auto result = realtime_tdvp_step(alpha_z_list, basis, dt_actual,
                                          terms, config, perms);
        basis = result.basis;
    }
}

void apply_kick(std::vector<BasisParams>& basis,
                const std::vector<AlphaIndex>& alpha_z_list,
                double dt_step, double T_pulse,
                const HamiltonianTerms& terms_with_kick,
                const SolverConfig& config,
                const PermutationSet* perms) {

    if (T_pulse <= 0.0) {
        // Ideal delta kick: single step with small dt
        auto result = realtime_tdvp_step(alpha_z_list, basis, dt_step,
                                          terms_with_kick, config, perms);
        basis = result.basis;
    } else {
        // Finite pulse: evolve through the pulse duration
        free_evolve(basis, alpha_z_list, dt_step, T_pulse,
                    terms_with_kick, config, perms);
    }
}

std::vector<Observable> kicked_evolution(
    std::vector<BasisParams>& basis,
    const std::vector<AlphaIndex>& alpha_z_list,
    const KickParams& kick_params,
    double dt_step,
    const HamiltonianTerms& terms_free,
    const HamiltonianTerms& terms_kick,
    const SolverConfig& config,
    const PermutationSet* perms) {

    std::vector<Observable> observables;
    double current_time = 0.0;

    // Record initial state
    {
        Cd E = compute_total_energy(basis, terms_free);
        Cd S = overlap(basis);
        Cd T = kinetic_energy_functional(basis);
        Observable obs;
        obs.time = current_time;
        obs.energy = E.real();
        obs.kinetic_energy = (T / S).real();
        obs.overlap_norm = S.real();
        observables.push_back(obs);

        std::cout << "Kick  0: t=" << std::setprecision(4) << current_time
                  << ", E=" << std::setprecision(8) << obs.energy
                  << ", T=" << std::setprecision(8) << obs.kinetic_energy
                  << ", S=" << std::setprecision(6) << obs.overlap_norm
                  << std::endl;
    }

    for (int n = 0; n < kick_params.n_kicks; n++) {
        // 1. Free evolution for one period
        double T_free = kick_params.T_period;
        if (kick_params.T_pulse > 0.0) {
            T_free -= kick_params.T_pulse;
        }
        if (T_free > 0.0) {
            free_evolve(basis, alpha_z_list, dt_step, T_free,
                        terms_free, config, perms);
        }
        current_time += T_free;

        // 2. Apply kick
        apply_kick(basis, alpha_z_list, dt_step, kick_params.T_pulse,
                   terms_kick, config, perms);
        current_time += kick_params.T_pulse;

        // 3. Record observables
        Cd E = compute_total_energy(basis, terms_free);
        Cd S = overlap(basis);
        Cd T = kinetic_energy_functional(basis);
        Observable obs;
        obs.time = current_time;
        obs.energy = E.real();
        obs.kinetic_energy = (T / S).real();
        obs.overlap_norm = S.real();
        observables.push_back(obs);

        std::cout << "Kick " << std::setw(2) << (n + 1)
                  << ": t=" << std::setprecision(4) << current_time
                  << ", E=" << std::setprecision(8) << obs.energy
                  << ", T=" << std::setprecision(8) << obs.kinetic_energy
                  << ", S=" << std::setprecision(6) << obs.overlap_norm
                  << std::endl;
    }

    return observables;
}

} // namespace ecg1d
