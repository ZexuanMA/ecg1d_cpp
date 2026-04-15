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

// Helper: compute dz = C^{-1}(-ig) at a given basis configuration
static VectorXcd compute_realtime_dz(
    const std::vector<AlphaIndex>& alpha_z_list,
    const std::vector<BasisParams>& basis,
    const HamiltonianTerms& terms,
    const SolverConfig& config,
    int updata_constant,
    TdvpDiagnostics* diag_out = nullptr) {

    int d_update = static_cast<int>(alpha_z_list.size()) - updata_constant;

    MatrixXcd C_1 = assemble_C(alpha_z_list, basis);
    MatrixXcd C_bar = C_1.conjugate().transpose();
    VectorXcd g_bar = assemble_grad(alpha_z_list, false, basis, terms);

    MatrixXcd C_bar_update = C_bar.bottomRightCorner(d_update, d_update);
    VectorXcd g_bar_update = g_bar.tail(d_update);

    if (config.lambda_C > 0.0) {
        C_bar_update += config.lambda_C * MatrixXcd::Identity(d_update, d_update);
    }

    Eigen::JacobiSVD<MatrixXcd> svd(C_bar_update, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double max_sv = svd.singularValues()(0);
    double threshold = config.rcond * max_sv;
    Eigen::VectorXd sv = svd.singularValues();

    if (diag_out) {
        diag_out->max_sv_C = sv(0);
        diag_out->min_sv_C = sv(sv.size() - 1);
        diag_out->cond_C = (diag_out->min_sv_C > 0) ? diag_out->max_sv_C / diag_out->min_sv_C : 1e30;
        diag_out->effective_rank = 0;
        for (int i = 0; i < sv.size(); i++) {
            if (sv(i) > threshold) diag_out->effective_rank++;
        }
    }

    VectorXcd rhs = Cd(0, -1) * g_bar_update;
    VectorXcd Utrhs = svd.matrixU().adjoint() * rhs;
    VectorXcd sinv_utrhs(sv.size());
    for (int i = 0; i < sv.size(); i++) {
        sinv_utrhs(i) = (sv(i) > threshold) ? Utrhs(i) / sv(i) : Cd(0, 0);
    }
    return svd.matrixV() * sinv_utrhs;
}

RealtimeStepResult realtime_tdvp_step(
    const std::vector<AlphaIndex>& alpha_z_list,
    const std::vector<BasisParams>& basis,
    double dt,
    const HamiltonianTerms& terms,
    const SolverConfig& config,
    const PermutationSet* perms) {

    // For real-time evolution, include u in the TDVP update (updata_constant=0).
    // This ensures u evolves consistently with z (A/B/R).
    // For imaginary time, u is excluded and re-solved as eigenstate.
    int updata_constant = 0;  // include u in the update

    // --- Classical RK4 method ---
    // k1 = f(t, y)
    TdvpDiagnostics diag;
    VectorXcd k1 = compute_realtime_dz(alpha_z_list, basis, terms, config,
                                        updata_constant, &diag);

    // k2 = f(t + dt/2, y + dt/2 * k1)
    std::vector<BasisParams> basis_s2 = basis;
    update_basis_function(basis_s2, k1, dt * 0.5, alpha_z_list, updata_constant);
    VectorXcd k2 = compute_realtime_dz(alpha_z_list, basis_s2, terms, config,
                                        updata_constant);

    // k3 = f(t + dt/2, y + dt/2 * k2)
    std::vector<BasisParams> basis_s3 = basis;
    update_basis_function(basis_s3, k2, dt * 0.5, alpha_z_list, updata_constant);
    VectorXcd k3 = compute_realtime_dz(alpha_z_list, basis_s3, terms, config,
                                        updata_constant);

    // k4 = f(t + dt, y + dt * k3)
    std::vector<BasisParams> basis_s4 = basis;
    update_basis_function(basis_s4, k3, dt, alpha_z_list, updata_constant);
    VectorXcd k4 = compute_realtime_dz(alpha_z_list, basis_s4, terms, config,
                                        updata_constant);

    // y_{n+1} = y_n + (dt/6)(k1 + 2*k2 + 2*k3 + k4)
    VectorXcd dz_rk4 = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    std::vector<BasisParams> basis_new = basis;
    update_basis_function(basis_new, dz_rk4, dt, alpha_z_list, updata_constant);

    double norm_dz = dz_rk4.norm();

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
