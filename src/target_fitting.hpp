#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include "tdvp_solver.hpp"
#include <vector>

namespace ecg1d {

// Target state = Jacobi-Anger-truncated kicked state, stored as a list of
// fixed ECG components. Each target_basis[m].u stores the complex coefficient
// c_m; everything else (A, B, R, name) defines |chi_m>.
using TargetBasis = std::vector<BasisParams>;

// Redistribute (A, B) within each diagonal axis so that B_aa has real part
// >= b_min, keeping A+B invariant (so the wave function is unchanged).
// The ECG kick shift R -> R + i*n*k_L/B_aa is only well-defined when B is
// well away from zero; a typical harmonic-oscillator ground state comes out
// with B=0 and needs this rebalance before target construction.
std::vector<BasisParams> rebalance_AB_for_kick(const std::vector<BasisParams>& basis,
                                               double b_min = 0.5);

// Build |psi_tar> = exp(-i*kappa*cos(2*k_L*x)) |psi_pre_kick> truncated to |n|<=n_max.
// For each pre_kick_basis[j] and each n in [-n_max, n_max] with |J_n(kappa)| above threshold,
// emit one target component:
//   c_m  = u_j^(0) * (-i)^n * J_n(kappa)
//   A_m  = A_j, B_m = B_j
//   R_m  = R_j + i * n * k_L / B_j (diagonal, per particle axis in N>=1)
//   name = target_name_offset + j*(2*n_max+1) + (n + n_max)
// N=1 is the primary use; for N>=2 the momentum shift is broadcast to every particle axis
// (consistent with a global cos(2*k_L*x_a) kick factored over particles).
TargetBasis build_kicked_target(const std::vector<BasisParams>& pre_kick_basis,
                                double kappa, double k_L, int n_max,
                                int target_name_offset = 10000);

// A(z) = <psi_tar|psi(z)> = sum_{m,i} conj(target[m].u) * trial[i].u *
//                           sum_perm sign(perm) * M_G(target[m], trial[i], perm)
Cd cross_overlap(const TargetBasis& target,
                 const std::vector<BasisParams>& trial,
                 const PermutationSet& perms);

// N_tar = <psi_tar|psi_tar> (convenience wrapper around overlap()).
Cd target_norm(const TargetBasis& target, const PermutationSet& perms);

// b_alpha(z) = <d_alpha psi(z) | psi_tar> where d_alpha is the BRA-side
// anti-holomorphic Wirtinger derivative d/d z_alpha^*. Returns a vector
// aligned with alpha_z_list.
//
// - alpha.a1 == 1 (u_i):          b = sum_m target[m].u * sum_perm sign * M_G(trial[i], target[m], perm)
// - alpha.a1 in {2,3,4} (B,R,A):  b = conj(trial[i].u) * sum_m target[m].u *
//                                       partial_z_addsn(a1, Real=false,
//                                                       pi=trial[i], pj=target[m],
//                                                       alpha_2=trial[i].name, alpha_3, alpha_4, perms)
// Real=false is the bra side (see derivatives.cpp resolve_case_first).
VectorXcd assemble_b_target(const std::vector<AlphaIndex>& alpha_z_list,
                            const std::vector<BasisParams>& trial,
                            const TargetBasis& target,
                            const PermutationSet& perms);

// g_alpha = -(A/N_tar) * b_alpha. Also returns A and N_tar via output params
// for the caller's fidelity monitor.
VectorXcd assemble_grad_target_fit(const std::vector<AlphaIndex>& alpha_z_list,
                                   const std::vector<BasisParams>& trial,
                                   const TargetBasis& target,
                                   const PermutationSet& perms,
                                   Cd& out_A, Cd& out_N_tar);

// Analytic optimum u given fixed {A,B,R}: solve S * u = b_target where
// b_target[i] = sum_m target[m].u * sum_perm sign * M_G(trial[i], target[m], perm).
// The returned u is renormalized so that <psi_new|S|psi_new> = <psi_old|S|psi_old>.
VectorXcd solve_optimal_u(const std::vector<BasisParams>& trial,
                          const TargetBasis& target,
                          const PermutationSet& perms);

// Target-fitting TDVP configuration (orthogonal to SolverConfig).
struct TargetFitConfig {
    enum UMode {
        UPDATE_LINEAR,   // include u in alpha_z_list, solve C dz = -g for everything
        UPDATE_RESOLVE   // strip u from alpha_z_list, step {A,B,R}, then re-solve u analytically
    };
    UMode u_mode = UPDATE_RESOLVE;

    double fidelity_tol = 1e-10;   // stop when 1 - F < fidelity_tol
    int    max_steps    = 2000;
    double dtao_init    = 1e-3;
    double dtao_max     = 0.0;     // 0 => 100 * dtao_init
    double dtao_grow    = 1.5;
    double lambda_C     = 1e-8;    // Tikhonov on C
    double rcond        = 1e-4;    // SVD truncation threshold
    bool   verbose      = true;
    int    print_every  = 50;      // verbose: print F every N steps
};

struct TargetFitResult {
    std::vector<BasisParams> basis;
    double              F_final;
    int                 n_steps;
    std::vector<double> F_history;     // per-step fidelity |A|^2/(N_tar * S)
    std::vector<double> A_abs_history; // per-step |A(z)|
    int                 line_search_fails;
};

TargetFitResult target_fitting_evolution(const std::vector<AlphaIndex>& alpha_z_list,
                                         std::vector<BasisParams> trial_init,
                                         const TargetBasis& target,
                                         const TargetFitConfig& fit_cfg = {});

} // namespace ecg1d
