#include "target_fitting.hpp"
#include "hamiltonian.hpp"
#include "pair_cache.hpp"
#include "derivatives.hpp"
#include "bessel.hpp"
#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace ecg1d {

// -----------------------------------------------------------------------------
// (A, B) redistribution to give B enough headroom for kick encoding
// -----------------------------------------------------------------------------

std::vector<BasisParams> rebalance_AB_for_kick(const std::vector<BasisParams>& basis,
                                               double b_min) {
    std::vector<BasisParams> out = basis;
    if (out.empty()) return out;
    const int N = out[0].N();
    for (auto& bp : out) {
        for (int a = 0; a < N; ++a) {
            double b_re = bp.B(a, a).real();
            if (b_re >= b_min) continue;
            // shift so Re(B_aa) = b_min, keep A+B invariant on diagonal
            double shift = b_min - b_re;
            bp.B(a, a) += Cd(shift, 0.0);
            bp.A(a, a) -= Cd(shift, 0.0);
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
// Target state construction (Jacobi-Anger expansion)
// -----------------------------------------------------------------------------

TargetBasis build_kicked_target(const std::vector<BasisParams>& pre_kick_basis,
                                double kappa, double k_L, int n_max,
                                int target_name_offset) {
    if (pre_kick_basis.empty()) return {};
    const int N = pre_kick_basis[0].N();
    if (N != 1) {
        throw std::runtime_error(
            "build_kicked_target: only N=1 supported for now. For N>=2 the "
            "Jacobi-Anger expansion factors per particle and requires a product "
            "of (2n_max+1) modes per axis; implement explicitly.");
    }

    const int per_basis = 2 * n_max + 1;
    TargetBasis out;
    out.reserve(pre_kick_basis.size() * per_basis);

    const double k_threshold = 1e-30;

    for (size_t j = 0; j < pre_kick_basis.size(); ++j) {
        const BasisParams& phi_j = pre_kick_basis[j];

        // precompute sums over particle axes: sum_a R_a and sum_a 1/B_aa
        Cd sum_R(0.0, 0.0);
        Cd sum_inv_B(0.0, 0.0);
        for (int a = 0; a < N; ++a) {
            sum_R += phi_j.R(a);
            sum_inv_B += 1.0 / phi_j.B(a, a);
        }

        for (int n = -n_max; n <= n_max; ++n) {
            // bessel_j handles negative n internally via J_{-n}=(-1)^n J_n
            double Jn = bessel_j(n, kappa);
            if (std::abs(Jn) < k_threshold) continue;

            // (-i)^n for n in Z via cycle: n mod 4 -> {1, -i, -1, i}
            // C++ modulo on negatives: ((n%4)+4)%4 gives Python-style residue.
            int nm = ((n % 4) + 4) % 4;
            Cd neg_i_n;
            switch (nm) {
                case 0: neg_i_n = Cd(1, 0);  break;
                case 1: neg_i_n = Cd(0, -1); break;
                case 2: neg_i_n = Cd(-1, 0); break;
                default: neg_i_n = Cd(0, 1); break;   // case 3
            }

            // ECG constant-term correction (see plan derivation):
            //   phi_kicked = phi_canon(R') * exp(R'^T B R' - R^T B R)
            //   for diagonal B: correction = exp(2 i n k_L sum_a R_a - n^2 k_L^2 sum_a 1/B_aa)
            Cd nk = static_cast<double>(n) * k_L;
            Cd correction = std::exp(Cd(0.0, 2.0) * nk * sum_R
                                     - nk * nk * sum_inv_B);

            // coefficient c_m = u_j^(0) * (-i)^n * J_n(kappa) * correction
            Cd c_m = phi_j.u * neg_i_n * Jn * correction;

            // shifted ECG: same A, B; R_a -> R_a + i*n*k_L/B_aa
            BasisParams chi;
            chi.u = c_m;
            chi.A = phi_j.A;
            chi.B = phi_j.B;
            chi.R = phi_j.R;
            for (int a = 0; a < N; ++a) {
                chi.R(a) += Cd(0.0, 1.0) * nk / phi_j.B(a, a);
            }
            chi.name = target_name_offset
                     + static_cast<int>(j) * per_basis + (n + n_max);

            out.push_back(chi);
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
// Cross-basis overlap  A(z) = <psi_tar | psi(z)>
// -----------------------------------------------------------------------------

Cd cross_overlap(const TargetBasis& target,
                 const std::vector<BasisParams>& trial,
                 const PermutationSet& perms) {
    Cd result(0.0, 0.0);
    for (const auto& chi : target) {
        Cd con_cm = std::conj(chi.u);
        for (const auto& phi : trial) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; ++p) {
                PairCache c = PairCache::build(chi, phi, perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G;
            }
            result += con_cm * phi.u * sum_p;
        }
    }
    return result;
}

Cd target_norm(const TargetBasis& target, const PermutationSet& perms) {
    // N_tar = <psi_tar|psi_tar>. overlap() regenerates its own PermutationSet
    // so `perms` is unused here; kept in signature for caller consistency.
    (void)perms;
    return overlap(target);
}

// -----------------------------------------------------------------------------
// b_alpha(z) = <d_alpha psi(z) | psi_tar>
// -----------------------------------------------------------------------------

VectorXcd assemble_b_target(const std::vector<AlphaIndex>& alpha_z_list,
                            const std::vector<BasisParams>& trial,
                            const TargetBasis& target,
                            const PermutationSet& perms) {
    const int d = static_cast<int>(alpha_z_list.size());
    VectorXcd b = VectorXcd::Zero(d);

    for (int a = 0; a < d; ++a) {
        const AlphaIndex& alpha = alpha_z_list[a];
        const BasisParams& phi_i = trial[alpha.a2];

        if (alpha.a1 == 1) {
            // d/d u_i^* of <psi(z)|psi_tar> = <phi_i|psi_tar>
            //   = sum_m c_m * sum_perm sign * M_G(phi_i, chi_m, perm)
            Cd val(0.0, 0.0);
            for (const auto& chi : target) {
                Cd sum_p(0.0, 0.0);
                for (int p = 0; p < perms.SN; ++p) {
                    PairCache c = PairCache::build(phi_i, chi, perms.matrices[p]);
                    sum_p += static_cast<double>(perms.signs[p]) * c.M_G;
                }
                val += chi.u * sum_p;
            }
            b(a) = val;
        } else {
            // alpha.a1 in {2, 3, 4}: geometric parameter of trial[i].
            // b_alpha = conj(u_i) * sum_m c_m * partial_z_addsn(
            //              alpha.a1, Real=false, pi=trial[i], pj=target[m],
            //              alpha_2=trial[i].name, alpha_3, alpha_4, perms)
            // Real=false is the BRA-side anti-holomorphic derivative d/d z_alpha^*
            // (see resolve_case_first in derivatives.cpp).
            Cd con_ui = std::conj(phi_i.u);
            Cd val(0.0, 0.0);
            for (const auto& chi : target) {
                Cd term = partial_z_addsn(alpha.a1, /*Real=*/false,
                                          phi_i, chi,
                                          phi_i.name, alpha.a3, alpha.a4,
                                          perms);
                val += chi.u * term;
            }
            b(a) = con_ui * val;
        }
    }
    return b;
}

// -----------------------------------------------------------------------------
// g_alpha = -(A / N_tar) * b_alpha
// -----------------------------------------------------------------------------

VectorXcd assemble_grad_target_fit(const std::vector<AlphaIndex>& alpha_z_list,
                                   const std::vector<BasisParams>& trial,
                                   const TargetBasis& target,
                                   const PermutationSet& perms,
                                   Cd& out_A, Cd& out_N_tar) {
    out_A = cross_overlap(target, trial, perms);
    out_N_tar = target_norm(target, perms);
    VectorXcd b = assemble_b_target(alpha_z_list, trial, target, perms);
    VectorXcd g = -(out_A / out_N_tar) * b;
    return g;
}

// -----------------------------------------------------------------------------
// Analytic u re-solve: u = S^{-1} b_target (bra-side target projection vector)
// -----------------------------------------------------------------------------

static MatrixXcd build_S_matrix(const std::vector<BasisParams>& basis,
                                const PermutationSet& perms) {
    const int K = static_cast<int>(basis.size());
    MatrixXcd S = MatrixXcd::Zero(K, K);
    for (int i = 0; i < K; ++i) {
        for (int j = i; j < K; ++j) {
            Cd s_ij(0.0, 0.0);
            for (int p = 0; p < perms.SN; ++p) {
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                s_ij += static_cast<double>(perms.signs[p]) * c.M_G;
            }
            S(i, j) = s_ij;
            S(j, i) = std::conj(s_ij);
        }
    }
    return S;
}

VectorXcd solve_optimal_u(const std::vector<BasisParams>& trial,
                          const TargetBasis& target,
                          const PermutationSet& perms) {
    const int K = static_cast<int>(trial.size());
    MatrixXcd S = build_S_matrix(trial, perms);

    // b_target[i] = sum_m c_m * sum_perm sign * M_G(phi_i, chi_m, perm)
    // (same object as assemble_b_target for a1==1, without conj(u_i) factor)
    VectorXcd b = VectorXcd::Zero(K);
    for (int i = 0; i < K; ++i) {
        Cd val(0.0, 0.0);
        for (const auto& chi : target) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; ++p) {
                PairCache c = PairCache::build(trial[i], chi, perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G;
            }
            val += chi.u * sum_p;
        }
        b(i) = val;
    }

    // Solve S u = b via SVD (robust to singular S, same as apply_analytic_kick)
    Eigen::JacobiSVD<MatrixXcd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXcd u_new = svd.solve(b);

    // Rescale u so that <psi|psi> == N_tar (the target's own norm). This is the
    // physically-correct norm to match — since |psi_tar> has norm N_tar and a
    // unitary kick preserves the pre-kick ground-state norm sqrt(pi), our trial
    // should too. Using N_tar as the fixed target avoids the previous drift
    // where norm_old = u_old^dagger * S_new * u_old shifted arbitrarily as TDVP
    // moved A/R during target fitting.
    Cd N_tar = target_norm(target, perms);
    Cd norm_new = (u_new.adjoint() * S * u_new)(0);
    if (norm_new.real() > 1e-15 && N_tar.real() > 1e-15) {
        double scale = std::sqrt(N_tar.real() / norm_new.real());
        u_new *= scale;
    }
    return u_new;
}

// -----------------------------------------------------------------------------
// Main loop: target-fitting TDVP (imaginary-time descent toward |psi_tar>)
// -----------------------------------------------------------------------------

static double compute_fidelity(const std::vector<BasisParams>& trial,
                               const TargetBasis& target,
                               Cd N_tar,
                               const PermutationSet& perms,
                               Cd* out_A = nullptr) {
    Cd A = cross_overlap(target, trial, perms);
    Cd S = overlap(trial);
    double F = (std::norm(A) / (N_tar.real() * S.real()));
    if (out_A) *out_A = A;
    return F;
}

TargetFitResult target_fitting_evolution(const std::vector<AlphaIndex>& alpha_z_list,
                                         std::vector<BasisParams> trial_init,
                                         const TargetBasis& target,
                                         const TargetFitConfig& fit_cfg) {
    const int N = trial_init[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    // Split alpha list by u_mode: active_list is what we actually step through TDVP.
    std::vector<AlphaIndex> active_list;
    active_list.reserve(alpha_z_list.size());
    for (const auto& a : alpha_z_list) {
        if (fit_cfg.u_mode == TargetFitConfig::UPDATE_RESOLVE && a.a1 == 1) continue;
        active_list.push_back(a);
    }
    const int d = static_cast<int>(active_list.size());

    // Precompute target norm (invariant across steps)
    Cd N_tar = target_norm(target, perms);
    if (fit_cfg.verbose) {
        std::cout << "[target_fit] N_tar = " << N_tar
                  << ", |target| = " << target.size()
                  << ", |alpha_active| = " << d
                  << " (mode=" << (fit_cfg.u_mode == TargetFitConfig::UPDATE_LINEAR
                                   ? "LINEAR" : "RESOLVE") << ")\n";
    }

    std::vector<BasisParams> basis = std::move(trial_init);

    // If RESOLVE, start from optimal u given the initial {A,B,R}
    if (fit_cfg.u_mode == TargetFitConfig::UPDATE_RESOLVE) {
        VectorXcd u_opt = solve_optimal_u(basis, target, perms);
        for (int i = 0; i < static_cast<int>(basis.size()); ++i) basis[i].u = u_opt(i);
    }

    TargetFitResult result;
    result.line_search_fails = 0;

    double F_pre = compute_fidelity(basis, target, N_tar, perms);
    result.F_history.push_back(F_pre);
    result.A_abs_history.push_back(std::sqrt(F_pre * N_tar.real() * overlap(basis).real()));

    double dtao = fit_cfg.dtao_init;
    double dtao_max = (fit_cfg.dtao_max > 0.0) ? fit_cfg.dtao_max : 100.0 * fit_cfg.dtao_init;

    int step = 0;
    for (; step < fit_cfg.max_steps; ++step) {
        // ---- Assemble C and g on the active parameter list ----
        MatrixXcd C_mat = assemble_C(active_list, basis);
        MatrixXcd C_bar = C_mat.conjugate().transpose();
        if (fit_cfg.lambda_C > 0.0) {
            C_bar += fit_cfg.lambda_C * MatrixXcd::Identity(d, d);
        }

        Cd A_now, N_chk;
        VectorXcd g = assemble_grad_target_fit(active_list, basis, target, perms,
                                               A_now, N_chk);

        // ---- Solve C_bar * dz = -g via SVD ----
        Eigen::JacobiSVD<MatrixXcd> svd(C_bar, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::VectorXd sv = svd.singularValues();
        double thr = fit_cfg.rcond * sv(0);
        VectorXcd Utrhs = svd.matrixU().adjoint() * (-g);
        VectorXcd sinv_utrhs(sv.size());
        for (int i = 0; i < sv.size(); ++i) {
            sinv_utrhs(i) = (sv(i) > thr) ? Utrhs(i) / sv(i) : Cd(0, 0);
        }
        VectorXcd dz = svd.matrixV() * sinv_utrhs;

        // ---- Line search: halve dtao until F improves ----
        bool accepted = false;
        double used_dtao = dtao;
        std::vector<BasisParams> trial_accept = basis;
        double F_new = F_pre;

        for (int k = 0; k < 100; ++k) {
            std::vector<BasisParams> trial = basis;
            update_basis_function(trial, dz, used_dtao, active_list, /*updata_constant=*/0);

            // If RESOLVE, re-solve u after stepping {A,B,R}
            if (fit_cfg.u_mode == TargetFitConfig::UPDATE_RESOLVE) {
                VectorXcd u_opt = solve_optimal_u(trial, target, perms);
                for (int i = 0; i < static_cast<int>(trial.size()); ++i) trial[i].u = u_opt(i);
            }

            double F_try = compute_fidelity(trial, target, N_tar, perms);
            if (F_try > F_pre) {
                trial_accept = std::move(trial);
                F_new = F_try;
                accepted = true;
                break;
            }
            used_dtao *= 0.5;
            if (used_dtao < 1e-20) break;
        }

        if (!accepted) {
            // Exploratory tiny step (mirrors tdvp_step fallback)
            double tiny = 1e-6;
            std::vector<BasisParams> trial = basis;
            update_basis_function(trial, dz, tiny, active_list, 0);
            if (fit_cfg.u_mode == TargetFitConfig::UPDATE_RESOLVE) {
                VectorXcd u_opt = solve_optimal_u(trial, target, perms);
                for (int i = 0; i < static_cast<int>(trial.size()); ++i) trial[i].u = u_opt(i);
            }
            double F_try = compute_fidelity(trial, target, N_tar, perms);
            if (F_try > F_pre) {
                trial_accept = std::move(trial);
                F_new = F_try;
                used_dtao = tiny;
                accepted = true;
            }
        }

        if (!accepted) {
            result.line_search_fails++;
            if (fit_cfg.verbose) {
                std::cout << "[target_fit] step " << step
                          << ": line search failed at F=" << std::setprecision(10) << F_pre
                          << " -> stopping\n";
            }
            break;
        }

        basis = std::move(trial_accept);
        F_pre = F_new;
        result.F_history.push_back(F_new);
        result.A_abs_history.push_back(std::abs(A_now));

        if (fit_cfg.verbose && (step % fit_cfg.print_every == 0)) {
            std::cout << "[target_fit] step " << step
                      << ": F=" << std::setprecision(10) << F_new
                      << ", 1-F=" << std::scientific << (1.0 - F_new)
                      << ", dtao=" << used_dtao
                      << std::defaultfloat << "\n";
        }

        if (1.0 - F_new < fit_cfg.fidelity_tol) {
            if (fit_cfg.verbose) {
                std::cout << "[target_fit] converged: 1-F=" << std::scientific
                          << (1.0 - F_new) << " < " << fit_cfg.fidelity_tol
                          << std::defaultfloat << " at step " << step << "\n";
            }
            step++;   // count the final accepted step
            break;
        }

        // Step size growth
        dtao = std::min(used_dtao * fit_cfg.dtao_grow, dtao_max);
    }

    result.basis = std::move(basis);
    result.F_final = F_pre;
    result.n_steps = step;

    if (fit_cfg.verbose) {
        std::cout << "[target_fit] summary: F_final=" << std::setprecision(12) << F_pre
                  << ", 1-F=" << std::scientific << (1.0 - F_pre)
                  << ", steps=" << result.n_steps
                  << ", line_search_fails=" << result.line_search_fails
                  << std::defaultfloat << "\n";
    }
    return result;
}

} // namespace ecg1d
