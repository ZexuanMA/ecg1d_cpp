#include "svm.hpp"
#include "pair_cache.hpp"
#include "interaction_kernels.hpp"
#include "physical_constants.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <Eigen/Dense>

namespace ecg1d {

std::pair<Cd, Cd> compute_HS_ij(const BasisParams& bi, const BasisParams& bj,
                                 const PermutationSet& perms,
                                 const HamiltonianTerms& terms) {
    int N = bi.N();
    Cd s_ij(0.0, 0.0);
    Cd h_ij(0.0, 0.0);

    for (int p = 0; p < perms.SN; p++) {
        double sign = static_cast<double>(perms.signs[p]);
        PairCache c = PairCache::build(bi, bj, perms.matrices[p]);

        s_ij += sign * c.M_G;

        Cd kernel(0.0, 0.0);

        if (terms.kinetic) {
            kernel += (-hbar * hbar / (2.0 * mass)) * compute_P_Mij(c);
        }
        if (terms.harmonic) {
            kernel += (mass * omega * omega / 2.0) * compute_rTr_Mij(c);
        }
        if (terms.delta && N >= 2) {
            Cd d_sum(0.0, 0.0);
            for (int a = 0; a < N; a++)
                for (int b = a + 1; b < N; b++)
                    d_sum += compute_G_Mijab(c, a, b);
            kernel += g_contact * d_sum;
        }
        if (terms.gaussian && N >= 2) {
            Cd g_sum(0.0, 0.0);
            for (int a = 0; a < N; a++)
                for (int b = a + 1; b < N; b++)
                    g_sum += compute_H_Mijab(c, a, b);
            kernel += g_gauss * g_sum;
        }
        if (terms.kicking) {
            Cd q_sum(0.0, 0.0);
            for (int a = 0; a < N; a++)
                q_sum += compute_Q_Mija(c, a);
            kernel += hbar * kappa * q_sum;
        }

        h_ij += sign * c.M_G * kernel;
    }

    return {h_ij, s_ij};
}

std::pair<MatrixXcd, MatrixXcd> build_HS(const std::vector<BasisParams>& basis,
                                          const PermutationSet& perms,
                                          const HamiltonianTerms& terms) {
    int n = static_cast<int>(basis.size());
    MatrixXcd H = MatrixXcd::Zero(n, n);
    MatrixXcd S = MatrixXcd::Zero(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            auto [h, s] = compute_HS_ij(basis[i], basis[j], perms, terms);
            H(i, j) = h;  H(j, i) = std::conj(h);
            S(i, j) = s;  S(j, i) = std::conj(s);
        }
    }
    return {H, S};
}

// Check if column k of overlap matrix S has excessive overlap with any other column.
// Returns true if max |S(i,k)| / sqrt(|S(i,i)*S(k,k)|) > threshold.
bool has_excessive_overlap(const MatrixXcd& S, int k, double threshold) {
    int n = S.rows();
    double s_kk = std::abs(S(k, k).real());
    if (s_kk < 1e-15) return true;  // degenerate basis function

    for (int i = 0; i < n; i++) {
        if (i == k) continue;
        double s_ii = std::abs(S(i, i).real());
        if (s_ii < 1e-15) continue;
        double overlap_norm = std::abs(S(i, k)) / std::sqrt(s_ii * s_kk);
        if (overlap_norm > threshold) return true;
    }
    return false;
}

// Core eigenvalue solver with truncation. Returns full result with variance.
EigenResult lowest_energy_full(const MatrixXcd& H, const MatrixXcd& S,
                                double max_cond, double E_lower_bound,
                                double rcond_trunc) {
    const double INF = std::numeric_limits<double>::infinity();
    EigenResult bad;
    bad.energy = INF;

    int n = H.rows();

    // Symmetrize and take real parts
    Eigen::MatrixXd Hs = (0.5 * (H + H.adjoint())).real();
    Eigen::MatrixXd Ss = (0.5 * (S + S.adjoint())).real();

    // Eigendecompose S
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ss);
    if (es.info() != Eigen::Success) return bad;

    Eigen::VectorXd w = es.eigenvalues();
    Eigen::MatrixXd v = es.eigenvectors();
    double w_max = w(n - 1);

    // S must have at least some positive eigenvalues
    if (w_max < 1e-15) return bad;

    // Truncation: keep only eigenvectors with w_i > w_max * rcond_trunc
    double w_cutoff = w_max * rcond_trunc;
    int n_keep = 0;
    for (int i = 0; i < n; i++) {
        if (w(i) > w_cutoff) n_keep++;
    }
    if (n_keep == 0) return bad;

    // Check condition number of kept subspace
    int i_start = n - n_keep;  // eigenvalues are sorted ascending
    double w_min_kept = w(i_start);
    if (w_max / w_min_kept > max_cond) return bad;

    // Build S^{-1/2} in the truncated subspace
    // V_keep: n x n_keep matrix of kept eigenvectors
    Eigen::MatrixXd V_keep = v.rightCols(n_keep);
    Eigen::VectorXd w_keep = w.tail(n_keep);
    Eigen::VectorXd w_inv_sqrt = w_keep.array().rsqrt();
    Eigen::MatrixXd S_inv_half = V_keep * w_inv_sqrt.asDiagonal() * V_keep.transpose();

    // H_tilde = S^{-1/2} H S^{-1/2} in truncated space
    Eigen::MatrixXd H_tilde = S_inv_half * Hs * S_inv_half;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eh(H_tilde);
    if (eh.info() != Eigen::Success) return bad;

    double E0 = eh.eigenvalues()(0);

    // Reject unphysical energies
    if (E0 < E_lower_bound) return bad;

    // Cross-validation: back-transform eigenvector and check Rayleigh quotient
    Eigen::VectorXd c_tilde = eh.eigenvectors().col(0);
    Eigen::VectorXd c = S_inv_half * c_tilde;
    double num = c.transpose() * Hs * c;
    double den = c.transpose() * Ss * c;
    if (den < 1e-15) return bad;
    double E_check = num / den;

    // If eigenvalue and Rayleigh quotient disagree, untrustworthy
    if (std::abs(E_check - E0) > 1e-6 * std::abs(E0) + 1e-12)
        return bad;

    // Compute energy variance: σ² = c^T H_tilde^2 c - E0^2
    // (c_tilde is normalized in the transformed space)
    Eigen::VectorXd Hc = H_tilde * c_tilde;
    double H2_expect = Hc.squaredNorm();
    double variance = H2_expect - E0 * E0;

    EigenResult result;
    result.energy = E0;
    result.variance = variance;
    result.n_kept = n_keep;
    return result;
}

double lowest_energy(const MatrixXcd& H, const MatrixXcd& S,
                     double max_cond, double E_lower_bound,
                     double rcond_trunc) {
    return lowest_energy_full(H, S, max_cond, E_lower_bound, rcond_trunc).energy;
}

void set_u_from_eigenvector(std::vector<BasisParams>& basis,
                            const MatrixXcd& H, const MatrixXcd& S,
                            double max_cond, double rcond_trunc) {
    int n = H.rows();

    Eigen::MatrixXd Hs = (0.5 * (H + H.adjoint())).real();
    Eigen::MatrixXd Ss = (0.5 * (S + S.adjoint())).real();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ss);
    Eigen::VectorXd w = es.eigenvalues();
    Eigen::MatrixXd v = es.eigenvectors();
    double w_max = w(n - 1);
    if (w_max < 1e-15) return;

    // Truncation: keep eigenvectors with w_i > w_max * rcond_trunc
    double w_cutoff = w_max * rcond_trunc;
    int n_keep = 0;
    for (int i = 0; i < n; i++) {
        if (w(i) > w_cutoff) n_keep++;
    }
    if (n_keep == 0) return;

    int i_start = n - n_keep;
    double w_min_kept = w(i_start);
    if (w_max / w_min_kept > max_cond) return;

    Eigen::MatrixXd V_keep = v.rightCols(n_keep);
    Eigen::VectorXd w_keep = w.tail(n_keep);
    Eigen::VectorXd w_inv_sqrt = w_keep.array().rsqrt();
    Eigen::MatrixXd S_inv_half = V_keep * w_inv_sqrt.asDiagonal() * V_keep.transpose();

    Eigen::MatrixXd H_tilde = S_inv_half * Hs * S_inv_half;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eh(H_tilde);
    Eigen::VectorXd c = S_inv_half * eh.eigenvectors().col(0);

    // Normalize so max |c_i| = 1
    double max_c = c.array().abs().maxCoeff();
    if (max_c > 0) c /= max_c;

    for (int i = 0; i < n; i++) {
        basis[i].u = Cd(c(i), 0.0);
    }
}

BasisParams random_basis_2particle(std::mt19937_64& rng, int name) {
    int N = 2;
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    std::normal_distribution<double> normal01(0.0, 1.0);

    // CM and relative widths (log-uniform sampling)
    double a_plus = std::exp(std::log(0.05) + uniform01(rng) * (std::log(8.0) - std::log(0.05)));
    double a_minus = std::exp(std::log(0.005) + uniform01(rng) * (std::log(100.0) - std::log(0.005)));

    double a = 0.5 * (a_plus + a_minus);
    double c = 0.5 * (a_plus - a_minus);

    MatrixXcd A = MatrixXcd::Zero(N, N);

    // 30% chance of asymmetric A
    if (uniform01(rng) < 0.3) {
        double asym = 0.7 + uniform01(rng) * 0.6;  // uniform(0.7, 1.3)
        A(0, 0) = Cd(a * asym, 0.0);
        A(1, 1) = Cd(a / asym, 0.0);
        A(0, 1) = Cd(c, 0.0);
        A(1, 0) = Cd(c, 0.0);
    } else {
        A(0, 0) = Cd(a, 0.0);
        A(1, 1) = Cd(a, 0.0);
        A(0, 1) = Cd(c, 0.0);
        A(1, 0) = Cd(c, 0.0);
    }

    // Ensure A is positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_a(A.real());
    double w_a_min = es_a.eigenvalues()(0);
    if (w_a_min < 0.01) {
        double shift = 0.01 - w_a_min + 0.005;
        A(0, 0) += Cd(shift, 0.0);
        A(1, 1) += Cd(shift, 0.0);
    }

    // B: small positive diagonal (log-uniform)
    MatrixXcd B = MatrixXcd::Zero(N, N);
    for (int i = 0; i < N; i++) {
        double b_val = std::exp(std::log(0.005) + uniform01(rng) * (std::log(3.0) - std::log(0.005)));
        B(i, i) = Cd(b_val, 0.0);
    }

    // Ensure A+B positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_ab((A + B).real());
    double w_ab_min = es_ab.eigenvalues()(0);
    if (w_ab_min < 0.05) {
        double shift = 0.05 - w_ab_min + 0.01;
        A(0, 0) += Cd(shift, 0.0);
        A(1, 1) += Cd(shift, 0.0);
    }

    // R: Gaussian centers
    VectorXcd R(N);
    for (int i = 0; i < N; i++) {
        R(i) = Cd(0.5 * normal01(rng), 0.0);
    }

    return BasisParams::from_arrays(Cd(1.0, 0.0), A, B, R, name);
}

BasisParams perturb_basis(const BasisParams& base, std::mt19937_64& rng,
                          double scale, int N) {
    std::normal_distribution<double> normal01(0.0, 1.0);

    // Perturb A
    MatrixXcd A_new = base.A;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A_new(i, j) += Cd(scale * normal01(rng), 0.0);
    // Symmetrize
    A_new = 0.5 * (A_new + A_new.transpose());

    // Ensure A positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_a(A_new.real());
    double w_a_min = es_a.eigenvalues()(0);
    if (w_a_min < 0.01) {
        double shift = 0.01 - w_a_min + 0.005;
        for (int i = 0; i < N; i++)
            A_new(i, i) += Cd(shift, 0.0);
    }

    // Perturb B
    MatrixXcd B_new = base.B;
    for (int i = 0; i < N; i++) {
        double b_val = std::max(0.005, base.B(i, i).real() + scale * normal01(rng));
        B_new(i, i) = Cd(b_val, 0.0);
    }

    // Perturb R
    VectorXcd R_new = base.R;
    for (int i = 0; i < N; i++) {
        R_new(i) += Cd(scale * 0.3 * normal01(rng), 0.0);
    }

    // Ensure A+B positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_ab((A_new + B_new).real());
    double w_ab_min = es_ab.eigenvalues()(0);
    if (w_ab_min < 0.05) {
        double shift = 0.05 - w_ab_min + 0.01;
        for (int i = 0; i < N; i++)
            A_new(i, i) += Cd(shift, 0.0);
    }

    return BasisParams::from_arrays(Cd(1.0, 0.0), A_new, B_new, R_new, base.name);
}

SvmResult svm_build_basis(int N, int K_max, int n_trials,
                           const HamiltonianTerms& terms,
                           int seed,
                           double E_lower_bound) {
    std::mt19937_64 rng(seed);
    PermutationSet perms = PermutationSet::generate(N);

    double E_lower = E_lower_bound;

    // Initial basis: harmonic oscillator ground state
    MatrixXcd A0 = MatrixXcd::Zero(N, N);
    for (int i = 0; i < N; i++) A0(i, i) = Cd(0.5, 0.0);
    MatrixXcd B0 = MatrixXcd::Zero(N, N);
    VectorXcd R0 = VectorXcd::Zero(N);

    BasisParams b0 = BasisParams::from_arrays(Cd(1.0, 0.0), A0, B0, R0, 0);
    std::vector<BasisParams> basis = {b0};

    auto [H, S] = build_HS(basis, perms, terms);
    auto res0 = lowest_energy_full(H, S, 1e8, E_lower);
    std::cout << "  K= 1: E = " << std::setprecision(10) << res0.energy
              << "  sigma=" << std::scientific << std::sqrt(std::max(0.0, res0.variance))
              << std::defaultfloat << std::endl;

    for (int k = 1; k < K_max; k++) {
        double best_E = std::numeric_limits<double>::infinity();
        BasisParams best_basis = b0;
        MatrixXcd best_H, best_S;
        int n = static_cast<int>(basis.size());

        for (int trial = 0; trial < n_trials; trial++) {
            BasisParams tb = random_basis_2particle(rng, k);

            // Extend H and S by one row/column
            MatrixXcd H_ext = MatrixXcd::Zero(n + 1, n + 1);
            MatrixXcd S_ext = MatrixXcd::Zero(n + 1, n + 1);
            H_ext.topLeftCorner(n, n) = H;
            S_ext.topLeftCorner(n, n) = S;

            for (int i = 0; i < n; i++) {
                auto [h, s] = compute_HS_ij(basis[i], tb, perms, terms);
                H_ext(i, n) = h;  H_ext(n, i) = std::conj(h);
                S_ext(i, n) = s;  S_ext(n, i) = std::conj(s);
            }

            auto [h_nn, s_nn] = compute_HS_ij(tb, tb, perms, terms);
            H_ext(n, n) = h_nn;
            S_ext(n, n) = s_nn;

            // Overlap screening: reject trial if too similar to existing basis
            if (has_excessive_overlap(S_ext, n, 0.99)) continue;

            double E_trial = lowest_energy(H_ext, S_ext, 1e8, E_lower);

            if (E_trial < best_E) {
                best_E = E_trial;
                best_basis = tb;
                best_H = H_ext;
                best_S = S_ext;
            }
        }

        if (best_E == std::numeric_limits<double>::infinity()) {
            std::cout << "  K=" << std::setw(2) << (k + 1)
                      << ": no valid trial found, stopping." << std::endl;
            break;
        }

        basis.push_back(best_basis);
        H = best_H;
        S = best_S;

        auto res = lowest_energy_full(H, S, 1e8, E_lower);
        std::cout << "  K=" << std::setw(2) << (k + 1)
                  << ": E = " << std::setprecision(10) << best_E
                  << "  sigma=" << std::scientific << std::sqrt(std::max(0.0, res.variance))
                  << "  n_kept=" << res.n_kept
                  << std::defaultfloat << std::endl;
    }

    return {basis, H, S};
}

SvmResult stochastic_refine(std::vector<BasisParams> basis,
                             MatrixXcd H, MatrixXcd S,
                             const PermutationSet& perms,
                             int N, const HamiltonianTerms& terms,
                             int n_trials, int max_rounds,
                             int seed,
                             double E_lower_bound) {
    std::mt19937_64 rng(seed);
    int n = static_cast<int>(basis.size());

    double E_current = lowest_energy(H, S, 1e8, E_lower_bound);
    double scale0 = 0.5;

    // Max total energy drop allowed per round (across all K replacements).
    // Anchored at start of round to prevent cascading drift.
    // Legitimate total improvement per round is typically 1e-5 to 1e-7.
    // Setting tight to prevent ghost states from sneaking through.
    double max_round_drop = 1e-4;

    std::cout << "  Refine: K=" << n << ", initial E=" << std::setprecision(10)
              << E_current << ", E_lower_bound=" << E_lower_bound
              << ", max_round_drop=" << max_round_drop << std::endl;

    for (int round_idx = 0; round_idx < max_rounds; round_idx++) {
        bool improved = false;
        double scale = scale0 / (1.0 + round_idx * 0.1);
        int n_replaced = 0;

        // Anchor: energy at start of this round. No replacement in this round
        // can push energy below E_round_start - max_round_drop.
        double E_round_start = E_current;
        double E_floor = std::max(E_lower_bound, E_round_start - max_round_drop);

        for (int k = 0; k < n; k++) {
            double best_E = E_current;
            BasisParams best_basis_k = basis[k];
            MatrixXcd best_H = H;
            MatrixXcd best_S = S;
            bool found = false;

            std::uniform_real_distribution<double> uniform01(0.0, 1.0);

            for (int trial = 0; trial < n_trials; trial++) {
                BasisParams tb = (uniform01(rng) < 0.5)
                    ? random_basis_2particle(rng, k)
                    : perturb_basis(basis[k], rng, scale, N);

                // Rebuild row/col k
                MatrixXcd H_new = H;
                MatrixXcd S_new = S;
                for (int i = 0; i < n; i++) {
                    if (i == k) {
                        auto [h, s] = compute_HS_ij(tb, tb, perms, terms);
                        H_new(k, k) = h;
                        S_new(k, k) = s;
                    } else {
                        auto [h, s] = compute_HS_ij(basis[i], tb, perms, terms);
                        H_new(i, k) = h;  H_new(k, i) = std::conj(h);
                        S_new(i, k) = s;  S_new(k, i) = std::conj(s);
                    }
                }

                // Overlap screening: reject if trial too similar to existing basis
                if (has_excessive_overlap(S_new, k, 0.95)) continue;

                // Stricter eigenvalue truncation for refinement (rcond=1e-8)
                double E_trial = lowest_energy(H_new, S_new, 1e8, E_floor, 1e-8);

                if (E_trial < best_E) {
                    best_E = E_trial;
                    best_basis_k = tb;
                    best_H = H_new;
                    best_S = S_new;
                    found = true;
                }
            }

            if (found) {
                basis[k] = best_basis_k;
                H = best_H;
                S = best_S;
                improved = true;
                E_current = best_E;
                n_replaced++;
            }
        }

        auto res = lowest_energy_full(H, S, 1e8, E_lower_bound);
        std::cout << "  Round " << (round_idx + 1) << " done: E = "
                  << std::setprecision(10) << E_current
                  << "  sigma=" << std::scientific << std::sqrt(std::max(0.0, res.variance))
                  << "  replaced=" << n_replaced
                  << std::defaultfloat << std::endl;

        if (!improved) {
            std::cout << "  No improvement in round " << (round_idx + 1)
                      << ". Stopping." << std::endl;
            break;
        }
    }

    return {basis, H, S};
}

} // namespace ecg1d
