#include "kick_operator.hpp"
#include "pair_cache.hpp"
#include "svm.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace ecg1d {

// Single-particle kick kernel for particle a:
// kick_a = sum_{n=-N}^{N} (-i)^n J_n(kappa) exp(-n^2 k_L^2 K_inv(a,a) + i*n*2*k_L*mu(a))
//
// Using symmetry J_{-n}(x) = (-1)^n J_n(x), rewrite as:
//   kick_a = J_0(kappa) * exp(0) + sum_{n=1}^{N} J_n(kappa) * [(-i)^n e^{...n} + (-1)^n (-i)^{-n} e^{...-n}]
//          = J_0(kappa) + sum_{n=1}^{N} J_n(kappa) * [(-i)^n e^{inq-n^2s} + i^n e^{-inq-n^2s}]
// where q = 2*k_L*mu(a), s = k_L^2*K_inv(a,a)
static Cd single_particle_kick_kernel(const PairCache& c, int a,
                                       double kappa, double k_L,
                                       int n_max) {
    Cd s = k_L * k_L * c.K_inv(a, a);   // gaussian damping factor
    Cd q = 2.0 * k_L * c.mu(a);          // phase factor

    // n=0 term: (-i)^0 J_0(kappa) exp(0) = J_0(kappa)
    Cd result = std::cyl_bessel_j(0, kappa);

    // Precompute (-i)^n for n=1,2,3,4 and cycle: (-i)^1=-i, (-i)^2=-1, (-i)^3=i, (-i)^4=1
    Cd neg_i(0.0, -1.0);
    Cd neg_i_power(1.0, 0.0);  // will be (-i)^n

    for (int n = 1; n <= n_max; n++) {
        neg_i_power *= neg_i;  // (-i)^n

        double Jn = std::cyl_bessel_j(n, kappa);
        if (std::abs(Jn) < 1e-30) continue;  // skip negligible terms

        // Gaussian damping: exp(-n^2 * s)
        Cd damping = std::exp(-static_cast<double>(n * n) * s);

        // Positive n term: (-i)^n J_n exp(-n^2 s + i n q)
        Cd pos_term = neg_i_power * Jn * damping * std::exp(Cd(0, 1) * static_cast<double>(n) * q);

        // Negative n term: (-i)^{-n} J_{-n} exp(-n^2 s - i n q)
        // J_{-n}(x) = (-1)^n J_n(x), (-i)^{-n} = i^n = conj((-i)^n) for real n
        Cd neg_term = std::conj(neg_i_power) * std::pow(-1.0, n) * Jn *
                      damping * std::exp(Cd(0, -1) * static_cast<double>(n) * q);

        result += pos_term + neg_term;
    }
    return result;
}

MatrixXcd build_kick_matrix(const std::vector<BasisParams>& basis,
                            const PermutationSet& perms,
                            double kappa, double k_L,
                            int n_bessel) {
    int K = static_cast<int>(basis.size());
    int N = basis[0].N();
    MatrixXcd Kick = MatrixXcd::Zero(K, K);

    // Kick operator is UNITARY, not Hermitian: K(j,i) != conj(K(i,j))
    // Must compute all K×K elements separately
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            Cd kick_ij(0.0, 0.0);

            for (int p = 0; p < perms.SN; p++) {
                double sign = static_cast<double>(perms.signs[p]);
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);

                // Total kick kernel = product over all particles
                Cd kick_kernel(1.0, 0.0);
                for (int a = 0; a < N; a++) {
                    kick_kernel *= single_particle_kick_kernel(c, a, kappa, k_L, n_bessel);
                }

                kick_ij += sign * c.M_G * kick_kernel;
            }

            Kick(i, j) = kick_ij;
        }
    }

    return Kick;
}

void apply_analytic_kick(std::vector<BasisParams>& basis,
                         const PermutationSet& perms,
                         double kappa, double k_L,
                         int n_bessel) {
    int K = static_cast<int>(basis.size());

    // Build kick matrix and overlap matrix
    MatrixXcd Kick = build_kick_matrix(basis, perms, kappa, k_L, n_bessel);

    // Build overlap matrix S
    MatrixXcd S = MatrixXcd::Zero(K, K);
    for (int i = 0; i < K; i++) {
        for (int j = i; j < K; j++) {
            Cd s_ij(0.0, 0.0);
            for (int p = 0; p < perms.SN; p++) {
                double sign = static_cast<double>(perms.signs[p]);
                PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                s_ij += sign * c.M_G;
            }
            S(i, j) = s_ij;
            S(j, i) = std::conj(s_ij);
        }
    }

    // Current u vector
    VectorXcd u(K);
    for (int i = 0; i < K; i++) u(i) = basis[i].u;

    // Diagnostics: norm before kick
    Cd norm_before = (u.adjoint() * S * u)(0);

    // Apply kick: u' = S^{-1} K u
    // Solve S u' = K u via SVD for robustness
    VectorXcd Ku = Kick * u;
    Eigen::JacobiSVD<MatrixXcd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXcd u_new = svd.solve(Ku);

    // Diagnostics: norm after kick
    Cd norm_after = (u_new.adjoint() * S * u_new)(0);

    // Check unitarity: kick should preserve norm
    // |<ψ'|ψ'>| should equal |<ψ|ψ>| since U is unitary
    std::cout << "  [kick] norm before=" << std::setprecision(6) << norm_before.real()
              << ", norm after=" << norm_after.real()
              << ", ratio=" << norm_after.real() / norm_before.real()
              << ", |u|=" << u.norm() << " -> " << u_new.norm()
              << ", S cond=" << std::scientific
              << svd.singularValues()(0) / svd.singularValues()(K-1)
              << std::defaultfloat << std::endl;

    // Renormalize u so that u'†S u' = u†S u (preserve norm)
    // This doesn't change the physical wave function (just a global factor),
    // but prevents gradient explosion in subsequent TDVP evolution.
    if (norm_after.real() > 1e-15 && norm_before.real() > 1e-15) {
        double scale = std::sqrt(norm_before.real() / norm_after.real());
        u_new *= scale;
    }

    // Update basis coefficients
    for (int i = 0; i < K; i++) {
        basis[i].u = u_new(i);
    }
}

void free_evolve_fixed_basis(std::vector<BasisParams>& basis,
                              const MatrixXcd& H, const MatrixXcd& S,
                              double T_duration) {
    int K = static_cast<int>(basis.size());

    // Solve generalized eigenvalue problem: H v = E S v
    // via S^{-1/2} transformation
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es_S(S);
    if (es_S.info() != Eigen::Success) {
        std::cerr << "free_evolve_fixed_basis: S eigendecomposition failed" << std::endl;
        return;
    }

    Eigen::VectorXd w = es_S.eigenvalues();
    double w_max = w(K - 1);
    double w_cutoff = w_max * 1e-10;

    // Count kept eigenvalues
    int n_keep = 0;
    for (int i = 0; i < K; i++) {
        if (w(i) > w_cutoff) n_keep++;
    }

    int i_start = K - n_keep;
    MatrixXcd V_keep = es_S.eigenvectors().rightCols(n_keep);
    Eigen::VectorXd w_keep = w.tail(n_keep);

    // S^{-1/2} in kept subspace
    Eigen::VectorXd w_inv_sqrt = w_keep.array().rsqrt();
    MatrixXcd S_inv_half = V_keep * w_inv_sqrt.asDiagonal() * V_keep.adjoint();

    // Transform H: H_tilde = S^{-1/2} H S^{-1/2}
    MatrixXcd H_tilde = S_inv_half * H * S_inv_half;
    // Symmetrize
    H_tilde = (0.5 * (H_tilde + H_tilde.adjoint())).eval();

    // Diagonalize H_tilde
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es_H(H_tilde);
    if (es_H.info() != Eigen::Success) {
        std::cerr << "free_evolve_fixed_basis: H_tilde eigendecomposition failed" << std::endl;
        return;
    }

    Eigen::VectorXd eigenvalues = es_H.eigenvalues();   // real
    MatrixXcd eigenvectors = es_H.eigenvectors();         // in S^{-1/2} space

    // Current u in original space
    VectorXcd u(K);
    for (int i = 0; i < K; i++) u(i) = basis[i].u;

    // Transform u to S^{-1/2} space: u_tilde = S^{1/2} u
    MatrixXcd S_half = V_keep * w_keep.array().sqrt().matrix().asDiagonal() * V_keep.adjoint();
    VectorXcd u_tilde = S_half * u;

    // Expand in eigenbasis of H_tilde: u_tilde = Σ c_n v_n
    VectorXcd c = eigenvectors.adjoint() * u_tilde;

    // Time evolve: c_n(t) = c_n(0) * exp(-i E_n t)
    for (int n = 0; n < n_keep; n++) {
        c(n) *= std::exp(Cd(0, -1) * eigenvalues(n) * T_duration);
    }

    // Back-transform: u_tilde(t) = Σ c_n(t) v_n
    VectorXcd u_tilde_new = eigenvectors * c;

    // Back to original space: u = S^{-1/2} u_tilde
    VectorXcd u_new = S_inv_half * u_tilde_new;

    // Update basis
    for (int i = 0; i < K; i++) {
        basis[i].u = u_new(i);
    }
}

std::vector<BasisParams> augment_basis_with_momentum(
    const std::vector<BasisParams>& basis,
    double k_L, int n_mom, double b_val) {

    int N = basis[0].N();
    std::vector<BasisParams> augmented = basis;

    // Candidate momentum functions: various widths × various momenta
    std::vector<double> widths = {0.15, 0.3, 0.5, 0.8, 1.2, 2.0, 3.5, 6.0, 10.0};
    std::vector<BasisParams> candidates;

    for (double w : widths) {
        for (int n = -n_mom; n <= n_mom; n++) {
            if (n == 0) continue;

            MatrixXcd A_new = MatrixXcd::Zero(N, N);
            MatrixXcd B_new = MatrixXcd::Zero(N, N);
            VectorXcd R_new = VectorXcd::Zero(N);

            for (int a = 0; a < N; a++) {
                A_new(a, a) = Cd(w, 0.0);
                B_new(a, a) = Cd(b_val, 0.0);
                double mom = static_cast<double>(n) * k_L;
                R_new(a) = Cd(0.0, mom / b_val);
            }
            if (N >= 2) {
                A_new(0, 1) = Cd(0.1 * w, 0.0);
                A_new(1, 0) = Cd(0.1 * w, 0.0);
            }

            candidates.push_back(BasisParams::from_arrays(
                Cd(0.0, 0.0), A_new, B_new, R_new, 1000 * n));
        }
    }

    // Greedy selection: add candidates one-by-one, checking S condition number.
    // Only keep those that don't make S too ill-conditioned.
    double max_cond = 1e6;  // keep S well-conditioned for accurate kick projection
    PermutationSet perms = PermutationSet::generate(N);

    int added = 0;
    for (auto& cand : candidates) {
        std::vector<BasisParams> trial = augmented;
        trial.push_back(cand);

        // Quick S conditioning check
        int Kt = static_cast<int>(trial.size());
        MatrixXcd St = MatrixXcd::Zero(Kt, Kt);
        for (int i = 0; i < Kt; i++) {
            for (int j = i; j < Kt; j++) {
                Cd s_ij(0.0, 0.0);
                for (int p = 0; p < perms.SN; p++) {
                    double sign = static_cast<double>(perms.signs[p]);
                    PairCache c = PairCache::build(trial[i], trial[j], perms.matrices[p]);
                    s_ij += sign * c.M_G;
                }
                St(i, j) = s_ij;
                St(j, i) = std::conj(s_ij);
            }
        }

        Eigen::SelfAdjointEigenSolver<MatrixXcd> es(St);
        double w_min = es.eigenvalues()(0);
        double w_max = es.eigenvalues()(Kt - 1);
        if (w_min > 0 && w_max / w_min < max_cond) {
            augmented.push_back(cand);
            added++;
        }
    }

    std::cout << "  Augmented basis: " << basis.size() << " original + "
              << added << " momentum (from " << candidates.size() << " candidates, cond<"
              << std::scientific << max_cond << std::defaultfloat << ") = "
              << augmented.size() << " total" << std::endl;

    return augmented;
}

} // namespace ecg1d
