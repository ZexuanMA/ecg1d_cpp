// Step B: Orthogonalized subspace prototype for kicked gamma=0
// This is a standalone test to answer: "Is the Gaussian subspace itself adequate,
// or is the greedy basis selection the bottleneck?"
//
// Strategy:
// 1. Build ground state basis (same as --kicked-gamma0)
// 2. Generate ALL momentum candidates (no greedy filtering)
// 3. Combine into one big basis, build full S, H, K matrices
// 4. SVD-truncate S to get orthonormal subspace
// 5. Transform H, K into this orthonormal subspace
// 6. Do kicked evolution in orthonormal subspace (no cond(S) problem)
// 7. Compare with exact

#include "svm.hpp"
#include "tdvp_solver.hpp"
#include "hamiltonian.hpp"
#include "kicked_exact.hpp"
#include "kick_operator.hpp"
#include "observables.hpp"
#include "physical_constants.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace ecg1d;
using Eigen::VectorXd;

// Silence cout during setup
class NullBuffer : public std::streambuf {
protected:
    int overflow(int c) override { return c; }
};

class ScopedSilence {
public:
    explicit ScopedSilence(bool enabled) : enabled_(enabled) {
        if (enabled_) { old_ = std::cout.rdbuf(&null_); }
    }
    ~ScopedSilence() {
        if (enabled_) { std::cout.rdbuf(old_); }
    }
private:
    bool enabled_ = false;
    NullBuffer null_;
    std::streambuf* old_ = nullptr;
};

// Generate ALL momentum candidates (no filtering)
static std::vector<BasisParams> generate_all_momentum_candidates(
    const std::vector<BasisParams>& gs_basis,
    double k_L_val, int n_mom, double b_val)
{
    int N = gs_basis[0].N();

    // Collect width candidates
    std::vector<double> widths = {0.3, 0.8, 2.0, 6.0};
    for (const auto& bp : gs_basis) {
        for (int a = 0; a < N; a++) {
            double a_diag = bp.A(a, a).real();
            if (a_diag > 0.05 && a_diag < 20.0) {
                bool dup = false;
                for (double w : widths) {
                    if (std::abs(a_diag - w) / w < 0.1) { dup = true; break; }
                }
                if (!dup) widths.push_back(a_diag);
            }
        }
    }
    std::sort(widths.begin(), widths.end());

    std::vector<BasisParams> candidates;

    if (N == 1) {
        for (double w : widths) {
            for (int n = -n_mom; n <= n_mom; n++) {
                if (n == 0) continue;
                MatrixXcd A_new = MatrixXcd::Zero(1, 1);
                MatrixXcd B_new = MatrixXcd::Zero(1, 1);
                VectorXcd R_new = VectorXcd::Zero(1);
                A_new(0, 0) = Cd(w, 0.0);
                B_new(0, 0) = Cd(b_val, 0.0);
                R_new(0) = Cd(0.0, static_cast<double>(n) * k_L_val / b_val);
                candidates.push_back(BasisParams::from_arrays(
                    Cd(0.0, 0.0), A_new, B_new, R_new, 1000 * n));
            }
        }
    } else {
        for (double w : widths) {
            for (int m1 = -n_mom; m1 <= n_mom; m1++) {
                for (int m2 = m1; m2 <= n_mom; m2++) {
                    if (m1 == 0 && m2 == 0) continue;

                    MatrixXcd A_new = MatrixXcd::Zero(N, N);
                    MatrixXcd B_new = MatrixXcd::Zero(N, N);
                    VectorXcd R_new = VectorXcd::Zero(N);

                    for (int a = 0; a < N; a++) {
                        A_new(a, a) = Cd(w, 0.0);
                        B_new(a, a) = Cd(b_val, 0.0);
                    }
                    if (N >= 2) {
                        A_new(0, 1) = Cd(0.1 * w, 0.0);
                        A_new(1, 0) = Cd(0.1 * w, 0.0);
                    }

                    R_new(0) = Cd(0.0, static_cast<double>(m1) * k_L_val / b_val);
                    R_new(1) = Cd(0.0, static_cast<double>(m2) * k_L_val / b_val);

                    int label = 1000 * m1 + m2;
                    candidates.push_back(BasisParams::from_arrays(
                        Cd(0.0, 0.0), A_new, B_new, R_new, label));
                }
            }
        }
    }

    return candidates;
}

int main(int argc, char* argv[]) {
    std::cout << std::setprecision(10);

    // --- Parameters ---
    int N = 2;
    int K_gs = 10;
    int n_mom = 4;
    double b_val = 0.5;
    double kappa_val = kappa;
    double k_L_val = k_L;
    double T_period = 2.0 * M_PI / omega;
    int n_kicks = 50;
    double svd_rcond = 1e-8;  // SVD truncation threshold

    // Parse optional args
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--n-kicks" && i + 1 < argc) n_kicks = std::stoi(argv[++i]);
        if (arg == "--n-mom" && i + 1 < argc) n_mom = std::stoi(argv[++i]);
        if (arg == "--svd-rcond" && i + 1 < argc) svd_rcond = std::stod(argv[++i]);
        if (arg == "--k-gs" && i + 1 < argc) K_gs = std::stoi(argv[++i]);
    }

    std::cout << "=== Step B: Orthogonalized Subspace Test ===" << std::endl;
    std::cout << "N=" << N << ", K_gs=" << K_gs << ", n_mom=" << n_mom
              << ", b_val=" << b_val << ", svd_rcond=" << std::scientific << svd_rcond
              << std::defaultfloat << std::endl;
    std::cout << "kappa=" << kappa_val << ", k_L=" << k_L_val
              << ", T=" << T_period << ", n_kicks=" << n_kicks << std::endl;

    // --- Phase 1: Ground state basis ---
    HamiltonianTerms terms_free;
    terms_free.kinetic  = true;
    terms_free.harmonic = true;
    terms_free.delta    = false;
    terms_free.gaussian = false;
    terms_free.kicking  = false;

    double E_lower = N * 0.5;
    PermutationSet perms = PermutationSet::generate(N);

    auto svm = [&]() {
        ScopedSilence silence(true);
        return svm_build_basis(N, K_gs, 8000, terms_free, 42, E_lower);
    }();

    auto refined = [&]() {
        ScopedSilence silence(true);
        return stochastic_refine(svm.basis, svm.H, svm.S,
                                 perms, N, terms_free, 500, 5, 123, E_lower);
    }();

    // Skip TDVP polish (build_alpha_z_list is in main.cpp); refined basis is good enough
    set_u_from_eigenvector(refined.basis, refined.H, refined.S);
    double E_polished = lowest_energy(refined.H, refined.S);
    std::cout << "\nGround state E = " << E_polished
              << " (exact = " << N * 0.5 << ")" << std::endl;

    // --- Phase 2: Generate ALL candidates ---
    auto candidates = generate_all_momentum_candidates(refined.basis, k_L_val, n_mom, b_val);
    std::cout << "Candidate pool: " << candidates.size()
              << " momentum candidates" << std::endl;

    // Combine: ground state basis + ALL candidates
    std::vector<BasisParams> full_basis = refined.basis;
    for (auto& c : candidates) {
        full_basis.push_back(c);
    }
    int M = static_cast<int>(full_basis.size());
    std::cout << "Full (unfiltered) basis size: " << M << std::endl;

    // --- Phase 3: Build full S, H, K matrices ---
    std::cout << "\nBuilding full S, H, K matrices (" << M << "x" << M << ")..." << std::flush;

    MatrixXcd S_full = MatrixXcd::Zero(M, M);
    MatrixXcd H_full = MatrixXcd::Zero(M, M);

    for (int i = 0; i < M; i++) {
        for (int j = i; j < M; j++) {
            auto [h_ij, s_ij] = compute_HS_ij(full_basis[i], full_basis[j], perms, terms_free);
            S_full(i, j) = s_ij;
            S_full(j, i) = std::conj(s_ij);
            H_full(i, j) = h_ij;
            H_full(j, i) = std::conj(h_ij);
        }
    }

    MatrixXcd K_full = build_kick_matrix(full_basis, perms, kappa_val, k_L_val);

    std::cout << " done." << std::endl;

    // --- Phase 3.5: Diagonal preconditioning ---
    // Normalize so S has unit diagonal: S_scaled = D S D, where D_ii = 1/sqrt(S_ii)
    // This removes the huge dynamic range from different Gaussian normalizations.
    VectorXcd D_inv(M);
    for (int i = 0; i < M; i++) {
        double s_ii = S_full(i, i).real();
        if (s_ii > 1e-15) {
            D_inv(i) = 1.0 / std::sqrt(s_ii);
        } else {
            D_inv(i) = 0.0;  // degenerate basis function
        }
    }

    // Apply: S_scaled(i,j) = D_i * S(i,j) * D_j, same for H, K
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            Cd scale = D_inv(i) * D_inv(j);
            S_full(i, j) *= scale;
            H_full(i, j) *= scale;
            K_full(i, j) *= scale;
        }
    }

    // Verify diagonal is 1
    double max_diag_err = 0.0;
    for (int i = 0; i < M; i++) {
        double err = std::abs(S_full(i, i).real() - 1.0);
        if (err > max_diag_err) max_diag_err = err;
    }
    std::cout << "After preconditioning: max |S_ii - 1| = " << std::scientific
              << max_diag_err << std::defaultfloat << std::endl;

    // --- Phase 4: SVD truncation of S ---
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(S_full);
    VectorXd evals = es.eigenvalues();     // ascending order
    MatrixXcd evecs = es.eigenvectors();

    double w_max = evals(M - 1);
    double threshold = w_max * svd_rcond;

    // Count how many to keep
    int r = 0;
    for (int i = 0; i < M; i++) {
        if (evals(i) > threshold) r++;
    }

    std::cout << "\nS eigenvalue spectrum:" << std::endl;
    std::cout << "  w_max = " << std::scientific << w_max << std::endl;
    std::cout << "  w_min = " << evals(0) << std::endl;
    std::cout << "  threshold = " << threshold << std::endl;
    std::cout << "  kept = " << r << " / " << M
              << " (discarded " << (M - r) << " near-null directions)"
              << std::defaultfloat << std::endl;

    // Show eigenvalue distribution
    std::cout << "  Eigenvalue distribution (log10):" << std::endl;
    int n_neg = 0, n_tiny = 0;
    for (int i = 0; i < M; i++) {
        if (evals(i) < 0) n_neg++;
        else if (evals(i) < threshold) n_tiny++;
    }
    std::cout << "    negative: " << n_neg
              << ", below threshold: " << n_tiny
              << ", kept: " << r << std::endl;

    // Build transformation matrix X = V_r * Lambda_r^{-1/2}
    // Maps from r-dim orthonormal space to M-dim original space
    MatrixXcd X(M, r);
    int col = 0;
    for (int i = 0; i < M; i++) {
        if (evals(i) > threshold) {
            X.col(col) = evecs.col(i) / std::sqrt(evals(i));
            col++;
        }
    }

    // --- Phase 5: Transform to orthonormal subspace ---
    // H_orth = X† H X,  K_orth = X† K X
    // In this basis, S_orth = I (by construction)
    MatrixXcd H_orth = X.adjoint() * H_full * X;
    MatrixXcd K_orth = X.adjoint() * K_full * X;

    // Verify: S_orth should be identity
    MatrixXcd S_orth_check = X.adjoint() * S_full * X;
    double s_orth_err = (S_orth_check - MatrixXcd::Identity(r, r)).norm();
    std::cout << "  ||S_orth - I|| = " << std::scientific << s_orth_err
              << std::defaultfloat << std::endl;

    // --- Phase 6: Ground state in orthonormal subspace ---
    // Standard eigenvalue problem: H_orth c = E c
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es_h(H_orth);
    VectorXd energies = es_h.eigenvalues();
    MatrixXcd states = es_h.eigenvectors();

    double E_orth = energies(0);
    VectorXcd c0 = states.col(0);  // ground state in orthonormal basis

    std::cout << "\nOrthonormal subspace ground state:" << std::endl;
    std::cout << "  E = " << E_orth << " (exact = " << N * 0.5 << ")"
              << ", error = " << std::scientific << std::abs(E_orth - N * 0.5)
              << std::defaultfloat << std::endl;

    // Compute <x^2> and <p^2> for ground state
    // <T> = c† (X† H_kin X) c,  <V> = c† (X† H_harm X) c
    // We need kinetic and harmonic matrices separately
    HamiltonianTerms terms_kin_only;
    terms_kin_only.kinetic = true;
    terms_kin_only.harmonic = false;
    terms_kin_only.delta = false;
    terms_kin_only.gaussian = false;
    terms_kin_only.kicking = false;

    HamiltonianTerms terms_harm_only;
    terms_harm_only.kinetic = false;
    terms_harm_only.harmonic = true;
    terms_harm_only.delta = false;
    terms_harm_only.gaussian = false;
    terms_harm_only.kicking = false;

    MatrixXcd T_full = MatrixXcd::Zero(M, M);
    MatrixXcd V_full = MatrixXcd::Zero(M, M);
    for (int i = 0; i < M; i++) {
        for (int j = i; j < M; j++) {
            Cd scale = D_inv(i) * D_inv(j);  // same preconditioning

            auto [t_ij, s_ij_t] = compute_HS_ij(full_basis[i], full_basis[j], perms, terms_kin_only);
            (void)s_ij_t;
            T_full(i, j) = t_ij * scale;
            T_full(j, i) = std::conj(t_ij * scale);

            auto [v_ij, s_ij_v] = compute_HS_ij(full_basis[i], full_basis[j], perms, terms_harm_only);
            (void)s_ij_v;
            V_full(i, j) = v_ij * scale;
            V_full(j, i) = std::conj(v_ij * scale);
        }
    }

    MatrixXcd T_orth = X.adjoint() * T_full * X;
    MatrixXcd V_orth = X.adjoint() * V_full * X;

    auto compute_observables = [&](const VectorXcd& c) {
        // In orthonormal basis, <O> = c† O_orth c
        double E = (c.adjoint() * H_orth * c)(0).real();
        double T = (c.adjoint() * T_orth * c)(0).real();
        double V = (c.adjoint() * V_orth * c)(0).real();
        double p2 = 2.0 * mass * T;
        double x2 = (mass * omega * omega > 0.0) ? (2.0 * V / (mass * omega * omega)) : 0.0;
        return std::make_tuple(E, x2, p2);
    };

    {
        auto [E0, x2_0, p2_0] = compute_observables(c0);
        std::cout << "  <x^2> = " << x2_0 << " (exact = 1)" << std::endl;
        std::cout << "  <p^2> = " << p2_0 << " (exact = 1)" << std::endl;
    }

    // --- Phase 7: Kicked evolution in orthonormal subspace ---
    // In orthonormal basis:
    //   Kick: c_new = K_orth * c_old  (since S_orth = I)
    //   Free evolution: c(t) = U exp(-i D t) U† c(0) where H_orth = U D U†
    //
    // Free evolution uses the eigenbasis of H_orth (already computed above).

    // Get exact reference
    auto exact = kicked_exact_1particle(1.0, 1.0, k_L_val, kappa_val, T_period, n_kicks);

    // First-kick instantaneous diagnostic
    {
        VectorXcd c_kicked = K_orth * c0;
        double norm_before = c0.squaredNorm();  // should be 1 in orthonormal basis
        double norm_after = c_kicked.squaredNorm();
        double fidelity = norm_after / norm_before;

        // Renormalize
        c_kicked /= std::sqrt(norm_after / norm_before);

        auto [E_kick, x2_kick, p2_kick] = compute_observables(c_kicked);
        double exact_E1 = 2.0 * exact.kick_energies[1];
        double exact_x2_1 = static_cast<double>(N) * 0.5;
        double exact_V1 = 0.5 * mass * omega * omega * exact_x2_1;
        double exact_p2_1 = 2.0 * mass * (exact_E1 - exact_V1);

        std::cout << "\n--- First Kick Instantaneous Diagnostic ---" << std::endl;
        std::cout << "Projection fidelity = " << fidelity
                  << "  (|1-fid|=" << std::scientific << std::abs(1.0 - fidelity)
                  << std::defaultfloat << ")" << std::endl;
        std::cout << "E after kick:     orth=" << E_kick
                  << ", exact=" << exact_E1
                  << ", dE=" << std::scientific << (E_kick - exact_E1)
                  << std::defaultfloat << std::endl;
        std::cout << "<x^2> after kick: orth=" << x2_kick
                  << ", exact=" << exact_x2_1
                  << ", dx2=" << std::scientific << (x2_kick - exact_x2_1)
                  << std::defaultfloat << std::endl;
        std::cout << "<p^2> after kick: orth=" << p2_kick
                  << ", exact=" << exact_p2_1
                  << ", dp2=" << std::scientific << (p2_kick - exact_p2_1)
                  << std::defaultfloat << std::endl;
    }

    // --- Full kicked evolution ---
    std::cout << "\n--- Kicked Evolution (orthonormal subspace, r=" << r << ") ---" << std::endl;

    VectorXcd c = c0;

    // Pre-compute free evolution propagator
    // H_orth eigenbasis already in es_h
    // Free propagator for one period:
    //   c(T) = states * diag(exp(-i E_n T)) * states† * c(0)
    VectorXcd phases(r);
    for (int i = 0; i < r; i++) {
        phases(i) = std::exp(Cd(0.0, -energies(i) * T_period));
    }

    // Print header
    {
        double E_ecg = (c.adjoint() * H_orth * c)(0).real();
        double E_ex = 2.0 * exact.kick_energies[0];
        double rel = std::abs(E_ecg - E_ex) / std::abs(E_ex);
        std::cout << "Kick  0: E_orth=" << std::setprecision(8) << E_ecg
                  << ", E_exact=" << E_ex
                  << ", rel_error=" << std::scientific << rel
                  << std::defaultfloat << std::endl;
    }

    double max_rel_err = 0.0;
    std::vector<double> fidelities;

    for (int n = 0; n < n_kicks; n++) {
        // 1. Apply kick: c_new = K_orth * c_old
        VectorXcd c_kicked = K_orth * c;
        double norm_before = c.squaredNorm();
        double norm_after = c_kicked.squaredNorm();
        double fid = norm_after / norm_before;
        fidelities.push_back(fid);

        // Renormalize to preserve norm
        c_kicked *= std::sqrt(norm_before / norm_after);

        // 2. Free evolution for one period
        // c(T) = states * diag(phases) * states† * c_kicked
        VectorXcd coeffs = states.adjoint() * c_kicked;
        for (int i = 0; i < r; i++) {
            coeffs(i) *= phases(i);
        }
        c = states * coeffs;

        // 3. Measure
        double E_ecg = (c.adjoint() * H_orth * c)(0).real();
        double E_ex = 2.0 * exact.kick_energies[n + 1];
        double rel = std::abs(E_ecg - E_ex) / std::abs(E_ex);
        if (rel > max_rel_err) max_rel_err = rel;

        const char* status = (std::abs(fid - 1.0) < 0.05) ? "PASS" : "FAIL";
        std::cout << "Kick " << std::setw(2) << (n + 1)
                  << ": fidelity=" << std::setprecision(6) << fid
                  << " (" << status << ")"
                  << ", E_orth=" << std::setprecision(8) << E_ecg
                  << ", E_exact=" << E_ex
                  << ", rel_error=" << std::scientific << rel
                  << std::defaultfloat << std::endl;
    }

    // --- Summary ---
    double fid_min = *std::min_element(fidelities.begin(), fidelities.end());
    double fid_max = *std::max_element(fidelities.begin(), fidelities.end());
    double fid_mean = 0.0;
    for (double f : fidelities) fid_mean += f;
    fid_mean /= fidelities.size();

    std::cout << "\n--- Summary ---" << std::endl;
    std::cout << "Subspace dimension: " << r << " / " << M
              << " (from " << static_cast<int>(refined.basis.size()) << " gs + "
              << candidates.size() << " candidates)" << std::endl;
    std::cout << "SVD rcond: " << std::scientific << svd_rcond << std::defaultfloat << std::endl;
    std::cout << "Kick fidelity: min=" << std::setprecision(6) << fid_min
              << ", max=" << fid_max
              << ", mean=" << fid_mean << std::endl;
    std::cout << "Max relative energy error: " << std::scientific << max_rel_err
              << std::defaultfloat << std::endl;

    // Compare with greedy baseline
    std::cout << "\n--- Comparison with greedy baseline ---" << std::endl;
    std::cout << "Greedy (14 basis, cond=1e10): kick 50 rel_error ~ 869%" << std::endl;
    std::cout << "Orthogonal (r=" << r << "): max rel_error = "
              << std::scientific << max_rel_err << std::defaultfloat << std::endl;

    return 0;
}
