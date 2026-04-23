#include "step1_imag_time_ecg.hpp"
#include "common.hpp"
#include "ecg_wavefunction.hpp"
#include "permutation.hpp"
#include "physical_constants.hpp"
#include "svm.hpp"
#include "observables.hpp"
#include "hamiltonian.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

namespace ecg1d {
namespace verify {

static HamiltonianTerms make_terms_full() {
    HamiltonianTerms t;
    t.kinetic = true;
    t.harmonic = true;
    t.delta = false;
    t.gaussian = false;
    t.kicking = true;
    t.kicking_scale = 1.0;
    return t;
}

static std::vector<AlphaIndex> make_alpha_list(int N, int K) {
    std::vector<AlphaIndex> a;
    for (int i = 0; i < K; i++) a.push_back({1, i, 0, 0});                    // u
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) a.push_back({2, i, j, 0});                // B diagonal
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) a.push_back({3, i, j, 0});                // R
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
            for (int k = j; k < N; k++) a.push_back({4, i, j, k});            // A symmetric
    return a;
}

Step1Result step1_imag_time_ecg(int N, int K,
                                const Eigen::VectorXd& x_grid,
                                const Eigen::VectorXd& k_grid) {
    std::cout << "\n=== Step 1 — ECG imaginary-time ground state (N=" << N
              << ", K=" << K << ") ===\n";

    HamiltonianTerms terms = make_terms_full();
    PermutationSet perms = PermutationSet::generate(N);

    // E_floor: V(x) = 1/2 x^2 + kappa*cos(2 k_L x) has min V = V(0) = 1 for default
    // parameters, so E > 1 per particle. Use N*1.0 to reject numerical ghosts.
    double E_floor = static_cast<double>(N) * 1.0;
    double E_svm = 0.0;
    std::vector<BasisParams> basis;

    if (N == 1) {
        // svm_build_basis / stochastic_refine both use `random_basis_2particle`
        // internally, which is hardcoded for N=2 and produces invalid N=1
        // trial basis members (ghost states with indefinite S). For N=1 we
        // therefore seed a small basis manually (HO ground state + random
        // perturbations via BasisParams::randoman) and let TDVP do the work.
        std::cout << "[step1] N=1 init: HO ground state + "
                  << (K - 1) << " randoman perturbations" << std::endl;
        MatrixXcd A0 = MatrixXcd::Zero(1, 1); A0(0, 0) = Cd(0.5, 0.0);
        MatrixXcd B0 = MatrixXcd::Zero(1, 1);
        VectorXcd R0 = VectorXcd::Zero(1);
        basis.push_back(BasisParams::from_arrays(Cd(1.0, 0.0), A0, B0, R0, 0));
        for (int k = 1; k < K; k++) {
            basis.push_back(BasisParams::randoman(N, /*seed=*/42 + k,
                                                  /*name=*/k));
        }
        auto [H0, S0] = build_HS(basis, perms, terms);
        E_svm = lowest_energy(H0, S0, 1e8, E_floor);
        set_u_from_eigenvector(basis, H0, S0);
        std::cout << "[step1] init E = " << std::setprecision(12) << E_svm << std::endl;
    } else {
        // N>=2: follow main.cpp run_svm_tdvp pattern (SVM + refine).
        std::cout << "[step1] Phase 1: SVM greedy build (K_max=" << K
                  << ", n_trials=5000, E_floor=" << E_floor << ")..." << std::endl;
        SvmResult svm = svm_build_basis(N, /*K_max=*/K, /*n_trials=*/5000,
                                        terms, /*seed=*/42, /*E_lower_bound=*/E_floor);
        E_svm = lowest_energy(svm.H, svm.S);
        std::cout << "[step1] SVM E = " << std::setprecision(12) << E_svm << std::endl;

        std::cout << "[step1] Phase 1.5: stochastic refinement..." << std::endl;
        SvmResult refined = stochastic_refine(svm.basis, svm.H, svm.S,
                                              perms, N, terms,
                                              /*n_trials=*/500, /*max_rounds=*/10,
                                              /*seed=*/123, E_floor);
        double E_refine = lowest_energy(refined.H, refined.S);
        std::cout << "[step1] refined E = " << std::setprecision(12) << E_refine << std::endl;

        basis = refined.basis;
        set_u_from_eigenvector(basis, refined.H, refined.S);
    }

    // 2. TDVP imaginary-time refinement
    std::vector<AlphaIndex> alpha = make_alpha_list(N, K);

    SolverConfig cfg;
    cfg.lambda_C          = 1e-8;
    cfg.rcond             = 1e-4;
    cfg.resolve_u         = true;
    cfg.dtao_grow         = 1.5;
    cfg.dtao_max          = 10.0;
    cfg.energy_tol        = 1e-12;
    cfg.stagnation_window = 50;
    cfg.adaptive_lambda   = true;
    cfg.lambda_max        = 1e-4;
    cfg.optimize_A        = true;
    cfg.optimize_B        = true;
    cfg.optimize_R        = true;

    std::cout << "[step1] Phase 2: TDVP imaginary-time polish (max 300 steps)..."
              << std::endl;
    evolution(alpha, basis, /*dtao=*/1.0, /*max_steps=*/300,
              /*tol=*/1e-10, terms, cfg, &perms);

    // 3. Final energy via generalized eigenvalue problem (optimal u)
    auto [H_final, S_final] = build_HS(basis, perms, terms);
    double E_final = lowest_energy(H_final, S_final);
    set_u_from_eigenvector(basis, H_final, S_final);

    EnergyComponents ec = compute_energy_components(basis, terms);

    std::cout << "[step1] Final energy: " << std::setprecision(12) << E_final << std::endl;
    std::cout << "        components:  T=" << ec.kinetic
              << "  V_ho=" << ec.harmonic
              << "  V_kick=" << ec.kicking << std::endl;

    // 4. Write CSVs
    std::ostringstream suffix;
    suffix << "_N" << N << "_K" << K << ".csv";

    {
        std::ofstream f(out_path("step1_ecg_energy" + suffix.str()));
        f << "quantity,value\n";
        f << std::setprecision(15);
        f << "N," << N << std::endl;
        f << "K," << K << std::endl;
        f << "E_ECG," << E_final << std::endl;
        f << "E_SVM_init," << E_svm << std::endl;
        f << "E_kinetic," << ec.kinetic << std::endl;
        f << "E_harmonic," << ec.harmonic << std::endl;
        f << "E_kicking," << ec.kicking << std::endl;
    }

    save_basis_csv(out_path("step1_ecg_basis" + suffix.str()), basis);

    Eigen::VectorXd n_x = ecg_single_particle_density(basis, perms, x_grid, 0);
    write_two_column(out_path("step1_ecg_density" + suffix.str()),
                     "x,n_x", x_grid, n_x);

    Eigen::VectorXd n_k = momentum_distribution(basis, perms, k_grid, 0);
    write_two_column(out_path("step1_ecg_nk" + suffix.str()),
                     "k,n_k", k_grid, n_k);

    // N=1: write the full complex ECG wavefunction on the x grid so step3 /
    // plotting can compare Re(psi), Im(psi), and compute |<psi_ECG|psi_ref>|^2.
    if (N == 1) {
        Eigen::VectorXcd psi = ecg_wavefunction_1p(basis, x_grid);
        std::ofstream fpsi(out_path("step1_ecg_psi" + suffix.str()));
        fpsi << "x,Re_psi,Im_psi\n";
        fpsi << std::setprecision(15);
        for (int j = 0; j < x_grid.size(); j++) {
            fpsi << x_grid(j) << "," << psi(j).real() << "," << psi(j).imag() << "\n";
        }
    }

    Step1Result r;
    r.basis      = basis;
    r.E_ecg      = E_final;
    r.E_kinetic  = ec.kinetic;
    r.E_harmonic = ec.harmonic;
    r.E_kicking  = ec.kicking;
    r.n_steps    = 0;
    return r;
}

} // namespace verify
} // namespace ecg1d
