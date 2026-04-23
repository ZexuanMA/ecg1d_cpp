#include "realtime_tdvp.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_gradient.hpp"
#include "physical_constants.hpp"
#include "kick_operator.hpp"  // free_evolve_fixed_basis, build_HS-like S
#include "svm.hpp"            // build_HS: (H, S) on arbitrary basis
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <iomanip>

namespace ecg1d {

// Compute dz = -i * C_bar^{-1} * g at the current basis, with SVD-based
// pseudoinverse and (optional) Tikhonov regularization.
// Returns dz (size = alpha_z_list.size()) and cond_C / effective_rank.
static VectorXcd compute_rhs_dz(const std::vector<AlphaIndex>& alpha_z_list,
                                const std::vector<BasisParams>& basis,
                                const HamiltonianTerms& terms,
                                const SolverConfig& config,
                                double* out_cond_C = nullptr,
                                int* out_rank = nullptr,
                                double out_sv_small[3] = nullptr) {
    const int d = static_cast<int>(alpha_z_list.size());
    MatrixXcd C_mat = assemble_C(alpha_z_list, basis);
    MatrixXcd C_bar = C_mat.conjugate().transpose();
    if (config.lambda_C > 0.0) {
        C_bar += config.lambda_C * MatrixXcd::Identity(d, d);
    }
    VectorXcd g = assemble_grad(alpha_z_list, /*Real=*/false, basis, terms);
    VectorXcd rhs = Cd(0.0, -1.0) * g;   // real-time: i C z_dot = g  ->  dz = -i C^{-1} g

    Eigen::JacobiSVD<MatrixXcd> svd(C_bar, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd sv = svd.singularValues();
    double thr = config.rcond * sv(0);
    if (out_cond_C) *out_cond_C = sv(0) / std::max(sv(sv.size() - 1), 1e-300);
    if (out_rank) {
        int r = 0;
        for (int i = 0; i < sv.size(); ++i) if (sv(i) > thr) r++;
        *out_rank = r;
    }
    // sv is returned from JacobiSVD in DESCENDING order, so last entries are
    // the smallest. Pull the bottom 3 (or fewer if d < 3).
    if (out_sv_small) {
        int n = sv.size();
        for (int k = 0; k < 3; ++k) {
            int idx = n - 1 - k;
            out_sv_small[k] = (idx >= 0) ? sv(idx)
                                         : std::numeric_limits<double>::quiet_NaN();
        }
    }

    VectorXcd Utrhs = svd.matrixU().adjoint() * rhs;
    VectorXcd sinv(sv.size());
    if (config.wiener_smooth) {
        // Wiener-style smooth pseudoinverse: sinv = sv / (sv^2 + lam^2).
        // lam = rcond * sv_max (reuse rcond as the soft-shoulder scale).
        // Continuous in sv — no discontinuous truncation when sv drifts across
        // the "rcond threshold". Amplification peaks at 1/(2*lam) when sv=lam.
        double lam = thr;   // thr is already rcond * sv(0)
        double lam2 = lam * lam;
        for (int i = 0; i < sv.size(); ++i) {
            double s = sv(i);
            sinv(i) = Utrhs(i) * s / (s * s + lam2);
        }
    } else {
        // Hard rcond truncation (original behavior).
        for (int i = 0; i < sv.size(); ++i) {
            sinv(i) = (sv(i) > thr) ? Utrhs(i) / sv(i) : Cd(0.0, 0.0);
        }
    }
    return svd.matrixV() * sinv;
}

RealtimeStepResult realtime_tdvp_step(const std::vector<AlphaIndex>& alpha_z_list,
                                      const std::vector<BasisParams>& basis,
                                      double dt,
                                      const HamiltonianTerms& terms,
                                      const SolverConfig& config,
                                      RtIntegrator integrator) {
    double cond_C = 0.0;
    int rank = 0;
    double sv_small[3] = {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};

    VectorXcd dz_combined;
    if (integrator == RtIntegrator::Euler) {
        dz_combined = compute_rhs_dz(alpha_z_list, basis, terms, config,
                                     &cond_C, &rank, sv_small);
    } else {
        // Classical RK4 on z: z' = f(z) = -i C^{-1} g
        VectorXcd k1 = compute_rhs_dz(alpha_z_list, basis, terms, config,
                                      &cond_C, &rank, sv_small);

        std::vector<BasisParams> b1 = basis;
        update_basis_function(b1, k1, 0.5 * dt, alpha_z_list, 0);
        VectorXcd k2 = compute_rhs_dz(alpha_z_list, b1, terms, config);

        std::vector<BasisParams> b2 = basis;
        update_basis_function(b2, k2, 0.5 * dt, alpha_z_list, 0);
        VectorXcd k3 = compute_rhs_dz(alpha_z_list, b2, terms, config);

        std::vector<BasisParams> b3 = basis;
        update_basis_function(b3, k3, dt, alpha_z_list, 0);
        VectorXcd k4 = compute_rhs_dz(alpha_z_list, b3, terms, config);

        dz_combined = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }

    std::vector<BasisParams> basis_new = basis;
    update_basis_function(basis_new, dz_combined, dt, alpha_z_list, /*updata_constant=*/0);

    RealtimeStepResult result;
    result.basis = std::move(basis_new);
    result.used_dt = dt;
    result.cond_C = cond_C;
    result.effective_rank = rank;
    result.sv_small[0] = sv_small[0];
    result.sv_small[1] = sv_small[1];
    result.sv_small[2] = sv_small[2];
    result.dz_norm = dz_combined.norm();
    return result;
}

// -----------------------------------------------------------------------------
// Main loop: evolve from t=0 to T_total sampling observables into a trace
// -----------------------------------------------------------------------------

static void sample_observables(const std::vector<BasisParams>& basis,
                               double t,
                               const HamiltonianTerms& terms,
                               const std::vector<BasisParams>& basis0,
                               Cd norm0,
                               RealtimeTrace& trace) {
    // <psi|psi>
    Cd S = overlap(basis);
    // Energies (kinetic + harmonic by default; other terms if enabled)
    Cd T_kin = kinetic_energy_functional(basis);
    Cd V_ho  = Harmonic_functional(basis);
    Cd E = (T_kin + V_ho) / S;   // normalized

    double x2 = 2.0 * (V_ho / S).real() / (mass * omega * omega);
    double p2 = 2.0 * mass * (T_kin / S).real();

    // <psi(0)|psi(t)>  — a bit ugly because overlap() wants one basis; build 2K basis
    // Simpler: manually compute via PairCache.
    const int Kt = static_cast<int>(basis.size());
    const int K0 = static_cast<int>(basis0.size());
    const int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);
    Cd amp(0.0, 0.0);
    for (int i = 0; i < K0; ++i) {
        Cd conj_ui = std::conj(basis0[i].u);
        for (int j = 0; j < Kt; ++j) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; ++p) {
                PairCache c = PairCache::build(basis0[i], basis[j], perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G;
            }
            amp += conj_ui * basis[j].u * sum_p;
        }
    }

    trace.t.push_back(t);
    trace.E.push_back(E.real());
    trace.norm.push_back(S.real());
    trace.x2.push_back(x2);
    trace.p2.push_back(p2);
    trace.overlap0.push_back(amp);
    (void)norm0;
}

RealtimeEvolutionResult realtime_tdvp_evolution(
    const std::vector<AlphaIndex>& alpha_z_list,
    std::vector<BasisParams> basis_init,
    double T_total,
    const HamiltonianTerms& terms,
    const SolverConfig& solver_cfg,
    const RealtimeEvolutionConfig& rt_cfg) {

    RealtimeEvolutionResult result;
    std::vector<BasisParams> basis = std::move(basis_init);
    const std::vector<BasisParams> basis0 = basis;   // t=0 snapshot for fidelity
    Cd norm0 = overlap(basis0);

    sample_observables(basis, 0.0, terms, basis0, norm0, result.trace);

    const int n_steps = static_cast<int>(std::round(T_total / rt_cfg.dt));
    // Blowup guard on cond(C_bar).
    // - Hard pinv mode:   cond is bounded by Tikhonov floor (lambda_C),
    //                     so a genuine blowup (>1e14) means basis has drifted
    //                     into pathological region.
    // - Wiener mode:      Wiener filter handles sv=0 natively (sinv=0), so
    //                     cond(C) = sv_max / numerical_zero ~ 1e16 is the
    //                     EXPECTED value for healthy state (true gauge nullspace
    //                     uncovered). Effectively disable the guard there.
    const double cond_abort = solver_cfg.wiener_smooth ? 1e30 : 1e14;
    double t = 0.0;
    // For u-split mode we need permutations to build (H, S) each step
    const int N = basis0[0].N();
    PermutationSet perms_step = PermutationSet::generate(N);

    // Fixed initial norm — used as post-correction target if enforce_norm is on
    const double norm_init = overlap(basis0).real();

    for (int step = 0; step < n_steps; ++step) {
        RealtimeStepResult r = realtime_tdvp_step(alpha_z_list, basis, rt_cfg.dt,
                                                  terms, solver_cfg, rt_cfg.integrator);
        basis = std::move(r.basis);

        // Lie-Trotter splitting for u: after TDVP moved {A,R} (and whatever
        // else was in alpha_z_list), evolve u exactly on the updated basis.
        if (rt_cfg.u_split_trotter) {
            auto [H_mat, S_mat] = build_HS(basis, perms_step, terms);
            free_evolve_fixed_basis(basis, H_mat, S_mat, rt_cfg.dt);
        }

        // Post-correct norm: rescale u so <psi|psi> == norm_init. Pragmatic
        // unitarity enforcement; does nothing for E conservation.
        if (rt_cfg.enforce_norm) {
            double norm_now = overlap(basis).real();
            if (norm_now > 1e-15 && std::isfinite(norm_now)) {
                double scale = std::sqrt(norm_init / norm_now);
                for (auto& bp : basis) bp.u *= scale;
            }
        }

        t += rt_cfg.dt;

        if ((step + 1) % rt_cfg.sample_every == 0) {
            sample_observables(basis, t, terms, basis0, norm0, result.trace);
        }

        // Per-step diagnostic print (coarsened by print_every to avoid flooding)
        if (rt_cfg.verbose && ((step + 1) % rt_cfg.print_every == 0)) {
            Cd S = overlap(basis);
            Cd T_kin = kinetic_energy_functional(basis);
            Cd V_ho  = Harmonic_functional(basis);
            double E = ((T_kin + V_ho) / S).real();

            // Per-basis min Re(A+B) — if this approaches 0, the Gaussian is
            // losing normalizability (pathology suspected for step-348 blowup).
            double min_rab = 1e99;
            int    argmin_i = -1;
            double min_reA = 1e99;
            int    argmin_A = -1;
            for (size_t i = 0; i < basis.size(); ++i) {
                double rab = (basis[i].A(0, 0) + basis[i].B(0, 0)).real();
                if (rab < min_rab) { min_rab = rab; argmin_i = (int)i; }
                double reA = basis[i].A(0, 0).real();
                if (reA < min_reA) { min_reA = reA; argmin_A = (int)i; }
            }

            std::cout << "[rt_tdvp] step " << (step + 1) << "/" << n_steps
                      << " t=" << std::fixed << std::setprecision(4) << t
                      << " E=" << std::setprecision(8) << E
                      << " norm=" << S.real()
                      << " cond=" << std::scientific << std::setprecision(2)
                      << r.cond_C
                      << " sv_min3=(" << r.sv_small[2] << "," << r.sv_small[1]
                      << "," << r.sv_small[0] << ")"
                      << " |dz|=" << std::setprecision(2) << r.dz_norm
                      << " min(Re A+B)=" << std::fixed << std::setprecision(4) << min_rab
                      << "@" << argmin_i
                      << " min(Re A)=" << min_reA << "@" << argmin_A
                      << std::defaultfloat << "\n";
        }

        // Blowup guard: cond(C) >> 1 means the variational metric has a near-zero
        // singular value, i.e. two parameters have become degenerate. Dump basis
        // so we can see which one drifted, then stop.
        if (r.cond_C > cond_abort || !std::isfinite(r.cond_C)) {
            std::cout << "\n[rt_tdvp] BLOWUP at step " << (step + 1)
                      << " t=" << std::fixed << std::setprecision(4) << t
                      << " cond(C)=" << std::scientific << r.cond_C
                      << std::defaultfloat << " (threshold " << cond_abort << ")\n";
            std::cout << "  basis dump:\n";
            for (size_t i = 0; i < basis.size(); ++i) {
                std::cout << "    [" << i << "] u=" << basis[i].u
                          << "  A=" << basis[i].A(0, 0)
                          << "  B=" << basis[i].B(0, 0)
                          << "  R=" << basis[i].R(0)
                          << "  A+B=" << (basis[i].A(0, 0) + basis[i].B(0, 0))
                          << "\n";
            }
            result.basis_final = std::move(basis);
            result.n_steps = step + 1;
            return result;
        }
    }

    // Final sample (in case we don't land on sample_every)
    if (result.trace.t.empty() || result.trace.t.back() < t) {
        sample_observables(basis, t, terms, basis0, norm0, result.trace);
    }

    result.basis_final = std::move(basis);
    result.n_steps = n_steps;
    return result;
}

} // namespace ecg1d
