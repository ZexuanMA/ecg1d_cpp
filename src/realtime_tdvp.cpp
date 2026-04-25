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
#include <limits>

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
// Dormand-Prince 5(4) adaptive step (used when integrator == RK45)
// 7 stages (last stage k7 = f at 5th-order result = FSAL, but we don't cache).
// Returns proposed dt_next and accepted flag.
// -----------------------------------------------------------------------------

struct DP5StepResult {
    std::vector<BasisParams> basis_new;
    double used_dt;
    double dt_next;
    bool   accepted;
    double err_norm;
    double cond_C;
    int    effective_rank;
    double sv_small[3];
    double dz_norm;
};

static DP5StepResult realtime_tdvp_step_dp5(
    const std::vector<AlphaIndex>& alpha_z_list,
    const std::vector<BasisParams>& basis,
    double dt,
    const HamiltonianTerms& terms,
    const SolverConfig& config,
    double abs_tol, double rel_tol) {

    // Dormand-Prince 5(4) Butcher tableau
    const double a21 = 1.0/5.0;
    const double a31 = 3.0/40.0,     a32 = 9.0/40.0;
    const double a41 = 44.0/45.0,    a42 = -56.0/15.0,  a43 = 32.0/9.0;
    const double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0,
                 a53 = 64448.0/6561.0, a54 = -212.0/729.0;
    const double a61 = 9017.0/3168.0, a62 = -355.0/33.0,
                 a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
    // 5th-order weights (also = last-row a7i, FSAL)
    const double w1 = 35.0/384.0, w3 = 500.0/1113.0, w4 = 125.0/192.0,
                 w5 = -2187.0/6784.0, w6 = 11.0/84.0;
    // Error coefficients e_i = b5_i - b4_i
    const double e1 = 71.0/57600.0,  e3 = -71.0/16695.0, e4 = 71.0/1920.0,
                 e5 = -17253.0/339200.0, e6 = 22.0/525.0, e7 = -1.0/40.0;

    double cond_C = 0.0;
    int rank = 0;
    double sv_small[3] = {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};

    // k1
    VectorXcd k1 = compute_rhs_dz(alpha_z_list, basis, terms, config,
                                   &cond_C, &rank, sv_small);

    auto scaled_step = [&](const VectorXcd& dz_combined, double scale) {
        std::vector<BasisParams> b_out = basis;
        update_basis_function(b_out, dz_combined, scale, alpha_z_list, 0);
        return b_out;
    };

    // k2 = f(basis + dt*a21*k1)
    std::vector<BasisParams> bs2 = scaled_step(k1, a21 * dt);
    VectorXcd k2 = compute_rhs_dz(alpha_z_list, bs2, terms, config);

    // k3 = f(basis + dt*(a31 k1 + a32 k2))
    std::vector<BasisParams> bs3 = scaled_step(a31 * k1 + a32 * k2, dt);
    VectorXcd k3 = compute_rhs_dz(alpha_z_list, bs3, terms, config);

    // k4
    std::vector<BasisParams> bs4 = scaled_step(a41 * k1 + a42 * k2 + a43 * k3, dt);
    VectorXcd k4 = compute_rhs_dz(alpha_z_list, bs4, terms, config);

    // k5
    std::vector<BasisParams> bs5 = scaled_step(
        a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4, dt);
    VectorXcd k5 = compute_rhs_dz(alpha_z_list, bs5, terms, config);

    // k6
    std::vector<BasisParams> bs6 = scaled_step(
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5, dt);
    VectorXcd k6 = compute_rhs_dz(alpha_z_list, bs6, terms, config);

    // 5th-order combined dz
    VectorXcd dz_5 = w1 * k1 + w3 * k3 + w4 * k4 + w5 * k5 + w6 * k6;

    // k7 = f(basis + dt*dz_5) (FSAL endpoint, used for error estimate)
    std::vector<BasisParams> bs7 = scaled_step(dz_5, dt);
    VectorXcd k7 = compute_rhs_dz(alpha_z_list, bs7, terms, config);

    // Error estimate: err_vec = e1 k1 + e3 k3 + e4 k4 + e5 k5 + e6 k6 + e7 k7
    VectorXcd err_vec = e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7;

    // Scaled RMS error norm
    const int n_alpha = static_cast<int>(err_vec.size());
    double sum_sq = 0.0;
    for (int i = 0; i < n_alpha; ++i) {
        double y_mag = std::abs(dz_5(i)) * dt;
        double scale = abs_tol + rel_tol * y_mag;
        double e_i = std::abs(err_vec(i)) * dt;
        double r = e_i / scale;
        sum_sq += r * r;
    }
    double err_norm = (n_alpha > 0) ? std::sqrt(sum_sq / n_alpha) : 0.0;

    // Update basis with 5th-order dz
    std::vector<BasisParams> basis_new = basis;
    update_basis_function(basis_new, dz_5, dt, alpha_z_list, 0);

    DP5StepResult result;
    result.basis_new = std::move(basis_new);
    result.used_dt = dt;
    result.err_norm = err_norm;
    result.cond_C = cond_C;
    result.effective_rank = rank;
    for (int i = 0; i < 3; ++i) result.sv_small[i] = sv_small[i];
    result.dz_norm = dz_5.norm();

    // Check basis health: reject step if the 5th-order solution produced
    // NaN/Inf OR absurdly-large finite values (> 1e15). The latter happens
    // when the adaptive controller overshoots and the basis drifts into
    // numerically unphysical regions (e.g. Re(A+B) = 1e96 = delta function).
    bool basis_healthy = true;
    const double mag_abort = 1e15;
    auto unhealthy = [&](Cd v) {
        return !std::isfinite(v.real()) || !std::isfinite(v.imag())
             || std::abs(v.real()) > mag_abort || std::abs(v.imag()) > mag_abort;
    };
    for (const auto& bp : result.basis_new) {
        if (unhealthy(bp.u)) { basis_healthy = false; break; }
        for (int a = 0; a < bp.N(); ++a) {
            if (unhealthy(bp.A(a, a) + bp.B(a, a))) { basis_healthy = false; break; }
            if (unhealthy(bp.R(a)))                 { basis_healthy = false; break; }
        }
        if (!basis_healthy) break;
    }

    // Accept only if (err_norm <= 1) AND basis is numerically healthy.
    // NaN err with healthy basis (shouldn't happen but safe): reject + shrink.
    if (!basis_healthy || !std::isfinite(err_norm)) {
        result.accepted = false;
    } else {
        result.accepted = (err_norm <= 1.0);
    }

    // Step-size control: dt_next = dt * 0.9 * err_norm^{-1/5}, clamped to [0.2, 5]
    double fac;
    if (!std::isfinite(err_norm) || !basis_healthy) {
        fac = 0.2;   // aggressive shrink when unhealthy
    } else if (err_norm > 1e-30) {
        fac = 0.9 * std::pow(err_norm, -0.2);
    } else {
        fac = 5.0;
    }
    fac = std::max(0.2, std::min(5.0, fac));
    result.dt_next = dt * fac;

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

    // <x>, <p> from analytic ECG matrix elements. Single-particle (N=1) only;
    // for N>=2 we record NaN (the marginal <x_a> would require a kernel that
    // sums over particle index a with permutation symmetrization — added later
    // if Step-4 is extended to N>=2).
    //   Per pair (i,j), single perm:
    //     <phi_i| x |phi_j>      = M_G * mu                    (mu = c.mu(0))
    //     <phi_i| -i d/dx |phi_j> = -i * M_G * ((-2 alpha_j) mu + beta_j)
    //   where alpha_j = c.K_Mj(0,0) = A_j+B_j and beta_j = c.g_Mj(0) = 2 R_j B_j.
    double x_mean = std::numeric_limits<double>::quiet_NaN();
    double p_mean = std::numeric_limits<double>::quiet_NaN();
    if (N == 1) {
        Cd amp_x(0.0, 0.0);
        Cd amp_p(0.0, 0.0);
        for (int i = 0; i < Kt; ++i) {
            Cd conj_ui = std::conj(basis[i].u);
            for (int j = 0; j < Kt; ++j) {
                Cd term_x(0.0, 0.0);
                Cd term_p(0.0, 0.0);
                for (int p = 0; p < perms.SN; ++p) {
                    PairCache c = PairCache::build(basis[i], basis[j], perms.matrices[p]);
                    double sign = static_cast<double>(perms.signs[p]);
                    Cd mu     = c.mu(0);
                    Cd alpha  = c.K_Mj(0, 0);
                    Cd beta   = c.g_Mj(0);
                    term_x += sign * c.M_G * mu;
                    term_p += sign * Cd(0.0, -1.0) * c.M_G * ((-2.0 * alpha) * mu + beta);
                }
                amp_x += conj_ui * basis[j].u * term_x;
                amp_p += conj_ui * basis[j].u * term_p;
            }
        }
        x_mean = (amp_x / S).real();
        p_mean = (amp_p / S).real();
    }

    trace.t.push_back(t);
    trace.E.push_back(E.real());
    trace.norm.push_back(S.real());
    trace.x_mean.push_back(x_mean);
    trace.p_mean.push_back(p_mean);
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

    // For fixed-step Euler/RK4 use n_steps; for RK45 adaptive, the loop is
    // controlled by (t < T_total) and dt adjusts itself each iteration.
    const bool adaptive = (rt_cfg.integrator == RtIntegrator::RK45);
    const int n_steps_fixed = static_cast<int>(std::round(T_total / rt_cfg.dt));
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

    // Current dt (for adaptive: updated each step; for fixed: never changes)
    double dt_current = rt_cfg.dt;
    int step = 0;
    int rejected_count = 0;
    while (true) {
        if (!adaptive && step >= n_steps_fixed) break;
        if (adaptive && t >= T_total - 1e-15) break;

        // Cap dt so we don't overshoot T_total (adaptive only)
        double dt_try = adaptive
            ? std::min(dt_current, T_total - t)
            : rt_cfg.dt;

        RealtimeStepResult r;
        if (rt_cfg.integrator == RtIntegrator::RK45) {
            // Adaptive Dormand-Prince loop: retry with smaller dt on rejection
            int retry = 0;
            const int max_retry = 20;
            while (true) {
                dt_try = std::min(dt_try, T_total - t);
                dt_try = std::max(dt_try, rt_cfg.rk45_dt_min);
                dt_try = std::min(dt_try, rt_cfg.rk45_dt_max);
                DP5StepResult dp = realtime_tdvp_step_dp5(
                    alpha_z_list, basis, dt_try, terms, solver_cfg,
                    rt_cfg.rk45_abs_tol, rt_cfg.rk45_rel_tol);
                if (dp.accepted || retry >= max_retry || dt_try <= rt_cfg.rk45_dt_min * 1.001) {
                    r.basis = std::move(dp.basis_new);
                    r.used_dt = dp.used_dt;
                    r.cond_C = dp.cond_C;
                    r.effective_rank = dp.effective_rank;
                    for (int i = 0; i < 3; ++i) r.sv_small[i] = dp.sv_small[i];
                    r.dz_norm = dp.dz_norm;
                    dt_current = std::max(rt_cfg.rk45_dt_min,
                                          std::min(rt_cfg.rk45_dt_max, dp.dt_next));
                    if (!dp.accepted) rejected_count++;
                    break;
                }
                // Reject: shrink dt_try and retry
                rejected_count++;
                dt_try = std::max(rt_cfg.rk45_dt_min, dp.dt_next);
                retry++;
            }
        } else {
            r = realtime_tdvp_step(alpha_z_list, basis, dt_try,
                                   terms, solver_cfg, rt_cfg.integrator);
        }
        basis = std::move(r.basis);

        // Lie-Trotter splitting for u: after TDVP moved {A,R} (and whatever
        // else was in alpha_z_list), evolve u exactly on the updated basis.
        if (rt_cfg.u_split_trotter) {
            auto [H_mat, S_mat] = build_HS(basis, perms_step, terms);
            free_evolve_fixed_basis(basis, H_mat, S_mat, r.used_dt);
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

        t += r.used_dt;
        step++;

        if (step % rt_cfg.sample_every == 0) {
            sample_observables(basis, t, terms, basis0, norm0, result.trace);
        }

        // Per-step diagnostic print (coarsened by print_every to avoid flooding)
        if (rt_cfg.verbose && (step % rt_cfg.print_every == 0)) {
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

            std::cout << "[rt_tdvp] step " << step;
            if (!adaptive) std::cout << "/" << n_steps_fixed;
            std::cout << " t=" << std::fixed << std::setprecision(4) << t
                      << " dt=" << std::scientific << std::setprecision(2) << r.used_dt
                      << " E=" << std::fixed << std::setprecision(8) << E
                      << " norm=" << S.real()
                      << " cond=" << std::scientific << std::setprecision(2)
                      << r.cond_C
                      << " sv_min3=(" << r.sv_small[2] << "," << r.sv_small[1]
                      << "," << r.sv_small[0] << ")"
                      << " |dz|=" << std::setprecision(2) << r.dz_norm
                      << " min(Re A+B)=" << std::fixed << std::setprecision(4) << min_rab
                      << "@" << argmin_i
                      << " min(Re A)=" << min_reA << "@" << argmin_A;
            if (adaptive) std::cout << " rej=" << rejected_count;
            std::cout << std::defaultfloat << "\n";
        }

        // Blowup guard: cond(C) >> 1 means the variational metric has a near-zero
        // singular value, i.e. two parameters have become degenerate. Dump basis
        // so we can see which one drifted, then stop.
        if (r.cond_C > cond_abort || !std::isfinite(r.cond_C)) {
            std::cout << "\n[rt_tdvp] BLOWUP at step " << step
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
            result.n_steps = step;
            return result;
        }
    }

    // Final sample (in case we don't land on sample_every)
    if (result.trace.t.empty() || result.trace.t.back() < t) {
        sample_observables(basis, t, terms, basis0, norm0, result.trace);
    }

    result.basis_final = std::move(basis);
    result.n_steps = step;
    if (adaptive && rt_cfg.verbose) {
        std::cout << "[rt_tdvp/RK45] total steps=" << step
                  << ", rejected=" << rejected_count
                  << ", final dt=" << std::scientific << dt_current
                  << std::defaultfloat << "\n";
    }
    return result;
}

} // namespace ecg1d
