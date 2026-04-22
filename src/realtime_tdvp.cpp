#include "realtime_tdvp.hpp"
#include "hamiltonian.hpp"
#include "hamiltonian_gradient.hpp"
#include "physical_constants.hpp"
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
                                int* out_rank = nullptr) {
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

    VectorXcd Utrhs = svd.matrixU().adjoint() * rhs;
    VectorXcd sinv(sv.size());
    for (int i = 0; i < sv.size(); ++i) {
        sinv(i) = (sv(i) > thr) ? Utrhs(i) / sv(i) : Cd(0.0, 0.0);
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

    VectorXcd dz_combined;
    if (integrator == RtIntegrator::Euler) {
        dz_combined = compute_rhs_dz(alpha_z_list, basis, terms, config,
                                     &cond_C, &rank);
    } else {
        // Classical RK4 on z: z' = f(z) = -i C^{-1} g
        VectorXcd k1 = compute_rhs_dz(alpha_z_list, basis, terms, config,
                                      &cond_C, &rank);

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
    double t = 0.0;
    for (int step = 0; step < n_steps; ++step) {
        RealtimeStepResult r = realtime_tdvp_step(alpha_z_list, basis, rt_cfg.dt,
                                                  terms, solver_cfg, rt_cfg.integrator);
        basis = std::move(r.basis);
        t += rt_cfg.dt;

        if ((step + 1) % rt_cfg.sample_every == 0) {
            sample_observables(basis, t, terms, basis0, norm0, result.trace);
        }

        if (rt_cfg.verbose && ((step + 1) % rt_cfg.print_every == 0)) {
            Cd S = overlap(basis);
            Cd T_kin = kinetic_energy_functional(basis);
            Cd V_ho  = Harmonic_functional(basis);
            double E = ((T_kin + V_ho) / S).real();
            std::cout << "[rt_tdvp] step " << (step + 1) << "/" << n_steps
                      << " t=" << std::fixed << std::setprecision(4) << t
                      << " E=" << std::setprecision(8) << E
                      << " norm=" << S.real()
                      << " cond(C)=" << std::scientific << std::setprecision(2)
                      << r.cond_C
                      << std::defaultfloat << "\n";
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
