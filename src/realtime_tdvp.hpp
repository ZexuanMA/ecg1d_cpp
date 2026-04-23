#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include "tdvp_solver.hpp"
#include <vector>
#include <limits>

namespace ecg1d {

// Real-time TDVP step for Schrodinger equation i*d|psi>/dt = H|psi>:
//   i * C_{alpha,beta} * z_dot_beta = g_alpha
//   -> dz = -i * C^{-1} * g    (C here is the tdvp metric with P_perp)
// vs imaginary-time version  (C * z_dot = -g  ->  dz = -C^{-1} * g)
//
// `dt` is a real positive step; integrator options (order 1 Euler or order 4 RK4).
//
// Returns dz as VectorXcd of size = alpha_z_list.size(); caller applies to basis.
enum class RtIntegrator { Euler, RK4, RK45 };

struct RealtimeStepResult {
    std::vector<BasisParams> basis;   // updated basis
    double used_dt;
    double cond_C;                    // diagnostic
    int    effective_rank;
    // Diagnostic: smallest 3 singular values of C_bar at step entry (k1 stage).
    // sv_small[0]  = smallest, sv_small[1]  = second smallest, sv_small[2] = third.
    // NaN if fewer than that many singular values exist.
    double sv_small[3] = {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN()};
    // Diagnostic: norm squared of dz returned by the full step (combined RK4 dz).
    double dz_norm = 0.0;
};

// Apply one real-time TDVP step. Does NOT adaptively rescale dt.
RealtimeStepResult realtime_tdvp_step(const std::vector<AlphaIndex>& alpha_z_list,
                                      const std::vector<BasisParams>& basis,
                                      double dt,
                                      const HamiltonianTerms& terms,
                                      const SolverConfig& config,
                                      RtIntegrator integrator = RtIntegrator::RK4);

struct RealtimeTrace {
    std::vector<double> t;
    std::vector<double> E;        // <H> / <psi|psi>  (real, conserved)
    std::vector<double> norm;     // <psi|psi> (real, conserved)
    std::vector<double> x2;       // <x^2>
    std::vector<double> p2;       // <p^2>
    std::vector<Cd>     overlap0; // <psi(0)|psi(t)> — fidelity |<psi(0)|psi(t)>|^2/(n0*n(t))
};

struct RealtimeEvolutionConfig {
    double dt          = 1e-3;
    RtIntegrator integrator = RtIntegrator::RK4;
    int    sample_every = 10;        // record trace every N steps
    bool   verbose      = true;
    int    print_every  = 100;

    // If true, per step do Lie-Trotter splitting:
    //   (1) TDVP moves only the parameters listed in alpha_z_list (caller's job
    //       to have excluded u/B from the list)
    //   (2) u is then propagated by exp(-i H dt) on the updated {A,B,R} basis
    //       via free_evolve_fixed_basis (exact, eigendecomposition-based).
    // This decouples the non-linear {A,R} flow from the linear u flow; avoids
    // the pathology where TDVP tries to absorb u dynamics into {A,R} drift.
    bool u_split_trotter = false;

    // RK45 (Dormand-Prince 5(4)) adaptive stepping — ignored unless integrator=RK45
    //   rk45_abs_tol: absolute tolerance
    //   rk45_rel_tol: relative tolerance
    //   rk45_dt_min/max: clamp on step size
    //   The initial dt is taken from rt_cfg.dt.
    double rk45_abs_tol = 1e-6;
    double rk45_rel_tol = 1e-6;
    double rk45_dt_min  = 1e-8;
    double rk45_dt_max  = 1.0;

    // If true, after each step rescale u so that <psi|psi> == initial_norm
    // (where initial_norm is computed at t=0 from basis_init). This is a
    // pragmatic post-correction that enforces physical unitarity when the raw
    // TDVP+SVD step drifts norm due to rcond truncation or Tikhonov λ not
    // being exactly consistent with ψ̇ ⊥ ψ. E conservation is NOT enforced,
    // only norm.
    bool enforce_norm = false;
};

struct RealtimeEvolutionResult {
    std::vector<BasisParams> basis_final;
    RealtimeTrace            trace;
    int    n_steps;
};

// Evolve `basis_init` from t=0 to t=T_total using real-time TDVP.
// At every sample_every'th step, record observables into RealtimeTrace.
RealtimeEvolutionResult realtime_tdvp_evolution(
    const std::vector<AlphaIndex>& alpha_z_list,
    std::vector<BasisParams> basis_init,
    double T_total,
    const HamiltonianTerms& terms,
    const SolverConfig& solver_cfg = SolverConfig::defaults(),
    const RealtimeEvolutionConfig& rt_cfg = {});

} // namespace ecg1d
