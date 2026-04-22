#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include "tdvp_solver.hpp"
#include <vector>

namespace ecg1d {

// Real-time TDVP step for Schrodinger equation i*d|psi>/dt = H|psi>:
//   i * C_{alpha,beta} * z_dot_beta = g_alpha
//   -> dz = -i * C^{-1} * g    (C here is the tdvp metric with P_perp)
// vs imaginary-time version  (C * z_dot = -g  ->  dz = -C^{-1} * g)
//
// `dt` is a real positive step; integrator options (order 1 Euler or order 4 RK4).
//
// Returns dz as VectorXcd of size = alpha_z_list.size(); caller applies to basis.
enum class RtIntegrator { Euler, RK4 };

struct RealtimeStepResult {
    std::vector<BasisParams> basis;   // updated basis
    double used_dt;
    double cond_C;                    // diagnostic
    int    effective_rank;
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
