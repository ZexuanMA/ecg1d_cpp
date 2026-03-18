#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "tdvp_solver.hpp"
#include "permutation.hpp"
#include <vector>

namespace ecg1d {

// Parameters for kicked evolution
struct KickParams {
    double T_period  = 1.0;   // kick period
    double T_pulse   = 0.0;   // pulse duration (0 = ideal delta kick)
    double kappa_kick = 1.0;  // kick strength
    double k_L_kick  = 0.5;   // lattice wavevector
    int    n_kicks   = 10;    // number of kicks
};

// Observable snapshot at a given time
struct Observable {
    double time;
    double energy;
    double kinetic_energy;
    double overlap_norm;       // <psi|psi> for normalization check
};

// Real-time TDVP step result
struct RealtimeStepResult {
    std::vector<BasisParams> basis;
    Cd E_new;
    double norm_dz;
    double used_dt;
    TdvpDiagnostics diag;
};

// Real-time TDVP step: solves C dz = -i g (no energy-monotonic line search)
RealtimeStepResult realtime_tdvp_step(
    const std::vector<AlphaIndex>& alpha_z_list,
    const std::vector<BasisParams>& basis,
    double dt,
    const HamiltonianTerms& terms,
    const SolverConfig& config = SolverConfig::defaults(),
    const PermutationSet* perms = nullptr);

// Free evolution for duration T using real-time TDVP steps
void free_evolve(std::vector<BasisParams>& basis,
                 const std::vector<AlphaIndex>& alpha_z_list,
                 double dt_step, double T_duration,
                 const HamiltonianTerms& terms,
                 const SolverConfig& config = SolverConfig::defaults(),
                 const PermutationSet* perms = nullptr);

// Apply one kick (short pulse of duration T_pulse, or instantaneous if T_pulse ~ 0)
void apply_kick(std::vector<BasisParams>& basis,
                const std::vector<AlphaIndex>& alpha_z_list,
                double dt_step, double T_pulse,
                const HamiltonianTerms& terms_with_kick,
                const SolverConfig& config = SolverConfig::defaults(),
                const PermutationSet* perms = nullptr);

// Full kicked evolution: alternating free + kick for n_kicks periods
// Returns observables recorded after each kick
std::vector<Observable> kicked_evolution(
    std::vector<BasisParams>& basis,
    const std::vector<AlphaIndex>& alpha_z_list,
    const KickParams& kick_params,
    double dt_step,
    const HamiltonianTerms& terms_free,
    const HamiltonianTerms& terms_kick,
    const SolverConfig& config = SolverConfig::defaults(),
    const PermutationSet* perms = nullptr);

} // namespace ecg1d
