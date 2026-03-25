#pragma once
#include <vector>
#include <Eigen/Dense>

namespace ecg1d {

// Result of exact kicked rotor simulation for one particle
struct KickedExactResult {
    std::vector<double> kick_energies;    // E after each kick (single particle)
    std::vector<double> kick_kinetic;     // <p^2>/2 after each kick
    std::vector<double> kick_norm;        // |<psi|psi>| (should stay 1)
};

// Exact single-particle kicked rotor simulation via finite-difference method.
// H_free = -1/(2m) d²/dz² + 1/2 m ω² z²
// Kick = κ cos(2 k_L z) applied instantaneously each period T.
//
// Method:
//   1. Diagonalize H_free on a grid: H_free = U Λ U†
//   2. Propagator: exp(-i H_free T) = U exp(-iΛT) U†
//   3. Kick operator: ψ(z) → exp(-i κ cos(2 k_L z)) ψ(z) (pointwise)
//   4. Repeat for n_kicks periods
//
// N=2 non-interacting result is simply 2× the single-particle energy.
KickedExactResult kicked_exact_1particle(
    double omega, double mass,
    double k_L, double kappa,
    double T_period,
    int n_kicks,
    int N_grid = 512,
    double z_max = 15.0);

// Run N=2 non-interacting exact kicked dynamics and print results
void run_kicked_exact_test(int n_kicks = 50);

} // namespace ecg1d
