#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include "permutation.hpp"
#include <vector>
#include <tuple>

namespace ecg1d {

// Alpha parameter index: (a1, a2, a3, a4)
// a1=1: u, a1=2: B, a1=3: R, a1=4: A
struct AlphaIndex {
    int a1, a2, a3, a4;
};

// Configuration: which Hamiltonian terms to include
struct HamiltonianTerms {
    bool kinetic  = true;
    bool harmonic = true;
    bool delta    = true;
    bool gaussian = true;
    bool kicking  = true;

    static HamiltonianTerms all() { return {true, true, true, true, true}; }
    static HamiltonianTerms kinetic_harmonic() { return {true, true, false, false, false}; }
};

// Solver configuration
struct SolverConfig {
    double lambda_C   = 1e-8;    // Tikhonov regularization for C matrix
    double rcond      = 1e-4;    // SVD truncation threshold
    bool   resolve_u  = false;   // Re-solve u via eigenvalue after each step
    double dtao_grow  = 1.5;     // Step size growth factor after successful step
    double dtao_max   = 0.0;     // Max step size (0 = 100*initial dtao)

    // Energy-based convergence
    double energy_tol        = 1e-12;  // Stop when |dE| < energy_tol for consecutive steps
    int    stagnation_window = 50;     // Steps without improvement before recovery attempt

    // Adaptive regularization
    bool   adaptive_lambda = false;    // Auto-increase lambda_C when condition number too high
    double lambda_max      = 1e-4;    // Upper bound for adaptive regularization

    // Parameter selection
    bool optimize_A = true;   // Include A matrix parameters
    bool optimize_B = false;  // Include B diagonal parameters (only for dynamics)
    bool optimize_R = false;  // Include R vector parameters (only for dynamics)
    // Note: B is redundant with A diagonal at B=0; R derivative is zero at B=0.
    // For static ground state, only A matters. For real-time dynamics, enable all three.

    static SolverConfig defaults() { return {}; }
    static SolverConfig dynamics() {
        SolverConfig c;
        c.optimize_B = true;
        c.optimize_R = true;
        return c;
    }
};

// Diagnostics from a TDVP step
struct TdvpDiagnostics {
    double cond_C       = 0.0;   // Condition number of C (max_sv / min_sv)
    int    effective_rank = 0;    // Number of singular values above threshold
    double min_sv_S     = 0.0;   // Smallest singular value of overlap matrix S
    double max_sv_C     = 0.0;   // Largest singular value of C_bar_update
    double min_sv_C     = 0.0;   // Smallest singular value of C_bar_update
};

// Assemble the full C metric tensor matrix
MatrixXcd assemble_C(const std::vector<AlphaIndex>& alpha_z_list,
                     const std::vector<BasisParams>& basis);

// Compute gradient of total Hamiltonian for a single alpha parameter
Cd grad_H_for_alpha(const AlphaIndex& alpha, bool Real,
                    const std::vector<BasisParams>& basis,
                    const HamiltonianTerms& terms = HamiltonianTerms::all());

// Assemble gradient vector
VectorXcd assemble_grad(const std::vector<AlphaIndex>& alpha_z_list,
                        bool Real,
                        const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms = HamiltonianTerms::all());

// Update basis function parameters with dz step
void update_basis_function(std::vector<BasisParams>& basis,
                           const VectorXcd& dz, double dtao,
                           const std::vector<AlphaIndex>& alpha_z_list,
                           int updata_constant);

// Compute total energy
Cd compute_total_energy(const std::vector<BasisParams>& basis,
                        const HamiltonianTerms& terms = HamiltonianTerms::all());

// Single TDVP step result
struct TdvpStepResult {
    std::vector<BasisParams> basis;
    Cd E_new, E_pre;
    double norm_dz;
    double used_dtao;
    TdvpDiagnostics diag;
};

TdvpStepResult tdvp_step(const std::vector<AlphaIndex>& alpha_z_list,
                          const std::vector<BasisParams>& basis,
                          double dtao,
                          const HamiltonianTerms& terms = HamiltonianTerms::all(),
                          const SolverConfig& config = SolverConfig::defaults(),
                          const PermutationSet* perms = nullptr);

// Main evolution loop
void evolution(const std::vector<AlphaIndex>& alpha_z_list,
               std::vector<BasisParams>& basis,
               double dtao = 1e-3,
               int max_steps = 10000,
               double tol = 1e-12,
               const HamiltonianTerms& terms = HamiltonianTerms::all(),
               const SolverConfig& config = SolverConfig::defaults(),
               const PermutationSet* perms = nullptr);

} // namespace ecg1d
