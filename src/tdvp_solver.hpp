#pragma once
#include "types.hpp"
#include "basis_params.hpp"
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
};

TdvpStepResult tdvp_step(const std::vector<AlphaIndex>& alpha_z_list,
                          const std::vector<BasisParams>& basis,
                          double dtao,
                          const HamiltonianTerms& terms = HamiltonianTerms::all());

// Main evolution loop
void evolution(const std::vector<AlphaIndex>& alpha_z_list,
               std::vector<BasisParams>& basis,
               double dtao = 1e-3,
               int max_steps = 10000,
               double tol = 1e-12,
               const HamiltonianTerms& terms = HamiltonianTerms::all());

} // namespace ecg1d
