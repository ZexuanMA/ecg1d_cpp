#pragma once
#include "basis_params.hpp"
#include "permutation.hpp"
#include "pair_cache.hpp"
#include "tdvp_solver.hpp"  // HamiltonianTerms
#include <vector>

namespace ecg1d {

// Single-pair kernels, N=1 only (throw otherwise).
// Return ⟨φ_i|·|φ_j⟩ / M_G (i.e. kernel divided by overlap M_G).
Cd compute_T2_Mij(const PairCache& c);           // ⟨·|¼∂⁴|·⟩ / M_G
Cd compute_V2_Mij(const PairCache& c);           // ⟨·|¼x⁴|·⟩ / M_G
Cd compute_TV_plus_VT_Mij(const PairCache& c);   // ⟨·|-¼(∂²x²+x²∂²)|·⟩ / M_G

// Single matrix element ⟨φ_i|Ĥ²|φ_j⟩ summed over permutations.
// terms selects which H² components to include; only kinetic+harmonic supported in M1.
Cd compute_H2_ij(const BasisParams& bi, const BasisParams& bj,
                 const PermutationSet& perms,
                 const HamiltonianTerms& terms);

// K×K Hermitian assembly (upper-triangular loop + conj mirror).
MatrixXcd build_H2(const std::vector<BasisParams>& basis,
                   const PermutationSet& perms,
                   const HamiltonianTerms& terms);

} // namespace ecg1d
