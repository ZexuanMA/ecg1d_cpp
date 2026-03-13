#pragma once
#include "types.hpp"
#include "pair_cache.hpp"
#include "permutation.hpp"
#include "basis_params.hpp"
#include <vector>
#include <array>

namespace ecg1d {

// 12-case derivative array
using Case12 = std::array<Cd, 12>;
using DerivTable = std::vector<Case12>;

int upper_index_inclusive(int N, int i, int j);

// --- Log-overlap first derivatives (12-case tables) ---
// Returns N entries for B/R, N*(N+1)/2 entries for A (upper-triangular order)

DerivTable dln_m_g_b(const PairCache& c, const VectorXi& sigma,
                     const BasisParams& pi, const BasisParams& pj);

DerivTable dln_m_g_r(const PairCache& c, const VectorXi& sigma,
                     const BasisParams& pi, const BasisParams& pj);

DerivTable dln_m_g_a(const PairCache& c, const VectorXi& sigma);

// --- K-inverse derivatives (12-case tables) ---

DerivTable dK_fu_B(int t, int o, const PairCache& c);
DerivTable dK_fu_R(int N);  // all zeros
DerivTable dK_fu_A(int t, int o, const PairCache& c);

// --- (K_inv @ b) derivatives (12-case tables) ---

DerivTable dK_fu_b_B(int t, const PairCache& c, const VectorXi& sigma,
                     const BasisParams& pi);
DerivTable dK_fu_b_R(int t, const PairCache& c,
                     const BasisParams& pi);
DerivTable dK_fu_b_A(int t, const PairCache& c);

// --- Case index resolution ---
// Maps (alpha_2, name_i, name_j, Real) -> (last_index_1, last_index_2)
// For first derivatives (holomorphic Real=True / anti-holomorphic Real=False)
std::pair<int,int> resolve_case_first(int alpha_2, int name_i, int name_j, bool Real);

// For second derivatives (beta side, always anti-holomorphic)
std::pair<int,int> resolve_case_second(int beta_2, int name_i, int name_j);

// --- partial_z_lnmg: per-permutation log-overlap derivative ---
// Returns SN-length vector (one value per permutation)
std::vector<Cd> partial_z_lnmg(int alpha_1, bool Real,
                                const BasisParams& pi, const BasisParams& pj,
                                int alpha_2, int alpha_3, int alpha_4,
                                const PermutationSet& perms);

// --- partial_z_lnmg_second: per-permutation second log-overlap derivative ---
std::vector<Cd> partial_z_lnmg_second(int alpha_1, int beta_1,
                                       const BasisParams& pi, const BasisParams& pj,
                                       int alpha_2, int alpha_3, int alpha_4,
                                       int beta_2, int beta_3, int beta_4,
                                       const PermutationSet& perms);

// --- partial_z_addsn: sum over permutations for first derivative ---
Cd partial_z_addsn(int alpha_1, bool Real,
                   const BasisParams& pi, const BasisParams& pj,
                   int alpha_2, int alpha_3, int alpha_4,
                   const PermutationSet& perms);

// --- partial_z_addsn_second: sum over permutations for second derivative ---
Cd partial_z_addsn_second(int alpha_1, int beta_1,
                          const BasisParams& pi, const BasisParams& pj,
                          int alpha_2, int alpha_3, int alpha_4,
                          int beta_2, int beta_3, int beta_4,
                          const PermutationSet& perms);

// --- partial_z_first: full first derivative over all basis functions ---
Cd partial_z_first(int alpha_1, bool Real,
                   const std::vector<BasisParams>& basis,
                   int alpha_2, int alpha_3, int alpha_4);

// --- partial_z_second: full second derivative over all basis functions ---
Cd partial_z_second(int alpha_1, int beta_1,
                    const std::vector<BasisParams>& basis,
                    int alpha_2, int alpha_3, int alpha_4,
                    int beta_2, int beta_3, int beta_4);

// --- calculate_C: single C matrix element ---
Cd calculate_C(int alpha_1, int alpha_2, int alpha_3, int alpha_4,
               int beta_1, int beta_2, int beta_3, int beta_4,
               const std::vector<BasisParams>& basis);

} // namespace ecg1d
