#pragma once
#include "types.hpp"
#include "pair_cache.hpp"
#include "basis_params.hpp"
#include "physical_constants.hpp"
#include "derivatives.hpp"
#include <vector>
#include <array>

namespace ecg1d {

// ============================================================
// Observable kernel derivative tables
//
// Each function returns a DerivTable (vector of Case12) where:
//   - B and R functions return N entries (one per particle)
//   - A functions return N*(N+1)/2 entries in upper-triangular order:
//       diagonal (m,m) then off-diagonal (m,n) for n>m, for m=0..N-1
//
// Case12 layout:
//   For B/R (12-case first-derivative indexing):
//     [0]=case_1, [1]=case_2, [2]=case_3, [3]=case_4,
//     [4]=case_5, [5]=case_6, [6]=case_7, [7]=case_8,
//     [8]=case_9, [9]=case_10,[10]=case_11,[11]=case_12
//
//   For A diagonal (m,m):
//     [0]=case_5, [1]=case_6, [2]=case_7, [3]=case_8,
//     [4]=case_13,[5]=case_14,[6]=case_15,[7]=case_16,
//     [8]=case_19,[9]=case_20,[10]=case_23,[11]=case_24
//
//   For A off-diagonal (m,n), n>m:
//     [0]=case_1, [1]=case_2, [2]=case_3, [3]=case_4,
//     [4]=case_9, [5]=case_10,[6]=case_11,[7]=case_12,
//     [8]=case_17,[9]=case_18,[10]=case_21,[11]=case_22
// ============================================================

// --- P_Mij (kinetic energy kernel) derivatives ---

DerivTable dP_Mij_b(const PairCache& c, const VectorXi& sigma,
                    const MatrixXcd& perm_matrix_cd,
                    const BasisParams& pi, const BasisParams& pj);

DerivTable dP_Mij_r(const PairCache& c, const VectorXi& sigma,
                    const BasisParams& pi, const BasisParams& pj);

DerivTable dP_Mij_a(const PairCache& c, const VectorXi& sigma,
                    const MatrixXcd& perm_matrix_cd);

// --- rTr_Mij (harmonic potential kernel) derivatives ---

DerivTable drTr_Mij_b(const PairCache& c, const VectorXi& sigma,
                      const BasisParams& pi, const BasisParams& pj);

DerivTable drTr_Mij_r(const PairCache& c, const VectorXi& sigma,
                      const BasisParams& pi, const BasisParams& pj);

DerivTable drTr_Mij_a(const PairCache& c, const VectorXi& sigma);

// --- G_Mijab (delta contact kernel) derivatives (per pair a,b) ---

DerivTable dG_Mijab_b(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj);

DerivTable dG_Mijab_r(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj);

DerivTable dG_Mijab_a(const PairCache& c, const VectorXi& sigma, int a, int b);

// --- H_Mijab (Gaussian interaction kernel) derivatives (per pair a,b) ---

DerivTable dH_Mijab_b(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj);

DerivTable dH_Mijab_r(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj);

DerivTable dH_Mijab_a(const PairCache& c, const VectorXi& sigma, int a, int b);

// --- Q_Mija (kicking term kernel) derivatives (per particle a) ---

DerivTable dQ_Mija_b(const PairCache& c, const VectorXi& sigma, int a,
                     const BasisParams& pi, const BasisParams& pj);

DerivTable dQ_Mija_r(const PairCache& c, const VectorXi& sigma, int a,
                     const BasisParams& pi, const BasisParams& pj);

DerivTable dQ_Mija_a(const PairCache& c, const VectorXi& sigma, int a);

} // namespace ecg1d
