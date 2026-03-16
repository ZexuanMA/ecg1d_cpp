#pragma once
#include "types.hpp"
#include "basis_params.hpp"

namespace ecg1d {

struct PairCache {
    int N;
    MatrixXcd K;       // (N,N) full K matrix
    MatrixXcd K_inv;   // K^{-1}
    VectorXcd b;       // (N,) linear term
    Cd C_val;          // scalar constant
    Cd det_K;          // det(K)
    Cd M_G;            // overlap matrix element
    VectorXcd mu;      // 0.5 * K_inv * b (mean position)
    MatrixXcd K_Mj;    // P^T A_j P + P^T B_j P (ket-side K)
    VectorXcd g_Mj;    // 2 * R_j^T B_j P (ket-side linear term, as column vector)

    static PairCache build(const BasisParams& pi, const BasisParams& pj,
                           const MatrixXi& perm_matrix);
};

} // namespace ecg1d
