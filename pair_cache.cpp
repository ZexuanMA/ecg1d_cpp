#include "pair_cache.hpp"
#include <cmath>

namespace ecg1d {

PairCache PairCache::build(const BasisParams& pi, const BasisParams& pj,
                           const MatrixXi& perm_matrix) {
    PairCache c;
    c.N = pi.N();
    int N = c.N;

    BasisParams pi_con = pi.conj_params();

    // Permutation matrix as double for arithmetic
    MatrixXcd P = perm_matrix.cast<double>().cast<Cd>();
    MatrixXcd PT = P.transpose();

    // K = conj(A_i) + conj(B_i) + P^T A_j P + P^T B_j P
    MatrixXcd AA = PT * pj.A * P;
    MatrixXcd BB = PT * pj.B * P;
    c.K = pi_con.A + pi_con.B + AA + BB;

    // K_Mj = P^T A_j P + P^T B_j P (ket side only)
    c.K_Mj = AA + BB;

    // b: bT = 2*conj(R_i)^T*conj(B_i) + 2*R_j^T*B_j*P
    // bT is a row vector; b = bT^T is column
    // conj(R_i)^T * conj(B_i) is a (1,N) row vector
    // R_j^T * B_j * P is a (1,N) row vector
    Eigen::RowVectorXcd ii = pi_con.R.transpose() * pi_con.B;  // (1,N)
    Eigen::RowVectorXcd jj = pj.R.transpose() * pj.B * P;      // (1,N)
    Eigen::RowVectorXcd bT = 2.0 * ii + 2.0 * jj;
    c.b = bT.transpose();  // column vector

    // g_Mj: g_Mj_T = 2 * R_j^T B_j P (row vector); g_Mj = column
    c.g_Mj = (2.0 * jj).transpose();

    // C = -conj(R_i)^T conj(B_i) conj(R_i) - R_j^T B_j R_j
    Cd C_ii = (pi_con.R.transpose() * pi_con.B * pi_con.R)(0);
    Cd C_jj = (pj.R.transpose() * pj.B * pj.R)(0);
    c.C_val = -C_ii - C_jj;

    // K_inv and det via LU
    Eigen::PartialPivLU<MatrixXcd> lu(c.K);
    c.K_inv = lu.inverse();
    c.det_K = lu.determinant();

    // mu = 0.5 * K_inv * b
    c.mu = 0.5 * c.K_inv * c.b;

    // M_G = pi^(N/2) * det(K)^(-1/2) * exp(C + 0.25 * bT * K_inv * b)
    double PI = std::pow(M_PI, N / 2.0);
    Cd sqrt_detK_inv = std::pow(c.det_K, -0.5);
    Cd exponent = c.C_val + 0.25 * (bT * c.K_inv * c.b)(0);
    c.M_G = PI * sqrt_detK_inv * std::exp(exponent);

    return c;
}

} // namespace ecg1d
