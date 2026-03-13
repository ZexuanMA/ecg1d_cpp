#include "observable_derivatives.hpp"
#include "interaction_kernels.hpp"
#include <cmath>
#include <complex>

namespace ecg1d {

// ============================================================
// Internal helpers
// ============================================================

static Case12 zeros12() {
    Case12 c{};
    for (auto& v : c) v = Cd(0, 0);
    return c;
}

// Find sigma_fu = sigma^{-1}(i): index z such that sigma(z) == i
static int sigma_inv(const VectorXi& sigma, int i) {
    int N = sigma.size();
    for (int z = 0; z < N; z++)
        if (sigma(z) == i) return z;
    return -1; // should never happen
}

// ============================================================
// 1. P_Mij (kinetic energy kernel) derivatives
// ============================================================

DerivTable dP_Mij_b(const PairCache& c, const VectorXi& sigma,
                    const MatrixXcd& perm_matrix_cd,
                    const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    // Precompute shared matrices
    // K_inv @ b = 2 * mu
    VectorXcd K_inv_b = 2.0 * c.mu;

    // g_Mj_T @ K_Mj @ K_inv: row vector, index i
    // = (c.K_Mj.transpose() * c.K_Mj * c.K_inv).row? No:
    // g_Mj_T is a row vector: g_Mj.transpose()
    // (g_Mj_T @ K_Mj @ K_inv)[i] = g_Mj.transpose() * (K_Mj * K_inv), col i
    // As a row vector: g_Mj.transpose() * K_Mj * K_inv
    Eigen::RowVectorXcd gMj_T_KMj_Kinv = c.g_Mj.transpose() * c.K_Mj * c.K_inv;

    // mu_Mj_T @ K_Mj @ K_Mj @ K_inv: row vector
    Eigen::RowVectorXcd muMj_T_KMj2_Kinv = c.mu.transpose() * c.K_Mj * c.K_Mj * c.K_inv;

    // K_inv @ K_Mj @ K_Mj @ K_inv: matrix, need (i,i) diagonal
    MatrixXcd KinvKMj2Kinv = c.K_inv * c.K_Mj * c.K_Mj * c.K_inv;

    // K_Mj @ mu: vector
    VectorXcd KMj_mu = c.K_Mj * c.mu;

    // P^T @ ((K_Mj @ K_inv) + (K_inv @ K_Mj)) @ P: matrix
    // perm_matrix_cd is P (the permutation matrix as complex)
    MatrixXcd PT = perm_matrix_cd.transpose();
    MatrixXcd KMj_Kinv = c.K_Mj * c.K_inv;
    MatrixXcd Kinv_KMj = c.K_inv * c.K_Mj;
    MatrixXcd sym_mat = PT * (KMj_Kinv + Kinv_KMj) * perm_matrix_cd;

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        Cd R_i_con_i = pi_con.R(i);
        Cd Kinv_b_i  = K_inv_b(i);
        Cd gMj_KMj_Kinv_i = gMj_T_KMj_Kinv(i);
        Cd muMj_KMj2_Kinv_i = muMj_T_KMj2_Kinv(i);
        Cd KinvKMj2Kinv_ii = KinvKMj2Kinv(i, i);
        Cd KMj_mu_i = KMj_mu(i);
        Cd gMj_i = c.g_Mj(i);
        Cd mu_i = c.mu(i);
        Cd R_j_sf = pj.R(sigma_fu);

        // case_1 = 0
        // case_2: bra Real=False
        Cd case_2 = (-2.0 * KinvKMj2Kinv_ii
                     - 2.0 * gMj_KMj_Kinv_i * (R_i_con_i + R_i_con_i - Kinv_b_i)
                     + 4.0 * muMj_KMj2_Kinv_i * (R_i_con_i + R_i_con_i - Kinv_b_i));

        // case_3: ket Real=True (direct, uses P^T sym_mat P)
        Cd case_3 = -2.0 + 2.0 * sym_mat(i, i);

        // case_4,5,6 = 0
        // case_7: ket permuted
        Cd case_7 = (-2.0 * KinvKMj2Kinv_ii
                     + 2.0 * gMj_i * (R_j_sf + R_j_sf)
                     - 4.0 * KMj_mu_i * (R_j_sf + R_j_sf)
                     - 2.0 * gMj_KMj_Kinv_i * (R_j_sf + R_j_sf - Kinv_b_i)
                     - 4.0 * gMj_i * mu_i
                     + 8.0 * mu_i * KMj_mu_i
                     + 4.0 * muMj_KMj2_Kinv_i * (R_j_sf + R_j_sf - Kinv_b_i));

        // case_8 = 0
        // Combined
        Cd case_9  = Cd(0) + case_3;   // case_1 + case_3
        Cd case_10 = case_2 + Cd(0);   // case_2 + case_4
        Cd case_11 = Cd(0) + case_7;   // case_5 + case_7
        Cd case_12 = Cd(0);            // case_6 + case_8

        result[i] = {Cd(0), case_2, case_3, Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dP_Mij_r(const PairCache& c, const VectorXi& sigma,
                    const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    // g_Mj_T @ K_Mj @ K_inv: row vector
    Eigen::RowVectorXcd gMj_T_KMj_Kinv = c.g_Mj.transpose() * c.K_Mj * c.K_inv;
    // mu_Mj_T @ K_Mj @ K_Mj @ K_inv: row vector
    Eigen::RowVectorXcd muMj_T_KMj2_Kinv = c.mu.transpose() * c.K_Mj * c.K_Mj * c.K_inv;
    // K_Mj @ mu
    VectorXcd KMj_mu = c.K_Mj * c.mu;

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        Cd B_i_con_ii = pi_con.B(i, i);
        Cd B_j_sf_sf  = pj.B(sigma_fu, sigma_fu);
        Cd gMj_KMj_Kinv_i   = gMj_T_KMj_Kinv(i);
        Cd muMj_KMj2_Kinv_i = muMj_T_KMj2_Kinv(i);
        Cd gMj_i   = c.g_Mj(i);
        Cd KMj_mu_i = KMj_mu(i);

        // case_2: bra Real=False
        Cd case_2 = (-4.0 * gMj_KMj_Kinv_i * B_i_con_ii
                     + 8.0 * muMj_KMj2_Kinv_i * B_i_con_ii);

        // case_3,4,5,6 = 0
        // case_7: ket permuted
        Cd case_7 = (-4.0 * gMj_KMj_Kinv_i * B_j_sf_sf
                     + 8.0 * muMj_KMj2_Kinv_i * B_j_sf_sf
                     + 4.0 * gMj_i * B_j_sf_sf
                     - 8.0 * KMj_mu_i * B_j_sf_sf);

        // case_8 = 0
        Cd case_9  = Cd(0);           // case_1(=0) + case_3(=0)
        Cd case_10 = case_2;          // case_2 + case_4(=0)
        Cd case_11 = case_7;          // case_5(=0) + case_7
        Cd case_12 = Cd(0);           // case_6(=0) + case_8(=0)

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dP_Mij_a(const PairCache& c, const VectorXi& sigma,
                    const MatrixXcd& perm_matrix_cd) {
    int N = c.N;
    DerivTable result;

    // Precompute shared matrices
    VectorXcd K_inv_b = 2.0 * c.mu;
    Eigen::RowVectorXcd gMj_T_KMj_Kinv = c.g_Mj.transpose() * c.K_Mj * c.K_inv;
    Eigen::RowVectorXcd muMj_T_KMj2_Kinv = c.mu.transpose() * c.K_Mj * c.K_Mj * c.K_inv;
    MatrixXcd KinvKMj2Kinv = c.K_inv * c.K_Mj * c.K_Mj * c.K_inv;
    VectorXcd KMj_mu = c.K_Mj * c.mu;

    // P^T @ ((K_Mj @ K_inv) + (K_inv @ K_Mj)) @ P
    MatrixXcd PT = perm_matrix_cd.transpose();
    MatrixXcd sym_mat = PT * (c.K_Mj * c.K_inv + c.K_inv * c.K_Mj) * perm_matrix_cd;

    for (int m = 0; m < N; m++) {
        Cd Kinv_b_m = K_inv_b(m);
        Cd gMj_KMj_Kinv_m   = gMj_T_KMj_Kinv(m);
        Cd muMj_KMj2_Kinv_m = muMj_T_KMj2_Kinv(m);
        Cd KinvKMj2Kinv_mm  = KinvKMj2Kinv(m, m);
        Cd gMj_m   = c.g_Mj(m);
        Cd mu_m    = c.mu(m);
        Cd KMj_mu_m = KMj_mu(m);

        // Diagonal (m, m)
        // case_5=0, case_6 (index[1]), case_7 (index[2]), case_8=0
        // case_13=0, case_14=0, case_15 (index[6]), case_16=0
        Cd case_6 = (-2.0 * KinvKMj2Kinv_mm
                     + 2.0 * gMj_KMj_Kinv_m * Kinv_b_m
                     - 4.0 * muMj_KMj2_Kinv_m * Kinv_b_m);

        Cd case_7 = -2.0 + 2.0 * sym_mat(m, m);

        Cd case_15 = (case_6
                      - 4.0 * gMj_m * mu_m
                      + 8.0 * mu_m * KMj_mu_m);

        // Combined: case_19=case_5+case_7, case_20=case_6+case_8, case_23=case_13+case_15, case_24=case_14+case_16
        Cd case_19 = Cd(0) + case_7;    // case_5(=0) + case_7
        Cd case_20 = case_6 + Cd(0);    // case_6 + case_8(=0)
        Cd case_23 = Cd(0) + case_15;   // case_13(=0) + case_15
        Cd case_24 = Cd(0);             // case_14(=0) + case_16(=0)

        // Store as [case_5, case_6, case_7, case_8, case_13, case_14, case_15, case_16, case_19, case_20, case_23, case_24]
        result.push_back({Cd(0), case_6, case_7, Cd(0),
                          Cd(0), Cd(0),  case_15, Cd(0),
                          case_19, case_20, case_23, case_24});

        // Off-diagonal (m, n), n > m
        for (int n = m + 1; n < N; n++) {
            Cd Kinv_b_n = K_inv_b(n);
            Cd gMj_KMj_Kinv_n   = gMj_T_KMj_Kinv(n);
            Cd muMj_KMj2_Kinv_n = muMj_T_KMj2_Kinv(n);
            Cd KinvKMj2Kinv_nm  = KinvKMj2Kinv(n, m);
            Cd KinvKMj2Kinv_mn  = KinvKMj2Kinv(m, n);
            Cd gMj_n   = c.g_Mj(n);
            Cd mu_n    = c.mu(n);
            Cd KMj_mu_n = KMj_mu(n);

            // case_1=0, case_2, case_3, case_4=0
            Cd case_2 = (-2.0 * KinvKMj2Kinv_nm
                         - 2.0 * KinvKMj2Kinv_mn
                         + 2.0 * gMj_KMj_Kinv_m * Kinv_b_n
                         + 2.0 * gMj_KMj_Kinv_n * Kinv_b_m
                         - 4.0 * muMj_KMj2_Kinv_m * Kinv_b_n
                         - 4.0 * muMj_KMj2_Kinv_n * Kinv_b_m);

            Cd case_3 = 2.0 * sym_mat(m, n) + 2.0 * sym_mat(n, m);

            // case_9=0, case_10=0, case_11, case_12=0
            Cd case_11 = (case_2
                          - 4.0 * gMj_m * mu_n
                          - 4.0 * gMj_n * mu_m
                          + 8.0 * mu_m * KMj_mu_n
                          + 8.0 * mu_n * KMj_mu_m);

            // Combined: case_17=case_1+case_3, case_18=case_2+case_4, case_21=case_9+case_11, case_22=case_10+case_12
            Cd case_17 = Cd(0) + case_3;    // case_1(=0) + case_3
            Cd case_18 = case_2 + Cd(0);    // case_2 + case_4(=0)
            Cd case_21 = Cd(0) + case_11;   // case_9(=0) + case_11
            Cd case_22 = Cd(0);             // case_10(=0) + case_12(=0)

            // Store as [case_1, case_2, case_3, case_4, case_9, case_10, case_11, case_12, case_17, case_18, case_21, case_22]
            result.push_back({Cd(0), case_2, case_3, Cd(0),
                               Cd(0), Cd(0),  case_11, Cd(0),
                               case_17, case_18, case_21, case_22});
        }
    }
    return result;
}

// ============================================================
// 2. rTr_Mij (harmonic potential kernel) derivatives
// ============================================================

DerivTable drTr_Mij_b(const PairCache& c, const VectorXi& sigma,
                      const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    // mu_Mj_T @ K_inv: row vector
    Eigen::RowVectorXcd muMj_T_Kinv = c.mu.transpose() * c.K_inv;
    // K_inv @ b = 2 * mu
    VectorXcd K_inv_b = 2.0 * c.mu;
    // K_inv @ K_inv: matrix, need diagonal
    MatrixXcd Kinv2 = c.K_inv * c.K_inv;

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        Cd R_i_con_i = pi_con.R(i);
        Cd R_j_sf    = pj.R(sigma_fu);
        Cd Kinv_b_i  = K_inv_b(i);
        Cd muMj_Kinv_i = muMj_T_Kinv(i);
        Cd Kinv2_ii  = Kinv2(i, i);

        // case_2: bra Real=False
        Cd case_2 = (muMj_Kinv_i * (R_i_con_i + R_i_con_i)
                     - muMj_Kinv_i * Kinv_b_i
                     - 0.5 * Kinv2_ii);

        // case_3,4,5,6 = 0
        // case_7: ket permuted
        Cd case_7 = (muMj_Kinv_i * (R_j_sf + R_j_sf)
                     - muMj_Kinv_i * Kinv_b_i
                     - 0.5 * Kinv2_ii);

        // case_8 = 0
        Cd case_9  = Cd(0);     // case_1(=0) + case_3(=0)
        Cd case_10 = case_2;    // case_2 + case_4(=0)
        Cd case_11 = case_7;    // case_5(=0) + case_7
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable drTr_Mij_r(const PairCache& c, const VectorXi& sigma,
                      const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    // mu_Mj_T @ K_inv: row vector
    Eigen::RowVectorXcd muMj_T_Kinv = c.mu.transpose() * c.K_inv;

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        Cd B_i_con_ii = pi_con.B(i, i);
        Cd B_j_sf_sf  = pj.B(sigma_fu, sigma_fu);
        Cd muMj_Kinv_i = muMj_T_Kinv(i);

        // case_2: bra Real=False
        Cd case_2 = 2.0 * muMj_Kinv_i * B_i_con_ii;

        // case_7: ket permuted
        Cd case_7 = 2.0 * muMj_Kinv_i * B_j_sf_sf;

        Cd case_9  = Cd(0);
        Cd case_10 = case_2;
        Cd case_11 = case_7;
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable drTr_Mij_a(const PairCache& c, const VectorXi& sigma) {
    int N = c.N;
    DerivTable result;

    VectorXcd K_inv_b = 2.0 * c.mu;
    Eigen::RowVectorXcd muMj_T_Kinv = c.mu.transpose() * c.K_inv;
    MatrixXcd Kinv2 = c.K_inv * c.K_inv;

    for (int m = 0; m < N; m++) {
        Cd Kinv_b_m  = K_inv_b(m);
        Cd muMj_Kinv_m = muMj_T_Kinv(m);
        Cd Kinv2_mm  = Kinv2(m, m);

        // Diagonal (m, m)
        // case_6 (index[1])
        Cd case_6 = (-muMj_Kinv_m * Kinv_b_m - 0.5 * Kinv2_mm);
        // case_15 = same as case_6 (index[6])
        Cd case_15 = case_6;

        Cd case_19 = Cd(0);        // case_5(=0) + case_7(=0)
        Cd case_20 = case_6;       // case_6 + case_8(=0)
        Cd case_23 = case_15;      // case_13(=0) + case_15
        Cd case_24 = Cd(0);

        result.push_back({Cd(0), case_6, Cd(0), Cd(0),
                          Cd(0), Cd(0),  case_15, Cd(0),
                          case_19, case_20, case_23, case_24});

        // Off-diagonal (m, n), n > m
        for (int n = m + 1; n < N; n++) {
            Cd Kinv_b_n  = K_inv_b(n);
            Cd muMj_Kinv_n = muMj_T_Kinv(n);
            Cd Kinv2_mn  = Kinv2(m, n);
            Cd Kinv2_nm  = Kinv2(n, m);

            // case_2 (index[1])
            Cd case_2 = (-muMj_Kinv_m * Kinv_b_n
                         - muMj_Kinv_n * Kinv_b_m
                         - 0.5 * Kinv2_mn
                         - 0.5 * Kinv2_nm);
            // case_11 = same as case_2 (index[6])
            Cd case_11 = case_2;

            Cd case_17 = Cd(0);       // case_1(=0) + case_3(=0)
            Cd case_18 = case_2;      // case_2 + case_4(=0)
            Cd case_21 = case_11;     // case_9(=0) + case_11
            Cd case_22 = Cd(0);

            result.push_back({Cd(0), case_2, Cd(0), Cd(0),
                               Cd(0), Cd(0),  case_11, Cd(0),
                               case_17, case_18, case_21, case_22});
        }
    }
    return result;
}

// ============================================================
// Internal helper: compute G_Mijab from PairCache
// ============================================================

static Cd G_Mijab_from_cache(const PairCache& c, int a, int b) {
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);
    return 1.0 / std::sqrt(M_PI * h) * std::exp(-p * p / h);
}

// ============================================================
// 3. G_Mijab (delta contact kernel) derivatives
// ============================================================

DerivTable dG_Mijab_b(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    BasisParams pj_con = pj.conj_params();
    DerivTable result(N);

    Cd G = G_Mijab_from_cache(c, a, b);
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);

    Cd factor_h = G * (-0.5 / h + p * p / (h * h));
    Cd factor_p = G * (-2.0 * p / h);

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        // dh/dK_ii: derivative of h w.r.t. K[i][i]
        Cd dh_dKii = (-c.K_inv(a, i) * c.K_inv(i, a)
                      - c.K_inv(b, i) * c.K_inv(i, b)
                      + 2.0 * c.K_inv(a, i) * c.K_inv(i, b));

        // case_2: bra Real=False
        // dp/dB_bra[i]: R_i_con enters b via 2*R_i_con*B_i_con, so
        // dp[a]/dB_bra_ii = K_inv[a][i]*(R_i_con[i] - mu[i])
        // dp[b]/dB_bra_ii = K_inv[b][i]*(R_i_con[i] - mu[i]) (negated in dp)
        Cd dp_bra = (c.K_inv(a, i) * pi_con.R(i)
                     - c.K_inv(a, i) * c.mu(i)
                     - c.K_inv(b, i) * pi_con.R(i)
                     + c.K_inv(b, i) * c.mu(i));
        Cd case_2 = factor_h * dh_dKii + factor_p * dp_bra;

        // case_7: ket permuted (Real=True)
        // In Python, uses params_j_con.R[sigma_fu]
        Cd dp_ket = (c.K_inv(a, i) * pj_con.R(sigma_fu)
                     - c.K_inv(a, i) * c.mu(i)
                     - c.K_inv(b, i) * pj_con.R(sigma_fu)
                     + c.K_inv(b, i) * c.mu(i));
        Cd case_7 = factor_h * dh_dKii + factor_p * dp_ket;

        Cd case_9  = Cd(0);
        Cd case_10 = case_2;
        Cd case_11 = case_7;
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dG_Mijab_r(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    Cd G = G_Mijab_from_cache(c, a, b);
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);

    Cd factor_p = G * (-2.0 * p / h);

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        // case_2: bra Real=False
        Cd case_2 = factor_p * (c.K_inv(a, i) * pi_con.B(i, i)
                                 - c.K_inv(b, i) * pi_con.B(i, i));

        // case_7: ket permuted (Real=True)
        Cd case_7 = factor_p * (c.K_inv(a, i) * pj.B(sigma_fu, sigma_fu)
                                 - c.K_inv(b, i) * pj.B(sigma_fu, sigma_fu));

        Cd case_9  = Cd(0);
        Cd case_10 = case_2;
        Cd case_11 = case_7;
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dG_Mijab_a(const PairCache& c, const VectorXi& sigma, int a, int b) {
    int N = c.N;
    DerivTable result;

    Cd G = G_Mijab_from_cache(c, a, b);
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);

    Cd factor_h = G * (-0.5 / h + p * p / (h * h));
    Cd factor_p = G * (-2.0 * p / h);

    for (int m = 0; m < N; m++) {
        // Diagonal (m, m)
        Cd dh_diag = (-c.K_inv(a, m) * c.K_inv(m, a)
                      - c.K_inv(b, m) * c.K_inv(m, b)
                      + 2.0 * c.K_inv(a, m) * c.K_inv(m, b));
        Cd dp_diag = (-c.K_inv(a, m) * c.mu(m)
                      + c.K_inv(b, m) * c.mu(m));

        // case_6 (index[1])
        Cd case_6 = factor_h * dh_diag + factor_p * dp_diag;
        // case_15 = same (index[6])
        Cd case_15 = case_6;

        Cd case_19 = Cd(0);      // case_5(=0) + case_7(=0)
        Cd case_20 = case_6;     // case_6 + case_8(=0)
        Cd case_23 = case_15;    // case_13(=0) + case_15
        Cd case_24 = Cd(0);

        result.push_back({Cd(0), case_6, Cd(0), Cd(0),
                          Cd(0), Cd(0),  case_15, Cd(0),
                          case_19, case_20, case_23, case_24});

        // Off-diagonal (m, n), n > m
        for (int n = m + 1; n < N; n++) {
            Cd dh_off = (-c.K_inv(a, m) * c.K_inv(n, a)
                         - c.K_inv(a, n) * c.K_inv(m, a)
                         - c.K_inv(b, m) * c.K_inv(n, b)
                         - c.K_inv(b, n) * c.K_inv(m, b)
                         + 2.0 * c.K_inv(a, m) * c.K_inv(n, b)
                         + 2.0 * c.K_inv(a, n) * c.K_inv(m, b));
            Cd dp_off = (-c.K_inv(a, m) * c.mu(n)
                         - c.K_inv(a, n) * c.mu(m)
                         + c.K_inv(b, m) * c.mu(n)
                         + c.K_inv(b, n) * c.mu(m));

            // case_2 (index[1])
            Cd case_2 = factor_h * dh_off + factor_p * dp_off;
            // case_11 = same (index[6])
            Cd case_11 = case_2;

            Cd case_17 = Cd(0);
            Cd case_18 = case_2;
            Cd case_21 = case_11;
            Cd case_22 = Cd(0);

            result.push_back({Cd(0), case_2, Cd(0), Cd(0),
                               Cd(0), Cd(0),  case_11, Cd(0),
                               case_17, case_18, case_21, case_22});
        }
    }
    return result;
}

// ============================================================
// Internal helper: compute H_Mijab from PairCache
// ============================================================

static Cd H_Mijab_from_cache(const PairCache& c, int a, int b) {
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);
    Cd h_eff = sigma_gauss * sigma_gauss + h;
    return sigma_gauss / std::sqrt(h_eff) * std::exp(-p * p / h_eff);
}

// ============================================================
// 4. H_Mijab (Gaussian interaction kernel) derivatives
// ============================================================

DerivTable dH_Mijab_b(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    BasisParams pj_con = pj.conj_params();
    DerivTable result(N);

    Cd H = H_Mijab_from_cache(c, a, b);
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);
    Cd h_eff = sigma_gauss * sigma_gauss + h;

    Cd factor_h = H * (-0.5 / h_eff + p * p / (h_eff * h_eff));
    Cd factor_p = H * (-2.0 * p / h_eff);

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        // dh/dK_ii (same formula as G case)
        Cd dh_dKii = (-c.K_inv(a, i) * c.K_inv(i, a)
                      - c.K_inv(b, i) * c.K_inv(i, b)
                      + 2.0 * c.K_inv(a, i) * c.K_inv(i, b));

        // case_2: bra Real=False (uses conj of i params)
        Cd dp_bra = (c.K_inv(a, i) * pi_con.R(i)
                     - c.K_inv(a, i) * c.mu(i)
                     - c.K_inv(b, i) * pi_con.R(i)
                     + c.K_inv(b, i) * c.mu(i));
        Cd case_2 = factor_h * dh_dKii + factor_p * dp_bra;

        // case_7: ket permuted (Real=True, uses conj of j params)
        Cd dp_ket = (c.K_inv(a, i) * pj_con.R(sigma_fu)
                     - c.K_inv(a, i) * c.mu(i)
                     - c.K_inv(b, i) * pj_con.R(sigma_fu)
                     + c.K_inv(b, i) * c.mu(i));
        Cd case_7 = factor_h * dh_dKii + factor_p * dp_ket;

        Cd case_9  = Cd(0);
        Cd case_10 = case_2;
        Cd case_11 = case_7;
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dH_Mijab_r(const PairCache& c, const VectorXi& sigma, int a, int b,
                      const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    Cd H = H_Mijab_from_cache(c, a, b);
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);
    Cd h_eff = sigma_gauss * sigma_gauss + h;

    Cd factor_p = H * (-2.0 * p / h_eff);

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        // case_2: bra Real=False (uses conj of i params)
        Cd case_2 = factor_p * (c.K_inv(a, i) * pi_con.B(i, i)
                                 - c.K_inv(b, i) * pi_con.B(i, i));

        // case_7: ket permuted (Real=True, uses j params directly)
        Cd case_7 = factor_p * (c.K_inv(a, i) * pj.B(sigma_fu, sigma_fu)
                                 - c.K_inv(b, i) * pj.B(sigma_fu, sigma_fu));

        Cd case_9  = Cd(0);
        Cd case_10 = case_2;
        Cd case_11 = case_7;
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dH_Mijab_a(const PairCache& c, const VectorXi& sigma, int a, int b) {
    int N = c.N;
    DerivTable result;

    Cd H = H_Mijab_from_cache(c, a, b);
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);
    Cd h_eff = sigma_gauss * sigma_gauss + h;

    Cd factor_h = H * (-0.5 / h_eff + p * p / (h_eff * h_eff));
    Cd factor_p = H * (-2.0 * p / h_eff);

    for (int m = 0; m < N; m++) {
        // Diagonal (m, m)
        Cd dh_diag = (-c.K_inv(a, m) * c.K_inv(m, a)
                      - c.K_inv(b, m) * c.K_inv(m, b)
                      + 2.0 * c.K_inv(a, m) * c.K_inv(m, b));
        Cd dp_diag = (-c.K_inv(a, m) * c.mu(m)
                      + c.K_inv(b, m) * c.mu(m));

        // case_6 (index[1])
        Cd case_6 = factor_h * dh_diag + factor_p * dp_diag;
        // case_15 = same (index[6])
        Cd case_15 = case_6;

        Cd case_19 = Cd(0);
        Cd case_20 = case_6;
        Cd case_23 = case_15;
        Cd case_24 = Cd(0);

        result.push_back({Cd(0), case_6, Cd(0), Cd(0),
                          Cd(0), Cd(0),  case_15, Cd(0),
                          case_19, case_20, case_23, case_24});

        // Off-diagonal (m, n), n > m
        for (int n = m + 1; n < N; n++) {
            Cd dh_off = (-c.K_inv(a, m) * c.K_inv(n, a)
                         - c.K_inv(a, n) * c.K_inv(m, a)
                         - c.K_inv(b, m) * c.K_inv(n, b)
                         - c.K_inv(b, n) * c.K_inv(m, b)
                         + 2.0 * c.K_inv(a, m) * c.K_inv(n, b)
                         + 2.0 * c.K_inv(a, n) * c.K_inv(m, b));
            Cd dp_off = (-c.K_inv(a, m) * c.mu(n)
                         - c.K_inv(a, n) * c.mu(m)
                         + c.K_inv(b, m) * c.mu(n)
                         + c.K_inv(b, n) * c.mu(m));

            // case_2 (index[1])
            Cd case_2 = factor_h * dh_off + factor_p * dp_off;
            // case_11 = same (index[6])
            Cd case_11 = case_2;

            Cd case_17 = Cd(0);
            Cd case_18 = case_2;
            Cd case_21 = case_11;
            Cd case_22 = Cd(0);

            result.push_back({Cd(0), case_2, Cd(0), Cd(0),
                               Cd(0), Cd(0),  case_11, Cd(0),
                               case_17, case_18, case_21, case_22});
        }
    }
    return result;
}

// ============================================================
// 5. Q_Mija (kicking term kernel) derivatives
// ============================================================

DerivTable dQ_Mija_b(const PairCache& c, const VectorXi& sigma, int a,
                     const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    Cd le = c.K_inv(a, a);
    Cd lc = 2.0 * c.mu(a);
    double k = k_L;

    Cd exp_factor = std::exp(-k * k * le);
    Cd cos_factor = std::cos(k * lc);
    Cd sin_factor = std::sin(k * lc);

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        // d(le)/dB_bra[i] = -K_inv[a][i]*K_inv[i][a]
        Cd d_le_bra = -c.K_inv(a, i) * c.K_inv(i, a);
        // d(lc)/dB_bra[i] = 2*K_inv[a][i]*(R_i_con[i] - mu[i])
        Cd d_lc_bra = 2.0 * c.K_inv(a, i) * (pi_con.R(i) - c.mu(i));

        // case_2: bra Real=False
        Cd case_2 = exp_factor * (-k * k * cos_factor * d_le_bra
                                   - k * sin_factor * d_lc_bra);

        // d(le)/dB_ket = same as bra (same K structure)
        Cd d_le_ket = -c.K_inv(a, i) * c.K_inv(i, a);
        // d(lc)/dB_ket[sigma_fu] = 2*K_inv[a][i]*(R_j[sigma_fu] - mu[i])
        Cd d_lc_ket = 2.0 * c.K_inv(a, i) * (pj.R(sigma_fu) - c.mu(i));

        // case_7: ket permuted (Real=True)
        Cd case_7 = exp_factor * (-k * k * cos_factor * d_le_ket
                                   - k * sin_factor * d_lc_ket);

        Cd case_9  = Cd(0);
        Cd case_10 = case_2;
        Cd case_11 = case_7;
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dQ_Mija_r(const PairCache& c, const VectorXi& sigma, int a,
                     const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    Cd le = c.K_inv(a, a);
    Cd lc = 2.0 * c.mu(a);
    double k = k_L;

    Cd exp_factor = std::exp(-k * k * le);
    Cd sin_factor = std::sin(k * lc);

    for (int i = 0; i < N; i++) {
        int sigma_fu = sigma_inv(sigma, i);

        // d(lc)/dR_i_con[i] = 2*K_inv[a][i]*B_i_con[i][i]
        Cd d_lc_bra = 2.0 * c.K_inv(a, i) * pi_con.B(i, i);

        // case_2: bra Real=False
        Cd case_2 = exp_factor * (-k * sin_factor * d_lc_bra);

        // d(lc)/dR_j[sigma_fu] = 2*K_inv[a][i]*B_j[sigma_fu][sigma_fu]
        Cd d_lc_ket = 2.0 * c.K_inv(a, i) * pj.B(sigma_fu, sigma_fu);

        // case_7: ket permuted (Real=True)
        Cd case_7 = exp_factor * (-k * sin_factor * d_lc_ket);

        Cd case_9  = Cd(0);
        Cd case_10 = case_2;
        Cd case_11 = case_7;
        Cd case_12 = Cd(0);

        result[i] = {Cd(0), case_2, Cd(0), Cd(0),
                     Cd(0), Cd(0),  case_7, Cd(0),
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dQ_Mija_a(const PairCache& c, const VectorXi& sigma, int a) {
    int N = c.N;
    DerivTable result;

    Cd le = c.K_inv(a, a);
    Cd lc = 2.0 * c.mu(a);
    double k = k_L;

    Cd exp_factor = std::exp(-k * k * le);
    Cd cos_factor = std::cos(k * lc);
    Cd sin_factor = std::sin(k * lc);

    for (int m = 0; m < N; m++) {
        // Diagonal (m, m)
        // d(le)/dA[m][m] = -K_inv[a][m]*K_inv[m][a]
        Cd d_le_diag = -c.K_inv(a, m) * c.K_inv(m, a);
        // d(lc)/dA[m][m] = -2*K_inv[a][m]*mu[m]
        Cd d_lc_diag = -2.0 * c.K_inv(a, m) * c.mu(m);

        // case_6 (index[1])
        Cd case_6 = exp_factor * (-k * k * cos_factor * d_le_diag
                                   - k * sin_factor * d_lc_diag);
        // case_15 = same (index[6])
        Cd case_15 = case_6;

        Cd case_19 = Cd(0);
        Cd case_20 = case_6;
        Cd case_23 = case_15;
        Cd case_24 = Cd(0);

        result.push_back({Cd(0), case_6, Cd(0), Cd(0),
                          Cd(0), Cd(0),  case_15, Cd(0),
                          case_19, case_20, case_23, case_24});

        // Off-diagonal (m, n), n > m
        for (int n = m + 1; n < N; n++) {
            // d(le)/dA[m][n] = -(K_inv[a][m]*K_inv[n][a] + K_inv[a][n]*K_inv[m][a])
            Cd d_le_off = -(c.K_inv(a, m) * c.K_inv(n, a)
                            + c.K_inv(a, n) * c.K_inv(m, a));
            // d(lc)/dA[m][n] = -2*(K_inv[a][m]*mu[n] + K_inv[a][n]*mu[m])
            Cd d_lc_off = -2.0 * (c.K_inv(a, m) * c.mu(n)
                                   + c.K_inv(a, n) * c.mu(m));

            // case_2 (index[1])
            Cd case_2 = exp_factor * (-k * k * cos_factor * d_le_off
                                       - k * sin_factor * d_lc_off);
            // case_11 = same (index[6])
            Cd case_11 = case_2;

            Cd case_17 = Cd(0);
            Cd case_18 = case_2;
            Cd case_21 = case_11;
            Cd case_22 = Cd(0);

            result.push_back({Cd(0), case_2, Cd(0), Cd(0),
                               Cd(0), Cd(0),  case_11, Cd(0),
                               case_17, case_18, case_21, case_22});
        }
    }
    return result;
}

} // namespace ecg1d
