#include "derivatives.hpp"
#include "hamiltonian.hpp"
#include <stdexcept>

namespace ecg1d {

int upper_index_inclusive(int N, int i, int j) {
    int S_i = i * (2 * N - i + 1) / 2;
    return S_i + (j - i);
}

static Case12 zeros12() {
    Case12 c{};
    for (auto& v : c) v = Cd(0, 0);
    return c;
}

// ============================================================
// Log-overlap derivatives
// ============================================================

DerivTable dln_m_g_b(const PairCache& c, const VectorXi& sigma,
                     const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    for (int i = 0; i < N; i++) {
        int sigma_zheng = sigma(i);
        int sigma_fu = -1;
        for (int z = 0; z < N; z++) if (sigma(z) == i) { sigma_fu = z; break; }

        Cd NINA_i = (c.K_inv.row(i) * c.b)(0);
        // NINA_sigma_zheng not used directly but kept for clarity

        Cd case_1(0), case_2(0), case_3(0), case_4(0);
        Cd case_5(0), case_6(0), case_7(0), case_8(0);

        // case_2: bra direct Real
        case_2 = -0.5 * c.K_inv(i, i) - pi_con.R(i) * pi_con.R(i)
                 - 0.25 * NINA_i * NINA_i + NINA_i * pi_con.R(i);

        // case_3: ket direct Real
        case_3 = -0.5 * c.K_inv(sigma_zheng, sigma_zheng) - pj.R(i) * pj.R(i);

        // case_7: ket permuted Real
        case_7 = -0.25 * NINA_i * NINA_i + NINA_i * pj.R(sigma_fu);

        // Combined cases
        Cd case_9  = case_1 + case_3;
        Cd case_10 = case_2 + case_4;
        Cd case_11 = case_5 + case_7;
        Cd case_12 = case_6 + case_8;

        result[i] = {case_1, case_2, case_3, case_4,
                     case_5, case_6, case_7, case_8,
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dln_m_g_r(const PairCache& c, const VectorXi& sigma,
                     const BasisParams& pi, const BasisParams& pj) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    for (int i = 0; i < N; i++) {
        int sigma_fu = -1;
        for (int z = 0; z < N; z++) if (sigma(z) == i) { sigma_fu = z; break; }

        Cd NINA_i = (c.K_inv.row(i) * c.b)(0);

        Cd case_1(0), case_2(0), case_3(0), case_4(0);
        Cd case_5(0), case_6(0), case_7(0), case_8(0);

        case_2 = -2.0 * pi_con.R(i) * pi_con.B(i, i) + NINA_i * pi_con.B(i, i);
        case_3 = -2.0 * pj.R(i) * pj.B(i, i);
        case_7 = NINA_i * pj.B(sigma_fu, sigma_fu);

        Cd case_9  = case_1 + case_3;
        Cd case_10 = case_2 + case_4;
        Cd case_11 = case_5 + case_7;
        Cd case_12 = case_6 + case_8;

        result[i] = {case_1, case_2, case_3, case_4,
                     case_5, case_6, case_7, case_8,
                     case_9, case_10, case_11, case_12};
    }
    return result;
}

DerivTable dln_m_g_a(const PairCache& c, const VectorXi& sigma) {
    int N = c.N;
    DerivTable result;

    for (int m = 0; m < N; m++) {
        int sigma_zheng_m = sigma(m);
        Cd NINA_m = (c.K_inv.row(m) * c.b)(0);

        // Diagonal entry (m,m)
        Cd c5(0), c6(0), c7(0), c8(0);
        Cd c13(0), c14(0), c15(0), c16(0);

        c6 = -0.5 * c.K_inv(m, m) - 0.25 * NINA_m * NINA_m;
        c7 = -0.5 * c.K_inv(sigma_zheng_m, sigma_zheng_m);
        c15 = -0.25 * NINA_m * NINA_m;

        Cd c19 = c5 + c7;
        Cd c20 = c6 + c8;
        Cd c23 = c13 + c15;
        Cd c24 = c14 + c16;

        result.push_back({c5, c6, c7, c8, c13, c14, c15, c16, c19, c20, c23, c24});

        // Off-diagonal entries (m, n) for n > m
        for (int n = m + 1; n < N; n++) {
            int sigma_zheng_n = sigma(n);
            Cd NINA_n = (c.K_inv.row(n) * c.b)(0);

            Cd o1(0), o2(0), o3(0), o4(0);
            Cd o9(0), o10(0), o11(0), o12(0);

            o2 = -1.0 * c.K_inv(n, m) - 0.5 * NINA_n * NINA_m;
            o3 = -1.0 * c.K_inv(sigma_zheng_n, sigma_zheng_m);
            o11 = -0.5 * NINA_n * NINA_m;

            Cd o17 = o1 + o3;
            Cd o18 = o2 + o4;
            Cd o21 = o9 + o11;
            Cd o22 = o10 + o12;

            result.push_back({o1, o2, o3, o4, o9, o10, o11, o12, o17, o18, o21, o22});
        }
    }
    return result;
}

// ============================================================
// K-inverse derivatives
// ============================================================

DerivTable dK_fu_B(int t, int o, const PairCache& c) {
    int N = c.N;
    DerivTable result(N);
    for (int i = 0; i < N; i++) {
        Case12 cs = zeros12();
        cs[1] = -1.0 * c.K_inv(t, i) * c.K_inv(i, o);  // case_2
        cs[9] = cs[1];  // case_10 = case_2 + case_4(=0)
        result[i] = cs;
    }
    return result;
}

DerivTable dK_fu_R(int N) {
    DerivTable result(N);
    for (int i = 0; i < N; i++) result[i] = zeros12();
    return result;
}

DerivTable dK_fu_A(int t, int o, const PairCache& c) {
    int N = c.N;
    DerivTable result;

    for (int m = 0; m < N; m++) {
        // Diagonal (m,m)
        Case12 cs = zeros12();
        cs[1] = -1.0 * c.K_inv(t, m) * c.K_inv(m, o);  // case_6
        cs[9] = cs[1];  // case_20 = case_6 + case_8(=0)
        result.push_back(cs);

        // Off-diagonal (m,n)
        for (int n = m + 1; n < N; n++) {
            Case12 os = zeros12();
            os[1] = -c.K_inv(t, m) * c.K_inv(n, o) - c.K_inv(t, n) * c.K_inv(m, o);  // case_2
            os[9] = os[1];  // case_18 = case_2 + case_4(=0)
            result.push_back(os);
        }
    }
    return result;
}

// ============================================================
// (K_inv @ b) derivatives
// ============================================================

DerivTable dK_fu_b_B(int t, const PairCache& c, const VectorXi& sigma,
                     const BasisParams& pi) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    for (int i = 0; i < N; i++) {
        Cd NINA_i = (c.K_inv.row(i) * c.b)(0);
        Case12 cs = zeros12();
        // case_2: 2*K_inv[t,i]*conj(R_i[i]) - K_inv[t,i]*NINA_i
        cs[1] = c.K_inv(t, i) * pi_con.R(i) + c.K_inv(t, i) * pi_con.R(i)
                - c.K_inv(t, i) * NINA_i;
        cs[9] = cs[1];  // case_10 = case_2
        result[i] = cs;
    }
    return result;
}

DerivTable dK_fu_b_R(int t, const PairCache& c,
                     const BasisParams& pi) {
    int N = c.N;
    BasisParams pi_con = pi.conj_params();
    DerivTable result(N);

    for (int i = 0; i < N; i++) {
        Case12 cs = zeros12();
        cs[1] = 2.0 * c.K_inv(t, i) * pi_con.B(i, i);  // case_2
        cs[9] = cs[1];  // case_10
        result[i] = cs;
    }
    return result;
}

DerivTable dK_fu_b_A(int t, const PairCache& c) {
    int N = c.N;
    DerivTable result;

    for (int m = 0; m < N; m++) {
        Cd NINA_m = (c.K_inv.row(m) * c.b)(0);
        Case12 cs = zeros12();
        cs[1] = -1.0 * c.K_inv(t, m) * NINA_m;  // case_6
        cs[9] = cs[1];  // case_20
        result.push_back(cs);

        for (int n = m + 1; n < N; n++) {
            Cd NINA_n = (c.K_inv.row(n) * c.b)(0);
            Case12 os = zeros12();
            os[1] = -c.K_inv(t, m) * NINA_n - c.K_inv(t, n) * NINA_m;  // case_2
            os[9] = os[1];  // case_18
            result.push_back(os);
        }
    }
    return result;
}

// ============================================================
// Case index resolution
// ============================================================

std::pair<int,int> resolve_case_first(int alpha_2, int name_i, int name_j, bool Real) {
    bool is_i = (alpha_2 == name_i);
    bool is_j = (alpha_2 == name_j);

    if (is_i && !is_j && Real)  return {0, 4};
    if (!is_i && is_j && Real)  return {2, 6};
    if (is_i && is_j && Real)   return {8, 10};
    if (is_i && !is_j && !Real) return {1, 5};
    if (!is_i && is_j && !Real) return {3, 7};
    if (is_i && is_j && !Real)  return {9, 11};

    return {0, 0}; // should not reach
}

std::pair<int,int> resolve_case_second(int beta_2, int name_i, int name_j) {
    bool is_i = (beta_2 == name_i);
    bool is_j = (beta_2 == name_j);

    if (is_i && !is_j) return {1, 5};
    if (!is_i && is_j) return {3, 7};
    if (is_i && is_j)  return {9, 11};

    return {0, 0};
}

// ============================================================
// partial_z_lnmg
// ============================================================

std::vector<Cd> partial_z_lnmg(int alpha_1, bool Real,
                                const BasisParams& pi, const BasisParams& pj,
                                int alpha_2, int alpha_3, int alpha_4,
                                const PermutationSet& perms) {
    int SN = perms.SN;
    int N = perms.N;

    if (alpha_2 != pi.name && alpha_2 != pj.name) {
        return std::vector<Cd>(SN, Cd(0, 0));
    }

    auto [idx1, idx2] = resolve_case_first(alpha_2, pi.name, pj.name, Real);

    if (alpha_1 == 2) {
        std::vector<Cd> result(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            DerivTable dt = dln_m_g_b(c, perms.sigmas[p], pi, pj);
            Cd val = dt[alpha_3][idx1];
            int a3_special = perms.sigmas[p](alpha_3);
            val += dt[a3_special][idx2];
            result[p] = val;
        }
        return result;
    }

    if (alpha_1 == 3) {
        std::vector<Cd> result(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            DerivTable dt = dln_m_g_r(c, perms.sigmas[p], pi, pj);
            Cd val = dt[alpha_3][idx1];
            int a3_special = perms.sigmas[p](alpha_3);
            val += dt[a3_special][idx2];
            result[p] = val;
        }
        return result;
    }

    if (alpha_1 == 4) {
        int index = upper_index_inclusive(N, alpha_3, alpha_4);
        std::vector<Cd> result(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            DerivTable dt = dln_m_g_a(c, perms.sigmas[p]);
            Cd val = dt[index][idx1];
            int a3s = perms.sigmas[p](alpha_3);
            int a4s = perms.sigmas[p](alpha_4);
            if (a4s < a3s) std::swap(a3s, a4s);
            int index_special = upper_index_inclusive(N, a3s, a4s);
            val += dt[index_special][idx2];
            result[p] = val;
        }
        return result;
    }

    return std::vector<Cd>(SN, Cd(0, 0));
}

// ============================================================
// partial_z_lnmg_second - second derivative of log-overlap
// ============================================================

// Helper: dsecond_b
static std::vector<Cd> dsecond_b_impl(int alpha_3, int /*alpha_4*/,
                                       const BasisParams& pi, const BasisParams& pj,
                                       int beta_1, int beta_2, int beta_3, int beta_4,
                                       const PermutationSet& perms) {
    int SN = perms.SN;
    int N = perms.N;

    if (beta_2 != pi.name && beta_2 != pj.name)
        return std::vector<Cd>(SN, Cd(0, 0));

    auto [idx1, idx2] = resolve_case_second(beta_2, pi.name, pj.name);

    std::vector<Cd> result(SN);
    for (int p = 0; p < SN; p++) {
        PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
        int a3s = perms.sigmas[p](alpha_3);
        int b3s = perms.sigmas[p](beta_3);
        Cd NINA_sigma = (c.K_inv.row(a3s) * c.b)(0);

        Cd val(0, 0);

        if (beta_1 == 2) {
            DerivTable dk = dK_fu_B(a3s, a3s, c);
            DerivTable dkb = dK_fu_b_B(a3s, c, perms.sigmas[p], pi);
            val += -0.5 * dk[beta_3][idx1];
            val += -0.5 * dk[b3s][idx2];
            val += -0.5 * NINA_sigma * dkb[beta_3][idx1];
            val += -0.5 * NINA_sigma * dkb[b3s][idx2];
            val += pj.R(alpha_3) * dkb[beta_3][idx1];
            val += pj.R(alpha_3) * dkb[b3s][idx2];
        } else if (beta_1 == 3) {
            DerivTable dk = dK_fu_R(N);
            DerivTable dkb = dK_fu_b_R(a3s, c, pi);
            val += -0.5 * dk[beta_3][idx1];
            val += -0.5 * dk[b3s][idx2];
            val += -0.5 * NINA_sigma * dkb[beta_3][idx1];
            val += -0.5 * NINA_sigma * dkb[b3s][idx2];
            val += pj.R(alpha_3) * dkb[beta_3][idx1];
            val += pj.R(alpha_3) * dkb[b3s][idx2];
        } else if (beta_1 == 4) {
            int b_index = upper_index_inclusive(N, beta_3, beta_4);
            int b4s = perms.sigmas[p](beta_4);
            if (b4s < b3s) std::swap(b3s, b4s);
            int b_index_s = upper_index_inclusive(N, b3s, b4s);

            DerivTable dk = dK_fu_A(a3s, a3s, c);
            DerivTable dkb = dK_fu_b_A(a3s, c);
            val += -0.5 * dk[b_index][idx1];
            val += -0.5 * dk[b_index_s][idx2];
            val += -0.5 * NINA_sigma * dkb[b_index][idx1];
            val += -0.5 * NINA_sigma * dkb[b_index_s][idx2];
            val += pj.R(alpha_3) * dkb[b_index][idx1];
            val += pj.R(alpha_3) * dkb[b_index_s][idx2];
        }
        result[p] = val;
    }
    return result;
}

// Helper: dsecond_r_1 - always zero (bra R with Real=True)
static std::vector<Cd> dsecond_r_1_impl(int SN) {
    return std::vector<Cd>(SN, Cd(0, 0));
}

// Helper: dsecond_r_2 - ket-side R second derivative
static std::vector<Cd> dsecond_r_2_impl(int alpha_3, int /*alpha_4*/,
                                         const BasisParams& pi, const BasisParams& pj,
                                         int beta_1, int beta_2, int beta_3, int beta_4,
                                         const PermutationSet& perms) {
    int SN = perms.SN;
    int N = perms.N;

    if (beta_2 != pi.name && beta_2 != pj.name)
        return std::vector<Cd>(SN, Cd(0, 0));

    auto [idx1, idx2] = resolve_case_second(beta_2, pi.name, pj.name);

    std::vector<Cd> result(SN);
    for (int p = 0; p < SN; p++) {
        PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
        int a3s = perms.sigmas[p](alpha_3);
        int b3s = perms.sigmas[p](beta_3);

        Cd val(0, 0);

        if (beta_1 == 2) {
            DerivTable dkb = dK_fu_b_B(a3s, c, perms.sigmas[p], pi);
            val += pj.B(alpha_3, alpha_3) * dkb[beta_3][idx1];
            val += pj.B(alpha_3, alpha_3) * dkb[b3s][idx2];
        } else if (beta_1 == 3) {
            DerivTable dkb = dK_fu_b_R(a3s, c, pi);
            val += pj.B(alpha_3, alpha_3) * dkb[beta_3][idx1];
            val += pj.B(alpha_3, alpha_3) * dkb[b3s][idx2];
        } else if (beta_1 == 4) {
            int b_index = upper_index_inclusive(N, beta_3, beta_4);
            int b4s = perms.sigmas[p](beta_4);
            if (b4s < b3s) std::swap(b3s, b4s);
            int b_index_s = upper_index_inclusive(N, b3s, b4s);

            DerivTable dkb = dK_fu_b_A(a3s, c);
            val += pj.B(alpha_3, alpha_3) * dkb[b_index][idx1];
            val += pj.B(alpha_3, alpha_3) * dkb[b_index_s][idx2];
        }
        result[p] = val;
    }
    return result;
}

// Helper: dsecond_a_3 (off-diagonal A, alpha_3 != alpha_4)
static std::vector<Cd> dsecond_a_3_impl(int alpha_3, int alpha_4,
                                         const BasisParams& pi, const BasisParams& pj,
                                         int beta_1, int beta_2, int beta_3, int beta_4,
                                         const PermutationSet& perms) {
    int SN = perms.SN;
    int N = perms.N;

    if (beta_2 != pi.name && beta_2 != pj.name)
        return std::vector<Cd>(SN, Cd(0, 0));

    auto [idx1, idx2] = resolve_case_second(beta_2, pi.name, pj.name);

    std::vector<Cd> result(SN);
    for (int p = 0; p < SN; p++) {
        PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
        int a3s = perms.sigmas[p](alpha_3);
        int a4s = perms.sigmas[p](alpha_4);
        int b3s = perms.sigmas[p](beta_3);

        Cd NINA_sigma = (c.K_inv.row(a3s) * c.b)(0);
        Cd NINA4_sigma = (c.K_inv.row(a4s) * c.b)(0);

        Cd val(0, 0);

        if (beta_1 == 2) {
            int b3s_perm = perms.sigmas[p](beta_3);
            DerivTable dk = dK_fu_B(a4s, a3s, c);
            DerivTable dkb3 = dK_fu_b_B(a3s, c, perms.sigmas[p], pi);
            DerivTable dkb4 = dK_fu_b_B(a4s, c, perms.sigmas[p], pi);
            val += -1.0 * dk[beta_3][idx1];
            val += -1.0 * dk[b3s_perm][idx2];
            val += -0.5 * NINA4_sigma * dkb3[beta_3][idx1];
            val += -0.5 * NINA4_sigma * dkb3[b3s_perm][idx2];
            val += -0.5 * NINA_sigma * dkb4[beta_3][idx1];
            val += -0.5 * NINA_sigma * dkb4[b3s_perm][idx2];
        } else if (beta_1 == 3) {
            int b3s_perm = perms.sigmas[p](beta_3);
            DerivTable dk = dK_fu_R(N);
            DerivTable dkb3 = dK_fu_b_R(a3s, c, pi);
            DerivTable dkb4 = dK_fu_b_R(a4s, c, pi);
            val += -1.0 * dk[beta_3][idx1];
            val += -1.0 * dk[b3s_perm][idx2];
            val += -0.5 * NINA4_sigma * dkb3[beta_3][idx1];
            val += -0.5 * NINA4_sigma * dkb3[b3s_perm][idx2];
            val += -0.5 * NINA_sigma * dkb4[beta_3][idx1];
            val += -0.5 * NINA_sigma * dkb4[b3s_perm][idx2];
        } else if (beta_1 == 4) {
            int b_index = upper_index_inclusive(N, beta_3, beta_4);
            int b3s_perm = perms.sigmas[p](beta_3);
            int b4s_perm = perms.sigmas[p](beta_4);
            if (b4s_perm < b3s_perm) std::swap(b3s_perm, b4s_perm);
            int b_index_s = upper_index_inclusive(N, b3s_perm, b4s_perm);

            DerivTable dk = dK_fu_A(a4s, a3s, c);
            DerivTable dkb3 = dK_fu_b_A(a3s, c);
            DerivTable dkb4 = dK_fu_b_A(a4s, c);
            val += -1.0 * dk[b_index][idx1];
            val += -1.0 * dk[b_index_s][idx2];
            val += -0.5 * NINA4_sigma * dkb3[b_index][idx1];
            val += -0.5 * NINA4_sigma * dkb3[b_index_s][idx2];
            val += -0.5 * NINA_sigma * dkb4[b_index][idx1];
            val += -0.5 * NINA_sigma * dkb4[b_index_s][idx2];
        }
        result[p] = val;
    }
    return result;
}

// Helper: dsecond_a_7 (diagonal A, alpha_3 == alpha_4)
static std::vector<Cd> dsecond_a_7_impl(int alpha_3, int alpha_4,
                                         const BasisParams& pi, const BasisParams& pj,
                                         int beta_1, int beta_2, int beta_3, int beta_4,
                                         const PermutationSet& perms) {
    int SN = perms.SN;
    int N = perms.N;

    if (beta_2 != pi.name && beta_2 != pj.name)
        return std::vector<Cd>(SN, Cd(0, 0));

    auto [idx1, idx2] = resolve_case_second(beta_2, pi.name, pj.name);

    std::vector<Cd> result(SN);
    for (int p = 0; p < SN; p++) {
        PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
        int a3s = perms.sigmas[p](alpha_3);
        int a4s = perms.sigmas[p](alpha_4);
        int b3s = perms.sigmas[p](beta_3);

        Cd NINA_sigma = (c.K_inv.row(a3s) * c.b)(0);
        Cd NINA4_sigma = (c.K_inv.row(a4s) * c.b)(0);

        Cd val(0, 0);

        if (beta_1 == 2) {
            int b3s_perm = perms.sigmas[p](beta_3);
            DerivTable dk = dK_fu_B(a4s, a3s, c);
            DerivTable dkb3 = dK_fu_b_B(a3s, c, perms.sigmas[p], pi);
            DerivTable dkb4 = dK_fu_b_B(a4s, c, perms.sigmas[p], pi);
            val += -0.5 * dk[beta_3][idx1];
            val += -0.5 * dk[b3s_perm][idx2];
            val += -0.25 * NINA4_sigma * dkb3[beta_3][idx1];
            val += -0.25 * NINA4_sigma * dkb3[b3s_perm][idx2];
            val += -0.25 * NINA_sigma * dkb4[beta_3][idx1];
            val += -0.25 * NINA_sigma * dkb4[b3s_perm][idx2];
        } else if (beta_1 == 3) {
            int b3s_perm = perms.sigmas[p](beta_3);
            DerivTable dk = dK_fu_R(N);
            DerivTable dkb3 = dK_fu_b_R(a3s, c, pi);
            DerivTable dkb4 = dK_fu_b_R(a4s, c, pi);
            val += -0.5 * dk[beta_3][idx1];
            val += -0.5 * dk[b3s_perm][idx2];
            val += -0.25 * NINA4_sigma * dkb3[beta_3][idx1];
            val += -0.25 * NINA4_sigma * dkb3[b3s_perm][idx2];
            val += -0.25 * NINA_sigma * dkb4[beta_3][idx1];
            val += -0.25 * NINA_sigma * dkb4[b3s_perm][idx2];
        } else if (beta_1 == 4) {
            int b_index = upper_index_inclusive(N, beta_3, beta_4);
            int b3s_perm = perms.sigmas[p](beta_3);
            int b4s_perm = perms.sigmas[p](beta_4);
            if (b4s_perm < b3s_perm) std::swap(b3s_perm, b4s_perm);
            int b_index_s = upper_index_inclusive(N, b3s_perm, b4s_perm);

            DerivTable dk = dK_fu_A(a4s, a3s, c);
            DerivTable dkb3 = dK_fu_b_A(a3s, c);
            DerivTable dkb4 = dK_fu_b_A(a4s, c);
            val += -0.5 * dk[b_index][idx1];
            val += -0.5 * dk[b_index_s][idx2];
            val += -0.25 * NINA4_sigma * dkb3[b_index][idx1];
            val += -0.25 * NINA4_sigma * dkb3[b_index_s][idx2];
            val += -0.25 * NINA_sigma * dkb4[b_index][idx1];
            val += -0.25 * NINA_sigma * dkb4[b_index_s][idx2];
        }
        result[p] = val;
    }
    return result;
}

std::vector<Cd> partial_z_lnmg_second(int alpha_1, int beta_1,
                                       const BasisParams& pi, const BasisParams& pj,
                                       int alpha_2, int alpha_3, int alpha_4,
                                       int beta_2, int beta_3, int beta_4,
                                       const PermutationSet& perms) {
    int SN = perms.SN;

    if (alpha_2 != pi.name && alpha_2 != pj.name)
        return std::vector<Cd>(SN, Cd(0, 0));
    if (beta_2 != pi.name && beta_2 != pj.name)
        return std::vector<Cd>(SN, Cd(0, 0));

    if (alpha_1 == 2) {
        if (alpha_2 != pj.name)
            return std::vector<Cd>(SN, Cd(0, 0));
        return dsecond_b_impl(alpha_3, alpha_4, pi, pj, beta_1, beta_2, beta_3, beta_4, perms);
    }

    if (alpha_1 == 3) {
        if (alpha_2 == pi.name && alpha_2 != pj.name) {
            return dsecond_r_1_impl(SN);
        }
        if (alpha_2 != pi.name && alpha_2 == pj.name) {
            return dsecond_r_2_impl(alpha_3, alpha_4, pi, pj, beta_1, beta_2, beta_3, beta_4, perms);
        }
        if (alpha_2 == pi.name && alpha_2 == pj.name) {
            auto r1 = dsecond_r_1_impl(SN);
            auto r2 = dsecond_r_2_impl(alpha_3, alpha_4, pi, pj, beta_1, beta_2, beta_3, beta_4, perms);
            for (int p = 0; p < SN; p++) r1[p] += r2[p];
            return r1;
        }
    }

    if (alpha_1 == 4) {
        if (alpha_2 != pj.name)
            return std::vector<Cd>(SN, Cd(0, 0));
        if (alpha_3 != alpha_4)
            return dsecond_a_3_impl(alpha_3, alpha_4, pi, pj, beta_1, beta_2, beta_3, beta_4, perms);
        else
            return dsecond_a_7_impl(alpha_3, alpha_4, pi, pj, beta_1, beta_2, beta_3, beta_4, perms);
    }

    return std::vector<Cd>(SN, Cd(0, 0));
}

// ============================================================
// partial_z_addsn
// ============================================================

Cd partial_z_addsn(int alpha_1, bool Real,
                   const BasisParams& pi, const BasisParams& pj,
                   int alpha_2, int alpha_3, int alpha_4,
                   const PermutationSet& perms) {
    int SN = perms.SN;

    // Compute M_G for all permutations
    std::vector<Cd> M_G(SN);
    for (int p = 0; p < SN; p++) {
        PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
        M_G[p] = c.M_G;
    }

    if (alpha_1 == 1) {
        Cd result(0, 0);
        for (int p = 0; p < SN; p++)
            result += static_cast<double>(perms.signs[p]) * M_G[p];
        return result;
    }

    // alpha_1 == 2, 3, 4
    std::vector<Cd> LNMGB = partial_z_lnmg(alpha_1, Real, pi, pj,
                                             alpha_2, alpha_3, alpha_4, perms);
    Cd result(0, 0);
    for (int p = 0; p < SN; p++)
        result += static_cast<double>(perms.signs[p]) * M_G[p] * LNMGB[p];
    return result;
}

// ============================================================
// partial_z_addsn_second
// ============================================================

Cd partial_z_addsn_second(int alpha_1, int beta_1,
                          const BasisParams& pi, const BasisParams& pj,
                          int alpha_2, int alpha_3, int alpha_4,
                          int beta_2, int beta_3, int beta_4,
                          const PermutationSet& perms) {
    int SN = perms.SN;

    std::vector<Cd> M_G(SN);
    for (int p = 0; p < SN; p++) {
        PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
        M_G[p] = c.M_G;
    }

    if (alpha_1 == 1 && beta_1 == 1) {
        Cd result(0, 0);
        for (int p = 0; p < SN; p++)
            result += static_cast<double>(perms.signs[p]) * M_G[p];
        return result;
    }

    if (alpha_1 == 1 && beta_1 != 1) {
        std::vector<Cd> LNMGB = partial_z_lnmg(beta_1, false, pi, pj,
                                                 beta_2, beta_3, beta_4, perms);
        Cd result(0, 0);
        for (int p = 0; p < SN; p++)
            result += static_cast<double>(perms.signs[p]) * M_G[p] * LNMGB[p];
        return result;
    }

    if (alpha_1 != 1 && beta_1 == 1) {
        std::vector<Cd> LNMGB = partial_z_lnmg(alpha_1, true, pi, pj,
                                                 alpha_2, alpha_3, alpha_4, perms);
        Cd result(0, 0);
        for (int p = 0; p < SN; p++)
            result += static_cast<double>(perms.signs[p]) * M_G[p] * LNMGB[p];
        return result;
    }

    // alpha_1 != 1 && beta_1 != 1
    std::vector<Cd> LNMGB_alpha = partial_z_lnmg(alpha_1, true, pi, pj,
                                                   alpha_2, alpha_3, alpha_4, perms);
    std::vector<Cd> LNMGB_beta = partial_z_lnmg(beta_1, false, pi, pj,
                                                  beta_2, beta_3, beta_4, perms);
    std::vector<Cd> LNMGB_ab = partial_z_lnmg_second(alpha_1, beta_1, pi, pj,
                                                       alpha_2, alpha_3, alpha_4,
                                                       beta_2, beta_3, beta_4, perms);

    Cd result(0, 0);
    for (int p = 0; p < SN; p++)
        result += static_cast<double>(perms.signs[p]) * M_G[p] *
                  (LNMGB_alpha[p] * LNMGB_beta[p] + LNMGB_ab[p]);
    return result;
}

// ============================================================
// partial_z_first
// ============================================================

Cd partial_z_first(int alpha_1, bool Real,
                   const std::vector<BasisParams>& basis,
                   int alpha_2, int alpha_3, int alpha_4) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    if (alpha_1 == 1 && Real) {
        Cd result(0, 0);
        for (int t = 0; t < basis_n; t++) {
            Cd con_ui = std::conj(basis[t].u);
            result += con_ui * partial_z_addsn(alpha_1, Real, basis[t], basis[alpha_2],
                                                alpha_2, alpha_3, alpha_4, perms);
        }
        return result;
    }

    if (alpha_1 == 1 && !Real) {
        Cd result(0, 0);
        for (int t = 0; t < basis_n; t++) {
            Cd uj = basis[t].u;
            result += uj * partial_z_addsn(alpha_1, Real, basis[alpha_2], basis[t],
                                            alpha_2, alpha_3, alpha_4, perms);
        }
        return result;
    }

    // alpha_1 != 1
    Cd result(0, 0);
    for (int i = 0; i < basis_n; i++) {
        Cd con_ui = std::conj(basis[i].u);
        for (int j = 0; j < basis_n; j++) {
            Cd uj = basis[j].u;
            result += con_ui * uj * partial_z_addsn(alpha_1, Real, basis[i], basis[j],
                                                     alpha_2, alpha_3, alpha_4, perms);
        }
    }
    return result;
}

// ============================================================
// partial_z_second
// ============================================================

Cd partial_z_second(int alpha_1, int beta_1,
                    const std::vector<BasisParams>& basis,
                    int alpha_2, int alpha_3, int alpha_4,
                    int beta_2, int beta_3, int beta_4) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);

    if (alpha_1 == 1 && beta_1 == 1) {
        return partial_z_addsn_second(alpha_1, beta_1, basis[beta_2], basis[alpha_2],
                                       alpha_2, alpha_3, alpha_4,
                                       beta_2, beta_3, beta_4, perms);
    }

    if (alpha_1 == 1 && beta_1 != 1) {
        Cd result(0, 0);
        for (int t = 0; t < basis_n; t++) {
            Cd con_ui = std::conj(basis[t].u);
            result += con_ui * partial_z_addsn_second(alpha_1, beta_1,
                                                       basis[t], basis[alpha_2],
                                                       alpha_2, alpha_3, alpha_4,
                                                       beta_2, beta_3, beta_4, perms);
        }
        return result;
    }

    if (alpha_1 != 1 && beta_1 == 1) {
        Cd result(0, 0);
        for (int t = 0; t < basis_n; t++) {
            Cd uj = basis[t].u;
            result += uj * partial_z_addsn_second(alpha_1, beta_1,
                                                   basis[beta_2], basis[t],
                                                   alpha_2, alpha_3, alpha_4,
                                                   beta_2, beta_3, beta_4, perms);
        }
        return result;
    }

    // alpha_1 != 1 && beta_1 != 1
    Cd result(0, 0);
    for (int i = 0; i < basis_n; i++) {
        Cd con_ui = std::conj(basis[i].u);
        for (int j = 0; j < basis_n; j++) {
            Cd uj = basis[j].u;
            result += con_ui * uj * partial_z_addsn_second(alpha_1, beta_1,
                                                            basis[i], basis[j],
                                                            alpha_2, alpha_3, alpha_4,
                                                            beta_2, beta_3, beta_4, perms);
        }
    }
    return result;
}

// ============================================================
// calculate_C
// ============================================================

Cd calculate_C(int alpha_1, int alpha_2, int alpha_3, int alpha_4,
               int beta_1, int beta_2, int beta_3, int beta_4,
               const std::vector<BasisParams>& basis) {
    Cd first_alpha = partial_z_first(alpha_1, true, basis, alpha_2, alpha_3, alpha_4);
    Cd first_beta = partial_z_first(beta_1, false, basis, beta_2, beta_3, beta_4);
    Cd S = overlap(basis);
    Cd second = partial_z_second(alpha_1, beta_1, basis,
                                  alpha_2, alpha_3, alpha_4,
                                  beta_2, beta_3, beta_4);

    return -1.0 / (S * S) * first_alpha * first_beta + second / S;
}

} // namespace ecg1d
