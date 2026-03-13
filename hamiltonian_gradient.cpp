#include "hamiltonian_gradient.hpp"
#include "hamiltonian.hpp"
#include "derivatives.hpp"
#include "observable_derivatives.hpp"
#include "pair_cache.hpp"
#include "permutation.hpp"
#include "interaction_kernels.hpp"
#include "physical_constants.hpp"

namespace ecg1d {

// ============================================================
// Generic observable gradient pattern
//
// Each Hamiltonian term has the form:
//   <H> = sum_{i,j} u_i* u_j * sum_p sign[p] * M_G[p] * kernel[p]
//
// d(<H>/S)/dz = (1/S)*d<H>/dz - (<H>/S²)*dS/dz
//
// d<H>/dz = d<H>_dMG + d<H>_dKernel
// where:
//   d<H>_dMG = sum over basis pairs * sum_p sign[p]*M_G[p]*kernel[p]*dln_MG[p]
//   d<H>_dKernel = sum over basis pairs * sum_p sign[p]*M_G[p]*d(kernel)[p]
//
// For alpha_1==1 (u derivative), dKernel = 0 (kernel doesn't depend on u)
// ============================================================

// Helper: compute per-perm kernel derivative dispatch (partial_z_kernel)
// Returns an SN-length vector of d(kernel)/dz for each permutation
// This is analogous to partial_z_lnmg but for the kernel
using KernelDerivFn_B = DerivTable(*)(const PairCache&, const VectorXi&,
                                       const MatrixXcd&,
                                       const BasisParams&, const BasisParams&);
using KernelDerivFn_R = DerivTable(*)(const PairCache&, const VectorXi&,
                                       const BasisParams&, const BasisParams&);
using KernelDerivFn_A_perm = DerivTable(*)(const PairCache&, const VectorXi&,
                                            const MatrixXcd&);

// For kernels that don't need perm_matrix (rTr, G, H, Q)
using KernelDerivFn_B_noperm = DerivTable(*)(const PairCache&, const VectorXi&,
                                              const BasisParams&, const BasisParams&);
using KernelDerivFn_A_noperm = DerivTable(*)(const PairCache&, const VectorXi&);

// Generic function for computing the "dMG" part of the Hamiltonian gradient
// This is: sum_p sign[p] * M_G[p] * kernel[p] * dln_MG[p]
// Has the same structure as partial_z_addsn but with kernel multiplication
static Cd addsn_dMG(int alpha_1, bool Real,
                     const BasisParams& pi, const BasisParams& pj,
                     int alpha_2, int alpha_3, int alpha_4,
                     const PermutationSet& perms,
                     const std::vector<Cd>& kernel_vals,
                     const std::vector<Cd>& MG_vals) {
    int SN = perms.SN;

    if (alpha_1 == 1) {
        // dln_MG is 1 for u derivative; just sum sign*M_G*kernel
        Cd result(0, 0);
        for (int p = 0; p < SN; p++)
            result += static_cast<double>(perms.signs[p]) * MG_vals[p] * kernel_vals[p];
        return result;
    }

    // alpha_1 == 2,3,4: need dln_MG tables
    std::vector<Cd> lnmg = partial_z_lnmg(alpha_1, Real, pi, pj,
                                            alpha_2, alpha_3, alpha_4, perms);
    Cd result(0, 0);
    for (int p = 0; p < SN; p++)
        result += static_cast<double>(perms.signs[p]) * MG_vals[p] * kernel_vals[p] * lnmg[p];
    return result;
}

// Generic function for computing the "dKernel" part
// Returns sum_p sign[p] * M_G[p] * d(kernel)/dz[p]
// For alpha_1==1, returns 0.
static std::vector<Cd> partial_z_kernel_generic(
    int alpha_1, bool Real,
    const BasisParams& pi, const BasisParams& pj,
    int alpha_2, int alpha_3, int alpha_4,
    const PermutationSet& perms,
    // Callbacks for each parameter type
    std::function<DerivTable(const PairCache&, const VectorXi&, const MatrixXcd&)> dk_b,
    std::function<DerivTable(const PairCache&, const VectorXi&)> dk_r,
    std::function<DerivTable(const PairCache&, const VectorXi&, const MatrixXcd&)> dk_a)
{
    int SN = perms.SN;
    int N = perms.N;

    if (alpha_1 == 1 || (alpha_2 != pi.name && alpha_2 != pj.name)) {
        return std::vector<Cd>(SN, Cd(0, 0));
    }

    auto [idx1, idx2] = resolve_case_first(alpha_2, pi.name, pj.name, Real);

    if (alpha_1 == 2) {
        std::vector<Cd> result(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            MatrixXcd perm_cd = perms.matrices[p].cast<Cd>();
            DerivTable dt = dk_b(c, perms.sigmas[p], perm_cd);
            Cd val = dt[alpha_3][idx1];
            int a3s = perms.sigmas[p](alpha_3);
            val += dt[a3s][idx2];
            result[p] = val;
        }
        return result;
    }

    if (alpha_1 == 3) {
        std::vector<Cd> result(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            DerivTable dt = dk_r(c, perms.sigmas[p]);
            Cd val = dt[alpha_3][idx1];
            int a3s = perms.sigmas[p](alpha_3);
            val += dt[a3s][idx2];
            result[p] = val;
        }
        return result;
    }

    if (alpha_1 == 4) {
        int index = upper_index_inclusive(N, alpha_3, alpha_4);
        std::vector<Cd> result(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            MatrixXcd perm_cd = perms.matrices[p].cast<Cd>();
            DerivTable dt = dk_a(c, perms.sigmas[p], perm_cd);
            Cd val = dt[index][idx1];
            int a3s = perms.sigmas[p](alpha_3);
            int a4s = perms.sigmas[p](alpha_4);
            if (a4s < a3s) std::swap(a3s, a4s);
            int index_s = upper_index_inclusive(N, a3s, a4s);
            val += dt[index_s][idx2];
            result[p] = val;
        }
        return result;
    }

    return std::vector<Cd>(SN, Cd(0, 0));
}

// ============================================================
// Kinetic energy gradient
// ============================================================

static Cd addsn_kinetic_dKernel(int alpha_1, bool Real,
                                 const BasisParams& pi, const BasisParams& pj,
                                 int alpha_2, int alpha_3, int alpha_4,
                                 const PermutationSet& perms,
                                 const std::vector<Cd>& MG_vals) {
    if (alpha_1 == 1) return Cd(0, 0);

    auto dk = partial_z_kernel_generic(alpha_1, Real, pi, pj,
        alpha_2, alpha_3, alpha_4, perms,
        [&](const PairCache& c, const VectorXi& sig, const MatrixXcd& pm) {
            return dP_Mij_b(c, sig, pm, pi, pj);
        },
        [&](const PairCache& c, const VectorXi& sig) {
            return dP_Mij_r(c, sig, pi, pj);
        },
        [&](const PairCache& c, const VectorXi& sig, const MatrixXcd& pm) {
            return dP_Mij_a(c, sig, pm);
        });

    Cd result(0, 0);
    int SN = perms.SN;
    for (int p = 0; p < SN; p++)
        result += static_cast<double>(perms.signs[p]) * MG_vals[p] * dk[p];
    return result;
}

Cd calculate_Hamiltonian_kinetic_partial(int a1, int a2, int a3, int a4,
                                          bool Real,
                                          const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);
    double coeff = -(hbar * hbar) / (2.0 * mass);

    Cd S = overlap(basis);
    Cd H_val = kinetic_energy_functional(basis);
    Cd dS = partial_z_first(a1, Real, basis, a2, a3, a4);

    // Compute d<H>/dz (both dMG and dKernel parts)
    Cd dH(0, 0);

    auto loop = [&](auto fn) {
        if (a1 == 1 && Real) {
            for (int t = 0; t < basis_n; t++) {
                Cd con_ui = std::conj(basis[t].u);
                dH += con_ui * fn(basis[t], basis[a2]);
            }
        } else if (a1 == 1 && !Real) {
            for (int t = 0; t < basis_n; t++) {
                Cd uj = basis[t].u;
                dH += uj * fn(basis[a2], basis[t]);
            }
        } else {
            for (int i = 0; i < basis_n; i++) {
                Cd con_ui = std::conj(basis[i].u);
                for (int j = 0; j < basis_n; j++) {
                    Cd uj = basis[j].u;
                    dH += con_ui * uj * fn(basis[i], basis[j]);
                }
            }
        }
    };

    loop([&](const BasisParams& pi, const BasisParams& pj) -> Cd {
        int SN = perms.SN;
        std::vector<Cd> MG_vals(SN), kernel_vals(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            MG_vals[p] = c.M_G;
            kernel_vals[p] = compute_P_Mij(c);
        }

        Cd dMG = addsn_dMG(a1, Real, pi, pj, a2, a3, a4, perms, kernel_vals, MG_vals);
        Cd dK = addsn_kinetic_dKernel(a1, Real, pi, pj, a2, a3, a4, perms, MG_vals);
        return dMG + dK;
    });

    dH *= coeff;

    return (1.0 / S) * dH - (H_val / (S * S)) * dS;
}

// ============================================================
// Harmonic gradient
// ============================================================

static Cd addsn_harmonic_dKernel(int alpha_1, bool Real,
                                  const BasisParams& pi, const BasisParams& pj,
                                  int alpha_2, int alpha_3, int alpha_4,
                                  const PermutationSet& perms,
                                  const std::vector<Cd>& MG_vals) {
    if (alpha_1 == 1) return Cd(0, 0);

    auto dk = partial_z_kernel_generic(alpha_1, Real, pi, pj,
        alpha_2, alpha_3, alpha_4, perms,
        [&](const PairCache& c, const VectorXi& sig, const MatrixXcd&) {
            return drTr_Mij_b(c, sig, pi, pj);
        },
        [&](const PairCache& c, const VectorXi& sig) {
            return drTr_Mij_r(c, sig, pi, pj);
        },
        [&](const PairCache& c, const VectorXi& sig, const MatrixXcd&) {
            return drTr_Mij_a(c, sig);
        });

    Cd result(0, 0);
    int SN = perms.SN;
    for (int p = 0; p < SN; p++)
        result += static_cast<double>(perms.signs[p]) * MG_vals[p] * dk[p];
    return result;
}

Cd calculate_Hamiltonian_harmonic_partial(int a1, int a2, int a3, int a4,
                                           bool Real,
                                           const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);
    double coeff = (mass * omega * omega) / 2.0;

    Cd S = overlap(basis);
    Cd H_val = Harmonic_functional(basis);
    Cd dS = partial_z_first(a1, Real, basis, a2, a3, a4);

    Cd dH(0, 0);

    auto loop = [&](auto fn) {
        if (a1 == 1 && Real) {
            for (int t = 0; t < basis_n; t++)
                dH += std::conj(basis[t].u) * fn(basis[t], basis[a2]);
        } else if (a1 == 1 && !Real) {
            for (int t = 0; t < basis_n; t++)
                dH += basis[t].u * fn(basis[a2], basis[t]);
        } else {
            for (int i = 0; i < basis_n; i++)
                for (int j = 0; j < basis_n; j++)
                    dH += std::conj(basis[i].u) * basis[j].u * fn(basis[i], basis[j]);
        }
    };

    loop([&](const BasisParams& pi, const BasisParams& pj) -> Cd {
        int SN = perms.SN;
        std::vector<Cd> MG_vals(SN), kernel_vals(SN);
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            MG_vals[p] = c.M_G;
            kernel_vals[p] = compute_rTr_Mij(c);
        }
        return addsn_dMG(a1, Real, pi, pj, a2, a3, a4, perms, kernel_vals, MG_vals)
             + addsn_harmonic_dKernel(a1, Real, pi, pj, a2, a3, a4, perms, MG_vals);
    });

    dH *= coeff;
    return (1.0 / S) * dH - (H_val / (S * S)) * dS;
}

// ============================================================
// Delta contact gradient
// ============================================================

Cd calculate_Hamiltonian_delta_partial(int a1, int a2, int a3, int a4,
                                        bool Real,
                                        const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    if (N < 2) return Cd(0, 0);

    PermutationSet perms = PermutationSet::generate(N);
    double coeff = g_contact;

    Cd S = overlap(basis);
    Cd H_val = Delta_contact_functional(basis);
    Cd dS = partial_z_first(a1, Real, basis, a2, a3, a4);

    Cd dH(0, 0);

    auto loop = [&](auto fn) {
        if (a1 == 1 && Real) {
            for (int t = 0; t < basis_n; t++)
                dH += std::conj(basis[t].u) * fn(basis[t], basis[a2]);
        } else if (a1 == 1 && !Real) {
            for (int t = 0; t < basis_n; t++)
                dH += basis[t].u * fn(basis[a2], basis[t]);
        } else {
            for (int i = 0; i < basis_n; i++)
                for (int j = 0; j < basis_n; j++)
                    dH += std::conj(basis[i].u) * basis[j].u * fn(basis[i], basis[j]);
        }
    };

    loop([&](const BasisParams& pi, const BasisParams& pj) -> Cd {
        int SN = perms.SN;
        std::vector<Cd> MG_vals(SN);
        // Kernel = sum over pairs a<b of G_Mijab
        std::vector<Cd> kernel_vals(SN, Cd(0, 0));
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            MG_vals[p] = c.M_G;
            kernel_vals[p] = compute_G_Mij(c);
        }

        Cd dMG = addsn_dMG(a1, Real, pi, pj, a2, a3, a4, perms, kernel_vals, MG_vals);

        // dKernel part: sum over pairs
        Cd dK(0, 0);
        if (a1 != 1) {
            if (a2 == pi.name || a2 == pj.name) {
                auto [idx1, idx2] = resolve_case_first(a2, pi.name, pj.name, Real);

                for (int p = 0; p < SN; p++) {
                    PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
                    Cd dk_val(0, 0);

                    for (int pa = 0; pa < N; pa++) {
                        for (int pb = pa + 1; pb < N; pb++) {
                            DerivTable dt;
                            if (a1 == 2) dt = dG_Mijab_b(c, perms.sigmas[p], pa, pb, pi, pj);
                            else if (a1 == 3) dt = dG_Mijab_r(c, perms.sigmas[p], pa, pb, pi, pj);
                            else dt = dG_Mijab_a(c, perms.sigmas[p], pa, pb);

                            int idx, idx_s;
                            if (a1 == 4) {
                                idx = upper_index_inclusive(N, a3, a4);
                                int a3s = perms.sigmas[p](a3), a4s = perms.sigmas[p](a4);
                                if (a4s < a3s) std::swap(a3s, a4s);
                                idx_s = upper_index_inclusive(N, a3s, a4s);
                            } else {
                                idx = a3;
                                idx_s = perms.sigmas[p](a3);
                            }
                            dk_val += dt[idx][idx1] + dt[idx_s][idx2];
                        }
                    }
                    dK += static_cast<double>(perms.signs[p]) * MG_vals[p] * dk_val;
                }
            }
        }

        return dMG + dK;
    });

    dH *= coeff;
    return (1.0 / S) * dH - (H_val / (S * S)) * dS;
}

// ============================================================
// Gaussian interaction gradient
// ============================================================

Cd calculate_Hamiltonian_gaussian_partial(int a1, int a2, int a3, int a4,
                                           bool Real,
                                           const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    if (N < 2) return Cd(0, 0);

    PermutationSet perms = PermutationSet::generate(N);
    double coeff = g_gauss;

    Cd S = overlap(basis);
    Cd H_val = Gaussian_interaction_functional(basis);
    Cd dS = partial_z_first(a1, Real, basis, a2, a3, a4);

    Cd dH(0, 0);

    auto loop = [&](auto fn) {
        if (a1 == 1 && Real) {
            for (int t = 0; t < basis_n; t++)
                dH += std::conj(basis[t].u) * fn(basis[t], basis[a2]);
        } else if (a1 == 1 && !Real) {
            for (int t = 0; t < basis_n; t++)
                dH += basis[t].u * fn(basis[a2], basis[t]);
        } else {
            for (int i = 0; i < basis_n; i++)
                for (int j = 0; j < basis_n; j++)
                    dH += std::conj(basis[i].u) * basis[j].u * fn(basis[i], basis[j]);
        }
    };

    loop([&](const BasisParams& pi, const BasisParams& pj) -> Cd {
        int SN = perms.SN;
        std::vector<Cd> MG_vals(SN);
        std::vector<Cd> kernel_vals(SN, Cd(0, 0));
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            MG_vals[p] = c.M_G;
            kernel_vals[p] = compute_H_Mij(c);
        }

        Cd dMG = addsn_dMG(a1, Real, pi, pj, a2, a3, a4, perms, kernel_vals, MG_vals);

        Cd dK(0, 0);
        if (a1 != 1) {
            if (a2 == pi.name || a2 == pj.name) {
                auto [idx1, idx2] = resolve_case_first(a2, pi.name, pj.name, Real);

                for (int p = 0; p < SN; p++) {
                    PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
                    Cd dk_val(0, 0);
                    for (int pa = 0; pa < N; pa++) {
                        for (int pb = pa + 1; pb < N; pb++) {
                            DerivTable dt;
                            if (a1 == 2) dt = dH_Mijab_b(c, perms.sigmas[p], pa, pb, pi, pj);
                            else if (a1 == 3) dt = dH_Mijab_r(c, perms.sigmas[p], pa, pb, pi, pj);
                            else dt = dH_Mijab_a(c, perms.sigmas[p], pa, pb);

                            int idx, idx_s;
                            if (a1 == 4) {
                                idx = upper_index_inclusive(N, a3, a4);
                                int a3s = perms.sigmas[p](a3), a4s = perms.sigmas[p](a4);
                                if (a4s < a3s) std::swap(a3s, a4s);
                                idx_s = upper_index_inclusive(N, a3s, a4s);
                            } else {
                                idx = a3;
                                idx_s = perms.sigmas[p](a3);
                            }
                            dk_val += dt[idx][idx1] + dt[idx_s][idx2];
                        }
                    }
                    dK += static_cast<double>(perms.signs[p]) * MG_vals[p] * dk_val;
                }
            }
        }

        return dMG + dK;
    });

    dH *= coeff;
    return (1.0 / S) * dH - (H_val / (S * S)) * dS;
}

// ============================================================
// Kicking term gradient
// ============================================================

Cd calculate_Hamiltonian_kicking_partial(int a1, int a2, int a3, int a4,
                                          bool Real,
                                          const std::vector<BasisParams>& basis) {
    int basis_n = static_cast<int>(basis.size());
    int N = basis[0].N();
    PermutationSet perms = PermutationSet::generate(N);
    double coeff = hbar * kappa;

    Cd S = overlap(basis);
    Cd H_val = kicking_term_functional(basis);
    Cd dS = partial_z_first(a1, Real, basis, a2, a3, a4);

    Cd dH(0, 0);

    auto loop = [&](auto fn) {
        if (a1 == 1 && Real) {
            for (int t = 0; t < basis_n; t++)
                dH += std::conj(basis[t].u) * fn(basis[t], basis[a2]);
        } else if (a1 == 1 && !Real) {
            for (int t = 0; t < basis_n; t++)
                dH += basis[t].u * fn(basis[a2], basis[t]);
        } else {
            for (int i = 0; i < basis_n; i++)
                for (int j = 0; j < basis_n; j++)
                    dH += std::conj(basis[i].u) * basis[j].u * fn(basis[i], basis[j]);
        }
    };

    loop([&](const BasisParams& pi, const BasisParams& pj) -> Cd {
        int SN = perms.SN;
        std::vector<Cd> MG_vals(SN);
        std::vector<Cd> kernel_vals(SN, Cd(0, 0));
        for (int p = 0; p < SN; p++) {
            PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
            MG_vals[p] = c.M_G;
            kernel_vals[p] = compute_Q_Mij(c);
        }

        Cd dMG = addsn_dMG(a1, Real, pi, pj, a2, a3, a4, perms, kernel_vals, MG_vals);

        Cd dK(0, 0);
        if (a1 != 1) {
            if (a2 == pi.name || a2 == pj.name) {
                auto [idx1, idx2] = resolve_case_first(a2, pi.name, pj.name, Real);

                for (int p = 0; p < SN; p++) {
                    PairCache c = PairCache::build(pi, pj, perms.matrices[p]);
                    Cd dk_val(0, 0);
                    for (int pa = 0; pa < N; pa++) {
                        DerivTable dt;
                        if (a1 == 2) dt = dQ_Mija_b(c, perms.sigmas[p], pa, pi, pj);
                        else if (a1 == 3) dt = dQ_Mija_r(c, perms.sigmas[p], pa, pi, pj);
                        else dt = dQ_Mija_a(c, perms.sigmas[p], pa);

                        int idx, idx_s;
                        if (a1 == 4) {
                            idx = upper_index_inclusive(N, a3, a4);
                            int a3s = perms.sigmas[p](a3), a4s = perms.sigmas[p](a4);
                            if (a4s < a3s) std::swap(a3s, a4s);
                            idx_s = upper_index_inclusive(N, a3s, a4s);
                        } else {
                            idx = a3;
                            idx_s = perms.sigmas[p](a3);
                        }
                        dk_val += dt[idx][idx1] + dt[idx_s][idx2];
                    }
                    dK += static_cast<double>(perms.signs[p]) * MG_vals[p] * dk_val;
                }
            }
        }

        return dMG + dK;
    });

    dH *= coeff;
    return (1.0 / S) * dH - (H_val / (S * S)) * dS;
}

} // namespace ecg1d
