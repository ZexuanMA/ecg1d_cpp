#include "h_squared.hpp"
#include "physical_constants.hpp"
#include <stdexcept>
#include <string>

namespace ecg1d {

/* =====================================================================
   DERIVATION (N=1, ℏ=m=ω=1; H = −½∂² + ½x²)

   PairCache (N=1) scalar fields:
     k ≡ K_Mj(0,0) = Aⱼ + Bⱼ           (ket-side quadratic coef)
     g ≡ g_Mj(0)  = 2·Rⱼ·Bⱼ            (ket-side linear coef)
     κ ≡ K_inv(0,0)                    (inverse of combined quadratic)
     μ ≡ mu(0)    = ½·K_inv·b           (mean position under combined Gaussian)

   Moments of the combined bra·ket Gaussian (normalized by M_G):
     m1 = μ
     m2 = κ/2 + μ²
     m3 = 3μ·κ/2 + μ³
     m4 = 3κ²/4 + 3κ·μ² + μ⁴
   (Isserlis; variance of combined Gaussian is κ/2.)

   Let L(x) = −2k·x + g   (log-derivative of φⱼ)
       P(x) = L² + L' = 4k²·x² − 4kg·x + (g² − 2k)   (so ∂²φⱼ = P·φⱼ)

   T² = ¼ ∂⁴:
     ∂⁴φⱼ/φⱼ = P² + 2L·P' + P''
             = 16k⁴·x⁴ − 32k³g·x³ + (24k²g² − 48k³)·x²
               + (48k²g − 8kg³)·x + (g⁴ − 12kg² + 12k²)
     ⟨∂⁴⟩/M_G = 16k⁴·m4 − 32k³g·m3 + (24k²g² − 48k³)·m2
                + (48k²g − 8kg³)·m1 + (g⁴ − 12kg² + 12k²)
     compute_T2_Mij = ¼·⟨∂⁴⟩/M_G

   V² = ¼ x⁴:
     compute_V2_Mij = ¼·m4

   TV + VT = −¼ (∂²x² + x²∂²):
     x²·∂²φⱼ/φⱼ = x²·P
     ∂²(x²·φⱼ)/φⱼ = 2 + 4g·x + (g²−10k)·x² − 4kg·x³ + 4k²·x⁴
     Sum/M_G       = 8k²·m4 − 8kg·m3 + (2g² − 12k)·m2 + 4g·m1 + 2
     compute_TV_plus_VT_Mij = −¼·Sum

   Sanity (K=1 HO ground state: k=0.5, g=0, κ=1, μ=0; m4=0.75):
     T² = 0.1875,  V² = 0.1875,  TV+VT = −0.125   → ⟨H²⟩/M_G = 0.25 = E² ✓
   ===================================================================== */

namespace {

struct Moments {
    Cd m1, m2, m3, m4;
};

Moments compute_moments_N1(const PairCache& c) {
    Cd kappa = c.K_inv(0, 0);
    Cd mu    = c.mu(0);
    Cd mu2   = mu * mu;
    Cd kh    = 0.5 * kappa;

    Moments m;
    m.m1 = mu;
    m.m2 = kh + mu2;
    m.m3 = 3.0 * mu * kh + mu * mu2;
    m.m4 = 0.75 * kappa * kappa + 3.0 * kappa * mu2 + mu2 * mu2;
    return m;
}

void check_N1(const PairCache& c) {
    if (c.N != 1) {
        throw std::runtime_error(
            "h_squared: N=1 only in M1 scope (got N=" + std::to_string(c.N) +
            "); multi-particle H² deferred to M5");
    }
}

}  // anonymous namespace

Cd compute_T2_Mij(const PairCache& c) {
    check_N1(c);
    Cd k = c.K_Mj(0, 0);
    Cd g = c.g_Mj(0);
    Moments m = compute_moments_N1(c);

    Cd k2 = k * k,  k3 = k2 * k,  k4 = k2 * k2;
    Cd g2 = g * g,  g3 = g2 * g,  g4 = g2 * g2;

    Cd d4 =  16.0 * k4 * m.m4
          -  32.0 * k3 * g  * m.m3
          + (24.0 * k2 * g2 - 48.0 * k3) * m.m2
          + (48.0 * k2 * g  -  8.0 * k  * g3) * m.m1
          + (g4 - 12.0 * k * g2 + 12.0 * k2);

    return 0.25 * d4;
}

Cd compute_V2_Mij(const PairCache& c) {
    check_N1(c);
    Moments m = compute_moments_N1(c);
    return 0.25 * m.m4;
}

Cd compute_TV_plus_VT_Mij(const PairCache& c) {
    check_N1(c);
    Cd k = c.K_Mj(0, 0);
    Cd g = c.g_Mj(0);
    Moments m = compute_moments_N1(c);

    Cd k2 = k * k;
    Cd g2 = g * g;

    Cd sum =  8.0 * k2 * m.m4
           -  8.0 * k  * g * m.m3
           + (2.0 * g2 - 12.0 * k) * m.m2
           +  4.0 * g  * m.m1
           +  2.0;

    return -0.25 * sum;
}

Cd compute_H2_ij(const BasisParams& bi, const BasisParams& bj,
                 const PermutationSet& perms,
                 const HamiltonianTerms& terms) {
    if (terms.delta || terms.gaussian || terms.kicking) {
        throw std::runtime_error(
            "compute_H2_ij: delta/gaussian/kicking H² not implemented (deferred to M5). "
            "Use HamiltonianTerms::kinetic_harmonic() for M1.");
    }
    if (!(terms.kinetic || terms.harmonic)) {
        return Cd(0.0, 0.0);
    }

    Cd h2_ij(0.0, 0.0);
    for (int p = 0; p < perms.SN; p++) {
        double sign = static_cast<double>(perms.signs[p]);
        PairCache c = PairCache::build(bi, bj, perms.matrices[p]);

        Cd kernel(0.0, 0.0);
        if (terms.kinetic)  kernel += compute_T2_Mij(c);
        if (terms.harmonic) kernel += compute_V2_Mij(c);
        if (terms.kinetic && terms.harmonic) {
            kernel += compute_TV_plus_VT_Mij(c);
        }

        h2_ij += sign * c.M_G * kernel;
    }
    return h2_ij;
}

MatrixXcd build_H2(const std::vector<BasisParams>& basis,
                   const PermutationSet& perms,
                   const HamiltonianTerms& terms) {
    int n = static_cast<int>(basis.size());
    if (n > 0 && basis[0].N() != 1) {
        throw std::runtime_error(
            "build_H2: N=1 only in M1 scope (got N=" +
            std::to_string(basis[0].N()) + ")");
    }
    MatrixXcd H2 = MatrixXcd::Zero(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            Cd h2 = compute_H2_ij(basis[i], basis[j], perms, terms);
            H2(i, j) = h2;
            H2(j, i) = std::conj(h2);
        }
    }
    return H2;
}

}  // namespace ecg1d
