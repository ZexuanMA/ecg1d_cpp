#include "interaction_kernels.hpp"
#include "physical_constants.hpp"
#include <cmath>

namespace ecg1d {

Cd compute_P_Mij(const PairCache& c) {
    // term1: -2 * tr(K_Mj)
    Cd term1 = -2.0 * c.K_Mj.trace();

    // term2: g_Mj^T * g_Mj
    Cd term2 = (c.g_Mj.transpose() * c.g_Mj)(0);

    // term3: -4 * g_Mj^T * K_Mj * mu
    Cd term3 = -4.0 * (c.g_Mj.transpose() * c.K_Mj * c.mu)(0);

    // term4: 4 * mu^T * K_Mj * K_Mj * mu
    Cd term4 = 4.0 * (c.mu.transpose() * c.K_Mj * c.K_Mj * c.mu)(0);

    // term5: 2 * tr(K_Mj * K_Mj * K_inv)
    Cd term5 = 2.0 * (c.K_Mj * c.K_Mj * c.K_inv).trace();

    return term1 + term2 + term3 + term4 + term5;
}

Cd compute_rTr_Mij(const PairCache& c) {
    Cd trace_term = 0.5 * c.K_inv.trace();
    Cd mu_term = (c.mu.transpose() * c.mu)(0);
    return trace_term + mu_term;
}

Cd compute_G_Mijab(const PairCache& c, int a, int b) {
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);
    return 1.0 / std::sqrt(M_PI * h) * std::exp(-p * p / h);
}

Cd compute_G_Mij(const PairCache& c) {
    Cd sum(0.0, 0.0);
    for (int a = 0; a < c.N; a++)
        for (int b = a + 1; b < c.N; b++)
            sum += compute_G_Mijab(c, a, b);
    return sum;
}

Cd compute_H_Mijab(const PairCache& c, int a, int b) {
    Cd h = c.K_inv(a, a) + c.K_inv(b, b) - 2.0 * c.K_inv(a, b);
    Cd p = c.mu(a) - c.mu(b);
    Cd h_eff = sigma_gauss * sigma_gauss + h;
    return sigma_gauss / std::sqrt(h_eff) * std::exp(-p * p / h_eff);
}

Cd compute_H_Mij(const PairCache& c) {
    Cd sum(0.0, 0.0);
    for (int a = 0; a < c.N; a++)
        for (int b = a + 1; b < c.N; b++)
            sum += compute_H_Mijab(c, a, b);
    return sum;
}

Cd compute_Q_Mija(const PairCache& c, int a) {
    Cd le = c.K_inv(a, a);
    Cd lc = 2.0 * c.mu(a);
    return std::exp(-k_L * k_L * le) * std::cos(k_L * lc);
}

Cd compute_Q_Mij(const PairCache& c) {
    Cd sum(0.0, 0.0);
    for (int a = 0; a < c.N; a++)
        sum += compute_Q_Mija(c, a);
    return sum;
}

} // namespace ecg1d
