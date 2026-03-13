#pragma once
#include "pair_cache.hpp"

namespace ecg1d {

// Kinetic energy kernel: P_Mij
// = -2*tr(K_Mj) + g^T g - 4*g^T K_Mj mu + 4*mu^T K_Mj^2 mu + 2*tr(K_Mj^2 K_inv)
Cd compute_P_Mij(const PairCache& c);

// Harmonic potential kernel: rTr_Mij = 0.5*tr(K_inv) + mu^T mu
Cd compute_rTr_Mij(const PairCache& c);

// Delta contact interaction kernel for particle pair (a,b)
// G_Mijab = 1/sqrt(pi*h) * exp(-p^2/h)
// where h = K_inv[a,a] + K_inv[b,b] - 2*K_inv[a,b], p = mu[a] - mu[b]
Cd compute_G_Mijab(const PairCache& c, int a, int b);

// Sum of G_Mijab over all pairs a < b
Cd compute_G_Mij(const PairCache& c);

// Gaussian interaction kernel for particle pair (a,b)
// H_Mijab = sigma / sqrt(sigma^2 + h) * exp(-p^2 / (sigma^2 + h))
Cd compute_H_Mijab(const PairCache& c, int a, int b);

// Sum of H_Mijab over all pairs a < b
Cd compute_H_Mij(const PairCache& c);

// Kicking term kernel for particle a
// Q_Mija = exp(-k_L^2 * K_inv[a,a]) * cos(k_L * 2*mu[a])
Cd compute_Q_Mija(const PairCache& c, int a);

// Sum of Q_Mija over all particles
Cd compute_Q_Mij(const PairCache& c);

} // namespace ecg1d
