#pragma once
#include "types.hpp"

namespace ecg1d {

// Rothe VarPro matrices for fixed α (basis_trial = basis_psi):
//   Â = 1 + i(dt/2)·H,   Â† = 1 − i(dt/2)·H
//   Â†·Â  = 1 + (dt/2)²·H²
//   (Â†)² = 1 − i·dt·H − (dt/2)²·H²
//
// S̃_mn = ⟨Â φ_m | Â φ_n⟩ = S_mn + (dt/2)² · (H²)_mn
// ρ̃_m  = Σ_k M_mk · u_k,  where M = S − i·dt·H − (dt/2)²·H²
//        (⟨Â φ_m | Â† ψ⟩ with ψ = Σ_k u_k |φ_k⟩ and basis_trial = basis_psi)

MatrixXcd build_rothe_Stilde(const MatrixXcd& S, const MatrixXcd& H2, double dt);

MatrixXcd build_rothe_M(const MatrixXcd& S, const MatrixXcd& H,
                        const MatrixXcd& H2, double dt);

} // namespace ecg1d
