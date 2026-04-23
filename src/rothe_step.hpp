#pragma once
#include "types.hpp"

namespace ecg1d {

// VarPro linear solve c = S̃⁻¹·(M·u_old).
// Uses Eigen PartialPivLU (works whether or not S̃ is explicitly Hermitian PD).
VectorXcd rothe_varpro_solve_c(const MatrixXcd& S_tilde, const MatrixXcd& M,
                               const VectorXcd& u_old);

// Rothe residual squared:  r² = ‖Â·c_new·|φ⟩ − Â†·u_old·|φ⟩‖²
// Formula (basis_trial = basis_psi):
//   r² = c*·S̃·c − 2·Re(c*·M·u) + u*·S̃·u
// On an eigenstate at VarPro optimum, r² ≡ 0 to machine precision.
double compute_rothe_residual_sq(const MatrixXcd& S_tilde, const MatrixXcd& M,
                                 const VectorXcd& u_old,
                                 const VectorXcd& c_new);

} // namespace ecg1d
