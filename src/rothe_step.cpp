#include "rothe_step.hpp"
#include <Eigen/Dense>

namespace ecg1d {

VectorXcd rothe_varpro_solve_c(const MatrixXcd& S_tilde, const MatrixXcd& M,
                               const VectorXcd& u_old) {
    VectorXcd rho_tilde = M * u_old;
    Eigen::PartialPivLU<MatrixXcd> lu(S_tilde);
    return lu.solve(rho_tilde);
}

double compute_rothe_residual_sq(const MatrixXcd& S_tilde, const MatrixXcd& M,
                                 const VectorXcd& u_old,
                                 const VectorXcd& c_new) {
    VectorXcd rho_tilde = M * u_old;
    Cd term1 = c_new.dot(S_tilde * c_new);   // c*·S̃·c
    Cd term2 = c_new.dot(rho_tilde);          // c*·ρ̃
    Cd term3 = u_old.dot(S_tilde * u_old);    // u*·S̃·u  (=‖Â†ψ‖²)
    Cd r2 = term1 - 2.0 * term2.real() + term3;
    return r2.real();
}

} // namespace ecg1d
