#include "rothe_operator.hpp"

namespace ecg1d {

MatrixXcd build_rothe_Stilde(const MatrixXcd& S, const MatrixXcd& H2, double dt) {
    const double half_dt = 0.5 * dt;
    return S + (half_dt * half_dt) * H2;
}

MatrixXcd build_rothe_M(const MatrixXcd& S, const MatrixXcd& H,
                        const MatrixXcd& H2, double dt) {
    const double half_dt = 0.5 * dt;
    const Cd i_dt(0.0, dt);
    return S - i_dt * H - (half_dt * half_dt) * H2;
}

} // namespace ecg1d
