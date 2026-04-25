#pragma once
#include "basis_params.hpp"
#include "types.hpp"
#include <Eigen/Dense>
#include <vector>

namespace ecg1d {
namespace verify {

// Result of a Step-4 reference real-time evolution under H_evolve = T + V_ho.
//
// Reference solvers (FD grid + sinc-DVR) propagate the same physical initial
// state used by the ECG side -- the Step-1 ECG ground state evaluated and
// renormalized on the discretization grid. Time evolution is done by
// eigendecomposing H_evolve once and applying psi(t) = U exp(-i Λ t) U^T psi0,
// which is *exact in time*; the only discretization error is spatial.
//
// Trace fields are sampled at every t in `t_trace`.
// Snapshot fields hold density / n(k) on the user grids at every t in `t_snap`.
struct Step4RefResult {
    // Trace (dense in time, for E(t), <x>(t), <p>(t), norm(t) plots).
    Eigen::VectorXd t_trace;
    Eigen::VectorXd E;        // <H_evolve> at trace time -- conserved by construction
    Eigen::VectorXd x_mean;
    Eigen::VectorXd p_mean;
    Eigen::VectorXd norm;     // integrated probability (should stay 1)
    Eigen::VectorXcd overlap0; // <psi(0) | psi(t)>  (for fidelity)

    // Snapshot grids (one row per snapshot time).
    Eigen::VectorXd t_snap;
    Eigen::VectorXd x;          // discretization x-grid
    Eigen::VectorXd k;          // user-supplied k-grid (n(k) evaluation)
    Eigen::MatrixXd density_snaps;  // [n_snap × Nx]   |psi(x,t)|^2
    Eigen::MatrixXd nk_snaps;       // [n_snap × Nk]   |psi_tilde(k,t)|^2

    // Initial energy (sanity reference for drift)
    double E_init = 0.0;
    double max_rel_E_drift = 0.0;
};

// Real-time spectral propagation on the FD 3-point grid.
// Initial state = ECG ground state from Step 1, evaluated on the grid via
// `ecg_wavefunction_1p` and trapezoid-renormalized.
Step4RefResult step4_realtime_grid(const std::vector<BasisParams>& ecg_psi0_basis,
                                   int n_grid, double x_max,
                                   const Eigen::VectorXd& k_grid,
                                   double T_total,
                                   const std::vector<double>& t_snap,
                                   const std::vector<double>& t_trace);

// Real-time spectral propagation on the sinc-DVR grid. Same signature.
Step4RefResult step4_realtime_dvr(const std::vector<BasisParams>& ecg_psi0_basis,
                                  int n_dvr, double x_max,
                                  const Eigen::VectorXd& k_grid,
                                  double T_total,
                                  const std::vector<double>& t_snap,
                                  const std::vector<double>& t_trace);

// Write a Step4RefResult to CSVs:
//   step4_{tag}_trace_N{N}.csv     : t,E,norm,x_mean,p_mean,fidelity
//   step4_{tag}_density_N{N}.csv   : t,x,n_x  (long format)
//   step4_{tag}_nk_N{N}.csv        : t,k,n_k  (long format)
void write_step4_ref_csvs(const std::string& tag, int N,
                          const Step4RefResult& r);

} // namespace verify
} // namespace ecg1d
