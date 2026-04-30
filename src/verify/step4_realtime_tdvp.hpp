#pragma once
#include "basis_params.hpp"
#include "realtime_tdvp.hpp"
#include "types.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace ecg1d {
namespace verify {

// Result of the Step-4 ECG-side real-time TDVP run under H_evolve = T + V_ho
// (cos quench off). Trace fields are sampled densely; snapshot fields hold
// density / n(k) on the user grids at every t in `t_snap`. Snapshot bases are
// also stashed so cross-check code can compute exact ECG-vs-grid fidelity.
struct Step4EcgResult {
    // Trace (dense in time). x_mean, p_mean, x2, p2 are normalized: <A>/<psi|psi>.
    // The *_raw versions are <psi|A|psi> without dividing by <psi|psi>.
    Eigen::VectorXd t_trace;
    Eigen::VectorXd E;
    Eigen::VectorXd norm;
    Eigen::VectorXd x_mean;
    Eigen::VectorXd p_mean;
    Eigen::VectorXd x2;
    Eigen::VectorXd p2;
    Eigen::VectorXd x_mean_raw;
    Eigen::VectorXd p_mean_raw;
    Eigen::VectorXd x2_raw;
    Eigen::VectorXd p2_raw;

    // Snapshots
    Eigen::VectorXd t_snap;
    Eigen::VectorXd x;
    Eigen::VectorXd k;
    Eigen::MatrixXd density_snaps; // [n_snap × Nx]
    Eigen::MatrixXd nk_snaps;      // [n_snap × Nk]
    Eigen::VectorXd E_snap;        // <H> at each snapshot
    Eigen::VectorXd fidelity_snap; // |<psi(0)|psi(t)>|^2 / (norm0*norm_t), global

    // Snapshot bases (kept for downstream cross-check — small, K basis params/snap)
    std::vector<std::vector<BasisParams>> basis_snaps;

    double E_init = 0.0;
    double max_rel_E_drift = 0.0;
};

// Which form of the four moments (x_mean, p_mean, x2, p2) to write to the
// trace CSV. Both forms are always computed; this only affects CSV output.
enum class MomentForm { Normalized, Raw, Both };

struct Step4RealtimeOptions {
    double lambda_C = 1e-8;
    double rcond = 1e-4;
    bool wiener_smooth = true;
    bool enforce_norm = false;
    bool u_split_trotter = false;
    MomentForm moment_form = MomentForm::Both;
    RtIntegrator integrator = RtIntegrator::RK4;
};

// Step 4 — ECG real-time TDVP under H_evolve = kinetic + harmonic trap (cos
// turned off). Starts from `basis_init` (Step-1 ECG ground state of H_full).
//
// `t_snap` are snapshot times where density and n(k) are written; t=0 and
// T_total are added if missing. `dt_trace` controls trace sample density.
//
// Snapshots are obtained by chaining `realtime_tdvp_evolution` calls between
// consecutive snapshot times — no re-runs from t=0.
//
// Writes:
//   step4_ecg_trace_N{N}_K{K}.csv     : t,E,norm,x_mean,p_mean
//   step4_ecg_density_N{N}_K{K}.csv   : t,x,n_x  (long format)
//   step4_ecg_nk_N{N}_K{K}.csv        : t,k,n_k  (long format)
//   step4_ecg_snap_N{N}_K{K}.csv      : t,E,fidelity (per snapshot)
Step4EcgResult step4_realtime_tdvp(int N, int K,
                                   const std::vector<BasisParams>& basis_init,
                                   double T_total, double dt,
                                   const Eigen::VectorXd& x_grid,
                                   const Eigen::VectorXd& k_grid,
                                   const std::vector<double>& t_snap,
                                   double dt_trace,
                                   bool rt_verbose = false,
                                   const Step4RealtimeOptions& rt_options = Step4RealtimeOptions{});

} // namespace verify
} // namespace ecg1d
