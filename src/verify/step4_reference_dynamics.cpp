#include "step4_reference_dynamics.hpp"
#include "common.hpp"
#include "ecg_wavefunction.hpp"
#include "step2_reference_gs.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace ecg1d {
namespace verify {

namespace {

// Direct DFT psi(x) -> psi_tilde(k) on user-supplied k-grid:
//   psi_tilde(k) = (2 pi)^{-1/2} integral psi(x) exp(-i k x) dx.
// Returns the FULL complex psi_tilde (not |psi_tilde|^2, unlike compute_nk).
Eigen::VectorXcd dft_to_kgrid(const Eigen::VectorXcd& psi,
                              const Eigen::VectorXd& x,
                              const Eigen::VectorXd& k_grid) {
    double dx = x(1) - x(0);
    Eigen::VectorXcd psi_k(k_grid.size());
    double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    for (int ki = 0; ki < k_grid.size(); ki++) {
        Cd acc(0.0, 0.0);
        double kv = k_grid(ki);
        for (int j = 0; j < x.size(); j++) {
            acc += psi(j) * std::exp(Cd(0.0, -kv * x(j)));
        }
        psi_k(ki) = acc * dx * inv_sqrt_2pi;
    }
    return psi_k;
}

// Trapezoid integral on a uniform grid (dx constant).
double trapezoid(const Eigen::VectorXd& integrand, double dx) {
    double s = 0.0;
    const int n = static_cast<int>(integrand.size());
    if (n == 0) return 0.0;
    s = 0.5 * (integrand(0) + integrand(n - 1));
    for (int j = 1; j < n - 1; j++) s += integrand(j);
    return s * dx;
}

// Core spectral-propagation engine shared by grid and DVR. T_kinetic is the
// caller-supplied real-symmetric kinetic-energy matrix (FD or sinc-DVR).
Step4RefResult propagate_spectral(const Eigen::VectorXd& x,
                                  const Eigen::MatrixXd& T_kinetic,
                                  const std::vector<BasisParams>& ecg_psi0_basis,
                                  const Eigen::VectorXd& k_grid,
                                  double T_total,
                                  const std::vector<double>& t_snap,
                                  const std::vector<double>& t_trace) {
    (void)T_total;  // diagnostics only; trace times themselves bound t

    const int Nx = static_cast<int>(x.size());
    const double dx = x(1) - x(0);

    // Build H_evolve = T + V_ho (no cos)
    Eigen::VectorXd V = detail::build_V_ho_only(x);
    Eigen::MatrixXd H = T_kinetic;
    for (int j = 0; j < Nx; j++) H(j, j) += V(j);

    // Diagonalize once
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("step4_reference_dynamics: eigensolver failed");
    }
    const Eigen::VectorXd lambda = es.eigenvalues();    // (Nx,)
    const Eigen::MatrixXd U      = es.eigenvectors();   // (Nx, Nx)

    // Initial state: ECG ground state evaluated on x; renormalize.
    Eigen::VectorXcd psi0 = ecg_wavefunction_1p(ecg_psi0_basis, x);
    detail::normalize_grid_state(psi0, dx);

    // Spectral coefficients c = U^T psi0  (U is real)
    Eigen::VectorXcd c = U.transpose().cast<Cd>() * psi0;

    // E_init via <psi0|H|psi0> dx (should equal sum_n |c_n|^2 lambda_n; we use
    // both as a self-consistency check).
    Cd E_init_psi  = (psi0.adjoint() * (H.cast<Cd>() * psi0))(0) * dx;
    double E_init  = E_init_psi.real();
    {
        // Cross-check: spectral energy = sum_n |c_n|^2 lambda_n * dx
        // (factor of dx because c are L2-discrete coefficients of an L2-continuous psi)
        double E_spec = 0.0;
        for (int n = 0; n < Nx; n++) E_spec += std::norm(c(n)) * lambda(n);
        E_spec *= dx;
        if (std::abs(E_spec - E_init) / std::max(std::abs(E_init), 1e-14) > 1e-10) {
            std::cerr << "[step4_ref] warning: E_init mismatch psi=" << E_init
                      << " spectral=" << E_spec << "\n";
        }
    }

    // Helper: propagate to time t and return psi(t) on x.
    auto propagate_to = [&](double t) -> Eigen::VectorXcd {
        Eigen::VectorXcd phases(Nx);
        for (int n = 0; n < Nx; n++) {
            phases(n) = std::exp(Cd(0.0, -lambda(n) * t)) * c(n);
        }
        return U.cast<Cd>() * phases;
    };

    // Trace observables
    const int n_tr = static_cast<int>(t_trace.size());
    Step4RefResult r;
    r.t_trace = Eigen::VectorXd(n_tr);
    r.E       = Eigen::VectorXd(n_tr);
    r.x_mean  = Eigen::VectorXd(n_tr);
    r.p_mean  = Eigen::VectorXd(n_tr);
    r.x2      = Eigen::VectorXd(n_tr);
    r.p2      = Eigen::VectorXd(n_tr);
    r.norm    = Eigen::VectorXd(n_tr);
    r.overlap0 = Eigen::VectorXcd(n_tr);
    r.x = x;
    r.k = k_grid;
    r.E_init = E_init;
    double max_dE = 0.0;
    for (int i = 0; i < n_tr; i++) {
        double t = t_trace[i];
        r.t_trace(i) = t;

        Eigen::VectorXcd psi = propagate_to(t);

        // norm = sum |psi|^2 dx (local renamed from `p2` to `prob` to avoid
        // shadowing the new r.p2 trace field).
        Eigen::VectorXd prob(Nx);
        for (int j = 0; j < Nx; j++) prob(j) = std::norm(psi(j));
        double n_t = trapezoid(prob, dx);
        r.norm(i) = n_t;
        const double inv_n_t = (n_t > 0.0) ? 1.0 / n_t : 0.0;

        // <x>, <x^2> in position space (normalized by n_t for safety even though
        // unitary spectral propagation keeps n_t = 1 to round-off).
        Eigen::VectorXd x_p(Nx), x2_p(Nx);
        for (int j = 0; j < Nx; j++) {
            x_p(j)  = x(j) * prob(j);
            x2_p(j) = x(j) * x(j) * prob(j);
        }
        r.x_mean(i) = trapezoid(x_p,  dx) * inv_n_t;
        r.x2(i)     = trapezoid(x2_p, dx) * inv_n_t;

        // <p>, <p^2> via k-space integral on k_grid (one DFT, two moments).
        Eigen::VectorXcd psi_k = dft_to_kgrid(psi, x, k_grid);
        double dk = k_grid(1) - k_grid(0);
        Eigen::VectorXd kpk(k_grid.size()), k2pk(k_grid.size());
        for (int ki = 0; ki < k_grid.size(); ki++) {
            double pk = std::norm(psi_k(ki));
            kpk(ki)  = k_grid(ki)               * pk;
            k2pk(ki) = k_grid(ki) * k_grid(ki)  * pk;
        }
        r.p_mean(i) = trapezoid(kpk,  dk) * inv_n_t;
        r.p2(i)     = trapezoid(k2pk, dk) * inv_n_t;

        // <H>: should be conserved exactly to roundoff
        Cd Et = (psi.adjoint() * (H.cast<Cd>() * psi))(0) * dx;
        r.E(i) = Et.real();
        max_dE = std::max(max_dE, std::abs(r.E(i) - E_init));

        // <psi(0)|psi(t)>
        Cd ov(0.0, 0.0);
        for (int j = 0; j < Nx; j++) ov += std::conj(psi0(j)) * psi(j);
        ov *= dx;
        r.overlap0(i) = ov;
    }
    r.max_rel_E_drift = max_dE / std::max(std::abs(E_init), 1e-14);

    // Snapshots
    const int n_sn = static_cast<int>(t_snap.size());
    r.t_snap = Eigen::VectorXd(n_sn);
    r.density_snaps = Eigen::MatrixXd(n_sn, Nx);
    r.nk_snaps      = Eigen::MatrixXd(n_sn, k_grid.size());
    for (int i = 0; i < n_sn; i++) {
        double t = t_snap[i];
        r.t_snap(i) = t;

        Eigen::VectorXcd psi = propagate_to(t);
        for (int j = 0; j < Nx; j++) r.density_snaps(i, j) = std::norm(psi(j));
        Eigen::VectorXd nk = detail::compute_nk(psi, x, k_grid);
        for (int kj = 0; kj < k_grid.size(); kj++) r.nk_snaps(i, kj) = nk(kj);
    }

    return r;
}

} // anonymous namespace

Step4RefResult step4_realtime_grid(const std::vector<BasisParams>& ecg_psi0_basis,
                                   int n_grid, double x_max,
                                   const Eigen::VectorXd& k_grid,
                                   double T_total,
                                   const std::vector<double>& t_snap,
                                   const std::vector<double>& t_trace) {
    std::cout << "[step4-ref/grid] n_grid=" << n_grid << " x_max=" << x_max
              << " trace_pts=" << t_trace.size() << " snap_pts=" << t_snap.size() << "\n";
    double dx = 2.0 * x_max / (n_grid - 1);
    Eigen::VectorXd x(n_grid);
    for (int j = 0; j < n_grid; j++) x(j) = -x_max + j * dx;
    Eigen::MatrixXd T = detail::build_T_grid_fd(x);
    Step4RefResult r = propagate_spectral(x, T, ecg_psi0_basis, k_grid,
                                          T_total, t_snap, t_trace);
    std::cout << "[step4-ref/grid] E_init=" << std::setprecision(12) << r.E_init
              << " max_rel_E_drift=" << std::scientific << r.max_rel_E_drift
              << std::defaultfloat << "\n";
    return r;
}

Step4RefResult step4_realtime_dvr(const std::vector<BasisParams>& ecg_psi0_basis,
                                  int n_dvr, double x_max,
                                  const Eigen::VectorXd& k_grid,
                                  double T_total,
                                  const std::vector<double>& t_snap,
                                  const std::vector<double>& t_trace) {
    std::cout << "[step4-ref/dvr] n_dvr=" << n_dvr << " x_max=" << x_max
              << " trace_pts=" << t_trace.size() << " snap_pts=" << t_snap.size() << "\n";
    double dx = 2.0 * x_max / (n_dvr - 1);
    Eigen::VectorXd x(n_dvr);
    for (int j = 0; j < n_dvr; j++) x(j) = -x_max + j * dx;
    Eigen::MatrixXd T = detail::build_T_dvr_sinc(x);
    Step4RefResult r = propagate_spectral(x, T, ecg_psi0_basis, k_grid,
                                          T_total, t_snap, t_trace);
    std::cout << "[step4-ref/dvr] E_init=" << std::setprecision(12) << r.E_init
              << " max_rel_E_drift=" << std::scientific << r.max_rel_E_drift
              << std::defaultfloat << "\n";
    return r;
}

void write_step4_ref_csvs(const std::string& tag, int N,
                          const Step4RefResult& r) {
    std::ostringstream suf;
    suf << "_N" << N << ".csv";

    // Trace. x2/p2 appended at end so existing positional readers (which
    // expect cols 0-5 = t,E,norm,x_mean,p_mean,fidelity) still work.
    {
        std::ofstream f(out_path("step4_" + tag + "_trace" + suf.str()));
        f << "t,E,norm,x_mean,p_mean,fidelity,x2,p2\n";
        f << std::setprecision(15);
        const double n0 = (r.norm.size() > 0) ? r.norm(0) : 1.0;
        for (int i = 0; i < r.t_trace.size(); i++) {
            double fid = 0.0;
            if (i < r.overlap0.size()) {
                double denom = n0 * r.norm(i);
                fid = (denom > 0) ? std::norm(r.overlap0(i)) / denom : 0.0;
            }
            f << r.t_trace(i) << "," << r.E(i) << "," << r.norm(i) << ","
              << r.x_mean(i) << "," << r.p_mean(i) << "," << fid << ","
              << r.x2(i)     << "," << r.p2(i)     << "\n";
        }
    }
    // Density (long format)
    {
        std::ofstream f(out_path("step4_" + tag + "_density" + suf.str()));
        f << "t,x,n_x\n";
        f << std::setprecision(15);
        for (int i = 0; i < r.t_snap.size(); i++) {
            for (int j = 0; j < r.x.size(); j++) {
                f << r.t_snap(i) << "," << r.x(j) << ","
                  << r.density_snaps(i, j) << "\n";
            }
        }
    }
    // n(k) (long format)
    {
        std::ofstream f(out_path("step4_" + tag + "_nk" + suf.str()));
        f << "t,k,n_k\n";
        f << std::setprecision(15);
        for (int i = 0; i < r.t_snap.size(); i++) {
            for (int kj = 0; kj < r.k.size(); kj++) {
                f << r.t_snap(i) << "," << r.k(kj) << ","
                  << r.nk_snaps(i, kj) << "\n";
            }
        }
    }
}

} // namespace verify
} // namespace ecg1d
