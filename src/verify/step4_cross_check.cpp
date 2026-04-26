#include "step4_cross_check.hpp"
#include "common.hpp"
#include "ecg_wavefunction.hpp"
#include <Eigen/Dense>
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

// Linear interpolation of (xs_in -> ys_in) onto xs_out. Assumes xs_in is sorted.
// Out-of-range values clamp to the endpoints.
Eigen::VectorXd interp1_linear(const Eigen::VectorXd& xs_in,
                               const Eigen::VectorXd& ys_in,
                               const Eigen::VectorXd& xs_out) {
    const int n_in = static_cast<int>(xs_in.size());
    const int n_out = static_cast<int>(xs_out.size());
    Eigen::VectorXd ys_out(n_out);
    int i = 0;
    for (int j = 0; j < n_out; j++) {
        double x = xs_out(j);
        if (x <= xs_in(0))            { ys_out(j) = ys_in(0); continue; }
        if (x >= xs_in(n_in - 1))     { ys_out(j) = ys_in(n_in - 1); continue; }
        while (i + 1 < n_in - 1 && xs_in(i + 1) < x) i++;
        // After the loop: xs_in(i) <= x and (i+1==n_in-1 OR xs_in(i+1) >= x)
        if (xs_in(i) > x) i = std::max(0, i - 1);
        double xa = xs_in(i), xb = xs_in(i + 1);
        double t = (x - xa) / (xb - xa);
        ys_out(j) = (1.0 - t) * ys_in(i) + t * ys_in(i + 1);
    }
    return ys_out;
}

double trap_integral(const Eigen::VectorXd& f, double dx) {
    const int n = static_cast<int>(f.size());
    if (n == 0) return 0.0;
    double s = 0.5 * (f(0) + f(n - 1));
    for (int j = 1; j < n - 1; j++) s += f(j);
    return s * dx;
}

// Resample a snapshot density (nx points x_ref(0..nx-1) -> n_x_ref values) onto
// the user x_grid via linear interpolation.
Eigen::VectorXd resample_density(const Eigen::VectorXd& x_ref,
                                 const Eigen::VectorXd& n_ref,
                                 const Eigen::VectorXd& x_user) {
    return interp1_linear(x_ref, n_ref, x_user);
}

double L2_diff(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double dx) {
    if (a.size() != b.size())
        throw std::runtime_error("L2_diff: size mismatch");
    Eigen::VectorXd d2(a.size());
    for (int j = 0; j < a.size(); j++) {
        double v = a(j) - b(j);
        d2(j) = v * v;
    }
    return std::sqrt(trap_integral(d2, dx));
}

double bhattacharyya(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double dx) {
    Eigen::VectorXd s(a.size());
    for (int j = 0; j < a.size(); j++) {
        double va = std::max(0.0, a(j));
        double vb = std::max(0.0, b(j));
        s(j) = std::sqrt(va * vb);
    }
    return trap_integral(s, dx);
}

// Find the trace sample whose t is closest to t_target.
int nearest_index(const Eigen::VectorXd& ts, double t_target) {
    int best = 0;
    double best_d = std::abs(ts(0) - t_target);
    for (int i = 1; i < ts.size(); i++) {
        double d = std::abs(ts(i) - t_target);
        if (d < best_d) { best_d = d; best = i; }
    }
    return best;
}

// Trace value at t_target via nearest sample. Sufficient when traces are dense.
double trace_at(const Eigen::VectorXd& ts, const Eigen::VectorXd& vals, double t_target) {
    return vals(nearest_index(ts, t_target));
}

// Build the ECG single-particle wavefunction at snapshot s (N=1) on x_grid,
// then renormalize to unit L2.
Eigen::VectorXcd ecg_psi_on_grid_norm(const std::vector<BasisParams>& basis,
                                      const Eigen::VectorXd& x_grid) {
    Eigen::VectorXcd psi = ecg_wavefunction_1p(basis, x_grid);
    double dx = x_grid(1) - x_grid(0);
    double s = 0.0;
    for (int j = 0; j < x_grid.size(); j++) s += std::norm(psi(j));
    s *= dx;
    if (s > 0) psi /= std::sqrt(s);
    return psi;
}

// Resample a complex grid wavefunction (real, imag separately) onto x_grid,
// then renormalize to unit L2.
Eigen::VectorXcd resample_psi(const Eigen::VectorXd& x_in,
                              const Eigen::VectorXd& re_in,
                              const Eigen::VectorXd& im_in,
                              const Eigen::VectorXd& x_user) {
    Eigen::VectorXd re_out = interp1_linear(x_in, re_in, x_user);
    Eigen::VectorXd im_out = interp1_linear(x_in, im_in, x_user);
    Eigen::VectorXcd psi(x_user.size());
    for (int j = 0; j < x_user.size(); j++) psi(j) = Cd(re_out(j), im_out(j));
    double dx = x_user(1) - x_user(0);
    double s = 0.0;
    for (int j = 0; j < x_user.size(); j++) s += std::norm(psi(j));
    s *= dx;
    if (s > 0) psi /= std::sqrt(s);
    return psi;
}

double fidelity(const Eigen::VectorXcd& psi_a, const Eigen::VectorXcd& psi_b, double dx) {
    Cd inner(0.0, 0.0);
    for (int j = 0; j < psi_a.size(); j++) inner += std::conj(psi_a(j)) * psi_b(j);
    inner *= dx;
    return std::norm(inner);  // both are L2-normalized, so this IS the fidelity
}

} // anonymous namespace

void step4_cross_check(int N,
                       const Step4EcgResult& ecg,
                       const Step4RefResult& grid_res,
                       const Step4RefResult& dvr_res,
                       const Eigen::VectorXd& x_grid,
                       const Eigen::VectorXd& k_grid) {
    std::cout << "\n=== Step 4 — cross-check (N=" << N << ") ===\n";

    if (ecg.t_snap.size() != grid_res.t_snap.size() ||
        ecg.t_snap.size() != dvr_res.t_snap.size()) {
        throw std::runtime_error(
            "step4_cross_check: snapshot count mismatch between ECG/grid/DVR");
    }
    const int n_snap = static_cast<int>(ecg.t_snap.size());
    const double dx = x_grid(1) - x_grid(0);
    const double dk = k_grid(1) - k_grid(0);

    std::ostringstream suf;
    suf << "_N" << N << ".csv";
    std::ofstream f(out_path("step4_crosscheck" + suf.str()));
    f << "pair,t,quantity,value\n";
    f << std::setprecision(15);

    auto write = [&](const std::string& pair, double t,
                     const std::string& q, double v) {
        f << pair << "," << t << "," << q << "," << v << "\n";
    };

    // For propagating the reference's complex psi to a snapshot time we'd need
    // it stored — we don't have it (only |psi|^2). Fidelity ECG-vs-ref is
    // computed in N=1 by using the ECG basis-snapshot to build |psi_ECG(t)>
    // on x_grid, and the reference |psi_ref(t)> from sqrt(n_ref(x,t)) (real,
    // assuming the propagator preserves real-valuedness — only valid for
    // real H_evolve and real psi0; the cos-quench breaks this mildly because
    // psi0 itself is real for the harmonic ground state but the cos-on
    // ground state may have small imaginary parts from chirp).
    // Bhattacharyya is the fully-faithful density-level fidelity; we keep it
    // as the headline density check, and report psi-fidelity for N=1 with
    // a note that it uses sqrt(n_ref).

    for (int s = 0; s < n_snap; s++) {
        double t = ecg.t_snap(s);

        // Energies at this snapshot
        double E_ecg  = ecg.E_snap(s);
        double E_grid = trace_at(grid_res.t_trace, grid_res.E, t);
        double E_dvr  = trace_at(dvr_res.t_trace,  dvr_res.E,  t);

        // <x>, <p> at this snapshot (interpolate from trace)
        double x_ecg  = trace_at(ecg.t_trace, ecg.x_mean, t);
        double p_ecg  = trace_at(ecg.t_trace, ecg.p_mean, t);
        double x_grid_v = trace_at(grid_res.t_trace, grid_res.x_mean, t);
        double p_grid_v = trace_at(grid_res.t_trace, grid_res.p_mean, t);
        double x_dvr_v  = trace_at(dvr_res.t_trace,  dvr_res.x_mean,  t);
        double p_dvr_v  = trace_at(dvr_res.t_trace,  dvr_res.p_mean,  t);

        // Densities resampled to user x_grid (ECG already is)
        Eigen::VectorXd nx_ecg = ecg.density_snaps.row(s);
        Eigen::VectorXd nx_grid = resample_density(grid_res.x,
                                                   Eigen::VectorXd(grid_res.density_snaps.row(s)),
                                                   x_grid);
        Eigen::VectorXd nx_dvr  = resample_density(dvr_res.x,
                                                   Eigen::VectorXd(dvr_res.density_snaps.row(s)),
                                                   x_grid);

        // n(k) -- already on user k_grid for all three
        Eigen::VectorXd nk_ecg = ecg.nk_snaps.row(s);
        Eigen::VectorXd nk_grid = grid_res.nk_snaps.row(s);
        Eigen::VectorXd nk_dvr  = dvr_res.nk_snaps.row(s);

        // For fidelity (N=1): build complex psi on x_grid for each.
        Eigen::VectorXcd psi_ecg, psi_grid, psi_dvr;
        if (N == 1) {
            psi_ecg = ecg_psi_on_grid_norm(ecg.basis_snaps[s], x_grid);
            // Reference: take sqrt(n_ref(x,t)) (real, approximate — assumes
            // psi_ref is real-valued at snapshot, which holds for the HO+cos
            // ground state propagated under HO-only since both Hamiltonians
            // are real-symmetric and the initial wavefunction is real).
            Eigen::VectorXd phi_grid_re(x_grid.size());
            Eigen::VectorXd phi_grid_im = Eigen::VectorXd::Zero(x_grid.size());
            for (int j = 0; j < x_grid.size(); j++) phi_grid_re(j) = std::sqrt(std::max(0.0, nx_grid(j)));
            psi_grid = resample_psi(x_grid, phi_grid_re, phi_grid_im, x_grid);
            Eigen::VectorXd phi_dvr_re(x_grid.size());
            Eigen::VectorXd phi_dvr_im = Eigen::VectorXd::Zero(x_grid.size());
            for (int j = 0; j < x_grid.size(); j++) phi_dvr_re(j) = std::sqrt(std::max(0.0, nx_dvr(j)));
            psi_dvr  = resample_psi(x_grid, phi_dvr_re, phi_dvr_im, x_grid);
        }

        auto pair_diag = [&](const std::string& pair,
                             double Ea, double Eb,
                             double xa, double xb,
                             double pa, double pb,
                             const Eigen::VectorXd& nxa, const Eigen::VectorXd& nxb,
                             const Eigen::VectorXd& nka, const Eigen::VectorXd& nkb,
                             const Eigen::VectorXcd* psi_a,
                             const Eigen::VectorXcd* psi_b) {
            write(pair, t, "delta_E",       Ea - Eb);
            write(pair, t, "abs_delta_x",   std::abs(xa - xb));
            write(pair, t, "abs_delta_p",   std::abs(pa - pb));
            write(pair, t, "L2_nx",         L2_diff(nxa, nxb, dx));
            write(pair, t, "L2_nk",         L2_diff(nka, nkb, dk));
            write(pair, t, "bhattacharyya", bhattacharyya(nxa, nxb, dx));
            if (psi_a && psi_b)
                write(pair, t, "fidelity",   fidelity(*psi_a, *psi_b, dx));
        };

        pair_diag("ECG_vs_grid",
                  E_ecg, E_grid, x_ecg, x_grid_v, p_ecg, p_grid_v,
                  nx_ecg, nx_grid, nk_ecg, nk_grid,
                  (N == 1) ? &psi_ecg : nullptr,
                  (N == 1) ? &psi_grid : nullptr);
        pair_diag("ECG_vs_DVR",
                  E_ecg, E_dvr, x_ecg, x_dvr_v, p_ecg, p_dvr_v,
                  nx_ecg, nx_dvr, nk_ecg, nk_dvr,
                  (N == 1) ? &psi_ecg : nullptr,
                  (N == 1) ? &psi_dvr : nullptr);
        pair_diag("grid_vs_DVR",
                  E_grid, E_dvr, x_grid_v, x_dvr_v, p_grid_v, p_dvr_v,
                  nx_grid, nx_dvr, nk_grid, nk_dvr,
                  (N == 1) ? &psi_grid : nullptr,
                  (N == 1) ? &psi_dvr : nullptr);
    }

    std::cout << "[step4-xc] wrote step4_crosscheck" << suf.str()
              << " (" << 3 * n_snap << " rows × per-quantity)\n";
}

} // namespace verify
} // namespace ecg1d
