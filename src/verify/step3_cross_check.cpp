#include "step3_cross_check.hpp"
#include "common.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ecg1d {
namespace verify {

struct PairStats {
    std::string label;
    double dE_abs, dE_rel;
    double L2_nx, L2_nk;
    double overlap_nx;     // Bhattacharyya on n(x): sum sqrt(n1 n2) dx
    double fidelity_psi;   // |<psi_1|psi_2>|^2 / (||psi_1||^2 ||psi_2||^2), NaN if not computed
};

// Linearly interpolate samples (x_src, y_src) onto x_tgt. x_src is assumed sorted.
static Eigen::VectorXd interp1d(const Eigen::VectorXd& x_src,
                                const Eigen::VectorXd& y_src,
                                const Eigen::VectorXd& x_tgt) {
    Eigen::VectorXd out(x_tgt.size());
    int n = x_src.size();
    for (int i = 0; i < x_tgt.size(); i++) {
        double xt = x_tgt(i);
        if (xt <= x_src(0))        { out(i) = y_src(0);     continue; }
        if (xt >= x_src(n - 1))    { out(i) = y_src(n - 1); continue; }
        // binary search
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x_src(mid) <= xt) lo = mid; else hi = mid;
        }
        double t = (xt - x_src(lo)) / (x_src(hi) - x_src(lo));
        out(i) = (1 - t) * y_src(lo) + t * y_src(hi);
    }
    return out;
}

static double trapz(const Eigen::VectorXd& y, const Eigen::VectorXd& x) {
    double s = 0;
    for (int i = 0; i + 1 < x.size(); i++) {
        s += 0.5 * (y(i) + y(i + 1)) * (x(i + 1) - x(i));
    }
    return s;
}

// Linear interpolation of a complex-valued function.
static Eigen::VectorXcd interp1d_complex(const Eigen::VectorXd& x_src,
                                         const Eigen::VectorXcd& y_src,
                                         const Eigen::VectorXd& x_tgt) {
    Eigen::VectorXcd out(x_tgt.size());
    int n = x_src.size();
    for (int i = 0; i < x_tgt.size(); i++) {
        double xt = x_tgt(i);
        if (xt <= x_src(0))     { out(i) = y_src(0);     continue; }
        if (xt >= x_src(n - 1)) { out(i) = y_src(n - 1); continue; }
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x_src(mid) <= xt) lo = mid; else hi = mid;
        }
        double t = (xt - x_src(lo)) / (x_src(hi) - x_src(lo));
        out(i) = (1 - t) * y_src(lo) + t * y_src(hi);
    }
    return out;
}

// Complex wavefunction fidelity: |<psi_1|psi_2>|^2 / (||psi_1||^2 * ||psi_2||^2)
// using trapezoidal quadrature on x_tgt. Does not require psi to be normalized.
static double wavefunction_fidelity(const Eigen::VectorXd& x1, const Eigen::VectorXcd& p1,
                                    const Eigen::VectorXd& x2, const Eigen::VectorXcd& p2,
                                    const Eigen::VectorXd& x_tgt) {
    Eigen::VectorXcd a = interp1d_complex(x1, p1, x_tgt);
    Eigen::VectorXcd b = interp1d_complex(x2, p2, x_tgt);
    int n = x_tgt.size();
    Cd     ov(0.0, 0.0);
    double na = 0.0, nb = 0.0;
    for (int i = 0; i + 1 < n; i++) {
        double dx = x_tgt(i + 1) - x_tgt(i);
        Cd     contribA = std::conj(a(i))     * b(i);
        Cd     contribB = std::conj(a(i + 1)) * b(i + 1);
        ov += 0.5 * (contribA + contribB) * dx;
        na += 0.5 * (std::norm(a(i)) + std::norm(a(i + 1))) * dx;
        nb += 0.5 * (std::norm(b(i)) + std::norm(b(i + 1))) * dx;
    }
    double denom = na * nb;
    if (denom <= 0) return std::numeric_limits<double>::quiet_NaN();
    return std::norm(ov) / denom;
}

// Resample n1 (on grid x1) to x_tgt, compare L2 difference and Bhattacharyya overlap.
static void compare_on_grid(const Eigen::VectorXd& x1, const Eigen::VectorXd& n1,
                            const Eigen::VectorXd& x2, const Eigen::VectorXd& n2,
                            const Eigen::VectorXd& x_tgt,
                            double& l2, double& bh) {
    Eigen::VectorXd a = interp1d(x1, n1, x_tgt);
    Eigen::VectorXd b = interp1d(x2, n2, x_tgt);
    Eigen::VectorXd d = (a - b).array().square();
    l2 = std::sqrt(trapz(d, x_tgt));
    Eigen::VectorXd prod(a.size());
    for (int i = 0; i < a.size(); i++) {
        double p = std::max(a(i), 0.0) * std::max(b(i), 0.0);
        prod(i) = std::sqrt(p);
    }
    bh = trapz(prod, x_tgt);
}

void step3_cross_check(int N,
                       const Step1Result& ecg_res,
                       const Eigen::VectorXd& ecg_x_grid,
                       const Eigen::VectorXd& ecg_n_x,
                       const Eigen::VectorXd& ecg_k_grid,
                       const Eigen::VectorXd& ecg_n_k,
                       const Eigen::VectorXcd& ecg_psi_x,
                       const RefSolverResult& grid_res,
                       const RefSolverResult& dvr_res) {
    std::cout << "\n=== Step 3 — cross-check (N=" << N << ") ===\n";

    bool have_psi = (ecg_psi_x.size() == ecg_x_grid.size());

    // For density comparison, choose a common x-grid: use ecg_x_grid (usually densest/widest).
    // For momentum: use ecg_k_grid (which is what ECG and both ref solvers were sampled on).
    std::vector<PairStats> rows;
    auto add_pair = [&](const std::string& label,
                        double E1, double E2,
                        const Eigen::VectorXd& x1, const Eigen::VectorXd& nx1,
                        const Eigen::VectorXd& x2, const Eigen::VectorXd& nx2,
                        const Eigen::VectorXd& kgrid,
                        const Eigen::VectorXd& nk1, const Eigen::VectorXd& nk2,
                        const Eigen::VectorXcd* psi1, const Eigen::VectorXd* xpsi1,
                        const Eigen::VectorXcd* psi2, const Eigen::VectorXd* xpsi2) {
        PairStats s;
        s.label = label;
        s.dE_abs = E1 - E2;
        s.dE_rel = (std::abs(E2) > 0) ? s.dE_abs / E2 : 0.0;
        double l2x, bhx;
        compare_on_grid(x1, nx1, x2, nx2, ecg_x_grid, l2x, bhx);
        s.L2_nx = l2x;
        s.overlap_nx = bhx;
        // momentum: both on kgrid (same), direct compare
        Eigen::VectorXd dk = nk1 - nk2;
        double l2k_sq = 0.0;
        for (int i = 0; i + 1 < kgrid.size(); i++) {
            double d1 = dk(i) * dk(i);
            double d2 = dk(i + 1) * dk(i + 1);
            l2k_sq += 0.5 * (d1 + d2) * (kgrid(i + 1) - kgrid(i));
        }
        s.L2_nk = std::sqrt(l2k_sq);
        s.fidelity_psi = std::numeric_limits<double>::quiet_NaN();
        if (psi1 && psi2 && xpsi1 && xpsi2 &&
            psi1->size() == xpsi1->size() && psi2->size() == xpsi2->size()) {
            s.fidelity_psi = wavefunction_fidelity(*xpsi1, *psi1, *xpsi2, *psi2, ecg_x_grid);
        }
        rows.push_back(s);
    };

    const Eigen::VectorXcd* p_ecg_psi = have_psi ? &ecg_psi_x      : nullptr;
    const Eigen::VectorXd*  p_ecg_x   = have_psi ? &ecg_x_grid     : nullptr;
    const Eigen::VectorXcd* p_grid_psi = (N == 1) ? &grid_res.psi  : nullptr;
    const Eigen::VectorXd*  p_grid_x   = (N == 1) ? &grid_res.x    : nullptr;
    const Eigen::VectorXcd* p_dvr_psi  = (N == 1) ? &dvr_res.psi   : nullptr;
    const Eigen::VectorXd*  p_dvr_x    = (N == 1) ? &dvr_res.x     : nullptr;

    add_pair("ECG_vs_grid", ecg_res.E_ecg, grid_res.E,
             ecg_x_grid, ecg_n_x, grid_res.x, grid_res.n_x,
             ecg_k_grid, ecg_n_k, grid_res.n_k,
             p_ecg_psi, p_ecg_x, p_grid_psi, p_grid_x);
    add_pair("ECG_vs_dvr",  ecg_res.E_ecg, dvr_res.E,
             ecg_x_grid, ecg_n_x, dvr_res.x, dvr_res.n_x,
             ecg_k_grid, ecg_n_k, dvr_res.n_k,
             p_ecg_psi, p_ecg_x, p_dvr_psi, p_dvr_x);
    add_pair("grid_vs_dvr", grid_res.E,    dvr_res.E,
             grid_res.x, grid_res.n_x, dvr_res.x, dvr_res.n_x,
             ecg_k_grid, grid_res.n_k, dvr_res.n_k,
             p_grid_psi, p_grid_x, p_dvr_psi, p_dvr_x);

    // Print + write CSV
    std::cout << std::setw(14) << "pair"
              << std::setw(16) << "dE_abs"
              << std::setw(16) << "dE_rel"
              << std::setw(14) << "L2_n(x)"
              << std::setw(14) << "L2_n(k)"
              << std::setw(14) << "BH_n(x)"
              << std::setw(16) << "|<psi1|psi2>|^2" << "\n";
    std::cout << std::string(104, '-') << "\n";
    for (const auto& r : rows) {
        std::cout << std::setw(14) << r.label
                  << std::setw(16) << std::scientific << std::setprecision(4) << r.dE_abs
                  << std::setw(16) << r.dE_rel
                  << std::setw(14) << r.L2_nx
                  << std::setw(14) << r.L2_nk
                  << std::setw(14) << std::fixed << std::setprecision(8) << r.overlap_nx
                  << std::setw(16);
        if (std::isnan(r.fidelity_psi)) std::cout << "-";
        else std::cout << std::fixed << std::setprecision(10) << r.fidelity_psi;
        std::cout << "\n";
    }
    std::cout << std::defaultfloat;

    std::ostringstream suffix;
    suffix << "_N" << N << ".csv";
    std::ofstream f(out_path("step3_crosscheck" + suffix.str()));
    f << "pair,dE_abs,dE_rel,L2_n_x,L2_n_k,overlap_BH_n_x,fidelity_psi\n";
    f << std::setprecision(15);
    for (const auto& r : rows) {
        f << r.label << "," << r.dE_abs << "," << r.dE_rel
          << "," << r.L2_nx << "," << r.L2_nk << "," << r.overlap_nx
          << "," << r.fidelity_psi << "\n";
    }
    // Also log each individual energy for easy reading in Python.
    f << "E_ECG," << ecg_res.E_ecg << ",,,,,\n";
    f << "E_grid," << grid_res.E   << ",,,,,\n";
    f << "E_dvr," << dvr_res.E     << ",,,,,\n";
}

} // namespace verify
} // namespace ecg1d
