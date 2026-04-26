#include "step4_realtime_tdvp.hpp"
#include "common.hpp"
#include "ecg_wavefunction.hpp"
#include "realtime_tdvp.hpp"
#include "observables.hpp"
#include "pair_cache.hpp"
#include "permutation.hpp"
#include "tdvp_solver.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace ecg1d {
namespace verify {

namespace {

std::vector<AlphaIndex> make_alpha_list(int N, int K) {
    std::vector<AlphaIndex> a;
    for (int i = 0; i < K; i++) a.push_back({1, i, 0, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) a.push_back({2, i, j, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++) a.push_back({3, i, j, 0});
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
            for (int k = j; k < N; k++) a.push_back({4, i, j, k});
    return a;
}

// Append `src` onto `dst`, skipping its first element (which duplicates the
// previous leg's final sample).
template <class Vec>
void append_skip_first(Vec& dst, const Vec& src) {
    if (src.size() <= 1) return;
    dst.insert(dst.end(), src.begin() + 1, src.end());
}

// Compute <bra|ket> (unnormalized) by summing PairCache.M_G across (i, j, perm),
// using the same convention as sample_observables in realtime_tdvp.cpp.
Cd ecg_basis_overlap(const std::vector<BasisParams>& bra,
                     const std::vector<BasisParams>& ket,
                     const PermutationSet& perms) {
    Cd amp(0.0, 0.0);
    const int Kb = static_cast<int>(bra.size());
    const int Kk = static_cast<int>(ket.size());
    for (int i = 0; i < Kb; ++i) {
        Cd conj_ui = std::conj(bra[i].u);
        for (int j = 0; j < Kk; ++j) {
            Cd sum_p(0.0, 0.0);
            for (int p = 0; p < perms.SN; ++p) {
                PairCache c = PairCache::build(bra[i], ket[j], perms.matrices[p]);
                sum_p += static_cast<double>(perms.signs[p]) * c.M_G;
            }
            amp += conj_ui * ket[j].u * sum_p;
        }
    }
    return amp;
}

} // anonymous namespace

Step4EcgResult step4_realtime_tdvp(int N, int K,
                                   const std::vector<BasisParams>& basis_init,
                                   double T_total, double dt,
                                   const Eigen::VectorXd& x_grid,
                                   const Eigen::VectorXd& k_grid,
                                   const std::vector<double>& t_snap_in,
                                   double dt_trace,
                                   bool rt_verbose) {
    std::cout << "\n=== Step 4 — ECG real-time TDVP (N=" << N << ", K=" << K
              << ", T=" << T_total << ", dt=" << dt << ") ===\n";

    HamiltonianTerms terms = HamiltonianTerms::kinetic_harmonic();
    PermutationSet perms   = PermutationSet::generate(N);
    std::vector<AlphaIndex> alpha = make_alpha_list(N, K);

    SolverConfig solver_cfg = SolverConfig::dynamics();
    solver_cfg.lambda_C = 1e-8;
    solver_cfg.rcond    = 1e-4;
    solver_cfg.wiener_smooth = true;

    // Sort, dedup, and bracket snapshot times by 0 and T_total.
    std::vector<double> t_snap = t_snap_in;
    if (t_snap.empty() || t_snap.front() > 1e-12) t_snap.insert(t_snap.begin(), 0.0);
    std::sort(t_snap.begin(), t_snap.end());
    t_snap.erase(std::unique(t_snap.begin(), t_snap.end(),
                             [](double a, double b) { return std::abs(a - b) < 1e-12; }),
                 t_snap.end());
    if (t_snap.back() < T_total - 1e-12) t_snap.push_back(T_total);
    const int n_snap = static_cast<int>(t_snap.size());

    Step4EcgResult R;
    R.x = x_grid;
    R.k = k_grid;
    R.t_snap = Eigen::VectorXd(n_snap);
    R.density_snaps = Eigen::MatrixXd(n_snap, x_grid.size());
    R.nk_snaps      = Eigen::MatrixXd(n_snap, k_grid.size());
    R.E_snap        = Eigen::VectorXd(n_snap);
    R.fidelity_snap = Eigen::VectorXd(n_snap);
    R.basis_snaps.reserve(n_snap);

    auto record_snapshot = [&](int idx, const std::vector<BasisParams>& b, double leg_E_at_b) {
        Eigen::VectorXd nx = ecg_single_particle_density(b, perms, x_grid, 0);
        Eigen::VectorXd nk;
        if (N == 1) {
            nk = ecg_momentum_density_1p(b, k_grid);    // |psi_tilde(k)|^2 (correct)
        } else {
            nk = momentum_density_from_density(x_grid, nx, k_grid);
        }
        for (int j = 0; j < x_grid.size(); j++) R.density_snaps(idx, j) = nx(j);
        for (int kj = 0; kj < k_grid.size(); kj++) R.nk_snaps(idx, kj) = nk(kj);
        R.E_snap(idx) = leg_E_at_b;
        // Global fidelity: |<basis_init|b>|^2 / (norm_init * norm_b)
        Cd amp = ecg_basis_overlap(basis_init, b, perms);
        Cd s_init = ecg_basis_overlap(basis_init, basis_init, perms);
        Cd s_t    = ecg_basis_overlap(b, b, perms);
        double denom = (s_init * s_t).real();
        R.fidelity_snap(idx) = (denom > 0) ? std::norm(amp) / denom : 0.0;
        R.basis_snaps.push_back(b);
    };

    std::vector<BasisParams> basis = basis_init;

    // Trace accumulators
    std::vector<double> t_all, E_all, norm_all, xm_all, pm_all, x2_all, p2_all;
    bool first_leg = true;

    // Initial snapshot
    {
        // E_init from a one-step "leg" of length zero is awkward; compute it
        // directly via the same observable pathway used inside realtime_tdvp.
        // Cheapest: do a tiny no-op evolution to get a one-sample trace,
        // OR we can wait for the first leg's t=0 sample (which IS basis_init).
    }
    R.t_snap(0) = 0.0;

    double last_E_seen = 0.0;
    int next_snap_to_record = 0;  // record snapshot 0 once we have its E

    for (int s = 1; s < n_snap; s++) {
        double t_a = t_snap[s - 1];
        double t_b = t_snap[s];
        double dT  = t_b - t_a;
        if (dT <= 1e-14) continue;

        RealtimeEvolutionConfig rt_cfg;
        rt_cfg.dt           = dt;
        rt_cfg.integrator   = RtIntegrator::RK4;
        rt_cfg.sample_every = std::max(1, static_cast<int>(std::round(dt_trace / dt)));
        rt_cfg.verbose      = rt_verbose;     // CLI-controlled: --rt-verbose
        rt_cfg.print_every  = rt_verbose ? 1 : 0;
        rt_cfg.enforce_norm = false;

        std::cout << "[step4] leg t=[" << t_a << ", " << t_b << "]\n";
        RealtimeEvolutionResult rt = realtime_tdvp_evolution(
            alpha, basis, dT, terms, solver_cfg, rt_cfg);
        const RealtimeTrace& tr = rt.trace;

        // Record initial snapshot now that we have its E (from the first leg's
        // t=0 sample).
        if (next_snap_to_record == 0 && !tr.E.empty()) {
            record_snapshot(0, basis, tr.E.front());
            next_snap_to_record = 1;
        }

        basis = rt.basis_final;
        last_E_seen = tr.E.empty() ? last_E_seen : tr.E.back();

        // Shift trace times to absolute, append.
        std::vector<double> t_shift = tr.t;
        for (auto& v : t_shift) v += t_a;
        if (first_leg) {
            t_all   = t_shift;
            E_all   = tr.E;
            norm_all= tr.norm;
            xm_all  = tr.x_mean;
            pm_all  = tr.p_mean;
            x2_all  = tr.x2;
            p2_all  = tr.p2;
            first_leg = false;
        } else {
            append_skip_first(t_all,    t_shift);
            append_skip_first(E_all,    tr.E);
            append_skip_first(norm_all, tr.norm);
            append_skip_first(xm_all,   tr.x_mean);
            append_skip_first(pm_all,   tr.p_mean);
            append_skip_first(x2_all,   tr.x2);
            append_skip_first(p2_all,   tr.p2);
        }

        R.t_snap(s) = t_b;
        record_snapshot(s, basis, last_E_seen);
    }

    // Pack trace
    const int n_tr = static_cast<int>(t_all.size());
    R.t_trace = Eigen::Map<Eigen::VectorXd>(t_all.data(), n_tr);
    R.E       = Eigen::Map<Eigen::VectorXd>(E_all.data(), n_tr);
    R.norm    = Eigen::Map<Eigen::VectorXd>(norm_all.data(), n_tr);
    R.x_mean  = Eigen::Map<Eigen::VectorXd>(xm_all.data(), n_tr);
    R.p_mean  = Eigen::Map<Eigen::VectorXd>(pm_all.data(), n_tr);
    R.x2      = Eigen::Map<Eigen::VectorXd>(x2_all.data(), n_tr);
    R.p2      = Eigen::Map<Eigen::VectorXd>(p2_all.data(), n_tr);

    // Energy drift summary
    R.E_init = (n_tr > 0) ? R.E(0) : 0.0;
    double maxdE = 0.0;
    for (int i = 0; i < n_tr; i++) maxdE = std::max(maxdE, std::abs(R.E(i) - R.E_init));
    R.max_rel_E_drift = (std::abs(R.E_init) > 0) ? maxdE / std::abs(R.E_init) : maxdE;
    std::cout << "[step4] E0=" << std::setprecision(10) << R.E_init
              << "  max|dE|=" << std::scientific << maxdE
              << "  rel=" << R.max_rel_E_drift
              << (R.max_rel_E_drift < 1e-4 ? "  [PASS]" : "  [WARN]")
              << std::defaultfloat << "\n";

    // Write CSVs
    std::ostringstream suf;
    suf << "_N" << N << "_K" << K << ".csv";

    // Trace
    {
        std::ofstream f(out_path("step4_ecg_trace" + suf.str()));
        f << "t,E,norm,x_mean,p_mean\n";
        f << std::setprecision(15);
        for (int i = 0; i < n_tr; i++) {
            f << R.t_trace(i) << "," << R.E(i) << "," << R.norm(i) << ","
              << R.x_mean(i) << "," << R.p_mean(i) << "\n";
        }
    }
    // Density snapshots (long format)
    {
        std::ofstream f(out_path("step4_ecg_density" + suf.str()));
        f << "t,x,n_x\n";
        f << std::setprecision(15);
        for (int s = 0; s < n_snap; s++) {
            for (int j = 0; j < x_grid.size(); j++) {
                f << R.t_snap(s) << "," << x_grid(j) << ","
                  << R.density_snaps(s, j) << "\n";
            }
        }
    }
    // Momentum-density snapshots (long format)
    {
        std::ofstream f(out_path("step4_ecg_nk" + suf.str()));
        f << "t,k,n_k\n";
        f << std::setprecision(15);
        for (int s = 0; s < n_snap; s++) {
            for (int kj = 0; kj < k_grid.size(); kj++) {
                f << R.t_snap(s) << "," << k_grid(kj) << ","
                  << R.nk_snaps(s, kj) << "\n";
            }
        }
    }
    // Snapshot meta (one row per snapshot: E, fidelity)
    {
        std::ofstream f(out_path("step4_ecg_snap" + suf.str()));
        f << "t,E,fidelity\n";
        f << std::setprecision(15);
        for (int s = 0; s < n_snap; s++) {
            f << R.t_snap(s) << "," << R.E_snap(s) << "," << R.fidelity_snap(s) << "\n";
        }
    }

    std::cout << "[step4] wrote " << n_tr << " trace rows, " << n_snap << " snapshots\n";
    return R;
}

} // namespace verify
} // namespace ecg1d
