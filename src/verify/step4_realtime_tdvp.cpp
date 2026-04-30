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
                                   bool rt_verbose,
                                   const Step4RealtimeOptions& rt_options) {
    std::cout << "\n=== Step 4 — ECG real-time TDVP (N=" << N << ", K=" << K
              << ", T=" << T_total << ", dt=" << dt << ") ===\n";

    HamiltonianTerms terms = HamiltonianTerms::kinetic_harmonic();
    PermutationSet perms   = PermutationSet::generate(N);
    std::vector<AlphaIndex> alpha = make_alpha_list(N, K);
    if (rt_options.u_split_trotter) {
        alpha.erase(std::remove_if(alpha.begin(), alpha.end(),
                                   [](const AlphaIndex& a) { return a.a1 == 1; }),
                    alpha.end());
    }

    SolverConfig solver_cfg = SolverConfig::dynamics();
    solver_cfg.lambda_C = rt_options.lambda_C;
    solver_cfg.rcond    = rt_options.rcond;
    solver_cfg.wiener_smooth = rt_options.wiener_smooth;

    std::cout << "[step4] solver lambda_C=" << solver_cfg.lambda_C
              << " rcond=" << solver_cfg.rcond
              << " wiener=" << (solver_cfg.wiener_smooth ? "on" : "off")
              << " enforce_norm=" << (rt_options.enforce_norm ? "on" : "off")
              << " u_split_trotter=" << (rt_options.u_split_trotter ? "on" : "off")
              << "\n";

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
    std::vector<double> xmraw_all, pmraw_all, x2raw_all, p2raw_all;
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
        rt_cfg.integrator   = rt_options.integrator;
        rt_cfg.sample_every = std::max(1, static_cast<int>(std::round(dt_trace / dt)));
        rt_cfg.verbose      = rt_verbose;     // CLI-controlled: --rt-verbose
        rt_cfg.print_every  = rt_verbose ? 1 : 0;
        rt_cfg.enforce_norm = rt_options.enforce_norm;
        rt_cfg.u_split_trotter = rt_options.u_split_trotter;

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
            t_all     = t_shift;
            E_all     = tr.E;
            norm_all  = tr.norm;
            xm_all    = tr.x_mean;
            pm_all    = tr.p_mean;
            x2_all    = tr.x2;
            p2_all    = tr.p2;
            xmraw_all = tr.x_mean_raw;
            pmraw_all = tr.p_mean_raw;
            x2raw_all = tr.x2_raw;
            p2raw_all = tr.p2_raw;
            first_leg = false;
        } else {
            append_skip_first(t_all,     t_shift);
            append_skip_first(E_all,     tr.E);
            append_skip_first(norm_all,  tr.norm);
            append_skip_first(xm_all,    tr.x_mean);
            append_skip_first(pm_all,    tr.p_mean);
            append_skip_first(x2_all,    tr.x2);
            append_skip_first(p2_all,    tr.p2);
            append_skip_first(xmraw_all, tr.x_mean_raw);
            append_skip_first(pmraw_all, tr.p_mean_raw);
            append_skip_first(x2raw_all, tr.x2_raw);
            append_skip_first(p2raw_all, tr.p2_raw);
        }

        R.t_snap(s) = t_b;
        record_snapshot(s, basis, last_E_seen);
    }

    // Pack trace
    const int n_tr = static_cast<int>(t_all.size());
    R.t_trace = Eigen::Map<Eigen::VectorXd>(t_all.data(), n_tr);
    R.E       = Eigen::Map<Eigen::VectorXd>(E_all.data(), n_tr);
    R.norm    = Eigen::Map<Eigen::VectorXd>(norm_all.data(), n_tr);
    R.x_mean     = Eigen::Map<Eigen::VectorXd>(xm_all.data(),    n_tr);
    R.p_mean     = Eigen::Map<Eigen::VectorXd>(pm_all.data(),    n_tr);
    R.x2         = Eigen::Map<Eigen::VectorXd>(x2_all.data(),    n_tr);
    R.p2         = Eigen::Map<Eigen::VectorXd>(p2_all.data(),    n_tr);
    R.x_mean_raw = Eigen::Map<Eigen::VectorXd>(xmraw_all.data(), n_tr);
    R.p_mean_raw = Eigen::Map<Eigen::VectorXd>(pmraw_all.data(), n_tr);
    R.x2_raw     = Eigen::Map<Eigen::VectorXd>(x2raw_all.data(), n_tr);
    R.p2_raw     = Eigen::Map<Eigen::VectorXd>(p2raw_all.data(), n_tr);

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

    // Raw <psi|H|psi> = E*norm should be flat to round-off under any Hermitian
    // projector, regardless of norm leak (Q9 / Q11 conservation law).
    if (n_tr > 0) {
        double Eraw0 = R.E(0) * R.norm(0);
        double maxdEraw = 0.0;
        for (int i = 0; i < n_tr; i++) {
            double Eraw_i = R.E(i) * R.norm(i);
            maxdEraw = std::max(maxdEraw, std::abs(Eraw_i - Eraw0));
        }
        double rel_raw = (std::abs(Eraw0) > 0) ? maxdEraw / std::abs(Eraw0) : maxdEraw;
        std::cout << "[step4] <psi|H|psi>_0=" << std::setprecision(10) << Eraw0
                  << "  max|d(E*norm)|=" << std::scientific << maxdEraw
                  << "  rel=" << rel_raw
                  << (rel_raw < 1e-10 ? "  [structural]" : "  [drift]")
                  << std::defaultfloat << "\n";
    }

    // Write CSVs
    std::ostringstream suf;
    suf << "_N" << N << "_K" << K << ".csv";

    // Trace. Column policy (preserves backward compat with plot_verification.py
    // which reads positional columns 3=x_mean, 4=p_mean for `normalized` and
    // `both`; `raw` writes NaN at columns 3-6 and is plotter-incompatible):
    //   normalized: t,E,norm,x_mean,p_mean,x2,p2                         (7 cols)
    //   both:       t,E,norm,x_mean,p_mean,x2,p2,x_mean_raw,p_mean_raw,x2_raw,p2_raw (11)
    //   raw:        t,E,norm,NaN,NaN,NaN,NaN,x_mean_raw,p_mean_raw,x2_raw,p2_raw     (11)
    {
        std::ofstream f(out_path("step4_ecg_trace" + suf.str()));
        const MomentForm mf = rt_options.moment_form;
        const bool emit_norm = (mf != MomentForm::Raw);
        const bool emit_raw  = (mf != MomentForm::Normalized);
        if (mf == MomentForm::Normalized) {
            f << "t,E,norm,x_mean,p_mean,x2,p2\n";
        } else if (mf == MomentForm::Raw) {
            f << "t,E,norm,x_mean,p_mean,x2,p2,x_mean_raw,p_mean_raw,x2_raw,p2_raw\n";
        } else {
            f << "t,E,norm,x_mean,p_mean,x2,p2,x_mean_raw,p_mean_raw,x2_raw,p2_raw\n";
        }
        f << std::setprecision(15);
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        for (int i = 0; i < n_tr; i++) {
            f << R.t_trace(i) << "," << R.E(i) << "," << R.norm(i);
            if (emit_norm) {
                f << "," << R.x_mean(i) << "," << R.p_mean(i)
                  << "," << R.x2(i)     << "," << R.p2(i);
            } else {
                f << "," << NaN << "," << NaN << "," << NaN << "," << NaN;
            }
            if (emit_raw) {
                f << "," << R.x_mean_raw(i) << "," << R.p_mean_raw(i)
                  << "," << R.x2_raw(i)     << "," << R.p2_raw(i);
            }
            f << "\n";
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
