#include "step4_realtime_tdvp.hpp"
#include "common.hpp"
#include "ecg_wavefunction.hpp"
#include "realtime_tdvp.hpp"
#include "observables.hpp"
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

static std::vector<AlphaIndex> make_alpha_list(int N, int K) {
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

void step4_realtime_tdvp(int N, int K,
                         const std::vector<BasisParams>& basis_init,
                         double T_total, double dt,
                         const Eigen::VectorXd& x_grid,
                         const Eigen::VectorXd& k_grid,
                         int n_snapshots) {
    std::cout << "\n=== Step 4 — real-time TDVP (N=" << N << ", K=" << K
              << ", T=" << T_total << ", dt=" << dt << ") ===\n";

    HamiltonianTerms terms = HamiltonianTerms::kinetic_harmonic();
    PermutationSet perms   = PermutationSet::generate(N);
    std::vector<AlphaIndex> alpha = make_alpha_list(N, K);

    SolverConfig solver_cfg = SolverConfig::dynamics();
    solver_cfg.lambda_C = 1e-8;
    solver_cfg.rcond    = 1e-4;
    solver_cfg.wiener_smooth = true;

    RealtimeEvolutionConfig rt_cfg;
    rt_cfg.dt           = dt;
    rt_cfg.integrator   = RtIntegrator::RK4;
    rt_cfg.sample_every = std::max(1, static_cast<int>(std::round(0.05 / dt)));
    rt_cfg.verbose      = false;
    rt_cfg.print_every  = 0;
    rt_cfg.enforce_norm = false;

    std::cout << "[step4] running trace-only evolution (sample_every="
              << rt_cfg.sample_every << ")...\n";
    RealtimeEvolutionResult rt = realtime_tdvp_evolution(
        alpha, basis_init, T_total, terms, solver_cfg, rt_cfg);

    const RealtimeTrace& tr = rt.trace;
    std::cout << "[step4] trace samples: " << tr.t.size() << "\n";

    // Energy conservation check
    double E0 = tr.E.front();
    double Emax_dev = 0.0;
    for (double e : tr.E) Emax_dev = std::max(Emax_dev, std::abs(e - E0));
    double rel_drift = (std::abs(E0) > 0) ? Emax_dev / std::abs(E0) : Emax_dev;
    std::cout << "[step4] E0=" << std::setprecision(10) << E0
              << "  max|dE|=" << std::scientific << Emax_dev
              << "  rel=" << rel_drift
              << (rel_drift < 1e-4 ? "  [PASS]" : "  [WARN]")
              << std::defaultfloat << "\n";

    // Write trace
    std::ostringstream suffix;
    suffix << "_N" << N << "_K" << K << ".csv";

    {
        std::ofstream f(out_path("step4_trace" + suffix.str()));
        f << "t,E,norm,x2,p2,fidelity\n";
        f << std::setprecision(15);
        for (size_t i = 0; i < tr.t.size(); i++) {
            double fid = 0.0;
            if (i < tr.overlap0.size() && tr.norm.size() > 0) {
                double n0 = tr.norm.front();
                double nt = tr.norm[i];
                double denom = n0 * nt;
                fid = (denom > 0) ? std::norm(tr.overlap0[i]) / denom : 0.0;
            }
            f << tr.t[i] << "," << tr.E[i] << "," << tr.norm[i]
              << "," << tr.x2[i] << "," << tr.p2[i]
              << "," << fid << "\n";
        }
    }

    // Snapshots: to save basis at intermediate times, re-run evolution to each
    // snapshot time independently. Cost is modest for small N, K.
    std::vector<double> snap_times;
    for (int s = 0; s < n_snapshots; s++) {
        double t_s = T_total * static_cast<double>(s) / (n_snapshots - 1);
        snap_times.push_back(t_s);
    }

    {
        std::ofstream f(out_path("step4_snapshots" + suffix.str()));
        f << "t,kind,grid,value\n";
        f << std::setprecision(15);

        RealtimeEvolutionConfig snap_cfg = rt_cfg;
        snap_cfg.sample_every = std::numeric_limits<int>::max(); // no trace overhead
        snap_cfg.verbose      = false;

        for (double t_s : snap_times) {
            std::vector<BasisParams> b;
            if (t_s <= 0.0) {
                b = basis_init;
            } else {
                RealtimeEvolutionResult r = realtime_tdvp_evolution(
                    alpha, basis_init, t_s, terms, solver_cfg, snap_cfg);
                b = r.basis_final;
            }
            Eigen::VectorXd nx = ecg_single_particle_density(b, perms, x_grid, 0);
            Eigen::VectorXd nk = momentum_distribution(b, perms, k_grid, 0);
            for (int j = 0; j < x_grid.size(); j++)
                f << t_s << ",x," << x_grid(j) << "," << nx(j) << "\n";
            for (int j = 0; j < k_grid.size(); j++)
                f << t_s << ",k," << k_grid(j) << "," << nk(j) << "\n";
            std::cout << "[step4] snapshot t=" << std::setprecision(4) << t_s
                      << " written\n";
        }
    }
}

} // namespace verify
} // namespace ecg1d
