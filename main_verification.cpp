// main_verification.cpp — end-to-end verification of ECG+TDVP against reference solvers.
//
// Steps:
//   1. Imaginary-time TDVP ground state of H_full = T + V_ho + kappa*cos(2 k_L x) in ECG.
//   2. Reference ground states: FD grid + sinc-DVR on the same single-particle H.
//      (N=2 is handled via product state since there is no inter-particle interaction.)
//   3. Cross-check: energy, n(x), n(k) pairwise comparison.
//   4. Real-time TDVP with cos potential turned off. Check energy conservation.
//
// Usage:
//   ./build/ecg1d_verify --N 1 --K 3
//   ./build/ecg1d_verify --N 2 --K 5 --run-step4
//   ./build/ecg1d_verify --help

#include "verify/common.hpp"
#include "verify/ecg_wavefunction.hpp"
#include "verify/step1_imag_time_ecg.hpp"
#include "verify/step2_reference_gs.hpp"
#include "verify/step3_cross_check.hpp"
#include "verify/step4_realtime_tdvp.hpp"
#include "verify/step4_reference_dynamics.hpp"
#include "verify/step4_cross_check.hpp"
#include "observables.hpp"
#include "permutation.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace ecg1d;
using namespace ecg1d::verify;

struct Args {
    int    N       = 1;
    int    K       = 0;       // 0 -> default (N=1:3, N=2:5)
    bool   run_step1 = false;
    bool   run_step2 = false;
    bool   run_step3 = false;
    bool   run_step4 = false;
    double dt      = 1e-4;
    double T_total = 2.0 * M_PI;   // one HO period (omega = 1)
    double dt_trace = 5e-2;
    int    n_grid  = 512;
    double x_max   = 15.0;
    int    n_dvr   = 256;
    int    n_snap  = 11;
    bool   rt_verbose = false;
    double rt_lambda_C = 1e-8;
    double rt_rcond = 1e-4;
    bool   rt_wiener_smooth = true;
    bool   rt_enforce_norm = false;
    bool   rt_u_split_trotter = false;
    bool   help    = false;
};

static void print_help() {
    std::cout <<
R"(ecg1d_verify - verification harness for ECG + TDVP.

Options:
  --N <int>         particle number (1 or 2)                    [default 1]
  --K <int>         ECG basis size                              [default 3 for N=1, 5 for N=2]
  --run-step1       run Step 1 (ECG imaginary-time ground state)
  --run-step2       run Step 2 (FD grid + sinc-DVR reference)
  --run-step3       run Step 3 (cross-check; requires step1+2 in same session)
  --run-step4       run Step 4 (real-time TDVP; requires step1 basis)
  (if no --run-stepN flag is given, ALL four steps run)

  --dt <float>      real-time step size                         [default 1e-4]
  --T  <float>      total real-time evolution time              [default 2*pi]
  --dt-trace <float> trace sampling interval                    [default 5e-2]
  --n-grid <int>    FD grid size                                [default 512]
  --x-max <float>   FD/DVR grid half-width                      [default 15.0]
  --n-dvr <int>     sinc-DVR grid size                          [default 256]
  --n-snap <int>    number of real-time density snapshots       [default 11]
  --rt-verbose      print Step-4 ECG E(t) at every RK4 step      [default off]
  --enforce-norm    rescale u after each realtime step           [default off]
  --u-split-trotter update geometry by TDVP, then evolve u       [default off]
  --rcond <float>   realtime pseudoinverse rcond                 [default 1e-4]
  --lambda-C <float> realtime C-matrix Tikhonov lambda           [default 1e-8]
  --wiener          use smooth Wiener pseudoinverse filter       [default on]
  --no-wiener       use hard rcond pseudoinverse truncation
  --help            show this message
)";
}

static Args parse(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string f = argv[i];
        auto next_int   = [&]() { return std::atoi(argv[++i]); };
        auto next_float = [&]() { return std::atof(argv[++i]); };
        if      (f == "--N")          a.N = next_int();
        else if (f == "--K")          a.K = next_int();
        else if (f == "--run-step1")  a.run_step1 = true;
        else if (f == "--run-step2")  a.run_step2 = true;
        else if (f == "--run-step3")  a.run_step3 = true;
        else if (f == "--run-step4")  a.run_step4 = true;
        else if (f == "--dt")         a.dt = next_float();
        else if (f == "--T")          a.T_total = next_float();
        else if (f == "--dt-trace")   a.dt_trace = next_float();
        else if (f == "--n-grid")     a.n_grid = next_int();
        else if (f == "--x-max")      a.x_max = next_float();
        else if (f == "--n-dvr")      a.n_dvr = next_int();
        else if (f == "--n-snap")     a.n_snap = next_int();
        else if (f == "--rt-verbose") a.rt_verbose = true;
        else if (f == "--enforce-norm") a.rt_enforce_norm = true;
        else if (f == "--u-split-trotter") a.rt_u_split_trotter = true;
        else if (f == "--rcond")      a.rt_rcond = next_float();
        else if (f == "--lambda-C")   a.rt_lambda_C = next_float();
        else if (f == "--wiener")     a.rt_wiener_smooth = true;
        else if (f == "--no-wiener")  a.rt_wiener_smooth = false;
        else if (f == "--help" || f == "-h") a.help = true;
        else {
            std::cerr << "unknown flag: " << f << "\n";
            a.help = true;
        }
    }
    if (!(a.run_step1 || a.run_step2 || a.run_step3 || a.run_step4)) {
        a.run_step1 = a.run_step2 = a.run_step3 = a.run_step4 = true;
    }
    if (a.K == 0) a.K = (a.N == 1) ? 3 : 5;
    return a;
}

int main_inner(int argc, char** argv);

int main(int argc, char** argv) {
    try {
        return main_inner(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << std::endl;
        return 2;
    } catch (...) {
        std::cerr << "FATAL: unknown exception" << std::endl;
        return 3;
    }
}

int main_inner(int argc, char** argv) {
    Args args = parse(argc, argv);
    if (args.help) { print_help(); return 0; }

    std::cout << "=== ecg1d_verify (N=" << args.N << ", K=" << args.K << ") ===\n";
    std::cout << "output directory: " << out_dir() << "\n";

    // Shared grids
    Eigen::VectorXd x_grid = linspace(-args.x_max, args.x_max, 401);
    Eigen::VectorXd k_grid = linspace(-10.0, 10.0, 401);

    // Step 1
    Step1Result step1;
    Eigen::VectorXd ecg_n_x, ecg_n_k;
    Eigen::VectorXcd ecg_psi_x;                       // N=1 only; empty otherwise
    bool have_step1 = false;
    if (args.run_step1) {
        step1 = step1_imag_time_ecg(args.N, args.K, x_grid, k_grid);
        PermutationSet perms = PermutationSet::generate(args.N);
        ecg_n_x = ecg_single_particle_density(step1.basis, perms, x_grid, 0);
        // Momentum probability density n(k) = |phi_tilde(k)|^2.
        //   N=1: analytic FT of each basis function.
        //   N>=2 non-interacting bosonic GS: phi = sqrt(n(x)), then FT.
        if (args.N == 1) {
            ecg_psi_x = ecg_wavefunction_1p(step1.basis, x_grid);
            ecg_n_k   = ecg_momentum_density_1p(step1.basis, k_grid);
        } else {
            ecg_n_k   = momentum_density_from_density(x_grid, ecg_n_x, k_grid);
        }
        have_step1 = true;
    }

    // Step 2
    RefSolverResult grid_res, dvr_res;
    bool have_step2 = false;
    if (args.run_step2) {
        std::cout << "\n=== Step 2 — reference ground states (N=" << args.N << ") ===\n";
        grid_res = reference_gs_grid(args.N, args.n_grid, args.x_max, k_grid);
        dvr_res  = reference_gs_dvr (args.N, args.n_dvr,  args.x_max, k_grid);
        write_reference_csvs("grid", args.N, grid_res);
        write_reference_csvs("dvr",  args.N, dvr_res);
        have_step2 = true;
    }

    // Step 3
    if (args.run_step3) {
        if (!(have_step1 && have_step2)) {
            std::cerr << "[step3] requires Step 1 AND Step 2 in the same run; skipping.\n";
        } else {
            step3_cross_check(args.N, step1, x_grid, ecg_n_x, k_grid, ecg_n_k,
                              ecg_psi_x, grid_res, dvr_res);
        }
    }

    // Step 4
    if (args.run_step4) {
        std::vector<BasisParams> basis_for_rt;
        if (have_step1) {
            basis_for_rt = step1.basis;
        } else {
            std::ostringstream p;
            p << "step1_ecg_basis_N" << args.N << "_K" << args.K << ".csv";
            std::string path = out_path(p.str());
            std::cout << "\n[step4] loading cached Step-1 basis from " << path << "\n";
            basis_for_rt = load_basis_csv(path, args.N);
            if (static_cast<int>(basis_for_rt.size()) != args.K) {
                std::cerr << "[step4] warning: basis size " << basis_for_rt.size()
                          << " != requested K=" << args.K << "\n";
            }
        }

        // Snapshot times: evenly spaced in [0, T_total], n_snap points.
        std::vector<double> t_snap(args.n_snap);
        for (int s = 0; s < args.n_snap; s++) {
            t_snap[s] = args.T_total * static_cast<double>(s)
                        / static_cast<double>(args.n_snap - 1);
        }
        // Trace times: dense, every dt_trace.
        std::vector<double> t_trace_ref;
        for (double t = 0.0; t <= args.T_total + 1e-12; t += args.dt_trace) {
            t_trace_ref.push_back(std::min(t, args.T_total));
        }

        // ECG side
        Step4RealtimeOptions rt_options;
        rt_options.lambda_C = args.rt_lambda_C;
        rt_options.rcond = args.rt_rcond;
        rt_options.wiener_smooth = args.rt_wiener_smooth;
        rt_options.enforce_norm = args.rt_enforce_norm;
        rt_options.u_split_trotter = args.rt_u_split_trotter;

        Step4EcgResult ecg = step4_realtime_tdvp(
            args.N, args.K, basis_for_rt, args.T_total, args.dt,
            x_grid, k_grid, t_snap, args.dt_trace, args.rt_verbose, rt_options);

        // N=1 only: the reference dynamics solvers below propagate a single-
        // particle wavefunction. For N>=2 we'd need a different reference.
        if (args.N == 1) {
            std::cout << "\n=== Step 4 — reference dynamics (FD grid) ===\n";
            Step4RefResult grid = step4_realtime_grid(
                basis_for_rt, args.n_grid, args.x_max, k_grid,
                args.T_total, t_snap, t_trace_ref);
            write_step4_ref_csvs("grid", args.N, grid);

            std::cout << "\n=== Step 4 — reference dynamics (sinc-DVR) ===\n";
            Step4RefResult dvr = step4_realtime_dvr(
                basis_for_rt, args.n_dvr, args.x_max, k_grid,
                args.T_total, t_snap, t_trace_ref);
            write_step4_ref_csvs("dvr", args.N, dvr);

            step4_cross_check(args.N, ecg, grid, dvr, x_grid, k_grid);
        } else {
            std::cout << "[step4] reference dynamics solvers are N=1 only; "
                      << "skipping for N=" << args.N << ".\n";
        }
    }

    std::cout << "\n=== done ===\n";
    return 0;
}
