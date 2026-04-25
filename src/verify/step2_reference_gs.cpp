#include "step2_reference_gs.hpp"
#include "common.hpp"
#include "physical_constants.hpp"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace ecg1d {
namespace verify {

namespace detail {

Eigen::MatrixXd build_T_grid_fd(const Eigen::VectorXd& x) {
    const int n = static_cast<int>(x.size());
    const double dx = x(1) - x(0);
    const double coeff = 1.0 / (2.0 * mass * dx * dx);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n, n);
    for (int j = 0; j < n; j++) {
        T(j, j) = 2.0 * coeff;
        if (j > 0)         T(j, j - 1) = -coeff;
        if (j < n - 1)     T(j, j + 1) = -coeff;
    }
    return T;
}

Eigen::MatrixXd build_T_dvr_sinc(const Eigen::VectorXd& x) {
    const int n = static_cast<int>(x.size());
    const double dx = x(1) - x(0);
    const double inv_mdx2 = 1.0 / (mass * dx * dx);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                T(i, j) = (M_PI * M_PI / 6.0) * inv_mdx2;
            } else {
                double sign = ((i - j) % 2 == 0) ? 1.0 : -1.0;
                double denom = static_cast<double>(i - j);
                T(i, j) = sign * inv_mdx2 / (denom * denom);
            }
        }
    }
    return T;
}

Eigen::VectorXd build_V_full(const Eigen::VectorXd& x) {
    Eigen::VectorXd V(x.size());
    for (int j = 0; j < x.size(); j++) {
        V(j) = 0.5 * mass * omega * omega * x(j) * x(j)
             + kappa * std::cos(2.0 * k_L * x(j));
    }
    return V;
}

Eigen::VectorXd build_V_ho_only(const Eigen::VectorXd& x) {
    Eigen::VectorXd V(x.size());
    for (int j = 0; j < x.size(); j++) {
        V(j) = 0.5 * mass * omega * omega * x(j) * x(j);
    }
    return V;
}

void normalize_grid_state(Eigen::VectorXcd& psi, double dx) {
    double s = 0.0;
    for (int j = 0; j < psi.size(); j++) s += std::norm(psi(j));
    s *= dx;
    if (s > 0) psi /= std::sqrt(s);
}

Eigen::VectorXd compute_nk(const Eigen::VectorXcd& psi,
                           const Eigen::VectorXd& x,
                           const Eigen::VectorXd& k_grid) {
    double dx = x(1) - x(0);
    Eigen::VectorXd nk(k_grid.size());
    double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
    for (int ki = 0; ki < k_grid.size(); ki++) {
        Cd acc(0.0, 0.0);
        double kv = k_grid(ki);
        for (int j = 0; j < x.size(); j++) {
            Cd phase = std::exp(Cd(0.0, -kv * x(j)));
            acc += psi(j) * phase;
        }
        acc *= dx * inv_sqrt_2pi;
        nk(ki) = std::norm(acc);
    }
    return nk;
}

} // namespace detail

RefSolverResult reference_gs_grid(int N, int n_grid, double x_max,
                                  const Eigen::VectorXd& k_grid) {
    std::cout << "[step2] grid diag: N=" << N << ", n_grid=" << n_grid
              << ", x_max=" << x_max << "\n";

    double dx = 2.0 * x_max / (n_grid - 1);
    Eigen::VectorXd x(n_grid);
    for (int j = 0; j < n_grid; j++) x(j) = -x_max + j * dx;

    Eigen::MatrixXd H = detail::build_T_grid_fd(x);
    Eigen::VectorXd V = detail::build_V_full(x);
    for (int j = 0; j < n_grid; j++) H(j, j) += V(j);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    double E_sp = es.eigenvalues()(0);
    Eigen::VectorXcd psi = es.eigenvectors().col(0).cast<Cd>();
    detail::normalize_grid_state(psi, dx);

    RefSolverResult r;
    r.E_sp = E_sp;
    r.E   = static_cast<double>(N) * E_sp;
    r.x   = x;
    r.psi = psi;
    r.n_x = Eigen::VectorXd(n_grid);
    for (int j = 0; j < n_grid; j++) r.n_x(j) = std::norm(psi(j));
    r.k   = k_grid;
    r.n_k = detail::compute_nk(psi, x, k_grid);

    std::cout << "        E_sp=" << std::setprecision(12) << E_sp
              << "  E_total=" << r.E << "\n";
    return r;
}

RefSolverResult reference_gs_dvr(int N, int n_dvr, double x_max,
                                 const Eigen::VectorXd& k_grid) {
    std::cout << "[step2] sinc-DVR: N=" << N << ", n_dvr=" << n_dvr
              << ", x_max=" << x_max << "\n";

    double dx = 2.0 * x_max / (n_dvr - 1);
    Eigen::VectorXd x(n_dvr);
    for (int j = 0; j < n_dvr; j++) x(j) = -x_max + j * dx;

    Eigen::MatrixXd H = detail::build_T_dvr_sinc(x);
    Eigen::VectorXd V = detail::build_V_full(x);
    for (int j = 0; j < n_dvr; j++) H(j, j) += V(j);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    double E_sp = es.eigenvalues()(0);
    Eigen::VectorXcd psi = es.eigenvectors().col(0).cast<Cd>();
    detail::normalize_grid_state(psi, dx);

    RefSolverResult r;
    r.E_sp = E_sp;
    r.E   = static_cast<double>(N) * E_sp;
    r.x   = x;
    r.psi = psi;
    r.n_x = Eigen::VectorXd(n_dvr);
    for (int j = 0; j < n_dvr; j++) r.n_x(j) = std::norm(psi(j));
    r.k   = k_grid;
    r.n_k = detail::compute_nk(psi, x, k_grid);

    std::cout << "        E_sp=" << std::setprecision(12) << E_sp
              << "  E_total=" << r.E << "\n";
    return r;
}

void write_reference_csvs(const std::string& method_tag, int N,
                          const RefSolverResult& res) {
    std::ostringstream suffix;
    suffix << "_N" << N << ".csv";

    {
        std::ofstream f(out_path("step2_" + method_tag + "_energy" + suffix.str()));
        f << "quantity,value\n";
        f << std::setprecision(15);
        f << "N," << N << "\n";
        f << "E_ref," << res.E << "\n";
        f << "E_single_particle," << res.E_sp << "\n";
    }
    write_two_column(out_path("step2_" + method_tag + "_density" + suffix.str()),
                     "x,n_x", res.x, res.n_x);
    write_two_column(out_path("step2_" + method_tag + "_nk" + suffix.str()),
                     "k,n_k", res.k, res.n_k);
    // Save real-space wavefunction (complex) so step3 can compute grid<->DVR overlap
    std::ofstream fpsi(out_path("step2_" + method_tag + "_psi" + suffix.str()));
    fpsi << "x,Re_psi,Im_psi\n";
    fpsi << std::setprecision(15);
    for (int j = 0; j < res.x.size(); j++) {
        fpsi << res.x(j) << "," << res.psi(j).real() << "," << res.psi(j).imag() << "\n";
    }
}

} // namespace verify
} // namespace ecg1d
