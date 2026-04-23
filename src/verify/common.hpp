#pragma once
#include "types.hpp"
#include "basis_params.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace ecg1d {
namespace verify {

std::string out_dir();
std::string out_path(const std::string& filename);

Eigen::VectorXd linspace(double a, double b, int n);

void save_basis_csv(const std::string& path,
                    const std::vector<BasisParams>& basis);

std::vector<BasisParams> load_basis_csv(const std::string& path, int N);

void write_two_column(const std::string& path,
                      const std::string& header,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y);

void write_trace_csv(const std::string& path,
                     const std::string& header,
                     const std::vector<std::vector<double>>& columns);

} // namespace verify
} // namespace ecg1d
