#include "common.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ecg1d {
namespace verify {

std::string out_dir() {
    return "out/verify";
}

std::string out_path(const std::string& filename) {
    std::filesystem::path p(out_dir());
    std::filesystem::create_directories(p);
    return (p / filename).string();
}

Eigen::VectorXd linspace(double a, double b, int n) {
    Eigen::VectorXd v(n);
    if (n == 1) { v(0) = a; return v; }
    double dx = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) v(i) = a + i * dx;
    return v;
}

static void write_cd(std::ostream& os, Cd z) {
    os << std::setprecision(17) << z.real() << "," << z.imag() << "\n";
}

void save_basis_csv(const std::string& path,
                    const std::vector<BasisParams>& basis) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    int K = static_cast<int>(basis.size());
    int N = basis[0].N();
    for (int k = 0; k < K; k++) {
        write_cd(f, basis[k].u);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                write_cd(f, basis[k].A(i, j));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                write_cd(f, basis[k].B(i, j));
        for (int i = 0; i < N; i++)
            write_cd(f, basis[k].R(i));
        f << basis[k].name << "\n";
    }
}

std::vector<BasisParams> load_basis_csv(const std::string& path, int N) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);

    auto read_cd = [&]() -> Cd {
        std::string line;
        if (!std::getline(f, line)) throw std::runtime_error("short CSV");
        std::stringstream ss(line);
        double re = 0, im = 0;
        char comma = ',';
        ss >> re >> comma >> im;
        return Cd(re, im);
    };
    auto read_int = [&]() -> int {
        std::string line;
        if (!std::getline(f, line)) throw std::runtime_error("short CSV");
        return std::stoi(line);
    };

    std::vector<BasisParams> basis;
    while (f.peek() != EOF) {
        BasisParams b;
        b.u = read_cd();
        b.A = MatrixXcd(N, N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) b.A(i, j) = read_cd();
        b.B = MatrixXcd(N, N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) b.B(i, j) = read_cd();
        b.R = VectorXcd(N);
        for (int i = 0; i < N; i++) b.R(i) = read_cd();
        b.name = read_int();
        basis.push_back(b);
        // skip possible trailing blanks
        while (f.peek() == '\n' || f.peek() == '\r') f.get();
    }
    return basis;
}

void write_two_column(const std::string& path,
                      const std::string& header,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    f << header << "\n";
    f << std::setprecision(15);
    for (int i = 0; i < x.size(); i++) f << x(i) << "," << y(i) << "\n";
}

void write_trace_csv(const std::string& path,
                     const std::string& header,
                     const std::vector<std::vector<double>>& columns) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("cannot open " + path);
    f << header << "\n";
    f << std::setprecision(15);
    if (columns.empty()) return;
    size_t n_rows = columns[0].size();
    for (size_t r = 0; r < n_rows; r++) {
        for (size_t c = 0; c < columns.size(); c++) {
            f << columns[c][r];
            if (c + 1 < columns.size()) f << ",";
        }
        f << "\n";
    }
}

} // namespace verify
} // namespace ecg1d
