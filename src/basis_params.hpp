#pragma once
#include "types.hpp"
#include <optional>

namespace ecg1d {

struct BasisParams {
    Cd u;
    MatrixXcd A;  // (N,N)
    MatrixXcd B;  // (N,N) diagonal
    VectorXcd R;  // (N,)
    int name;

    int N() const { return static_cast<int>(R.size()); }

    BasisParams conj_params() const;

    static BasisParams from_arrays(Cd u, const MatrixXcd& A, const MatrixXcd& B,
                                   const VectorXcd& R, int name = 0);

    static BasisParams randoman(int N, std::optional<int> seed = std::nullopt,
                                int name = 0, double margin = 1e-6,
                                double take_diag_fraction = 0.34);
};

} // namespace ecg1d
