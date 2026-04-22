#pragma once
#include <cmath>

namespace ecg1d {

// Portable cylindrical Bessel function J_n(x) for integer order n.
// Motivation: std::cyl_bessel_j is C++17 special-math and is NOT in libc++
// (macOS default). This ascending-series implementation covers our regime
// (|x| up to ~30, |n| up to ~30) with double precision.
//
// J_n(x) = sum_{m=0}^inf (-1)^m / (m! (n+m)!) * (x/2)^(2m+n)
// J_{-n}(x) = (-1)^n J_n(x)
//
// Series converges rapidly for |x/2|^2 / (m*(n+m)) < 1, i.e. small x or large m.
// For x = 1 and n<=10, ~15 terms suffice.
inline double bessel_j(int n, double x) {
    if (n < 0) {
        double v = bessel_j(-n, x);
        return (((-n) & 1) ? -v : v);
    }
    const double half_x = 0.5 * x;
    // Leading term: (x/2)^n / n!
    double term = 1.0;
    for (int k = 1; k <= n; ++k) term *= half_x / static_cast<double>(k);
    double sum = term;
    // Recurrence: term_{m+1} = term_m * -(x/2)^2 / ((m+1)(n+m+1))
    const double minus_half_x_sq = -half_x * half_x;
    for (int m = 0; m < 500; ++m) {
        term *= minus_half_x_sq / (static_cast<double>(m + 1) * static_cast<double>(n + m + 1));
        sum += term;
        // Termination: |term| < 1e-18 * |sum| gives ~15 correct digits
        if (std::abs(term) < 1e-18 * std::abs(sum)) break;
    }
    return sum;
}

} // namespace ecg1d
