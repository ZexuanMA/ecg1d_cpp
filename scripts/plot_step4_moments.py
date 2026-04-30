"""Plot the four moments <x>, <p>, <x^2>, <p^2> from step 4 traces.

Reads:
  out/verify/step4_ecg_trace_N{N}_K{K}.csv    (11 cols when moment-form=both)
  out/verify/step4_grid_trace_N{N}.csv
  out/verify/step4_dvr_trace_N{N}.csv

Writes:
  out/verify/fig_step4_moments_N{N}_K{K}.png      (2x2 panel: <x>, <p>, <x^2>, <p^2>)
  out/verify/fig_step4_norm_energy_N{N}_K{K}.png  (norm + E_norm + E*norm vs t)
  out/verify/fig_step4_raw_vs_norm_N{N}_K{K}.png  (ECG raw vs normalized x^2, p^2)

Usage:
  python scripts/plot_step4_moments.py [N=1] [K=3]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "out/verify"


def load_csv(path):
    arr = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    return arr


def main(N=1, K=3):
    ecg_p  = os.path.join(OUT, f"step4_ecg_trace_N{N}_K{K}.csv")
    grid_p = os.path.join(OUT, f"step4_grid_trace_N{N}.csv")
    dvr_p  = os.path.join(OUT, f"step4_dvr_trace_N{N}.csv")
    if not os.path.exists(ecg_p):
        print(f"missing {ecg_p}; run ./build/ecg1d_verify --run-step4 first", file=sys.stderr)
        sys.exit(1)

    e = load_csv(ecg_p)
    have_ref = os.path.exists(grid_p) and os.path.exists(dvr_p)
    g = load_csv(grid_p) if have_ref else None
    d = load_csv(dvr_p)  if have_ref else None
    if not have_ref:
        print(f"[note] no reference traces for N={N} (grid/DVR solvers are N=1 only); ECG-only plots.")

    # --- 2x2 panel: <x>, <p>, <x^2>, <p^2> ---
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    series = [(e, "ECG", "-", "C0", 1.6)]
    if have_ref:
        series += [
            (g, "FD grid",  "--", "C1", 1.3),
            (d, "sinc-DVR", ":",  "C2", 1.3),
        ]
    for arr, name, ls, c, lw in series:
        t = arr["t"]
        axes[0, 0].plot(t, arr["x_mean"], ls, color=c, label=name, lw=lw)
        axes[0, 1].plot(t, arr["p_mean"], ls, color=c, label=name, lw=lw)
        axes[1, 0].plot(t, arr["x2"],     ls, color=c, label=name, lw=lw)
        axes[1, 1].plot(t, arr["p2"],     ls, color=c, label=name, lw=lw)
    axes[0, 0].set_ylabel(r"$\langle x \rangle$")
    axes[0, 1].set_ylabel(r"$\langle p \rangle$")
    axes[1, 0].set_ylabel(r"$\langle x^2 \rangle$")
    axes[1, 1].set_ylabel(r"$\langle p^2 \rangle$")
    for ax in axes.flat:
        ax.set_xlabel(r"$t$")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle(rf"Step 4 — moments  (N={N}, K={K})  [normalized form]", fontsize=12)
    fig.tight_layout()
    out_moments = os.path.join(OUT, f"fig_step4_moments_N{N}_K{K}.png")
    fig.savefig(out_moments, dpi=130)
    plt.close(fig)
    print(f"wrote {out_moments}")

    # --- norm leak + energy diagnostic ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].plot(e["t"], e["norm"] / e["norm"][0], "-",  color="C0", label="ECG", lw=1.6)
    if have_ref:
        axes[0].plot(g["t"], g["norm"] / g["norm"][0], "--", color="C1", label="FD grid", lw=1.3)
        axes[0].plot(d["t"], d["norm"] / d["norm"][0], ":",  color="C2", label="sinc-DVR", lw=1.3)
    axes[0].set_xlabel(r"$t$"); axes[0].set_ylabel(r"$\langle\psi|\psi\rangle / \langle\psi(0)|\psi(0)\rangle$")
    axes[0].set_title("Norm (relative)")
    axes[0].axhline(1.0, color="k", lw=0.5, alpha=0.3)
    axes[0].grid(alpha=0.3); axes[0].legend(fontsize=8)

    # Normalized energy <H>/<psi|psi>
    for arr, name, ls, c, lw in series:
        rel = np.abs(arr["E"] - arr["E"][0]) / max(abs(arr["E"][0]), 1e-30)
        axes[1].semilogy(arr["t"], np.maximum(rel, 1e-18), ls, color=c, label=name, lw=lw)
    axes[1].set_xlabel(r"$t$"); axes[1].set_ylabel(r"$|E_{\rm norm}(t)-E_{\rm norm}(0)|/|E_{\rm norm}(0)|$")
    axes[1].set_title("Normalized energy drift")
    axes[1].grid(alpha=0.3, which="both"); axes[1].legend(fontsize=8)

    # Raw <psi|H|psi> = E_norm * norm — Q9 structural conservation under any Hermitian P
    Eraw_e = e["E"] * e["norm"]
    axes[2].plot(e["t"], Eraw_e / Eraw_e[0], "-",  color="C0", label="ECG", lw=1.6)
    if have_ref:
        Eraw_g = g["E"] * g["norm"]
        Eraw_d = d["E"] * d["norm"]
        axes[2].plot(g["t"], Eraw_g / Eraw_g[0], "--", color="C1", label="FD grid", lw=1.3)
        axes[2].plot(d["t"], Eraw_d / Eraw_d[0], ":",  color="C2", label="sinc-DVR", lw=1.3)
    axes[2].set_xlabel(r"$t$"); axes[2].set_ylabel(r"$\langle\psi|H|\psi\rangle(t)/\langle\psi|H|\psi\rangle(0)$")
    axes[2].set_title(r"Raw $\langle\psi|H|\psi\rangle$ — Q9 diagnostic")
    axes[2].axhline(1.0, color="k", lw=0.5, alpha=0.3)
    axes[2].grid(alpha=0.3); axes[2].legend(fontsize=8)

    fig.suptitle(rf"Step 4 — norm + energy  (N={N}, K={K})", fontsize=12)
    fig.tight_layout()
    out_ne = os.path.join(OUT, f"fig_step4_norm_energy_N{N}_K{K}.png")
    fig.savefig(out_ne, dpi=130)
    plt.close(fig)
    print(f"wrote {out_ne}")

    # --- ECG raw vs normalized comparison (only meaningful on ECG side) ---
    if "x2_raw" in e.dtype.names:
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
        axes[0, 0].plot(e["t"], e["x_mean"],     "-", color="C0", label=r"$\langle x\rangle$ norm", lw=1.6)
        axes[0, 0].plot(e["t"], e["x_mean_raw"], "--", color="C3", label=r"$\langle\psi|x|\psi\rangle$ raw", lw=1.3)
        axes[0, 1].plot(e["t"], e["p_mean"],     "-", color="C0", label=r"$\langle p\rangle$ norm", lw=1.6)
        axes[0, 1].plot(e["t"], e["p_mean_raw"], "--", color="C3", label=r"$\langle\psi|p|\psi\rangle$ raw", lw=1.3)
        axes[1, 0].plot(e["t"], e["x2"],     "-", color="C0", label=r"$\langle x^2\rangle$ norm", lw=1.6)
        axes[1, 0].plot(e["t"], e["x2_raw"], "--", color="C3", label=r"$\langle\psi|x^2|\psi\rangle$ raw", lw=1.3)
        axes[1, 1].plot(e["t"], e["p2"],     "-", color="C0", label=r"$\langle p^2\rangle$ norm", lw=1.6)
        axes[1, 1].plot(e["t"], e["p2_raw"], "--", color="C3", label=r"$\langle\psi|p^2|\psi\rangle$ raw", lw=1.3)
        for ax in axes.flat:
            ax.set_xlabel(r"$t$"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
        fig.suptitle(rf"Step 4 — ECG raw vs normalized moments  (N={N}, K={K})", fontsize=12)
        fig.tight_layout()
        out_rn = os.path.join(OUT, f"fig_step4_raw_vs_norm_N{N}_K{K}.png")
        fig.savefig(out_rn, dpi=130)
        plt.close(fig)
        print(f"wrote {out_rn}")
    else:
        print("skipping raw-vs-norm plot: ECG trace has no *_raw columns "
              "(rerun with --moment-form both)")


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    main(N, K)
