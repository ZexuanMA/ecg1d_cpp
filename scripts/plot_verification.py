#!/usr/bin/env python3
"""Plot the output of ecg1d_verify. Reads out/verify/*.csv, writes PNGs into the
same directory. No CLI arguments; scans for every (N, K) combination present.

Run after:
    ./build/ecg1d_verify --N 1 --K 3
    ./build/ecg1d_verify --N 2 --K 5
    python scripts/plot_verification.py
"""
import os
import re
import sys
import glob
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "out/verify"


def load_key_value(path):
    d = {}
    with open(path) as f:
        rdr = csv.reader(f)
        next(rdr, None)  # header
        for row in rdr:
            if len(row) < 2 or not row[0]:
                continue
            try:
                d[row[0]] = float(row[1])
            except ValueError:
                d[row[0]] = row[1]
    return d


def load_two_col(path):
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    return arr[:, 0], arr[:, 1]


def discover():
    """Return dict: N -> {K, grid, dvr, step3, step4_trace, step4_snap}"""
    found = {}
    for p in glob.glob(os.path.join(OUT, "step1_ecg_energy_N*_K*.csv")):
        m = re.match(r".*step1_ecg_energy_N(\d+)_K(\d+)\.csv$", p)
        if not m: continue
        N, K = int(m.group(1)), int(m.group(2))
        found.setdefault(N, {})["K"] = K
    return found


def plot_energy_bars(runs):
    fig, ax = plt.subplots(figsize=(6, 4))
    labels, ecg_vals, grid_vals, dvr_vals = [], [], [], []
    for N, info in sorted(runs.items()):
        K = info["K"]
        step1 = load_key_value(os.path.join(OUT, f"step1_ecg_energy_N{N}_K{K}.csv"))
        step2g = load_key_value(os.path.join(OUT, f"step2_grid_energy_N{N}.csv"))
        step2d = load_key_value(os.path.join(OUT, f"step2_dvr_energy_N{N}.csv"))
        labels.append(f"N={N}")
        ecg_vals.append(step1["E_ECG"])
        grid_vals.append(step2g["E_ref"])
        dvr_vals.append(step2d["E_ref"])
    x = np.arange(len(labels))
    w = 0.27
    ax.bar(x - w, ecg_vals,  w, label="ECG")
    ax.bar(x,     grid_vals, w, label="FD grid")
    ax.bar(x + w, dvr_vals,  w, label="sinc-DVR")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("E_0")
    ax.set_title("Ground-state energy: ECG vs reference solvers")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step3_energy.png"), dpi=120)
    plt.close(fig)


def plot_density_overlay(runs):
    n = len(runs)
    fig, axes = plt.subplots(1, max(1, n), figsize=(5 * max(1, n), 4), squeeze=False)
    for i, (N, info) in enumerate(sorted(runs.items())):
        K = info["K"]
        ax = axes[0, i]
        x_e, n_e = load_two_col(os.path.join(OUT, f"step1_ecg_density_N{N}_K{K}.csv"))
        x_g, n_g = load_two_col(os.path.join(OUT, f"step2_grid_density_N{N}.csv"))
        x_d, n_d = load_two_col(os.path.join(OUT, f"step2_dvr_density_N{N}.csv"))
        ax.plot(x_e, n_e, "-",  label="ECG",      lw=2)
        ax.plot(x_g, n_g, "--", label="FD grid",  lw=1.2)
        ax.plot(x_d, n_d, ":",  label="sinc-DVR", lw=1.2)
        ax.set_xlabel("x"); ax.set_ylabel("n(x)")
        ax.set_title(f"N={N}, K={K}")
        ax.set_xlim(-6, 6)
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step3_density.png"), dpi=120)
    plt.close(fig)


def plot_momentum_overlay(runs):
    n = len(runs)
    fig, axes = plt.subplots(1, max(1, n), figsize=(5 * max(1, n), 4), squeeze=False)
    for i, (N, info) in enumerate(sorted(runs.items())):
        K = info["K"]
        ax = axes[0, i]
        k_e, nk_e = load_two_col(os.path.join(OUT, f"step1_ecg_nk_N{N}_K{K}.csv"))
        k_g, nk_g = load_two_col(os.path.join(OUT, f"step2_grid_nk_N{N}.csv"))
        k_d, nk_d = load_two_col(os.path.join(OUT, f"step2_dvr_nk_N{N}.csv"))
        ax.plot(k_e, nk_e, "-",  label="ECG",      lw=2)
        ax.plot(k_g, nk_g, "--", label="FD grid",  lw=1.2)
        ax.plot(k_d, nk_d, ":",  label="sinc-DVR", lw=1.2)
        ax.set_xlabel("k"); ax.set_ylabel("n(k)")
        ax.set_title(f"N={N}, K={K}")
        ax.set_xlim(-6, 6)
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step3_momentum.png"), dpi=120)
    plt.close(fig)


def _load_psi_csv(path):
    """Returns (x, psi_complex) or None if not found."""
    if not os.path.exists(path):
        return None
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    x = arr[:, 0]
    psi = arr[:, 1] + 1j * arr[:, 2]
    return x, psi


def _fix_global_phase(psi, ref):
    """Rotate psi so that <ref|psi> is real positive. Global phase only."""
    ov = np.vdot(ref, psi)
    if abs(ov) == 0:
        return psi
    return psi * (np.conj(ov) / abs(ov))


def plot_wavefunction_overlay(runs):
    """N=1 only: overlay Re(psi), Im(psi), |psi| for ECG, grid, DVR."""
    n1_runs = {N: info for N, info in runs.items() if N == 1}
    if not n1_runs:
        return
    for N, info in sorted(n1_runs.items()):
        K = info["K"]
        ecg = _load_psi_csv(os.path.join(OUT, f"step1_ecg_psi_N{N}_K{K}.csv"))
        grd = _load_psi_csv(os.path.join(OUT, f"step2_grid_psi_N{N}.csv"))
        dvr = _load_psi_csv(os.path.join(OUT, f"step2_dvr_psi_N{N}.csv"))
        if ecg is None or grd is None or dvr is None:
            continue

        # Fix global phase: make grid psi real-positive-peak, align others to it.
        xg, pg = grd
        # flip sign so peak is positive, imag ~ 0
        pg = pg * (1.0 if pg.real[np.argmax(np.abs(pg))] >= 0 else -1.0)
        # Normalize each to unit L2 (trap integration) on its own grid.
        def normalize(x, p):
            dx = x[1] - x[0]
            nrm = np.sqrt(np.sum(np.abs(p) ** 2) * dx)
            return p / nrm if nrm > 0 else p
        pg = normalize(xg, pg)
        xd, pd = dvr
        # interpolate pg onto xd for phase alignment
        pg_on_xd = np.interp(xd, xg, pg.real) + 1j * np.interp(xd, xg, pg.imag)
        pd = normalize(xd, _fix_global_phase(pd, pg_on_xd))
        xe, pe = ecg
        pg_on_xe = np.interp(xe, xg, pg.real) + 1j * np.interp(xe, xg, pg.imag)
        pe = normalize(xe, _fix_global_phase(pe, pg_on_xe))

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        for ax, field, label in zip(
            axes,
            [lambda p: p.real, lambda p: p.imag, lambda p: np.abs(p)],
            ["Re ψ(x)", "Im ψ(x)", "|ψ(x)|"]
        ):
            ax.plot(xe, field(pe), "-",  label="ECG",      lw=2.0)
            ax.plot(xg, field(pg), "--", label="FD grid",  lw=1.2)
            ax.plot(xd, field(pd), ":",  label="sinc-DVR", lw=1.2)
            ax.set_xlabel("x"); ax.set_ylabel(label)
            ax.set_xlim(-6, 6)
            ax.set_title(f"{label}  (N={N}, K={K})")
            ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"fig_step3_wavefunction_N{N}.png"), dpi=120)
        plt.close(fig)

        # ψ diff plot: |ψ_ECG - ψ_ref| on the common ECG grid.
        pg_on_xe2 = np.interp(xe, xg, pg.real) + 1j * np.interp(xe, xg, pg.imag)
        pd_on_xe2 = np.interp(xe, xd, pd.real) + 1j * np.interp(xe, xd, pd.imag)
        pe_aligned_grid = normalize(xe, _fix_global_phase(pe, pg_on_xe2))
        pe_aligned_dvr  = normalize(xe, _fix_global_phase(pe, pd_on_xe2))
        fig2, ax2 = plt.subplots(figsize=(6, 4))
        ax2.semilogy(xe, np.maximum(np.abs(pe_aligned_grid - pg_on_xe2), 1e-16),
                     label="|ψ_ECG - ψ_grid|")
        ax2.semilogy(xe, np.maximum(np.abs(pe_aligned_dvr - pd_on_xe2), 1e-16),
                     label="|ψ_ECG - ψ_dvr|")
        ax2.semilogy(xe, np.maximum(np.abs(pg_on_xe2 - pd_on_xe2), 1e-16),
                     label="|ψ_grid - ψ_dvr|")
        ax2.set_xlabel("x"); ax2.set_ylabel("point-wise |Δψ|")
        ax2.set_xlim(-6, 6)
        ax2.set_title(f"Wavefunction pointwise error (N={N}, K={K})")
        ax2.legend()
        fig2.tight_layout()
        fig2.savefig(os.path.join(OUT, f"fig_step3_wavefunction_error_N{N}.png"), dpi=120)
        plt.close(fig2)


def plot_step4_energy(runs):
    n = len(runs)
    fig, axes = plt.subplots(2, max(1, n), figsize=(6 * max(1, n), 6), squeeze=False)
    for i, (N, info) in enumerate(sorted(runs.items())):
        K = info["K"]
        path = os.path.join(OUT, f"step4_trace_N{N}_K{K}.csv")
        if not os.path.exists(path): continue
        data = np.loadtxt(path, delimiter=",", skiprows=1)
        t, E, nrm, x2, p2, fid = data.T
        axes[0, i].plot(t, E)
        axes[0, i].set_title(f"<H>(t), N={N}, K={K}")
        axes[0, i].set_xlabel("t"); axes[0, i].set_ylabel("<H>")
        dE = np.abs(E - E[0])
        axes[1, i].semilogy(t, np.maximum(dE, 1e-16))
        axes[1, i].set_title("|E(t) - E(0)|")
        axes[1, i].set_xlabel("t"); axes[1, i].set_ylabel("|dE|")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step4_energy_trace.png"), dpi=120)
    plt.close(fig)


def plot_step4_observables(runs):
    n = len(runs)
    fig, axes = plt.subplots(1, max(1, n), figsize=(6 * max(1, n), 4), squeeze=False)
    for i, (N, info) in enumerate(sorted(runs.items())):
        K = info["K"]
        path = os.path.join(OUT, f"step4_trace_N{N}_K{K}.csv")
        if not os.path.exists(path): continue
        data = np.loadtxt(path, delimiter=",", skiprows=1)
        t, E, nrm, x2, p2, fid = data.T
        ax = axes[0, i]
        ax.plot(t, x2, label="<x²>")
        ax.plot(t, p2, label="<p²>")
        ax.plot(t, fid, label="|<ψ₀|ψ(t)>|²")
        ax.set_xlabel("t")
        ax.set_title(f"N={N}, K={K}")
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step4_observables.png"), dpi=120)
    plt.close(fig)


def plot_step4_snapshots(runs):
    n = len(runs)
    fig, axes = plt.subplots(2, max(1, n), figsize=(6 * max(1, n), 8), squeeze=False)
    for i, (N, info) in enumerate(sorted(runs.items())):
        K = info["K"]
        path = os.path.join(OUT, f"step4_snapshots_N{N}_K{K}.csv")
        if not os.path.exists(path): continue
        with open(path) as f:
            rdr = csv.reader(f); next(rdr)
            rows = [row for row in rdr]
        times = sorted({float(r[0]) for r in rows})
        for kind, row_ax in [("x", axes[0, i]), ("k", axes[1, i])]:
            for t in times:
                sel = [r for r in rows if r[1] == kind and float(r[0]) == t]
                g = np.array([float(r[2]) for r in sel])
                v = np.array([float(r[3]) for r in sel])
                order = np.argsort(g)
                row_ax.plot(g[order], v[order], label=f"t={t:.2f}")
            row_ax.set_xlabel(kind); row_ax.set_ylabel(f"n({kind},t)")
            row_ax.set_title(f"N={N}, K={K}, n({kind},t)")
            row_ax.legend(fontsize=8)
            if kind == "x": row_ax.set_xlim(-6, 6)
            if kind == "k": row_ax.set_xlim(-6, 6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step4_snapshots.png"), dpi=120)
    plt.close(fig)


def main():
    if not os.path.isdir(OUT):
        print(f"no output directory {OUT}; run ecg1d_verify first", file=sys.stderr)
        sys.exit(1)
    runs = discover()
    if not runs:
        print("no Step-1 results found; run ecg1d_verify first", file=sys.stderr)
        sys.exit(1)

    print(f"found runs: {runs}")
    plot_energy_bars(runs)
    plot_density_overlay(runs)
    plot_momentum_overlay(runs)
    plot_wavefunction_overlay(runs)
    plot_step4_energy(runs)
    plot_step4_observables(runs)
    plot_step4_snapshots(runs)
    print(f"wrote figures to {OUT}/fig_*.png")


if __name__ == "__main__":
    main()
