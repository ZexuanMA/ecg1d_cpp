#!/usr/bin/env python3
"""Plot the output of ecg1d_verify. Reads out/verify/*.csv, writes PNGs into the
same directory. No CLI arguments; scans for every (N, K) combination present.

Every figure is written per-N (filename ends in _N{N}.png) so running a new N
does NOT overwrite the figures from a previous N.

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
    """Return dict: N -> {K}"""
    found = {}
    for p in glob.glob(os.path.join(OUT, "step1_ecg_energy_N*_K*.csv")):
        m = re.match(r".*step1_ecg_energy_N(\d+)_K(\d+)\.csv$", p)
        if not m: continue
        N, K = int(m.group(1)), int(m.group(2))
        found.setdefault(N, {})["K"] = K
    return found


def plot_energy_bars_single(N, K):
    step1 = load_key_value(os.path.join(OUT, f"step1_ecg_energy_N{N}_K{K}.csv"))
    step2g = load_key_value(os.path.join(OUT, f"step2_grid_energy_N{N}.csv"))
    step2d = load_key_value(os.path.join(OUT, f"step2_dvr_energy_N{N}.csv"))
    fig, ax = plt.subplots(figsize=(5, 4))
    methods = ["ECG", "FD grid", "sinc-DVR"]
    vals    = [step1["E_ECG"], step2g["E_ref"], step2d["E_ref"]]
    ax.bar(methods, vals)
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.8f}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("E_0")
    ax.set_title(f"Ground-state energy (N={N}, K={K})")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step3_energy_N{N}.png"), dpi=120)
    plt.close(fig)


def plot_density_overlay_single(N, K):
    x_e, n_e = load_two_col(os.path.join(OUT, f"step1_ecg_density_N{N}_K{K}.csv"))
    x_g, n_g = load_two_col(os.path.join(OUT, f"step2_grid_density_N{N}.csv"))
    x_d, n_d = load_two_col(os.path.join(OUT, f"step2_dvr_density_N{N}.csv"))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x_e, n_e, "-",  label="ECG",      lw=2)
    ax.plot(x_g, n_g, "--", label="FD grid",  lw=1.2)
    ax.plot(x_d, n_d, ":",  label="sinc-DVR", lw=1.2)
    ax.set_xlabel("x"); ax.set_ylabel("n(x)")
    ax.set_title(f"Position density n(x)  (N={N}, K={K})")
    ax.set_xlim(-6, 6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step3_density_N{N}.png"), dpi=120)
    plt.close(fig)


def plot_momentum_overlay_single(N, K):
    k_e, nk_e = load_two_col(os.path.join(OUT, f"step1_ecg_nk_N{N}_K{K}.csv"))
    k_g, nk_g = load_two_col(os.path.join(OUT, f"step2_grid_nk_N{N}.csv"))
    k_d, nk_d = load_two_col(os.path.join(OUT, f"step2_dvr_nk_N{N}.csv"))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(k_e, nk_e, "-",  label="ECG",      lw=2)
    ax.plot(k_g, nk_g, "--", label="FD grid",  lw=1.2)
    ax.plot(k_d, nk_d, ":",  label="sinc-DVR", lw=1.2)
    ax.set_xlabel("k"); ax.set_ylabel("n(k)")
    ax.set_title(f"Momentum distribution n(k)  (N={N}, K={K})")
    ax.set_xlim(-6, 6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step3_momentum_N{N}.png"), dpi=120)
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


def _load_ecg_psi(N, K):
    """Load the ECG single-particle orbital.
      N=1: complex psi from step1_ecg_psi_N1_K{K}.csv (written by step1).
      N>=2: build phi(x) = sqrt(n(x)) from the density (exact for non-
            interacting bosonic ground states); Im(phi) is identically 0.
    Returns (x, psi_complex) or None if the inputs are missing.
    """
    direct = os.path.join(OUT, f"step1_ecg_psi_N{N}_K{K}.csv")
    if os.path.exists(direct):
        return _load_psi_csv(direct)
    dens = os.path.join(OUT, f"step1_ecg_density_N{N}_K{K}.csv")
    if not os.path.exists(dens):
        return None
    arr = np.loadtxt(dens, delimiter=",", skiprows=1)
    x = arr[:, 0]
    phi = np.sqrt(np.maximum(arr[:, 1], 0.0)).astype(complex)
    return x, phi


def plot_wavefunction_single(N, K):
    """Overlay Re(psi), Im(psi), |psi| for ECG, FD grid, sinc-DVR and
    write pointwise |Δψ|. For N=1 uses the complex ECG psi written by
    step1; for N>=2 uses phi = sqrt(n(x)) (valid for non-interacting GS)."""
    ecg = _load_ecg_psi(N, K)
    grd = _load_psi_csv(os.path.join(OUT, f"step2_grid_psi_N{N}.csv"))
    dvr = _load_psi_csv(os.path.join(OUT, f"step2_dvr_psi_N{N}.csv"))
    if ecg is None or grd is None or dvr is None:
        return

    xg, pg = grd
    pg = pg * (1.0 if pg.real[np.argmax(np.abs(pg))] >= 0 else -1.0)

    def normalize(x, p):
        dx = x[1] - x[0]
        nrm = np.sqrt(np.sum(np.abs(p) ** 2) * dx)
        return p / nrm if nrm > 0 else p

    pg = normalize(xg, pg)
    xd, pd = dvr
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


def plot_step4_energy_single(N, K):
    path = os.path.join(OUT, f"step4_trace_N{N}_K{K}.csv")
    if not os.path.exists(path):
        return
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    t, E, nrm, x2, p2, fid = data.T
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(t, E)
    axes[0].set_title(f"<H>(t)  (N={N}, K={K})")
    axes[0].set_xlabel("t"); axes[0].set_ylabel("<H>")
    dE = np.abs(E - E[0])
    axes[1].semilogy(t, np.maximum(dE, 1e-16))
    axes[1].set_title(f"|E(t) - E(0)|  (N={N}, K={K})")
    axes[1].set_xlabel("t"); axes[1].set_ylabel("|dE|")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_energy_trace_N{N}.png"), dpi=120)
    plt.close(fig)


def plot_step4_observables_single(N, K):
    path = os.path.join(OUT, f"step4_trace_N{N}_K{K}.csv")
    if not os.path.exists(path):
        return
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    t, E, nrm, x2, p2, fid = data.T
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, x2, label="<x²>")
    ax.plot(t, p2, label="<p²>")
    ax.plot(t, fid, label="|<ψ₀|ψ(t)>|²")
    ax.set_xlabel("t")
    ax.set_title(f"Real-time observables (N={N}, K={K})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_observables_N{N}.png"), dpi=120)
    plt.close(fig)


def plot_step4_snapshots_single(N, K):
    path = os.path.join(OUT, f"step4_snapshots_N{N}_K{K}.csv")
    if not os.path.exists(path):
        return
    with open(path) as f:
        rdr = csv.reader(f); next(rdr)
        rows = [row for row in rdr]
    times = sorted({float(r[0]) for r in rows})
    fig, axes = plt.subplots(2, 1, figsize=(7, 8))
    for kind, ax in [("x", axes[0]), ("k", axes[1])]:
        for t in times:
            sel = [r for r in rows if r[1] == kind and float(r[0]) == t]
            g = np.array([float(r[2]) for r in sel])
            v = np.array([float(r[3]) for r in sel])
            order = np.argsort(g)
            ax.plot(g[order], v[order], label=f"t={t:.2f}")
        ax.set_xlabel(kind); ax.set_ylabel(f"n({kind},t)")
        ax.set_title(f"n({kind},t)  (N={N}, K={K})")
        ax.legend(fontsize=8)
        ax.set_xlim(-6, 6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_snapshots_N{N}.png"), dpi=120)
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
    for N, info in sorted(runs.items()):
        K = info["K"]
        plot_energy_bars_single(N, K)
        plot_density_overlay_single(N, K)
        plot_momentum_overlay_single(N, K)
        plot_wavefunction_single(N, K)
        plot_step4_energy_single(N, K)
        plot_step4_observables_single(N, K)
        plot_step4_snapshots_single(N, K)
    print(f"wrote figures to {OUT}/fig_*_N*.png")


if __name__ == "__main__":
    main()
