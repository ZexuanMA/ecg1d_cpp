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


# --- Step 4 (real-time TDVP) plotting ---------------------------------------
# CSV layout produced by ecg1d_verify --run-step4:
#   step4_ecg_trace_N{N}_K{K}.csv     : t,E,norm,x_mean,p_mean
#   step4_ecg_density_N{N}_K{K}.csv   : t,x,n_x  (long format)
#   step4_ecg_nk_N{N}_K{K}.csv        : t,k,n_k  (long format)
#   step4_ecg_snap_N{N}_K{K}.csv      : t,E,fidelity   (one row per snapshot)
#   step4_grid_trace_N{N}.csv         : t,E,norm,x_mean,p_mean,fidelity
#   step4_grid_density_N{N}.csv       : t,x,n_x  (long)
#   step4_grid_nk_N{N}.csv            : t,k,n_k  (long)
#   step4_dvr_*                       : same as grid_, on DVR discretization
#   step4_crosscheck_N{N}.csv         : pair,t,quantity,value (long)


def _step4_paths(N, K):
    return {
        "ecg_trace":  os.path.join(OUT, f"step4_ecg_trace_N{N}_K{K}.csv"),
        "ecg_snap":   os.path.join(OUT, f"step4_ecg_snap_N{N}_K{K}.csv"),
        "ecg_dens":   os.path.join(OUT, f"step4_ecg_density_N{N}_K{K}.csv"),
        "ecg_nk":     os.path.join(OUT, f"step4_ecg_nk_N{N}_K{K}.csv"),
        "grid_trace": os.path.join(OUT, f"step4_grid_trace_N{N}.csv"),
        "grid_dens":  os.path.join(OUT, f"step4_grid_density_N{N}.csv"),
        "grid_nk":    os.path.join(OUT, f"step4_grid_nk_N{N}.csv"),
        "dvr_trace":  os.path.join(OUT, f"step4_dvr_trace_N{N}.csv"),
        "dvr_dens":   os.path.join(OUT, f"step4_dvr_density_N{N}.csv"),
        "dvr_nk":     os.path.join(OUT, f"step4_dvr_nk_N{N}.csv"),
        "xc":         os.path.join(OUT, f"step4_crosscheck_N{N}.csv"),
    }


def _have(p):
    return os.path.exists(p)


def _load_long_grid(path, value_col):
    """Read a 't,grid,value' long-format CSV. Returns (times, grid, M[t,grid])."""
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    ts_all = arr[:, 0]
    g_all  = arr[:, 1]
    v_all  = arr[:, 2]
    times = np.unique(ts_all)
    grid  = np.unique(g_all)
    M = np.zeros((len(times), len(grid)))
    t_idx = {t: i for i, t in enumerate(times)}
    g_idx = {g: i for i, g in enumerate(grid)}
    for ti, gi, vi in zip(ts_all, g_all, v_all):
        M[t_idx[ti], g_idx[gi]] = vi
    return times, grid, M


def plot_step4_energy_drift_single(N, K):
    p = _step4_paths(N, K)
    if not _have(p["ecg_trace"]) or not _have(p["grid_trace"]) or not _have(p["dvr_trace"]):
        return
    e = np.loadtxt(p["ecg_trace"],  delimiter=",", skiprows=1)
    g = np.loadtxt(p["grid_trace"], delimiter=",", skiprows=1)
    d = np.loadtxt(p["dvr_trace"],  delimiter=",", skiprows=1)
    fig, ax = plt.subplots(figsize=(7, 4))
    for arr, name, ls, c in [
        (e, "ECG",      "-",  "C0"),
        (g, "FD grid",  "--", "C1"),
        (d, "sinc-DVR", ":",  "C2"),
    ]:
        t = arr[:, 0]; E = arr[:, 1]
        E0 = E[0]
        denom = max(abs(E0), 1e-30)
        rel = np.abs(E - E0) / denom
        ax.semilogy(t, np.maximum(rel, 1e-18), ls, color=c, label=name, lw=1.5)
    ax.set_xlabel("t"); ax.set_ylabel("|E(t) - E(0)| / |E(0)|")
    ax.set_title(f"Step 4 — energy-drift  (N={N}, K={K})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_energy_drift_N{N}.png"), dpi=120)
    plt.close(fig)


def plot_step4_observables_single(N, K):
    p = _step4_paths(N, K)
    if not _have(p["ecg_trace"]) or not _have(p["grid_trace"]) or not _have(p["dvr_trace"]):
        return
    e = np.loadtxt(p["ecg_trace"],  delimiter=",", skiprows=1)  # t,E,norm,x,p
    g = np.loadtxt(p["grid_trace"], delimiter=",", skiprows=1)  # +fidelity
    d = np.loadtxt(p["dvr_trace"],  delimiter=",", skiprows=1)
    snap_e = np.loadtxt(p["ecg_snap"], delimiter=",", skiprows=1) if _have(p["ecg_snap"]) else None

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for arr, name, ls, c in [
        (e, "ECG",      "-",  "C0"),
        (g, "FD grid",  "--", "C1"),
        (d, "sinc-DVR", ":",  "C2"),
    ]:
        t = arr[:, 0]
        x_mean = arr[:, 3]
        p_mean = arr[:, 4]
        axes[0].plot(t, x_mean, ls, color=c, label=name, lw=1.5)
        axes[1].plot(t, p_mean, ls, color=c, label=name, lw=1.5)
        if arr.shape[1] >= 6:  # has fidelity column (grid/DVR)
            axes[2].plot(t, arr[:, 5], ls, color=c, label=name, lw=1.5)
    if snap_e is not None and snap_e.size:
        if snap_e.ndim == 1:
            snap_e = snap_e.reshape(1, -1)
        axes[2].plot(snap_e[:, 0], snap_e[:, 2], "o", color="C0",
                     label="ECG (snap)")
    axes[0].set_xlabel("t"); axes[0].set_ylabel("<x>")
    axes[0].set_title(f"<x>(t)  (N={N}, K={K})")
    axes[1].set_xlabel("t"); axes[1].set_ylabel("<p>")
    axes[1].set_title(f"<p>(t)  (N={N}, K={K})")
    axes[2].set_xlabel("t"); axes[2].set_ylabel("|<ψ(0)|ψ(t)>|²")
    axes[2].set_title(f"Fidelity  (N={N}, K={K})")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_observables_N{N}.png"), dpi=120)
    plt.close(fig)


def _heatmap(ax, times, grid, M, title, cmap="viridis", vmin=None, vmax=None):
    """times along y-axis, grid along x-axis."""
    extent = [grid.min(), grid.max(), times.min(), times.max()]
    im = ax.imshow(M, aspect="auto", origin="lower", extent=extent,
                   cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    return im


def plot_step4_density_heatmap_single(N, K):
    p = _step4_paths(N, K)
    if not _have(p["ecg_dens"]) or not _have(p["grid_dens"]):
        return
    te, xe, Me = _load_long_grid(p["ecg_dens"],  "n_x")
    tg, xg, Mg = _load_long_grid(p["grid_dens"], "n_x")

    # Resample grid M onto ECG x-grid for difference panel
    Mg_re = np.zeros_like(Me)
    for i in range(len(tg)):
        Mg_re[i, :] = np.interp(xe, xg, Mg[i, :])
    diff = Me - Mg_re

    vmax_d = max(Me.max(), Mg.max())
    vmax_diff = np.max(np.abs(diff))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    _heatmap(axes[0], te, xe, Me, f"ECG  n(x,t)  (N={N})", vmin=0, vmax=vmax_d)
    _heatmap(axes[1], tg, xg, Mg, "FD grid  n(x,t)",       vmin=0, vmax=vmax_d)
    im = _heatmap(axes[2], te, xe, diff, "ECG - grid", cmap="RdBu_r",
                  vmin=-vmax_diff, vmax=vmax_diff)
    for ax in axes:
        ax.set_xlabel("x"); ax.set_ylabel("t")
        ax.set_xlim(-6, 6)
    fig.colorbar(im, ax=axes[2], fraction=0.04)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_density_heatmap_N{N}.png"), dpi=120)
    plt.close(fig)


def plot_step4_momentum_heatmap_single(N, K):
    p = _step4_paths(N, K)
    if not _have(p["ecg_nk"]) or not _have(p["grid_nk"]):
        return
    te, ke, Me = _load_long_grid(p["ecg_nk"],  "n_k")
    tg, kg, Mg = _load_long_grid(p["grid_nk"], "n_k")
    Mg_re = np.zeros_like(Me)
    for i in range(len(tg)):
        Mg_re[i, :] = np.interp(ke, kg, Mg[i, :])
    diff = Me - Mg_re

    vmax_d = max(Me.max(), Mg.max())
    vmax_diff = np.max(np.abs(diff))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    _heatmap(axes[0], te, ke, Me, f"ECG  n(k,t)  (N={N})", vmin=0, vmax=vmax_d)
    _heatmap(axes[1], tg, kg, Mg, "FD grid  n(k,t)",       vmin=0, vmax=vmax_d)
    im = _heatmap(axes[2], te, ke, diff, "ECG - grid", cmap="RdBu_r",
                  vmin=-vmax_diff, vmax=vmax_diff)
    for ax in axes:
        ax.set_xlabel("k"); ax.set_ylabel("t")
        ax.set_xlim(-6, 6)
    fig.colorbar(im, ax=axes[2], fraction=0.04)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_momentum_heatmap_N{N}.png"), dpi=120)
    plt.close(fig)


def plot_step4_snapshots_overlay_single(N, K):
    p = _step4_paths(N, K)
    if not (_have(p["ecg_dens"]) and _have(p["grid_dens"]) and _have(p["dvr_dens"])):
        return
    te, xe, Me_dens = _load_long_grid(p["ecg_dens"],  "n_x")
    tg, xg, Mg_dens = _load_long_grid(p["grid_dens"], "n_x")
    td, xd, Md_dens = _load_long_grid(p["dvr_dens"],  "n_x")
    te2, ke, Me_nk  = _load_long_grid(p["ecg_nk"],    "n_k")
    tg2, kg, Mg_nk  = _load_long_grid(p["grid_nk"],   "n_k")
    td2, kd, Md_nk  = _load_long_grid(p["dvr_nk"],    "n_k")

    # Pick 5 snapshot times evenly spaced
    n_t = len(te)
    if n_t < 5:
        chosen = list(range(n_t))
    else:
        chosen = [int(round(i * (n_t - 1) / 4)) for i in range(5)]
    fig, axes = plt.subplots(2, len(chosen), figsize=(3 * len(chosen), 6),
                             squeeze=False)
    for col, idx in enumerate(chosen):
        t = te[idx]
        # n(x)
        ax = axes[0, col]
        ax.plot(xe, Me_dens[idx],            "-",  color="C0", label="ECG", lw=1.5)
        ax.plot(xg, Mg_dens[idx],            "--", color="C1", label="grid", lw=1.0)
        ax.plot(xd, Md_dens[idx],            ":",  color="C2", label="DVR",  lw=1.0)
        ax.set_xlim(-6, 6)
        ax.set_title(f"t = {t:.3f}")
        ax.set_xlabel("x"); ax.set_ylabel("n(x)")
        if col == 0:
            ax.legend(fontsize=8)
        # n(k)
        ax = axes[1, col]
        ax.plot(ke, Me_nk[idx],            "-",  color="C0", label="ECG", lw=1.5)
        ax.plot(kg, Mg_nk[idx],            "--", color="C1", label="grid", lw=1.0)
        ax.plot(kd, Md_nk[idx],            ":",  color="C2", label="DVR",  lw=1.0)
        ax.set_xlim(-6, 6)
        ax.set_xlabel("k"); ax.set_ylabel("n(k)")
        if col == 0:
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"fig_step4_snapshots_overlay_N{N}.png"), dpi=120)
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
        plot_step4_energy_drift_single(N, K)
        plot_step4_observables_single(N, K)
        plot_step4_density_heatmap_single(N, K)
        plot_step4_momentum_heatmap_single(N, K)
        plot_step4_snapshots_overlay_single(N, K)
    print(f"wrote figures to {OUT}/fig_*_N*.png")


if __name__ == "__main__":
    main()
