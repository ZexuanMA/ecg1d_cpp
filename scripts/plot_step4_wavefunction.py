"""Plot wavefunction-density dynamics from Step 4 traces.

Reads (long format: t,x,n_x   or   t,k,n_k):
  out/verify/step4_ecg_density_N{N}_K{K}.csv
  out/verify/step4_grid_density_N{N}.csv
  out/verify/step4_dvr_density_N{N}.csv
  out/verify/step4_ecg_nk_N{N}_K{K}.csv
  out/verify/step4_grid_nk_N{N}.csv
  out/verify/step4_dvr_nk_N{N}.csv

Writes:
  fig_step4_density_overlay_N{N}_K{K}.png   density |psi(x,t)|^2 vs x at multiple t
  fig_step4_nk_overlay_N{N}_K{K}.png        n(k,t) at multiple t
  fig_step4_density_heat_N{N}_K{K}.png      ECG density heatmap (x,t)
  fig_step4_nk_heat_N{N}_K{K}.png           ECG n(k) heatmap (k,t)
  fig_step4_density_diff_N{N}_K{K}.png      pointwise |n_ECG - n_grid|(x,t) heatmap

Usage:
  python scripts/plot_step4_wavefunction.py [N=1] [K=3]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

OUT = "out/verify"


def load_long(path):
    """Read t,grid,value long-format CSV. Returns (times, grid, M[t_idx, g_idx])."""
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    ts = np.unique(arr[:, 0])
    gs = np.unique(arr[:, 1])
    M = np.zeros((len(ts), len(gs)))
    t_idx = {t: i for i, t in enumerate(ts)}
    g_idx = {g: i for i, g in enumerate(gs)}
    for ti, gi, vi in arr:
        M[t_idx[ti], g_idx[gi]] = vi
    return ts, gs, M


def small_multiples(times_e, xs_e, M_e, times_g, xs_g, M_g, times_d, xs_d, M_d,
                     xlabel, ylabel, suptitle, outpath, xlim):
    """One panel per snapshot time. ECG (solid) vs FD grid (dashed) vs DVR (dotted).
    Grid/DVR are interpolated onto the ECG x-grid for direct overlay; the
    nearest snapshot in the reference's own time list is used.
    If grid/DVR data is None (e.g. N>=2), only ECG is shown.
    """
    n = len(times_e)
    n_cols = 4
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.0 * n_cols, 2.8 * n_rows),
                              sharex=True, sharey=True)
    axes_flat = axes.flatten() if n_rows > 1 else axes

    # Find global y-max for consistent scale.
    candidates = [M_e.max()]
    if M_g is not None: candidates.append(M_g.max())
    if M_d is not None: candidates.append(M_d.max())
    ymax = max(candidates) * 1.05

    for i, t in enumerate(times_e):
        ax = axes_flat[i]
        ax.plot(xs_e, M_e[i], "-", color="C0", lw=1.8, label="ECG")
        if M_g is not None:
            jg = int(np.argmin(np.abs(times_g - t)))
            ax.plot(xs_e, np.interp(xs_e, xs_g, M_g[jg]), "--", color="C1", lw=1.4,
                    label="FD grid")
        if M_d is not None:
            jd = int(np.argmin(np.abs(times_d - t)))
            ax.plot(xs_e, np.interp(xs_e, xs_d, M_d[jd]), ":",  color="C2", lw=1.4,
                    label="sinc-DVR")
        ax.set_title(f"t = {t:.3f}", fontsize=10)
        ax.grid(alpha=0.3)
        ax.set_xlim(xlim)
        ax.set_ylim(0, ymax)

    # Hide unused panels.
    for j in range(n, len(axes_flat)):
        axes_flat[j].axis("off")

    # Single legend on the first panel.
    axes_flat[0].legend(fontsize=9, loc="upper right")

    # Common axis labels.
    for ax in axes[:, 0] if n_rows > 1 else [axes_flat[0]]:
        ax.set_ylabel(ylabel)
    bottom_axes = axes[-1, :] if n_rows > 1 else axes_flat
    for ax in bottom_axes:
        ax.set_xlabel(xlabel)

    fig.suptitle(suptitle, fontsize=12)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print(f"wrote {outpath}")


def comparison_residuals(times_e, xs_e, M_e, times_g, xs_g, M_g,
                          xlabel, ylabel, suptitle, outpath, xlim):
    """Stacked residual plot: M_ECG - M_grid at every snapshot, one line per t,
    with a colormap encoding time so the eye can pick out where the ansatz
    starts to deviate."""
    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = cm.get_cmap("plasma")
    n = len(times_e)
    for i, t in enumerate(times_e):
        jg = int(np.argmin(np.abs(times_g - t)))
        diff = M_e[i] - np.interp(xs_e, xs_g, M_g[jg])
        ax.plot(xs_e, diff, color=cmap(i / max(n - 1, 1)),
                lw=1.4, label=f"t={t:.2f}")
    ax.axhline(0.0, color="k", lw=0.5, alpha=0.4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    ax.set_title(suptitle)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print(f"wrote {outpath}")


def main(N=1, K=3):
    paths = {
        "ecg_dx":  os.path.join(OUT, f"step4_ecg_density_N{N}_K{K}.csv"),
        "grid_dx": os.path.join(OUT, f"step4_grid_density_N{N}.csv"),
        "dvr_dx":  os.path.join(OUT, f"step4_dvr_density_N{N}.csv"),
        "ecg_nk":  os.path.join(OUT, f"step4_ecg_nk_N{N}_K{K}.csv"),
        "grid_nk": os.path.join(OUT, f"step4_grid_nk_N{N}.csv"),
        "dvr_nk":  os.path.join(OUT, f"step4_dvr_nk_N{N}.csv"),
    }
    if not (os.path.exists(paths["ecg_dx"]) and os.path.exists(paths["ecg_nk"])):
        print(f"missing ECG snapshot CSVs for N={N} K={K}; run step 4 first",
              file=sys.stderr)
        sys.exit(1)
    have_ref = all(os.path.exists(paths[k]) for k in
                   ("grid_dx", "dvr_dx", "grid_nk", "dvr_nk"))
    if not have_ref:
        print(f"[note] no reference snapshots for N={N} (grid/DVR are N=1 only); "
              "skipping residual plots, small-multiples will show ECG only.")

    # --- small multiples: |psi(x,t)|^2 at each snapshot ---
    te, xe, De = load_long(paths["ecg_dx"])
    tg, xg, Dg = load_long(paths["grid_dx"]) if have_ref else (None, None, None)
    td, xd, Dd = load_long(paths["dvr_dx"])  if have_ref else (None, None, None)
    title_compare = "ECG vs FD grid vs sinc-DVR" if have_ref else "ECG only"
    small_multiples(
        te, xe, De, tg, xg, Dg, td, xd, Dd,
        xlabel=r"$x$",
        ylabel=r"$|\psi(x,t)|^2$",
        suptitle=rf"Step 4 — density snapshots: {title_compare}  (N={N}, K={K})",
        outpath=os.path.join(OUT, f"fig_step4_density_overlay_N{N}_K{K}.png"),
        xlim=(-6, 6),
    )
    if have_ref:
        comparison_residuals(
            te, xe, De, tg, xg, Dg,
            xlabel=r"$x$",
            ylabel=r"$|\psi_{\rm ECG}|^2 - |\psi_{\rm grid}|^2$",
            suptitle=rf"Step 4 — ECG−grid density residual at each snapshot  (N={N}, K={K})",
            outpath=os.path.join(OUT, f"fig_step4_density_residual_N{N}_K{K}.png"),
            xlim=(-6, 6),
        )

    # --- small multiples: n(k,t) ---
    te, ke, Ne = load_long(paths["ecg_nk"])
    tg, kg, Ng = load_long(paths["grid_nk"]) if have_ref else (None, None, None)
    td, kd, Nd = load_long(paths["dvr_nk"])  if have_ref else (None, None, None)
    small_multiples(
        te, ke, Ne, tg, kg, Ng, td, kd, Nd,
        xlabel=r"$k$",
        ylabel=r"$n(k,t)$",
        suptitle=rf"Step 4 — momentum density snapshots: {title_compare}  (N={N}, K={K})",
        outpath=os.path.join(OUT, f"fig_step4_nk_overlay_N{N}_K{K}.png"),
        xlim=(-5, 5),
    )
    if have_ref:
        comparison_residuals(
            te, ke, Ne, tg, kg, Ng,
            xlabel=r"$k$",
            ylabel=r"$n_{\rm ECG}(k) - n_{\rm grid}(k)$",
            suptitle=rf"Step 4 — ECG−grid momentum-density residual at each snapshot  (N={N}, K={K})",
            outpath=os.path.join(OUT, f"fig_step4_nk_residual_N{N}_K{K}.png"),
            xlim=(-5, 5),
        )

    # --- ECG-only heatmaps in (x,t) and (k,t) ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    im0 = axes[0].imshow(De.T, aspect="auto", origin="lower",
                         extent=[te[0], te[-1], xe[0], xe[-1]],
                         cmap="viridis")
    axes[0].set_xlabel(r"$t$"); axes[0].set_ylabel(r"$x$")
    axes[0].set_ylim(-6, 6)
    axes[0].set_title(rf"ECG  $|\psi(x,t)|^2$")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(Ne.T, aspect="auto", origin="lower",
                         extent=[te[0], te[-1], ke[0], ke[-1]],
                         cmap="viridis")
    axes[1].set_xlabel(r"$t$"); axes[1].set_ylabel(r"$k$")
    axes[1].set_ylim(-5, 5)
    axes[1].set_title(rf"ECG  $n(k,t)$")
    fig.colorbar(im1, ax=axes[1])
    fig.suptitle(rf"Step 4 — ECG wavefunction heatmaps  (N={N}, K={K})", fontsize=12)
    fig.tight_layout()
    out = os.path.join(OUT, f"fig_step4_density_heat_N{N}_K{K}.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")

    # --- difference heatmap ECG vs grid (resampled to ECG x-grid) ---
    if not have_ref:
        return  # ECG-only mode (N>=2): no reference grid to diff against.
    # Reload position-density data here (the n(k) load above clobbered te/De)
    te, xe, De = load_long(paths["ecg_dx"])
    tg, xg, Dg = load_long(paths["grid_dx"])
    # Grid is on a much denser grid; resample to the ECG x grid by interpolation.
    Dg_on_xe = np.zeros_like(De)
    for i, t_ecg in enumerate(te):
        # find nearest grid time
        j = int(np.argmin(np.abs(tg - t_ecg)))
        Dg_on_xe[i] = np.interp(xe, xg, Dg[j])
    diff = De - Dg_on_xe

    fig, ax = plt.subplots(figsize=(7, 4.5))
    vmax = np.max(np.abs(diff))
    im = ax.imshow(diff.T, aspect="auto", origin="lower",
                   extent=[te[0], te[-1], xe[0], xe[-1]],
                   cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xlabel(r"$t$"); ax.set_ylabel(r"$x$")
    ax.set_ylim(-6, 6)
    ax.set_title(rf"$|\psi_{{\rm ECG}}|^2 - |\psi_{{\rm grid}}|^2$  (N={N}, K={K})")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    out = os.path.join(OUT, f"fig_step4_density_diff_N{N}_K{K}.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    main(N, K)
