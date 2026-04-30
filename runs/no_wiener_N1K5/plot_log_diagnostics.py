"""Parse the verbose [rt_tdvp] step lines from the --no-wiener log and plot
per-step quantities not available in the trace CSV: cond(C), |dz|, min Re(A+B),
min Re(A). Also plots snapshot fidelity from step4_ecg_snap_N1_K5.csv,
and observables (E, norm, <x>, <p>, <x^2>, <p^2>) parsed straight from the log.

Writes:
  out/verify/fig_step4_log_diagnostics_N1_K5.png
  out/verify/fig_step4_fidelity_N1_K5.png
  out/verify/fig_step4_log_observables_N1_K5.png
"""
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOG  = "/home/gyqyan/zexuan/ecg1d_cpp/logs/realtime_N1_K5_no_wiener.log"
OUT  = "/home/gyqyan/zexuan/ecg1d_cpp/runs/no_wiener_N1K5/out/verify"
SNAP = os.path.join(OUT, "step4_ecg_snap_N1_K5.csv")

LINE_RE = re.compile(
    r"\[rt_tdvp\]\s+step\s+\d+/\d+\s+t=([\d\.eE+-]+)\s+dt=[\d\.eE+-]+\s+"
    r"E=([\d\.eE+-]+)\s+norm=([\d\.eE+-]+)\s+"
    r"<x>=([\d\.eE+-]+)\s+<p>=([\d\.eE+-]+)\s+"
    r"<x2>=([\d\.eE+-]+)\s+<p2>=([\d\.eE+-]+)\s+"
    r"cond=([\d\.eE+-]+).*?\|dz\|=([\d\.eE+-]+)\s+"
    r"min\(Re A\+B\)=([\d\.eE+-]+)@\d+\s+min\(Re A\)=([\d\.eE+-]+)@\d+"
)
LEG_RE = re.compile(r"\[step4\]\s+leg\s+t=\[([\d\.eE+-]+),\s*([\d\.eE+-]+)\]")

def parse_log():
    rows = []
    leg_start = 0.0
    with open(LOG) as f:
        for line in f:
            m_leg = LEG_RE.search(line)
            if m_leg:
                leg_start = float(m_leg.group(1))
                continue
            m = LINE_RE.search(line)
            if m:
                vals = [float(g) for g in m.groups()]
                vals[0] += leg_start  # promote in-leg t to global t
                rows.append(vals)
    a = np.array(rows)
    return {
        "t":    a[:, 0],
        "E":    a[:, 1],
        "norm": a[:, 2],
        "x":    a[:, 3],
        "p":    a[:, 4],
        "x2":   a[:, 5],
        "p2":   a[:, 6],
        "cond": a[:, 7],
        "dz":   a[:, 8],
        "minABreal": a[:, 9],
        "minAreal":  a[:,10],
    }

def main():
    d = parse_log()
    t = d["t"]
    print(f"parsed {len(t)} step lines  t in [{t[0]:.3f}, {t[-1]:.3f}]")

    # ------ 4-panel solver-health diagnostic ------
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

    ax = axes[0, 0]
    ax.semilogy(t, d["cond"], "-", color="C3", lw=1.0)
    ax.set_ylabel(r"cond$(C)$")
    ax.set_title("C-matrix conditioning")
    ax.grid(alpha=0.3, which="both")

    ax = axes[0, 1]
    ax.semilogy(t, d["dz"], "-", color="C0", lw=1.0)
    ax.set_ylabel(r"$|dz|$ per step")
    ax.set_title(r"per-step parameter increment")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1, 0]
    ax.plot(t, d["minABreal"], "-", color="C2", lw=1.0, label=r"$\min\,\mathrm{Re}(A_k+B_k)$")
    ax.plot(t, d["minAreal"],  "-", color="C4", lw=1.0, label=r"$\min\,\mathrm{Re}(A_k)$")
    ax.axhline(0.0, color="k", lw=0.4, alpha=0.4)
    ax.set_ylabel(r"min real part")
    ax.set_xlabel(r"$t$")
    ax.set_title("Gaussian width health (positivity)")
    ax.legend(fontsize=9, loc="best")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.plot(t, d["norm"], "-", color="C1", lw=1.0)
    ax.axhline(d["norm"][0], color="k", lw=0.4, alpha=0.4,
               label=fr"$\langle\psi|\psi\rangle_0={d['norm'][0]:.3f}$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle\psi|\psi\rangle(t)$")
    ax.set_title("Raw norm trajectory  (--no-wiener, no enforce)")
    ax.legend(fontsize=9, loc="best")
    ax.grid(alpha=0.3)

    fig.suptitle("Step 4 — solver diagnostics (parsed from verbose log)  N=1, K=5, --no-wiener",
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(OUT, "fig_step4_log_diagnostics_N1_K5.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")

    # ------ snapshot fidelity ------
    snap = np.genfromtxt(SNAP, delimiter=",", names=True, dtype=float)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(snap["t"], snap["fidelity"], "o-", color="C0", lw=1.6, ms=6)
    ax.axhline(1.0, color="k", lw=0.4, alpha=0.4)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|\langle\psi_{\rm ECG}(t)|\psi_{\rm grid}(t)\rangle|^2$")
    ax.set_title("Step 4 — snapshot fidelity vs FD-grid reference  (N=1, K=5, --no-wiener)")
    ax.set_ylim(0.85, 1.02)
    ax.grid(alpha=0.3)
    for ti, fi in zip(snap["t"], snap["fidelity"]):
        ax.annotate(f"{fi:.4f}", (ti, fi), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=8)
    fig.tight_layout()
    out = os.path.join(OUT, "fig_step4_fidelity_N1_K5.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")

    # ------ observables traced straight from the log ------
    # Convention reminder: the log's E= field IS the rescaled energy
    # E_resc = <psi|H|psi> / <psi|psi> (see realtime_tdvp.cpp:489).
    # Raw matrix element: E_raw = E_resc * norm.
    E_resc = d["E"]
    E_raw  = d["E"] * d["norm"]

    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)

    ax = axes[0, 0]
    ax.plot(t, E_resc, "-", color="C0", lw=1.0)
    ax.axhline(E_resc[0], color="k", lw=0.4, alpha=0.4,
               label=fr"$E_0={E_resc[0]:.6f}$")
    ax.set_ylabel(r"$\langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$")
    ax.set_title("Energy (norm-rescaled)  — log's E= field")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.plot(t, E_raw, "-", color="C3", lw=1.0)
    ax.axhline(E_raw[0], color="k", lw=0.4, alpha=0.4,
               label=fr"$\langle\psi|H|\psi\rangle_0={E_raw[0]:.4f}$")
    ax.set_ylabel(r"$\langle\psi|H|\psi\rangle$")
    ax.set_title("Energy (raw matrix element = E·norm)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(t, d["x"], "-", color="C2", lw=1.0)
    ax.axhline(0.0, color="k", lw=0.4, alpha=0.4)
    ax.set_ylabel(r"$\langle x\rangle$")
    ax.set_title(r"position centroid")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.plot(t, d["p"], "-", color="C4", lw=1.0)
    ax.axhline(0.0, color="k", lw=0.4, alpha=0.4)
    ax.set_ylabel(r"$\langle p\rangle$")
    ax.set_title(r"momentum centroid")
    ax.grid(alpha=0.3)

    ax = axes[2, 0]
    ax.plot(t, d["x2"], "-", color="C5", lw=1.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle x^2\rangle$")
    ax.set_title(r"position variance proxy")
    ax.grid(alpha=0.3)

    ax = axes[2, 1]
    ax.plot(t, d["p2"], "-", color="C6", lw=1.0)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle p^2\rangle$")
    ax.set_title(r"momentum variance proxy")
    ax.grid(alpha=0.3)

    fig.suptitle("Step 4 — observables parsed from log  N=1, K=5, --no-wiener",
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(OUT, "fig_step4_log_observables_N1_K5.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"wrote {out}")

if __name__ == "__main__":
    main()
