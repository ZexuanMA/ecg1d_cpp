#!/usr/bin/env python3
"""
Grid-based fluorescent spot intensity analysis v3.

Key improvements over v2:
- Per-image rotation detection (chip may be repositioned between measurements)
- Contrast enhancement for reliable angle detection on quenched images
- Periodic grid fitting: finds ALL well positions using known grid period,
  not just bright peaks. Works even when most spots are quenched.
- Center-of-mass refinement for quantum dot drift within wells.
"""

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage, signal
import json, os

DIR = os.path.dirname(os.path.abspath(__file__))

# Grid parameters (measured from 0ugL autocorrelation)
ROW_PERIOD = 42   # pixels between rows
COL_PERIOD = 49   # pixels between columns within a row


def load_red_channel(name):
    path = os.path.join(DIR, f"{name}.bmp")
    img = np.array(Image.open(path))
    return img[:, :, 0].astype(np.float64)


# ──────────────────────────────────────────────
# Rotation detection with contrast enhancement
# ──────────────────────────────────────────────
def detect_rotation_angle(R):
    """Two-stage rotation detection with contrast enhancement.
    Works for both bright (0ugL) and quenched (20ugL) images."""
    scale = 2
    R_s = R[::scale, ::scale]
    H, W = R_s.shape
    margin = int(min(H, W) * 0.12)
    R_crop = R_s[margin:H - margin, margin:W - margin]

    # Contrast enhancement (critical for dim images)
    p2, p98 = np.percentile(R_crop, [2, 98])
    R_e = np.clip((R_crop - p2) / max(p98 - p2, 1), 0, 1)

    def score(a):
        Rr = ndimage.rotate(R_e, a, reshape=False, order=1,
                            mode='reflect') if abs(a) > 0.01 else R_e
        v = ndimage.gaussian_filter1d(Rr.mean(axis=1), sigma=1.5)
        return np.sum(np.diff(v) ** 2)

    # Coarse search
    coarse = np.arange(-8, 8.01, 0.2)
    cs = np.array([score(a) for a in coarse])
    best_c = coarse[np.argmax(cs)]

    # Fine search
    fine = np.arange(best_c - 0.5, best_c + 0.501, 0.02)
    fs = np.array([score(a) for a in fine])
    best_f = fine[np.argmax(fs)]

    quality = cs.max() / max(np.median(cs), 1e-10)
    return best_f, quality


# ──────────────────────────────────────────────
# Row detection
# ──────────────────────────────────────────────
def detect_rows(R, expected_spacing=ROW_PERIOD):
    """Find row y-positions from vertical projection."""
    v_proj = R.mean(axis=1)
    v_smooth = ndimage.gaussian_filter1d(v_proj, sigma=3)
    row_peaks, _ = signal.find_peaks(
        v_smooth,
        distance=int(expected_spacing * 0.7),
        prominence=0.5,  # low prominence to catch faint bands
    )

    if len(row_peaks) < 10:
        # Fallback: use periodic grid with phase search
        row_peaks = _periodic_fit_1d(v_smooth, expected_spacing)

    return row_peaks


def _periodic_fit_1d(profile, period):
    """Find best periodic grid that matches the profile."""
    best_score = -np.inf
    best_offset = 0
    for offset in np.arange(0, period, 0.5):
        positions = np.arange(offset, len(profile), period).astype(int)
        positions = positions[(positions >= 0) & (positions < len(profile))]
        if len(positions) < 3:
            continue
        score = profile[positions].sum()
        if score > best_score:
            best_score = score
            best_offset = offset
    positions = np.arange(best_offset, len(profile), period).astype(int)
    positions = positions[(positions >= 5) & (positions < len(profile) - 5)]
    return positions


# ──────────────────────────────────────────────
# Column detection within each row (periodic fitting)
# ──────────────────────────────────────────────
def detect_columns_in_row(R, row_y, col_period=COL_PERIOD, half_width=5):
    """Find ALL column (well) positions using periodic grid fitting.
    Always uses periodicity to ensure quenched wells are included."""
    y0 = max(0, row_y - half_width)
    y1 = min(R.shape[0], row_y + half_width + 1)
    strip = R[y0:y1, :].mean(axis=0)
    strip_smooth = ndimage.gaussian_filter1d(strip, sigma=2)

    # Estimate local period from autocorrelation
    s = strip_smooth - strip_smooth.mean()
    auto = np.correlate(s, s, mode='full')
    auto = auto[len(auto) // 2:]
    auto_peaks, _ = signal.find_peaks(auto, distance=int(col_period * 0.5))

    # Pick the autocorrelation peak closest to expected period
    local_period = col_period
    if len(auto_peaks) > 0:
        diffs = np.abs(auto_peaks - col_period)
        best_idx = np.argmin(diffs)
        if diffs[best_idx] < col_period * 0.3:
            local_period = auto_peaks[best_idx]

    # Periodic fitting: find best phase offset
    best_score = -np.inf
    best_offset = 0
    for offset in np.arange(0, local_period, 0.5):
        positions = np.arange(offset, len(strip_smooth), local_period).astype(int)
        positions = positions[(positions >= 3) & (positions < len(strip_smooth) - 3)]
        if len(positions) < 3:
            continue
        score = strip_smooth[positions].sum()
        if score > best_score:
            best_score = score
            best_offset = offset

    columns = np.arange(best_offset, len(strip_smooth), local_period).astype(int)
    columns = columns[(columns >= 3) & (columns < len(strip_smooth) - 3)]
    return columns


# ──────────────────────────────────────────────
# Center-of-mass refinement
# ──────────────────────────────────────────────
def refine_center(R, y, x, radius=4):
    """Sub-pixel position refinement using center-of-mass."""
    y0 = max(0, y - radius)
    y1 = min(R.shape[0], y + radius + 1)
    x0 = max(0, x - radius)
    x1 = min(R.shape[1], x + radius + 1)
    patch = R[y0:y1, x0:x1]
    patch_bg = patch - patch.min()
    total = patch_bg.sum()
    if total < 1e-6:
        return float(y), float(x)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    return (yy * patch_bg).sum() / total, (xx * patch_bg).sum() / total


# ──────────────────────────────────────────────
# Full grid detection
# ──────────────────────────────────────────────
def detect_grid(R):
    """Detect full grid on one image. Returns list of (row_idx, cy, cx)."""
    row_ys = detect_rows(R)
    print(f"    {len(row_ys)} rows, spacing ~{np.median(np.diff(row_ys)):.0f} px")

    all_spots = []
    counts = []
    for ri, ry in enumerate(row_ys):
        cols = detect_columns_in_row(R, ry)
        for cx in cols:
            cy, cx_ref = refine_center(R, ry, cx)
            all_spots.append((ri, cy, cx_ref))
        counts.append(len(cols))

    counts = np.array(counts)
    print(f"    Spots/row: mean={counts.mean():.1f}, min={counts.min()}, max={counts.max()}")
    print(f"    Total: {len(all_spots)}")
    return all_spots, row_ys, counts


# ──────────────────────────────────────────────
# Intensity measurement
# ──────────────────────────────────────────────
def measure_intensity(R, cy, cx, spot_radius=3, bg_inner=6, bg_outer=10):
    """Net signal = spot mean (7x7) - background mean (annulus r=6..10)."""
    H, W = R.shape
    iy, ix = int(round(cy)), int(round(cx))

    sy0 = max(0, iy - spot_radius)
    sy1 = min(H, iy + spot_radius + 1)
    sx0 = max(0, ix - spot_radius)
    sx1 = min(W, ix + spot_radius + 1)
    spot_val = R[sy0:sy1, sx0:sx1].mean()

    by0 = max(0, iy - bg_outer)
    by1 = min(H, iy + bg_outer + 1)
    bx0 = max(0, ix - bg_outer)
    bx1 = min(W, ix + bg_outer + 1)
    local_yy = np.arange(by1 - by0)[:, None]
    local_xx = np.arange(bx1 - bx0)[None, :]
    local_dist_sq = (local_yy - (cy - by0))**2 + (local_xx - (cx - bx0))**2
    local_mask = (local_dist_sq >= bg_inner**2) & (local_dist_sq <= bg_outer**2)
    if local_mask.sum() > 5:
        bg_val = R[by0:by1, bx0:bx1][local_mask].mean()
    else:
        bg_val = np.median(R[by0:by1, bx0:bx1])

    return spot_val, bg_val, spot_val - bg_val


def measure_all_spots(R, spots):
    """Measure intensity for all spots in one image."""
    data = [measure_intensity(R, cy, cx) for _, cy, cx in spots]
    arr = np.array(data)
    nets = arr[:, 2]
    print(f"    net signal: mean={nets.mean():.1f}, std={nets.std():.1f}, "
          f"median={np.median(nets):.1f}")
    return arr


# ──────────────────────────────────────────────
# Block partitioning and quality filtering
# ──────────────────────────────────────────────
def partition_and_filter(spots, measurements, n_br=7, n_bc=7, min_spots=5,
                         bg_outlier_sigma=2.5):
    """Partition into blocks and filter out artifact regions."""
    arr = np.array([(cy, cx) for _, cy, cx in spots])
    y_edges = np.linspace(arr[:, 0].min() - 1, arr[:, 0].max() + 1, n_br + 1)
    x_edges = np.linspace(arr[:, 1].min() - 1, arr[:, 1].max() + 1, n_bc + 1)

    block_ids = np.zeros(len(spots), dtype=int)
    for i, (_, cy, cx) in enumerate(spots):
        by = np.clip(np.searchsorted(y_edges, cy) - 1, 0, n_br - 1)
        bx = np.clip(np.searchsorted(x_edges, cx) - 1, 0, n_bc - 1)
        block_ids[i] = by * n_bc + bx

    n_blocks = n_br * n_bc
    keep = np.ones(len(spots), dtype=bool)
    bg_means = []
    for b in range(n_blocks):
        mask = block_ids == b
        if mask.sum() < min_spots:
            keep[mask] = False
            bg_means.append(np.nan)
            continue
        bg_means.append(measurements[mask, 1].mean())

    bg_means = np.array(bg_means)
    valid = bg_means[~np.isnan(bg_means)]
    med = np.median(valid)
    mad = np.median(np.abs(valid - med))
    thresh = med + bg_outlier_sigma * mad * 1.4826

    n_excl = 0
    for b in range(n_blocks):
        if np.isnan(bg_means[b]):
            continue
        if bg_means[b] > thresh:
            keep[block_ids == b] = False
            n_excl += 1

    print(f"    {n_excl} blocks excluded, {keep.sum()}/{len(keep)} spots kept")
    return keep, block_ids, y_edges, x_edges


# ──────────────────────────────────────────────
# Intensity discretization
# ──────────────────────────────────────────────
def discretize_intensity(ref_net):
    """5-level discretization based on noise estimate and signal range."""
    p5 = np.percentile(ref_net, 5)
    p10 = np.percentile(ref_net, 10)
    quench = max(abs(p5), abs(p10), 5.0)

    bright = ref_net[ref_net > quench]
    if len(bright) < 20:
        return np.array([-np.inf, quench, np.inf]), 2

    p90 = np.percentile(bright, 90)
    r = p90 - quench
    edges = np.array([-np.inf, quench, quench + r * 0.22, quench + r * 0.50,
                      quench + r * 0.80, np.inf])
    n_levels = len(edges) - 1
    print(f"    {n_levels} levels, quench={quench:.1f}")
    print(f"    edges: {[f'{e:.1f}' for e in edges[1:-1]]}")
    return edges, n_levels


def compute_distribution(net_signals, bin_edges):
    n = len(bin_edges) - 1
    counts = np.zeros(n)
    for i in range(n):
        if i == 0:
            counts[i] = (net_signals <= bin_edges[1]).sum()
        elif i == n - 1:
            counts[i] = (net_signals > bin_edges[-2]).sum()
        else:
            counts[i] = ((net_signals > bin_edges[i]) &
                         (net_signals <= bin_edges[i + 1])).sum()
    total = counts.sum()
    return (counts / total if total > 0 else counts), counts


# ──────────────────────────────────────────────
# Visualization
# ──────────────────────────────────────────────
def save_crop(R, spots, keep, net_signals, bin_edges, filename,
              crop_center, crop_size=300, spot_radius=4):
    """Save a crop region with color-coded spot circles."""
    cy_c, cx_c = crop_center
    H, W = R.shape
    y0 = max(0, cy_c - crop_size)
    y1 = min(H, cy_c + crop_size)
    x0 = max(0, cx_c - crop_size)
    x1 = min(W, cx_c + crop_size)

    crop = R[y0:y1, x0:x1]
    vmin, vmax = crop.min(), crop.max()
    norm = ((crop - vmin) / max(vmax - vmin, 1) * 220).astype(np.uint8)
    rgb = np.stack([norm, norm // 4, norm // 6], axis=2)
    img = Image.fromarray(rgb)
    draw = ImageDraw.Draw(img)

    colors = [
        (80, 80, 255),    # 0: quenched
        (0, 200, 200),    # 1: dim
        (0, 255, 0),      # 2: medium
        (255, 200, 0),    # 3: bright
        (255, 50, 50),    # 4: very bright
    ]

    for i, (_, cy, cx) in enumerate(spots):
        if not (y0 <= cy < y1 and x0 <= cx < x1):
            continue
        iy = int(round(cy - y0))
        ix = int(round(cx - x0))
        if not keep[i]:
            c = (60, 60, 60)
        else:
            lv = np.clip(np.searchsorted(bin_edges, net_signals[i]) - 1,
                         0, len(colors) - 1)
            c = colors[lv]
        draw.ellipse([ix - spot_radius, iy - spot_radius,
                      ix + spot_radius, iy + spot_radius],
                     outline=c, width=1)

    img.save(os.path.join(DIR, filename))
    print(f"    Saved {filename}")


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    print("=" * 60)
    print("Grid-based Fluorescent Spot Intensity Analysis v3")
    print("  per-image rotation + periodic grid fitting")
    print("=" * 60)

    # Load images
    print("\n[1] Loading images...")
    raw_images = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        raw_images[name] = load_red_channel(name)
        print(f"  {name}: {raw_images[name].shape}")

    # Per-image rotation detection and correction
    print("\n[2] Detecting rotation per image...")
    images = {}
    angles = {}
    for name, R in raw_images.items():
        angle, quality = detect_rotation_angle(R)
        angles[name] = angle
        if abs(angle) > 0.01:
            images[name] = ndimage.rotate(R, angle, reshape=False, order=1,
                                          mode='nearest')
        else:
            images[name] = R
        print(f"  {name}: angle={angle:+.2f} deg, quality={quality:.1f}")

    # Per-image grid detection
    print("\n[3] Detecting grid on each image...")
    all_spots = {}
    all_meas = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        print(f"\n  {name}:")
        spots, row_ys, spr = detect_grid(images[name])
        all_spots[name] = spots

        print(f"  Measuring intensities...")
        meas = measure_all_spots(images[name], spots)
        all_meas[name] = meas

    # Block filtering (per image)
    print("\n[4] Block filtering...")
    all_keep = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        print(f"  {name}:")
        keep, bids, ye, xe = partition_and_filter(
            all_spots[name], all_meas[name])
        all_keep[name] = keep

    # Discretization based on 0ugL
    print("\n[5] Discretization (from 0ugL reference)...")
    ref_net = all_meas["0ugL"][all_keep["0ugL"], 2]
    bin_edges, n_levels = discretize_intensity(ref_net)

    # Distributions
    print("\n[6] Computing distributions...")
    distributions = {}
    level_names = ["quench", "dim", "medium", "bright", "v.bright"][:n_levels]
    for name in ["0ugL", "5ugL", "20ugL"]:
        net = all_meas[name][all_keep[name], 2]
        dist, counts = compute_distribution(net, bin_edges)
        distributions[name] = dist
        parts = [f"{level_names[i]}={dist[i]:.4f}" for i in range(n_levels)]
        print(f"  {name} ({int(all_keep[name].sum())} spots): {', '.join(parts)}")

    # Stochastic dominance
    print("\n[7] Stochastic dominance check...")
    all_ok = True
    for k in range(n_levels):
        p0 = distributions["0ugL"][k:].sum()
        p5 = distributions["5ugL"][k:].sum()
        p20 = distributions["20ugL"][k:].sum()
        ok = (p0 >= p5 - 1e-10) and (p0 >= p20 - 1e-10) and (p5 >= p20 - 1e-10)
        if not ok:
            all_ok = False
        print(f"  k>={k}: 0ugL={p0:.4f}  5ugL={p5:.4f}  20ugL={p20:.4f}  "
              f"{'OK' if ok else 'FAIL'}")
    print(f"  Overall: {'PASS' if all_ok else 'FAIL'}")

    # Save crops: 2 per concentration
    print("\n[8] Saving visualizations...")
    for name in ["0ugL", "5ugL", "20ugL"]:
        R = images[name]
        H, W = R.shape
        spots = all_spots[name]
        keep = all_keep[name]
        net = all_meas[name][:, 2]

        crop_centers = [
            (int(H * 0.33), int(W * 0.40)),
            (int(H * 0.67), int(W * 0.65)),
        ]
        for ci, center in enumerate(crop_centers, 1):
            save_crop(R, spots, keep, net, bin_edges,
                      f"crop{ci}_{name}.png", center, crop_size=300)

    # Distribution chart and stochastic dominance plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
        labels = ["Quenched", "Dim", "Medium", "Bright", "V.Bright"][:n_levels]
        bar_colors = ["#4444ff", "#00c8c8", "#00ff00", "#ffc800", "#ff3232"][:n_levels]
        for ax, name in zip(axes, ["0ugL", "5ugL", "20ugL"]):
            d = distributions[name]
            ax.bar(range(n_levels), d, color=bar_colors)
            ax.set_xticks(range(n_levels))
            ax.set_xticklabels(labels, fontsize=9)
            ax.set_title(name)
            for i, v in enumerate(d):
                ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=8)
        axes[0].set_ylabel("Probability")
        fig.suptitle("Discrete Intensity Distribution (v3)", fontsize=13)
        fig.tight_layout()
        fig.savefig(os.path.join(DIR, "distributions_bar.png"), dpi=150)
        plt.close()
        print("  Saved distributions_bar.png")

        fig, ax = plt.subplots(figsize=(7, 5))
        for name, mk, color in [("0ugL", "o-", "#4488ff"),
                                 ("5ugL", "s-", "#ff9900"),
                                 ("20ugL", "D-", "#ff4444")]:
            surv = [distributions[name][k:].sum() for k in range(n_levels)]
            ax.plot(range(n_levels), surv, mk, label=name, color=color,
                    markersize=8, linewidth=2)
        ax.set_xticks(range(n_levels))
        ax.set_xticklabels([f">={k}" for k in range(n_levels)])
        ax.set_xlabel("Level threshold k")
        ax.set_ylabel("P(Level >= k)")
        ax.set_title("Stochastic Dominance Verification")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(DIR, "stochastic_dominance.png"), dpi=150)
        plt.close()
        print("  Saved stochastic_dominance.png")
    except ImportError:
        print("  (matplotlib not available, skipping charts)")

    # Ridge-Projection summary
    print("\n" + "=" * 60)
    print("Ridge-Projection Input")
    print("=" * 60)
    for x_name, y_name in [("0ugL", "5ugL"), ("0ugL", "20ugL")]:
        print(f"\n  {x_name} -> {y_name}:")
        print(f"  X = {np.array2string(distributions[x_name], precision=6, separator=', ')}")
        print(f"  Y = {np.array2string(distributions[y_name], precision=6, separator=', ')}")

    # Save JSON
    output = {
        "method": "v3: per-image rotation + periodic grid fitting",
        "rotation_angles_deg": {n: round(a, 3) for n, a in angles.items()},
        "grid_periods_px": {"row": ROW_PERIOD, "col": COL_PERIOD},
        "per_image_stats": {},
        "bin_edges": [round(e, 1) if abs(e) < 1e10 else str(e) for e in bin_edges],
        "n_levels": n_levels,
        "distributions": {n: d.tolist() for n, d in distributions.items()},
    }
    for name in ["0ugL", "5ugL", "20ugL"]:
        output["per_image_stats"][name] = {
            "total_spots": len(all_spots[name]),
            "spots_after_filter": int(all_keep[name].sum()),
            "net_signal_mean": round(all_meas[name][all_keep[name], 2].mean(), 2),
        }

    with open(os.path.join(DIR, "grid_analysis_results.json"), "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to grid_analysis_results.json")


if __name__ == "__main__":
    main()
