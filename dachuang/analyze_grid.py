#!/usr/bin/env python3
"""
Grid-based fluorescent spot intensity analysis v3.

Each image is processed independently:
1. Detect and correct rotation (contrast-enhanced, works even on dim images)
2. Detect rows from vertical projection
3. Detect luminous spots via peak detection within each row
4. Measure intensity and classify

For distribution calculation, 0ugL spot count is used as total well count
to estimate quenched fraction in other images.
"""

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage, signal
import json, os

DIR = os.path.dirname(os.path.abspath(__file__))


def load_red_channel(name):
    path = os.path.join(DIR, f"{name}.bmp")
    img = np.array(Image.open(path))
    return img[:, :, 0].astype(np.float64)


# ──────────────────────────────────────────────
# Rotation detection
# ──────────────────────────────────────────────
def detect_rotation_angle(R):
    """Contrast-enhanced two-stage rotation detection."""
    scale = 2
    R_s = R[::scale, ::scale]
    H, W = R_s.shape
    margin = int(min(H, W) * 0.12)
    R_crop = R_s[margin:H - margin, margin:W - margin]

    p2, p98 = np.percentile(R_crop, [2, 98])
    R_e = np.clip((R_crop - p2) / max(p98 - p2, 1), 0, 1)

    def score(a):
        Rr = ndimage.rotate(R_e, a, reshape=False, order=1,
                            mode='reflect') if abs(a) > 0.01 else R_e
        v = ndimage.gaussian_filter1d(Rr.mean(axis=1), sigma=1.5)
        return np.sum(np.diff(v) ** 2)

    coarse = np.arange(-8, 8.01, 0.2)
    cs = np.array([score(a) for a in coarse])
    best_c = coarse[np.argmax(cs)]

    fine = np.arange(best_c - 0.5, best_c + 0.501, 0.02)
    fs = np.array([score(a) for a in fine])
    return fine[np.argmax(fs)]


# ──────────────────────────────────────────────
# Grid detection (peak-based)
# ──────────────────────────────────────────────
def detect_rows(R, expected_spacing=42):
    """Find rows from vertical projection. Works even on dim images
    because the chip's physical bands create contrast."""
    v_proj = R.mean(axis=1)
    v_smooth = ndimage.gaussian_filter1d(v_proj, sigma=3)
    rows, _ = signal.find_peaks(
        v_smooth, distance=int(expected_spacing * 0.7), prominence=0.5)
    return rows


def detect_spots_in_row(R, row_y, half_width=5):
    """Find luminous spot positions within a row using peak detection.
    min_distance=30 prevents detecting noise between wells (~45px apart)."""
    y0 = max(0, row_y - half_width)
    y1 = min(R.shape[0], row_y + half_width + 1)
    strip = R[y0:y1, :].mean(axis=0)
    strip_smooth = ndimage.gaussian_filter1d(strip, sigma=2)
    bg = np.median(strip_smooth)
    spots, _ = signal.find_peaks(
        strip_smooth, height=bg + 5, distance=15, prominence=2)
    return spots


def refine_center(R, y, x, radius=4):
    """Sub-pixel center refinement via center-of-mass."""
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


def fill_row_gaps(peaks, image_width, expected_period=46, margin=10):
    """Fill gaps between detected peaks to recover ALL well positions.
    Uses detected peaks as anchors, estimates well spacing from their
    closest-neighbor distances, then interpolates missing wells.
    expected_period: approximate well-to-well spacing (pixels)."""
    if len(peaks) < 2:
        return peaks

    peaks = np.sort(peaks)
    spacings = np.diff(peaks)

    # Only consider spacings that are plausible single-well distances
    # (0.6x to 1.5x expected period)
    lo, hi = expected_period * 0.6, expected_period * 1.5
    candidates = spacings[(spacings >= lo) & (spacings <= hi)]
    if len(candidates) >= 3:
        period = np.median(candidates)
    else:
        period = expected_period

    if period < 15:  # sanity check
        return peaks

    # Fill gaps between consecutive peaks
    all_pos = []
    for i in range(len(peaks)):
        all_pos.append(float(peaks[i]))
        if i < len(peaks) - 1:
            gap = peaks[i + 1] - peaks[i]
            n_missing = round(gap / period) - 1
            if n_missing > 0:
                step = gap / (n_missing + 1)
                for j in range(1, n_missing + 1):
                    all_pos.append(peaks[i] + j * step)

    # No edge extension — only fill between detected anchors

    return np.array(sorted(all_pos))


def merge_close_peaks(peaks, strip_values, min_sep=28):
    """Merge double-detections of the same well. Keeps the brightest peak
    in each cluster of peaks closer than min_sep."""
    if len(peaks) < 2:
        return peaks
    peaks = np.sort(peaks)
    merged = []
    i = 0
    while i < len(peaks):
        cluster = [peaks[i]]
        j = i + 1
        while j < len(peaks) and peaks[j] - cluster[0] < min_sep:
            cluster.append(peaks[j])
            j += 1
        cluster = np.array(cluster)
        best = cluster[np.argmax(strip_values[cluster])]
        merged.append(best)
        i = j
    return np.array(merged)


def detect_grid(R):
    """Detect all well positions on one image.
    1. Peak detection finds bright anchors (may include double-peaks)
    2. Merge double-peaks (< 28px apart)
    3. Fill gaps to recover quenched wells"""
    row_ys = detect_rows(R)
    H, W = R.shape
    all_spots = []
    counts = []
    for ri, ry in enumerate(row_ys):
        # Step 1: raw peak detection
        raw_xs = detect_spots_in_row(R, ry)

        # Step 2: merge double-peaks
        y0 = max(0, ry - 5)
        y1 = min(H, ry + 6)
        strip = R[y0:y1, :].mean(axis=0)
        strip_smooth = ndimage.gaussian_filter1d(strip, sigma=2)
        clean_xs = merge_close_peaks(raw_xs, strip_smooth)

        # Step 3: fill gaps between anchors
        all_xs = fill_row_gaps(clean_xs, W, expected_period=46)

        for sx in all_xs:
            cy, cx = refine_center(R, ry, int(round(sx)))
            all_spots.append((ri, cy, cx))
        counts.append(len(all_xs))

    counts = np.array(counts) if len(counts) > 0 else np.array([0])
    return all_spots, row_ys, counts


# ──────────────────────────────────────────────
# Intensity measurement
# ──────────────────────────────────────────────
def measure_intensity(R, cy, cx, spot_radius=3, bg_inner=6, bg_outer=10):
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


# ──────────────────────────────────────────────
# Block filtering
# ──────────────────────────────────────────────
def block_filter(spots, measurements, n_br=7, n_bc=7):
    if len(spots) == 0:
        return np.array([], dtype=bool)
    arr = np.array([(cy, cx) for _, cy, cx in spots])
    y_edges = np.linspace(arr[:, 0].min() - 1, arr[:, 0].max() + 1, n_br + 1)
    x_edges = np.linspace(arr[:, 1].min() - 1, arr[:, 1].max() + 1, n_bc + 1)

    block_ids = np.zeros(len(spots), dtype=int)
    for i, (_, cy, cx) in enumerate(spots):
        by = np.clip(np.searchsorted(y_edges, cy) - 1, 0, n_br - 1)
        bx = np.clip(np.searchsorted(x_edges, cx) - 1, 0, n_bc - 1)
        block_ids[i] = by * n_bc + bx

    keep = np.ones(len(spots), dtype=bool)
    bg_means = []
    for b in range(n_br * n_bc):
        mask = block_ids == b
        if mask.sum() < 5:
            keep[mask] = False
            bg_means.append(np.nan)
            continue
        bg_means.append(measurements[mask, 1].mean())

    bg_means = np.array(bg_means)
    valid = bg_means[~np.isnan(bg_means)]
    if len(valid) == 0:
        return keep
    med = np.median(valid)
    mad = np.median(np.abs(valid - med))
    thresh = med + 2.5 * mad * 1.4826

    for b in range(n_br * n_bc):
        if not np.isnan(bg_means[b]) and bg_means[b] > thresh:
            keep[block_ids == b] = False

    return keep


# ──────────────────────────────────────────────
# Coordinate mapping (rotated → original)
# ──────────────────────────────────────────────
def unrotate_spots(spots, angle_deg, img_shape):
    """Map spot coordinates from rotated image back to original image."""
    if abs(angle_deg) < 0.01:
        return spots
    theta = np.radians(angle_deg)
    cy_img, cx_img = img_shape[0] / 2.0, img_shape[1] / 2.0
    out = []
    for ri, y, x in spots:
        dy, dx = y - cy_img, x - cx_img
        y_orig = cy_img + dy * np.cos(theta) + dx * np.sin(theta)
        x_orig = cx_img - dy * np.sin(theta) + dx * np.cos(theta)
        out.append((ri, y_orig, x_orig))
    return out


# ──────────────────────────────────────────────
# Visualization (on original unrotated image)
# ──────────────────────────────────────────────
def save_crop(R, spots, keep, net_signals, bin_edges, filename,
              crop_center, crop_size=300, spot_radius=4):
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
    print("  independent per-image processing")
    print("=" * 60)

    # Process each image independently
    results = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        print(f"\n{'─'*40}")
        print(f"  Processing {name}")
        print(f"{'─'*40}")

        R_raw = load_red_channel(name)
        print(f"  Shape: {R_raw.shape}")

        # Rotation correction
        angle = detect_rotation_angle(R_raw)
        if abs(angle) > 0.01:
            R = ndimage.rotate(R_raw, angle, reshape=False, order=1,
                               mode='nearest')
        else:
            R = R_raw
        print(f"  Rotation: {angle:+.2f} deg")

        # Grid detection
        spots, row_ys, counts = detect_grid(R)
        n_rows = len(row_ys)
        n_spots = len(spots)
        print(f"  Detected: {n_rows} rows, {n_spots} luminous spots")
        if len(counts) > 0:
            print(f"  Spots/row: mean={counts.mean():.1f}, "
                  f"min={counts.min()}, max={counts.max()}")

        # Measure intensities
        meas = np.array([measure_intensity(R, cy, cx)
                         for _, cy, cx in spots]) if spots else np.zeros((0, 3))
        if n_spots > 0:
            nets = meas[:, 2]
            print(f"  Net signal: mean={nets.mean():.1f}, "
                  f"std={nets.std():.1f}, median={np.median(nets):.1f}")

        # Block filtering
        keep = block_filter(spots, meas) if n_spots > 0 else np.array([], dtype=bool)
        n_kept = keep.sum() if len(keep) > 0 else 0
        print(f"  After block filter: {n_kept}/{n_spots}")

        results[name] = {
            'angle': angle,
            'R_raw': R_raw,
            'R': R,
            'spots': spots,
            'meas': meas,
            'keep': keep,
            'n_rows': n_rows,
        }

    # Discretization based on 0ugL
    print(f"\n{'─'*40}")
    print("  Discretization & Distributions")
    print(f"{'─'*40}")

    ref = results["0ugL"]
    ref_net = ref['meas'][ref['keep'], 2]
    p5 = np.percentile(ref_net, 5)
    p10 = np.percentile(ref_net, 10)
    quench = max(abs(p5), abs(p10), 5.0)
    bright = ref_net[ref_net > quench]
    p90 = np.percentile(bright, 90) if len(bright) >= 20 else quench + 50
    r = p90 - quench
    bin_edges = np.array([-np.inf, quench, quench + r * 0.22,
                          quench + r * 0.50, quench + r * 0.80, np.inf])
    n_levels = 5
    print(f"  Quench threshold: {quench:.1f}")
    print(f"  Bin edges: {[f'{e:.1f}' for e in bin_edges[1:-1]]}")

    # Compute distributions (each image has all wells via gap filling)
    distributions = {}
    level_names = ["quench", "dim", "medium", "bright", "v.bright"]

    for name in ["0ugL", "5ugL", "20ugL"]:
        res = results[name]
        kept_net = res['meas'][res['keep'], 2] if res['keep'].sum() > 0 else np.array([])
        n_wells = len(kept_net)

        counts = np.zeros(n_levels)
        for i in range(n_levels):
            if i == 0:
                counts[i] = (kept_net <= bin_edges[1]).sum()
            elif i == n_levels - 1:
                counts[i] = (kept_net > bin_edges[-2]).sum()
            else:
                counts[i] = ((kept_net > bin_edges[i]) &
                             (kept_net <= bin_edges[i + 1])).sum()

        dist = counts / max(counts.sum(), 1)
        distributions[name] = dist
        parts = [f"{level_names[i]}={dist[i]:.4f}" for i in range(n_levels)]
        print(f"  {name} ({n_wells} wells): {', '.join(parts)}")

    # Stochastic dominance check
    print(f"\n  Stochastic dominance:")
    all_ok = True
    for k in range(n_levels):
        p0 = distributions["0ugL"][k:].sum()
        p5 = distributions["5ugL"][k:].sum()
        p20 = distributions["20ugL"][k:].sum()
        ok = (p0 >= p5 - 1e-10) and (p0 >= p20 - 1e-10) and (p5 >= p20 - 1e-10)
        if not ok:
            all_ok = False
        print(f"    k>={k}: 0ugL={p0:.4f}  5ugL={p5:.4f}  20ugL={p20:.4f}  "
              f"{'OK' if ok else 'FAIL'}")
    print(f"    Overall: {'PASS' if all_ok else 'FAIL'}")

    # Save visualizations
    print(f"\n{'─'*40}")
    print("  Saving visualizations")
    print(f"{'─'*40}")

    for name in ["0ugL", "5ugL", "20ugL"]:
        res = results[name]
        R_raw = res['R_raw']  # original unrotated image
        H, W = R_raw.shape
        # Map detected spots back to original image coordinates
        spots_orig = unrotate_spots(res['spots'], res['angle'], (H, W))
        keep = res['keep']
        net = res['meas'][:, 2] if len(res['meas']) > 0 else np.array([])

        crop_centers = [
            (int(H * 0.33), int(W * 0.40)),
            (int(H * 0.67), int(W * 0.65)),
        ]
        for ci, center in enumerate(crop_centers, 1):
            save_crop(R_raw, spots_orig, keep, net, bin_edges,
                      f"crop{ci}_{name}.png", center, crop_size=300)

    # Charts
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
        labels = ["Quenched", "Dim", "Medium", "Bright", "V.Bright"]
        bar_colors = ["#4444ff", "#00c8c8", "#00ff00", "#ffc800", "#ff3232"]
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
        "method": "v3: independent per-image rotation + peak detection + gap filling",
        "per_image": {},
        "bin_edges": [round(e, 1) if abs(e) < 1e10 else str(e) for e in bin_edges],
        "n_levels": n_levels,
        "distributions": {n: d.tolist() for n, d in distributions.items()},
    }
    for name in ["0ugL", "5ugL", "20ugL"]:
        res = results[name]
        output["per_image"][name] = {
            "rotation_deg": round(res['angle'], 3),
            "n_rows": res['n_rows'],
            "spots_detected": len(res['spots']),
            "spots_after_filter": int(res['keep'].sum()),
        }

    with open(os.path.join(DIR, "grid_analysis_results.json"), "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to grid_analysis_results.json")

    # Clean up diagnostic files
    for f in ["_diag_0ugL.png", "_diag_5ugL.png", "_diag_20ugL.png"]:
        p = os.path.join(DIR, f)
        if os.path.exists(p):
            os.remove(p)


if __name__ == "__main__":
    main()
