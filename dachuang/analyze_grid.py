#!/usr/bin/env python3
"""
Grid-based fluorescent spot intensity analysis with automatic rotation correction.

Algorithm:
1. Detect image tilt from vertical projection contrast maximization
2. Correct rotation on all images
3. Detect horizontal row positions from vertical projection
4. Within each row, detect spot positions from horizontal profile
5. Refine spot centers via local center-of-mass
6. Measure spot intensity and local background at each grid node
7. Partition into rectangular blocks, assess quality per block
8. Exclude noisy/artifact blocks
9. Compute discretized intensity distributions
"""

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage, signal
import json, os

DIR = os.path.dirname(os.path.abspath(__file__))


# ──────────────────────────────────────────────
# 1. Load image, extract red channel
# ──────────────────────────────────────────────
def load_red_channel(name):
    path = os.path.join(DIR, f"{name}.bmp")
    img = np.array(Image.open(path))
    return img[:, :, 0].astype(np.float64)


# ──────────────────────────────────────────────
# 2. Rotation detection and correction
# ──────────────────────────────────────────────
def detect_rotation_angle(R):
    """Find tilt angle using two-stage search on gradient energy of vertical projection.
    Uses 2x downsampling and center crop to avoid edge artifacts."""
    scale = 2
    R_small = R[::scale, ::scale]
    H, W = R_small.shape
    # Crop center to avoid edge fill artifacts when rotating
    margin = int(min(H, W) * 0.12)
    R_crop = R_small[margin:H - margin, margin:W - margin]

    def score_angle(a):
        if abs(a) < 0.001:
            R_rot = R_crop
        else:
            R_rot = ndimage.rotate(R_crop, a, reshape=False, order=1, mode='reflect')
        v = R_rot.mean(axis=1)
        v_s = ndimage.gaussian_filter1d(v, sigma=1.5)
        # Gradient energy: sharper row peaks = larger squared gradients
        return np.sum(np.diff(v_s) ** 2)

    # Stage 1: coarse search -5° to +5° in 0.2° steps
    coarse = np.arange(-5, 5.01, 0.2)
    coarse_scores = np.array([score_angle(a) for a in coarse])
    best_coarse = coarse[np.argmax(coarse_scores)]

    # Stage 2: fine search around best ±0.5° in 0.02° steps
    fine = np.arange(best_coarse - 0.5, best_coarse + 0.501, 0.02)
    fine_scores = np.array([score_angle(a) for a in fine])
    best_fine = fine[np.argmax(fine_scores)]

    return best_fine


# ──────────────────────────────────────────────
# 3. Detect row positions
# ──────────────────────────────────────────────
def detect_rows(R, expected_spacing=42, smooth_sigma=3):
    """Find row y-positions from vertical projection of red channel."""
    v_proj = R.mean(axis=1)
    v_smooth = ndimage.gaussian_filter1d(v_proj, sigma=smooth_sigma)
    row_peaks, _ = signal.find_peaks(
        v_smooth,
        distance=int(expected_spacing * 0.7),
        prominence=1.5,
    )
    return row_peaks, v_smooth


# ──────────────────────────────────────────────
# 4. Detect spots within each row
# ──────────────────────────────────────────────
def detect_spots_in_row(R, row_y, half_width=5, smooth_sigma=2,
                        min_spacing=15, prominence=2):
    """Find spot x-positions within a horizontal strip around row_y."""
    y0 = max(0, row_y - half_width)
    y1 = min(R.shape[0], row_y + half_width + 1)
    strip = R[y0:y1, :].mean(axis=0)
    strip_smooth = ndimage.gaussian_filter1d(strip, sigma=smooth_sigma)
    bg = np.median(strip_smooth)
    spots, _ = signal.find_peaks(
        strip_smooth, height=bg + 5, distance=min_spacing, prominence=prominence,
    )
    return spots, strip_smooth


# ──────────────────────────────────────────────
# 5. Refine spot center via center-of-mass
# ──────────────────────────────────────────────
def refine_center(R, y, x, radius=4):
    """Sub-pixel center refinement using center-of-mass."""
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
# 6. Measure spot intensity and local background
# ──────────────────────────────────────────────
def measure_intensity(R, cy, cx, spot_radius=3, bg_inner=6, bg_outer=10):
    """Net signal = spot mean (7x7) - background mean (annulus r=6..10)."""
    H, W = R.shape
    iy, ix = int(round(cy)), int(round(cx))

    # Spot intensity
    sy0 = max(0, iy - spot_radius)
    sy1 = min(H, iy + spot_radius + 1)
    sx0 = max(0, ix - spot_radius)
    sx1 = min(W, ix + spot_radius + 1)
    spot_val = R[sy0:sy1, sx0:sx1].mean()

    # Background annulus
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
# 7. Full grid detection
# ──────────────────────────────────────────────
def detect_grid(R_ref):
    """Detect all spots on the reference image. Returns list of (row_idx, cy, cx)."""
    row_ys, _ = detect_rows(R_ref)
    print(f"  Detected {len(row_ys)} rows, spacing ~{np.median(np.diff(row_ys)):.0f} px")

    all_spots = []
    spots_per_row = []
    for ri, ry in enumerate(row_ys):
        spot_xs, _ = detect_spots_in_row(R_ref, ry)
        for sx in spot_xs:
            cy, cx = refine_center(R_ref, ry, sx)
            all_spots.append((ri, cy, cx))
        spots_per_row.append(len(spot_xs))

    spots_per_row = np.array(spots_per_row)
    print(f"  Spots/row: mean={spots_per_row.mean():.1f}, "
          f"min={spots_per_row.min()}, max={spots_per_row.max()}")
    print(f"  Total: {len(all_spots)}")
    return all_spots, row_ys, spots_per_row


# ──────────────────────────────────────────────
# 8. Measure all spots on all images
# ──────────────────────────────────────────────
def measure_all_spots(images, spots):
    """Returns dict: name -> array of (spot_val, bg_val, net_signal)."""
    results = {}
    for name, R in images.items():
        data = [measure_intensity(R, cy, cx) for _, cy, cx in spots]
        results[name] = np.array(data)
        nets = results[name][:, 2]
        print(f"  {name}: mean={nets.mean():.1f}, std={nets.std():.1f}, "
              f"median={np.median(nets):.1f}")
    return results


# ──────────────────────────────────────────────
# 9. Block partitioning and quality filtering
# ──────────────────────────────────────────────
def partition_into_blocks(spots, row_ys, n_block_rows=7, n_block_cols=7):
    spots_arr = np.array([(cy, cx) for _, cy, cx in spots])
    y_edges = np.linspace(spots_arr[:, 0].min() - 1, spots_arr[:, 0].max() + 1,
                          n_block_rows + 1)
    x_edges = np.linspace(spots_arr[:, 1].min() - 1, spots_arr[:, 1].max() + 1,
                          n_block_cols + 1)
    block_ids = np.zeros(len(spots), dtype=int)
    for i, (_, cy, cx) in enumerate(spots):
        by = np.clip(np.searchsorted(y_edges, cy) - 1, 0, n_block_rows - 1)
        bx = np.clip(np.searchsorted(x_edges, cx) - 1, 0, n_block_cols - 1)
        block_ids[i] = by * n_block_cols + bx
    return block_ids, y_edges, x_edges


def filter_blocks(block_ids, ref_data, n_blocks, min_spots=5,
                  bg_outlier_sigma=2.5):
    """Exclude blocks with too few spots, outlier background, or high bg variance."""
    keep = np.ones(len(block_ids), dtype=bool)
    bg_means = []
    bg_stds = []
    for b in range(n_blocks):
        mask = block_ids == b
        if mask.sum() < min_spots:
            keep[mask] = False
            bg_means.append(np.nan)
            bg_stds.append(np.nan)
            continue
        bg = ref_data[mask, 1]
        bg_means.append(bg.mean())
        bg_stds.append(bg.std())

    bg_means = np.array(bg_means)
    bg_stds = np.array(bg_stds)
    valid = bg_means[~np.isnan(bg_means)]
    med = np.median(valid)
    mad = np.median(np.abs(valid - med))
    thresh = med + bg_outlier_sigma * mad * 1.4826

    n_excl = 0
    for b in range(n_blocks):
        if np.isnan(bg_means[b]):
            continue
        mask = block_ids == b
        if bg_means[b] > thresh or bg_stds[b] > mad * 1.4826 * 3:
            keep[mask] = False
            n_excl += 1

    print(f"  {n_excl} blocks excluded, {keep.sum()}/{len(keep)} spots kept "
          f"({100*keep.sum()/len(keep):.1f}%)")
    return keep


# ──────────────────────────────────────────────
# 10. Intensity discretization
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
    edges = np.array([-np.inf, quench, quench + r*0.22, quench + r*0.50,
                      quench + r*0.80, np.inf])
    n_levels = len(edges) - 1

    print(f"  {n_levels} levels, quench={quench:.1f}")
    print(f"  Bin edges: {[f'{e:.1f}' for e in edges[1:-1]]}")
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
                         (net_signals <= bin_edges[i+1])).sum()
    total = counts.sum()
    return (counts / total if total > 0 else counts), counts


# ──────────────────────────────────────────────
# 11. Visualization — crop images only
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
        (80, 80, 255),    # 0: quenched (blue)
        (0, 200, 200),    # 1: dim (cyan)
        (0, 255, 0),      # 2: medium (green)
        (255, 200, 0),    # 3: bright (yellow)
        (255, 50, 50),    # 4: very bright (red)
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
    print(f"  Saved {filename}")


def save_block_map(spots, block_ids, keep, shape, filename, y_edges, x_edges):
    n_br = len(y_edges) - 1
    n_bc = len(x_edges) - 1
    canvas = np.zeros((*shape, 3), dtype=np.uint8)
    img = Image.fromarray(canvas)
    draw = ImageDraw.Draw(img)
    for br in range(n_br):
        for bc in range(n_bc):
            bid = br * n_bc + bc
            mask = block_ids == bid
            n_kept = (mask & keep).sum()
            n_total = mask.sum()
            y0, y1 = int(y_edges[br]), int(y_edges[br+1])
            x0, x1 = int(x_edges[bc]), int(x_edges[bc+1])
            if n_total == 0:
                color = (40, 40, 40)
            elif n_kept == 0:
                color = (100, 0, 0)
            else:
                color = (0, int(200 * n_kept / n_total), 0)
            draw.rectangle([x0, y0, x1, y1], fill=color, outline=(80, 80, 80))
            draw.text((x0 + 5, y0 + 5), f"{n_kept}/{n_total}", fill=(255, 255, 255))
    img.save(os.path.join(DIR, filename))
    print(f"  Saved {filename}")


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    print("=" * 60)
    print("Grid-based Fluorescent Spot Intensity Analysis")
    print("  with automatic rotation correction")
    print("=" * 60)

    # Load images
    print("\n[1] Loading images...")
    raw_images = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        raw_images[name] = load_red_channel(name)
        print(f"  {name}: {raw_images[name].shape}")

    # Detect and correct rotation
    print("\n[2] Detecting rotation angle...")
    angle = detect_rotation_angle(raw_images["0ugL"])
    print(f"  Detected tilt: {angle:.2f} deg")

    images = {}
    if abs(angle) > 0.01:
        for name, R in raw_images.items():
            images[name] = ndimage.rotate(R, angle, reshape=False, order=1,
                                          mode='nearest')
        print(f"  All images rotated by {angle:.2f} deg")
    else:
        images = raw_images
        print("  No rotation needed")

    R_ref = images["0ugL"]

    # Grid detection
    print("\n[3] Detecting grid structure on 0ugL...")
    spots, row_ys, spots_per_row = detect_grid(R_ref)

    # Measure intensities
    print("\n[4] Measuring spot intensities...")
    measurements = measure_all_spots(images, spots)

    # Block filtering
    print("\n[5] Block partitioning and quality filtering...")
    n_br, n_bc = 7, 7
    block_ids, y_edges, x_edges = partition_into_blocks(
        spots, row_ys, n_block_rows=n_br, n_block_cols=n_bc)
    keep = filter_blocks(block_ids, measurements["0ugL"], n_br * n_bc)

    # Excluded blocks detail
    for b in range(n_br * n_bc):
        mask = block_ids == b
        if mask.sum() > 0 and (mask & keep).sum() == 0:
            br, bc = b // n_bc, b % n_bc
            bg = measurements["0ugL"][mask, 1].mean()
            print(f"    Excluded block ({br},{bc}): {mask.sum()} spots, bg={bg:.1f}")

    # Discretize
    ref_net_filtered = measurements["0ugL"][keep, 2]
    print("\n[6] Discretizing intensity levels...")
    bin_edges, n_levels = discretize_intensity(ref_net_filtered)

    # Distributions
    print("\n[7] Computing distributions...")
    distributions = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        net = measurements[name][keep, 2]
        dist, counts = compute_distribution(net, bin_edges)
        distributions[name] = dist
        level_names = ["quench", "dim", "medium", "bright", "v.bright"][:n_levels]
        parts = [f"{level_names[i]}={dist[i]:.4f}" for i in range(n_levels)]
        print(f"  {name}: {', '.join(parts)}")

    # Stochastic dominance
    print("\n[8] Stochastic dominance check...")
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

    # Save 2 crops per concentration
    print("\n[9] Saving crop visualizations (2 per concentration)...")
    H, W = R_ref.shape
    crop_centers = [
        (int(H * 0.33), int(W * 0.40)),   # upper-left area
        (int(H * 0.67), int(W * 0.65)),   # lower-right area
    ]
    for name in ["0ugL", "5ugL", "20ugL"]:
        net = measurements[name][:, 2]
        for ci, center in enumerate(crop_centers, 1):
            save_crop(images[name], spots, keep, net, bin_edges,
                      f"crop{ci}_{name}.png", center, crop_size=300)

    # Block map
    save_block_map(spots, block_ids, keep, R_ref.shape, "block_map.png",
                   y_edges, x_edges)

    # Distribution bar chart and stochastic dominance plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Bar chart
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
        fig.suptitle("Discrete Intensity Distribution", fontsize=13)
        fig.tight_layout()
        fig.savefig(os.path.join(DIR, "distributions_bar.png"), dpi=150)
        plt.close()
        print("  Saved distributions_bar.png")

        # Stochastic dominance
        fig, ax = plt.subplots(figsize=(7, 5))
        styles = [("0ugL", "o-", "#4488ff"), ("5ugL", "s-", "#ff9900"),
                  ("20ugL", "D-", "#ff4444")]
        for name, marker, color in styles:
            surv = [distributions[name][k:].sum() for k in range(n_levels)]
            ax.plot(range(n_levels), surv, marker, label=name, color=color,
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

    # Save JSON results
    output = {
        "method": "Grid-based detection with rotation correction and block filtering",
        "rotation_angle_deg": round(angle, 3),
        "grid": {
            "n_rows": int(len(row_ys)),
            "row_spacing_px": float(np.median(np.diff(row_ys))),
            "total_spots": len(spots),
            "spots_after_filter": int(keep.sum()),
            "block_grid": f"{n_br}x{n_bc}",
        },
        "bin_edges": [round(e, 1) if abs(e) < 1e10 else str(e) for e in bin_edges],
        "n_levels": n_levels,
        "distributions": {n: d.tolist() for n, d in distributions.items()},
    }
    with open(os.path.join(DIR, "grid_analysis_results.json"), "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to grid_analysis_results.json")


if __name__ == "__main__":
    main()
