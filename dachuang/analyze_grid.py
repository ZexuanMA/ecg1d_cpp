#!/usr/bin/env python3
"""
Grid-based fluorescent spot intensity analysis.

Algorithm:
1. Detect horizontal row positions from vertical projection
2. Within each row, detect spot positions from horizontal profile
3. Refine spot centers via local center-of-mass
4. Measure spot intensity and local background at each grid node
5. Partition into rectangular blocks, assess quality per block
6. Exclude noisy/artifact blocks
7. Compute discretized intensity distributions (adaptive binning)
"""

import numpy as np
from PIL import Image
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
# 2. Detect row positions
# ──────────────────────────────────────────────
def detect_rows(R, expected_spacing=42, smooth_sigma=3):
    """Find row y-positions from vertical projection of red channel."""
    v_proj = R.mean(axis=1)
    v_smooth = ndimage.gaussian_filter1d(v_proj, sigma=smooth_sigma)
    row_peaks, props = signal.find_peaks(
        v_smooth,
        distance=int(expected_spacing * 0.7),
        prominence=1.5,
    )
    return row_peaks, v_smooth


# ──────────────────────────────────────────────
# 3. Detect spots within each row
# ──────────────────────────────────────────────
def detect_spots_in_row(R, row_y, half_width=3, smooth_sigma=2,
                        min_spacing=15, prominence=3):
    """Find spot x-positions within a horizontal strip around row_y."""
    y0 = max(0, row_y - half_width)
    y1 = min(R.shape[0], row_y + half_width + 1)
    strip = R[y0:y1, :].mean(axis=0)
    strip_smooth = ndimage.gaussian_filter1d(strip, sigma=smooth_sigma)

    bg = np.median(strip_smooth)
    spots, _ = signal.find_peaks(
        strip_smooth,
        height=bg + 5,
        distance=min_spacing,
        prominence=prominence,
    )
    return spots, strip_smooth


# ──────────────────────────────────────────────
# 4. Refine spot center via center-of-mass
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
    cy = (yy * patch_bg).sum() / total
    cx = (xx * patch_bg).sum() / total
    return cy, cx


# ──────────────────────────────────────────────
# 5. Measure spot intensity and local background
# ──────────────────────────────────────────────
def measure_intensity(R, cy, cx, spot_radius=3, bg_inner=6, bg_outer=10):
    """
    Spot intensity = mean in spot_radius circle.
    Background = mean in annulus [bg_inner, bg_outer].
    Net signal = spot - background.
    """
    H, W = R.shape
    iy, ix = int(round(cy)), int(round(cx))

    # Spot intensity (square patch for speed)
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
    bg_patch = R[by0:by1, bx0:bx1]
    yy, xx = np.mgrid[by0:by1, bx0:bx1]
    dist_sq = (yy - cy) ** 2 + (xx - cx) ** 2
    mask = (dist_sq >= bg_inner**2) & (dist_sq <= bg_outer**2)
    if mask.sum() > 0:
        bg_val = bg_patch[mask[(yy - by0), (xx - bx0)]].mean() if mask.any() else spot_val
        # Simpler approach
        bg_pixels = R[by0:by1, bx0:bx1]
        # Build mask in local coords
        local_yy = np.arange(by1 - by0)[:, None]
        local_xx = np.arange(bx1 - bx0)[None, :]
        local_dist_sq = (local_yy - (cy - by0))**2 + (local_xx - (cx - bx0))**2
        local_mask = (local_dist_sq >= bg_inner**2) & (local_dist_sq <= bg_outer**2)
        if local_mask.sum() > 5:
            bg_val = bg_pixels[local_mask].mean()
        else:
            bg_val = np.median(bg_pixels)
    else:
        bg_val = np.median(R[by0:by1, bx0:bx1])

    return spot_val, bg_val, spot_val - bg_val


# ──────────────────────────────────────────────
# 6. Full grid detection on reference image
# ──────────────────────────────────────────────
def detect_grid(R_ref):
    """Detect all spots on the reference (0ugL) image. Returns list of (row_idx, cy, cx)."""
    row_ys, v_smooth = detect_rows(R_ref)
    print(f"  Detected {len(row_ys)} rows, spacing ~{np.median(np.diff(row_ys)):.0f} px")

    all_spots = []
    spots_per_row = []
    for ri, ry in enumerate(row_ys):
        spot_xs, _ = detect_spots_in_row(R_ref, ry)
        row_spots = []
        for sx in spot_xs:
            cy, cx = refine_center(R_ref, ry, sx)
            row_spots.append((ri, cy, cx))
        all_spots.extend(row_spots)
        spots_per_row.append(len(row_spots))

    spots_per_row = np.array(spots_per_row)
    print(f"  Spots per row: mean={spots_per_row.mean():.1f}, "
          f"min={spots_per_row.min()}, max={spots_per_row.max()}")
    print(f"  Total spots: {len(all_spots)}")
    return all_spots, row_ys, spots_per_row


# ──────────────────────────────────────────────
# 7. Measure all spots on all images
# ──────────────────────────────────────────────
def measure_all_spots(images, spots):
    """Measure intensity for each spot in each image.
    Returns dict: name -> array of (spot_val, bg_val, net_signal)."""
    results = {}
    for name, R in images.items():
        data = []
        for ri, cy, cx in spots:
            sv, bv, net = measure_intensity(R, cy, cx)
            data.append((sv, bv, net))
        results[name] = np.array(data)
        nets = results[name][:, 2]
        print(f"  {name}: net signal mean={nets.mean():.1f}, "
              f"std={nets.std():.1f}, median={np.median(nets):.1f}")
    return results


# ──────────────────────────────────────────────
# 8. Block partitioning and quality filtering
# ──────────────────────────────────────────────
def partition_into_blocks(spots, row_ys, n_block_rows=7, n_block_cols=7):
    """Divide spots into rectangular blocks for quality assessment."""
    spots_arr = np.array([(cy, cx) for _, cy, cx in spots])
    y_min, y_max = spots_arr[:, 0].min(), spots_arr[:, 0].max()
    x_min, x_max = spots_arr[:, 1].min(), spots_arr[:, 1].max()

    y_edges = np.linspace(y_min - 1, y_max + 1, n_block_rows + 1)
    x_edges = np.linspace(x_min - 1, x_max + 1, n_block_cols + 1)

    block_ids = np.zeros(len(spots), dtype=int)
    for i, (_, cy, cx) in enumerate(spots):
        by = np.searchsorted(y_edges, cy) - 1
        bx = np.searchsorted(x_edges, cx) - 1
        by = np.clip(by, 0, n_block_rows - 1)
        bx = np.clip(bx, 0, n_block_cols - 1)
        block_ids[i] = by * n_block_cols + bx

    return block_ids, y_edges, x_edges


def filter_blocks(block_ids, ref_data, n_blocks, min_spots=5,
                  bg_outlier_sigma=2.5, net_cv_max=1.5):
    """
    Exclude blocks that are noisy or have artifacts.
    Criteria:
      - Too few spots
      - Background too high or too variable (artifact region)
      - Net signal coefficient-of-variation too high
    Returns: boolean mask (True = keep).
    """
    keep = np.ones(len(block_ids), dtype=bool)

    # Compute per-block statistics
    block_bg_means = []
    block_bg_stds = []
    block_counts = []
    for b in range(n_blocks):
        mask = block_ids == b
        count = mask.sum()
        block_counts.append(count)
        if count < min_spots:
            keep[mask] = False
            block_bg_means.append(np.nan)
            block_bg_stds.append(np.nan)
            continue
        bg = ref_data[mask, 1]  # background values
        block_bg_means.append(bg.mean())
        block_bg_stds.append(bg.std())

    block_bg_means = np.array(block_bg_means)
    block_bg_stds = np.array(block_bg_stds)
    block_counts = np.array(block_counts)

    # Remove blocks with outlier background (too bright = artifact/contamination)
    valid_bgs = block_bg_means[~np.isnan(block_bg_means)]
    bg_median = np.median(valid_bgs)
    bg_mad = np.median(np.abs(valid_bgs - bg_median))
    bg_threshold = bg_median + bg_outlier_sigma * bg_mad * 1.4826

    n_excluded = 0
    for b in range(n_blocks):
        mask = block_ids == b
        if np.isnan(block_bg_means[b]):
            continue
        if block_bg_means[b] > bg_threshold:
            keep[mask] = False
            n_excluded += 1
        # Also check if background is too variable within the block
        if block_bg_stds[b] > bg_mad * 1.4826 * 3:
            keep[mask] = False
            n_excluded += 1

    print(f"  Block filtering: {n_excluded} blocks excluded, "
          f"{keep.sum()}/{len(keep)} spots retained "
          f"({100*keep.sum()/len(keep):.1f}%)")
    return keep


# ──────────────────────────────────────────────
# 9. Adaptive intensity discretization
# ──────────────────────────────────────────────
def discretize_intensity(ref_net, method="auto"):
    """
    Determine bin edges from the reference (0ugL) net signal distribution.
    Uses physically motivated thresholds based on distribution shape analysis.
    Returns: bin_edges, n_levels

    The distribution typically has:
    - A noise/quenched peak near zero (net ≤ ~8)
    - A broad dim region (8-25)
    - A medium region (25-45)
    - A bright tail (45-70, >70)
    """
    # Estimate noise level from the lower tail
    p5 = np.percentile(ref_net, 5)
    p10 = np.percentile(ref_net, 10)
    noise_estimate = max(abs(p5), abs(p10), 5.0)

    # Quench threshold: where signal is indistinguishable from noise
    quench_threshold = noise_estimate

    bright = ref_net[ref_net > quench_threshold]
    if len(bright) < 20:
        return np.array([-np.inf, quench_threshold, np.inf]), 2

    # Use physically motivated bin edges for 5 levels:
    # Level 0: quenched  (≤ quench_threshold)
    # Level 1: dim       (quench_threshold, edge1]
    # Level 2: medium    (edge1, edge2]
    # Level 3: bright    (edge2, edge3]
    # Level 4: very bright (> edge3)
    #
    # Edges based on signal range, roughly equal-width in log-ish space
    p90 = np.percentile(bright, 90)

    # Divide the range [quench_threshold, p90] into intervals
    # that capture the decay shape of the distribution
    sig_range = p90 - quench_threshold
    edge1 = quench_threshold + sig_range * 0.22  # dim → medium boundary
    edge2 = quench_threshold + sig_range * 0.50  # medium → bright boundary
    edge3 = quench_threshold + sig_range * 0.80  # bright → very bright boundary

    bin_edges = np.array([-np.inf, quench_threshold, edge1, edge2, edge3, np.inf])
    n_levels = len(bin_edges) - 1

    # Verify no level is too sparse (< 3%)
    for i in range(n_levels):
        if i == 0:
            frac = (ref_net <= bin_edges[1]).sum() / len(ref_net)
        elif i == n_levels - 1:
            frac = (ref_net > bin_edges[-2]).sum() / len(ref_net)
        else:
            frac = ((ref_net > bin_edges[i]) & (ref_net <= bin_edges[i+1])).sum() / len(ref_net)
        if frac < 0.03:
            # Merge with adjacent level — fall back to 4 levels
            bin_edges = np.array([-np.inf, quench_threshold, p50, p90, np.inf])
            n_levels = 4
            break

    print(f"  Discretization: {n_levels} levels, quench threshold = {quench_threshold:.1f}")
    print(f"  Bin edges: {[f'{e:.1f}' if abs(e) < 1e10 else str(e) for e in bin_edges]}")
    return bin_edges, n_levels


def compute_distribution(net_signals, bin_edges):
    """Compute discrete probability distribution from net signals."""
    counts = np.zeros(len(bin_edges) - 1)
    for i in range(len(bin_edges) - 1):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        counts[i] = ((net_signals > lo) & (net_signals <= hi)).sum()
    # Handle edge: signals exactly at -inf edge
    counts[0] += (net_signals <= bin_edges[1]).sum() - counts[0]
    # Recalculate properly
    counts = np.zeros(len(bin_edges) - 1)
    for i in range(len(bin_edges) - 1):
        if i == 0:
            mask = net_signals <= bin_edges[1]
        elif i == len(bin_edges) - 2:
            mask = net_signals > bin_edges[-2]
        else:
            mask = (net_signals > bin_edges[i]) & (net_signals <= bin_edges[i + 1])
        counts[i] = mask.sum()

    total = counts.sum()
    if total > 0:
        dist = counts / total
    else:
        dist = counts
    return dist, counts


# ──────────────────────────────────────────────
# 10. Visualization
# ──────────────────────────────────────────────
def save_marked_image(R, spots, keep_mask, filename, spot_radius=5,
                      net_signals=None, bin_edges=None):
    """Save image with circles marking detected spots.
    If net_signals and bin_edges provided, color-code by intensity level.
    Otherwise: Green=kept, Red=excluded."""
    from PIL import ImageDraw

    # Normalize to 0-255 and convert to RGB (warm tone for fluorescence)
    R_norm = ((R - R.min()) / max(R.max() - R.min(), 1) * 220).astype(np.uint8)
    img_rgb = np.stack([R_norm, R_norm // 4, R_norm // 6], axis=2)
    pil_img = Image.fromarray(img_rgb)
    draw = ImageDraw.Draw(pil_img)

    # Color palette for levels (level 0=blue/quenched, up to level 4=bright yellow)
    level_colors = [
        (80, 80, 255),    # 0: blue (quenched)
        (0, 200, 200),    # 1: cyan (dim)
        (0, 255, 0),      # 2: green (medium)
        (255, 200, 0),    # 3: yellow (bright)
        (255, 50, 50),    # 4: red (very bright)
        (255, 0, 255),    # 5+: magenta
    ]

    for i, (_, cy, cx) in enumerate(spots):
        iy, ix = int(round(cy)), int(round(cx))
        if not keep_mask[i]:
            color = (60, 60, 60)  # gray for excluded
        elif net_signals is not None and bin_edges is not None:
            val = net_signals[i]
            level = np.searchsorted(bin_edges, val) - 1
            level = np.clip(level, 0, len(level_colors) - 1)
            color = level_colors[level]
        else:
            color = (0, 255, 0)
        draw.ellipse(
            [ix - spot_radius, iy - spot_radius,
             ix + spot_radius, iy + spot_radius],
            outline=color, width=1,
        )
    pil_img.save(os.path.join(DIR, filename))
    print(f"  Saved {filename}")


def save_block_map(spots, block_ids, keep_mask, R_shape, filename,
                   y_edges, x_edges):
    """Save block partition map showing which blocks are kept/excluded."""
    from PIL import ImageDraw

    canvas = np.zeros((*R_shape, 3), dtype=np.uint8)
    pil_img = Image.fromarray(canvas)
    draw = ImageDraw.Draw(pil_img)

    n_block_rows = len(y_edges) - 1
    n_block_cols = len(x_edges) - 1

    for br in range(n_block_rows):
        for bc in range(n_block_cols):
            bid = br * n_block_cols + bc
            mask = block_ids == bid
            n_kept = (mask & keep_mask).sum()
            n_total = mask.sum()

            y0, y1 = int(y_edges[br]), int(y_edges[br + 1])
            x0, x1 = int(x_edges[bc]), int(x_edges[bc + 1])

            if n_total == 0:
                color = (40, 40, 40)
            elif n_kept == 0:
                color = (100, 0, 0)
            else:
                green = int(200 * n_kept / n_total)
                color = (0, green, 0)

            draw.rectangle([x0, y0, x1, y1], fill=color, outline=(80, 80, 80))
            # Label
            draw.text((x0 + 5, y0 + 5), f"{n_kept}/{n_total}", fill=(255, 255, 255))

    pil_img.save(os.path.join(DIR, filename))
    print(f"  Saved {filename}")


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    print("=" * 60)
    print("Grid-based Fluorescent Spot Intensity Analysis")
    print("=" * 60)

    # Load images
    print("\n[1] Loading images...")
    images = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        images[name] = load_red_channel(name)
        print(f"  {name}: shape={images[name].shape}, "
              f"R range=[{images[name].min():.0f}, {images[name].max():.0f}]")

    R_ref = images["0ugL"]

    # Detect grid on reference image
    print("\n[2] Detecting grid structure on 0ugL...")
    spots, row_ys, spots_per_row = detect_grid(R_ref)

    # Measure intensities on all images
    print("\n[3] Measuring spot intensities...")
    measurements = measure_all_spots(images, spots)

    # Block partitioning
    print("\n[4] Block partitioning and quality filtering...")
    n_br, n_bc = 7, 7
    block_ids, y_edges, x_edges = partition_into_blocks(
        spots, row_ys, n_block_rows=n_br, n_block_cols=n_bc
    )
    n_blocks = n_br * n_bc
    keep = filter_blocks(block_ids, measurements["0ugL"], n_blocks)

    # Show per-block stats
    print("\n  Block quality summary:")
    for b in range(n_blocks):
        mask = (block_ids == b)
        n_total = mask.sum()
        n_kept = (mask & keep).sum()
        if n_total > 0:
            br, bc = b // n_bc, b % n_bc
            ref_net = measurements["0ugL"][mask, 2]
            bg_mean = measurements["0ugL"][mask, 1].mean()
            status = "KEEP" if n_kept > 0 else "EXCL"
            if n_kept == 0:
                print(f"    Block ({br},{bc}): {n_total:3d} spots, "
                      f"bg={bg_mean:.1f}, [{status}]")

    # Apply filter
    ref_net_filtered = measurements["0ugL"][keep, 2]

    # Intensity discretization
    print("\n[5] Discretizing intensity levels...")
    bin_edges, n_levels = discretize_intensity(ref_net_filtered)

    # Compute distributions for all images
    print("\n[6] Computing intensity distributions...")
    distributions = {}
    for name in ["0ugL", "5ugL", "20ugL"]:
        net = measurements[name][keep, 2]
        dist, counts = compute_distribution(net, bin_edges)
        distributions[name] = dist
        print(f"\n  {name}:")
        for i in range(n_levels):
            lo = bin_edges[i] if bin_edges[i] > -1e10 else "-inf"
            hi = bin_edges[i + 1] if bin_edges[i + 1] < 1e10 else "+inf"
            print(f"    Level {i} ({lo} ~ {hi}]: "
                  f"count={int(counts[i])}, prob={dist[i]:.4f}")

    # Stochastic dominance check
    print("\n[7] Stochastic dominance check (P(level >= k))...")
    print(f"  {'k':<4} {'0ugL':>8} {'5ugL':>8} {'20ugL':>8} "
          f"{'0>5?':>6} {'0>20?':>6} {'5>20?':>6}")
    for k in range(n_levels):
        p0 = distributions["0ugL"][k:].sum()
        p5 = distributions["5ugL"][k:].sum()
        p20 = distributions["20ugL"][k:].sum()
        c05 = "OK" if p0 >= p5 - 1e-10 else "FAIL"
        c020 = "OK" if p0 >= p20 - 1e-10 else "FAIL"
        c520 = "OK" if p5 >= p20 - 1e-10 else "FAIL"
        print(f"  {k:<4} {p0:>8.4f} {p5:>8.4f} {p20:>8.4f} "
              f"{c05:>6} {c020:>6} {c520:>6}")

    # Save visualizations
    print("\n[8] Saving visualizations...")

    # Full image with level-colored spots for each sample
    for name in ["0ugL", "5ugL", "20ugL"]:
        net = measurements[name][:, 2]
        save_marked_image(
            images[name], spots, keep,
            f"grid_marked_{name}.png", spot_radius=5,
            net_signals=net, bin_edges=bin_edges,
        )

    save_block_map(spots, block_ids, keep, R_ref.shape, "block_map.png",
                   y_edges, x_edges)

    # Save crop for inspection (center region, all three samples side by side)
    cy_center = R_ref.shape[0] // 2
    cx_center = R_ref.shape[1] // 2
    crop_size = 300
    for name in ["0ugL", "5ugL", "20ugL"]:
        R_img = images[name]
        R_crop = R_img[cy_center - crop_size:cy_center + crop_size,
                       cx_center - crop_size:cx_center + crop_size]
        crop_spots = []
        crop_nets = []
        crop_keep_list = []
        for idx, (ri, cy, cx) in enumerate(spots):
            if abs(cy - cy_center) < crop_size and abs(cx - cx_center) < crop_size:
                crop_spots.append((ri, cy - cy_center + crop_size,
                                   cx - cx_center + crop_size))
                crop_nets.append(measurements[name][idx, 2])
                crop_keep_list.append(keep[idx])
        crop_nets = np.array(crop_nets) if crop_nets else np.array([])
        crop_keep_arr = np.array(crop_keep_list) if crop_keep_list else np.array([], dtype=bool)
        save_marked_image(
            R_crop, crop_spots, crop_keep_arr,
            f"grid_crop_{name}.png", spot_radius=4,
            net_signals=crop_nets, bin_edges=bin_edges,
        )

    # Output for Ridge-Projection
    print("\n" + "=" * 60)
    print("Ridge-Projection Input")
    print("=" * 60)
    for pair in [("0ugL", "5ugL"), ("0ugL", "20ugL")]:
        X = distributions[pair[0]]
        Y = distributions[pair[1]]
        print(f"\n  {pair[0]} -> {pair[1]}:")
        print(f"  X = {np.array2string(X, precision=6, separator=', ')}")
        print(f"  Y = {np.array2string(Y, precision=6, separator=', ')}")

    # Save results
    output = {
        "method": "Grid-based detection with block quality filtering",
        "grid": {
            "n_rows": int(len(row_ys)),
            "row_spacing_px": float(np.median(np.diff(row_ys))),
            "total_spots_detected": len(spots),
            "spots_after_filtering": int(keep.sum()),
            "block_grid": f"{n_br}x{n_bc}",
        },
        "bin_edges": bin_edges.tolist(),
        "n_levels": n_levels,
        "distributions": {
            name: dist.tolist() for name, dist in distributions.items()
        },
    }
    out_path = os.path.join(DIR, "grid_analysis_results.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
