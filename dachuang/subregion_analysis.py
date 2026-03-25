#!/usr/bin/env python3
"""
Subregion distribution uniformity test.

Split the 0ugL image into NxM sub-regions and compute the intensity
distribution in each. If the detection algorithm is correct, distributions
should be consistent across sub-regions (0ugL has no quencher).

Also processes 20ugL for comparison.
"""

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage
import os, sys

# Reuse detection functions from analyze_grid
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_grid import (
    load_red_channel, detect_rotation_angle, detect_grid,
    measure_intensity, block_filter
)

DIR = os.path.dirname(os.path.abspath(__file__))


def analyze_subregions(name, n_rows_grid=2, n_cols_grid=4):
    """Analyze one image, split into sub-regions, return per-region stats."""
    print(f"\n{'='*60}")
    print(f"  Subregion analysis: {name}  ({n_rows_grid}x{n_cols_grid} grid)")
    print(f"{'='*60}")

    R_raw = load_red_channel(name)
    H, W = R_raw.shape

    # Rotation correction
    angle = detect_rotation_angle(R_raw)
    if abs(angle) > 0.01:
        R = ndimage.rotate(R_raw, angle, reshape=False, order=1, mode='nearest')
    else:
        R = R_raw
    print(f"  Rotation: {angle:+.2f} deg")

    # Detect all spots
    spots, row_ys, counts = detect_grid(R)
    n_spots = len(spots)
    print(f"  Total spots: {n_spots}")

    # Measure intensities
    meas = np.array([measure_intensity(R, cy, cx)
                     for _, cy, cx in spots]) if spots else np.zeros((0, 3))

    # Get spot coordinates
    coords = np.array([(cy, cx) for _, cy, cx in spots])

    # Define sub-region boundaries
    y_min, y_max = coords[:, 0].min(), coords[:, 0].max()
    x_min, x_max = coords[:, 1].min(), coords[:, 1].max()

    y_edges = np.linspace(y_min - 1, y_max + 1, n_rows_grid + 1)
    x_edges = np.linspace(x_min - 1, x_max + 1, n_cols_grid + 1)

    # Assign each spot to a sub-region
    region_ids = np.zeros(n_spots, dtype=int)
    for i in range(n_spots):
        ry = np.clip(np.searchsorted(y_edges, coords[i, 0]) - 1, 0, n_rows_grid - 1)
        rx = np.clip(np.searchsorted(x_edges, coords[i, 1]) - 1, 0, n_cols_grid - 1)
        region_ids[i] = ry * n_cols_grid + rx

    # Compute distribution per region using SAME bin edges
    # Use percentile-based bins from the full image
    all_net = meas[:, 2]
    p5 = np.percentile(all_net, 5)
    p10 = np.percentile(all_net, 10)
    quench = max(abs(p5), abs(p10), 5.0)
    bright = all_net[all_net > quench]
    p90 = np.percentile(bright, 90) if len(bright) >= 20 else quench + 50
    r = p90 - quench
    bin_edges = np.array([-np.inf, quench, quench + r * 0.22,
                          quench + r * 0.50, quench + r * 0.80, np.inf])
    n_levels = 5
    level_names = ["quenched", "dim", "medium", "bright", "v.bright"]

    print(f"  Bin edges: {[f'{e:.1f}' for e in bin_edges[1:-1]]}")

    # Per-region distributions
    region_data = []
    for rid in range(n_rows_grid * n_cols_grid):
        mask = region_ids == rid
        n_in = mask.sum()
        net = meas[mask, 2]

        dist = np.zeros(n_levels)
        for k in range(n_levels):
            if k == 0:
                dist[k] = (net <= bin_edges[1]).sum()
            elif k == n_levels - 1:
                dist[k] = (net > bin_edges[-2]).sum()
            else:
                dist[k] = ((net > bin_edges[k]) & (net <= bin_edges[k + 1])).sum()

        frac = dist / max(dist.sum(), 1)

        ry = rid // n_cols_grid
        rx = rid % n_cols_grid
        label = f"R{ry}C{rx}"

        region_data.append({
            'label': label,
            'n_spots': n_in,
            'counts': dist.astype(int),
            'frac': frac,
            'mean_net': net.mean() if n_in > 0 else 0,
            'std_net': net.std() if n_in > 0 else 0,
            'median_net': np.median(net) if n_in > 0 else 0,
        })

        parts = [f"{level_names[k]}={frac[k]:.3f}" for k in range(n_levels)]
        print(f"  {label} ({n_in:5d} spots): {', '.join(parts)}  "
              f"mean_net={net.mean():.1f}")

    # Summary statistics across regions
    fracs = np.array([rd['frac'] for rd in region_data])
    print(f"\n  Distribution consistency (std across regions):")
    for k in range(n_levels):
        vals = fracs[:, k]
        print(f"    {level_names[k]:10s}: mean={vals.mean():.4f}, "
              f"std={vals.std():.4f}, "
              f"CV={vals.std()/max(vals.mean(), 1e-10):.3f}, "
              f"range=[{vals.min():.4f}, {vals.max():.4f}]")

    return region_data, R_raw, R, spots, meas, coords, region_ids, \
           y_edges, x_edges, bin_edges, angle


def plot_results(name, region_data, R_raw, R, spots, meas, coords,
                 region_ids, y_edges, x_edges, bin_edges, angle,
                 n_rows_grid, n_cols_grid):
    """Generate visualization."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_levels = 5
    n_regions = n_rows_grid * n_cols_grid
    level_names = ["Quenched", "Dim", "Medium", "Bright", "V.Bright"]
    bar_colors = ["#4444ff", "#00c8c8", "#00ff00", "#ffc800", "#ff3232"]

    # --- Figure 1: Bar chart grid ---
    fig, axes = plt.subplots(n_rows_grid, n_cols_grid,
                             figsize=(4 * n_cols_grid, 3.5 * n_rows_grid),
                             sharey=True)
    if n_rows_grid == 1:
        axes = axes[np.newaxis, :]
    if n_cols_grid == 1:
        axes = axes[:, np.newaxis]

    for rid in range(n_regions):
        ry = rid // n_cols_grid
        rx = rid % n_cols_grid
        ax = axes[ry, rx]
        rd = region_data[rid]
        ax.bar(range(n_levels), rd['frac'], color=bar_colors)
        ax.set_xticks(range(n_levels))
        ax.set_xticklabels(level_names, fontsize=7, rotation=30)
        ax.set_title(f"{rd['label']} (n={rd['n_spots']})", fontsize=10)
        ax.set_ylim(0, 0.55)
        for k, v in enumerate(rd['frac']):
            ax.text(k, v + 0.01, f"{v:.3f}", ha="center", fontsize=7)

    fig.suptitle(f"{name}: Subregion Intensity Distributions "
                 f"({n_rows_grid}x{n_cols_grid})", fontsize=14)
    fig.tight_layout()
    outpath = os.path.join(DIR, f"_subregion_dist_{name}.png")
    fig.savefig(outpath, dpi=150)
    plt.close()
    print(f"  Saved {outpath}")

    # --- Figure 2: Overlay on image showing sub-region boundaries ---
    H, W = R.shape
    vmin, vmax = R.min(), R.max()
    norm = ((R - vmin) / max(vmax - vmin, 1) * 220).astype(np.uint8)
    rgb = np.stack([norm, norm // 4, norm // 6], axis=2)
    img = Image.fromarray(rgb)
    draw = ImageDraw.Draw(img)

    # Draw sub-region grid lines
    for ye in y_edges[1:-1]:
        draw.line([(0, int(ye)), (W - 1, int(ye))], fill=(255, 255, 0), width=2)
    for xe in x_edges[1:-1]:
        draw.line([(int(xe), 0), (int(xe), H - 1)], fill=(255, 255, 0), width=2)

    # Label each region with spot count
    for rid in range(n_regions):
        ry = rid // n_cols_grid
        rx = rid % n_cols_grid
        cy = (y_edges[ry] + y_edges[ry + 1]) / 2
        cx = (x_edges[rx] + x_edges[rx + 1]) / 2
        rd = region_data[rid]
        text = f"{rd['label']}\nn={rd['n_spots']}"
        draw.text((int(cx) - 30, int(cy) - 15), text, fill=(255, 255, 0))

    outpath2 = os.path.join(DIR, f"_subregion_grid_{name}.png")
    img.save(outpath2)
    print(f"  Saved {outpath2}")

    # --- Figure 3: Heatmap of quenched fraction per region ---
    fig, axes = plt.subplots(1, n_levels, figsize=(3.5 * n_levels, 3))
    for k in range(n_levels):
        ax = axes[k]
        grid = np.zeros((n_rows_grid, n_cols_grid))
        for rid in range(n_regions):
            ry = rid // n_cols_grid
            rx = rid % n_cols_grid
            grid[ry, rx] = region_data[rid]['frac'][k]
        im = ax.imshow(grid, cmap='YlOrRd', vmin=0,
                       vmax=max(grid.max(), 0.01), aspect='auto')
        ax.set_title(level_names[k], fontsize=10)
        for iy in range(n_rows_grid):
            for ix in range(n_cols_grid):
                ax.text(ix, iy, f"{grid[iy, ix]:.3f}",
                        ha='center', va='center', fontsize=8)
        ax.set_xticks(range(n_cols_grid))
        ax.set_yticks(range(n_rows_grid))
        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle(f"{name}: Fraction heatmap per subregion", fontsize=13)
    fig.tight_layout()
    outpath3 = os.path.join(DIR, f"_subregion_heatmap_{name}.png")
    fig.savefig(outpath3, dpi=150)
    plt.close()
    print(f"  Saved {outpath3}")


def main():
    n_rows_grid = 2
    n_cols_grid = 4

    for name in ["0ugL", "20ugL"]:
        result = analyze_subregions(name, n_rows_grid, n_cols_grid)
        region_data = result[0]
        try:
            plot_results(name, *result, n_rows_grid, n_cols_grid)
        except ImportError:
            print("  (matplotlib not available, skipping plots)")


if __name__ == "__main__":
    main()
