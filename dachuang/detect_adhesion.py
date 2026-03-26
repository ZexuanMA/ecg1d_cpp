#!/usr/bin/env python3
"""
Adhesion (粘连) detection and microsphere-level intensity measurement.

Physical model:
- Each microsphere has ~5 branches that can independently emit light
- When multiple branches emit, their combined PSF forms a large blob
- This blob can span multiple well positions

Algorithm:
1. Detect wells as before (peak detection + gap filling)
2. Binary threshold the smoothed image to find connected bright regions
3. For each connected region, count how many well positions fall inside
4. Multi-well regions = adhesion → merge into one microsphere
"""

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_grid import (
    load_red_channel, detect_rotation_angle, detect_grid,
    measure_intensity
)

DIR = os.path.dirname(os.path.abspath(__file__))


def find_adhesion_clusters(R, spots, min_prominence=5.0):
    """Find connected bright regions and map wells to them."""
    H, W = R.shape

    R_smooth = ndimage.gaussian_filter(R, sigma=1.5)
    bg = ndimage.uniform_filter(R_smooth, size=50)
    above = R_smooth > (bg + min_prominence)

    struct = ndimage.generate_binary_structure(2, 2)  # 8-connectivity
    labels, n_components = ndimage.label(above, structure=struct)

    # Map each well to its component
    n_spots = len(spots)
    spot_comp = np.zeros(n_spots, dtype=int)
    for i, (ri, cy, cx) in enumerate(spots):
        iy, ix = int(round(cy)), int(round(cx))
        if 0 <= iy < H and 0 <= ix < W:
            spot_comp[i] = labels[iy, ix]

    # Group by component → clusters
    cluster_labels = np.zeros(n_spots, dtype=int)
    comp_to_cluster = {}
    next_id = 1
    for i in range(n_spots):
        c = spot_comp[i]
        if c == 0:
            cluster_labels[i] = next_id
            next_id += 1
        else:
            if c not in comp_to_cluster:
                comp_to_cluster[c] = next_id
                next_id += 1
            cluster_labels[i] = comp_to_cluster[c]

    return cluster_labels, labels, spot_comp


def main():
    print("=" * 60)
    print("Adhesion Detection")
    print("=" * 60)

    R_raw = load_red_channel("0ugL")
    angle = detect_rotation_angle(R_raw)
    R = ndimage.rotate(R_raw, angle, reshape=False, order=1,
                       mode='nearest') if abs(angle) > 0.01 else R_raw
    H, W = R.shape
    print(f"  Image: {H}x{W}, rotation: {angle:+.2f} deg")

    spots, row_ys, counts = detect_grid(R)
    n_spots = len(spots)
    print(f"  Detected wells: {n_spots}")

    prom = 5.0
    print(f"  Prominence threshold: {prom}")

    cluster_labels, comp_map, spot_comp = \
        find_adhesion_clusters(R, spots, min_prominence=prom)

    # Cluster size statistics
    unique_clusters, cluster_counts = np.unique(cluster_labels, return_counts=True)
    n_clusters = len(unique_clusters)

    size_hist = {}
    for sz in cluster_counts:
        size_hist[sz] = size_hist.get(sz, 0) + 1

    n_single = (cluster_counts == 1).sum()
    n_adhesion = (cluster_counts > 1).sum()
    wells_in_adhesion = sum(c for c in cluster_counts if c > 1)
    n_quenched = (spot_comp == 0).sum()

    print(f"\n  Results:")
    print(f"    Total microspheres: {n_clusters}")
    print(f"    Single-well: {n_single}")
    print(f"    Adhesion (>1 well): {n_adhesion}")
    print(f"    Wells in adhesion: {wells_in_adhesion}/{n_spots} "
          f"({wells_in_adhesion/n_spots*100:.1f}%)")
    print(f"    Quenched (below threshold): {n_quenched}")

    print(f"\n  Cluster size distribution:")
    for sz in sorted(size_hist.keys())[:15]:
        print(f"    size={sz}: {size_hist[sz]} clusters")
    if max(size_hist.keys()) > 15:
        big = {k: v for k, v in size_hist.items() if k > 15}
        for sz in sorted(big.keys()):
            print(f"    size={sz}: {big[sz]} clusters")

    # Intensity comparison: single vs adhesion
    print(f"\n  Intensity comparison:")
    single_idx = []
    adhesion_clusters = {}  # cluster_id -> member indices

    for cid, cnt in zip(unique_clusters, cluster_counts):
        members = np.where(cluster_labels == cid)[0]
        if cnt == 1:
            single_idx.append(members[0])
        else:
            adhesion_clusters[cid] = members

    # Measure single wells
    single_nets = []
    for idx in single_idx[:500]:  # sample for speed
        _, cy, cx = spots[idx]
        _, _, net = measure_intensity(R, cy, cx)
        single_nets.append(net)
    single_nets = np.array(single_nets)
    print(f"    Single wells: mean={single_nets.mean():.1f}, "
          f"median={np.median(single_nets):.1f}")

    # Measure adhesion: integrated flux over component
    bg_map = ndimage.uniform_filter(R, size=50)
    adhesion_results = []
    count = 0
    for cid, members in adhesion_clusters.items():
        if count >= 200:
            break
        comp_id = spot_comp[members[0]]
        if comp_id == 0:
            continue

        comp_mask = (comp_map == comp_id)
        comp_pixels = R[comp_mask]

        # Background from border
        dilated = ndimage.binary_dilation(comp_mask, iterations=5)
        border = dilated & ~comp_mask
        bg_val = np.median(R[border]) if border.sum() > 10 else bg_map[
            int(spots[members[0]][1]), int(spots[members[0]][2])]

        total_flux = (comp_pixels - bg_val).clip(0).sum()
        n_members = len(members)

        adhesion_results.append({
            'n_members': n_members,
            'total_flux': total_flux,
            'per_well': total_flux / n_members,
            'comp_area': comp_mask.sum(),
        })
        count += 1

    if adhesion_results:
        per_wells = np.array([r['per_well'] for r in adhesion_results])
        sizes = np.array([r['n_members'] for r in adhesion_results])
        print(f"    Adhesion per-well: mean={per_wells.mean():.1f}, "
              f"median={np.median(per_wells):.1f}")
        print(f"    Adhesion sizes: mean={sizes.mean():.1f}, "
              f"median={np.median(sizes):.0f}")

    # Visualization
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 3, figsize=(16, 10))

        # 1. Cluster size histogram
        ax = axes[0, 0]
        max_show = min(12, max(size_hist.keys()))
        bar_x = list(range(1, max_show + 1))
        bar_y = [size_hist.get(s, 0) for s in bar_x]
        bar_y.append(sum(v for k, v in size_hist.items() if k > max_show))
        ax.bar(bar_x + [max_show + 1], bar_y, color='steelblue', alpha=0.8)
        ax.set_xlabel('Wells per microsphere')
        ax.set_ylabel('Count')
        ax.set_title('Cluster size distribution')
        ax.set_xticks(bar_x + [max_show + 1])
        ax.set_xticklabels([str(s) for s in bar_x] + [f'>{max_show}'])

        # 2. Single vs adhesion intensity
        ax = axes[0, 1]
        bins = np.linspace(-10, 100, 60)
        ax.hist(single_nets, bins=bins, alpha=0.6, label='Single',
                color='steelblue', density=True)
        if adhesion_results:
            ax.hist(per_wells, bins=bins, alpha=0.6, label='Adhesion/well',
                    color='coral', density=True)
        ax.set_xlabel('Net signal per well')
        ax.set_title('Single vs adhesion intensity')
        ax.legend()

        # 3. Per-well intensity vs cluster size
        ax = axes[0, 2]
        if adhesion_results:
            ax.scatter(sizes, per_wells, s=10, alpha=0.4, c='coral')
            for sz in range(2, min(10, int(sizes.max()) + 1)):
                mask = sizes == sz
                if mask.sum() >= 3:
                    ax.plot(sz, per_wells[mask].mean(), 'ko', ms=8)
            ax.set_xlabel('Cluster size')
            ax.set_ylabel('Net signal per well')
            ax.set_title('Intensity vs cluster size (dots=mean)')

        # 4-6. Example adhesion crops
        # Sort adhesion by size, show top examples
        sorted_adhesion = sorted(adhesion_clusters.items(),
                                 key=lambda x: len(x[1]), reverse=True)

        for panel_idx, (cid, members) in enumerate(sorted_adhesion[:3]):
            ax = axes[1, panel_idx]
            cy_mean = np.mean([spots[m][1] for m in members])
            cx_mean = np.mean([spots[m][2] for m in members])
            rad = max(25, len(members) * 10)
            y0, y1 = max(0, int(cy_mean) - rad), min(H, int(cy_mean) + rad)
            x0, x1 = max(0, int(cx_mean) - rad), min(W, int(cx_mean) + rad)

            patch = R[y0:y1, x0:x1]
            ax.imshow(patch, cmap='hot', vmin=0,
                      vmax=np.percentile(R, 99.5))
            for m in members:
                my = spots[m][1] - y0
                mx = spots[m][2] - x0
                ax.plot(mx, my, 'c+', markersize=10, markeredgewidth=2)
            ax.set_title(f'Adhesion: {len(members)} wells merged')
            ax.axis('off')

        fig.suptitle(f'Adhesion Detection (0ugL, prominence={prom})',
                     fontsize=14)
        fig.tight_layout()
        outpath = os.path.join(DIR, '_adhesion_detection.png')
        fig.savefig(outpath, dpi=150)
        plt.close()
        print(f"\n  Saved {outpath}")

    except ImportError:
        print("  (matplotlib not available)")


if __name__ == "__main__":
    main()
