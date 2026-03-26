#!/usr/bin/env python3
"""
Spot morphology analysis — understand the "adhesion" (粘连) phenomenon.

For each detected spot, measure:
- Gaussian-fit width (sigma_x, sigma_y)
- Peak intensity
- Whether the spot is elongated or round
- Whether there are sub-peaks within the spot region

Goal: characterize single-branch vs multi-branch (adhesion) spots.
"""

import numpy as np
from PIL import Image
from scipy import ndimage, optimize
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_grid import (
    load_red_channel, detect_rotation_angle, detect_grid,
    measure_intensity
)

DIR = os.path.dirname(os.path.abspath(__file__))


def gaussian_2d(coords, amp, x0, y0, sigma_x, sigma_y, bg):
    """2D Gaussian model."""
    y, x = coords
    return bg + amp * np.exp(-((x - x0)**2 / (2 * sigma_x**2) +
                               (y - y0)**2 / (2 * sigma_y**2)))


def fit_gaussian(R, cy, cx, radius=8):
    """Fit a 2D Gaussian to a spot. Returns (sigma_x, sigma_y, amplitude, residual)."""
    H, W = R.shape
    y0 = max(0, int(cy) - radius)
    y1 = min(H, int(cy) + radius + 1)
    x0 = max(0, int(cx) - radius)
    x1 = min(W, int(cx) + radius + 1)
    patch = R[y0:y1, x0:x1]

    yy, xx = np.mgrid[0:patch.shape[0], 0:patch.shape[1]]
    coords = (yy.ravel(), xx.ravel())
    data = patch.ravel()

    bg_est = np.percentile(patch, 10)
    amp_est = patch.max() - bg_est
    cy_local = cy - y0
    cx_local = cx - x0

    try:
        popt, _ = optimize.curve_fit(
            gaussian_2d, coords, data,
            p0=[amp_est, cx_local, cy_local, 2.0, 2.0, bg_est],
            bounds=([0, -1, -1, 0.5, 0.5, 0],
                    [500, patch.shape[1]+1, patch.shape[0]+1, 15, 15, 200]),
            maxfev=2000
        )
        amp, xc, yc, sx, sy, bg = popt
        fitted = gaussian_2d(coords, *popt)
        residual = np.sqrt(np.mean((data - fitted)**2))
        return sx, sy, amp, residual, bg
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, np.nan


def count_local_peaks(R, cy, cx, radius=10, min_prominence=3.0):
    """Count distinct peaks within a radius of the spot center."""
    H, W = R.shape
    y0 = max(0, int(cy) - radius)
    y1 = min(H, int(cy) + radius + 1)
    x0 = max(0, int(cx) - radius)
    x1 = min(W, int(cx) + radius + 1)
    patch = R[y0:y1, x0:x1]

    # Use local maxima detection
    max_filt = ndimage.maximum_filter(patch, size=3)
    local_max = (patch == max_filt)

    # Background level
    bg = np.percentile(patch, 20)
    prominent = local_max & (patch > bg + min_prominence)

    return prominent.sum()


def measure_spot_extent(R, cy, cx, threshold_frac=0.3, max_radius=15):
    """Measure the extent of the bright region around a spot.
    Returns the number of pixels above threshold and the effective radius."""
    H, W = R.shape
    y0 = max(0, int(cy) - max_radius)
    y1 = min(H, int(cy) + max_radius + 1)
    x0 = max(0, int(cx) - max_radius)
    x1 = min(W, int(cx) + max_radius + 1)
    patch = R[y0:y1, x0:x1]

    bg = np.percentile(patch, 10)
    peak = patch.max()
    threshold = bg + (peak - bg) * threshold_frac

    above = patch > threshold
    n_pixels = above.sum()
    # Effective radius: sqrt(area / pi)
    r_eff = np.sqrt(n_pixels / np.pi) if n_pixels > 0 else 0

    return n_pixels, r_eff


def main():
    print("=" * 60)
    print("Spot Morphology Analysis")
    print("=" * 60)

    R_raw = load_red_channel("0ugL")
    angle = detect_rotation_angle(R_raw)
    R = ndimage.rotate(R_raw, angle, reshape=False, order=1,
                       mode='nearest') if abs(angle) > 0.01 else R_raw
    print(f"  Rotation: {angle:+.2f} deg")

    spots, row_ys, counts = detect_grid(R)
    n_spots = len(spots)
    print(f"  Total spots: {n_spots}")

    # Analyze morphology — sample 2000 spots for speed
    np.random.seed(42)
    sample_n = min(2000, n_spots)
    sample_idx = np.random.choice(n_spots, sample_n, replace=False)
    print(f"\n  Fitting Gaussians on {sample_n}/{n_spots} sampled spots...")
    results = []
    for count, i in enumerate(sample_idx):
        if count % 500 == 0:
            print(f"    {count}/{sample_n}...")
        ri, cy, cx = spots[i]

        sx, sy, amp, resid, bg = fit_gaussian(R, cy, cx)
        n_peaks = count_local_peaks(R, cy, cx)
        n_pix, r_eff = measure_spot_extent(R, cy, cx)
        spot_val, bg_val, net = measure_intensity(R, cy, cx)

        results.append({
            'ri': ri, 'cy': cy, 'cx': cx,
            'sigma_x': sx, 'sigma_y': sy,
            'amplitude': amp, 'fit_residual': resid, 'fit_bg': bg,
            'n_peaks': n_peaks,
            'n_pix_above': n_pix, 'r_eff': r_eff,
            'net_signal': net,
            'elongation': max(sx, sy) / min(sx, sy) if (not np.isnan(sx) and min(sx, sy) > 0.1) else np.nan,
            'sigma_mean': (sx + sy) / 2 if not np.isnan(sx) else np.nan,
        })

    results = [r for r in results if not np.isnan(r['sigma_x'])]
    print(f"  Successfully fitted: {len(results)}/{n_spots}")

    # Statistics
    sigmas = np.array([r['sigma_mean'] for r in results])
    elongs = np.array([r['elongation'] for r in results])
    n_peaks_arr = np.array([r['n_peaks'] for r in results])
    r_effs = np.array([r['r_eff'] for r in results])
    amps = np.array([r['amplitude'] for r in results])
    nets = np.array([r['net_signal'] for r in results])

    print(f"\n  Gaussian sigma (mean of x,y):")
    print(f"    mean={sigmas.mean():.2f}, std={sigmas.std():.2f}, "
          f"median={np.median(sigmas):.2f}")
    for p in [10, 25, 50, 75, 90, 95, 99]:
        print(f"    P{p}={np.percentile(sigmas, p):.2f}")

    print(f"\n  Elongation (max_sigma / min_sigma):")
    print(f"    mean={elongs.mean():.2f}, median={np.median(elongs):.2f}")
    for p in [50, 75, 90, 95, 99]:
        print(f"    P{p}={np.percentile(elongs, p):.2f}")

    print(f"\n  Local peaks within r=10:")
    for n in range(1, 8):
        frac = (n_peaks_arr == n).sum() / len(n_peaks_arr)
        print(f"    {n} peaks: {(n_peaks_arr == n).sum()} ({frac:.3f})")
    print(f"    >=8 peaks: {(n_peaks_arr >= 8).sum()}")

    print(f"\n  Effective radius (pixels at 30% threshold):")
    print(f"    mean={r_effs.mean():.2f}, std={r_effs.std():.2f}, "
          f"median={np.median(r_effs):.2f}")
    for p in [50, 75, 90, 95, 99]:
        print(f"    P{p}={np.percentile(r_effs, p):.2f}")

    # Correlation: bright spots tend to be wider?
    print(f"\n  Correlation analysis:")
    bright = nets > np.percentile(nets, 75)
    dim = nets < np.percentile(nets, 25)
    print(f"    Bright spots (top 25%): sigma={sigmas[bright].mean():.2f}, "
          f"r_eff={r_effs[bright].mean():.2f}, "
          f"elongation={elongs[bright].mean():.2f}, "
          f"n_peaks={n_peaks_arr[bright].mean():.1f}")
    print(f"    Dim spots (bottom 25%): sigma={sigmas[dim].mean():.2f}, "
          f"r_eff={r_effs[dim].mean():.2f}, "
          f"elongation={elongs[dim].mean():.2f}, "
          f"n_peaks={n_peaks_arr[dim].mean():.1f}")

    # Identify candidate "adhesion" spots
    # Criteria: large sigma OR high elongation OR many sub-peaks
    adhesion_sigma = sigmas > np.percentile(sigmas, 90)
    adhesion_elong = elongs > np.percentile(elongs, 90)
    adhesion_peaks = n_peaks_arr >= 4
    adhesion_reff = r_effs > np.percentile(r_effs, 90)
    adhesion_any = adhesion_sigma | adhesion_elong | adhesion_peaks | adhesion_reff

    print(f"\n  Candidate adhesion spots:")
    print(f"    By sigma > P90 ({np.percentile(sigmas, 90):.2f}): "
          f"{adhesion_sigma.sum()}")
    print(f"    By elongation > P90 ({np.percentile(elongs, 90):.2f}): "
          f"{adhesion_elong.sum()}")
    print(f"    By n_peaks >= 4: {adhesion_peaks.sum()}")
    print(f"    By r_eff > P90 ({np.percentile(r_effs, 90):.2f}): "
          f"{adhesion_reff.sum()}")
    print(f"    Any criterion: {adhesion_any.sum()} "
          f"({adhesion_any.sum()/len(results)*100:.1f}%)")

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 3, figsize=(15, 9))

        # 1. Sigma distribution
        ax = axes[0, 0]
        ax.hist(sigmas, bins=80, range=(0, 10), color='steelblue', alpha=0.8)
        ax.axvline(np.percentile(sigmas, 90), color='red', ls='--',
                   label=f'P90={np.percentile(sigmas, 90):.1f}')
        ax.set_xlabel('Gaussian sigma (px)')
        ax.set_ylabel('Count')
        ax.set_title('Spot size distribution')
        ax.legend()

        # 2. Elongation distribution
        ax = axes[0, 1]
        ax.hist(elongs, bins=80, range=(1, 5), color='darkorange', alpha=0.8)
        ax.axvline(np.percentile(elongs, 90), color='red', ls='--',
                   label=f'P90={np.percentile(elongs, 90):.1f}')
        ax.set_xlabel('Elongation (max_σ / min_σ)')
        ax.set_title('Spot elongation distribution')
        ax.legend()

        # 3. Number of local peaks
        ax = axes[0, 2]
        peak_vals = list(range(1, 8))
        peak_counts = [(n_peaks_arr == n).sum() for n in peak_vals]
        peak_counts.append((n_peaks_arr >= 8).sum())
        ax.bar(peak_vals + [8], peak_counts, color='seagreen', alpha=0.8)
        ax.set_xlabel('Number of local peaks')
        ax.set_ylabel('Count')
        ax.set_xticks(peak_vals + [8])
        ax.set_xticklabels([str(n) for n in peak_vals] + ['≥8'])
        ax.set_title('Sub-peaks per spot')

        # 4. Effective radius distribution
        ax = axes[1, 0]
        ax.hist(r_effs, bins=80, range=(0, 10), color='mediumpurple', alpha=0.8)
        ax.axvline(np.percentile(r_effs, 90), color='red', ls='--',
                   label=f'P90={np.percentile(r_effs, 90):.1f}')
        ax.set_xlabel('Effective radius (px)')
        ax.set_title('Spot extent distribution')
        ax.legend()

        # 5. Sigma vs net signal scatter
        ax = axes[1, 1]
        idx = np.random.choice(len(results), min(3000, len(results)), replace=False)
        ax.scatter(nets[idx], sigmas[idx], s=2, alpha=0.3, c='steelblue')
        ax.set_xlabel('Net signal')
        ax.set_ylabel('Gaussian sigma (px)')
        ax.set_title('Size vs brightness')

        # 6. Elongation vs net signal scatter
        ax = axes[1, 2]
        ax.scatter(nets[idx], elongs[idx], s=2, alpha=0.3, c='darkorange')
        ax.set_xlabel('Net signal')
        ax.set_ylabel('Elongation')
        ax.set_title('Elongation vs brightness')

        fig.suptitle('0ugL Spot Morphology Analysis', fontsize=14)
        fig.tight_layout()
        outpath = os.path.join(DIR, '_spot_morphology.png')
        fig.savefig(outpath, dpi=150)
        plt.close()
        print(f"\n  Saved {outpath}")

        # Crop examples of adhesion spots
        fig, axes = plt.subplots(4, 8, figsize=(16, 8))
        fig.suptitle('Top 32 spots by Gaussian sigma (candidate adhesion)', fontsize=13)

        sorted_idx = np.argsort(sigmas)[::-1]
        for k in range(32):
            ax = axes[k // 8, k % 8]
            r = results[sorted_idx[k]]
            cy, cx = int(r['cy']), int(r['cx'])
            rad = 12
            y0 = max(0, cy - rad)
            y1 = min(R.shape[0], cy + rad + 1)
            x0 = max(0, cx - rad)
            x1 = min(R.shape[1], cx + rad + 1)
            patch = R[y0:y1, x0:x1]
            ax.imshow(patch, cmap='hot', vmin=0,
                      vmax=np.percentile(R, 99.5))
            ax.set_title(f'σ={r["sigma_mean"]:.1f}\nn={r["n_peaks"]}',
                         fontsize=7)
            ax.axis('off')

        fig.tight_layout()
        outpath2 = os.path.join(DIR, '_adhesion_examples.png')
        fig.savefig(outpath2, dpi=150)
        plt.close()
        print(f"  Saved {outpath2}")

        # Also show examples of "normal" small spots for comparison
        fig, axes = plt.subplots(4, 8, figsize=(16, 8))
        fig.suptitle('32 typical small spots (sigma near median)', fontsize=13)

        median_sigma = np.median(sigmas)
        near_median = np.abs(sigmas - median_sigma)
        median_idx = np.argsort(near_median)
        for k in range(32):
            ax = axes[k // 8, k % 8]
            r = results[median_idx[k]]
            cy, cx = int(r['cy']), int(r['cx'])
            rad = 12
            y0 = max(0, cy - rad)
            y1 = min(R.shape[0], cy + rad + 1)
            x0 = max(0, cx - rad)
            x1 = min(R.shape[1], cx + rad + 1)
            patch = R[y0:y1, x0:x1]
            ax.imshow(patch, cmap='hot', vmin=0,
                      vmax=np.percentile(R, 99.5))
            ax.set_title(f'σ={r["sigma_mean"]:.1f}\nn={r["n_peaks"]}',
                         fontsize=7)
            ax.axis('off')

        fig.tight_layout()
        outpath3 = os.path.join(DIR, '_normal_spot_examples.png')
        fig.savefig(outpath3, dpi=150)
        plt.close()
        print(f"  Saved {outpath3}")

    except ImportError:
        print("  (matplotlib not available)")


if __name__ == "__main__":
    main()
