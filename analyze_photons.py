#!/usr/bin/env python3
"""
NEXT Honeycomb Cell — Photon Analysis
=======================================
Reads ROOT ntuple from ELHolePhotonSim and computes:
  - Photon acceptance: fraction reaching SiPM
  - Reflection statistics
  - Loss channels: wall absorption, top escape
  - Per-electron detection efficiency

Usage:
    python3 analyze_photons.py output.root [--plot photon_analysis.png]
"""

import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def load_root_data(fname):
    """Load photon ntuple using uproot."""
    try:
        import uproot
        f = uproot.open(fname)
        tree = f["photons"]
        return tree.arrays(library="np")
    except ImportError:
        raise ImportError("Install uproot: pip install uproot awkward")


def analyze_and_plot(data, outpng):
    electron_id   = data["electron_id"]
    n_reflections = data["n_reflections"]
    reached_sipm  = data["reached_sipm"].astype(bool)
    absorbed_wall = data["absorbed_wall"].astype(bool)
    escaped_top   = data.get("escaped_top", np.zeros(len(electron_id))).astype(bool)
    src_x = data["src_x_cm"]
    src_y = data["src_y_cm"]
    src_z = data["src_z_cm"]
    hit_x = data["hit_x_cm"]
    hit_y = data["hit_y_cm"]

    total     = len(electron_id)
    n_sipm    = reached_sipm.sum()
    n_wall    = absorbed_wall.sum()
    n_escaped = escaped_top.sum()
    n_other   = total - n_sipm - n_wall - n_escaped
    eff       = n_sipm / total * 100 if total > 0 else 0

    print(f"\n=== Photon Statistics ===")
    print(f"  Total photons      : {total:,}")
    print(f"  Reached SiPM       : {n_sipm:,} ({eff:.1f}%)")
    print(f"  Absorbed by wall   : {n_wall:,} ({n_wall/total*100:.1f}%)")
    print(f"  Escaped top        : {n_escaped:,} ({n_escaped/total*100:.1f}%)")
    print(f"  Other/lost         : {n_other:,} ({n_other/total*100:.1f}%)")
    mean_r = n_reflections[reached_sipm].mean() if n_sipm > 0 else 0
    if n_sipm > 0:
        print(f"  Mean reflections (SiPM hits): {mean_r:.2f}")

    # ── Per-electron statistics ──
    eids = np.unique(electron_id)
    e_effs = []
    e_total_ph = []
    e_detected = []
    for eid in eids:
        mask_e = electron_id == eid
        n_ph = mask_e.sum()
        n_det = (mask_e & reached_sipm).sum()
        e_effs.append(n_det / n_ph * 100 if n_ph > 0 else 0)
        e_total_ph.append(n_ph)
        e_detected.append(n_det)

    print(f"\n=== Per-electron (N={len(eids)}) ===")
    print(f"  Mean photons/electron : {np.mean(e_total_ph):.1f}")
    print(f"  Mean detected/electron: {np.mean(e_detected):.1f}")
    print(f"  Mean efficiency       : {np.mean(e_effs):.1f}%")

    # ── Figure ──────────────────────────────────────────────────
    DARK  = "#0a0e1a"; PANEL = "#141824"; CYAN = "#00e5ff"
    AMBER = "#ffb300"; MINT  = "#69ffb4"; ROSE = "#ff5580"; WHITE = "#e8eaf6"

    fig = plt.figure(figsize=(18, 12), facecolor=DARK)
    gs  = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.36,
                   left=0.07, right=0.97, top=0.93, bottom=0.07)

    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)]
    for ax in axes:
        ax.set_facecolor(PANEL)
        for sp in ax.spines.values(): sp.set_edgecolor("#2a3050")
        ax.tick_params(colors=WHITE, labelsize=8)
        ax.xaxis.label.set_color(WHITE); ax.yaxis.label.set_color(WHITE)
        ax.title.set_color(CYAN)

    ax_hit, ax_ref, ax_src, ax_loss, ax_refvz, ax_eff = axes

    # 1. Hit position map on SiPM
    if n_sipm > 0:
        ax_hit.hist2d(hit_x[reached_sipm], hit_y[reached_sipm],
                      bins=40, cmap="plasma")
        ax_hit.set_xlabel("x at SiPM [cm]"); ax_hit.set_ylabel("y at SiPM [cm]")
        ax_hit.set_title("Photon Arrival on SiPM", fontsize=9, fontweight="bold")
        ax_hit.set_aspect("equal")

    # 2. Reflection count distribution
    if n_sipm > 0:
        max_ref = int(n_reflections[reached_sipm].max())
        bins_r = np.arange(0, min(max_ref + 2, 40))
        ax_ref.hist(n_reflections[reached_sipm], bins=bins_r,
                    color=CYAN, alpha=0.8, edgecolor=DARK)
        ax_ref.axvline(mean_r, color=AMBER, lw=1.5, ls="--",
                       label=f"mean={mean_r:.1f}")
        ax_ref.legend(fontsize=8, framealpha=0.3, labelcolor=WHITE, facecolor=PANEL)
    ax_ref.set_xlabel("Wall reflections before SiPM")
    ax_ref.set_ylabel("Photon count")
    ax_ref.set_title("Reflections Distribution", fontsize=9, fontweight="bold")

    # 3. Source positions (x-z)
    ax_src.hist2d(src_x, src_z, bins=40, cmap="viridis")
    ax_src.set_xlabel("src x [cm]"); ax_src.set_ylabel("src z [cm] (0=SiPM)")
    ax_src.set_title("Photon Source Positions", fontsize=9, fontweight="bold")

    # 4. Loss channel pie/bar chart
    labels = ['SiPM', 'Wall abs.', 'Escaped top', 'Other']
    counts = [n_sipm, n_wall, n_escaped, n_other]
    colors_pie = [MINT, ROSE, AMBER, '#555555']
    bars = ax_loss.barh(labels, counts, color=colors_pie, edgecolor=DARK)
    for bar, count in zip(bars, counts):
        ax_loss.text(bar.get_width() + total * 0.01, bar.get_y() + bar.get_height() / 2,
                     f'{count:,} ({count/total*100:.1f}%)',
                     va='center', color=WHITE, fontsize=8)
    ax_loss.set_xlabel("Photon count")
    ax_loss.set_title("Photon Fate", fontsize=9, fontweight="bold")

    # 5. Mean reflections vs source z
    if n_sipm > 0:
        z_bins = np.linspace(src_z.min(), src_z.max(), 16)
        z_mid = 0.5 * (z_bins[:-1] + z_bins[1:])
        r_mean, r_std = [], []
        for lo, hi in zip(z_bins[:-1], z_bins[1:]):
            m = reached_sipm & (src_z >= lo) & (src_z < hi)
            if m.sum() > 0:
                r_mean.append(n_reflections[m].mean())
                r_std.append(n_reflections[m].std())
            else:
                r_mean.append(np.nan); r_std.append(0)
        r_mean = np.array(r_mean); r_std = np.array(r_std)
        ax_refvz.plot(z_mid, r_mean, color=MINT, lw=2)
        ax_refvz.fill_between(z_mid, r_mean - r_std, r_mean + r_std,
                              color=MINT, alpha=0.25)
    ax_refvz.set_xlabel("Source z [cm] (0=SiPM)")
    ax_refvz.set_ylabel("Mean reflections")
    ax_refvz.set_title("Reflections vs Source Depth", fontsize=9, fontweight="bold")

    # 6. Per-electron detection efficiency
    ax_eff.hist(e_effs, bins=20, color=AMBER, alpha=0.8, edgecolor=DARK)
    ax_eff.axvline(np.mean(e_effs), color=ROSE, lw=1.5, ls="--",
                   label=f"mean={np.mean(e_effs):.1f}%")
    ax_eff.legend(fontsize=8, framealpha=0.3, labelcolor=WHITE, facecolor=PANEL)
    ax_eff.set_xlabel("Photon detection efficiency [%]")
    ax_eff.set_ylabel("Electrons")
    ax_eff.set_title("Per-Electron Acceptance", fontsize=9, fontweight="bold")

    fig.suptitle(
        f"NEXT Honeycomb Hole — Photon Simulation  ·  {total:,} photons  ·  "
        f"SiPM acceptance {eff:.1f}%  ·  Mean reflections {mean_r:.1f}",
        color=WHITE, fontsize=12, fontweight="bold", y=0.97)

    fig.savefig(outpng, dpi=150, facecolor=DARK, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved → {outpng}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile")
    parser.add_argument("--plot", default="photon_analysis.png")
    args = parser.parse_args()
    data = load_root_data(args.rootfile)
    analyze_and_plot(data, args.plot)
