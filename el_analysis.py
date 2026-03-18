#!/usr/bin/env python3
"""
NEXT Honeycomb Cell — Electron Trajectory Analysis
====================================================
Reads el_electrons.json and produces analysis plots:
  1. Transverse plane map (initial positions, fate)
  2. Collection efficiency vs distance to hole centre
  3. Mean trajectories grouped by initial distance
  4. Final positions at anode/hole level
  5. Transverse displacement distribution
  6. Photon yield per collected electron

Usage:
    python3 el_analysis.py [el_electrons.json]
"""

import json
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec

# ─── Load data ─────────────────────────────────────────────────────
indir = os.path.dirname(os.path.abspath(__file__))
inpath = sys.argv[1] if len(sys.argv) > 1 else os.path.join(indir, "el_electrons.json")
with open(inpath) as f:
    data = json.load(f)

meta = data["metadata"]
electrons = data["electrons"]
hole_centres = np.array(meta["hole_centres"])
R_hole = meta["hole_radius_cm"]
pitch = meta["hole_pitch_cm"]

print(f"Loaded {len(electrons)} electrons, {len(hole_centres)} holes")
print(f"Geometric acceptance: {meta['geometric_acceptance_pct']:.1f}%")
print(f"Collection efficiency (→ SiPM): {meta['collection_efficiency_pct']:.1f}%")

# ─── Derived arrays ────────────────────────────────────────────────
x0 = np.array([e["x0"] for e in electrons])
y0 = np.array([e["y0"] for e in electrons])
xf = np.array([e["x_final"] for e in electrons])
yf = np.array([e["y_final"] for e in electrons])
entered = np.array([e["entered_hole"] for e in electrons], dtype=bool)
hit_wall = np.array([e["hit_wall"] for e in electrons], dtype=bool)
reached_sipm = np.array([e["reached_sipm"] for e in electrons], dtype=bool)
lost = ~entered  # didn't enter hole

# Initial distance to nearest hole
def dist_to_nearest(xs, ys, centres):
    dists = []
    for x, y in zip(xs, ys):
        d2 = (centres[:, 0] - x)**2 + (centres[:, 1] - y)**2
        dists.append(np.sqrt(d2.min()))
    return np.array(dists)

d_init = dist_to_nearest(x0, y0, hole_centres)

# ─── Figure layout ─────────────────────────────────────────────────
fig = plt.figure(figsize=(20, 16), facecolor="#0a0e1a")
gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38,
              left=0.07, right=0.97, top=0.93, bottom=0.06)

ax_map   = fig.add_subplot(gs[0, :2])
ax_eff   = fig.add_subplot(gs[0, 2])
ax_traj  = fig.add_subplot(gs[1, :])
ax_final = fig.add_subplot(gs[2, 0])
ax_disp  = fig.add_subplot(gs[2, 1])
ax_phot  = fig.add_subplot(gs[2, 2])

DARK  = "#0a0e1a"
PANEL = "#141824"
CYAN  = "#00e5ff"
AMBER = "#ffb300"
MINT  = "#69ffb4"
ROSE  = "#ff5580"
WHITE = "#e8eaf6"
BLUE  = "#4488ff"

for ax in [ax_map, ax_eff, ax_traj, ax_final, ax_disp, ax_phot]:
    ax.set_facecolor(PANEL)
    for spine in ax.spines.values():
        spine.set_edgecolor("#2a3050")
    ax.tick_params(colors=WHITE, labelsize=8)
    ax.xaxis.label.set_color(WHITE)
    ax.yaxis.label.set_color(WHITE)
    ax.title.set_color(CYAN)

# ═══ 1. Transverse map: initial positions ══════════════════════════
disk_r = meta["disk_radius_cm"]
for hc in hole_centres:
    circ = Circle(hc, R_hole, color="#1a2540", zorder=1)
    ax_map.add_patch(circ)
    circ2 = Circle(hc, R_hole, fill=False, edgecolor=CYAN, lw=0.8, alpha=0.7, zorder=2)
    ax_map.add_patch(circ2)

# Generation disk
gen_circ = Circle((0, 0), disk_r, fill=False, edgecolor=WHITE, lw=1, ls="--", alpha=0.5, zorder=2)
ax_map.add_patch(gen_circ)

ax_map.scatter(x0[lost], y0[lost], c=ROSE, s=6, alpha=0.5, label="Lost (anode)", zorder=3)
ax_map.scatter(x0[hit_wall], y0[hit_wall], c=AMBER, s=8, alpha=0.7, label="Hit wall", zorder=4)
ax_map.scatter(x0[reached_sipm], y0[reached_sipm], c=MINT, s=10, alpha=0.9, label="→ SiPM", zorder=5)

ax_map.set_xlim(-disk_r * 1.1, disk_r * 1.1)
ax_map.set_ylim(-disk_r * 1.1, disk_r * 1.1)
ax_map.set_aspect("equal")
ax_map.set_xlabel("x [cm]")
ax_map.set_ylabel("y [cm]")
ax_map.set_title("Transverse Plane — Electron Initial Positions & Fate", fontsize=10, fontweight="bold")
ax_map.legend(fontsize=8, loc="upper right", framealpha=0.3, labelcolor=WHITE, facecolor=PANEL)

# ═══ 2. Collection efficiency vs initial distance ══════════════════
bins = np.linspace(0, disk_r, 25)
bin_centres = 0.5 * (bins[:-1] + bins[1:])
eff_entered = []
eff_sipm = []
for lo, hi in zip(bins[:-1], bins[1:]):
    mask = (d_init >= lo) & (d_init < hi)
    if mask.sum() > 0:
        eff_entered.append(entered[mask].mean() * 100)
        eff_sipm.append(reached_sipm[mask].mean() * 100)
    else:
        eff_entered.append(np.nan)
        eff_sipm.append(np.nan)

ax_eff.bar(bin_centres, eff_entered, width=np.diff(bins) * 0.85,
           color=CYAN, alpha=0.4, edgecolor=DARK, label="Entered hole")
ax_eff.bar(bin_centres, eff_sipm, width=np.diff(bins) * 0.85,
           color=MINT, alpha=0.7, edgecolor=DARK, label="Reached SiPM")
ax_eff.axvline(R_hole, color=AMBER, ls=":", lw=1.5, label=f"R_hole = {R_hole} cm")
ax_eff.set_xlabel("Initial dist to nearest hole [cm]")
ax_eff.set_ylabel("Efficiency [%]")
ax_eff.set_title("Geometric Acceptance vs Distance", fontsize=10, fontweight="bold")
ax_eff.set_ylim(0, 110)
ax_eff.legend(fontsize=7, framealpha=0.3, labelcolor=WHITE, facecolor=PANEL)

# ═══ 3. Mean trajectories: z vs r displacement ═════════════════════
# Group by initial distance rings
d_bins = [0, R_hole * 0.5, R_hole, R_hole * 2, disk_r]
colors_traj = [MINT, CYAN, AMBER, ROSE]
labels_traj = [
    f"d < {d_bins[1]:.2f} cm (inside hole)",
    f"{d_bins[1]:.2f}–{d_bins[2]:.2f} cm (rim)",
    f"{d_bins[2]:.2f}–{d_bins[3]:.2f} cm (near wall)",
    f"{d_bins[3]:.2f}–{d_bins[4]:.2f} cm (far)",
]

for lo, hi, col, lab in zip(d_bins[:-1], d_bins[1:], colors_traj, labels_traj):
    mask = (d_init >= lo) & (d_init < hi)
    group = [electrons[i] for i in np.where(mask)[0]]
    if not group:
        continue

    # Plot z vs radial displacement from start
    all_dr = []
    all_z = []
    for e in group:
        traj = e["trajectory"]
        z_arr = np.array(traj["z"])
        dx = np.array(traj["x"]) - e["x0"]
        dy = np.array(traj["y"]) - e["y0"]
        dr = np.sqrt(dx**2 + dy**2)
        all_dr.append(dr)
        all_z.append(z_arr)

    # Interpolate to common z grid
    z_grid = np.linspace(meta["z_electron_start_cm"], 0, 50)
    dr_interp = []
    for z_arr, dr in zip(all_z, all_dr):
        # z is decreasing, flip for interp
        dr_interp.append(np.interp(z_grid, z_arr[::-1], dr[::-1]))

    dr_mat = np.array(dr_interp)
    mean_dr = np.nanmean(dr_mat, axis=0)
    std_dr = np.nanstd(dr_mat, axis=0)

    ax_traj.plot(z_grid, mean_dr, color=col, lw=2, label=lab)
    ax_traj.fill_between(z_grid, mean_dr - std_dr, mean_dr + std_dr,
                         color=col, alpha=0.15)

ax_traj.axvline(meta["z_hole_top_cm"], color=WHITE, lw=1, ls="--", alpha=0.5)
ax_traj.text(meta["z_hole_top_cm"] + 0.02, ax_traj.get_ylim()[1] * 0.9,
             "hole top", color=WHITE, fontsize=8, alpha=0.7)
ax_traj.set_xlabel("z [cm] (2.0 = start, 0 = SiPM)")
ax_traj.set_ylabel("⟨Δr⟩ from start [cm]")
ax_traj.set_title("Mean Radial Displacement vs Depth", fontsize=10, fontweight="bold")
ax_traj.legend(fontsize=8, framealpha=0.3, labelcolor=WHITE, facecolor=PANEL, loc="upper left")

# ═══ 4. Final positions at anode plane ═════════════════════════════
for hc in hole_centres:
    circ = Circle(hc, R_hole, fill=False, edgecolor=CYAN, lw=0.6, alpha=0.7)
    ax_final.add_patch(circ)

ax_final.scatter(xf[lost], yf[lost], c=ROSE, s=4, alpha=0.4, label="Lost")
ax_final.scatter(xf[hit_wall], yf[hit_wall], c=AMBER, s=6, alpha=0.6, label="Wall")
ax_final.scatter(xf[reached_sipm], yf[reached_sipm], c=MINT, s=8, alpha=0.8, label="SiPM")
ax_final.set_xlim(-disk_r * 1.1, disk_r * 1.1)
ax_final.set_ylim(-disk_r * 1.1, disk_r * 1.1)
ax_final.set_aspect("equal")
ax_final.set_xlabel("x_final [cm]")
ax_final.set_ylabel("y_final [cm]")
ax_final.set_title("Final Positions", fontsize=10, fontweight="bold")
ax_final.legend(fontsize=7, framealpha=0.3, labelcolor=WHITE, facecolor=PANEL)

# ═══ 5. Transverse displacement ════════════════════════════════════
disp = np.sqrt((xf - x0)**2 + (yf - y0)**2)
if reached_sipm.sum() > 0:
    ax_disp.hist(disp[reached_sipm], bins=20, color=MINT, alpha=0.7,
                 label="→ SiPM", density=True)
if lost.sum() > 0:
    ax_disp.hist(disp[lost], bins=20, color=ROSE, alpha=0.5,
                 label="Lost", density=True)
ax_disp.set_xlabel("Total transverse displacement [cm]")
ax_disp.set_ylabel("Probability density")
ax_disp.set_title("Transverse Displacement", fontsize=10, fontweight="bold")
ax_disp.legend(fontsize=7, framealpha=0.3, labelcolor=WHITE, facecolor=PANEL)

# ═══ 6. Photon yield per collected electron ════════════════════════
photons = np.array([e["total_photons"] for e in electrons])
if reached_sipm.sum() > 0:
    ax_phot.hist(photons[reached_sipm], bins=20, color=AMBER, alpha=0.8, edgecolor=DARK)
    mean_ph = photons[reached_sipm].mean()
    ax_phot.axvline(mean_ph, color=ROSE, lw=1.5, ls="--",
                    label=f"mean = {mean_ph:.0f}")
    ax_phot.legend(fontsize=8, framealpha=0.3, labelcolor=WHITE, facecolor=PANEL)
ax_phot.set_xlabel("EL photons per electron")
ax_phot.set_ylabel("Electrons")
ax_phot.set_title("Photon Yield (collected e⁻)", fontsize=10, fontweight="bold")

# ─── Title ─────────────────────────────────────────────────────────
fig.suptitle(
    f"NEXT Honeycomb Cell  ·  {meta['pressure_bar']} bar  ·  "
    f"EL {meta['el_field_v_cm']:.0f} V/cm  ·  Drift {meta['drift_field_v_cm']:.0f} V/cm  ·  "
    f"Acceptance {meta['geometric_acceptance_pct']:.1f}%  ·  "
    f"Eff. {meta['collection_efficiency_pct']:.1f}%",
    color=WHITE, fontsize=12, fontweight="bold", y=0.97
)

outpath = os.path.join(indir, "el_analysis.png")
fig.savefig(outpath, dpi=160, facecolor=DARK, bbox_inches="tight")
plt.close(fig)
print(f"  Plot saved → {outpath}")
