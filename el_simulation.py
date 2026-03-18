#!/usr/bin/env python3
"""
NEXT Honeycomb Cell — Electron Drift & EL Simulation
=====================================================
Simulates electron drift through two regions:
  1. Drift region (above hole): uniform drift field, no EL
  2. EL region (inside hole): strong EL field, photon production

Electrons start uniformly distributed in a disk of radius = pitch (10 mm)
centred on the hole, at z = 2 cm from the SiPM (hole bottom).

Coordinate system:
  z = 0       : SiPM plane (hole bottom)
  z = 0.5 cm  : hole top (anode plane)
  z = 2.0 cm  : electron start

Output: el_electrons.json with trajectories and excitation points.
"""

import numpy as np
import json
import os

# ─── Pressure & Fields ───────────────────────────────────────────
PRESSURE_BAR       = 4.0            # bar
EL_REDUCED_FIELD   = 1.2            # kV/(cm·bar) — reduced EL field
DRIFT_FIELD_V_CM   = 100.0          # V/cm — drift field above hole
EL_FIELD_V_CM      = EL_REDUCED_FIELD * 1e3 * PRESSURE_BAR  # V/cm (4800 at 4 bar)

# ─── Hole geometry ────────────────────────────────────────────────
HOLE_RADIUS_CM     = 0.25           # 2.5 mm radius (5 mm diameter)
HOLE_HEIGHT_CM     = 0.50           # 5 mm height (EL gap)
HOLE_PITCH_CM      = 1.00           # 10 mm hexagonal pitch

# ─── Drift geometry ──────────────────────────────────────────────
Z_SIPM             = 0.0            # z of SiPM (hole bottom)
Z_HOLE_TOP         = HOLE_HEIGHT_CM # z of hole top = 0.5 cm
Z_ELECTRON_START   = 2.0            # z where electrons start = 2.0 cm
DRIFT_LENGTH_CM    = Z_ELECTRON_START - Z_HOLE_TOP  # 1.5 cm of drift

# ─── EL yield ────────────────────────────────────────────────────
# Y = 140 * (E/P - 0.83) * P  [photons/cm]  (E/P in kV/(cm·bar))
EL_THRESHOLD_RED   = 0.83           # kV/(cm·bar) — EL threshold in Xe
EL_YIELD_COEFF     = 140.0          # photons/(e·cm·bar) per unit (E/P - threshold)
PHOTONS_PER_CM     = EL_YIELD_COEFF * (EL_REDUCED_FIELD - EL_THRESHOLD_RED) * PRESSURE_BAR

# ─── Diffusion coefficients ──────────────────────────────────────
# σ_T per step = DIFF_COEFF * sqrt(Δz)
# Higher diffusion at low field (drift), lower at high field (EL)
DIFF_TRANS_DRIFT   = 0.05           # cm/√cm — drift region (100 V/cm)
DIFF_TRANS_EL      = 0.01           # cm/√cm — EL region (4800 V/cm)

# ─── Simulation parameters ───────────────────────────────────────
N_ELECTRONS        = 100
DISK_RADIUS_CM     = HOLE_PITCH_CM  # 1.0 cm — generation disk radius
STEP_SIZE_CM       = 0.002          # Δz per step [cm]
N_TRAJECTORY_SAVE  = 40             # trajectory points saved per electron
SEED               = 42


def make_hex_holes(pitch, n_rings=2):
    """Generate hexagonal lattice of hole centres around the central hole."""
    centres = [(0.0, 0.0)]
    for ring in range(1, n_rings + 1):
        for i in range(6):
            angle0 = np.pi / 3 * i
            for k in range(ring):
                angle1 = np.pi / 3 * (i + 2)
                cx = ring * pitch * np.cos(angle0) + k * pitch * np.cos(angle1)
                cy = ring * pitch * np.sin(angle0) + k * pitch * np.sin(angle1)
                centres.append((cx, cy))
    return np.array(centres)


def nearest_hole(x, y, hole_centres):
    """Return distance to nearest hole centre and its index."""
    dx = hole_centres[:, 0] - x
    dy = hole_centres[:, 1] - y
    d2 = dx**2 + dy**2
    idx = np.argmin(d2)
    return np.sqrt(d2[idx]), idx


def simulate_electron(x0, y0, hole_centres, rng):
    """
    Simulate one electron drifting from z=Z_ELECTRON_START downward.

    Two regions:
      - Drift (z > HOLE_HEIGHT): field=100 V/cm, diffusion only, no EL
      - Hole  (0 < z < HOLE_HEIGHT): field=4800 V/cm, EL production

    Electron stops if:
      - It reaches z=0 (SiPM) — collected
      - At z=HOLE_HEIGHT, it's outside the hole — hits anode plate
      - Inside hole, r > HOLE_RADIUS — hits wall
    """
    z = Z_ELECTRON_START
    x, y = x0, y0
    steps_x, steps_y, steps_z = [float(x)], [float(y)], [float(z)]
    excitations = []

    # Total steps
    n_steps_drift = int(DRIFT_LENGTH_CM / STEP_SIZE_CM)
    n_steps_hole = int(HOLE_HEIGHT_CM / STEP_SIZE_CM)
    n_steps_total = n_steps_drift + n_steps_hole
    save_every = max(1, n_steps_total // N_TRAJECTORY_SAVE)

    step_count = 0
    entered_hole = False
    hit_wall = False
    reached_sipm = False

    # ── Phase 1: Drift region (z from Z_ELECTRON_START down to Z_HOLE_TOP) ──
    for i in range(n_steps_drift):
        sigma_t = DIFF_TRANS_DRIFT * np.sqrt(STEP_SIZE_CM)
        x += rng.normal(0, sigma_t)
        y += rng.normal(0, sigma_t)
        z -= STEP_SIZE_CM

        step_count += 1
        if step_count % save_every == 0:
            steps_x.append(float(x))
            steps_y.append(float(y))
            steps_z.append(float(z))

    # ── Check if electron enters hole at z = HOLE_HEIGHT ──
    dist_to_center = np.sqrt(x**2 + y**2)  # distance to central hole
    # Find nearest hole
    dist_nearest, hole_idx = nearest_hole(x, y, hole_centres)

    if dist_nearest < HOLE_RADIUS_CM:
        entered_hole = True
        # Shift coordinates so hole centre is at (0,0) for excitation points
        hx, hy = hole_centres[hole_idx]
        x_rel = x - hx
        y_rel = y - hy

        # ── Phase 2: Inside hole (z from HOLE_HEIGHT down to 0) ──
        for i in range(n_steps_hole):
            sigma_t = DIFF_TRANS_EL * np.sqrt(STEP_SIZE_CM)
            x += rng.normal(0, sigma_t)
            y += rng.normal(0, sigma_t)
            x_rel += x - (steps_x[-1] if steps_x else x0)  # track relative motion
            y_rel = y - hy
            x_rel = x - hx
            z -= STEP_SIZE_CM

            # Check wall hit
            r_rel = np.sqrt((x - hx)**2 + (y - hy)**2)
            if r_rel >= HOLE_RADIUS_CM:
                hit_wall = True
                step_count += 1
                steps_x.append(float(x))
                steps_y.append(float(y))
                steps_z.append(float(z))
                break

            # EL excitation — photons proportional to dz
            n_ph = rng.poisson(PHOTONS_PER_CM * STEP_SIZE_CM)
            if n_ph > 0:
                excitations.append({
                    "x": float(x - hx),  # relative to hole centre
                    "y": float(y - hy),
                    "z": float(z),        # z from SiPM
                    "n_photons": int(n_ph)
                })

            step_count += 1
            if step_count % save_every == 0:
                steps_x.append(float(x))
                steps_y.append(float(y))
                steps_z.append(float(z))

        if not hit_wall:
            reached_sipm = True

    # Final position
    steps_x.append(float(x))
    steps_y.append(float(y))
    steps_z.append(float(z))

    total_photons = sum(e["n_photons"] for e in excitations)

    return {
        "x0": float(x0),
        "y0": float(y0),
        "z0": float(Z_ELECTRON_START),
        "x_final": float(x),
        "y_final": float(y),
        "z_final": float(z),
        "entered_hole": bool(entered_hole),
        "hit_wall": bool(hit_wall),
        "reached_sipm": bool(reached_sipm),
        "hole_index": int(hole_idx) if entered_hole else -1,
        "dist_to_nearest_hole_cm": float(dist_nearest),
        "total_photons": int(total_photons),
        "trajectory": {
            "x": steps_x,
            "y": steps_y,
            "z": steps_z,
        },
        "excitation_points": excitations,
    }


def main():
    rng = np.random.default_rng(SEED)
    hole_centres = make_hex_holes(HOLE_PITCH_CM, n_rings=2)

    # Random initial positions in disk of radius DISK_RADIUS_CM
    # Uniform in disk: r = R * sqrt(U), theta = 2π * V
    r = DISK_RADIUS_CM * np.sqrt(rng.uniform(0, 1, N_ELECTRONS))
    theta = rng.uniform(0, 2 * np.pi, N_ELECTRONS)
    x0s = r * np.cos(theta)
    y0s = r * np.sin(theta)

    electrons = []
    for i, (x0, y0) in enumerate(zip(x0s, y0s)):
        print(f"  Simulating electron {i+1}/{N_ELECTRONS}...", end="\r", flush=True)
        e = simulate_electron(x0, y0, hole_centres, rng)
        e["electron_id"] = i
        electrons.append(e)
    print(f"\n  Done. {N_ELECTRONS} electrons simulated.")

    # Summary
    n_entered = sum(1 for e in electrons if e["entered_hole"])
    n_wall = sum(1 for e in electrons if e["hit_wall"])
    n_sipm = sum(1 for e in electrons if e["reached_sipm"])
    n_lost = N_ELECTRONS - n_entered
    total_photons = sum(e["total_photons"] for e in electrons)

    print(f"\n  === Results ===")
    print(f"  Entered hole   : {n_entered}/{N_ELECTRONS} ({n_entered/N_ELECTRONS*100:.1f}%)")
    print(f"  Hit wall       : {n_wall}/{N_ELECTRONS} ({n_wall/N_ELECTRONS*100:.1f}%)")
    print(f"  Reached SiPM   : {n_sipm}/{N_ELECTRONS} ({n_sipm/N_ELECTRONS*100:.1f}%)")
    print(f"  Lost (anode)   : {n_lost}/{N_ELECTRONS} ({n_lost/N_ELECTRONS*100:.1f}%)")
    print(f"  Total EL photons: {total_photons}")
    if n_sipm > 0:
        avg_ph = total_photons / n_sipm
        print(f"  Avg photons/collected e⁻: {avg_ph:.1f}")

    metadata = {
        "description": "NEXT honeycomb cell EL simulation",
        "pressure_bar": PRESSURE_BAR,
        "el_reduced_field_kv_cm_bar": EL_REDUCED_FIELD,
        "el_field_v_cm": EL_FIELD_V_CM,
        "drift_field_v_cm": DRIFT_FIELD_V_CM,
        "hole_radius_cm": HOLE_RADIUS_CM,
        "hole_height_cm": HOLE_HEIGHT_CM,
        "hole_pitch_cm": HOLE_PITCH_CM,
        "z_electron_start_cm": Z_ELECTRON_START,
        "z_hole_top_cm": Z_HOLE_TOP,
        "drift_length_cm": DRIFT_LENGTH_CM,
        "photons_per_cm": round(PHOTONS_PER_CM, 2),
        "diff_trans_drift": DIFF_TRANS_DRIFT,
        "diff_trans_el": DIFF_TRANS_EL,
        "step_size_cm": STEP_SIZE_CM,
        "disk_radius_cm": DISK_RADIUS_CM,
        "n_holes": len(hole_centres),
        "hole_centres": hole_centres.tolist(),
        "n_electrons": N_ELECTRONS,
        "n_entered_hole": n_entered,
        "n_hit_wall": n_wall,
        "n_reached_sipm": n_sipm,
        "geometric_acceptance_pct": round(n_entered / N_ELECTRONS * 100, 2),
        "collection_efficiency_pct": round(n_sipm / N_ELECTRONS * 100, 2),
        "total_el_photons": total_photons,
        "random_seed": SEED,
    }

    output = {"metadata": metadata, "electrons": electrons}

    outdir = os.path.dirname(os.path.abspath(__file__))
    outpath = os.path.join(outdir, "el_electrons.json")
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n  Saved → {outpath}")
    return output


if __name__ == "__main__":
    main()
