#!/usr/bin/env python3
"""
NEXT Honeycomb Cell — Electron Drift & EL Simulation (COMSOL field)
====================================================================
Simulates electron drift using the interpolated COMSOL E-field map.

Coordinate system (hole-centered):
  y = 0       : anode mesh (SiPM below)
  y ~ 0.3 cm  : top electrode of MTHGEM plate
  y = Y_START : electron generation plane (default 1.0 cm = 10 mm)
  (x', z')   : transverse plane, central hole at origin

Electrons drift along -y (from Y_START toward y=0).
The E-field has 3 components (Ex, Ey, Ez) from the COMSOL map,
providing funneling of electrons into the hole.

Output: el_electrons.json with trajectories and excitation points.
"""

import numpy as np
import json
import os
from field_map import FieldMap

# ─── Pressure & EL yield ────────────────────────────────────────
PRESSURE_BAR       = 4.0            # bar
EL_THRESHOLD_RED   = 0.83           # kV/(cm·bar) — EL threshold in Xe
EL_YIELD_COEFF     = 140.0          # photons/(e·cm·bar) per unit (E/P - threshold)

# ─── Hole geometry (from COMSOL analysis) ───────────────────────
HOLE_RADIUS_CM     = 0.25           # 2.5 mm
HOLE_HEIGHT_CM     = 0.30           # 3 mm plate thickness (EL gap)
HOLE_PITCH_CM      = 1.00           # 10 mm hexagonal pitch

# ─── Drift geometry ────────────────────────────────────────────
Y_ANODE            = 0.0            # y of anode mesh
Y_PLATE_TOP        = HOLE_HEIGHT_CM # y of top electrode (~3 mm)
Y_ELECTRON_START   = 1.0            # y where electrons start (10 mm)

# ─── Diffusion coefficients ────────────────────────────────────
DIFF_TRANS_DRIFT   = 0.05           # cm/√cm — drift region (~100 V/cm)
DIFF_TRANS_EL      = 0.01           # cm/√cm — EL region (~4700 V/cm)

# ─── Simulation parameters ─────────────────────────────────────
N_ELECTRONS        = 10000
DISK_RADIUS_CM     = HOLE_PITCH_CM  # 1.0 cm — covers full Voronoi cell
STEP_SIZE_CM       = 0.002          # Δy per step [cm]
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
                cz = ring * pitch * np.sin(angle0) + k * pitch * np.sin(angle1)
                centres.append((cx, cz))
    return np.array(centres)


def nearest_hole(xp, zp, hole_centres):
    """Return distance to nearest hole centre and its index."""
    dx = hole_centres[:, 0] - xp
    dz = hole_centres[:, 1] - zp
    d2 = dx**2 + dz**2
    idx = np.argmin(d2)
    return np.sqrt(d2[idx]), idx


def simulate_electron(x0, z0, field_map, hole_centres, rng):
    """
    Simulate one electron drifting from y=Y_ELECTRON_START toward y=0.

    At each step:
      - Get E(x', y, z') from the COMSOL field map
      - Advance y by -STEP_SIZE_CM
      - Compute transverse displacement from field ratios: dx = (Ex/Ey)*dy
      - Add Gaussian diffusion in x' and z'
      - Check boundary conditions (hole wall, anode plate, SiPM)
    """
    y = Y_ELECTRON_START
    xp, zp = x0, z0
    steps_x, steps_z, steps_y = [float(xp)], [float(zp)], [float(y)]
    excitations = []

    n_steps_total = int(Y_ELECTRON_START / STEP_SIZE_CM)
    save_every = max(1, n_steps_total // N_TRAJECTORY_SAVE)

    step_count = 0
    entered_hole = False
    hit_wall = False
    reached_sipm = False
    hole_idx = -1

    dy = -STEP_SIZE_CM  # negative: drifting toward y=0

    while y > 0 and not hit_wall:
        # Get E-field
        Ex, Ey, Ez = field_map.get_E(xp, y, zp)

        # Avoid division by zero
        if abs(Ey) < 1.0:
            Ey_safe = -1.0  # fallback: small drift field
        else:
            Ey_safe = Ey

        # Transverse displacement from field
        dx_field = (Ex / Ey_safe) * dy
        dz_field = (Ez / Ey_safe) * dy

        # Diffusion
        in_hole = (y < Y_PLATE_TOP) and entered_hole
        sigma_t = DIFF_TRANS_EL if in_hole else DIFF_TRANS_DRIFT
        sigma = sigma_t * np.sqrt(STEP_SIZE_CM)
        dx_diff = rng.normal(0, sigma)
        dz_diff = rng.normal(0, sigma)

        # Update position
        xp += dx_field + dx_diff
        zp += dz_field + dz_diff
        y += dy

        if y < 0:
            y = 0.0

        step_count += 1

        # ── Check if entering hole ──
        if not entered_hole and y <= Y_PLATE_TOP:
            dist_nearest, hole_idx = nearest_hole(xp, zp, hole_centres)
            if dist_nearest < HOLE_RADIUS_CM:
                entered_hole = True
            else:
                # Hits anode plate outside hole
                break

        # ── Inside hole: check wall ──
        if entered_hole and y <= Y_PLATE_TOP:
            hx, hz = hole_centres[hole_idx]
            r_rel = np.sqrt((xp - hx)**2 + (zp - hz)**2)
            if r_rel >= HOLE_RADIUS_CM:
                hit_wall = True
                steps_x.append(float(xp))
                steps_z.append(float(zp))
                steps_y.append(float(y))
                break

            # EL photon production
            E_mag = np.sqrt(Ex**2 + Ey**2 + Ez**2)  # V/cm
            E_reduced = E_mag / (PRESSURE_BAR * 1e3)  # kV/(cm·bar)
            if E_reduced > EL_THRESHOLD_RED:
                photons_per_cm = EL_YIELD_COEFF * (E_reduced - EL_THRESHOLD_RED) * PRESSURE_BAR
                n_ph = rng.poisson(photons_per_cm * STEP_SIZE_CM)
                if n_ph > 0:
                    excitations.append({
                        "x": float(xp - hx),
                        "z": float(zp - hz),
                        "y": float(y),
                        "n_photons": int(n_ph)
                    })

        # Save trajectory point
        if step_count % save_every == 0:
            steps_x.append(float(xp))
            steps_z.append(float(zp))
            steps_y.append(float(y))

    if entered_hole and not hit_wall and y <= 0:
        reached_sipm = True

    # Final position
    steps_x.append(float(xp))
    steps_z.append(float(zp))
    steps_y.append(float(y))

    dist_nearest, nearest_idx = nearest_hole(xp, zp, hole_centres)

    total_photons = sum(e["n_photons"] for e in excitations)

    return {
        "x0": float(x0),
        "z0": float(z0),
        "y0": float(Y_ELECTRON_START),
        "x_final": float(xp),
        "z_final": float(zp),
        "y_final": float(y),
        "entered_hole": bool(entered_hole),
        "hit_wall": bool(hit_wall),
        "reached_sipm": bool(reached_sipm),
        "hole_index": int(hole_idx) if entered_hole else -1,
        "dist_to_nearest_hole_cm": float(dist_nearest),
        "total_photons": int(total_photons),
        "trajectory": {
            "x": steps_x,
            "z": steps_z,
            "y": steps_y,
        },
        "excitation_points": excitations,
    }


def main():
    print("Loading COMSOL field map...")
    field_map = FieldMap()

    rng = np.random.default_rng(SEED)
    hole_centres = make_hex_holes(HOLE_PITCH_CM, n_rings=2)
    print(f"Hexagonal lattice: {len(hole_centres)} holes")

    # Random initial positions in disk of radius DISK_RADIUS_CM
    r = DISK_RADIUS_CM * np.sqrt(rng.uniform(0, 1, N_ELECTRONS))
    theta = rng.uniform(0, 2 * np.pi, N_ELECTRONS)
    x0s = r * np.cos(theta)
    z0s = r * np.sin(theta)

    electrons = []
    for i, (x0, z0) in enumerate(zip(x0s, z0s)):
        if (i + 1) % 100 == 0 or i == 0:
            print(f"  Electron {i+1}/{N_ELECTRONS}...", end="\r", flush=True)
        e = simulate_electron(x0, z0, field_map, hole_centres, rng)
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
        "description": "NEXT honeycomb cell EL simulation (COMSOL field map)",
        "pressure_bar": PRESSURE_BAR,
        "hole_radius_cm": HOLE_RADIUS_CM,
        "hole_height_cm": HOLE_HEIGHT_CM,
        "hole_pitch_cm": HOLE_PITCH_CM,
        "y_electron_start_cm": Y_ELECTRON_START,
        "y_plate_top_cm": Y_PLATE_TOP,
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
