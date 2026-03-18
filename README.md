# NEXT Honeycomb Cell — EL Simulation
## Electroluminescence in cylindrical holes · Honeycomb geometry

---

## Overview

Full simulation chain for a single honeycomb cell in a NEXT-like TPC:
electron drift → EL photon generation inside hole → VUV photon propagation → SiPM acceptance.

```
Electron starts at z=2 cm (drift region)
         │  drift field: 100 V/cm
         ▼
   Hole top (z=0.5 cm) — enters hole?
         │  EL field: 4800 V/cm (1200 V/(cm·bar) × 4 bar)
         ▼
   Excitation points → VUV photons (172 nm)
         │
         ▼  el_electrons.json
   Geant4 propagates photons in hole
   (PTFE walls 95% reflective, SiPM at bottom)
         │
         ▼  output.root
   Analysis: acceptance, reflections, efficiency
```

---

## Geometry

```
              ┌─── Disk r=10 mm ───┐
              │  electron generation │
              └──────────┬──────────┘
                         │ drift 15 mm (100 V/cm)
    ╔═══════╤═══════════╧═══════════╤═══════╗  z = 5 mm (hole top)
    ║ PTFE  │     Xe gas (EL)       │ PTFE  ║
    ║ wall  │   Ø5 mm × 5 mm       │ wall  ║  EL: 4800 V/cm
    ╚═══════╧═══════════╤═══════════╧═══════╝  z = 0 (hole bottom)
                    ┌───┴───┐
                    │ SiPM  │
                    └───────┘
```

- **Hole**: Ø 5 mm, height 5 mm
- **Pitch**: 10 mm (hexagonal honeycomb)
- **Wall**: ~5 mm between hole edges

---

## Files

| File | Description |
|------|-------------|
| `el_simulation.py` | Step 1: electron drift + EL excitation |
| `el_analysis.py` | Step 2: trajectory analysis & efficiency plots |
| `ELHolePhotonSim.cc` | Step 3: Geant4 VUV photon propagation |
| `CMakeLists.txt` | Geant4 build configuration |
| `analyze_photons.py` | Step 4: photon statistics & acceptance |
| `cell_diagram.py` | Generates geometry diagram |

---

## Step 1 — Electron Drift

```bash
python3 el_simulation.py
```

- 100 electrons, uniformly distributed in disk (r=10 mm) at z=2 cm
- Drift region (z > 0.5 cm): 100 V/cm, transverse diffusion
- Hole (0 < z < 0.5 cm): 4800 V/cm, EL photon production
- Output: `el_electrons.json`

**EL yield** at 4 bar:
- Y = 140 × (E/P − 0.83) × P ≈ 207 photons/(e⁻·cm)
- Over 5 mm gap: ~104 photons per electron

---

## Step 2 — Analysis

```bash
python3 el_analysis.py
```

Output: `el_analysis.png` with 6 panels showing geometric acceptance,
trajectories, displacement, and photon yield.

---

## Step 3 — Geant4 Photon Simulation

```bash
mkdir build && cd build
cmake .. -DGeant4_DIR=/path/to/geant4/lib/Geant4-11.X.X
make -j$(nproc)
./ELHolePhotonSim ../el_electrons.json output.root
```

Propagates VUV photons (172 nm) through the hole:
- PTFE walls: 95% diffuse Lambertian reflectivity
- SiPM at bottom: perfect absorber
- Tracks: reflections, arrival position, fate (SiPM/wall/escape)

---

## Step 4 — Photon Analysis

```bash
python3 analyze_photons.py output.root
```

Output: `photon_analysis.png` — acceptance, reflections, per-electron efficiency.

---

## Parameters (baseline: 4 bar)

| Parameter | Value |
|-----------|-------|
| Hole diameter | 5 mm |
| Hole height | 5 mm |
| Pitch | 10 mm |
| EL reduced field | 1200 V/(cm·bar) |
| EL field (4 bar) | 4800 V/cm |
| Drift field | 100 V/cm |
| Pressure | 4 bar (also: 10, 13.5 bar) |
| VUV wavelength | 172 nm (Xe) |
| PTFE reflectivity | 95% |

---

## Dependencies

```bash
pip install numpy matplotlib scipy uproot awkward
# Geant4 >= 10.7 (for C++ photon simulation)
```
