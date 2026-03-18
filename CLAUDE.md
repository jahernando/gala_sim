# NEXT Honeycomb Cell EL Simulation

## Project Goal

Simulate the electroluminescence (EL) response of a single cell in a NEXT-like TPC
with honeycomb hole structure. The full chain: electron drift → EL photon generation
inside hole → Geant4 VUV photon propagation → SiPM acceptance.

## Geometry

- **Hole**: cylindrical, diameter 5 mm (radius 2.5 mm), height 5 mm
- **Pitch**: 10 mm hexagonal honeycomb
- **Wall thickness**: ~5 mm between hole edges (PTFE)
- **SiPM**: at the bottom of each hole (z=0)
- **Fill factor**: ~5.7%

## Physics Parameters

| Parameter              | Value                |
|------------------------|----------------------|
| Hole radius            | 2.5 mm               |
| Hole height (EL gap)   | 5 mm                 |
| Pitch                  | 10 mm                |
| EL reduced field       | 1200 V/(cm·bar)      |
| Drift field            | 100 V/cm             |
| Pressure (baseline)    | 4 bar                |
| Pressure (planned)     | 10 bar, 13.5 bar     |
| EL real field (4 bar)  | 4800 V/cm            |
| VUV wavelength (Xe)    | 172 nm               |
| PTFE reflectivity      | ~95%                 |

## Coordinate System

- z = 0: SiPM plane (hole bottom)
- z = 0.5 cm: hole top (anode plane)
- z = 2.0 cm: electron start position
- Electrons drift downward (z decreasing)
- Hole centered at (x, y) = (0, 0)

## Simulation Chain

1. `el_simulation.py` — Electron drift + EL excitation → `el_electrons.json`
2. `el_analysis.py` — Trajectory analysis + efficiency plots → `el_analysis.png`
3. `ELHolePhotonSim.cc` (Geant4) — VUV photon propagation → `output.root`
4. `analyze_photons.py` — Photon statistics + acceptance → `photon_analysis.png`

## Electron Generation

- Electrons generated uniformly in a disk of radius 10 mm (= pitch) centered on the hole
- This covers the full Voronoi cell + "no-man's land" triangles between hexagons
- Gives geometric acceptance directly: n_entering_hole / n_total

## Key Design Decisions

- EL happens INSIDE the hole (not in a gap above)
- Two drift regions: drift (100 V/cm above hole) + EL (4800 V/cm inside hole)
- Electron stops if it hits the anode plate (outside hole) or hole wall
- Future: user may provide detailed electric field map (potential) to replace uniform fields
- Pressure is parameterizable for scans at 4, 10, 13.5 bar

## Dependencies

```bash
pip install numpy matplotlib scipy uproot awkward
# Geant4 >= 10.7 (for C++ photon simulation)
```
