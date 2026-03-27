# NEXT Honeycomb Cell EL Simulation

## Project Goal

Simulate the electroluminescence (EL) response of a single cell in a NEXT-like TPC
with honeycomb hole structure. The full chain: electron drift → EL photon generation
inside hole → Geant4 VUV photon propagation → SiPM acceptance.

## Geometry (from COMSOL field analysis)

- **Hole**: cylindrical, diameter 5 mm (radius 2.5 mm), depth ~3 mm (EL gap)
- **Pitch**: 10 mm hexagonal honeycomb
- **Wall thickness**: ~5 mm between hole edges (PTFE)
- **SiPM**: below anode mesh at y = 0
- **Anode mesh electrode**: y ~ 0, V = 0
- **Top electrode**: y ~ 3 mm, V ~ -1500 V
- **Cathode**: y ~ 103 mm, V = -2500 V

## Physics Parameters

| Parameter              | Value                |
|------------------------|----------------------|
| Hole radius            | 2.5 mm               |
| Hole depth (EL gap)    | ~3 mm                |
| Pitch                  | 10 mm                |
| EL field (inside hole) | ~4700 V/cm           |
| Drift field            | ~105 V/cm            |
| Pressure (baseline)    | 4 bar                |
| Pressure (planned)     | 10 bar, 13.5 bar     |
| VUV wavelength (Xe)    | 172 nm               |
| PTFE reflectivity      | ~95%                 |

## Coordinate System (hole-centered)

- y = 0: anode mesh (SiPM below)
- y = 0.3 cm: top electrode of MTHGEM plate
- y = 1.0 cm: electron start position (default)
- Electrons drift along -y (y decreasing)
- (x', z') = transverse plane, central hole at origin (0, 0)

## COMSOL Field Map

- **File**: `data/Mesh-MTHGEM_fields_4bar.txt` (COMSOL 5.5, 544k nodes)
- **Domain**: rectangle from central hole (origin) to neighbor hole at 150°
  - x' in [-8.70, 0] mm, z' in [0, 5.00] mm
  - dx = pitch × cos(30°), dz = pitch × sin(30°)
- **Symmetry**: 30° sector (150°-180°), extended by hexagonal 6-fold symmetry + mirror
- **Interpolation**: `field_map.py` — 2D per-slab griddata with fold-back
- **Analysis**: `data/field_analysis.md`

## Simulation Chain

1. `field_map.py` — COMSOL field interpolator with hexagonal fold-back
2. `el_simulation.py` — Electron drift (COMSOL E-field) + EL excitation → `el_electrons.json`
3. `el_analysis.py` — Trajectory analysis + efficiency plots → `el_analysis.png`
4. `ELHolePhotonSim.cc` (Geant4) — VUV photon propagation → `output.root`
5. `analyze_photons.py` — Photon statistics + acceptance → `photon_analysis.png`

## Electron Generation

- Electrons generated uniformly in a disk of radius 10 mm (= pitch) centered on the hole
- Start at y = 10 mm (configurable via Y_ELECTRON_START)
- Covers the full Voronoi cell
- Gives geometric acceptance directly: n_entering_hole / n_total

## Key Design Decisions

- EL happens INSIDE the hole (not in a gap above)
- E-field from COMSOL interpolation (3 components: Ex, Ey, Ez)
- Transverse field components provide electron funneling into the hole
- Electron stops if it hits the anode plate (outside hole), hole wall, or reaches SiPM (y=0)
- EL yield computed from local |E|, not uniform field
- Pressure is parameterizable for scans at 4, 10, 13.5 bar

## Pending Tasks

- [x] Integrate COMSOL electric field map into `el_simulation.py`
- [ ] Pressure scans at 10 bar and 13.5 bar
- [ ] Explore different hole geometries (larger diameter, smaller pitch)

## Dependencies

```bash
pip install numpy matplotlib scipy uproot awkward
# Geant4 >= 10.7 (for C++ photon simulation)
```
