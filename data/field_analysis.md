# COMSOL Field Data Analysis: Mesh-MTHGEM_fields_4bar.txt

## File Structure

- **Model**: `Mesh-MTHGEM.mph` (COMSOL 5.5)
- **Dimension**: 3D, 544,638 nodes
- **Quantity**: Electric potential V (Volts)
- **Length unit**: mm
- **Format**: 4 columns: x, y, z, V(V) — whitespace-separated, 9-line header

## Coordinate Ranges

| Axis | Min     | Max     | Span    |
|------|---------|---------|---------|
| x    | -5.33   | +3.37   | 8.7 mm  |
| y    | -0.85   | +103.0  | ~104 mm |
| z    | 6.20    | 11.20   | 5.0 mm  |
| V    | -2500   | 0       | 2500 V  |

## Deduced Geometry

### Coordinate System

- **y-axis = drift direction** (parallel to hole axis, field points along -y)
- **(x, z) = transverse plane** perpendicular to drift
- y = 0: anode (mesh electrode, V = 0)
- y = 103: cathode (V = -2500 V)
- Electrons drift from cathode (high y) toward anode (y = 0)

### MTHGEM Plate

- **Bottom electrode (mesh/anode)**: y ~ -0.3 mm, V = 0
- **Top electrode**: y ~ 3 mm, V ~ -1500 V
- **Plate thickness between electrodes**: ~3 mm
- **Hole**: cylindrical, axis along y, radius 2.5 mm (5 mm diameter)
- **Hole depth (EL gap)**: ~3 mm (y = 0 to y ~ 3)

### Electric Fields

| Region              | y range (mm)  | E (V/cm)  | E/P at 4 bar       |
|---------------------|---------------|-----------|---------------------|
| Inside hole (EL)    | 0 to ~3       | ~4700     | ~1175 V/(cm·bar)    |
| Drift               | ~5 to 103     | ~105      | ~26 V/(cm·bar)      |

### Simulation Domain in (x, z) — Minimal Hexagonal Sector

The domain is a **rectangle with holes at two diagonally opposite corners**:

- **Corner A** (x_min, z_min): center of **hole A** (central hole)
- **Corner B** (x_max, z_max): center of **hole B** (neighbor at 30°)
- **Diagonal A→B**: 10.03 mm = pitch, angle = 29.9°
- **x span**: 8.70 mm ≈ pitch × cos(30°) = 8.66 mm
- **z span**: 5.00 mm = pitch × sin(30°) = 5.00 mm
- **Edges along x=0 and z=0** (in model coords) = symmetry planes

This rectangle is the **minimal domain** covering the 0°–30° sector. By mirror reflections in x and z, the full hexagonal cell around each hole is reconstructed.

Each hole appears as a **quarter-circle** at its corner of the rectangle. The PTFE walls occupy the remaining area.

### Coordinate Offset (data vs model)

The COMSOL data file has an offset from the model coordinates:
- x_model = x_data + 5.33 (hole A at x_model = 0)
- z_model = z_data - 6.20 (hole A at z_model = 0)

### Boundary Conditions

- **V = 0** at y ~ -0.3 (anode mesh electrode, all x, z)
- **V = -2500** at y ~ 103 (cathode plane, all x, z)
- **Lateral boundaries** (x = 0, z = 0 in model coords): **mirror symmetry**

### Voltage Budget

- Total: 2500 V
- Across EL gap (~3 mm): ~1500 V
- Across drift (~100 mm): ~1000 V

### Node Density (mesh refinement)

The COMSOL mesh is heavily refined near the electrodes:

| y range (mm) | Nodes   | Region                    |
|--------------|---------|---------------------------|
| [-1, 0]      | 322,680 | Anode/mesh electrode       |
| [0, 1]       | 30,758  | Inside hole (near anode)   |
| [1, 2]       | 5,818   | Inside hole (mid)          |
| [2, 3]       | 60,978  | Top electrode of plate     |
| [3, 4]       | 38,717  | Just above plate           |
| [4, 5]       | 3,351   | Transition to drift        |
| [5, 100]     | ~2,000  | Drift region (sparse)      |
| [101, 103]   | 75,507  | Cathode                    |
