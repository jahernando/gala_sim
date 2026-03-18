"""
Diagram of a single NEXT honeycomb cell for the EL hole simulation.
Shows both a top view (hexagonal cell with hole) and a side cross-section.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch

fig, axes = plt.subplots(1, 2, figsize=(14, 7), facecolor='white')

# ============================================================
# LEFT: Top view — hexagonal cell with central hole
# ============================================================
ax = axes[0]
ax.set_aspect('equal')
ax.set_title('Top View — Hexagonal Cell', fontsize=14, fontweight='bold')

pitch = 10.0  # mm
hole_diameter = 5.0  # mm
hole_radius = hole_diameter / 2.0

# Draw hexagonal cell (flat-top orientation)
# For flat-top hex with pitch = center-to-center distance:
# The hex circumradius (center to vertex) = pitch / sqrt(3) for pointy-top
# For a honeycomb with pitch=10mm center-to-center, the hex apothem = pitch/2 = 5mm
apothem = pitch / 2.0  # distance from center to flat side = 5 mm
circumradius = apothem / np.cos(np.pi / 6)  # = apothem / (sqrt(3)/2) ≈ 5.77 mm

# Hexagon vertices (flat-top: first vertex at 0°)
hex_angles = np.linspace(0, 2 * np.pi, 7) + np.pi / 6  # rotate 30° for flat-top
hex_x = circumradius * np.cos(hex_angles)
hex_y = circumradius * np.sin(hex_angles)
ax.fill(hex_x, hex_y, color='#D2B48C', alpha=0.4, label='PTFE/copper wall')
ax.plot(hex_x, hex_y, 'k-', linewidth=2)

# Draw the hole
circle = plt.Circle((0, 0), hole_radius, color='lightblue', alpha=0.6, label=f'Hole (Ø{hole_diameter} mm)')
ax.add_patch(circle)
circle_edge = plt.Circle((0, 0), hole_radius, fill=False, edgecolor='blue', linewidth=2)
ax.add_patch(circle_edge)

# Draw neighboring holes (to show the honeycomb pattern)
neighbor_angles = np.linspace(0, 2 * np.pi, 7)[:-1]  # 6 neighbors
for angle in neighbor_angles:
    nx = pitch * np.cos(angle + np.pi / 6)
    ny = pitch * np.sin(angle + np.pi / 6)
    # Draw partial hex
    nh_x = nx + circumradius * np.cos(hex_angles)
    nh_y = ny + circumradius * np.sin(hex_angles)
    ax.plot(nh_x, nh_y, 'k-', linewidth=1, alpha=0.3)
    # Draw neighbor hole
    nc = plt.Circle((nx, ny), hole_radius, color='lightblue', alpha=0.3)
    ax.add_patch(nc)
    nc_edge = plt.Circle((nx, ny), hole_radius, fill=False, edgecolor='blue', linewidth=1, alpha=0.5)
    ax.add_patch(nc_edge)

# Annotations
# Pitch arrow
ax.annotate('', xy=(pitch * np.cos(np.pi / 6), pitch * np.sin(np.pi / 6)),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle='<->', color='red', lw=2))
mid_x = pitch * np.cos(np.pi / 6) / 2
mid_y = pitch * np.sin(np.pi / 6) / 2
ax.text(mid_x - 2.5, mid_y + 0.5, f'pitch = {pitch} mm', color='red',
        fontsize=11, fontweight='bold', rotation=30)

# Diameter arrow
ax.annotate('', xy=(hole_radius, 0), xytext=(-hole_radius, 0),
            arrowprops=dict(arrowstyle='<->', color='darkblue', lw=2))
ax.text(-1.2, -1.2, f'Ø = {hole_diameter} mm', color='darkblue',
        fontsize=11, fontweight='bold')

# Wall thickness annotation
wall_thickness = pitch - hole_diameter  # approximate, between flat sides
ax.annotate('', xy=(hole_radius + 0.2, 4.5), xytext=(apothem - 0.2, 4.5),
            arrowprops=dict(arrowstyle='<->', color='brown', lw=1.5))
ax.text(hole_radius + 0.5, 5.0, f'wall ≈ {wall_thickness:.1f} mm', color='brown',
        fontsize=10)

# SiPM at center
ax.plot(0, 0, 'rs', markersize=8, label='SiPM (at bottom)')

ax.set_xlim(-14, 14)
ax.set_ylim(-14, 14)
ax.set_xlabel('x (mm)', fontsize=12)
ax.set_ylabel('y (mm)', fontsize=12)
ax.legend(loc='upper left', fontsize=9)
ax.grid(True, alpha=0.3)

# ============================================================
# RIGHT: Side cross-section
# ============================================================
ax = axes[1]
ax.set_title('Side Cross-Section — Single Cell', fontsize=14, fontweight='bold')

hole_height = 5.0    # mm
drift_length = 15.0   # mm (electron starts at 20 mm from bottom, hole top at 5 mm)
el_bottom = 0.0       # mm (bottom of hole = SiPM plane)
el_top = hole_height  # mm (top of hole)
drift_top = 20.0      # mm (electron start: 2 cm from bottom)

# Wall material (PTFE/copper) — left and right of hole
wall_width = 4.0  # mm to draw on each side
# Left wall
rect_left = patches.Rectangle((-hole_radius - wall_width, el_bottom),
                               wall_width, hole_height,
                               facecolor='#D2B48C', edgecolor='black', linewidth=2,
                               label='Wall (PTFE/Cu)')
ax.add_patch(rect_left)
# Right wall
rect_right = patches.Rectangle((hole_radius, el_bottom),
                                wall_width, hole_height,
                                facecolor='#D2B48C', edgecolor='black', linewidth=2)
ax.add_patch(rect_right)

# Anode plate above hole (extends beyond walls)
plate_width = hole_radius + wall_width
# Left plate extension above
rect_plate_left = patches.Rectangle((-plate_width, el_top),
                                     plate_width - hole_radius, 1.0,
                                     facecolor='#D2B48C', edgecolor='black', linewidth=1.5)
ax.add_patch(rect_plate_left)
# Right plate extension above
rect_plate_right = patches.Rectangle((hole_radius, el_top),
                                      plate_width - hole_radius, 1.0,
                                      facecolor='#D2B48C', edgecolor='black', linewidth=1.5)
ax.add_patch(rect_plate_right)

# Hole interior (Xe gas with EL field)
rect_hole = patches.Rectangle((-hole_radius, el_bottom),
                               hole_diameter, hole_height,
                               facecolor='#ADD8E6', edgecolor='blue', linewidth=2,
                               alpha=0.4, label=f'Hole (Xe, EL region)')
ax.add_patch(rect_hole)

# SiPM at bottom
rect_sipm = patches.Rectangle((-hole_radius, el_bottom - 0.8),
                                hole_diameter, 0.8,
                                facecolor='green', edgecolor='darkgreen', linewidth=2,
                                alpha=0.7, label='SiPM')
ax.add_patch(rect_sipm)

# Drift region above
rect_drift = patches.Rectangle((-plate_width, el_top + 1.0),
                                2 * plate_width, drift_top - el_top - 1.0,
                                facecolor='#E0E8FF', edgecolor='none', alpha=0.3,
                                label='Drift region (Xe)')
ax.add_patch(rect_drift)

# Electric field arrows inside hole (EL field, pointing down)
for y_arr in np.linspace(el_bottom + 0.8, el_top - 0.5, 4):
    ax.annotate('', xy=(0, y_arr - 0.6), xytext=(0, y_arr + 0.6),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))
ax.text(hole_radius + 0.5, hole_height / 2, 'E$_{EL}$\n4800\nV/cm',
        color='red', fontsize=10, fontweight='bold', ha='left', va='center')

# Drift field arrows (above hole)
for y_arr in [8, 11, 14, 17]:
    ax.annotate('', xy=(0, y_arr - 1.0), xytext=(0, y_arr + 1.0),
                arrowprops=dict(arrowstyle='->', color='orange', lw=1.5, alpha=0.7))
ax.text(3.5, 13, 'E$_{drift}$\n100 V/cm', color='orange', fontsize=10,
        fontweight='bold', ha='left', va='center')

# Electron trajectory (schematic)
e_x = [0.8, 0.6, 0.3, 0.1, -0.1, -0.2, 0.0, 0.1, -0.1]
e_z = [20.0, 17.0, 14.0, 11.0, 8.0, 5.0, 3.5, 2.0, 0.5]
ax.plot(e_x, e_z, 'b-', linewidth=2, alpha=0.7)
ax.plot(e_x[0], e_z[0], 'bo', markersize=10, label='e⁻ start (z=20 mm)')

# Excitation points (VUV photon emission) — only inside hole
exc_z = [4.5, 3.0, 1.5]
exc_x = [-0.15, 0.05, -0.05]
for ex, ez in zip(exc_x, exc_z):
    ax.plot(ex, ez, 'r*', markersize=14, zorder=5)
    # Photon rays
    for angle in [30, 120, 210, 300]:
        dx = 1.2 * np.cos(np.radians(angle))
        dz = 1.2 * np.sin(np.radians(angle))
        ax.annotate('', xy=(ex + dx, ez + dz), xytext=(ex, ez),
                    arrowprops=dict(arrowstyle='->', color='gold', lw=1, alpha=0.6))

ax.plot([], [], 'r*', markersize=12, label='Excitation → VUV γ')

# Dimension annotations
# Hole height
ax.annotate('', xy=(hole_radius + wall_width + 1.5, el_bottom),
            xytext=(hole_radius + wall_width + 1.5, el_top),
            arrowprops=dict(arrowstyle='<->', color='navy', lw=2))
ax.text(hole_radius + wall_width + 2.0, hole_height / 2, '5 mm\n(hole)',
        color='navy', fontsize=10, va='center')

# Drift distance
ax.annotate('', xy=(-plate_width - 1.5, el_top + 1.0),
            xytext=(-plate_width - 1.5, drift_top),
            arrowprops=dict(arrowstyle='<->', color='purple', lw=1.5))
ax.text(-plate_width - 1.3, (el_top + 1.0 + drift_top) / 2, '15 mm\n(drift)',
        color='purple', fontsize=9, ha='left', va='center')

# Total distance from bottom
ax.annotate('', xy=(plate_width + 3.5, el_bottom),
            xytext=(plate_width + 3.5, drift_top),
            arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
ax.text(plate_width + 4.0, drift_top / 2, '20 mm\n(e⁻ start)',
        color='black', fontsize=10, va='center')

# Hole diameter
ax.annotate('', xy=(-hole_radius, el_bottom - 1.5),
            xytext=(hole_radius, el_bottom - 1.5),
            arrowprops=dict(arrowstyle='<->', color='darkblue', lw=2))
ax.text(0, el_bottom - 2.5, f'Ø = {hole_diameter} mm', color='darkblue',
        fontsize=11, fontweight='bold', ha='center')

ax.set_xlim(-12, 16)
ax.set_ylim(-4, 23)
ax.set_xlabel('x (mm)', fontsize=12)
ax.set_ylabel('z (mm)', fontsize=12)
ax.legend(loc='upper left', fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/hernando/work/investigacion/NEXT/software/gala_sim/cell_diagram.png',
            dpi=150, bbox_inches='tight')
plt.show()
print("Diagram saved to cell_diagram.png")
