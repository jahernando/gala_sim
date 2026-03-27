"""
COMSOL Field Map Interpolator
==============================
Loads the MTHGEM potential map and provides E-field at any (x', y, z') point.

Coordinate system (hole-centered):
  - Central hole at origin (0, 0) in (x', z') plane
  - y = drift direction (electrons drift from high y toward y = 0)
  - COMSOL rectangle: x' in [-8.70, 0], z' in [0, 5.00]
  - Sector covered: 150° to 180° (30° unique sector of hexagonal symmetry)

Usage:
    from field_map import FieldMap
    fm = FieldMap("data/Mesh-MTHGEM_fields_4bar.txt")
    Ex, Ey, Ez = fm.get_E(x, y, z)  # all in V/cm, x/y/z in cm
"""

import numpy as np
from scipy.interpolate import griddata, RegularGridInterpolator
import os


class FieldMap:
    """3D electric field interpolator from COMSOL potential data.

    Builds a regular 3D grid of V, then computes E = -grad(V) on that grid.
    Lookups via RegularGridInterpolator are O(1) per point.
    """

    def __init__(self, filepath=None, y_transition_cm=0.5,
                 nx=80, ny=100, nz=50):
        """
        Parameters
        ----------
        filepath : str
            Path to COMSOL data file.
        y_transition_cm : float
            Above this y (in cm), use uniform drift field.
        nx, ny, nz : int
            Grid resolution for the regular interpolation grid.
        """
        if filepath is None:
            here = os.path.dirname(os.path.abspath(__file__))
            filepath = os.path.join(here, "data", "Mesh-MTHGEM_fields_4bar.txt")

        self.y_transition_cm = y_transition_cm
        self._load(filepath)
        self._build_grid(nx, ny, nz)

    def _load(self, filepath):
        """Load COMSOL data and shift to hole-centered coordinates."""
        data = np.loadtxt(filepath, skiprows=9)
        xr, yr, zr, V = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

        # Central hole at (x_max, z_min) → origin
        x0, z0 = xr.max(), zr.min()
        self._x_mm = xr - x0
        self._z_mm = zr - z0
        self._y_mm = yr
        self._V = V

        # Domain bounds in mm
        self.dom_x_min = self._x_mm.min()
        self.dom_x_max = self._x_mm.max()
        self.dom_z_min = self._z_mm.min()
        self.dom_z_max = self._z_mm.max()

        # Drift field from linear region
        mask_drift = self._y_mm > 50
        if mask_drift.sum() > 100:
            coeffs = np.polyfit(self._y_mm[mask_drift], self._V[mask_drift], 1)
            self.drift_Ey_V_per_mm = -coeffs[0]
        else:
            self.drift_Ey_V_per_mm = 1.05

        self.drift_Ey_V_per_cm = self.drift_Ey_V_per_mm * 10.0

        print(f"FieldMap loaded: {len(self._V)} nodes")
        print(f"  Domain (mm): x'=[{self.dom_x_min:.2f}, {self.dom_x_max:.2f}], "
              f"z'=[{self.dom_z_min:.2f}, {self.dom_z_max:.2f}]")
        print(f"  Drift Ey = {self.drift_Ey_V_per_cm:.1f} V/cm")

    def _build_grid(self, nx, ny, nz):
        """Build regular 3D grid of V and E, layer by layer (2D per y-slab)."""
        y_max_mm = self.y_transition_cm * 10.0

        self._xg = np.linspace(self.dom_x_min, self.dom_x_max, nx)
        self._yg = np.linspace(-0.5, y_max_mm, ny)
        self._zg = np.linspace(self.dom_z_min, self.dom_z_max, nz)

        XG, ZG = np.meshgrid(self._xg, self._zg, indexing='ij')
        query_xz = np.column_stack([XG.ravel(), ZG.ravel()])

        V_grid = np.zeros((nx, ny, nz))

        print(f"  Building {nx}x{ny}x{nz} grid (layer-by-layer)...", flush=True)
        for j, yv in enumerate(self._yg):
            y_tol = max(0.3, (self._yg[1] - self._yg[0]) * 1.5)
            mask = np.abs(self._y_mm - yv) < y_tol
            if mask.sum() < 20:
                mask = np.abs(self._y_mm - yv) < 2.0
            if mask.sum() < 5:
                V_grid[:, j, :] = V_grid[:, max(0, j-1), :]
                continue

            pts2d = np.column_stack([self._x_mm[mask], self._z_mm[mask]])
            vals = self._V[mask]

            Vi = griddata(pts2d, vals, query_xz, method='linear')
            nans = np.isnan(Vi)
            if nans.any():
                Vi[nans] = griddata(pts2d, vals, query_xz[nans], method='nearest')
            V_grid[:, j, :] = Vi.reshape(nx, nz)

        # Compute E = -grad(V), convert V/mm → V/cm
        dx = self._xg[1] - self._xg[0]
        dy = self._yg[1] - self._yg[0]
        dz = self._zg[1] - self._zg[0]

        Ex_grid = -np.gradient(V_grid, dx, axis=0) * 10.0
        Ey_grid = -np.gradient(V_grid, dy, axis=1) * 10.0
        Ez_grid = -np.gradient(V_grid, dz, axis=2) * 10.0

        axes = (self._xg, self._yg, self._zg)
        self._interp_Ex = RegularGridInterpolator(axes, Ex_grid,
                                                   bounds_error=False, fill_value=0.0)
        self._interp_Ey = RegularGridInterpolator(axes, Ey_grid,
                                                   bounds_error=False,
                                                   fill_value=self.drift_Ey_V_per_cm)
        self._interp_Ez = RegularGridInterpolator(axes, Ez_grid,
                                                   bounds_error=False, fill_value=0.0)
        print("  Grid ready.")

    def fold_to_sector(self, x_mm, z_mm):
        """Fold any (x', z') into the COMSOL 150°-180° sector.

        Uses hexagonal 6-fold symmetry + mirror to reduce any angle
        to the unique 30° sector, then maps to 150°-180°.

        Parameters & returns in mm.
        """
        x_mm = np.asarray(x_mm, dtype=float)
        z_mm = np.asarray(z_mm, dtype=float)

        r = np.sqrt(x_mm**2 + z_mm**2)
        theta = np.arctan2(z_mm, x_mm) % (2 * np.pi)

        # Reduce to [0, 60°) by 6-fold symmetry
        sector = theta % (np.pi / 3)

        # Mirror at 30°: fold [30°, 60°) → [0°, 30°]
        sector = np.where(sector > np.pi / 6, np.pi / 3 - sector, sector)

        # Map to COMSOL sector: phi = 180° - sector
        phi = np.pi - sector

        xf = r * np.cos(phi)
        zf = r * np.sin(phi)

        # Clip to domain
        xf = np.clip(xf, self.dom_x_min, self.dom_x_max)
        zf = np.clip(zf, self.dom_z_min, self.dom_z_max)

        return xf, zf

    def get_E(self, x_cm, y_cm, z_cm):
        """Compute E-field at position (x', y, z') in hole-centered coords.

        Parameters
        ----------
        x_cm, y_cm, z_cm : float or array
            Position in cm (hole-centered coordinates).

        Returns
        -------
        Ex, Ey, Ez : float or array
            Electric field components in V/cm.
            Ey > 0 means field points toward -y (toward anode).
        """
        scalar = np.ndim(x_cm) == 0
        x_cm = np.atleast_1d(np.asarray(x_cm, dtype=float))
        y_cm = np.atleast_1d(np.asarray(y_cm, dtype=float))
        z_cm = np.atleast_1d(np.asarray(z_cm, dtype=float))

        Ex = np.zeros_like(x_cm)
        Ey = np.full_like(y_cm, self.drift_Ey_V_per_cm)
        Ez = np.zeros_like(z_cm)

        # Points in the interpolated region
        near = y_cm < self.y_transition_cm
        if not near.any():
            return _squeeze(Ex, Ey, Ez) if scalar else (Ex, Ey, Ez)

        # Convert to mm and fold
        x_mm = x_cm[near] * 10.0
        y_mm = y_cm[near] * 10.0
        z_mm = z_cm[near] * 10.0
        xf, zf = self.fold_to_sector(x_mm, z_mm)

        # Query the pre-built interpolators (fast O(1) lookups)
        query = np.column_stack([xf, y_mm, zf])
        Ex[near] = self._interp_Ex(query)
        Ey[near] = self._interp_Ey(query)
        Ez[near] = self._interp_Ez(query)

        if scalar:
            return _squeeze(Ex, Ey, Ez)
        return Ex, Ey, Ez


def _squeeze(Ex, Ey, Ez):
    """Return scalars if input was scalar."""
    if Ex.size == 1:
        return float(Ex[0]), float(Ey[0]), float(Ez[0])
    return Ex, Ey, Ez
