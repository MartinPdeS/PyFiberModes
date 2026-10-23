"""Length-valued Cartesian and cylindrical coordinate grids."""

from dataclasses import dataclass
import numpy as np

from PyFiberModes import _coordinates
from PyFiberModes.units import Angle, Length, magnitude, quantity


@dataclass
class CylindricalCoordinates:
    """Cylindrical positions with length-valued rho/z and angle-valued phi."""

    rho: Length
    phi: Angle
    z: Length

    def __post_init__(self):
        """Validate and normalize coordinate quantities."""
        self.rho = quantity(magnitude(self.rho, "meter", "rho"), "meter")
        self.phi = quantity(magnitude(self.phi, "radian", "phi"), "radian")
        self.z = quantity(magnitude(self.z, "meter", "z"), "meter")

    def to_cartesian(self):
        """Convert cylindrical positions to Cartesian positions."""
        return CartesianCoordinates(
            x=self.rho * np.cos(self.phi),
            y=self.rho * np.sin(self.phi),
            z=self.z,
        )

    def to_cylindrical(self):
        """Return these cylindrical positions."""
        return self


@dataclass
class CartesianCoordinates:
    """Cartesian positions with explicit length units on every axis."""

    x: Length
    y: Length
    z: Length
    is_structured: bool = False
    is_3D: bool = False

    def __post_init__(self):
        """Validate lengths and normalize numerical coordinate arrays to meters."""
        for name in ("x", "y", "z"):
            setattr(
                self,
                name,
                quantity(
                    np.asarray(magnitude(getattr(self, name), "meter", name), dtype=float), "meter"
                ),
            )

    @classmethod
    def _from_si(cls, coordinates):
        """Attach length units to a numerical coordinate grid."""
        return cls(
            x=quantity(coordinates.x, "meter"),
            y=quantity(coordinates.y, "meter"),
            z=quantity(coordinates.z, "meter"),
            is_structured=coordinates.is_structured,
            is_3D=coordinates.is_3D,
        )

    @property
    def x_boundaries(self):
        """Minimum and maximum x positions."""
        return self.x.min(), self.x.max()

    @property
    def y_boundaries(self):
        """Minimum and maximum y positions."""
        return self.y.min(), self.y.max()

    @property
    def dx(self):
        """Structured grid spacing along x."""
        if not self.is_structured:
            raise ValueError("dx value cannot be inferred from an unstructured mesh")
        return quantity(np.diff(np.unique(self.x.magnitude))[0], self.x.units)

    @property
    def dy(self):
        """Structured grid spacing along y."""
        if not self.is_structured:
            raise ValueError("dy value cannot be inferred from an unstructured mesh")
        return quantity(np.diff(np.unique(self.y.magnitude))[0], self.y.units)

    @classmethod
    def generate_from_square(cls, length: Length, center=None, n_points: int = 100):
        """Generate a square grid with a length-valued side and center."""
        center_si = (
            (0.0, 0.0)
            if center is None
            else tuple(magnitude(item, "meter", "center") for item in center)
        )
        return cls._from_si(
            _coordinates.CartesianCoordinates.generate_from_square(
                length=magnitude(length, "meter", "length"),
                center=center_si,
                n_points=n_points,
            )
        )

    @classmethod
    def generate_from_boundaries(cls, x_limits, y_limits, z_limits, x_points, y_points, z_points):
        """Generate a structured grid from length-valued axis limits."""
        limits = [
            tuple(magnitude(item, "meter", "limits") for item in axis)
            for axis in (x_limits, y_limits, z_limits)
        ]
        return cls._from_si(
            _coordinates.CartesianCoordinates.generate_from_boundaries(
                x_limits=limits[0],
                y_limits=limits[1],
                z_limits=limits[2],
                x_points=x_points,
                y_points=y_points,
                z_points=z_points,
            )
        )

    def to_cylindrical(self):
        """Convert positions to cylindrical coordinates."""
        return CylindricalCoordinates(
            rho=np.sqrt(self.x**2 + self.y**2),
            phi=np.arctan2(self.y, self.x),
            z=self.z,
        )

    def centering(self):
        """Center each axis in place."""
        self.x -= self.x.mean()
        self.y -= self.y.mean()
        self.z -= self.z.mean()

    def shift_coordinates(self, shift):
        """Shift the grid by three length-valued displacements."""
        # Validate all inputs before mutating the grid.
        values = [quantity(magnitude(item, "meter", "shift"), "meter") for item in shift]
        self.x += values[0]
        self.y += values[1]
        self.z += values[2]

    def scale_coordinates(self, scale):
        """Scale axes by dimensionless factors."""
        self.x *= float(scale[0])
        self.y *= float(scale[1])
        self.z *= float(scale[2])
