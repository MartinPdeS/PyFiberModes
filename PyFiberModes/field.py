"""Unit-aware mode fields and integrated modal observables.

Analytical fields have an arbitrary overall amplitude. Attaching V/m and A/m
expresses the solver's electric/magnetic convention; absolute optical power
requires choosing an amplitude normalization.
"""

import math
from types import SimpleNamespace

from PyFiberModes import _field
from PyFiberModes.coordinates import CartesianCoordinates, CylindricalCoordinates  # noqa: F401
from PyFiberModes.units import Angle, Length, quantity, scalar, ureg


def radial_quantities(vectors):
    """Attach field units to cylindrical component containers (not positions)."""
    return tuple(
        SimpleNamespace(
            **{name: quantity(getattr(vector, name), unit) for name in ("rho", "phi", "z")}
        )
        for vector, unit in zip(vectors, ("volt / meter", "ampere / meter"), strict=True)
    )


class Field:
    """Sample an optical mode over a grid with a length-valued half-width."""

    def __init__(self, fiber, mode, limit: Length, n_point: int = 101):
        """Convert the sampling limit to SI and initialize numerical fields."""
        if n_point < 2:
            raise ValueError("n_point must be at least two")
        limit_si = scalar(limit, "meter", "limit")
        if not math.isfinite(limit_si) or limit_si <= 0:
            raise ValueError("limit must be positive and finite")
        self.fiber = fiber
        self.mode = mode
        self._si = _field.Field(fiber=fiber._si, mode=mode, limit=limit_si, n_point=n_point)

    @property
    def limit(self):
        """Grid half-width in meters."""
        return quantity(self._si.limit, "meter")

    @property
    def n_point(self):
        """Number of grid points along each axis."""
        return self._si.n_point

    @property
    def cartesian_coordinates(self):
        """Sample positions as length-valued Cartesian coordinates."""
        return CartesianCoordinates._from_si(self._si.cartesian_coordinates)

    @property
    def cylindrical_coordinates(self):
        """Sample positions in cylindrical coordinates."""
        return self.cartesian_coordinates.to_cylindrical()

    def get_components(
        self,
        components=("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"),
        phi: Angle = 0 * ureg.radian,
        theta: Angle = 0 * ureg.radian,
    ):
        """Return electric and magnetic component arrays with physical units."""
        values = self._si.get_components(
            components=components,
            phi=scalar(phi, "radian", "phi"),
            theta=scalar(theta, "radian", "theta"),
        )
        return {
            name: quantity(value, "volt / meter" if name.startswith("E") else "ampere / meter")
            for name, value in values.items()
        }

    def get_azimuthal_dependency(self, phi: Angle, dependency_type: str):
        """Evaluate dimensionless sine/cosine angular dependencies."""
        return quantity(
            self._si.get_azimuthal_dependency(scalar(phi, "radian", "phi"), dependency_type),
            "dimensionless",
        )

    def get_effective_area(self):
        """Return effective mode area in square meters."""
        return quantity(self._si.get_effective_area(), "meter ** 2")

    def get_integrale_square(self, array):
        """Integrate a quantity-valued squared field over the sample area."""
        if not hasattr(array, "units"):
            raise TypeError("array must be a pint Quantity")
        return quantity(
            self._si.get_integrale_square(array.magnitude), array.units**2 * ureg.meter**2
        )

    def get_intensity(self):
        """Return the solver's integrated electric intensity in square volts."""
        return quantity(self._si.get_intensity(), "volt ** 2")

    def get_normalization_constant(self):
        """Return the field's power normalization constant in watts."""
        return quantity(self._si.get_normalization_constant(), "watt")

    def get_poynting_vector(self, phi: Angle = 0 * ureg.radian, theta: Angle = 0 * ureg.radian):
        """Return Cartesian power flux components in watts per square meter."""
        return tuple(
            quantity(value, "watt / meter ** 2")
            for value in self._si.get_poynting_vector(
                phi=scalar(phi, "radian", "phi"),
                theta=scalar(theta, "radian", "theta"),
            )
        )

    def get_power(self, phi: Angle = 0 * ureg.radian, theta: Angle = 0 * ureg.radian):
        """Return integrated longitudinal power in watts."""
        return quantity(
            self._si.get_power(
                phi=scalar(phi, "radian", "phi"),
                theta=scalar(theta, "radian", "theta"),
            ),
            "watt",
        )

    def get_confinement_factor(self, radius: Length):
        """Return the dimensionless fraction of electric energy within radius."""
        return quantity(
            self._si.get_confinement_factor(scalar(radius, "meter", "radius")), "dimensionless"
        )

    def overlap(self, other):
        """Return the dimensionless normalized electric-field overlap."""
        from PyFiberModes.propagation import overlap

        return overlap(self, other)

    def plot(self, plot_type=(), show=True, save_filename=None):
        """Plot selected numerical field components using the optional backend."""
        return self._si.plot(plot_type=plot_type, show=show, save_filename=save_filename)


def _component_method(name):
    """Wrap a numerical field component with angle checks and output units."""
    original = getattr(_field.Field, name)

    def method(self, phi: Angle = 0 * ureg.radian, theta: Angle = 0 * ureg.radian):
        """Evaluate a field component at quantity-valued phase/rotation angles."""
        value = original(
            self._si, phi=scalar(phi, "radian", "phi"), theta=scalar(theta, "radian", "theta")
        )
        unit = (
            "radian"
            if name.endswith("pol")
            else "volt / meter"
            if name.startswith("E")
            else "ampere / meter"
        )
        return quantity(value, unit)

    method.__name__ = name
    return method


for _name in (
    "Ex",
    "Ey",
    "Ez",
    "Er",
    "Ephi",
    "Et",
    "Epol",
    "Emod",
    "Hx",
    "Hy",
    "Hz",
    "Hr",
    "Hphi",
    "Ht",
    "Hpol",
    "Hmod",
):
    setattr(Field, _name, _component_method(_name))
