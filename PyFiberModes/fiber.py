"""Unit-checked public fiber model backed by numerical SI solvers."""

import math
from dataclasses import replace

from PyFiberModes import _fiber
from PyFiberModes.materials import DEFAULT_MATERIALS, MaterialRegistry
from PyFiberModes.units import (
    Length,
    RefractiveIndex,
    Quantity,
    MODAL_UNITS,
    magnitude,
    quantity,
    scalar,
)


class Layer:
    """Unit-aware view of a numerical step-index layer.

    Radius and index assignments update the owning fiber and invalidate caches.
    """

    def __init__(self, fiber, layer):
        """Retain the owning fiber and its numerical layer."""
        self._fiber = fiber
        self._layer = layer

    @property
    def name(self):
        """Name assigned when the layer was added to the fiber."""
        return self._fiber.layer_names[self._fiber._si.layers.index(self._layer)]

    @property
    def position(self):
        """Zero-based radial order of this layer."""
        return self._fiber._si.layers.index(self._layer)

    @property
    def radius_in(self) -> Length:
        """Inner radius as a length quantity."""
        return quantity(self._layer.radius_in, "meter")

    @radius_in.setter
    def radius_in(self, value: Length):
        """Update the inner radius in meters."""
        radius = scalar(value, "meter", "radius_in")
        if not math.isfinite(radius) or radius < 0:
            raise ValueError("radius_in must be non-negative and finite")
        self._layer.radius_in = radius
        self._fiber.clear_caches()

    @property
    def radius_out(self) -> Length:
        """Outer radius as a length quantity."""
        return quantity(self._layer.radius_out, "meter")

    @radius_out.setter
    def radius_out(self, value: Length):
        """Update the outer radius in meters."""
        radius = scalar(value, "meter", "radius_out")
        if math.isnan(radius) or radius <= 0:
            raise ValueError("radius_out must be positive")
        self._layer.radius_out = radius
        self._fiber.clear_caches()

    @property
    def refractive_index(self) -> RefractiveIndex:
        """Dimensionless refractive index."""
        return quantity(self._layer.refractive_index, "dimensionless")

    @refractive_index.setter
    def refractive_index(self, value: RefractiveIndex):
        """Update the refractive index."""
        import numpy as np

        index = scalar(value, "dimensionless", "refractive_index")
        if not np.isfinite(index) or index <= 0:
            raise ValueError("refractive index must be positive and finite")
        self._layer.refractive_index = index
        self._fiber.clear_caches()

    @property
    def wavelength(self) -> Length:
        """Wavelength shared with the owning fiber."""
        return self._fiber.wavelength

    @property
    def is_first_layer(self):
        """Whether this is the innermost layer."""
        return self._layer.is_first_layer

    @property
    def is_last_layer(self):
        """Whether this is the unbounded outer layer."""
        return self._layer.is_last_layer


class Fiber:
    """Multilayer fiber with explicit quantities for all physical parameters.

    Parameters
    ----------
    wavelength : TypedUnit.Length
        Positive finite vacuum wavelength, for example ``1550 * ureg.nanometer``.
    """

    def __init__(self, wavelength: Length):
        """Validate wavelength and initialize an empty numerical fiber."""
        self._si = _fiber.Fiber(wavelength=scalar(wavelength, "meter", "wavelength"))

    @classmethod
    def _from_si(cls, fiber):
        """Wrap an existing SI fiber without reinterpreting numerical values."""
        instance = cls.__new__(cls)
        instance._si = fiber
        return instance

    @property
    def analysis(self):
        """Return modal analysis with unit-aware results and solver settings."""
        return ModalAnalysis(self)

    @property
    def fields(self):
        """Return a unit-aware field-construction service."""
        return FieldAnalysis(self)

    @property
    def geometry_signature(self):
        """Return an immutable numerical signature used for cache diagnostics."""
        return self._si.geometry_signature

    @property
    def wavelength(self) -> Length:
        """Vacuum wavelength in meters."""
        return quantity(self._si.wavelength, "meter")

    @wavelength.setter
    def wavelength(self, value: Length):
        """Update wavelength consistently across all layers."""
        self.update_wavelength(value)

    def update_wavelength(self, wavelength: Length) -> None:
        """Validate and update the wavelength, clearing numerical caches."""
        import numpy as np

        value = scalar(wavelength, "meter", "wavelength")
        if not np.isfinite(value) or value <= 0:
            raise ValueError("wavelength must be positive and finite")
        self._si.update_wavelength(value)

    def add_layer(self, name: str, radius: Length, index: RefractiveIndex) -> None:
        """Add a layer with a length radius and a dimensionless index quantity."""
        import numpy as np

        radius_si = scalar(radius, "meter", "radius")
        if np.isnan(radius_si) or radius_si <= 0:
            raise ValueError(
                "radius must be positive; use infinity * ureg.meter for the outer layer"
            )
        self._si.add_layer(
            name=name,
            radius=radius_si,
            index=scalar(index, "dimensionless", "index"),
        )

    def initialize_layers(self) -> None:
        """Mark inner and outer layers and finalize solver geometry."""
        self._si.initialize_layers()

    def clear_caches(self) -> None:
        """Invalidate cached numerical solutions."""
        self._si.clear_caches()

    @property
    def layers(self):
        """Ordered unit-aware views of the step-index layers."""
        return [Layer(self, layer) for layer in self._si.layers]

    def __getitem__(self, index):
        """Return a layer by its radial order."""
        return self.layers[index]

    @property
    def first_layer(self):
        """Innermost fiber layer."""
        return self[0]

    @property
    def last_layer(self):
        """Unbounded outermost fiber layer."""
        return self[-1]

    @property
    def penultimate_layer(self):
        """Layer immediately inside the outermost layer."""
        return self[-2]

    @property
    def n_layer(self):
        """Number of layers."""
        return self._si.n_layer

    @property
    def n_interface(self):
        """Number of interfaces."""
        return self._si.n_interface

    @property
    def layer_names(self):
        """Names of the fiber layers."""
        return self._si.layer_names

    @property
    def layer_radius(self):
        """Layer radii as a length array."""
        return quantity([layer.radius_out for layer in self._si.layers], "meter")

    @property
    def index_list(self):
        """Layer refractive indices as a dimensionless array."""
        return quantity([layer.refractive_index for layer in self._si.layers], "dimensionless")

    @property
    def radius(self) -> Length:
        """Largest finite interface radius."""
        return quantity(self._si.radius, "meter")

    @property
    def maximum_index(self) -> RefractiveIndex:
        """Largest material refractive index."""
        return quantity(self._si.maximum_index, "dimensionless")

    @property
    def minimum_index(self) -> RefractiveIndex:
        """Smallest material refractive index."""
        return quantity(self._si.minimum_index, "dimensionless")

    @property
    def V_number(self):
        """Dimensionless normalized frequency."""
        return quantity(self._si.V_number, "dimensionless")

    @property
    def M_number(self):
        """Dimensionless approximate guided-mode count."""
        return quantity(self._si.M_number, "dimensionless")

    def get_NA(self):
        """Return the dimensionless numerical aperture."""
        return quantity(self._si.get_NA(), "dimensionless")

    def scale(self, factor: float):
        """Return an independent fiber scaled by a dimensionless factor."""
        return self._from_si(self._si.scale(factor=float(factor)))

    def iterate_interfaces(self):
        """Yield adjacent unit-aware layer pairs."""
        from itertools import pairwise

        yield from pairwise(self.layers)

    def get_layer_at_radius(self, radius: Length):
        """Find the layer containing a length-valued radial position."""
        layer = self._si.get_layer_at_radius(scalar(radius, "meter", "radius"))
        return Layer(self, layer)

    def get_index_at_radius(self, radius: Length):
        """Return the dimensionless refractive index at a radial position."""
        return self.get_layer_at_radius(radius).refractive_index

    def get_index_profile(self, radius: Length):
        """Evaluate a dimensionless index profile at length-valued positions."""
        return quantity(
            self._si.get_index_profile(magnitude(radius, "meter", "radius")), "dimensionless"
        )

    def solve_effective_index(self, mode):
        """Return solver diagnostics with a quantity-valued effective index."""
        result = self._si.solve_effective_index(mode)
        return replace(
            result, value=None if result.value is None else quantity(result.value, "dimensionless")
        )

    def get_mode_field(self, mode, limit: Length = None, n_point: int = 101):
        """Sample a field on a grid with an optional length-valued half-width."""
        from PyFiberModes.field import Field

        return Field(
            fiber=self,
            mode=mode,
            limit=self.radius * 5.5 if limit is None else limit,
            n_point=n_point,
        )

    def get_radial_field(self, mode, radius: Length):
        """Return electric (V/m) and magnetic (A/m) cylindrical components."""
        from PyFiberModes.field import radial_quantities

        return radial_quantities(self._si.get_radial_field(mode, scalar(radius, "meter", "radius")))

    def get_radial_fields(self, mode, radius: Length):
        """Return electric and magnetic components over a length-valued array."""
        from PyFiberModes.field import radial_quantities

        return radial_quantities(
            self._si.get_radial_fields(mode, magnitude(radius, "meter", "radius"))
        )

    def get_radial_field_norm(self, mode, radius: Length):
        """Return the electric-field norm in volts per meter."""
        return quantity(
            self._si.get_radial_field_norm(mode, scalar(radius, "meter", "radius")), "volt / meter"
        )

    def does_mode_exist(self, *mode_list):
        """Return existence flags for the supplied modes."""
        return self._si.does_mode_exist(*mode_list)

    def find_modes(self, families=("LP",), max_nu=6, max_m=6):
        """Discover guided modes within the requested search bounds."""
        return self._si.find_modes(families=families, max_nu=max_nu, max_m=max_m)

    def print_data(self, data_type_list: list[str], mode_list: list):
        """Print modal quantities including their output units."""
        for data_type in data_type_list:
            print(f"{data_type} @ wavelength: {self.wavelength}")
            for mode in mode_list:
                from PyFiberModes.exceptions import ConvergenceError

                try:
                    value = getattr(self, "get_" + data_type.lower())(mode=mode)
                except ConvergenceError:
                    value = quantity(float("nan"), MODAL_UNITS[data_type.lower()])
                print(f"{mode}: {value}")

    def sweep(self, parameters, **kwargs):
        """Sweep length-valued wavelengths or a custom quantity parameter."""
        from PyFiberModes.analysis import sweep_modes

        return sweep_modes(self, parameters, **kwargs)


def _modal_method(name, unit):
    """Expose a numerical modal method with its physical output unit."""
    original = getattr(_fiber.Fiber, f"get_{name}")

    def method(self, mode) -> Quantity:
        """Calculate a modal quantity in its documented output units."""
        return quantity(original(self._si, mode), unit)

    method.__name__ = original.__name__
    method.__qualname__ = f"Fiber.{original.__name__}"
    method.__doc__ = (
        original.__doc__.replace("float\n", "pint.Quantity\n") + f"\n        Output unit: {unit}.\n"
    )
    return method


for _name, _unit in MODAL_UNITS.items():
    setattr(Fiber, f"get_{_name}", _modal_method(_name, _unit))


def load_fiber(
    fiber_name: str,
    wavelength: Length,
    add_air_layer: bool = False,
    materials: MaterialRegistry = DEFAULT_MATERIALS,
) -> Fiber:
    """Load YAML geometry (stored in SI) at a length-valued vacuum wavelength."""
    value = scalar(wavelength, "meter", "wavelength")
    return Fiber._from_si(
        _fiber.load_fiber(
            fiber_name=fiber_name,
            wavelength=value,
            add_air_layer=add_air_layer,
            materials=materials,
        )
    )


def get_fiber_from_delta_and_V0(delta, V0, wavelength: Length) -> Fiber:
    """Construct a fiber from dimensionless contrast, V-number and wavelength."""
    return Fiber._from_si(
        _fiber.get_fiber_from_delta_and_V0(
            delta=float(delta),
            V0=float(V0),
            wavelength=scalar(wavelength, "meter", "wavelength"),
        )
    )


class ModalAnalysis:
    """Unit-aware view of a fiber's cached numerical modal-analysis service."""

    def __init__(self, fiber):
        """Retain the owning fiber."""
        self.fiber = fiber

    @property
    def settings(self):
        """Numerical solver settings."""
        return self.fiber._si.analysis.settings

    @settings.setter
    def settings(self, value):
        """Update numerical solver settings and invalidate cached results."""
        self.fiber._si.analysis.settings = value
        self.clear()

    @property
    def cache_size(self):
        """Number of cached modal results."""
        return self.fiber._si.analysis.cache_size

    def clear(self):
        """Discard cached modal results."""
        self.fiber._si.analysis.clear()

    def effective_index(self, mode):
        """Return a dimensionless effective-index quantity."""
        return self.fiber.get_effective_index(mode)

    def effective_index_result(self, mode):
        """Return unit-bearing effective-index solver diagnostics."""
        return self.fiber.solve_effective_index(mode)

    def cutoff_result(self, mode):
        """Return unit-bearing normalized-cutoff solver diagnostics."""
        result = self.fiber._si.analysis.cutoff_result(mode)
        return replace(
            result, value=None if result.value is None else quantity(result.value, "dimensionless")
        )


class FieldAnalysis:
    """Unit-aware field-construction service."""

    def __init__(self, fiber):
        """Retain the owning fiber."""
        self.fiber = fiber

    def mode_field(self, mode, limit: Length = None, n_point: int = 101):
        """Construct a field with a length-valued sampling limit."""
        return self.fiber.get_mode_field(mode=mode, limit=limit, n_point=n_point)
