"""Core fiber model and modal-property calculations."""

import numpy
from typing import Generator
from functools import cache
from dataclasses import dataclass, field
from scipy import constants
from copy import deepcopy
from itertools import pairwise
from scipy.constants import c

from PyFiberModes.stepindex import StepIndex
from PyFiberModes import Mode
from PyFinitDiff.finite_difference_1D import get_function_derivative
from PyFiberModes.field import Field

from PyFiberModes.fundamentals import (
    get_radial_field,
    get_propagation_constant_from_omega
)

from PyFiberModes import loader
from PyFiberModes.coordinates import CylindricalCoordinates
from PyFiberModes.exceptions import ConvergenceError, ValidationError
from PyFiberModes.materials import DEFAULT_MATERIALS, MaterialRegistry
from PyFiberModes.services import FieldAnalysis, ModalAnalysis


@dataclass
class Fiber(object):
    """
    A class to represent an optical fiber with multiple layers.

    This class provides functionality to define and manipulate
    the properties of optical fibers, including their layers,
    refractive indices, and propagation characteristics.

    Parameters
    ----------
    wavelength : float
        The wavelength of light used for fiber simulations.
    layer_names : list of str
        Names of the fiber layers.
    layer_radius : list of float
        Radii of the fiber layers in meters.
    layer_types : list of str
        Types of layers (e.g., core, cladding).
    index_list : list of float
        Refractive indices for each layer.
    """
    wavelength: float
    layer_names: list = field(default_factory=list)
    layer_radius: list = field(default_factory=list)
    layer_types: list = field(default_factory=list)
    index_list: list = field(default_factory=list)

    def __post_init__(self):
        """Initialize mutable layer and solver state after construction."""
        if self.wavelength is not None and (
            not numpy.isfinite(self.wavelength) or self.wavelength <= 0
        ):
            raise ValidationError("wavelength must be positive and finite")
        self.layers_parameters = []
        self.radius_in = 0
        self.layers = []
        self._analysis_service = None
        self._field_service = None

    @property
    def geometry_signature(self) -> tuple:
        """Return an immutable signature of all solver-relevant fiber state.

        Returns
        -------
        tuple
            Wavelength and ordered layer geometry/index values.
        """
        layers = tuple(
            (layer.radius_in, layer.radius_out, layer.refractive_index)
            for layer in self.layers
        )
        return self.wavelength, layers

    @property
    def analysis(self) -> ModalAnalysis:
        """Return the focused modal-analysis service.

        Returns
        -------
        ModalAnalysis
            Geometry-aware cached analysis service.
        """
        if self._analysis_service is None:
            self._analysis_service = ModalAnalysis(self)
        return self._analysis_service

    @property
    def fields(self) -> FieldAnalysis:
        """Return the focused field-construction service.

        Returns
        -------
        FieldAnalysis
            Field orchestration service.
        """
        if self._field_service is None:
            self._field_service = FieldAnalysis(self)
        return self._field_service

    def clear_caches(self) -> None:
        """Clear modal and radial-field caches after a physical mutation."""
        if self._analysis_service is not None:
            self._analysis_service.clear()
        self.get_radial_field.cache_clear()

    def scale(self, factor: float) -> None:
        """
        Scale the fiber's geometry by a given factor.

        This method scales all radii in the fiber by the specified factor.

        Parameters
        ----------
        factor : float
            Scaling factor to apply to the radii.

        Returns
        -------
        Fiber
            A scaled copy of the fiber.

        Notes
        -----
        This operation does not modify the original fiber.
        """
        scaled_fiber = deepcopy(self)
        for layer in scaled_fiber.layers:
            layer.radius_in *= factor
            layer.radius_out *= factor

        return scaled_fiber

    @property
    def n_layer(self) -> int:
        """
        Number of layers in the fiber.

        Returns
        -------
        int
            The total number of layers.
        """
        return len(self.layers)

    @property
    def n_interface(self) -> int:
        """
        Number of interfaces in the fiber.

        Returns
        -------
        int
            The number of interfaces, calculated as `n_layer - 1`.
        """
        return len(self.layers) - 1

    @property
    def last_layer(self) -> StepIndex:
        """Returns the last layer

        Returns
        -------
        StepIndex
            The last layer.
        """
        return self.layers[-1]

    @property
    def penultimate_layer(self) -> StepIndex:
        """Returns the second to last layer

        Returns
        -------
        StepIndex
            The second to last layer.
        """
        return self.layers[-2]

    @property
    def first_layer(self) -> StepIndex:
        """Returns the first layer

        Returns
        -------
        StepIndex
            The first layer.
        """
        return self.layers[0]

    def __hash__(self):
        """Return a hash derived from the ordered fiber layers.

        Returns
        -------
        int
            Hash of the layer tuple.
        """
        return hash(tuple(self.layers))

    def __getitem__(self, index: int) -> StepIndex:
        """Returns the nth layer.

        Parameters
        ----------
        index : int
            The index

        Returns
        -------
        StepIndex
            The step index.
        """
        return self.layers[index]

    def iterate_interfaces(self) -> Generator:
        """Iterates through pair of layers that forms interfaces

        Yields
        ------
        tuple[StepIndex, StepIndex]
            The two layers that form the interfaces.
        """
        for layer_in, layer_out in pairwise(self.layers):
            yield layer_in, layer_out

    def update_wavelength(self, wavelength: float) -> None:
        """Update the wavelength of the fiber and all its layers

        Parameters
        ----------
        wavelength : float
            The wavelength

        Returns
        -------
        None
            No return
        """
        self.wavelength = wavelength
        for layer in self.layers:
            layer.wavelength = wavelength
        self.clear_caches()

    def add_layer(self, name: str, radius: float, index: float) -> None:
        """
        Add a layer to the fiber.

        Parameters
        ----------
        name : str
            Name of the layer (e.g., "core", "cladding").
        radius : float
            Outer radius of the layer in meters.
        index : float
            Refractive index of the layer.

        Notes
        -----
        Layers should be added in order, starting with the innermost layer (core).
        """
        if not name:
            raise ValidationError("layer name must not be empty")
        if not numpy.isfinite(index) or index <= 0:
            raise ValidationError("refractive index must be positive and finite")
        radius = float(radius)
        if radius <= 0:
            radius = numpy.inf
        if radius <= self.radius_in:
            raise ValidationError("layer radii must be strictly increasing")

        self.layer_names.append(name)
        self.index_list.append(index)

        self.layer_radius.append(radius)

        layer = StepIndex(
            radius_in=self.radius_in,
            radius_out=radius,
            index_list=[index],
        )
        layer.is_last_layer = False
        layer.is_first_layer = False
        layer.wavelength = self.wavelength

        self.layers.append(layer)

        self.radius_in = radius
        self.clear_caches()

    def initialize_layers(self) -> None:
        """Initializes the layers.

        Returns
        -------
        None
            No returns
        """
        if not self.layers:
            raise ValidationError("a fiber must contain at least one layer")
        self.layers[-1].is_last_layer = True
        self.layers[0].is_first_layer = True

        self.layers[-1].radius_out = numpy.inf

        for position, layer in enumerate(self.layers):
            layer.position = position
        self.clear_caches()

    def get_layer_at_radius(self, radius: float) -> StepIndex:
        """Gets the layer that is associated to a given radius.

        Parameters
        ----------
        radius : float
            The radius

        Returns
        -------
        StepIndex
            Layer whose radial interval contains ``radius``.
        """
        radius = abs(radius)
        for layer in self.layers:
            if layer.radius_in <= radius <= layer.radius_out:
                return layer
        raise ValidationError(f"radius {radius} lies outside the fiber geometry")

    def get_index_profile(self, radius) -> numpy.ndarray:
        """Evaluate the refractive-index profile for many radii at once.

        Parameters
        ----------
        radius : array-like
            Scalar or array of radial positions in meters.

        Returns
        -------
        numpy.ndarray
            Refractive indices with the same shape as ``radius``.
        """
        radii = numpy.abs(numpy.asarray(radius, dtype=float))
        boundaries = numpy.asarray([layer.radius_out for layer in self.layers])
        indices = numpy.asarray([layer.refractive_index for layer in self.layers])
        positions = numpy.searchsorted(boundaries, radii, side="left")
        positions = numpy.clip(positions, 0, len(indices) - 1)
        return indices[positions]

    @property
    def radius(self) -> float:
        """Gets the fiber total radius taking account for all layers.

        Returns
        -------
        float
            The fiber radius.
        """
        layer_radius = [
            layer.radius_out for layer in self.layers[:-1]
        ]

        largest_radius = numpy.max(layer_radius)

        return largest_radius

    def get_index_at_radius(self, radius: float) -> float:
        """Gets the refractive index at a given radius.

        Parameters
        ----------
        radius : float
            The radius

        Returns
        -------
        float
            The refractive index at given radius.
        """
        layer = self.get_layer_at_radius(radius)

        return layer.refractive_index

    @property
    def maximum_index(self) -> float:
        """Gets the maximum refractive index of the fiber.

        Returns
        -------
        float
            The minimum index.
        """
        layers_maximum_index = [
            layer.refractive_index for layer in self.layers
        ]

        return numpy.max(layers_maximum_index)

    @property
    def minimum_index(self) -> float:
        """Gets the minimum refractive index of the fiber.

        Returns
        -------
        float
            The minimum index.
        """
        layers_maximum_index = [
            layer.refractive_index for layer in self.layers
        ]

        return numpy.min(layers_maximum_index)

    def get_NA(self) -> float:
        r"""
        Compute the numerical aperture (NA) of the fiber.

        The numerical aperture is defined as:

        .. math::
            NA = \sqrt{n_{max}^2 - n_{min}^2}

        where:
            - :math:`n_{max}` is the maximum refractive index in the fiber.
            - :math:`n_{min}` is the refractive index of the outermost layer.

        Returns
        -------
        float
            The numerical aperture of the fiber.
        """
        n_max = self.maximum_index

        last_layer = self.layers[-1]

        n_min = last_layer.refractive_index

        return numpy.sqrt(n_max**2 - n_min**2)

    @property
    def M_number(self) -> float:
        r"""Gets the m number representing an approximation of the number of existing mode
        in the fiber. It's valide only for highly multimode fibers M number is defined as:

        .. math::
            M = \frac{V^2}{2}

        Returns
        -------
        float
            The M number.
        """
        return self.V_number**2 / 2

    @property
    def V_number(self) -> float:
        r"""
        Compute the V-number (normalized frequency) of the fiber.

        The V-number is given by:

        .. math::
            V = \frac{2 \pi a}{\lambda} \cdot NA

        where:
            - :math:`a` is the core radius.
            - :math:`\lambda` is the wavelength.
            - :math:`NA` is the numerical aperture.

        Returns
        -------
        float
            The V-number of the fiber.

        Notes
        -----
        The V-number determines whether the fiber supports single-mode
        or multimode operation.
        """
        NA = self.get_NA()

        inner_radius = self.last_layer.radius_in

        V0 = (2 * numpy.pi / self.wavelength) * inner_radius * NA

        return V0

    def get_mode_cutoff_v0(self, mode: Mode) -> float:
        """Gets the cutoff wavelength of the fiber.

        Parameters
        ----------
        mode : Mode
            The mode to consider

        Returns
        -------
        float
            The cutoff wavelength.
        """
        result = self.analysis.cutoff_result(mode)
        return result.value if result.converged else numpy.nan

    def get_mode_cutoff_wavelength(self, mode: Mode) -> float:
        """
        Compute the cutoff wavelength for a given mode.

        Parameters
        ----------
        mode : Mode
            The optical mode for which to calculate the cutoff wavelength.

        Returns
        -------
        float
            The cutoff wavelength of the specified mode.

        Notes
        -----
        A mode is no longer guided when the wavelength exceeds the cutoff.
        """
        cutoff_V0 = self.get_mode_cutoff_v0(mode=mode)

        if cutoff_V0 == 0:
            return numpy.inf

        if numpy.isinf(cutoff_V0):
            return 0

        NA = self.get_NA()

        inner_radius = self.last_layer.radius_in

        cutoff_wavelength = 2 * numpy.pi / cutoff_V0 * inner_radius * NA

        return cutoff_wavelength

    def get_effective_index(self, mode: Mode) -> float:
        """Gets the effective index.

        Parameters
        ----------
        mode : Mode
            The mode to consider

        Returns
        -------
        float
            The effective index.
        """
        return self.analysis.effective_index(mode)

    def solve_effective_index(self, mode: Mode):
        """Return a structured effective-index solver result.

        Parameters
        ----------
        mode : Mode
            Mode to solve.

        Returns
        -------
        SolverResult[float]
            Value and explicit convergence diagnostics.
        """
        return self.analysis.effective_index_result(mode)

    def get_normalized_beta(self, mode: Mode) -> float:
        """Gets the normalized propagation constant [beta].

        Parameters
        ----------
        mode : Mode
            The mode to consider

        Returns
        -------
        float
            The normalized propagation constant.
        """
        neff = self.get_effective_index(mode)

        n_max = self.maximum_index

        n_last_layer = self.last_layer.refractive_index

        numerator = neff**2 - n_last_layer**2

        denominator = n_max**2 - n_last_layer**2

        return numerator / denominator

    def get_propagation_constant(self, mode: Mode) -> float:
        r"""Gets the propagation constant [:math:`beta`].

        Parameters
        ----------
        mode : Mode
            The mode to consider

        Returns
        -------
        float
            The propagation constant [:math:`beta`].
        """
        neff = self.get_effective_index(mode)

        beta = neff * (2 * numpy.pi / self.wavelength)

        return beta

    def get_phase_velocity(self, mode: Mode) -> float:
        r"""
        Compute the phase velocity of a mode in the fiber.

        Parameters
        ----------
        mode : Mode
            The optical mode for which to compute the phase velocity.

        Returns
        -------
        float
            The phase velocity of the mode, in meters per second.

        Notes
        -----
        Phase velocity is given by:

        .. math::
            v_p = \frac{c}{n_{eff}}

        where :math:`n_{eff}` is the effective refractive index.
        """
        n_eff = self.get_effective_index(mode)

        return constants.c / n_eff

    def get_group_index(self, mode: Mode) -> float:
        """Gets the group index.

        Parameters
        ----------
        mode : Mode
            The mode to consider

        Returns
        -------
        float
            The group index.
        """
        omega = c * 2 * numpy.pi / self.wavelength

        derivative = get_function_derivative(
            function=get_propagation_constant_from_omega,
            x_eval=omega,
            derivative=1,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
        )

        value = derivative * constants.c
        if not numpy.isfinite(value):
            raise ConvergenceError(f"group-index derivative failed for {mode}")
        return value

    def get_group_velocity(self, mode: Mode) -> float:
        """Compute the correctly named modal group velocity.

        Parameters
        ----------
        mode : Mode
            Mode to evaluate.

        Returns
        -------
        float
            Group velocity in meters per second.
        """
        return constants.c / self.get_group_index(mode)

    def get_group_velocity_dispersion(self, mode: Mode) -> float:
        r"""Gets the fiber group velocity dispersion defined as:

        .. math::
            \frac{\partial^2 \beta}{\partial \omega^2}

        Parameters
        ----------
        mode : Mode
            The mode to consider

        Returns
        -------
        float
            The group_velocity dispersion
        """
        omega = c * 2 * numpy.pi / self.wavelength

        derivative = get_function_derivative(
            function=get_propagation_constant_from_omega,
            x_eval=omega,
            derivative=2,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
        )

        if not numpy.isfinite(derivative):
            raise ConvergenceError(
                f"group-velocity-dispersion derivative failed for {mode}"
            )
        return derivative

    def get_dispersion(self, mode: Mode) -> float:
        r"""
        Compute the modal dispersion in the fiber.

        Dispersion is defined as:

        .. math::
            D = - \frac{2 \pi c}{\lambda^2} \cdot \frac{\partial^2 \beta}{\partial \omega^2}

        where:
            - :math:`c` is the speed of light.
            - :math:`\lambda` is the wavelength.
            - :math:`\beta` is the propagation constant.

        Parameters
        ----------
        mode : Mode
            The optical mode to consider.

        Returns
        -------
        float
            The dispersion in units of ps/nm/km.

        Notes
        -----
        Dispersion affects pulse broadening in fiber optics, particularly for long-haul transmission.
        """
        gvd = self.get_group_velocity_dispersion(mode=mode)

        factor = - 2 * numpy.pi * constants.c / self.wavelength**2

        return 1e6 * factor * gvd

    def get_S_parameter(self, mode: Mode) -> float:
        r"""
        Compute the S-parameter (third-order dispersion).

        The S-parameter is defined as:

        .. math::
            S = 10^{-3} \cdot \left( \frac{2 \pi c}{\lambda^2} \right)^2 \cdot \frac{\partial^3 \beta}{\partial \omega^3}

        Parameters
        ----------
        mode : Mode
            The optical mode to consider.

        Returns
        -------
        float
            The S-parameter in appropriate units.
        """
        omega = c * 2 * numpy.pi / self.wavelength

        derivative = get_function_derivative(
            function=get_propagation_constant_from_omega,
            x_eval=omega,
            derivative=3,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
        )

        if not numpy.isfinite(derivative):
            raise ConvergenceError(f"dispersion-slope derivative failed for {mode}")

        factor = 2 * numpy.pi * constants.c / self.wavelength**2

        return 1e-3 * derivative * factor**2

    def get_mode_field(
            self,
            mode: Mode,
            limit: float = None,
            n_point: int = 101) -> Field:
        """Get field class

        Parameters
        ----------
        mode : Mode
            The mode to consider
        limit : float
            The limit boundary
        n_point : int
            The number of point for axis discreditization

        Returns
        -------
        Field
            The field instance of the mode.
        """
        return self.fields.mode_field(mode=mode, limit=limit, n_point=n_point)

    @cache
    def get_radial_field(
        self, mode: Mode, radius: float
    ) -> tuple[CylindricalCoordinates, CylindricalCoordinates]:
        r"""Gets the mode field without the azimuthal component.
        Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

        Parameters
        ----------
        mode : Mode
            The mode to consider
        radius : float
            The radius

        Returns
        -------
        tuple[CylindricalCoordinates, CylindricalCoordinates]
            Electric and magnetic radial fields.
        """
        radial_field = get_radial_field(
            fiber=self,
            mode=mode,
            wavelength=self.wavelength,
            radius=radius
        )

        return radial_field

    def get_radial_fields(self, mode: Mode, radius) -> tuple[CylindricalCoordinates, CylindricalCoordinates]:
        """Evaluate radial electric and magnetic fields over an array.

        Parameters
        ----------
        mode : Mode
            Mode to evaluate.
        radius : array-like
            Radial sampling positions.

        Returns
        -------
        tuple[CylindricalCoordinates, CylindricalCoordinates]
            Electric and magnetic cylindrical components as arrays.
        """
        radii = numpy.asarray(radius, dtype=float)
        flat = radii.ravel()
        solved = [self.get_radial_field(mode, float(value)) for value in flat]

        def assemble(position):
            """Assemble one electric or magnetic coordinate container."""
            values = [item[position] for item in solved]
            return CylindricalCoordinates(*(
                numpy.asarray([getattr(value, name) for value in values]).reshape(radii.shape)
                for name in ("rho", "phi", "z")
            ))

        return assemble(0), assemble(1)

    def get_radial_field_norm(self, mode: Mode, radius: float) -> float:
        r"""Gets the norm of the mode field without the azimuthal component.
        Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

        Parameters
        ----------
        mode : Mode
            The mode to consider
        radius : float
            The radius

        Returns
        -------
        float
            The radial field.
        """
        e_field, h_field = self.get_radial_field(mode=mode, radius=radius)

        norm = numpy.sqrt(e_field.rho**2 + e_field.phi**2 + e_field.z**2)

        return norm

    def does_mode_exist(self, *mode_list) -> list:
        """Test whether effective-index solutions exist for several modes.

        Parameters
        ----------
        *mode_list : Mode
            Modes to evaluate.

        Returns
        -------
        list of bool
            Existence flag for each supplied mode.
        """
        mode_exist = []
        for mode in mode_list:
            neff = self.get_effective_index(mode=mode)
            if not numpy.isfinite(neff):
                mode_exist.append(False)
            else:
                mode_exist.append(True)

        return mode_exist

    def find_modes(self, families=("LP",), max_nu: int = 6, max_m: int = 6):
        """Automatically discover guided modes at the current wavelength.

        Parameters
        ----------
        families : iterable of str, optional
            Mode families to search.
        max_nu : int, optional
            Largest azimuthal order to inspect.
        max_m : int, optional
            Largest radial order to inspect.

        Returns
        -------
        list of Mode
            Guided modes found in the requested search range.
        """
        from PyFiberModes.analysis import find_modes
        return find_modes(self, families=families, max_nu=max_nu, max_m=max_m)

    def sweep(self, parameters, **kwargs):
        """Evaluate modal metrics over wavelength or a custom parameter.

        Parameters
        ----------
        parameters : array-like or mapping
            Parameter values passed to :func:`PyFiberModes.analysis.sweep_modes`.
        **kwargs
            Additional sweep configuration.

        Returns
        -------
        ModeSweepResult
            Computed modal quantities across the parameter grid.
        """
        from PyFiberModes.analysis import sweep_modes
        return sweep_modes(self, parameters, **kwargs)

    def print_data(self, data_type_list: list[str], mode_list: list[Mode]) -> None:
        """Prints the given data for the given modes.

        Parameters
        ----------
        data_type_list : list[str]
            The data type list
        mode_list : list[Mode]
            The mode list

        Returns
        -------
        None
            No return
        """
        for data_type in data_type_list:
            first_line = f"{data_type} @ wavelength: {self.wavelength}:\n"
            print(first_line, "-" * len(first_line))
            for mode in mode_list:
                data_type_string = f"get_{data_type.lower()}"
                try:
                    data = getattr(self, data_type_string)(mode=mode)
                except ConvergenceError:
                    data = numpy.nan
                output_string = f"{mode = } \t {data_type}: {data}"
                print(output_string)

            print('\n\n')


def get_fiber_from_delta_and_V0(delta: float, V0: float, wavelength: float) -> Fiber:
    """Construct an SMF28-like fiber from relative contrast and V-number.

    Parameters
    ----------
    delta : float
        Relative core-cladding refractive-index difference.
    V0 : float
        Target normalized frequency.
    wavelength : float
        Vacuum wavelength in meters.

    Returns
    -------
    Fiber
        Two-layer fiber adjusted to the requested contrast and V-number.
    """
    fiber = load_fiber(fiber_name='SMF28', wavelength=wavelength)

    core, clad = fiber.layers

    n_clad = clad.refractive_index
    n_core = n_clad / numpy.sqrt(1 - 2 * delta)

    r_core = V0 / numpy.sqrt(n_core**2 - n_clad**2) * wavelength / (2 * numpy.pi)

    fiber.layers[0].radius_out = r_core
    fiber.layers[1].radius_in = r_core
    fiber.layers[0].refractive_index = n_core

    return fiber


def load_fiber(
    fiber_name: str,
    wavelength: float | None = None,
    add_air_layer: bool = False,
    materials: MaterialRegistry = DEFAULT_MATERIALS,
) -> Fiber:
    """Load a named YAML fiber definition.

    Parameters
    ----------
    fiber_name : str
        Configuration filename without the ``.yaml`` suffix.
    wavelength : float, optional
        Vacuum wavelength in meters used to evaluate material indices.
    add_air_layer : bool, optional
        Add an outer air layer when true.
    materials : MaterialRegistry, optional
        Registry used to resolve material names in the YAML definition.

    Returns
    -------
    Fiber
        Constructed fiber with initialized step-index layers.
    """
    fiber_dict = loader.load_fiber_as_dict(
        fiber_name=fiber_name,
        wavelength=wavelength,
        order='out-to-in',
        materials=materials,
    )

    fiber = Fiber(wavelength=wavelength)

    for _, layer in fiber_dict['layers'].items():
        fiber.add_layer(**layer)

    if add_air_layer:
        fiber.add_layer(name='air', radius=numpy.inf, index=1.0)

    fiber.initialize_layers()

    return fiber

# -
