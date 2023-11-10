#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from itertools import count
from functools import cache
from scipy.optimize import fixed_point
from dataclasses import dataclass, field
from scipy import constants

from PyFiberModes.fiber_geometry.stepindex import StepIndex
from PyFiberModes import Wavelength, Mode, ModeFamily
from PyFiberModes.functions import get_derivative
from PyFiberModes.field import Field
from PyFiberModes.fundamentals import get_effective_index, get_cutoff_V0, get_radial_field

from MPSTools.fiber_catalogue import loader


class NameSpace():
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


@dataclass
class Fiber(object):
    wavelength: Wavelength
    """ Wavelength to consider """
    layer_names: list = field(default_factory=list)
    """ Name of each layers """
    layer_radius: list = field(default_factory=list)
    """ Radius of each layers """
    layer_types: list = field(default_factory=list)
    """ Type of each layers """
    material_types: list = field(default_factory=list)
    """ Material of each layers """
    index_list: list = field(default_factory=list)
    """ Refractive index of each layers """

    logger = logging.getLogger(__name__)

    def __post_init__(self):
        self.wavelength = Wavelength(self.wavelength)
        self.layers_parameters = []
        self.radius_in = 0
        self.layers = []

    @property
    def n_layer(self) -> int:
        return len(self.layer_names)

    def __hash__(self):
        return hash(tuple(self.layer_radius))

    def __getitem__(self, index: int) -> object:
        return self.layers[index]

    def update_wavelength(self, wavelength: Wavelength) -> None:
        """
        Update the wavelength of the fiber and all its layers

        :param      wavelength:  The wavelength
        :type       wavelength:  Wavelength

        :returns:   No return
        :rtype:     None
        """
        self.wavelength = wavelength
        for layer in self.layers:
            layer.wavelength = wavelength

    def add_layer(self, name: str, radius: float, index: float, material_type: str, layer_type: str) -> None:
        self.layer_names.append(name)
        self.layer_types.append(layer_type)
        self.material_types.append(material_type)
        self.index_list.append(index)

        if name != 'cladding':
            self.layer_radius.append(radius)

        layer = StepIndex(
            radius_in=self.radius_in,
            radius_out=radius,
            material_type=material_type,
            index_list=[index],
        )
        layer.wavelength = self.wavelength

        self.layers.append(layer)

        self.radius_in = radius

    def initialize_layers(self) -> None:
        """
        Initializes the layers.

        :returns:   No returns
        :rtype:     None
        """
        self.layers[-1].radius_out = numpy.inf

    def get_layer_name(self, layer_index: int) -> str:
        """
        Gets the layer name.

        :param      layer_index:  The layer index
        :type       layer_index:  int

        :returns:   The layer name.
        :rtype:     str
        """
        return self.layer_names[layer_index]

    def get_layer_at_radius(self, radius: float):
        """
        Gets the layer that is associated to a given radius.

        :param      radius:  The radius
        :type       radius:  float
        """
        radius = abs(radius)
        for idx, layer_radius in enumerate(self.layer_radius):
            if radius < layer_radius:
                return self.layers[idx]

        return self.layers[-1]

    def get_inner_radius(self, layer_idx: int) -> float:
        """
        Gets the radius of the inner most layer

        :param      layer_idx:  The layer index
        :type       layer_idx:  int

        :returns:   The inner radius.
        :rtype:     float
        """
        if layer_idx < 0:
            layer_idx = len(self.layer_radius) + layer_idx + 1

        if layer_idx != 0:
            return self.layer_radius[layer_idx - 1]

        return 0

    def get_outer_radius(self, layer_idx: int) -> float:
        """
        Gets the radius of the outer most layer

        :param      layer_idx:  The layer index
        :type       layer_idx:  int

        :returns:   The inner radius.
        :rtype:     float
        """
        if layer_idx < len(self.layer_radius):
            return self.layer_radius[layer_idx]
        return float("inf")

    def get_layer_thickness(self, layer_idx: int) -> float:
        """
        Gets the thickness of a specific layer.

        :param      layer_idx:  The layer index
        :type       layer_idx:  int

        :returns:   The thickness.
        :rtype:     float
        """
        outer_radius = self.get_outer_radius(layer_idx=layer_idx)
        inner_radius = self.get_inner_radius(layer_idx=layer_idx)
        return outer_radius - inner_radius

    def get_fiber_radius(self) -> float:
        """
        Gets the fiber total radius taking account for all layers.

        :returns:   The fiber radius.
        :rtype:     float
        """
        layer_radius = [
            layer.radius_out for layer in self.layers
        ]

        largest_radius = numpy.max(layer_radius)

        return largest_radius

    def get_index_at_radius(self, radius: float) -> float:
        """
        Gets the refractive index at a given radius.

        :param      radius:      The radius
        :type       radius:      float

        :returns:   The refractive index at given radius.
        :rtype:     float
        """
        layer = self.get_layer_at_radius(radius)

        return layer.index(radius)

    def get_layer_minimum_index(self, layer_idx: int) -> float:
        """
        Gets the minimum refractive index of the layers.

        :param      layer_idx:   The layer index
        :type       layer_idx:   int

        :returns:   The minimum index.
        :rtype:     float
        """
        layer = self.layers[layer_idx]

        return layer.get_minimum_index()

    def get_maximum_index(self) -> float:
        """
        Gets the maximum refractive index of the fiber.

        :param      layer_idx:   The layer index
        :type       layer_idx:   int

        :returns:   The minimum index.
        :rtype:     float
        """
        layers_maximum_index = [
            layer.get_maximum_index() for layer in self.layers
        ]

        return numpy.max(layers_maximum_index)

    def get_minimum_index(self) -> float:
        """
        Gets the minimum refractive index of the fiber.

        :param      layer_idx:   The layer index
        :type       layer_idx:   int

        :returns:   The minimum index.
        :rtype:     float
        """
        layers_maximum_index = [
            layer.get_minimum_index() for layer in self.layers
        ]

        return numpy.min(layers_maximum_index)

    def get_layer_maximum_index(self, layer_idx: int) -> float:
        """
        Gets the maximum refractive index of the layers.

        :param      layer_idx:   The layer index
        :type       layer_idx:   int

        :returns:   The minimum index.
        :rtype:     float
        """
        layer = self.layers[layer_idx]

        return layer.get_maximum_index()

    def get_NA(self) -> float:
        """
        Gets the numerical aperture NA.

        :returns:   The numerical aperture.
        :rtype:     float
        """
        n_max = self.get_maximum_index()

        last_layer = self.layers[-1]

        n_min = last_layer.get_minimum_index()

        return numpy.sqrt(n_max**2 - n_min**2)

    def get_V0(self) -> float:
        """
        Gets the V0 parameter.

        :returns:   The parameter V0.
        :rtype:     float
        """
        NA = self.get_NA()

        inner_radius = self.layers[-1].radius_in

        V0 = self.wavelength.k0 * inner_radius * NA

        return V0

    def V0_to_wavelength(self, V0: float, maxiter: int = 500, tolerance: float = 1e-15) -> float:
        """
        Convert V0 number to wavelength.
        An iterative method is used, since the index can be wavelength dependant.

        :param      V0:       The V0
        :type       V0:       float
        :param      maxiter:  The maxiter
        :type       maxiter:  int
        :param      tol:      The tolerance
        :type       tol:      float

        :returns:   The associated wavelength
        :rtype:     float
        """
        if V0 == 0:
            return numpy.inf

        if numpy.isinf(V0):
            return 0

        def model(wl):
            last_layer = self.layers[-1]
            NA = self.get_NA()
            return 2 * numpy.pi / V0 * last_layer.radius_in * NA

        wavelength = model(wl=1.55e-6)

        if abs(wavelength - model(wl=wavelength)) > tolerance:
            for w in (1.55e-6, 5e-6, 10e-6):
                try:
                    wavelength = fixed_point(model, w, xtol=tolerance, maxiter=maxiter)
                except RuntimeError:
                    # FIXME: What should we do if it does not converge?
                    self.logger.info(
                        f"V0_to_wavelength: did not converged from {w * 1e6}µm for {V0=} ({wavelength=})"
                    )
                if wavelength > 0:
                    break

        if wavelength == 0:
            self.logger.error(
                f"V0_to_wavelength: did not converged for {V0=} {wavelength=})"
            )

        return Wavelength(wavelength)

    def get_cutoff_V0(self, mode: Mode) -> float:
        """
        Gets the cutoff wavelength of the fiber.

        :param      mode:  The mode to consider
        :type       mode:  Mode

        :returns:   The cutoff wavelength.
        :rtype:     float
        """
        cutoff_V0 = get_cutoff_V0(
            mode=mode,
            fiber=self,
            wavelength=self.wavelength
        )

        return cutoff_V0

    def get_cutoff_wavelength(self, mode: Mode) -> float:
        """
        Gets the cutoff wavelength.

        :param      mode:  The mode to consider
        :type       mode:  Mode

        :returns:   The cutoff wavelength.
        :rtype:     float
        """
        cutoff_V0 = self.get_cutoff_V0(mode=mode)
        wavelength = self.V0_to_wavelength(V0=cutoff_V0)
        return wavelength

    def get_effective_index(self,
            mode: Mode,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        """
        Gets the effective index.

        :param      mode:                   The mode to consider
        :type       mode:                   Mode
        :param      delta_neff:             The discretization for research of neff value
        :type       delta_neff:             float
        :param      lower_neff_boundary:    The minimum value neff can reach
        :type       lower_neff_boundary:    float

        :returns:   The effective index.
        :rtype:     float
        """
        neff = get_effective_index(
            fiber=self,
            wavelength=self.wavelength,
            mode=mode,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        return neff

    def get_derivative_beta_vs_omega(self,
            omega: float,
            mode: Mode,
            order: float = 0,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        """
        Gets the derivative of beta vs omega.

        :param      mode:                 The mode to consider
        :type       mode:                 Mode
        :param      order:                The order of the derivative
        :type       order:                float
        :param      delta_neff:           The discretization for research of neff value
        :type       delta_neff:           float
        :param      lower_neff_boundary:  The minimum value neff can reach
        :type       lower_neff_boundary:  float

        :returns:   The derivative beta vs omega.
        :rtype:     float
        """

        wl = Wavelength(omega=omega)
        if order == 0:
            neff = get_effective_index(fiber=self, wavelength=wl, mode=mode, delta_neff=delta_neff, lower_neff_boundary=lower_neff_boundary)
            return neff * wl.k0

        function_kwargs = dict(
            mode=mode,
            order=0,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        return get_derivative(
            function=self.get_derivative_beta_vs_omega,
            x=omega,
            order=order,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=function_kwargs
        )

    def get_normalized_beta(self,
            mode: Mode,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        """
        Gets the normalized propagation constant [beta].

        :param      mode:                   The mode to consider
        :type       mode:                   Mode
        :param      delta_neff:             The discretization for research of neff value
        :type       delta_neff:             float
        :param      lower_neff_boundary:    The minimum value neff can reach
        :type       lower_neff_boundary:    float

        :returns:   The normalized propagation constant.
        :rtype:     float
        """
        neff = get_effective_index(
            fiber=self,
            wavelength=self.wavelength,
            mode=mode,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        last_layer = self.layers[-1]

        n_max = self.get_maximum_index()

        n_min = last_layer.get_minimum_index()

        numerator = neff**2 - n_min**2

        denominator = n_max**2 - n_min**2

        return numerator / denominator

    def get_phase_velocity(self,
            mode: Mode,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        """
        Gets the phase velocity.

        :param      mode:                   The mode to consider
        :type       mode:                   Mode
        :param      delta_neff:             The discretization for research of neff value
        :type       delta_neff:             float
        :param      lower_neff_boundary:    The minimum value neff can reach
        :type       lower_neff_boundary:    float

        :returns:   The phase velocity.
        :rtype:     float
        """
        n_eff = get_effective_index(
            fiber=self,
            wavelength=self.wavelength,
            mode=mode,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        return constants.c / n_eff

    def get_group_index(self,
            mode: Mode,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        """
        Gets the group index.

        :param      mode:                   The mode to consider
        :type       mode:                   Mode
        :param      delta_neff:             The discretization for research of neff value
        :type       delta_neff:             float
        :param      lower_neff_boundary:    The minimum value neff can reach
        :type       lower_neff_boundary:    float

        :returns:   The group index.
        :rtype:     float
        """
        derivative = self.get_derivative_beta_vs_omega(
            omega=self.wavelength.omega,
            mode=mode,
            order=1,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        return derivative * constants.c

    def get_groupe_velocity(self,
            mode: Mode,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        """
        Gets the groupe velocity.

        :param      mode:                  The mode to consider
        :type       mode:                  Mode
        :param      delta_neff:            The discretization for research of neff value
        :type       delta_neff:            float
        :param      lower_neff_boundary:   The minimum value neff can reach
        :type       lower_neff_boundary:   float

        :returns:   The groupe velocity.
        :rtype:     float
        """
        derivative = self.get_derivative_beta_vs_omega(
            omega=self.wavelength.omega,
            mode=mode,
            order=1,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        return 1 / derivative

    def get_dispersion(self,
            mode: Mode,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        r"""
        Gets the modal dispersion :math:`\frac{2 * \pi * c}{\lambda^2} * \frac{\partial^2 beta}{\partial omega}`.

        :param      mode:                   The mode to consider
        :type       mode:                   Mode
        :param      delta_neff:             The discretization for research of neff value
        :type       delta_neff:             float
        :param      lower_neff_boundary:    The minimum value neff can reach
        :type       lower_neff_boundary:    float

        :returns:   The modal dispersion.
        :rtype:     float
        """
        derivative = self.get_derivative_beta_vs_omega(
            omega=self.wavelength.omega,
            mode=mode,
            order=2,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        factor = 2 * numpy.pi * constants.c / self.wavelength**2

        return 1e6 * -derivative * factor

    def get_S_parameter(self,
            mode: Mode,
            delta_neff: float = 1e-6,
            lower_neff_boundary: float = None) -> float:
        """
        Gets the s parameter.

        :param      mode:                   The mode to consider
        :type       mode:                   Mode
        :param      delta_neff:             The discretization for research of neff value
        :type       delta_neff:             float
        :param      lower_neff_boundary:    The minimum value neff can reach
        :type       lower_neff_boundary:    float

        :returns:   The s parameter.
        :rtype:     float
        """
        derivative = self.get_derivative_beta_vs_omega(
            omega=self.wavelength.omega,
            mode=mode,
            order=3,
            delta_neff=delta_neff,
            lower_neff_boundary=lower_neff_boundary
        )

        factor = 2 * numpy.pi * constants.c / self.wavelength**2

        return 1e-3 * derivative * factor**2

    def get_vectorial_modes(self,
            wavelength: float,
            nu_max=None,
            m_max=None,
            delta: float = 1e-6) -> set:
        """
        Gets the family of vectorial modes.

        :param      wavelength:  The wavelength to consider
        :type       wavelength:  float
        :param      nu_max:      The maximum value of nu parameter
        :type       nu_max:      int
        :param      m_max:       The maximum value of m parameter
        :type       m_max:       int
        :param      delta_neff:  The discretization for research of neff value
        :type       delta_neff:  float

        :returns:   The vectorial modes.
        :rtype:     set
        """
        families = (
            ModeFamily.HE,
            ModeFamily.EH,
            ModeFamily.TE,
            ModeFamily.TM
        )

        modes = self.get_modes_from_familly(
            families=families,
            wavelength=wavelength,
            nu_max=nu_max,
            m_max=m_max,
            delta=delta
        )

        return modes

    def get_LP_modes(self,
            wavelength: float,
            ellmax: int = None,
            m_max: int = None,
            delta_neff: float = 1e-6) -> set:
        """
        Gets the family of LP modes.

        :param      wavelength:  The wavelength to consider
        :type       wavelength:  float
        :param      ellmax:      The ellmax
        :type       ellmax:      int
        :param      mmax:        The maximum value of m parameter
        :type       mmax:        int
        :param      delta_neff:  The discretization for research of neff value
        :type       delta_neff:  float

        :returns:   The lp modes.
        :rtype:     set
        """
        families = (ModeFamily.LP,)

        modes = self.get_modes_from_familly(
            families=families,
            wavelength=wavelength,
            nu_max=ellmax,
            m_max=m_max,
            delta_neff=delta_neff
        )

        return modes

    def get_modes_from_familly(self,
            families,
            wavelength: float,
            nu_max: int = numpy.inf,
            m_max: int = numpy.inf,
            delta_neff: float = 1e-6) -> set:
        """
        Find all modes of given families, within given constraints

        :param      families:         The families
        :type       families:         object
        :param      wavelength:       The wavelength to consider
        :type       wavelength:       float
        :param      nu_max:           The radial number nu maximum to reach
        :type       nu_max:           int
        :param      m_max:            The azimuthal number m maximum to reach
        :type       m_max:            int
        :param      delta_neff:       The discretization for research of neff value
        :type       delta_neff:       float

        :returns:   The mode to considers from familly.
        :rtype:     set
        """
        modes = set()
        v0 = self.get_V0(wavelength=wavelength)

        for family in families:
            for nu in count(0):

                try:
                    _mmax = m_max[nu]
                except IndexError:
                    _mmax = m_max[-1]
                except TypeError:
                    _mmax = m_max

                if family in [ModeFamily.TE, ModeFamily.TM] and nu > 0:
                    break

                if family in [ModeFamily.HE, ModeFamily.EH] and nu == 0:
                    continue

                if nu > nu_max:
                    break

                for m in count(1):
                    if m > _mmax:
                        break

                    mode = Mode(family, nu, m)

                    try:
                        if self.get_cutoff_V0(mode=mode) > v0:
                            break

                    except (NotImplementedError, ValueError):
                        neff = get_effective_index(fiber=self, wavelength=self.wavelength, mode=mode, delta_neff=delta_neff)

                        if numpy.isnan(neff):
                            break

                    modes.add(mode)

                if m == 1:
                    break
        return modes

    def get_field(self,
            mode: Mode,
            limit: float = None,
            n_point: int = 101) -> Field:
        """
        Get field class

        :param      mode:        The mode to consider
        :type       mode:        Mode
        :param      wavelength:  The wavelength to consider
        :type       wavelength:  float
        :param      limit:       The limit boundary
        :type       limit:       float
        :param      n_point:     The number of point for axis discreditization
        :type       n_point:     int

        :returns:   The field instance of the mode.
        :rtype:     Field
        """
        if limit is None:
            limit = self.get_fiber_radius()

        field = Field(
            fiber=self,
            mode=mode,
            limit=limit,
            n_point=n_point
        )

        return field

    @cache
    def get_radial_field(self,
            mode: Mode,
            wavelength: float,
            radius: float) -> float:
        """
        Gets the mode field without the azimuthal component.

        :param      mode:        The mode to consider
        :type       mode:        Mode
        :param      wavelength:  The wavelength to consider
        :type       wavelength:  float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The radial field.
        :rtype:     float
        """
        radial_field = get_radial_field(
            fiber=self,
            mode=mode,
            wavelength=self.wavelength,
            radius=radius
        )

        return radial_field


def load_fiber(fiber_name: str, wavelength: float = None) -> Fiber:
    """
    Loads a fiber as type that suit PyFiberModes.

    :param      fiber_name:  The fiber name
    :type       fiber_name:  str
    :param      wavelength:  The wavelength to consider
    :type       wavelength:  float

    :returns:   The loaded fiber
    :rtype:     Fiber
    """
    fiber_dict = loader.load_fiber_as_dict(
        fiber_name=fiber_name,
        wavelength=wavelength,
        order='out-to-in'
    )

    fiber = Fiber(wavelength=wavelength)

    for _, layer in fiber_dict['layers'].items():
        if layer.get('name') in ['air']:
            continue

        fiber.add_layer(material_type='Fixed', layer_type='StepIndex', **layer)

    fiber.initialize_layers()

    return fiber

# -
