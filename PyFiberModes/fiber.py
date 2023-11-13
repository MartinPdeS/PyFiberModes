#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from itertools import count
from functools import cache
from dataclasses import dataclass, field
from scipy import constants

from PyFiberModes.stepindex import StepIndex
from PyFiberModes import Wavelength, Mode, ModeFamily
from PyFiberModes.functions import get_derivative
from PyFiberModes.field import Field

from PyFiberModes.fundamentals import (
    get_effective_index,
    get_cutoff_v0,
    get_radial_field,
    get_propagation_constant_from_omega
)

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
        return len(self.layers)

    @property
    def n_interface(self) -> int:
        return len(self.layers) - 1

    @property
    def last_layer(self):
        return self.layers[-1]

    @property
    def penultimate_layer(self):
        return self.layers[-2]

    @property
    def first_layer(self):
        return self.layers[0]

    def __hash__(self):
        return hash(tuple(self.layer_radius))

    def __getitem__(self, index: int) -> object:
        return self.layers[index]

    def iterate_interfaces(self) -> tuple[object, object]:
        """
        Iterates through pair of layers that forms interfaces

        :returns:   The two layers that form the interfaces.
        :rtype:     tuple[object, object]
        """
        for layer_idx in range(1, self.n_layer):
            layer_in = self.layers[layer_idx - 1]
            layer_out = self.layers[layer_idx]
            yield layer_in, layer_out

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

    def add_layer(self, name: str, radius: float, index: float) -> None:
        self.layer_names.append(name)
        self.index_list.append(index)

        if name != 'cladding':
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

    def initialize_layers(self) -> None:
        """
        Initializes the layers.

        :returns:   No returns
        :rtype:     None
        """
        self.layers[-1].is_last_layer = True
        self.layers[0].is_first_layer = True

        self.layers[-1].radius_out = numpy.inf

        for position, layer in enumerate(self.layers):
            layer.position = position

    def get_layer_at_radius(self, radius: float):
        """
        Gets the layer that is associated to a given radius.

        :param      radius:  The radius
        :type       radius:  float
        """
        radius = abs(radius)
        for layer in self.layers:
            if (radius > layer.radius_in) and (radius < layer.radius_out):
                return layer

    def get_fiber_radius(self) -> float:
        """
        Gets the fiber total radius taking account for all layers.

        :returns:   The fiber radius.
        :rtype:     float
        """
        layer_radius = [
            layer.radius_out for layer in self.layers[:-1]
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

        return layer.refractive_index

    def get_maximum_index(self) -> float:
        """
        Gets the maximum refractive index of the fiber.

        :param      layer_idx:   The layer index
        :type       layer_idx:   int

        :returns:   The minimum index.
        :rtype:     float
        """
        layers_maximum_index = [
            layer.refractive_index for layer in self.layers
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
            layer.refractive_index for layer in self.layers
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

        return layer.refractive_index

    def get_NA(self) -> float:
        """
        Gets the numerical aperture NA.

        :returns:   The numerical aperture.
        :rtype:     float
        """
        n_max = self.get_maximum_index()

        last_layer = self.layers[-1]

        n_min = last_layer.refractive_index

        return numpy.sqrt(n_max**2 - n_min**2)

    def get_M_number(self) -> float:
        r"""
        Gets the m number representing an approximation of the number of existing mode
        in the fiber. It's valide only for highly multimode fibers
        M number is defined as:

        .. math::
            M = \frac{V^2}{2}

        :returns:   The M number.
        :rtype:     float
        """
        pass

    def get_V0(self) -> float:
        r"""
        Gets the V number parameter defined as:

        .. math::
            V = \frac{2 * pi * a}{\lambda} * NA

        reference: https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=10417#:~:text=Single%20mode%20fibers%20are%20defined,fiber%20at%20the%20same%20wavelength.

        :returns:   The parameter V0.
        :rtype:     float
        """
        NA = self.get_NA()

        inner_radius = self.last_layer.radius_in

        V0 = self.wavelength.k0 * inner_radius * NA

        return V0

    def get_cutoff_v0(self, mode: Mode) -> float:
        """
        Gets the cutoff wavelength of the fiber.

        :param      mode:  The mode to consider
        :type       mode:  Mode

        :returns:   The cutoff wavelength.
        :rtype:     float
        """
        cutoff_V0 = get_cutoff_v0(
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
        cutoff_V0 = self.get_cutoff_v0(mode=mode)

        if cutoff_V0 == 0:
            return numpy.inf

        if numpy.isinf(cutoff_V0):
            return 0

        NA = self.get_NA()

        inner_radius = self.last_layer.radius_in

        cutoff_wavelength = 2 * numpy.pi / cutoff_V0 * inner_radius * NA

        return Wavelength(cutoff_wavelength)

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

    def get_normalized_beta(self, mode: Mode) -> float:
        """
        Gets the normalized propagation constant [beta].

        :param      mode:                   The mode to consider
        :type       mode:                   Mode

        :returns:   The normalized propagation constant.
        :rtype:     float
        """
        neff = get_effective_index(
            fiber=self,
            wavelength=self.wavelength,
            mode=mode,
        )

        n_max = self.get_maximum_index()

        n_last_layer = self.last_layer.refractive_index

        numerator = neff**2 - n_last_layer**2

        denominator = n_max**2 - n_last_layer**2

        return numerator / denominator

    def get_phase_velocity(self, mode: Mode) -> float:
        """
        Gets the phase velocity.

        :param      mode:                   The mode to consider
        :type       mode:                   Mode

        :returns:   The phase velocity.
        :rtype:     float
        """
        n_eff = get_effective_index(
            fiber=self,
            wavelength=self.wavelength,
            mode=mode,
        )

        return constants.c / n_eff

    def get_group_index(self, mode: Mode) -> float:
        """
        Gets the group index.

        :param      mode:                   The mode to consider
        :type       mode:                   Mode

        :returns:   The group index.
        :rtype:     float
        """
        derivative = get_derivative(
            function=get_propagation_constant_from_omega,
            x=self.wavelength.omega,
            order=1,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
        )

        return derivative * constants.c

    def get_groupe_velocity(self, mode: Mode) -> float:
        r"""
        Gets the groupe velocity defined as:

        .. math::
            \left( \frac{\partial \beta}{\partial \omega} \right)^{-1}

        :param      mode:                  The mode to consider
        :type       mode:                  Mode

        :returns:   The groupe velocity.
        :rtype:     float
        """
        derivative = get_derivative(
            function=get_propagation_constant_from_omega,
            x=self.wavelength.omega,
            order=1,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
        )

        return 1 / derivative

    def get_group_velocity_dispersion(self, mode: Mode) -> float:
        r"""
        Gets the fiber group velocity dispersion defined as:

        .. math::
            \frac{\partial^2 \beta}{\partial \omega^2}

        :param      mode:   The mode to consider
        :type       mode:   Mode

        :returns:   The group_velocity dispersion
        :rtype:     float
        """
        derivative = get_derivative(
            function=get_propagation_constant_from_omega,
            x=self.wavelength.omega,
            order=2,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
        )

        return derivative

    def get_dispersion(self, mode: Mode) -> float:
        r"""
        Gets the fiber dispersion defined as:

        .. math::
            10^6 * \frac{2 \pi c}{\lambda^2} * \frac{\partial^2 \beta}{\partial \omega}

        :param      mode:   The mode to consider
        :type       mode:   Mode

        :returns:   The modal dispersion in units of ps/nm/km.
        :rtype:     float
        """
        gvd = self.get_group_velocity_dispersion(mode=mode)

        factor = - 2 * numpy.pi * constants.c / self.wavelength**2

        return 1e6 * factor * gvd

    def get_S_parameter(self, mode: Mode) -> float:
        r"""
        Gets the s parameter defined as:

        .. math::
            10^{-3} * \left( \frac{2 \pi c}{\lambda^2} \right)^2 * \frac{\partial^3 \beta}{\partial \omega}

        :param      mode:                   The mode to consider
        :type       mode:                   Mode

        :returns:   The s parameter.
        :rtype:     float
        """
        derivative = get_derivative(
            function=get_propagation_constant_from_omega,
            x=self.wavelength.omega,
            order=3,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
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
        v0 = self.get_V0()

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
                        if self.get_cutoff_v0(mode=mode) > v0:
                            break

                    except (NotImplementedError, ValueError):
                        neff = get_effective_index(fiber=self, wavelength=self.wavelength, mode=mode, delta_neff=delta_neff)

                        if numpy.isnan(neff):
                            break

                    modes.add(mode)

                if m == 1:
                    break
        return modes

    def get_mode_field(self,
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
            limit = self.get_fiber_radius() * 1.2

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
        r"""
        Gets the mode field without the azimuthal component.
        Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

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

    def print_data(self, data_type_list: list[str], mode_list: list[Mode]) -> None:
        """
        Prints the given data for the given modes.

        :param      data_type_list:  The data type list
        :type       data_type_list:  list[str]
        :param      mode_list:       The mode list
        :type       mode_list:       list[Mode]

        :returns:   No return
        :rtype:     None
        """
        for data_type in data_type_list:
            print(f"{data_type} @ wavelength: {self.wavelength}", '\n')
            for mode in mode_list:
                data_type_string = f"get_{data_type.lower()}"
                data = getattr(self, data_type_string)(mode=mode)
                output_string = f"{mode = } \t {data_type}: {data}"
                print(output_string)

            print('\n\n')


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

        fiber.add_layer(**layer)

    fiber.initialize_layers()

    return fiber

# -
