#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from functools import cache
from dataclasses import dataclass, field
from scipy import constants
from copy import deepcopy

from PyFiberModes.stepindex import StepIndex
from PyFiberModes import Wavelength, Mode
from PyFinitDiff.tools.derivatives import get_derivative_1d
from PyFiberModes.field import Field

from PyFiberModes.fundamentals import (
    get_effective_index,
    get_cutoff_v0,
    get_radial_field,
    get_propagation_constant_from_omega
)

from MPSTools.fiber_catalogue import loader
from MPSTools.tools.coordinates import CylindricalCoordinates


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

    def __post_init__(self):
        self.wavelength = Wavelength(self.wavelength)
        self.layers_parameters = []
        self.radius_in = 0
        self.layers = []

    def scale(self, factor: float) -> None:
        scaled_fiber = deepcopy(self)
        for layer in scaled_fiber.layers:
            layer.radius_in *= factor
            layer.radius_out *= factor

        return scaled_fiber

    @property
    def n_layer(self) -> int:
        """
        Returns the number of layers in the fiber.

        :returns:   The number of layers
        :rtype:     int
        """
        return len(self.layers)

    @property
    def n_interface(self) -> int:
        """
        Returns the number of layer interfaces in the fiber.

        :returns:   The number of interface
        :rtype:     int
        """
        return len(self.layers) - 1

    @property
    def last_layer(self) -> StepIndex:
        """
        Returns the last layer

        :returns:   The last layer.
        :rtype:     StepIndex
        """
        return self.layers[-1]

    @property
    def penultimate_layer(self) -> StepIndex:
        """
        Returns the second to last layer

        :returns:   The second to last layer.
        :rtype:     StepIndex
        """
        return self.layers[-2]

    @property
    def first_layer(self):
        """
        Returns the first layer

        :returns:   The first layer.
        :rtype:     StepIndex
        """
        return self.layers[0]

    def __hash__(self):
        return hash(tuple(self.layers))

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
        """
        Adds a layer to the fiber the center-most layer have to be settled first.

        :param      name:    The name
        :type       name:    str
        :param      radius:  The radius
        :type       radius:  float
        :param      index:   The index
        :type       index:   float

        :returns:   No return
        :rtype:     None
        """
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

    def get_layer_at_radius(self, radius: float) -> StepIndex:
        """
        Gets the layer that is associated to a given radius.

        :param      radius:  The radius
        :type       radius:  float
        """
        radius = abs(radius)
        for layer in self.layers:
            if (radius > layer.radius_in) and (radius < layer.radius_out):
                return layer

    @property
    def radius(self) -> float:
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

    @property
    def maximum_index(self) -> float:
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

    @property
    def minimum_index(self) -> float:
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

    def get_NA(self) -> float:
        r"""
        Gets the numerical aperture NA defined as:

        .. math::
            NA = \sqrt{n_{max}^2-n_{min}^2}

        :returns:   The numerical aperture.
        :rtype:     float
        """
        n_max = self.maximum_index

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
            delta_neff: float = 1e-6) -> float:
        """
        Gets the effective index.

        :param      mode:    The mode to consider
        :type       mode:    Mode

        :returns:   The effective index.
        :rtype:     float
        """
        neff = get_effective_index(
            fiber=self,
            wavelength=self.wavelength,
            mode=mode
        )

        return neff

    def get_normalized_beta(self, mode: Mode) -> float:
        """
        Gets the normalized propagation constant [beta].

        :param      mode:    The mode to consider
        :type       mode:    Mode

        :returns:   The normalized propagation constant.
        :rtype:     float
        """
        neff = get_effective_index(
            fiber=self,
            wavelength=self.wavelength,
            mode=mode,
        )

        n_max = self.maximum_index

        n_last_layer = self.last_layer.refractive_index

        numerator = neff**2 - n_last_layer**2

        denominator = n_max**2 - n_last_layer**2

        return numerator / denominator

    def get_phase_velocity(self, mode: Mode) -> float:
        """
        Gets the phase velocity.

        :param      mode:    The mode to consider
        :type       mode:    Mode

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

        :param      mode:    The mode to consider
        :type       mode:    Mode

        :returns:   The group index.
        :rtype:     float
        """
        derivative = get_derivative_1d(
            function=get_propagation_constant_from_omega,
            x_eval=self.wavelength.omega,
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

        :param      mode:    The mode to consider
        :type       mode:    Mode

        :returns:   The groupe velocity.
        :rtype:     float
        """
        derivative = get_derivative_1d(
            function=get_propagation_constant_from_omega,
            x_eval=self.wavelength.omega,
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
        derivative = get_derivative_1d(
            function=get_propagation_constant_from_omega,
            x_eval=self.wavelength.omega,
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

        :param      mode:    The mode to consider
        :type       mode:    Mode

        :returns:   The s parameter.
        :rtype:     float
        """
        derivative = get_derivative_1d(
            function=get_propagation_constant_from_omega,
            x_eval=self.wavelength.omega,
            order=3,
            accuracy=4,
            delta=1e12,  # This value is critical for accurate computation
            function_kwargs=dict(fiber=self, mode=mode)
        )

        factor = 2 * numpy.pi * constants.c / self.wavelength**2

        return 1e-3 * derivative * factor**2

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
            limit = self.radius * 5.5

        field = Field(
            fiber=self,
            mode=mode,
            limit=limit,
            n_point=n_point
        )

        return field

    @cache
    def get_radial_field(self, mode: Mode, radius: float) -> CylindricalCoordinates:
        r"""
        Gets the mode field without the azimuthal component.
        Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

        :param      mode:        The mode to consider
        :type       mode:        Mode
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The radial field.
        :rtype:     CylindricalCoordinates
        """
        radial_field = get_radial_field(
            fiber=self,
            mode=mode,
            wavelength=self.wavelength,
            radius=radius
        )

        return radial_field

    def get_radial_field_norm(self, mode: Mode, radius: float) -> float:
        r"""
        Gets the norm of the mode field without the azimuthal component.
        Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

        :param      mode:        The mode to consider
        :type       mode:        Mode
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The radial field.
        :rtype:     float
        """
        e_field, h_field = get_radial_field(
            fiber=self,
            mode=mode,
            wavelength=self.wavelength,
            radius=radius
        )

        norm = numpy.sqrt(e_field.rho**2 + e_field.phi**2 + e_field.z**2)

        return norm

    def does_mode_exist(self, *mode_list) -> list:
        mode_exist = []
        for mode in mode_list:
            neff = self.get_effective_index(mode=mode)
            if neff is numpy.nan:
                mode_exist.append(False)
            else:
                mode_exist.append(True)

        return mode_exist

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
            first_line = f"{data_type} @ wavelength: {self.wavelength}:\n"
            print(first_line, "-" * len(first_line))
            for mode in mode_list:
                data_type_string = f"get_{data_type.lower()}"
                data = getattr(self, data_type_string)(mode=mode)
                output_string = f"{mode = } \t {data_type}: {data}"
                print(output_string)

            print('\n\n')


def get_fiber_from_delta_and_V0(delta: float, V0: float, wavelength: float) -> Fiber:
    fiber = load_fiber(fiber_name='SMF28', wavelength=wavelength)

    core, clad = fiber.layers

    n_clad = clad.refractive_index
    n_core = n_clad / numpy.sqrt(1 - 2 * delta)

    r_core = V0 / numpy.sqrt(n_core**2 - n_clad**2) * wavelength / (2 * numpy.pi)

    fiber.layers[0].radius_out = r_core
    fiber.layers[1].radius_in = r_core
    fiber.layers[0].refractive_index = n_core

    return fiber


def load_fiber(fiber_name: str, wavelength: float = None, add_air_layer: bool = False) -> Fiber:
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
        fiber.add_layer(**layer)

    if add_air_layer:
        fiber.add_layer(name='air', radius=numpy.inf, index=1.0)

    fiber.initialize_layers()

    return fiber

# -
