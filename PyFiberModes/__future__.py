#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

from PyFiberModes import Mode


def get_normalized_LP_coupling(fiber, mode_0: Mode, mode_1: Mode) -> float:
    r"""
    Gets the normalized coupling between two supermodes as defined in Equation 7.39 Jacques Bures.

    .. math::

        \tilde{C_{ij}} = \frac{0.5 k_0^2}{(\beta_0 - \beta_1) \sqrt{\beta_0 * \beta_1}} \sum_i r_i^2 \psi_0(r_i) * \psi_1(r_i)

    :param      mode_0:  The mode 0
    :type       mode_0:  Mode
    :param      mode_1:  The mode 1
    :type       mode_1:  Mode

    :returns:   The normalized coupling.
    :rtype:     float
    """
    assert mode_0.family == 'LP' and mode_1.family == 'LP', "The normalized coupling equation are only valid for scalar [LP] modes"

    beta_0 = fiber.get_propagation_constant(mode=mode_0)
    beta_1 = fiber.get_propagation_constant(mode=mode_1)

    norm_0 = get_LP_mode_norm(fiber=fiber, mode=mode_0)
    norm_1 = get_LP_mode_norm(fiber=fiber, mode=mode_1)

    coupling = 0
    for layer_in, layer_out in fiber.iterate_interfaces():
        delta_n = (layer_in.refractive_index**2 - layer_out.refractive_index**2)
        radius = layer_in.radius_out

        e_field_0, _ = fiber.get_radial_field(mode=mode_0, radius=radius)
        e_field_1, _ = fiber.get_radial_field(mode=mode_1, radius=radius)

        field_value_at_r_0 = e_field_0.rho / numpy.sqrt(norm_0)
        field_value_at_r_1 = e_field_1.rho / numpy.sqrt(norm_1)

        term = 2 * numpy.pi * delta_n * radius**2

        integral = - term * field_value_at_r_0 * field_value_at_r_1

    term_0 = abs(beta_0 - beta_1)
    term_1 = numpy.sqrt(beta_0 * beta_1)
    term_2 = 0.5 * (2 * numpy.pi / fiber.wavelength)**2
    term_3 = term_2 / (term_0 * term_1)

    coupling = integral * term_3

    return coupling


def get_LP_mode_norm(fiber, mode: Mode) -> float:
    radius_list = numpy.linspace(0, 2 * fiber.radius, 100)

    amplitudes = get_LP_mode_radial_field(
        fiber=fiber,
        mode=mode,
        radius_list=radius_list
    )

    norm = get_scalar_field_norm(rho_list=radius_list, field=amplitudes)

    return norm


def get_LP_mode_radial_normalized_field(fiber, mode: Mode, radius_list: numpy.ndarray = None) -> float:
    if radius_list is None:
        radius_list = numpy.linspace(0, 2 * fiber.radius, 200)

    scalar_field = get_LP_mode_radial_field(
        fiber=fiber,
        mode=mode,
        radius_list=radius_list
    )

    norm = get_scalar_field_norm(rho_list=radius_list, field=scalar_field)

    return scalar_field / numpy.sqrt(norm), radius_list


def get_LP_mode_radial_field(fiber, mode: Mode, radius_list: numpy.ndarray) -> float:
    amplitudes = numpy.empty(radius_list.size)

    for idx, rho in enumerate(radius_list):
        e_field, _ = fiber.get_radial_field(mode=mode, radius=rho)
        amplitudes[idx] = e_field.rho

    return amplitudes


def get_scalar_field_norm(rho_list: numpy.ndarray, field: numpy.ndarray) -> float:
    rho_list = numpy.asarray(rho_list)
    field = numpy.asarray(field)
    dr = rho_list[1] - rho_list[0]

    integral = numpy.trapz(field**2 * rho_list, dx=dr)

    norm = 2 * numpy.pi * integral

    return norm

# -
