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

    if mode_0.nu > 0 or mode_1.nu > 0:
        return 0

    norm_0 = get_LP_mode_norm(fiber=fiber, mode=mode_0)
    norm_1 = get_LP_mode_norm(fiber=fiber, mode=mode_1)

    beta_0 = fiber.get_propagation_constant(mode=mode_0)
    beta_1 = fiber.get_propagation_constant(mode=mode_1)
    # testing--------------------------------------------

    # radius_list = numpy.linspace(0, 1.2 * fiber.radius, 100)

    # normalized_field_0 = get_LP_mode_radial_normalized_field(
    #     fiber=fiber,
    #     mode=mode_0,
    #     radius_list=radius_list
    # )

    # normalized_field_1 = get_LP_mode_radial_normalized_field(
    #     fiber=fiber,
    #     mode=mode_1,
    #     radius_list=radius_list
    # )

    # idx_interface = numpy.argmin((radius_list - fiber.layers[0].radius_out)**2)

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(radius_list, normalized_field_0)
    # plt.plot(radius_list, normalized_field_1)
    # plt.vlines(x=radius_list[idx_interface], ymin=-40000, ymax=40000, color='red')
    # plt.vlines(x=fiber.first_layer.radius_out, ymin=-40000, ymax=40000, color='blue')
    # plt.show()

    # # testing--------------------------------------------

    coupling = 0
    for layer_in, layer_out in fiber.iterate_interfaces():
        radius = layer_in.radius_out

        index_in, index_out = layer_in.refractive_index, layer_out.refractive_index

        delta_index = (index_in**2 - index_out**2)

        e_field_0, _ = fiber.get_radial_field(mode=mode_0, radius=radius)
        e_field_1, _ = fiber.get_radial_field(mode=mode_1, radius=radius)

        field_0 = e_field_0.rho / numpy.sqrt(norm_0)
        field_1 = e_field_1.rho / numpy.sqrt(norm_1)

        interface_coupling = radius**2 * field_0 * field_1 * delta_index

        coupling += interface_coupling

    term_0 = abs(beta_0 - beta_1)
    term_1 = numpy.sqrt(beta_0 * beta_1)
    term_2 = 0.5 * fiber.wavelength.k0**2
    term_3 = term_2 / (term_0 * term_1)

    coupling = coupling * term_3

    return coupling


def get_LP_mode_norm(fiber, mode: Mode) -> float:
    radius_list = numpy.linspace(0, 2 * fiber.radius, 100)

    amplitudes = get_LP_mode_radial_field(
        fiber=fiber,
        mode=mode,
        radius_list=radius_list
    )

    norm = get_scalar_field_norm(x_list=radius_list, field=amplitudes)

    return norm


def get_LP_mode_radial_normalized_field(fiber, mode: Mode, radius_list: numpy.ndarray) -> float:
    scalar_field = get_LP_mode_radial_field(
        fiber=fiber,
        mode=mode,
        radius_list=radius_list
    )

    norm = get_scalar_field_norm(x_list=radius_list, field=scalar_field)

    return scalar_field / numpy.sqrt(norm)


def get_LP_mode_radial_field(fiber, mode: Mode, radius_list: numpy.ndarray) -> float:
    amplitudes = numpy.empty(radius_list.size)

    for idx, rho in enumerate(radius_list):
        e_field, _ = fiber.get_radial_field(mode=mode, radius=rho)
        amplitudes[idx] = e_field.rho

    return amplitudes


def get_scalar_field_norm(x_list: numpy.ndarray, field: numpy.ndarray) -> float:
    x_list = numpy.asarray(x_list)
    field = numpy.asarray(field)
    dx = x_list[1] - x_list[0]

    integral = numpy.trapz(field**2 * x_list, dx=dx)

    norm = 2 * numpy.pi * integral

    return norm

# -
