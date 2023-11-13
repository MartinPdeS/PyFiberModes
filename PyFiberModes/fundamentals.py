#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

from PyFiberModes import Wavelength, Mode, ModeFamily
from PyFiberModes.mode_instances import HE11, LP01


def get_wavelength_from_V0(fiber: object, V0: float) -> float:
    """
    Gets the wavelength associated to the V number V0.

    :param      fiber:  The fiber
    :type       fiber:  object
    :param      V0:     The V number
    :type       V0:     float

    :returns:   The wavelength from v 0.
    :rtype:     float
    """
    NA = fiber.get_NA()
    last_layer = fiber.last_layer

    wavelength = 2 * numpy.pi / V0 * last_layer.radius_in * NA

    return Wavelength(wavelength)


def get_propagation_constant_from_omega(
        omega: float,
        fiber: object,
        mode: Mode,
        delta_neff: float = 1e-6,
        lower_neff_boundary: float = None) -> float:
    """
    Gets the effective index of a given fiber and mode.

    :param      fiber:                The fiber to evaluate
    :type       fiber:                Fiber
    :param      wavelength:           The wavelength
    :type       wavelength:           Wavelength
    :param      mode:                 The mode
    :type       mode:                 Mode
    :param      delta_neff:           The delta neff
    :type       delta_neff:           float
    :param      lower_neff_boundary:  The lower neff boundary
    :type       lower_neff_boundary:  float

    :returns:   The effective index.
    :rtype:     float
    """
    wavelength = Wavelength(omega=omega)

    from PyFiberModes import solver

    if fiber.n_layer == 2:  # Standard Step-Index Fiber [SSIF]
        neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength)
    else:  # Multi-Layer Step-Index Fiber [MLSIF]
        neff_solver = solver.mlsif.NeffSolver(fiber=fiber, wavelength=wavelength)

    neff = neff_solver.solve(
        mode=mode,
        delta_neff=delta_neff,
        lower_neff_boundary=lower_neff_boundary
    )

    return neff * wavelength.k0


def get_effective_index(
        fiber,
        wavelength: Wavelength,
        mode: Mode,
        delta_neff: float = 1e-6,
        lower_neff_boundary: float = None) -> float:
    """
    Gets the effective index of a given fiber and mode.

    :param      fiber:                The fiber to evaluate
    :type       fiber:                Fiber
    :param      wavelength:           The wavelength
    :type       wavelength:           Wavelength
    :param      mode:                 The mode
    :type       mode:                 Mode
    :param      delta_neff:           The delta neff
    :type       delta_neff:           float
    :param      lower_neff_boundary:  The lower neff boundary
    :type       lower_neff_boundary:  float

    :returns:   The effective index.
    :rtype:     float
    """
    from PyFiberModes import solver

    if fiber.n_layer == 2:  # Standard Step-Index Fiber [SSIF]
        neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength)
    else:  # Multi-Layer Step-Index Fiber [MLSIF]
        neff_solver = solver.mlsif.NeffSolver(fiber=fiber, wavelength=wavelength)

    neff = neff_solver.solve(
        mode=mode,
        delta_neff=delta_neff,
        lower_neff_boundary=lower_neff_boundary
    )

    return neff


def get_cutoff_v0(
        fiber,
        wavelength: Wavelength,
        mode: Mode) -> float:
    """
    Gets the effective index of a given fiber and mode.

    :param      fiber:                The fiber to evaluate
    :type       fiber:                Fiber
    :param      wavelength:           The wavelength
    :type       wavelength:           Wavelength
    :param      mode:                 The mode
    :type       mode:                 Mode

    :returns:   The V0 value associated to cutoff.
    :rtype:     float
    """
    from PyFiberModes import solver

    if mode in [HE11, LP01]:
        return 0

    match fiber.n_layer:
        case 2:  # Standard Step-Index Fiber [SSIF|
            cutoff_solver = solver.ssif.CutoffSolver(fiber=fiber, wavelength=wavelength)
        case 3:  # Three-Layer Step-Index Fiber [TLSIF]
            cutoff_solver = solver.tlsif.CutoffSolver(fiber=fiber, wavelength=wavelength)
        case _:  # Multi-Layer Step-Index Fiber [MLSIF]
            cutoff_solver = solver.solver.FiberSolver(fiber=fiber, wavelength=wavelength)

    cutoff = cutoff_solver.solve(mode=mode)

    return cutoff


def get_radial_field(
        fiber,
        mode: Mode,
        wavelength: float,
        radius: float) -> float:
    r"""
    Gets the mode field without the azimuthal component.
    Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

    :param      fiber:       The fiber to evaluate
    :type       fiber:       Fiber
    :param      mode:        The mode to consider
    :type       mode:        Mode
    :param      wavelength:  The wavelength to consider
    :type       wavelength:  float
    :param      radius:      The radius
    :type       radius:      float

    :returns:   The radial field.
    :rtype:     float
    """
    from PyFiberModes import solver

    if fiber.n_layer == 2:  # Standard Step-Index Fiber [SSIF]
        neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength)
    else:  # Multi-Layer Step-Index Fiber [MLSIF]
        neff_solver = solver.mlsif.NeffSolver(fiber=fiber, wavelength=wavelength)

    neff = get_effective_index(
        fiber=fiber,
        wavelength=fiber.wavelength,
        mode=mode
    )

    kwargs = dict(
        nu=mode.nu,
        neff=neff,
        radius=radius
    )

    match mode.family:
        case ModeFamily.LP:
            return neff_solver.get_LP_field(**kwargs)
        case ModeFamily.TE:
            return neff_solver.get_TE_field(**kwargs)
        case ModeFamily.TM:
            return neff_solver.get_TM_field(**kwargs)
        case ModeFamily.EH:
            return neff_solver.get_EH_field(**kwargs)
        case ModeFamily.HE:
            return neff_solver.get_HE_field(**kwargs)

# -
