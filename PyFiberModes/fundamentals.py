#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

from PyFiberModes.mode import Mode
from PyFiberModes import Wavelength
from PyFiberModes.mode_instances import HE11, LP01

from MPSTools.tools.coordinates import CylindricalCoordinates


def get_delta_from_fiber(fiber) -> float:
    """
    Gets the delta from fiber as defined in reference Eq. 3.82 of Jacques Bures.

    :param      fiber:  The fiber
    :type       fiber:  Fiber

    :returns:   The delta from fiber.
    :rtype:     float
    """
    core, clad = fiber.layers
    n_ratio = clad.refractive_index**2 / core.refractive_index**2
    return 0.5 * (1 - n_ratio)


def get_wavelength_from_V0(fiber: object, V0: float) -> float:
    """
    Gets the wavelength associated to the V number V0.

    :param      fiber:  The fiber
    :type       fiber:  object
    :param      V0:     The V number
    :type       V0:     float

    :returns:   The wavelength from V number.
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
        delta_neff: float = 1e-6) -> float:
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
        delta_neff=delta_neff
    )

    return neff * wavelength.k0


def get_U_parameter(
        fiber,
        wavelength: Wavelength,
        mode: Mode,
        delta_neff: float = 1e-6) -> float:
    """
    Gets the U parameter for a given fiber and mode.

    :param      fiber:                The fiber to evaluate
    :type       fiber:                Fiber
    :param      wavelength:           The wavelength to consider
    :type       wavelength:           Wavelength
    :param      mode:                 The mode
    :type       mode:                 Mode
    :param      delta_neff:           The delta neff
    :type       delta_neff:           float

    :returns:   The U parameter.
    :rtype:     float
    """
    from PyFiberModes import solver
    assert fiber.n_layer == 2, "Cannot compute U number for more than two layers"

    neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength)

    neff = neff_solver.solve(
        mode=mode,
        delta_neff=delta_neff,
    )

    U, _, _ = neff_solver.get_U_W_V_parameter(neff=neff)

    return U


def get_effective_index(
        fiber,
        wavelength: Wavelength,
        mode: Mode,
        delta_neff: float = 1e-6) -> float:
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
    )

    return neff


def get_mode_cutoff_v0(
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
        radius: float) -> tuple:
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

    :returns:   The radial field in a tupler of CylindricalCoordinates.
    :rtype:     tuple
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
        case 'LP':
            (er, ephi, ez), (hr, hphi, hz) = neff_solver.get_LP_field(**kwargs)
        case 'TE':
            (er, ephi, ez), (hr, hphi, hz) = neff_solver.get_TE_field(**kwargs)
        case 'TM':
            (er, ephi, ez), (hr, hphi, hz) = neff_solver.get_TM_field(**kwargs)
        case 'EH':
            (er, ephi, ez), (hr, hphi, hz) = neff_solver.get_EH_field(**kwargs)
        case 'HE':
            (er, ephi, ez), (hr, hphi, hz) = neff_solver.get_HE_field(**kwargs)

    e_field = CylindricalCoordinates(rho=er, phi=ephi, z=ez)
    h_field = CylindricalCoordinates(rho=hr, phi=hphi, z=hz)

    return e_field, h_field

# -
