#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFiberModes import Wavelength, Mode, ModeFamily
from PyFiberModes.mode_instances import HE11, LP01
from PyFiberModes import solver


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
    n_layers = len(fiber.layers)

    if n_layers == 2:  # Standard Step-Index Fiber [SSIF]
        neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength)
    else:  # Multi-Layer Step-Index Fiber [MLSIF]
        neff_solver = solver.mlsif.NeffSolver(fiber=fiber, wavelength=wavelength)

    neff = neff_solver.solve(
        mode=mode,
        delta_neff=delta_neff,
        lower_neff_boundary=lower_neff_boundary
    )

    return neff


def get_cutoff_V0(
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
    if mode in [HE11, LP01]:
        return 0

    n_layers = len(fiber.layers)

    match n_layers:
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
    """
    Gets the mode field without the azimuthal component.

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
    n_layers = len(fiber.layers)

    if n_layers == 2:  # Standard Step-Index Fiber [SSIF]
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
