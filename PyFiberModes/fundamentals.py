#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import numpy as np
from PyFiberModes.mode import Mode
from PyFiberModes.mode_instances import HE11, LP01
from scipy.constants import c

from PyFiberModes.coordinates import CylindricalCoordinates


def get_delta_from_fiber(fiber) -> float:
    r"""
    Calculate the relative index difference, \(\Delta\), of the fiber.

    .. math::
        \Delta = \frac{1}{2} \left( 1 - \frac{n_{\text{clad}}^2}{n_{\text{core}}^2} \right)

    Parameters
    ----------
    fiber : Fiber
        The fiber object containing core and cladding properties.

    Returns
    -------
    float
        The relative index difference, \(\Delta\).
    """
    core, clad = fiber.layers
    n_ratio = clad.refractive_index**2 / core.refractive_index**2
    return 0.5 * (1 - n_ratio)


def get_wavelength_from_V0(fiber: object, V0: float) -> float:
    r"""
    Compute the wavelength corresponding to a given V-number, \(V_0\).

    .. math::
        \lambda = \frac{2 \pi a \cdot \text{NA}}{V_0}

    where:
    - \(a\) is the core radius.
    - \(\text{NA}\) is the numerical aperture.

    Parameters
    ----------
    fiber : Fiber
        The fiber object.
    V0 : float
        The V-number.

    Returns
    -------
    float
        The wavelength corresponding to V$_0$.
    """
    NA = fiber.get_NA()
    last_layer = fiber.last_layer
    wavelength = 2 * np.pi / V0 * last_layer.radius_in * NA
    return wavelength


def get_propagation_constant_from_omega(
        omega: float,
        fiber: object,
        mode: Mode,
        delta_neff: float = 1e-6) -> float:
    r"""
    Calculate the propagation constant, \(\beta\), for a given angular frequency, fiber, and mode.

    .. math::
        \beta = k_0 n_{\text{eff}}

    where:
    - \(k_0 = \frac{2 \pi}{\lambda}\) is the free-space wave number.
    - \(n_{\text{eff}}\) is the effective refractive index.

    Parameters
    ----------
    omega : float
        The angular frequency.
    fiber : Fiber
        The fiber object.
    mode : Mode
        The mode of interest.
    delta_neff : float, optional
        Convergence threshold for \(n_{\text{eff}}\) calculation, by default \(1 \times 10^{-6}\).

    Returns
    -------
    float
        The propagation constant, \(\beta\).
    """
    wavelength = c * 2 * np.pi / omega

    from PyFiberModes import solver

    neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength) if fiber.n_layer == 2 \
        else solver.mlsif.NeffSolver(fiber=fiber, wavelength=wavelength)

    neff = neff_solver.solve(mode=mode, delta_neff=delta_neff)
    return neff * (2 * numpy.pi / wavelength)


def get_U_parameter(
        fiber,
        wavelength: float,
        mode: Mode,
        delta_neff: float = 1e-6) -> float:
    r"""
    Calculate the \(U\) parameter for a given fiber and mode.

    .. math::
        U = k_0 a \sqrt{n_{\text{core}}^2 - n_{\text{eff}}^2}

    Parameters
    ----------
    fiber : Fiber
        The fiber object.
    wavelength : float
        The wavelength of interest.
    mode : Mode
        The mode of interest.
    delta_neff : float, optional
        Convergence threshold for \(n_{\text{eff}}\), by default \(1 \times 10^{-6}\).

    Returns
    -------
    float
        The \(U\) parameter.
    """
    from PyFiberModes import solver
    assert fiber.n_layer == 2, "U-parameter can only be calculated for two-layer fibers."

    neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength)
    neff = neff_solver.solve(mode=mode, delta_neff=delta_neff)
    U, _, _ = neff_solver.get_U_W_V_parameter(neff=neff)
    return U


def get_effective_index(
        fiber,
        wavelength: float,
        mode: Mode,
        delta_neff: float = 1e-6) -> float:
    r"""
    Compute the effective refractive index, \(n_{\text{eff}}\), for a given mode.

    Parameters
    ----------
    fiber : Fiber
        The fiber object.
    wavelength : float
        The wavelength of interest.
    mode : Mode
        The mode of interest.
    delta_neff : float, optional
        Convergence threshold for \(n_{\text{eff}}\), by default \(1 \times 10^{-6}\).

    Returns
    -------
    float
        The effective refractive index, \(n_{\text{eff}}\).
    """
    from PyFiberModes import solver

    if fiber.n_layer == 2:
        neff_solver = solver.ssif.NeffSolver(fiber=fiber, wavelength=wavelength)
    else:
        neff_solver = solver.mlsif.NeffSolver(fiber=fiber, wavelength=wavelength)

    return neff_solver.solve(mode=mode, delta_neff=delta_neff)


def get_mode_cutoff_v0(
        fiber,
        wavelength: float,
        mode: Mode) -> float:
    r"""
    Compute the cutoff V-number, \(V_0\), for a mode.

    Parameters
    ----------
    fiber : Fiber
        The fiber object.
    wavelength : float
        The wavelength of interest.
    mode : Mode
        The mode of interest.

    Returns
    -------
    float
        The cutoff V-number, \(V_0\).
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
    Compute the radial field components of a mode in cylindrical coordinates.

    The tuple structure contains:
    - \((E_r, E_\phi, E_z)\): Electric field components.
    - \((H_r, H_\phi, H_z)\): Magnetic field components.

    Parameters
    ----------
    fiber : Fiber
        The fiber object.
    mode : Mode
        The mode of interest.
    wavelength : float
        The wavelength of interest.
    radius : float
        The radial position.

    Returns
    -------
    tuple
        Radial field components as CylindricalCoordinates.
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
