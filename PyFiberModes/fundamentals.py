"""Fundamental normalized-frequency and propagation relations."""

from __future__ import annotations

import numpy
import numpy as np
from typing import TYPE_CHECKING
from PyFiberModes.mode import Mode
from PyFiberModes.mode_instances import HE11, LP01
from scipy.constants import c

from PyFiberModes.coordinates import CylindricalCoordinates
from PyFiberModes.exceptions import ConvergenceError, UnsupportedGeometryError

if TYPE_CHECKING:
    from PyFiberModes.solver.results import SolverResult


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


def get_wavelength_from_V0(fiber: object, normalized_frequency: float) -> float:
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
    normalized_frequency : float
        The V-number.

    Returns
    -------
    float
        The wavelength corresponding to V$_0$.
    """
    NA = fiber.get_NA()
    last_layer = fiber.last_layer
    wavelength = 2 * np.pi / normalized_frequency * last_layer.radius_in * NA
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

    effective_index_solver = (
        solver.two_layer.EffectiveIndexSolver(fiber=fiber, wavelength=wavelength)
        if fiber.n_layer == 2
        else solver.multilayer.EffectiveIndexSolver(fiber=fiber, wavelength=wavelength)
    )

    effective_index = effective_index_solver.solve(mode=mode, delta_neff=delta_neff)
    return effective_index * (2 * numpy.pi / wavelength)


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
    if fiber.n_layer != 2:
        raise UnsupportedGeometryError(
            "U-parameter calculations require exactly two fiber layers"
        )

    effective_index_solver = solver.two_layer.EffectiveIndexSolver(
        fiber=fiber, wavelength=wavelength
    )
    effective_index = effective_index_solver.solve(mode=mode, delta_neff=delta_neff)
    U, _, _ = effective_index_solver.get_U_W_V_parameter(neff=effective_index)
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
        effective_index_solver = solver.two_layer.EffectiveIndexSolver(
            fiber=fiber, wavelength=wavelength
        )
    else:
        effective_index_solver = solver.multilayer.EffectiveIndexSolver(
            fiber=fiber, wavelength=wavelength
        )

    return effective_index_solver.solve(mode=mode, delta_neff=delta_neff)


def get_effective_index_result(
        fiber,
        wavelength: float,
        mode: Mode,
        delta_neff: float = 1e-6) -> "SolverResult[float]":
    """Compute an effective index with structured convergence diagnostics.

    Parameters
    ----------
    fiber : Fiber
        Fiber to solve.
    wavelength : float
        Vacuum wavelength in meters.
    mode : Mode
        Mode to solve.
    delta_neff : float, optional
        Maximum effective-index search step.

    Returns
    -------
    SolverResult[float]
        Effective index, residual, bracket, iteration count, and explanation.
    """
    from PyFiberModes import solver

    solver_class = (
        solver.two_layer.EffectiveIndexSolver
        if fiber.n_layer == 2
        else solver.multilayer.EffectiveIndexSolver
    )
    effective_index_solver = solver_class(fiber=fiber, wavelength=wavelength)
    return effective_index_solver.solve_result(mode=mode, delta_neff=delta_neff)


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
        case 2:
            cutoff_solver = solver.two_layer.CutoffSolver(fiber=fiber, wavelength=wavelength)
        case 3:
            cutoff_solver = solver.three_layer.CutoffSolver(fiber=fiber, wavelength=wavelength)
        case _:
            raise UnsupportedGeometryError(
                "Cutoff calculations are currently available for two- and three-layer fibers"
            )

    cutoff = cutoff_solver.solve(mode=mode)

    return cutoff


def get_radial_field(
        fiber,
        mode: Mode,
        wavelength: float,
        radius: float) -> tuple:
    r"""
    Compute the radial field components of a mode in cylindrical coordinates.

    The returned pair contains the cylindrical electric components
    :math:`(E_r, E_\phi, E_z)` and magnetic components
    :math:`(H_r, H_\phi, H_z)`.

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

    if fiber.n_layer == 2:
        effective_index_solver = solver.two_layer.EffectiveIndexSolver(
            fiber=fiber, wavelength=wavelength
        )
    else:
        effective_index_solver = solver.multilayer.EffectiveIndexSolver(
            fiber=fiber, wavelength=wavelength
        )

    neff = get_effective_index(
        fiber=fiber,
        wavelength=fiber.wavelength,
        mode=mode
    )
    if not numpy.isfinite(neff):
        raise ConvergenceError(f"no guided effective-index solution for {mode}")

    kwargs = dict(
        nu=mode.nu,
        neff=neff,
        radius=radius
    )

    match mode.family:
        case 'LP':
            (er, ephi, ez), (hr, hphi, hz) = effective_index_solver.get_LP_field(**kwargs)
        case 'TE':
            (er, ephi, ez), (hr, hphi, hz) = effective_index_solver.get_TE_field(**kwargs)
        case 'TM':
            (er, ephi, ez), (hr, hphi, hz) = effective_index_solver.get_TM_field(**kwargs)
        case 'EH':
            (er, ephi, ez), (hr, hphi, hz) = effective_index_solver.get_EH_field(**kwargs)
        case 'HE':
            (er, ephi, ez), (hr, hphi, hz) = effective_index_solver.get_HE_field(**kwargs)

    e_field = CylindricalCoordinates(rho=er, phi=ephi, z=ez)
    h_field = CylindricalCoordinates(rho=hr, phi=hphi, z=hz)

    return e_field, h_field
