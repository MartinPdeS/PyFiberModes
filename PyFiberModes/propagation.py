"""Mode overlaps and coupled-mode longitudinal propagation."""

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.integrate import solve_ivp


def overlap(field_a, field_b, *, normalized: bool = True) -> complex:
    """Calculate the electric-field overlap on a shared Cartesian grid.

    Parameters
    ----------
    field_a, field_b : Field
        Field objects sampled on grids with identical shapes and spacing.
    normalized : bool, optional
        Divide the inner product by both field norms when true.

    Returns
    -------
    complex
        Complex overlap integral or normalized overlap coefficient.

    Raises
    ------
    ValueError
        If the two field grids have different shapes.
    """
    if field_a.cartesian_coordinates.x.shape != field_b.cartesian_coordinates.x.shape:
        raise ValueError("fields must use the same grid shape")
    components_a = (field_a.Ex(), field_a.Ey(), field_a.Ez())
    components_b = (field_b.Ex(), field_b.Ey(), field_b.Ez())
    inner = sum(np.vdot(a, b) for a, b in zip(components_a, components_b))
    inner *= field_a.cartesian_coordinates.dx * field_a.cartesian_coordinates.dy
    if not normalized:
        return inner
    norm_a = sum(field_a.get_integrale_square(component) for component in components_a)
    norm_b = sum(field_b.get_integrale_square(component) for component in components_b)
    denominator = np.sqrt(norm_a * norm_b)
    return inner / denominator if denominator else 0j


def coupling_matrix(fields: Sequence, perturbation=None) -> np.ndarray:
    """Build a Hermitian coupling matrix from sampled mode fields.

    Parameters
    ----------
    fields : sequence of Field
        Mode fields defined on a common Cartesian grid.
    perturbation : scalar or numpy.ndarray, optional
        Spatial weighting applied to each overlap integrand.

    Returns
    -------
    numpy.ndarray
        Complex Hermitian matrix with a zero diagonal.
    """
    count = len(fields)
    matrix = np.zeros((count, count), dtype=complex)
    weight = 1.0 if perturbation is None else np.asarray(perturbation)
    for i in range(count):
        for j in range(i + 1, count):
            components_a = (fields[i].Ex(), fields[i].Ey(), fields[i].Ez())
            components_b = (fields[j].Ex(), fields[j].Ey(), fields[j].Ez())
            value = sum(
                np.sum(np.conjugate(a) * weight * b)
                for a, b in zip(components_a, components_b)
            )
            value *= fields[i].cartesian_coordinates.dx * fields[i].cartesian_coordinates.dy
            norm_a = sum(fields[i].get_integrale_square(a) for a in components_a)
            norm_b = sum(fields[j].get_integrale_square(b) for b in components_b)
            value /= np.sqrt(norm_a * norm_b)
            matrix[i, j], matrix[j, i] = value, np.conjugate(value)
    return matrix


@dataclass(frozen=True)
class PropagationResult:
    """Store longitudinal coupled-mode propagation results.

    Parameters
    ----------
    z : numpy.ndarray
        Longitudinal sample positions.
    amplitudes : numpy.ndarray
        Complex modal amplitudes shaped ``(n_positions, n_modes)``.
    """
    z: np.ndarray
    amplitudes: np.ndarray

    @property
    def powers(self):
        """Return modal powers at every longitudinal position.

        Returns
        -------
        numpy.ndarray
            Squared amplitude magnitudes with the same shape as
            :attr:`amplitudes`.
        """
        return np.abs(self.amplitudes) ** 2


class CoupledModeSystem:
    r"""Represent a longitudinal coupled-mode system.

    Parameters
    ----------
    betas : array-like
        Propagation constants for each mode.
    coupling : array-like or callable
        Constant coupling matrix or function ``coupling(z)``.

    Notes
    -----
    Amplitudes satisfy :math:`dA/dz = -i(\beta + \kappa)A`.
    """

    def __init__(self, betas, coupling):
        """Initialize propagation constants and modal coupling.

        Parameters
        ----------
        betas : array-like
            Propagation constants for each mode.
        coupling : array-like or callable
            Constant or position-dependent coupling matrix.
        """
        self.betas = np.asarray(betas, dtype=float)
        self.coupling = coupling

    def propagate(self, initial, z, **solve_options) -> PropagationResult:
        """Propagate initial modal amplitudes over longitudinal positions.

        Parameters
        ----------
        initial : array-like
            Complex modal amplitudes at ``z[0]``.
        z : array-like
            Ordered longitudinal evaluation positions.
        **solve_options
            Additional options forwarded to :func:`scipy.integrate.solve_ivp`.

        Returns
        -------
        PropagationResult
            Positions and propagated modal amplitudes.

        Raises
        ------
        ValueError
            If initial amplitudes and propagation constants differ in length.
        RuntimeError
            If numerical integration fails.
        """
        z = np.asarray(z, dtype=float)
        initial = np.asarray(initial, dtype=complex)
        if initial.shape != self.betas.shape:
            raise ValueError("initial amplitudes and betas must have the same length")

        def derivative(position, amplitudes):
            """Evaluate the coupled-mode ordinary differential equation."""
            coupling = self.coupling(position) if callable(self.coupling) else self.coupling
            hamiltonian = np.diag(self.betas) + np.asarray(coupling, dtype=complex)
            return -1j * hamiltonian @ amplitudes

        solution = solve_ivp(
            derivative, (z[0], z[-1]), initial, t_eval=z,
            rtol=solve_options.pop("rtol", 1e-8),
            atol=solve_options.pop("atol", 1e-10), **solve_options,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        return PropagationResult(solution.t, solution.y.T)
