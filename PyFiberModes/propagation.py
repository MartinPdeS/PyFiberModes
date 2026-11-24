"""Mode overlaps and coupled-mode longitudinal propagation."""

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.integrate import solve_ivp
from PyFiberModes.units import magnitude, quantity


def _require_common_grid(*, reference, other):
    """Require equal dimensions and physical sampling intervals on centered grids."""
    first = reference.cartesian_coordinates
    second = other.cartesian_coordinates
    if first.x.shape != second.x.shape or not np.allclose(
        (first.dx, first.dy), (second.dx, second.dy), rtol=1e-12, atol=0
    ):
        raise ValueError("fields must use the same grid shape and spacing")


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
        If the two field grids have different shapes or spacings.
    """
    if hasattr(field_a, "_si"):
        result = overlap(field_a._si, field_b._si, normalized=normalized)
        return quantity(result, "dimensionless" if normalized else "volt ** 2")
    _require_common_grid(reference=field_a, other=field_b)
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
    if fields and hasattr(fields[0], "_si"):
        unit = perturbation.units if hasattr(perturbation, "units") else "dimensionless"
        weight = perturbation.magnitude if hasattr(perturbation, "units") else perturbation
        return quantity(coupling_matrix([field._si for field in fields], weight), unit)
    for field in fields[1:]:
        _require_common_grid(reference=fields[0], other=field)
    count = len(fields)
    matrix = np.zeros((count, count), dtype=complex)
    weight = 1.0 if perturbation is None else np.asarray(perturbation)
    for i in range(count):
        for j in range(i + 1, count):
            components_a = (fields[i].Ex(), fields[i].Ey(), fields[i].Ez())
            components_b = (fields[j].Ex(), fields[j].Ey(), fields[j].Ez())
            value = sum(
                np.sum(np.conjugate(a) * weight * b) for a, b in zip(components_a, components_b)
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
            Finite real inverse-length quantities, one per mode.
        coupling : array-like or callable
            Square inverse-length matrix, or a function mapping a length
            quantity to such a matrix. Mode order must match ``betas``.
        """
        constants = np.asarray(magnitude(betas, "1 / meter", "betas"))
        if constants.ndim != 1 or not constants.size or not np.all(np.isfinite(constants)):
            raise ValueError("betas must be a non-empty vector of finite propagation constants")
        if np.any(np.imag(constants)):
            raise ValueError("betas must be real")
        self.betas = quantity(np.asarray(constants.real, dtype=float), "1 / meter")
        self.coupling = coupling if callable(coupling) else quantity(
            self._coupling_matrix(value=coupling), "1 / meter"
        )

    def _coupling_matrix(self, *, value):
        """Return a finite square coupling matrix in inverse meters."""
        matrix = np.asarray(magnitude(value, "1 / meter", "coupling"), dtype=complex)
        if matrix.shape != (self.betas.size, self.betas.size):
            raise ValueError("coupling must have shape (n_modes, n_modes)")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("coupling must contain only finite values")
        return matrix

    def propagate(self, initial, z, **solve_options) -> PropagationResult:
        """Propagate initial modal amplitudes over longitudinal positions.

        Parameters
        ----------
        initial : array-like
            Complex modal amplitudes at ``z[0]``.
        z : array-like
            At least two finite length samples in strictly increasing or decreasing order.
        **solve_options
            Additional options forwarded to :func:`scipy.integrate.solve_ivp`.

        Returns
        -------
        PropagationResult
            Positions and propagated modal amplitudes.

        Raises
        ------
        ValueError
            If the position grid, amplitude vector, or coupling matrix is invalid.
        RuntimeError
            If numerical integration fails.
        """
        positions = np.asarray(magnitude(z, "meter", "z"))
        if np.any(np.imag(positions)):
            raise ValueError("z positions must be real")
        z = np.asarray(positions.real, dtype=float)
        if z.ndim != 1 or z.size < 2 or not np.all(np.isfinite(z)):
            raise ValueError("z must be a vector of at least two finite positions")
        if not (np.all(z[1:] > z[:-1]) or np.all(z[1:] < z[:-1])):
            raise ValueError("z positions must be strictly monotonic")
        initial = np.asarray(initial, dtype=complex)
        if initial.shape != self.betas.shape:
            raise ValueError("initial amplitudes and betas must have the same length")

        if not np.all(np.isfinite(initial)):
            raise ValueError("initial amplitudes must be finite")
        diagonal = np.diag(self.betas.magnitude)
        if callable(self.coupling):
            def derivative(position, amplitudes):
                """Evaluate the local coupling at the integrator's position."""
                local = self.coupling(quantity(position, "meter"))
                return -1j * (diagonal + self._coupling_matrix(value=local)) @ amplitudes
        else:
            generator = -1j * (diagonal + self._coupling_matrix(value=self.coupling))

            def derivative(position, amplitudes):
                """Apply the constant propagation operator."""
                return generator @ amplitudes

        solution = solve_ivp(
            derivative,
            (z[0], z[-1]),
            initial,
            t_eval=z,
            rtol=solve_options.pop("rtol", 1e-8),
            atol=solve_options.pop("atol", 1e-10),
            **solve_options,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        return PropagationResult(
            quantity(solution.t, "meter"), quantity(solution.y.T, "dimensionless")
        )
