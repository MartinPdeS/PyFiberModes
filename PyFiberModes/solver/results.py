"""Structured numerical results returned by solver service APIs."""

from dataclasses import dataclass, field
from typing import Generic, TypeVar

T = TypeVar("T")


@dataclass(frozen=True, slots=True)
class SolverResult(Generic[T]):
    """Describe a numerical solver outcome without ambiguous sentinel values.

    Parameters
    ----------
    value : T or None
        Computed value, or ``None`` when unsuccessful.
    converged : bool
        Whether a valid solution was obtained.
    residual : float or None, optional
        Residual of the characteristic equation when available.
    iterations : int or None, optional
        Number of iterations when reported by the backend.
    bracket : tuple of float or None, optional
        Root bracket used by the backend.
    message : str, optional
        Human-readable status information.
    metadata : dict, optional
        Additional backend-specific diagnostics.
    """

    value: T | None
    converged: bool
    residual: float | None = None
    iterations: int | None = None
    bracket: tuple[float, float] | None = None
    message: str = ""
    metadata: dict[str, object] = field(default_factory=dict)

    def unwrap(self) -> T:
        """Return the value or raise when the solve failed.

        Returns
        -------
        T
            Successfully computed value.

        Raises
        ------
        ConvergenceError
            If ``converged`` is false or no value is available.
        """
        from PyFiberModes.exceptions import ConvergenceError

        if not self.converged or self.value is None:
            raise ConvergenceError(self.message or "solver did not converge")
        return self.value
