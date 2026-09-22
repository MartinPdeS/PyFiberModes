"""Focused orchestration services used by the :class:`Fiber` façade."""

from __future__ import annotations

import numpy as np

from PyFiberModes.exceptions import ConvergenceError
from PyFiberModes.models import SolverSettings
from PyFiberModes.solver.results import SolverResult


class ModalAnalysis:
    """Solve and cache modal properties for one fiber.

    Parameters
    ----------
    fiber : Fiber
        Fiber model owned by this service.
    settings : SolverSettings, optional
        Numerical settings shared by solves.
    """

    def __init__(self, fiber, settings: SolverSettings | None = None):
        """Initialize an empty geometry-aware result cache."""
        self.fiber = fiber
        self.settings = settings or SolverSettings()
        self._cache: dict[tuple, SolverResult[float]] = {}

    def clear(self) -> None:
        """Discard all cached modal solutions."""
        self._cache.clear()

    def _key(self, quantity: str, mode) -> tuple:
        """Build a cache key including every mutable physical input."""
        return quantity, mode, self.fiber.geometry_signature

    def effective_index_result(self, mode) -> SolverResult[float]:
        """Return a diagnostic effective-index result for a mode.

        Parameters
        ----------
        mode : Mode
            Mode to solve.

        Returns
        -------
        SolverResult[float]
            Cached value and convergence status.
        """
        key = self._key("effective_index", mode)
        if key in self._cache:
            return self._cache[key]

        from PyFiberModes.fundamentals import get_effective_index

        try:
            value = float(get_effective_index(
                fiber=self.fiber,
                wavelength=self.fiber.wavelength,
                mode=mode,
                delta_neff=self.settings.effective_index_step,
            ))
            converged = bool(np.isfinite(value))
            result = SolverResult(
                value=value if converged else None,
                converged=converged,
                message="converged" if converged else f"no guided solution for {mode}",
            )
        except (ArithmeticError, RuntimeError, ValueError) as error:
            result = SolverResult(value=None, converged=False, message=str(error))

        self._cache[key] = result
        if self.settings.raise_on_failure and not result.converged:
            raise ConvergenceError(result.message)
        return result

    def effective_index(self, mode) -> float:
        """Return an effective index through the scalar convenience API."""
        result = self.effective_index_result(mode)
        return result.value if result.converged else float("nan")

    def cutoff_result(self, mode) -> SolverResult[float]:
        """Return a cached normalized-cutoff result.

        Parameters
        ----------
        mode : Mode
            Mode whose cutoff is requested.

        Returns
        -------
        SolverResult[float]
            Normalized cutoff and explicit status.
        """
        key = self._key("cutoff", mode)
        if key not in self._cache:
            from PyFiberModes.fundamentals import get_mode_cutoff_v0

            try:
                value = float(get_mode_cutoff_v0(
                    fiber=self.fiber,
                    wavelength=self.fiber.wavelength,
                    mode=mode,
                ))
                converged = not np.isnan(value)
                self._cache[key] = SolverResult(
                    value=value if converged else None,
                    converged=converged,
                    message="converged" if converged else f"no cutoff solution for {mode}",
                )
            except (ArithmeticError, RuntimeError, ValueError) as error:
                self._cache[key] = SolverResult(
                    value=None, converged=False, message=str(error)
                )
        return self._cache[key]

    @property
    def cache_size(self) -> int:
        """Return the number of cached modal quantities."""
        return len(self._cache)


class FieldAnalysis:
    """Construct sampled and radial fields for one fiber.

    Parameters
    ----------
    fiber : Fiber
        Fiber model owned by this service.
    """

    def __init__(self, fiber):
        """Store the owning fiber."""
        self.fiber = fiber

    def mode_field(self, mode, limit: float | None = None, n_point: int = 101):
        """Construct a two-dimensional field grid.

        Parameters
        ----------
        mode : Mode
            Mode to sample.
        limit : float, optional
            Half-width of the Cartesian domain.
        n_point : int, optional
            Samples along each Cartesian axis.

        Returns
        -------
        Field
            Sampled vector field.
        """
        from PyFiberModes.field import Field

        if n_point < 2:
            raise ValueError("n_point must be at least two")
        return Field(
            fiber=self.fiber,
            mode=mode,
            limit=self.fiber.radius * 5.5 if limit is None else limit,
            n_point=n_point,
        )
