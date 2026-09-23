"""Structural interfaces implemented by modal solver backends."""

from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from PyFiberModes.mode import Mode
from PyFiberModes.solver.results import SolverResult


@runtime_checkable
class EffectiveIndexSolver(Protocol):
    """Protocol for effective-index solver implementations."""

    def solve(self, mode: Mode, delta_neff: float) -> float:
        """Solve the effective index for ``mode``."""
        ...


@runtime_checkable
class DiagnosticEffectiveIndexSolver(EffectiveIndexSolver, Protocol):
    """Effective-index backend that exposes structured diagnostics."""

    def solve_result(self, mode: Mode, delta_neff: float) -> SolverResult[float]:
        """Solve ``mode`` and retain convergence diagnostics."""
        ...


@runtime_checkable
class CutoffSolver(Protocol):
    """Protocol for normalized-cutoff solver implementations."""

    def solve(self, mode: Mode) -> float:
        """Solve the normalized cutoff for ``mode``."""
        ...


@runtime_checkable
class RadialFieldSolver(Protocol):
    """Protocol for cylindrical radial-field solver implementations."""

    def get_LP_field(
        self, nu: int, neff: float, radius: float
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Evaluate an LP field at a radial position."""
        ...
