"""Validated immutable configuration models for fibers and solvers."""

from dataclasses import dataclass
import math

from PyFiberModes.exceptions import ValidationError
from PyFiberModes.units import Length, RefractiveIndex, scalar


@dataclass(frozen=True, slots=True)
class LayerSpec:
    """Describe one immutable concentric step-index layer.

    Parameters
    ----------
    name : str
        Human-readable layer name.
    radius : TypedUnit.Length
        Outer radius in meters. The outermost layer may use infinity.
    refractive_index : TypedUnit.RefractiveIndex
        Positive finite refractive index.
    """

    name: str
    radius: Length
    refractive_index: RefractiveIndex

    def __post_init__(self) -> None:
        """Validate physical layer values."""
        if not self.name:
            raise ValidationError("layer name must not be empty")
        radius = scalar(self.radius, "meter", "radius")
        index = scalar(self.refractive_index, "dimensionless", "refractive_index")
        if math.isnan(radius) or radius <= 0:
            raise ValidationError("layer radius must be positive")
        if not math.isfinite(index) or index <= 0:
            raise ValidationError("refractive index must be positive and finite")


@dataclass(frozen=True, slots=True)
class SolverSettings:
    """Configure root searches consistently across solver implementations.

    Parameters
    ----------
    tolerance : float, optional
        Absolute root tolerance.
    max_iterations : int, optional
        Maximum root-finding iterations.
    effective_index_step : float, optional
        Maximum effective-index bracketing step.
    raise_on_failure : bool, optional
        Raise :class:`ConvergenceError` instead of returning a failed result.
    """

    tolerance: float = 1e-8
    max_iterations: int = 100
    effective_index_step: float = 1e-6
    raise_on_failure: bool = False

    def __post_init__(self) -> None:
        """Validate solver tolerances and iteration limits."""
        if self.tolerance <= 0 or self.effective_index_step <= 0:
            raise ValidationError("solver tolerances must be positive")
        if self.max_iterations < 1:
            raise ValidationError("max_iterations must be at least one")
