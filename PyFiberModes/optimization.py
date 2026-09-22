"""Constrained inverse design for layered fibers."""

from dataclasses import dataclass
from typing import Callable, Sequence
from copy import deepcopy

import numpy as np
from scipy.optimize import differential_evolution, minimize


@dataclass(frozen=True)
class DesignParameter:
    """Describe one bounded scalar fiber design variable.

    Parameters
    ----------
    name : str
        Stable name used in optimization results.
    bounds : tuple of float
        Inclusive lower and upper bounds.
    setter : callable
        Function accepting ``(fiber, value)`` and mutating a candidate fiber.
    """
    name: str
    bounds: tuple[float, float]
    setter: Callable[[object, float], None]

    @classmethod
    def layer_radius(cls, layer: int, bounds: tuple[float, float]):
        """Create a parameter controlling a layer interface radius.

        Parameters
        ----------
        layer : int
            Inner layer index whose outer radius is varied.
        bounds : tuple of float
            Lower and upper radius bounds in meters.

        Returns
        -------
        DesignParameter
            Radius design variable that keeps adjacent interfaces continuous.
        """
        def set_radius(fiber, value):
            """Apply a shared interface radius to adjacent layers."""
            fiber.layers[layer].radius_out = value
            fiber.layers[layer + 1].radius_in = value
        return cls(f"layer_{layer}_radius", bounds, set_radius)

    @classmethod
    def layer_index(cls, layer: int, bounds: tuple[float, float]):
        """Create a parameter controlling one layer refractive index.

        Parameters
        ----------
        layer : int
            Layer index to modify.
        bounds : tuple of float
            Lower and upper refractive-index bounds.

        Returns
        -------
        DesignParameter
            Refractive-index design variable.
        """
        return cls(
            f"layer_{layer}_index", bounds,
            lambda fiber, value: setattr(fiber.layers[layer], "refractive_index", value),
        )


@dataclass(frozen=True)
class DesignResult:
    """Store the outcome of a fiber optimization.

    Parameters
    ----------
    fiber : Fiber
        Optimized copy of the input fiber.
    parameters : dict of str to float
        Best value for each named design parameter.
    objective : float
        Final penalized objective value.
    success : bool
        Whether the numerical optimizer reported convergence.
    message : str
        Optimizer termination message.
    evaluations : int
        Number of objective evaluations.
    """
    fiber: object
    parameters: dict[str, float]
    objective: float
    success: bool
    message: str
    evaluations: int


def optimize_fiber(
    fiber,
    parameters: Sequence[DesignParameter],
    objective: Callable[[object], float],
    *,
    constraints: Sequence[Callable[[object], float]] = (),
    method: str = "differential_evolution",
    penalty: float = 1e12,
    options: dict | None = None,
) -> DesignResult:
    """Minimize an arbitrary fiber objective with inequality constraints.

    Parameters
    ----------
    fiber : Fiber
        Fiber used as an immutable design template.
    parameters : sequence of DesignParameter
        Variables exposed to the optimizer.
    objective : callable
        Scalar function evaluated on each candidate fiber.
    constraints : sequence of callable, optional
        Inequality functions that are non-negative when satisfied.
    method : str, optional
        ``"differential_evolution"`` or a method accepted by
        :func:`scipy.optimize.minimize`.
    penalty : float, optional
        Quadratic penalty multiplier for constraint violations.
    options : dict, optional
        Options forwarded to the selected SciPy optimizer.

    Returns
    -------
    DesignResult
        Optimized fiber and numerical termination information.

    Raises
    ------
    ValueError
        If no design parameters are supplied.

    Notes
    -----
    The input fiber is never modified.
    """
    if not parameters:
        raise ValueError("at least one design parameter is required")
    template = deepcopy(fiber)

    def build(values):
        """Construct a candidate fiber from optimizer coordinates."""
        candidate = deepcopy(template)
        for parameter, value in zip(parameters, values):
            parameter.setter(candidate, float(value))
        return candidate

    def loss(values):
        """Evaluate the penalized objective for optimizer coordinates."""
        candidate = build(values)
        violation = sum(max(0.0, -float(rule(candidate))) ** 2 for rule in constraints)
        value = float(objective(candidate)) + penalty * violation
        return value if np.isfinite(value) else np.finfo(float).max

    bounds = [parameter.bounds for parameter in parameters]
    options = dict(options or {})
    if method == "differential_evolution":
        raw = differential_evolution(loss, bounds=bounds, **options)
    else:
        x0 = np.asarray([np.mean(bound) for bound in bounds])
        raw = minimize(loss, x0=x0, bounds=bounds, method=method, options=options)
    best = build(raw.x)
    return DesignResult(
        best,
        {parameter.name: float(value) for parameter, value in zip(parameters, raw.x)},
        float(raw.fun), bool(raw.success), str(raw.message), int(raw.nfev),
    )
