"""Constrained inverse design for layered fibers."""

from dataclasses import dataclass
from typing import Callable, Sequence
from copy import deepcopy

import numpy as np
from scipy.optimize import differential_evolution, minimize


@dataclass(frozen=True)
class DesignParameter:
    """A bounded scalar fiber parameter."""
    name: str
    bounds: tuple[float, float]
    setter: Callable[[object, float], None]

    @classmethod
    def layer_radius(cls, layer: int, bounds: tuple[float, float]):
        def set_radius(fiber, value):
            fiber.layers[layer].radius_out = value
            fiber.layers[layer + 1].radius_in = value
        return cls(f"layer_{layer}_radius", bounds, set_radius)

    @classmethod
    def layer_index(cls, layer: int, bounds: tuple[float, float]):
        return cls(
            f"layer_{layer}_index", bounds,
            lambda fiber, value: setattr(fiber.layers[layer], "refractive_index", value),
        )


@dataclass(frozen=True)
class DesignResult:
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

    Each constraint must be non-negative when satisfied. The input fiber is
    never modified; the optimized copy is returned in :class:`DesignResult`.
    """
    if not parameters:
        raise ValueError("at least one design parameter is required")
    template = deepcopy(fiber)

    def build(values):
        candidate = deepcopy(template)
        for parameter, value in zip(parameters, values):
            parameter.setter(candidate, float(value))
        return candidate

    def loss(values):
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
