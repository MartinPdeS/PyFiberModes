"""Mode discovery and parameter-sweep utilities."""

from dataclasses import dataclass
from copy import deepcopy
from typing import Callable, Iterable, Sequence

import numpy as np

from PyFiberModes.mode import Mode


@dataclass(frozen=True)
class ModeSweepResult:
    """Dense, dependency-free result returned by :func:`sweep_modes`.

    Rows correspond to ``parameters`` and columns to ``modes``.  Missing or
    unguided solutions are represented by NaN.
    """

    parameter_name: str
    parameters: np.ndarray
    modes: tuple[Mode, ...]
    values: dict[str, np.ndarray]

    def __getitem__(self, metric: str) -> np.ndarray:
        return self.values[metric]

    def for_mode(self, mode: Mode) -> dict[str, np.ndarray]:
        index = self.modes.index(mode)
        return {name: value[:, index] for name, value in self.values.items()}

    def as_records(self) -> list[dict]:
        return [
            {self.parameter_name: float(parameter), "mode": mode, **{
                name: float(array[i, j]) for name, array in self.values.items()
            }}
            for i, parameter in enumerate(self.parameters)
            for j, mode in enumerate(self.modes)
        ]


def candidate_modes(
    families: Sequence[str] = ("LP",), max_nu: int = 6, max_m: int = 6
) -> tuple[Mode, ...]:
    """Generate deterministic candidate modes for automatic discovery."""
    modes = []
    for family in families:
        nus = (0,) if family in ("TE", "TM") else range(max_nu + 1)
        for nu in nus:
            if family in ("HE", "EH") and nu == 0:
                continue
            modes.extend(Mode(family, nu, m) for m in range(1, max_m + 1))
    return tuple(modes)


def find_modes(fiber, families=("LP",), max_nu=6, max_m=6) -> tuple[Mode, ...]:
    """Return all guided modes among a bounded set of candidates."""
    found = []
    for mode in candidate_modes(families, max_nu, max_m):
        try:
            neff = fiber.get_effective_index(mode)
            guided = fiber.last_layer.refractive_index < neff <= fiber.maximum_index
            if np.isfinite(neff) and guided:
                found.append(mode)
        except (ArithmeticError, AssertionError, RuntimeError, ValueError):
            continue
    return tuple(found)


def sweep_modes(
    fiber,
    parameters: Iterable[float],
    *,
    modes: Sequence[Mode] | None = None,
    parameter_name: str = "wavelength",
    setter: Callable | None = None,
    metrics: Sequence[str] = ("effective_index", "propagation_constant"),
    families=("LP",),
    max_nu=6,
    max_m=6,
) -> ModeSweepResult:
    """Track modes over wavelength or any parameter handled by ``setter``."""
    parameters = np.asarray(tuple(parameters), dtype=float)
    if parameters.ndim != 1 or not len(parameters):
        raise ValueError("parameters must be a non-empty one-dimensional sequence")
    working_fiber = deepcopy(fiber)
    if setter is None:
        if parameter_name != "wavelength":
            raise ValueError("a setter is required for non-wavelength sweeps")

        def setter(current_fiber, value):
            current_fiber.update_wavelength(value)
    if modes is None:
        discovered = []
        for parameter in parameters:
            setter(working_fiber, float(parameter))
            for mode in find_modes(working_fiber, families, max_nu, max_m):
                if mode not in discovered:
                    discovered.append(mode)
        modes = discovered
    modes = tuple(modes)
    values = {metric: np.full((len(parameters), len(modes)), np.nan) for metric in metrics}
    for i, parameter in enumerate(parameters):
        setter(working_fiber, float(parameter))
        for j, mode in enumerate(modes):
            for metric in metrics:
                try:
                    value = getattr(working_fiber, f"get_{metric}")(mode)
                    values[metric][i, j] = value if np.isfinite(value) else np.nan
                except (ArithmeticError, AssertionError, RuntimeError, ValueError):
                    pass
    return ModeSweepResult(parameter_name, parameters, modes, values)
