"""Mode discovery and parameter-sweep utilities."""

from dataclasses import dataclass
from copy import deepcopy
from typing import Callable, Iterable, Sequence

import numpy as np

from PyFiberModes.exceptions import ConvergenceError, ValidationError
from PyFiberModes.mode import Mode


@dataclass(frozen=True)
class ModeSweepResult:
    """Store modal quantities evaluated over a parameter sweep.

    Rows correspond to ``parameters`` and columns to ``modes``.  Missing or
    unguided solutions are represented by NaN.

    Parameters
    ----------
    parameter_name : str
        Name of the swept parameter.
    parameters : numpy.ndarray
        One-dimensional parameter values.
    modes : tuple of Mode
        Modes represented by the columns of every value array.
    values : dict of str to numpy.ndarray
        Modal metrics shaped ``(n_parameters, n_modes)``.
    """

    parameter_name: str
    parameters: np.ndarray
    modes: tuple[Mode, ...]
    values: dict[str, np.ndarray]

    def __getitem__(self, metric: str) -> np.ndarray:
        """Return the array associated with a metric name.

        Parameters
        ----------
        metric : str
            Metric key stored in :attr:`values`.

        Returns
        -------
        numpy.ndarray
            Metric values shaped ``(n_parameters, n_modes)``.
        """
        return self.values[metric]

    def for_mode(self, mode: Mode) -> dict[str, np.ndarray]:
        """Extract every metric for one mode.

        Parameters
        ----------
        mode : Mode
            Mode whose column is requested.

        Returns
        -------
        dict of str to numpy.ndarray
            Mapping from metric names to one-dimensional sweep values.
        """
        index = self.modes.index(mode)
        return {name: value[:, index] for name, value in self.values.items()}

    def as_records(self) -> list[dict[str, float | Mode]]:
        """Convert the dense result into row-oriented records.

        Returns
        -------
        list of dict
            One record for every parameter and mode combination.
        """
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
    """Generate a deterministic bounded set of mode candidates.

    Parameters
    ----------
    families : sequence of str, optional
        Mode families to enumerate.
    max_nu : int, optional
        Largest azimuthal order to include.
    max_m : int, optional
        Largest radial order to include.

    Returns
    -------
    tuple of Mode
        Candidates ordered by family, azimuthal order, and radial order.
    """
    if max_nu < 0 or max_m < 1:
        raise ValidationError("max_nu must be non-negative and max_m must be positive")
    modes = []
    for family in families:
        nus = (0,) if family in ("TE", "TM") else range(max_nu + 1)
        for nu in nus:
            if family in ("HE", "EH") and nu == 0:
                continue
            modes.extend(Mode(family, nu, m) for m in range(1, max_m + 1))
    return tuple(modes)


def find_modes(fiber, families=("LP",), max_nu=6, max_m=6) -> tuple[Mode, ...]:
    """Find guided modes supported by a fiber.

    Parameters
    ----------
    fiber : Fiber
        Fiber on which effective indices are evaluated.
    families : sequence of str, optional
        Mode families to search.
    max_nu : int, optional
        Largest azimuthal order to test.
    max_m : int, optional
        Largest radial order to test.

    Returns
    -------
    tuple of Mode
        Candidates with finite effective indices between the cladding and
        maximum material indices.
    """
    found = []
    for mode in candidate_modes(families, max_nu, max_m):
        try:
            effective_index = fiber.get_effective_index(mode)
            guided = (
                fiber.last_layer.refractive_index
                < effective_index
                <= fiber.maximum_index
            )
            if np.isfinite(effective_index) and guided:
                found.append(mode)
        except ConvergenceError:
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
    """Evaluate modal metrics over wavelength or another scalar parameter.

    Parameters
    ----------
    fiber : Fiber
        Fiber used as the immutable sweep template.
    parameters : iterable of float
        Ordered parameter values.
    modes : sequence of Mode, optional
        Modes to evaluate. When omitted, modes are discovered automatically.
    parameter_name : str, optional
        Name recorded in the returned result.
    setter : callable, optional
        Function accepting ``(fiber, value)`` for non-wavelength sweeps.
    metrics : sequence of str, optional
        Suffixes of ``Fiber.get_<metric>`` methods to evaluate.
    families : sequence of str, optional
        Families used during automatic discovery.
    max_nu : int, optional
        Maximum discovered azimuthal order.
    max_m : int, optional
        Maximum discovered radial order.

    Returns
    -------
    ModeSweepResult
        Dense arrays for each requested metric.

    Raises
    ------
    ValueError
        If ``parameters`` is empty or not one-dimensional, or when a custom
        parameter is requested without a setter.
    """
    parameters = np.asarray(tuple(parameters), dtype=float)
    if parameters.ndim != 1 or not len(parameters):
        raise ValueError("parameters must be a non-empty one-dimensional sequence")
    if setter is None:
        if parameter_name != "wavelength":
            raise ValueError("a setter is required for non-wavelength sweeps")

        def setter(current_fiber, value):
            """Update wavelength for the default sweep behavior."""
            current_fiber.update_wavelength(value)
    if modes is None:
        discovered = []
        for parameter in parameters:
            parameter_fiber = deepcopy(fiber)
            setter(parameter_fiber, float(parameter))
            for mode in find_modes(parameter_fiber, families, max_nu, max_m):
                if mode not in discovered:
                    discovered.append(mode)
        modes = discovered
    modes = tuple(modes)
    values = {metric: np.full((len(parameters), len(modes)), np.nan) for metric in metrics}
    parameter_fibers = tuple(deepcopy(fiber) for _ in parameters)
    for parameter_fiber, parameter in zip(parameter_fibers, parameters, strict=True):
        setter(parameter_fiber, float(parameter))

    for i, parameter_fiber in enumerate(parameter_fibers):
        metric_functions = {
            metric: getattr(parameter_fiber, f"get_{metric}") for metric in metrics
        }
        for metric, metric_function in metric_functions.items():
            for j, mode in enumerate(modes):
                try:
                    value = metric_function(mode)
                    values[metric][i, j] = value if np.isfinite(value) else np.nan
                except ConvergenceError:
                    pass
    return ModeSweepResult(parameter_name, parameters, modes, values)
