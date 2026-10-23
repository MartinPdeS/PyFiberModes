"""Shared TypedUnit registry and conversions at the numerical SI boundary."""

import numpy as np
from pint import Quantity
from TypedUnit import (
    Angle,
    Area,
    Dimensionless,
    Frequency,
    Length,
    RefractiveIndex,
    Time,
    Velocity,
    ureg,
)

__all__ = [
    "ureg",
    "Quantity",
    "Angle",
    "Area",
    "Dimensionless",
    "Frequency",
    "Length",
    "RefractiveIndex",
    "Time",
    "Velocity",
]


def magnitude(value, unit, name):
    """Require a Pint quantity and return its magnitude in the requested unit.

    Dimensional checks use TypedUnit's shared registry. Bare numerical values
    are deliberately rejected, including zero and infinity.
    """
    if not isinstance(value, Quantity):
        raise TypeError(f"{name} must be a pint Quantity with units compatible with {unit}")
    return value.to(unit).magnitude


def quantity(value, unit):
    """Create a quantity in the same registry used by PyMieSim."""
    return ureg.Quantity(value, unit)


def scalar(value, unit, name):
    """Convert a scalar quantity, rejecting array-valued scalar parameters."""
    result = magnitude(value, unit, name)
    if np.ndim(result):
        raise ValueError(f"{name} must be scalar")
    return float(result)


MODAL_UNITS = {
    "effective_index": "dimensionless",
    "normalized_beta": "dimensionless",
    "mode_cutoff_v0": "dimensionless",
    "mode_cutoff_wavelength": "meter",
    "propagation_constant": "1 / meter",
    "phase_velocity": "meter / second",
    "group_index": "dimensionless",
    "group_velocity": "meter / second",
    "group_velocity_dispersion": "second ** 2 / meter",
    "dispersion": "picosecond / nanometer / kilometer",
    "S_parameter": "picosecond / nanometer ** 2 / kilometer",
}
