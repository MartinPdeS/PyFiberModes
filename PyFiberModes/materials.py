"""Pluggable refractive-index models used by fiber configuration loading."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class RefractiveIndexModel(Protocol):
    """Provide a refractive index at a vacuum wavelength."""

    def refractive_index(self, wavelength: float) -> float:
        """Return the refractive index at ``wavelength`` in meters."""
        ...


@dataclass(frozen=True, slots=True)
class ConstantIndex:
    """Represent a wavelength-independent refractive index."""

    value: float

    def refractive_index(self, wavelength: float) -> float:
        """Return the constant index after validating the wavelength."""
        if wavelength <= 0:
            raise ValueError("wavelength must be positive")
        return float(self.value)


@dataclass(frozen=True, slots=True)
class FusedSilica:
    """Evaluate fused silica with the Malitson Sellmeier coefficients."""

    def refractive_index(self, wavelength: float) -> float:
        """Return fused-silica index for a wavelength expressed in meters."""
        if wavelength <= 0:
            raise ValueError("wavelength must be positive")
        wavelength_micrometers = wavelength * 1e6
        wavelength_squared = wavelength_micrometers**2
        index_squared = 1.0
        for coefficient, resonance in (
            (0.6961663, 0.0684043),
            (0.4079426, 0.1162414),
            (0.8974794, 9.896161),
        ):
            index_squared += coefficient * wavelength_squared / (
                wavelength_squared - resonance**2
            )
        return float(np.sqrt(index_squared))


@dataclass(frozen=True, slots=True)
class CallableIndex:
    """Adapt a plain wavelength-to-index callable to the material protocol."""

    function: Callable[[float], float]

    def refractive_index(self, wavelength: float) -> float:
        """Evaluate the wrapped callable."""
        return float(self.function(wavelength))


@dataclass(slots=True)
class MaterialRegistry:
    """Resolve human-readable material names to refractive-index models."""

    models: dict[str, RefractiveIndexModel] = field(default_factory=dict)

    def __init__(
        self,
        models: Mapping[str, RefractiveIndexModel | Callable[[float], float]] | None = None,
    ) -> None:
        """Create a registry with built-in fused-silica aliases."""
        silica = FusedSilica()
        self.models = {"fused_silica": silica, "silica": silica}
        if models:
            for name, model in models.items():
                self.register(name, model)

    def register(
        self,
        name: str,
        model: RefractiveIndexModel | Callable[[float], float],
    ) -> None:
        """Register a model or callable under a normalized material name."""
        if not name.strip():
            raise ValueError("material name must not be empty")
        adapted = model if isinstance(model, RefractiveIndexModel) else CallableIndex(model)
        self.models[name.casefold()] = adapted

    def resolve(self, name: str) -> RefractiveIndexModel:
        """Return the model registered for ``name``."""
        try:
            return self.models[name.casefold()]
        except KeyError as error:
            available = ", ".join(sorted(self.models))
            raise KeyError(f"unknown material {name!r}; available materials: {available}") from error

    def refractive_index(self, name: str, wavelength: float) -> float:
        """Evaluate a named material at a vacuum wavelength."""
        return self.resolve(name).refractive_index(wavelength)


DEFAULT_MATERIALS = MaterialRegistry()


__all__ = [
    "CallableIndex",
    "ConstantIndex",
    "DEFAULT_MATERIALS",
    "FusedSilica",
    "MaterialRegistry",
    "RefractiveIndexModel",
]
