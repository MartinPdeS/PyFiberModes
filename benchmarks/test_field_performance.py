"""Microbenchmarks for modal caching and vectorized field evaluation."""

import numpy as np

from PyFiberModes import HE11, LP01
from PyFiberModes.fiber import load_fiber


def test_cached_effective_index(benchmark):
    """Benchmark the cached effective-index public API."""
    fiber = load_fiber("SMF28", wavelength=1550e-9)
    fiber.get_effective_index(LP01)
    benchmark(fiber.get_effective_index, LP01)


def test_vectorized_index_profile(benchmark):
    """Benchmark vectorized layer lookup over a dense radial grid."""
    fiber = load_fiber("SMF28", wavelength=1550e-9)
    radius = np.linspace(0, 100e-6, 100_000)
    benchmark(fiber.get_index_profile, radius)


def test_batched_field_components(benchmark):
    """Benchmark one shared radial pass for all Cartesian components."""
    fiber = load_fiber("SMF28", wavelength=1550e-9)

    def evaluate():
        """Build a fresh field and evaluate six components together."""
        field = fiber.get_mode_field(HE11, limit=12e-6, n_point=101)
        return field.get_components()

    benchmark(evaluate)
