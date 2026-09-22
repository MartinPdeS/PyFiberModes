"""Numerical reference and convergence tests for physical outputs."""

import numpy as np
import pytest

from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber


@pytest.mark.parametrize(
    ("wavelength", "reference"),
    [(1310e-9, 1.4509131222977205), (1550e-9, 1.447490504877299)],
)
def test_smf28_lp01_effective_index_reference(wavelength, reference):
    """Pin effective index values for the bundled SMF-28 definition."""
    fiber = load_fiber("SMF28", wavelength=wavelength)
    neff = fiber.get_effective_index(LP01)
    assert neff == pytest.approx(reference, abs=2e-10)
    assert fiber.last_layer.refractive_index < neff < fiber.maximum_index


def test_effective_area_converges_with_grid_resolution():
    """Require the LP01 effective area to converge under grid refinement."""
    fiber = load_fiber("SMF28", wavelength=1550e-9)
    areas = [
        fiber.get_mode_field(LP01, limit=20e-6, n_point=size).get_effective_area()
        for size in (81, 121, 161)
    ]
    assert np.ptp(areas) / np.mean(areas) < 2e-3


def test_modal_and_cutoff_solutions_are_cached():
    """Avoid repeated characteristic-root solves for identical state."""
    fiber = load_fiber("SMF28", wavelength=1550e-9)
    fiber.get_effective_index(LP01)
    fiber.get_mode_cutoff_v0(LP01)
    initial_size = fiber.analysis.cache_size
    assert initial_size >= 2
    fiber.get_effective_index(LP01)
    fiber.get_mode_cutoff_v0(LP01)
    assert fiber.analysis.cache_size == initial_size
