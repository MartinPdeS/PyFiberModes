"""Regression tests for the typed solver and service architecture."""

import numpy as np
import pytest

from PyFiberModes import LP01
from PyFiberModes.coordinates import CylindricalCoordinates
from PyFiberModes.exceptions import ConvergenceError, ValidationError
from PyFiberModes.fiber import load_fiber
from PyFiberModes.models import LayerSpec, SolverSettings
from PyFiberModes.solver.protocols import EffectiveIndexSolver
from PyFiberModes.solver.results import SolverResult


def test_field_uses_the_canonical_coordinate_model():
    """Ensure the field module does not define a second coordinate class."""
    from PyFiberModes import field

    assert field.CylindricalCoordinates is CylindricalCoordinates


def test_validated_immutable_models():
    """Reject nonphysical input and prevent configuration mutation."""
    layer = LayerSpec("core", 4.1e-6, 1.45)
    assert layer.refractive_index == 1.45
    with pytest.raises(ValidationError):
        LayerSpec("core", -1, 1.45)
    with pytest.raises(ValidationError):
        SolverSettings(tolerance=0)
    with pytest.raises(Exception):
        layer.radius = 5e-6


def test_solver_result_has_explicit_failure_semantics():
    """Expose failures without relying on a NaN sentinel."""
    failed = SolverResult[float](None, False, message="not guided")
    with pytest.raises(ConvergenceError, match="not guided"):
        failed.unwrap()


def test_solver_protocol_is_structural():
    """Allow third-party solvers without inheritance coupling."""
    class Solver:
        def solve(self, mode, delta_neff):
            """Return a deterministic test solution."""
            return 1.4

    assert isinstance(Solver(), EffectiveIndexSolver)


def test_modal_service_cache_tracks_geometry(monkeypatch):
    """Reuse solutions until wavelength or geometry changes."""
    import PyFiberModes.fundamentals as fundamentals

    fiber = load_fiber("SMF28", wavelength=1550e-9)
    calls = []

    def solve(**kwargs):
        """Record calls to the analytical backend."""
        calls.append(kwargs)
        return SolverResult(value=1.45, converged=True, message="converged")

    monkeypatch.setattr(fundamentals, "get_effective_index_result", solve)
    assert fiber.get_effective_index(LP01) == fiber.get_effective_index(LP01)
    assert len(calls) == 1

    fiber.layers[0].refractive_index += 1e-4
    fiber.get_effective_index(LP01)
    assert len(calls) == 2


def test_vectorized_index_profile_matches_scalar_lookup():
    """Evaluate a complete radial profile without a Python radius loop."""
    fiber = load_fiber("SMF28", wavelength=1550e-9)
    radii = np.array([0.0, fiber.first_layer.radius_out / 2, 2 * fiber.radius])
    profile = fiber.get_index_profile(radii)
    expected = np.array([fiber.get_index_at_radius(radius) for radius in radii])
    assert np.array_equal(profile, expected)


def test_structured_solve_result_and_scalar_value_agree():
    """Keep the float convenience API consistent with solver diagnostics."""
    fiber = load_fiber("SMF28", wavelength=1550e-9)
    result = fiber.solve_effective_index(LP01)
    assert result.converged
    assert result.value == pytest.approx(fiber.get_effective_index(LP01))
