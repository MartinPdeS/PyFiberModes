from dataclasses import dataclass

import numpy as np
import pytest

from PyFiberModes import (
    ConvergenceError,
    Family,
    Mode,
    ValidationError,
)
from PyFiberModes.analysis import candidate_modes, find_modes, sweep_modes
from PyFiberModes.optimization import DesignParameter, optimize_fiber
from PyFiberModes.propagation import CoupledModeSystem, coupling_matrix, overlap


@dataclass
class FakeLayer:
    refractive_index: float
    radius_in: float = 0.0
    radius_out: float = 1.0


class FakeFiber:
    def __init__(self):
        self.wavelength = 1.0
        self.layers = [FakeLayer(2.0), FakeLayer(1.0, 1.0, np.inf)]

    @property
    def last_layer(self):
        return self.layers[-1]

    @property
    def maximum_index(self):
        return max(layer.refractive_index for layer in self.layers)

    def update_wavelength(self, wavelength):
        self.wavelength = wavelength

    def get_effective_index(self, mode):
        return 1.5 + 0.01 * mode.nu - 0.1 * self.wavelength

    def get_propagation_constant(self, mode):
        return 2 * np.pi * self.get_effective_index(mode) / self.wavelength


def test_mode_sweep_is_structured_and_restores_fiber():
    fiber = FakeFiber()
    mode = Mode("LP", 0, 1)
    result = sweep_modes(fiber, [1.0, 2.0], modes=[mode])
    assert result["effective_index"].shape == (2, 1)
    assert result.for_mode(mode)["effective_index"][0] == pytest.approx(1.4)
    assert fiber.wavelength == 1.0


def test_candidate_modes_obey_family_rules():
    modes = candidate_modes(("LP", "HE", "TE"), max_nu=2, max_m=2)
    assert Mode("LP", 0, 1) in modes
    assert Mode("HE", 0, 1) not in modes
    assert {mode.nu for mode in modes if mode.family == "TE"} == {0}


def test_mode_accepts_family_enum_and_validates_orders():
    assert Mode(Family.LP, 0, 1) == Mode("LP", 0, 1)
    with pytest.raises(ValidationError, match="unknown mode family"):
        Mode("invalid", 0, 1)
    with pytest.raises(ValidationError, match="nu must"):
        Mode("LP", -1, 1)
    with pytest.raises(ValidationError, match="m must"):
        Mode("LP", 0, 0)


def test_find_modes_filters_unguided_and_failed_candidates():
    fiber = FakeFiber()

    def effective_index(mode):
        if mode.m == 2:
            raise ConvergenceError("solver failed")
        return 1.2 if mode.nu == 0 else 0.9

    fiber.get_effective_index = effective_index
    assert find_modes(fiber, max_nu=1, max_m=2) == (Mode("LP", 0, 1),)


def test_candidate_mode_bounds_are_validated():
    with pytest.raises(ValidationError, match="max_nu"):
        candidate_modes(max_nu=-1)
    with pytest.raises(ValidationError, match="max_nu"):
        candidate_modes(max_m=0)


def test_sweep_validates_parameters_and_custom_setter_is_non_mutating():
    fiber = FakeFiber()
    with pytest.raises(ValueError, match="non-empty"):
        sweep_modes(fiber, [], modes=[])
    with pytest.raises(ValueError, match="setter"):
        sweep_modes(fiber, [1], modes=[], parameter_name="core_index")

    result = sweep_modes(
        fiber,
        [1.7, 1.8],
        modes=[Mode("LP", 0, 1)],
        parameter_name="core_index",
        setter=lambda candidate, value: setattr(candidate.layers[0], "refractive_index", value),
        metrics=("effective_index",),
    )
    assert len(result.as_records()) == 2
    assert fiber.layers[0].refractive_index == 2.0


def test_inverse_design_does_not_mutate_source():
    fiber = FakeFiber()
    parameter = DesignParameter.layer_index(0, (1.1, 2.0))
    result = optimize_fiber(
        fiber,
        [parameter],
        lambda candidate: (candidate.layers[0].refractive_index - 1.4) ** 2,
        method="L-BFGS-B",
    )
    assert result.parameters[parameter.name] == pytest.approx(1.4, abs=1e-4)
    assert fiber.layers[0].refractive_index == 2.0


def test_design_parameter_radius_updates_both_interfaces():
    fiber = FakeFiber()
    parameter = DesignParameter.layer_radius(0, (0.5, 1.5))
    parameter.setter(fiber, 1.25)
    assert fiber.layers[0].radius_out == 1.25
    assert fiber.layers[1].radius_in == 1.25


def test_inverse_design_validates_parameters_and_enforces_constraint():
    fiber = FakeFiber()
    with pytest.raises(ValueError, match="at least one"):
        optimize_fiber(fiber, [], lambda candidate: 0)
    parameter = DesignParameter.layer_index(0, (1.1, 2.0))
    result = optimize_fiber(
        fiber,
        [parameter],
        lambda candidate: (candidate.layers[0].refractive_index - 1.2) ** 2,
        constraints=(lambda candidate: candidate.layers[0].refractive_index - 1.6,),
        method="L-BFGS-B",
    )
    assert result.parameters[parameter.name] >= 1.59


def test_lossless_coupled_mode_propagation_conserves_power():
    z = np.linspace(0, 5, 100)
    result = CoupledModeSystem([0, 0], [[0, 0.2], [0.2, 0]]).propagate([1, 0], z)
    assert result.amplitudes.shape == (100, 2)
    assert np.allclose(result.powers.sum(axis=1), 1, atol=1e-6)


def test_longitudinally_varying_coupling_and_input_validation():
    system = CoupledModeSystem([0, 0], lambda z: [[0, 0.1 * z], [0.1 * z, 0]])
    result = system.propagate([1, 0], np.linspace(0, 1, 10))
    assert result.powers[-1, 1] > 0
    with pytest.raises(ValueError, match="same length"):
        system.propagate([1], [0, 1])


class FakeCoordinates:
    x = np.zeros((2, 2))
    dx = 0.5
    dy = 0.5


class FakeField:
    def __init__(self, ex, ey=None, ez=None):
        self._components = (
            np.asarray(ex, dtype=complex),
            np.zeros((2, 2), dtype=complex) if ey is None else np.asarray(ey, dtype=complex),
            np.zeros((2, 2), dtype=complex) if ez is None else np.asarray(ez, dtype=complex),
        )
        self.cartesian_coordinates = FakeCoordinates()

    def Ex(self):
        return self._components[0]

    def Ey(self):
        return self._components[1]

    def Ez(self):
        return self._components[2]

    def get_integrale_square(self, array):
        return float(np.sum(np.abs(array) ** 2) * 0.25)


def test_vector_overlap_normalization_and_grid_validation():
    first = FakeField([[1, 1], [1, 1]])
    second = FakeField([[1j, 1j], [1j, 1j]])
    assert overlap(first, first) == pytest.approx(1)
    assert abs(overlap(first, second)) == pytest.approx(1)
    second.cartesian_coordinates.x = np.zeros((3, 3))
    with pytest.raises(ValueError, match="same grid"):
        overlap(first, second)


def test_coupling_matrix_is_hermitian_and_has_zero_diagonal():
    first = FakeField([[1, 0], [0, 1]])
    second = FakeField([[1j, 0], [0, 1j]])
    matrix = coupling_matrix([first, second])
    assert np.allclose(matrix, matrix.conjugate().T)
    assert np.allclose(np.diag(matrix), 0)
