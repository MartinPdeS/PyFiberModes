"""Analytical lossless coupled-mode evolution."""

import numpy as np
import pytest
from PyFiberModes import CoupledModeSystem
from PyFiberModes.units import ureg


def test_lossless_coupled_mode_propagation_conserves_power():
    z = np.linspace(0, 5, 100)
    result = CoupledModeSystem(np.array([0, 0]) / ureg.meter, np.array([[0, 0.2], [0.2, 0]]) / ureg.meter).propagate([1, 0], z * ureg.meter)
    assert result.amplitudes.shape == (100, 2)
    assert np.allclose(result.powers.sum(axis=1), 1, atol=1e-6)


@pytest.mark.parametrize("start, stop", [(2.0, 6.0), (6.0, 2.0)])
@pytest.mark.parametrize("varying", [False, True])
def test_phase_matched_mode_exchange(start, stop, varying):
    """Recover signed quadrature amplitudes for both propagation directions."""
    positions = np.linspace(start, stop, 13) * ureg.centimeter
    matrix = np.array([[0, 1], [1, 0]]) / ureg.centimeter

    def coupling(position):
        """Increase coupling linearly with displacement from launch."""
        return matrix * (1 + (position.to("centimeter").magnitude - start) / 10)

    system = CoupledModeSystem(np.zeros(2) / ureg.meter, coupling if varying else matrix)
    result = system.propagate([1, 0], positions)
    displacement = positions.magnitude - start
    phase = displacement + displacement**2 / 20 if varying else displacement
    expected = np.column_stack((np.cos(phase), -1j * np.sin(phase)))
    np.testing.assert_allclose(result.amplitudes.magnitude, expected, atol=3e-8, rtol=0)
