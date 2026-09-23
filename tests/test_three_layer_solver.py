#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pytest
from PyFiberModes.fiber import load_fiber
from PyFiberModes import HE11, TE01, TM01, EH11, LP01
from tests.helpers import get_mode_beta
from PyFiberModes.solver.multilayer import EffectiveIndexSolver


def test_three_layer_solver():
    fiber = load_fiber(fiber_name='SMF28', wavelength=1550e-9, add_air_layer=True)

    itr_list = numpy.linspace(1.0, 0.3, 2)

    get_mode_beta(fiber=fiber, mode_list=[LP01, HE11, TE01, TM01, EH11], itr_list=itr_list)


@pytest.mark.parametrize("mode", [TE01, TM01])
def test_multilayer_axisymmetric_fields_are_finite_in_every_region(mode):
    """Evaluate completed TE/TM branches across core, cladding, and air."""
    fiber = load_fiber(fiber_name="SMF28", wavelength=1550e-9, add_air_layer=True)
    solver = EffectiveIndexSolver(fiber=fiber, wavelength=fiber.wavelength)
    solution = solver.solve_result(mode=mode, delta_neff=1e-6)
    field_function = solver.get_TE_field if mode.family == "TE" else solver.get_TM_field

    assert solution.converged
    assert solution.value is not None
    for radius in (0.0, 2e-6, 5e-6, 70e-6):
        electric, magnetic = field_function(mode.nu, solution.value, radius)
        assert numpy.isfinite(electric).all()
        assert numpy.isfinite(magnetic).all()

    if mode.family == "TE":
        assert electric[2] == 0
    else:
        assert magnetic[2] == 0


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
