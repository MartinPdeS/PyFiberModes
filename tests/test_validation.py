#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pytest
from PyFiberModes.fiber import get_fiber_from_delta_and_V0
from PyFiberModes.fiber import load_fiber
from PyFiberModes.tools.utils import get_mode_beta
from PyFiberModes import HE11, TE01, TM01, HE21, EH11, HE31, HE12, LP01


def test_validation_cutoff_wavelength():
    """ Validation data are from Table 3.6 of Jacques Bures """
    fiber = get_fiber_from_delta_and_V0(delta=0.3, V0=5, wavelength=1310e-9)

    mode_list = [HE11, TE01, TM01, HE21, EH11, HE31, HE12]
    val_list = [0, 2.405, 2.405, 2.853, 3.832, 4.342, 3.832]

    for mode, val in zip(mode_list, val_list):
        cutoff_V0 = fiber.get_mode_cutoff_v0(mode=mode)
        assert numpy.isclose(cutoff_V0, val, atol=1e-3), f"Mode {mode} cutoff V0 do not match validation data."


def test_2_3_layer_solver():
    kwargs = dict(fiber_name='SMF28', wavelength=1550e-9)

    fiber_2L = load_fiber(**kwargs, add_air_layer=False)
    fiber_3L = load_fiber(**kwargs, add_air_layer=True)

    itr_list = numpy.linspace(1.0, 0.3, 50)

    kwargs = dict(mode_list=[LP01], itr_list=itr_list)

    data_3L = get_mode_beta(fiber=fiber_3L, **kwargs)
    data_2L = get_mode_beta(fiber=fiber_2L, **kwargs)

    beta_3L = data_3L['LP01']
    beta_2L = data_2L['LP01']

    discrepencies = numpy.isclose(beta_3L, beta_2L, atol=1e-4, equal_nan=True)

    assert discrepencies.all()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])

# -
