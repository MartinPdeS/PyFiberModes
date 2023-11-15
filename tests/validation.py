#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

from PyFiberModes.fiber import get_fiber_from_delta_and_V0
from PyFiberModes import HE11, TE01, TM01, HE21, EH11, HE31, HE12


def test_validation_cutoff_wavelength():
    """ Validation data are from Table 3.6 of Jacques Bures """
    fiber = get_fiber_from_delta_and_V0(delta=0.3, V0=5, wavelength=1310e-9)

    mode_list = [HE11, TE01, TM01, HE21, EH11, HE31, HE12]
    val_list = [0, 2.405, 2.405, 2.853, 3.832, 4.342, 3.832]

    for mode, val in zip(mode_list, val_list):
        cutoff_V0 = fiber.get_cutoff_v0(mode=mode)
        assert numpy.isclose(cutoff_V0, val, atol=1e-3), f"Mode {mode} cutoff V0 do not match validation data."

# -
