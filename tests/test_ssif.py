#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pytest
from PyFiberModes.fiber import load_fiber
from PyFiberModes.tools.utils import get_mode_beta
from PyFiberModes import HE11, TE01, TM01, EH11, LP01


def test_solver():
    fiber = load_fiber(fiber_name='SMF28', wavelength=1550e-9, add_air_layer=False)

    itr_list = numpy.linspace(1.0, 0.3, 10)

    get_mode_beta(fiber=fiber, mode_list=[HE11, TE01, TM01, EH11, LP01], itr_list=itr_list)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
