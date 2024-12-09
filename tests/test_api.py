#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy

from PyFiberModes import HE11, LP01, LP11, LP21, LP12
from PyFiberModes.fiber import load_fiber
from PyFiberModes import FiberFactory
from PyFiberModes.source import Source


function_list = [
    "get_dispersion",
    "get_effective_index",
    "get_normalized_beta",
    "get_phase_velocity",
    "get_group_index",
    "get_groupe_velocity",
    "get_S_parameter",
    "get_mode_field"
]

attribute_list = [
    "dispersion",
    "effective_index",
    "normalized_beta",
    "phase_velocity",
    "groupe_velocity",
]


@pytest.mark.parametrize('function_string', function_list, ids=function_list)
def test_get_attribute(function_string):
    factory = FiberFactory(wavelength=1550e-9)

    factory.add_layer(
        name="core",
        radius=4e-6,
        index=numpy.linspace(1.464, 1.494, 10)
    )

    factory.add_layer(name="cladding", index=1.4444)

    fiber = factory[0]

    for fiber in factory:
        function = getattr(fiber, function_string)

        _ = function(mode=HE11)


def test_print_data():
    source = Source(1310e-9)

    smf28 = load_fiber(fiber_name='SMF28', wavelength=source.wavelength)

    smf28.print_data(data_type_list=attribute_list, mode_list=[LP01, LP11, LP21, LP12])


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
