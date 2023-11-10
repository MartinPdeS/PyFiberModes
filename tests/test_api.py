#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy


def test_fiber_factory():
    from PyFiberModes import FiberFactory

    factory = FiberFactory(wavelength=1550e-9)

    factory.add_layer(
        name="core",
        radius=4e-6,
        index=numpy.linspace(1.464, 1.494, 10)
    )

    factory.add_layer(name="cladding", index=1.4444)


attribute_function = [
    "get_dispersion",
    "get_effective_index",
    "get_normalized_beta",
    "get_phase_velocity",
    "get_group_index",
    "get_groupe_velocity",
    "get_S_parameter",
    "get_field"
]


@pytest.mark.parametrize('function_string', attribute_function, ids=attribute_function)
def test_get_attribute(function_string):
    from PyFiberModes import FiberFactory, HE11
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


if __name__ == '__main__':
    test_get_attribute(function_string='get_dispersion')
# -
