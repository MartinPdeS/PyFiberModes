#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy


from PyFiberModes.fiber import load_fiber
from PyFiberModes import FiberFactory


def test_fiber_factory():
    factory = FiberFactory(wavelength=1550e-9)

    factory.add_layer(
        name="core",
        radius=4e-6,
        index=numpy.linspace(1.464, 1.494, 10)
    )

    factory.add_layer(name="cladding", index=1.4444)

    print(factory)


def test_fiber_loader():
    smf28 = load_fiber(fiber_name='SMF28', wavelength=1310e-9)

    print(smf28)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])

# -
