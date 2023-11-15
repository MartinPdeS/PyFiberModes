#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from itertools import product
from PyFiberModes.fiber import Fiber
from dataclasses import dataclass


@dataclass
class ProxyLayer():
    name: str
    radius: list
    index: list

    def __post_init__(self):
        self.name = [self.name]
        self.radius = numpy.atleast_1d(self.radius)
        self.index = numpy.atleast_1d(self.index)

    def get_generator(self):
        return product(self.name, self.radius, self.index)


class FiberFactory(object):
    """
    FiberFactory is used to instantiate a
    :py:class:`~PyFiberModes.fiber.fiber.Fiber` or a series of
    Fiber objects.

    It can read fiber definition from json file, and write it back.
    Convenient functions are available to set fiber parameters, and to
    iterate through fiber objects.

    All fibers build from a given factory share the same number
    of layers, the same kind of geometries, and the same
    materials. However, parameters can vary.

    Args:
        filename: Name of fiber file to load, or None to construct
                  empty Fiberfactory object.

    """
    def __init__(self, wavelength: float):
        self.layers_list = []
        self.neff_solver = None
        self.cutoff_solver = None
        self.wavelength = wavelength

    def add_layer(self,
            index: float,
            name: str = "",
            radius: float = 0) -> None:
        """
        Insert a new layer in the factory.

        :param      name:      Layer name.
        :type       name:      str
        :param      radius:    Radius of the layer (in meters).
        :type       radius:    float
        :param      kwargs:    The keywords arguments
        :type       kwargs:    dictionary
        """
        layer = ProxyLayer(
            name=name,
            radius=radius,
            index=index,
        )

        self.layers_list.append(layer)

    def get_overall_generator(self):
        """
        Return a generator of all combination of fibers.

        :returns:   Generator
        :rtype:     object
        """
        list_of_generator = []

        for layer in self.layers_list:
            generator = layer.get_generator()

            list_of_generator.append(generator)

        overall_generator = product(*list_of_generator)

        return overall_generator

    def __getitem__(self, index: int) -> Fiber:
        """
        Of all the fiber combination, returns the one associated to the given index.

        :returns:   Return the associated fiber
        :rtype:     Fiber
        """
        generator = self.get_overall_generator()

        structure = list(generator)[index]

        fiber = Fiber(wavelength=self.wavelength)

        for name, radius, index in structure:
            fiber.add_layer(
                name=name,
                radius=radius,
                index=index,
            )

        fiber.initialize_layers()

        return fiber

    def __iter__(self) -> Fiber:
        """
        Iterate through all combination of fibers.

        :returns:   Yield the next fiber
        :rtype:     Fiber
        """
        generator = self.get_overall_generator()
        for structure in generator:
            fiber = Fiber(wavelength=self.wavelength)

            for name, radius, index in structure:
                fiber.add_layer(
                    name=name,
                    radius=radius,
                    index=index,
                )

            fiber.initialize_layers()

            yield fiber
# -
