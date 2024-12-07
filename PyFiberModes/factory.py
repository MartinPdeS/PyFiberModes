#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import product
from dataclasses import dataclass, field
import numpy as np
from PyFiberModes.fiber import Fiber


@dataclass
class ProxyLayer:
    """
    Represents a layer configuration in a fiber, with name, radius, and refractive index.

    Attributes
    ----------
    name : str
        Name of the layer.
    radius : list[float]
        Radius values for the layer (in meters).
    index : list[float]
        Refractive index values for the layer.
    """
    name: str
    radius: list = field(default_factory=list)
    index: list = field(default_factory=list)

    def __post_init__(self):
        self.name = [self.name]
        self.radius = np.atleast_1d(self.radius)
        self.index = np.atleast_1d(self.index)

    def get_generator(self):
        """
        Create a generator for all combinations of name, radius, and index.

        Returns
        -------
        generator : itertools.product
            Generator yielding tuples of (name, radius, index).
        """
        return product(self.name, self.radius, self.index)


class FiberFactory:
    """
    Factory to create and manage multiple Fiber instances with various configurations.

    Parameters
    ----------
    wavelength : float
        Wavelength (in meters) used for fiber simulations.

    Attributes
    ----------
    layers_list : list[ProxyLayer]
        List of ProxyLayer objects representing the layers of the fiber.
    neff_solver : None
        Placeholder for a solver object (if required).
    cutoff_solver : None
        Placeholder for a cutoff solver object (if required).
    wavelength : float
        Wavelength used in simulations.
    """

    def __init__(self, wavelength: float):
        self.layers_list = []
        self.neff_solver = None
        self.cutoff_solver = None
        self.wavelength = wavelength

    def add_layer(self, index: float, name: str = "", radius: float = 0.0):
        """
        Add a new layer to the fiber factory.

        Parameters
        ----------
        index : float
            Refractive index of the layer.
        name : str, optional
            Name of the layer (default is an empty string).
        radius : float, optional
            Radius of the layer (default is 0.0 meters).
        """
        layer = ProxyLayer(name=name, radius=radius, index=index)
        self.layers_list.append(layer)

    def get_overall_generator(self):
        """
        Generate all possible combinations of layer configurations.

        Returns
        -------
        overall_generator : itertools.product
            Generator yielding all possible layer configurations.
        """
        generators = [layer.get_generator() for layer in self.layers_list]
        return product(*generators)

    def __getitem__(self, index: int) -> Fiber:
        """
        Retrieve a specific fiber configuration by index.

        Parameters
        ----------
        index : int
            Index of the desired fiber configuration.

        Returns
        -------
        fiber : Fiber
            A Fiber object with the selected configuration.
        """
        structure = list(self.get_overall_generator())[index]
        fiber = self._create_fiber_from_structure(structure)
        return fiber

    def __iter__(self):
        """
        Iterate through all possible fiber configurations.

        Yields
        ------
        fiber : Fiber
            A Fiber object with the next configuration in the sequence.
        """
        for structure in self.get_overall_generator():
            yield self._create_fiber_from_structure(structure)

    def _create_fiber_from_structure(self, structure) -> Fiber:
        """
        Create a Fiber object from a given structure.

        Parameters
        ----------
        structure : iterable
            An iterable containing tuples of (name, radius, index) for each layer.

        Returns
        -------
        fiber : Fiber
            A Fiber object initialized with the given structure.
        """
        fiber = Fiber(wavelength=self.wavelength)
        for name, radius, index in structure:
            fiber.add_layer(name=name, radius=radius, index=index)
        fiber.initialize_layers()
        return fiber
