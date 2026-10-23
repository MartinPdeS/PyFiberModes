"""Factories for constructing parameterized multilayer fibers."""

from itertools import product
from dataclasses import dataclass, field
import numpy as np
from PyFiberModes.fiber import Fiber
from PyFiberModes.units import Length, RefractiveIndex, magnitude, quantity, ureg


@dataclass
class ProxyLayer:
    """
    Represents a layer configuration in a fiber, with name, radius, and refractive index.

    Parameters
    ----------
    name : str
        Name of the layer.
    radius : pint.Quantity
        Radius values for the layer (in meters).
    index : pint.Quantity
        Refractive index values for the layer.
    """
    name: str
    radius: list = field(default_factory=list)
    index: list = field(default_factory=list)

    def __post_init__(self):
        """Normalize scalar layer parameters into one-dimensional arrays."""
        self.name = [self.name]
        self.radius = quantity(np.atleast_1d(magnitude(self.radius, "meter", "radius")), "meter")
        self.index = quantity(np.atleast_1d(magnitude(self.index, "dimensionless", "index")), "dimensionless")

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
    wavelength : TypedUnit.Length
        Wavelength (in meters) used for fiber simulations.

    Attributes
    ----------
    layers_list : list[ProxyLayer]
        List of ProxyLayer objects representing the layers of the fiber.
    wavelength : TypedUnit.Length
        Wavelength used in simulations.
    """

    def __init__(self, wavelength: Length):
        """Initialize an empty fiber-configuration factory.

        Parameters
        ----------
        wavelength : TypedUnit.Length
            Vacuum wavelength in meters for generated fibers.
        """
        self.layers_list = []
        self.wavelength = quantity(magnitude(wavelength, "meter", "wavelength"), "meter")

    def add_layer(self, index: RefractiveIndex, name: str = "", radius: Length = np.inf * ureg.meter):
        """
        Add a new layer to the fiber factory.

        Parameters
        ----------
        index : TypedUnit.RefractiveIndex
            Refractive index of the layer.
        name : str, optional
            Name of the layer (default is an empty string).
        radius : TypedUnit.Length, optional
            Radius of the layer (defaults to an unbounded outer layer).
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
