"""Optical-source wavelength and frequency conversions."""

import math

from scipy.constants import c, pi
from PyFiberModes.units import Length, scalar, quantity


class Source():
    """Represent a monochromatic optical source.

    Parameters
    ----------
    wavelength : TypedUnit.Length
        Vacuum wavelength in meters.
    """

    def __init__(self, wavelength: Length):
        """Initialize a source from its vacuum wavelength.

        Parameters
        ----------
        wavelength : TypedUnit.Length
            Vacuum wavelength in meters.
        """
        value = scalar(wavelength, "meter", "wavelength")
        if value <= 0 or not math.isfinite(value):
            raise ValueError("wavelength must be positive and finite")
        self._wavelength = quantity(value, "meter")

    @property
    def wavelength(self) -> Length:
        """Validated vacuum wavelength."""
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value: Length):
        """Validate wavelength mutations as well as construction."""
        wavelength = scalar(value, "meter", "wavelength")
        if not math.isfinite(wavelength) or wavelength <= 0:
            raise ValueError("wavelength must be positive and finite")
        self._wavelength = quantity(wavelength, "meter")

    @property
    def wavenumber(self):
        """
        Wave number (2π/λ).

        Returns
        -------
        pint.Quantity
            Wave number in radians per meter.
        """
        return 2 * pi * quantity(1, "radian") / self.wavelength

    @property
    def omega(self):
        """
        Angular frequency (in rad/s).

        Returns
        -------
        pint.Quantity
            Angular frequency.
        """
        return quantity(c, "meter / second") * 2 * pi * quantity(1, "radian") / self.wavelength

    @property
    def frequency(self):
        """
        Frequency (in Hertz).

        Returns
        -------
        pint.Quantity
            Frequency in Hz.
        """
        return quantity(c, "meter / second") / self.wavelength

    def __str__(self) -> str:
        """
        String representation of the wavelength in nanometers.

        Returns
        -------
        str
            Wavelength formatted as a string in nanometers with 2 decimal places.
        """
        return f"{self.wavelength.to('nanometer').magnitude:.2f} nm"

    def __repr__(self) -> str:
        """
        Debug representation of the Wavelength object.

        Returns
        -------
        str
            Debug string representation.
        """
        return f"Source(wavelength={self.wavelength!r})"
