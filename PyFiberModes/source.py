#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.constants import c, pi


class Source():
    def __init__(self, wavelength: float):
        self.wavelength = wavelength

    @property
    def wavenumber(self) -> float:
        """
        Wave number (2π/λ).

        Returns
        -------
        float
            Wave number in radians per meter.
        """
        return 2 * pi / self if self != 0 else float('inf')

    @property
    def omega(self) -> float:
        """
        Angular frequency (in rad/s).

        Returns
        -------
        float
            Angular frequency.
        """
        return c * 2 * pi / self if self != 0 else float('inf')

    @property
    def frequency(self) -> float:
        """
        Frequency (in Hertz).

        Returns
        -------
        float
            Frequency in Hz.
        """
        return c / self if self != 0 else float('inf')

    def __str__(self) -> str:
        """
        String representation of the wavelength in nanometers.

        Returns
        -------
        str
            Wavelength formatted as a string in nanometers with 2 decimal places.
        """
        return f"{1e9 * self.wavelength:.2f} nm"

    def __repr__(self) -> str:
        """
        Debug representation of the Wavelength object.

        Returns
        -------
        str
            Debug string representation.
        """
        return f"Wavelength({self.wavelength:.6e} meters)"
