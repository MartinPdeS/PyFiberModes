#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.constants import c, pi


class Wavelength(float):
    """
    Class for wavelength unit conversions.

    Inherits from :py:class:`float` to allow seamless use wherever floats are used.
    The wavelength is always stored in meters.

    Attributes can be accessed to convert wavelength into wave number, frequency, or angular frequency.
    """

    def __new__(cls, *args, **kwargs):
        """
        Create a new Wavelength object.

        Parameters
        ----------
        args : tuple
            Positional arguments to specify the wavelength (default in meters).
        kwargs : dict
            Keyword arguments for alternative unit-based construction:
                - 'k0': Wave number (2π/λ).
                - 'omega', 'w': Angular frequency (rad/s).
                - 'wl', 'wavelength': Wavelength (in meters).
                - 'frequency', 'v', 'f': Frequency (in Hertz).

        Returns
        -------
        Wavelength
            A new Wavelength instance in meters.

        Raises
        ------
        TypeError
            If multiple arguments are provided or if arguments are invalid.
        """
        if len(args) + len(kwargs) > 1:
            raise TypeError("Wavelength constructor requires exactly one parameter.")

        if args:
            wl = args[0]
        elif 'k0' in kwargs:
            wl = 2 * pi / kwargs['k0']
        elif 'omega' in kwargs or 'w' in kwargs:
            omega = kwargs.get('omega') or kwargs['w']
            wl = c * 2 * pi / omega
        elif 'wl' in kwargs or 'wavelength' in kwargs:
            wl = kwargs.get('wl') or kwargs['wavelength']
        elif 'frequency' in kwargs or 'v' in kwargs or 'f' in kwargs:
            frequency = kwargs.get('frequency') or kwargs.get('v') or kwargs['f']
            wl = c / frequency
        else:
            raise TypeError("Invalid argument for Wavelength constructor.")

        return super().__new__(cls, wl)

    @property
    def k0(self) -> float:
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

    w = omega

    @property
    def wavelength(self) -> float:
        """
        Wavelength (in meters).

        Returns
        -------
        float
            Wavelength in meters.
        """
        return self

    wl = wavelength

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

    v = frequency
    f = frequency

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
