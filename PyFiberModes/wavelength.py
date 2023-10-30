import numpy
import scipy


class Wavelength(float):
    """Easy wavelength units conversion class.

    This class is inherited from :py:class:`float`. Therefore, it can be
    used wherever you can use floats. Wavelength always is expressed in
    meters.

    Properties can be used to convert wavelength to frequency or wavenumber.

    """

    def __new__(cls, *args, **kwargs):
        """Construct a Wavelength object, using given value.

        You can pass to the constructor any keyword defined in properties
        (k0, omega, w, wl, wavelength, frequency, v, or f).
        If no keyword is given, value is considered to be wavelength.

        """
        nargs = len(args) + len(kwargs)
        if nargs > 1:
            raise TypeError("Wavelength constructor need exactly one parameter")
        if nargs == 0:
            wl = 0
        elif len(args) == 1:
            wl = args[0]
        elif 'k0' in kwargs:
            wl = 2 * numpy.pi / kwargs['k0']
        elif 'omega' in kwargs:
            wl = scipy.c * 2 * numpy.pi / kwargs['omega']
        elif 'w' in kwargs:
            wl = scipy.c * 2 * numpy.pi / kwargs['w']
        elif 'wl' in kwargs:
            wl = kwargs['wl']
        elif 'wavelength' in kwargs:
            wl = kwargs['wavelength']
        elif 'frequency' in kwargs:
            wl = scipy.c / kwargs['frequency']
        elif 'v' in kwargs:
            wl = scipy.c / kwargs['v']
        elif 'f' in kwargs:
            wl = scipy.c / kwargs['f']
        else:
            raise TypeError("Invalid argument")

        return float.__new__(cls, wl)

    @property
    def k0(self):
        """Wave number (:math:`2 \pi / \lambda`)."""
        return 2 * numpy.pi / self if self != 0 else float("inf")

    @property
    def omega(self):
        """Angular frequency (in rad/s)."""
        return scipy.c * 2 * numpy.pi / self if self != 0 else float("inf")

    w = omega

    @property
    def wavelength(self):
        """Wavelength (in meters)."""
        return self

    wl = wavelength

    @property
    def frequency(self):
        """Frequency (in Hertz)."""
        return scipy.c / self if self != 0 else float("inf")

    v = frequency
    f = frequency

    def __str__(self):
        """Format wavelength as string (nanometers, 2 digits)"""
        return "{:.2f} nm".format(1e9 * self.wavelength)

    def __repr__(self):
        return "wavelength({})".format(self.wavelength)
