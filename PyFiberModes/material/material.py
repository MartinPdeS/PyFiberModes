from scipy.optimize import brentq
import warnings


class OutOfRangeWarning(UserWarning):

    """Warning raised when a material parameter is out of the allowed range.

    If this warning is raised, it means that results are possibily
    innacurate.

    """
    pass


class Material(object):

    """Material abstract class.

    This gives the interface for the different materials, as well as
    some common functions. All methods in that class are class methods.

    """

    name = "Abstract Material"  # English name for the mateial
    nparams = 0  # Number of material parameters
    info = None  # Info about the material
    url = None  # URL for the reference about the material
    WLRANGE = None  # Acceptable range for the wavelength

    @classmethod
    def _testRange(cls, wavelength: float):
        if cls.WLRANGE is None:
            return
        if cls.WLRANGE[0] <= wavelength <= cls.WLRANGE[1]:
            return

        msg = f"""Wavelength {wavelength} out of supported range for material {cls.name}. 
        Wavelength should be in the range {cls.WLRANGE[0]} - {cls.WLRANGE[1]}.
        Results could be innacurate."""

        warnings.warn(msg, OutOfRangeWarning)

    @classmethod
    def n(cls, wavelength: float, *args, **kwargs):
        raise NotImplementedError(
            "This method must be implemented in derived class."
        )

    @classmethod
    def wlFromN(cls, n, *args, **kwargs):
        def f(wl):
            return n - cls.n(wl, *args, **kwargs)

        if cls.WLRANGE is None:
            raise NotImplementedError(
                "This method only works if WLRANGE is defined"
            )

        a = f(cls.WLRANGE[0])
        b = f(cls.WLRANGE[1])
        if a * b > 0:
            warnings.warn(f"Index {n} out of range.", OutOfRangeWarning)
            return None

        return brentq(f, cls.WLRANGE[0], cls.WLRANGE[1])

    @classmethod
    def str(cls, *args):
        return cls.name

    @classmethod
    def __str__(cls):
        return cls.name

    @classmethod
    def __repr__(cls):
        return f"{cls.__name__}()"
