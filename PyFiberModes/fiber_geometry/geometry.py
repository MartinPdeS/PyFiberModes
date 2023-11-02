from PyFiberModes import material


class Geometry(object):
    def __init__(self, radius_in: float, radius_out: float, *fp, m, mp, **kwargs):
        self._m = material.__dict__[m]()  # instantiate material object
        self._mp = mp
        cm = kwargs.get("cm", None)
        # dsa
        if cm:
            self._cm = material.__dict__[cm]()
            self._cmp = kwargs.get("cmp")

        self._fp = fp
        self.radius_in = radius_in
        self.radius_out = radius_out

    def __str__(self):
        return self.__class__.__name__ + ' ' + self._m.str(*self._mp)
