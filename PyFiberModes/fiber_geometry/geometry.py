from PyFiberModes import material
from PyFiberModes.material.fixed import Fixed


class Geometry(object):
    def __init__(self, radius_in: float, radius_out: float, *fp, material_type, material_parameter, **kwargs):
        if material_type == 'Fixed':
            self.material_type = Fixed()

        self.material_parameter = material_parameter
        cm = kwargs.get("cm", None)

        if cm:
            self._cm = material.__dict__[cm]()
            self._cmp = kwargs.get("cmp")

        self.radius_in = radius_in
        self.radius_out = radius_out

    def __str__(self):
        return self.__class__.__name__ + ' ' + self.material_type.str(*self.material_parameter)
