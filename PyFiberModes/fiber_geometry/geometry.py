from PyFiberModes.material.fixed import Fixed
from dataclasses import dataclass, field


@dataclass
class Geometry(object):
    radius_in: float
    """ Minimum radius of the structure """
    radius_out: float
    """ Maximum radius of the structure """
    index_list: str
    """ Refractive index of the structure """
    material_type: object = field(default_factory=None)
    """ Material type """

    def __post_init__(self):
        if self.material_type.lower() == 'fixed':
            self.material_type = Fixed()

# -
