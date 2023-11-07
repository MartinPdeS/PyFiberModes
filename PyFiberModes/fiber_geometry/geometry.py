from PyFiberModes.material.fixed import Fixed
from dataclasses import dataclass, field


@dataclass
class Geometry(object):
    radius_in: float
    radius_out: float
    material_parameter: str
    material_type: object = field(default_factory='Fixed')

    def __post_init__(self):
        if self.material_type.lower() == 'fixed':
            self.material_type = Fixed()

# -
