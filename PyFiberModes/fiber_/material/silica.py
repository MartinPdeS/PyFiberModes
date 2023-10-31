from .sellmeier import Sellmeier
from PyFiberModes.wavelength import Wavelength


class Silica(Sellmeier):

    """Silica material, based on Sellmeier forumla."""

    name = "Fused Silica"
    info = "Fused silica, at room temperature (20 °C)."
    url = "http://refractiveindex.info/legacy/?group=GLASSES&material=F_SILICA"
    WLRANGE = (Wavelength(0.21e-6), Wavelength(3.71e-6))

    B = (0.6961663, 0.4079426, 0.8974794)
    C = (0.0684043, 0.1162414, 9.896161)
