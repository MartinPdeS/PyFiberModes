from .sellmeier import Sellmeier
from PyFiberModes.wavelength import Wavelength


class Germania(Sellmeier):

    """Germania material, based on Sellmeier forumla."""

    name = "Germanium dioxide"
    info = "Fused germanium dioxide."
    url = "http://refractiveindex.info/?shelf=main&book=GeO2&page=Fleming"
    WLRANGE = (Wavelength(0.36e-6), Wavelength(4.3e-6))
    B = (0.80686642, 0.71815848, 0.85416831)
    C = (0.068972606, 0.15396605, 11.841931)
