from .material import Material


class Fixed(Material):

    """Fixed index material class.

    A material with a fixed index always have the same refractive index,
    whatever the wavelength is.

    """

    name = "Fixed index"
    info = "Fixed index"
    nparams = 1

    @classmethod
    def get_refractive_index(cls, wavelength: float, index: float):
        return index
