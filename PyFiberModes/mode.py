"""Mode-family definitions and immutable mode identifiers."""

from enum import Enum
from dataclasses import dataclass


class Family(str, Enum):
    """Supported scalar and vector mode families."""

    LP = "LP"
    HE = "HE"
    EH = "EH"
    TE = "TE"
    TM = "TM"


@dataclass(frozen=True, eq=True)
class Mode():
    """Identify one circular-fiber propagation mode.

    Parameters
    ----------
    family : {'LP', 'HE', 'EH', 'TE', 'TM'}
        Electromagnetic mode family.
    nu : int
        Non-negative azimuthal order.
    m : int
        Non-negative radial order.
    """
    family: str | Family
    """ Family of the mode """
    nu: int
    """ Parameter of the mode. It often corresponds to the parameter of the radial Bessel functions. """
    m: int
    """ Radial order of the mode (positive integer). It corresponds to the number of concentric rings in the mode fields. """

    def __post_init__(self):
        """Validate the mode family and modal orders.

        Raises
        ------
        ValidationError
            If the family is unknown, ``nu`` is negative, or ``m`` is not a
            positive integer.
        """
        from PyFiberModes.exceptions import ValidationError

        try:
            family = self.family if isinstance(self.family, Family) else Family(self.family)
        except (TypeError, ValueError) as error:
            choices = ", ".join(member.value for member in Family)
            raise ValidationError(
                f"unknown mode family {self.family!r}; expected one of {choices}"
            ) from error
        if not isinstance(self.nu, int) or isinstance(self.nu, bool) or self.nu < 0:
            raise ValidationError("nu must be a non-negative integer")
        if not isinstance(self.m, int) or isinstance(self.m, bool) or self.m < 1:
            raise ValidationError("m must be a positive integer")

        object.__setattr__(self, "family", family.value)

    def __repr__(self) -> str:
        """Return the compact family-and-order representation.

        Returns
        -------
        str
            Mode label such as ``"LP01"``.
        """
        return f"{self.family}{self.nu}{self.m}"

# -
