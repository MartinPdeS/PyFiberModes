"""Mode-family definitions and immutable mode identifiers."""

from enum import Enum
from dataclasses import dataclass

Family = Enum('Family', 'LP HE EH TE TM', module=__name__)


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
    family: str
    """ Family of the mode """
    nu: int
    """ Parameter of the mode. It often corresponds to the parameter of the radial Bessel functions. """
    m: int
    """ Radial order of the mode (positive integer). It corresponds to the number of concentric rings in the mode fields. """

    def __post_init__(self):
        """Validate the mode family and modal orders.

        Raises
        ------
        AssertionError
            If the family is unknown or either modal order is negative.
        """
        assert self.family in ['LP', 'HE', 'EH', 'TE', 'TM'], f'Unexpected mode family: {self.family}'

        assert self.nu >= 0, 'Unexpected negative nu value'

        assert self.m >= 0, 'Unexpected negative m value'

    def get_LP_equvalent_mode(self):  # previously lpEq
        """Return the weak-guidance LP equivalent of this mode.

        Returns
        -------
        Mode
            Corresponding linearly polarized mode.
        """
        if self.family is Family.LP:
            return self
        elif self.family is Family.HE:
            return Mode(Family.LP, self.nu - 1, self.m)
        else:
            return Mode(Family.LP, self.nu + 1, self.m)

    def __repr__(self) -> str:
        """Return the compact family-and-order representation.

        Returns
        -------
        str
            Mode label such as ``"LP01"``.
        """
        return f"{self.family}{self.nu}{self.m}"

    def get_lower_neff_mode(self):
        """Return the adjacent LP mode with lower effective index.

        Returns
        -------
        Mode or None
            Next LP azimuthal order, or ``None`` for non-LP families.
        """
        if self.family == 'LP':
            if self.nu == 0:
                return Mode(family='LP', nu=1, m=self.m)
            else:
                return Mode(family='LP', nu=self.nu + 1, m=self.m)


# -
