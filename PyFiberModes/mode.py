#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum
from collections import namedtuple


#: Constants for identifying mode family. Can be
#: LP, HE, EH, TE, or TM.
Family = Enum('Family', 'LP HE EH TE TM', module=__name__)


class Mode(namedtuple('Mode', 'family nu m')):
    """
    Fiber mode representation.

    The fiber mode consists of a mode family, and two mode parameters
    (*ν* and *m*). If is derived from namedtuple. Therefore, it is
    unmutable, and it can be used as dictionary key.

    .. py:attribute:: family

        Mode family.

        .. seealso:: :py:class:`Family`

    .. py:attribute:: nu

        *ν* parameter of the mode. It often corresponds to the parameter
        of the radial Bessel functions.

    .. py:attribute:: m

        (positive integer) Radial order of the mode.
        It corresponds to the number of concentric rings in the mode fields.

    """

    def __new__(cls, family, nu, m):
        if not isinstance(family, Family):
            family = Family[family]

        return super(Mode, cls).__new__(cls, family, nu, m)

    # def get_lower_neff_mode(self):
    #     if self.family is ModeFamily.HE:
    #         lower_neff_mode = Mode(ModeFamily.EH, mode.nu, mode.m - 1)

    def get_LP_equvalent_mode(self):  # previously lpEq
        """
        Gets the equivalent LP mode.
        """
        if self.family is Family.LP:
            return self
        elif self.family is Family.HE:
            return Mode(Family.LP, self.nu - 1, self.m)
        else:
            return Mode(Family.LP, self.nu + 1, self.m)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"{self.family.name}{self.nu}{self.m}"

# -
