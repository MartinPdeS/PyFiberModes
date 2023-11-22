#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum
from dataclasses import dataclass

Family = Enum('Family', 'LP HE EH TE TM', module=__name__)


@dataclass(frozen=True, eq=True)
class Mode():
    family: str
    """ Family of the mode """
    nu: int
    """ Parameter of the mode. It often corresponds to the parameter of the radial Bessel functions. """
    m: int
    """ Radial order of the mode (positive integer). It corresponds to the number of concentric rings in the mode fields. """

    def __post_init__(self):
        assert self.family in ['LP', 'HE', 'EH', 'TE', 'TM'], f'Unexpected mode family: {self.family}'

        assert self.nu >= 0, 'Unexpected negative nu value'

        assert self.m >= 0, 'Unexpected negative m value'

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

    def __repr__(self) -> str:
        return f"{self.family}{self.nu}{self.m}"


    def get_lower_neff_mode(self):
        if self.family == 'LP':
            if self.nu == 0:
                return Mode(family='LP', nu=1, m=self.m)
            else:
                return Mode(family='LP', nu=self.nu + 1, m=self.m)


# -
