#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFiberModes.mode import Mode, Family

#: Predefined HE modes
HE11 = Mode(Family.HE, 1, 1)
HE12 = Mode(Family.HE, 1, 2)
HE22 = Mode(Family.HE, 2, 2)

#: Predefined LP modes
LP01 = Mode(Family.LP, 0, 1)
LP11 = Mode(Family.LP, 1, 1)
LP02 = Mode(Family.LP, 0, 2)