#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFiberModes.mode import Mode
from PyFiberModes.wavelength import Wavelength
from PyFiberModes.mode_instances import HE11, LP01, LP11, TE01
from PyFiberModes.fundamentals import get_wavelength_from_V0
from PyFiberModes.solver.base_solver import BaseSolver

import numpy
from scipy.special import j0, y0, i0, k0
from scipy.special import j1, y1, i1, k1
from scipy.special import jn, yn, iv, kn
from scipy.special import jvp, ivp

"""
Solver for three layer step-index solver: TLSIF
"""


class NameSpace():
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


class CutoffSolver(BaseSolver):
    def get_lower_neff_mode(self, mode: Mode) -> Mode:
        lower_neff_mode = None

        if mode.family == 'HE':
            lower_neff_mode = Mode('EH', mode.nu, mode.m - 1)
        else:
            lower_neff_mode = Mode(mode.family, mode.nu, mode.m - 1)

        if lower_neff_mode == HE11:
            lower_neff_mode = TE01

        elif lower_neff_mode == LP01:
            lower_neff_mode = LP11

        elif mode.family == 'EH':
            lower_neff_mode = Mode('HE', mode.nu, mode.m)

        elif mode.nu >= 1:  # TE(0,1) is single-mode condition. Roots below TE(0,1) are false-positive
            lower_neff_mode = TE01

        return lower_neff_mode

    def solve(self, mode: Mode):
        lower_neff_mode = self.get_lower_neff_mode(mode=mode)

        if (mode.m >= 2 or mode.family == 'EH'):
            v0_lowbound = self.fiber.get_mode_cutoff_v0(mode=lower_neff_mode)
            delta = 0.05 / v0_lowbound if v0_lowbound > 4 else self._MCD
            v0_lowbound += delta / 100

        elif mode.nu >= 1:
            v0_lowbound = self.fiber.get_mode_cutoff_v0(mode=lower_neff_mode)
            delta = 0.05 / v0_lowbound
            v0_lowbound -= delta / 100
        else:
            v0_lowbound = delta = self._MCD

        if numpy.isnan(delta):
            print(v0_lowbound)

        match mode.family:
            case 'LP':
                function = self._lpcoeq
            case 'TE':
                function = self._tecoeq
            case 'TM':
                function = self._tmcoeq
            case 'HE':
                function = self._hecoeq
            case 'EH':
                function = self._ehcoeq

        return self.find_function_first_root(
            function=function,
            function_args=(mode.nu,),
            lowbound=v0_lowbound,
            delta=delta,
            maxiter=int(250 / delta)
        )

    def get_parameters(self, V0: float, nu: int) -> tuple:
        """
        { function_description }

        :param      V0:   The V0 parameter
        :type       V0:   float

        :returns:   The computed parameters
        :rtype:     tuple
        """
        r1 = self.fiber.layers[0].radius_out
        r2 = self.fiber.layers[1].radius_out

        wavelength = get_wavelength_from_V0(fiber=self.fiber, V0=V0)

        if numpy.isinf(wavelength):
            wavelength = Wavelength(k0=1)  # because it causes troubles if 0

        layers_minimum_index_squared = [
            layer.refractive_index**2 for layer in self.fiber.layers
        ]

        n1sq, n2sq, n3sq = layers_minimum_index_squared

        if wavelength == 0:  # Avoid floating point error. But there should be a way to do it better.
            Usq = [numpy.inf, numpy.inf, numpy.inf]
        else:
            Usq = [wavelength.k0**2 * (nsq - n3sq) for nsq in layers_minimum_index_squared]

        s1, s2, s3 = numpy.sign(Usq)
        u1, u2, u3 = numpy.sqrt(numpy.abs(Usq))

        data_structure = NameSpace(
            nu=nu,
            u1r1=u1 * r1,
            u2r1=u2 * r1,
            u2r2=u2 * r2,
            s1=s1,
            s2=s2,
            n1sq=n1sq,
            n2sq=n2sq,
            n3sq=n3sq,
            s3=s3
        )

        return data_structure

    def _get_delta_(self, p: NameSpace) -> float:
        """
        Gets the delta. s3 is sign of Delta

        :param      p:    { parameter_description }
        :type       p:    NameSpace

        :returns:   The delta.
        :rtype:     float
        """
        if p.s1 < 0:
            f = ivp(p.nu, p.u1r1) / (iv(p.nu, p.u1r1) * p.u1r1)  # c
        else:
            jnnuu1r1 = jn(p.nu, p.u1r1)
            if jnnuu1r1 == 0:  # Avoid zero division error
                return numpy.inf
            f = jvp(p.nu, p.u1r1) / (jnnuu1r1 * p.u1r1)  # a b d

        if p.s1 == p.s2:
            # b d
            kappa1 = -(p.n1sq + p.n2sq) * f / p.n2sq
            kappa2 = p.n1sq * f * f / p.n2sq - p.nu**2 * p.n3sq / p.n2sq * (1 / p.u1r1**2 - 1 / p.u2r1**2)**2
        else:
            # a c
            kappa1 = (p.n1sq + p.n2sq) * f / p.n2sq
            kappa2 = p.n1sq * f * f / p.n2sq - p.nu**2 * p.n3sq / p.n2sq * (1 / p.u1r1**2 + 1 / p.u2r1**2)**2

        d = kappa1**2 - 4 * kappa2
        if d < 0:
            return numpy.nan
        return p.u2r1 * (p.nu / p.u2r1**2 + (kappa1 + p.s3 * numpy.sqrt(d)) * 0.5)

    def _lpcoeq(self, v0: float, nu: int) -> float:
        """
        Not sure here.

        :param      v0:   The V parameter
        :type       v0:   float
        :param      nu:   The radial parameter of the mode.
        :type       nu:   int

        :returns:   Not to sure
        :rtype:     float
        """
        p = self.get_parameters(v0, nu)

        if p.s1 == 0:  # e
            return jn(nu + 1, p.u2r1) * yn(nu - 1, p.u2r2) - yn(nu + 1, p.u2r1) * jn(nu - 1, p.u2r2)

        if p.s1 > 0:
            f11a, f11b = jn(nu - 1, p.u1r1), jn(nu, p.u1r1)
        else:
            f11a, f11b = iv(nu - 1, p.u1r1), iv(nu, p.u1r1)

        if p.s2 > 0:
            f22a, f22b = jn(nu - 1, p.u2r2), yn(nu - 1, p.u2r2)
            f2a = jn(nu, p.u2r1) * f22b - yn(nu, p.u2r1) * f22a
            f2b = jn(nu - 1, p.u2r1) * f22b - yn(nu - 1, p.u2r1) * f22a
        else:  # a
            f22a, f22b = iv(nu - 1, p.u2r2), kn(nu - 1, p.u2r2)
            f2a = iv(nu, p.u2r1) * f22b + kn(nu, p.u2r1) * f22a
            f2b = iv(nu - 1, p.u2r1) * f22b - kn(nu - 1, p.u2r1) * f22a

        return f11a * f2a * p.u1r1 - f11b * f2b * p.u2r1

    def _tecoeq(self, v0: float, nu: int) -> float:
        """
        Not sure here.

        :param      v0:   The V parameter
        :type       v0:   float
        :param      nu:   The radial parameter of the mode.
        :type       nu:   int

        :returns:   Not to sure
        :rtype:     float
        """
        p = self.get_parameters(v0, nu)

        if p.s1 > 0:
            f11a, f11b = j0(p.u1r1), jn(2, p.u1r1)
        else:
            f11a, f11b = i0(p.u1r1), -iv(2, p.u1r1)
        if p.s2 > 0:
            f22a, f22b = j0(p.u2r2), y0(p.u2r2)
            f2a = jn(2, p.u2r1) * f22b - yn(2, p.u2r1) * f22a
            f2b = j0(p.u2r1) * f22b - y0(p.u2r1) * f22a
        else:  # a
            f22a, f22b = i0(p.u2r2), k0(p.u2r2)
            f2a = kn(2, p.u2r1) * f22a - iv(2, p.u2r1) * f22b
            f2b = i0(p.u2r1) * f22b - k0(p.u2r1) * f22a
        return f11a * f2a - f11b * f2b

    def _tmcoeq(self, v0: float, nu: int) -> float:
        """
        Not sure here.

        :param      v0:   The V parameter
        :type       v0:   float
        :param      nu:   The radial parameter of the mode.
        :type       nu:   int

        :returns:   Not to sure
        :rtype:     float
        """
        p = self.get_parameters(v0, nu)

        if p.s1 == 0:  # e
            f11a, f11b = 2, 1
        elif p.s1 > 0:  # a, b, d
            f11a, f11b = j0(p.u1r1) * p.u1r1, j1(p.u1r1)
        else:  # c
            f11a, f11b = i0(p.u1r1) * p.u1r1, i1(p.u1r1)
        if p.s2 > 0:
            f22a, f22b = j0(p.u2r2), y0(p.u2r2)
            f2a = j1(p.u2r1) * f22b - y1(p.u2r1) * f22a
            f2b = j0(p.u2r1) * f22b - y0(p.u2r1) * f22a
        else:  # a
            f22a, f22b = i0(p.u2r2), k0(p.u2r2)
            f2a = i1(p.u2r1) * f22b + k1(p.u2r1) * f22a
            f2b = i0(p.u2r1) * f22b - k0(p.u2r1) * f22a
        return f11a * p.n2sq * f2a - f11b * p.n1sq * f2b * p.u2r1

    def _ehcoeq(self, v0: float, nu: int) -> float:
        """
        Not sure here.

        :param      v0:   The V parameter
        :type       v0:   float
        :param      nu:   The radial parameter of the mode.
        :type       nu:   int

        :returns:   Not to sure
        :rtype:     float
        """
        p = self.get_parameters(v0, nu)

        if p.s1 == 0:
            value = self._function_3_(
                nu=nu,
                u2r1=p.u2r1,
                u2r2=p.u2r2,
                dn=2,
                n2sq=p.n2sq,
                n3sq=p.n3sq
            )
        else:
            p.s3 = 1 if p.s1 == p.s2 else -1

            value = self._function_1_(p=p)

        return value

    def _hecoeq(self, v0: float, nu: int):
        """
        Not sure here.

        :param      v0:   The V parameter
        :type       v0:   float
        :param      nu:   The radial parameter of the mode.
        :type       nu:   int

        :returns:   Not to sure
        :rtype:     float
        """
        p = self.get_parameters(v0, nu)

        if p.s1 == 0:
            p.dn = -2
            return self._function_3_(p=p)

        else:
            p.s3 = -1 if p.s1 == p.s2 else 1
            if p.n1sq > p.n2sq > p.n3sq:
                p.s3 = -1 if p.nu == 1 else 1

            return self._function_2_(p=p)

    def _function_1_(self, p: NameSpace) -> float:
        if p.s2 < 0:  # a
            b11 = iv(p.nu, p.u2r1)
            b12 = kn(p.nu, p.u2r1)
            b21 = iv(p.nu, p.u2r2)
            b22 = kn(p.nu, p.u2r2)
            b31 = iv(p.nu + 1, p.u2r1)
            b32 = kn(p.nu + 1, p.u2r1)
            f1 = b31 * b22 + b32 * b21
            f2 = b11 * b22 - b12 * b21
        else:
            b11 = jn(p.nu, p.u2r1)
            b12 = yn(p.nu, p.u2r1)
            b21 = jn(p.nu, p.u2r2)
            b22 = yn(p.nu, p.u2r2)
            if p.s1 == 0:
                f1 = 0
            else:
                b31 = jn(p.nu + 1, p.u2r1)
                b32 = yn(p.nu + 1, p.u2r1)
                f1 = b31 * b22 - b32 * b21
            f2 = b12 * b21 - b11 * b22
        if p.s1 == 0:
            delta = 1
        else:
            delta = self._get_delta_(p=p)
        return f1 + f2 * delta

    def _function_2_(self, p: NameSpace) -> float:

        with numpy.errstate(invalid='ignore'):
            delta = self._get_delta_(p=p)

            n0sq = (p.n3sq - p.n2sq) / (p.n2sq + p.n3sq)

            if p.s2 < 0:  # a
                b11 = iv(p.nu, p.u2r1)
                b12 = kn(p.nu, p.u2r1)
                b21 = iv(p.nu, p.u2r2)
                b22 = kn(p.nu, p.u2r2)
                b31 = iv(p.nu + 1, p.u2r1)
                b32 = kn(p.nu + 1, p.u2r1)
                b41 = iv(p.nu - 2, p.u2r2)
                b42 = kn(p.nu - 2, p.u2r2)
                g1 = b11 * delta + b31
                g2 = b12 * delta - b32
                f1 = b41 * g2 - b42 * g1
                f2 = b21 * g2 - b22 * g1
            else:
                b11 = jn(p.nu, p.u2r1)
                b12 = yn(p.nu, p.u2r1)
                b21 = jn(p.nu, p.u2r2)
                b22 = yn(p.nu, p.u2r2)
                b31 = jn(p.nu + 1, p.u2r1)
                b32 = yn(p.nu + 1, p.u2r1)
                b41 = jn(p.nu - 2, p.u2r2)
                b42 = yn(p.nu - 2, p.u2r2)
                g1 = b11 * delta - b31
                g2 = b12 * delta - b32
                f1 = b41 * g2 - b42 * g1
                f2 = b22 * g1 - b21 * g2
            return f1 + n0sq * f2

    def _function_3_(self, p: NameSpace) -> float:
        n0sq = (p.n3sq - p.n2sq) / (p.n2sq + p.n3sq)
        b11 = jn(p.nu, p.u2r1)
        b12 = yn(p.nu, p.u2r1)
        b21 = jn(p.nu, p.u2r2)
        b22 = yn(p.nu, p.u2r2)

        if p.dn > 0:
            b31 = jn(p.nu + p.dn, p.u2r1)
            b32 = yn(p.nu + p.dn, p.u2r1)
            f1 = b31 * b22 - b32 * b21
            f2 = b11 * b22 - b12 * b21
        else:
            b31 = jn(p.nu + p.dn, p.u2r2)
            b32 = yn(p.nu + p.dn, p.u2r2)
            f1 = b31 * b12 - b32 * b11
            f2 = b12 * b21 - b11 * b22

        return f1 - n0sq * f2
