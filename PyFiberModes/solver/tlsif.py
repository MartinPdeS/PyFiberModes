#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFiberModes import solver
from PyFiberModes import Mode, ModeFamily, Wavelength, HE11
from math import sqrt, isinf, isnan
import numpy
from scipy.special import j0, y0, i0, k0
from scipy.special import j1, y1, i1, k1
from scipy.special import jn, yn, iv, kn
from scipy.special import jvp, ivp


"""
Solver for three layer step-index solver: TLSIF
"""


class CutoffSolver(solver.solver.FiberSolver):

    def solve(self, mode: Mode):
        if mode.m > 1:
            if mode.family is ModeFamily.HE:
                pm = Mode(ModeFamily.EH, mode.nu, mode.m - 1)
            else:
                pm = Mode(mode.family, mode.nu, mode.m - 1)
            if pm == HE11:
                pm = Mode(ModeFamily.TE, 0, 1)
            elif pm == Mode(ModeFamily.LP, 0, 1):
                pm = Mode(ModeFamily.LP, 1, 1)
            lowbound = self.fiber.get_cutoff(mode=pm)
            delta = 0.05 / lowbound if lowbound > 4 else self._MCD
            lowbound += delta / 100

        elif mode.family is ModeFamily.EH:
            pm = Mode(ModeFamily.HE, mode.nu, mode.m)
            lowbound = self.fiber.get_cutoff(pm)
            delta = 0.05 / lowbound if lowbound > 4 else self._MCD
            lowbound += delta / 100

        elif mode.nu > 0:
            # TE(0,1) is single-mode condition. Roots below TE(0,1) are false-positive
            pm = Mode(ModeFamily.TE, 0, 1)
            lowbound = self.fiber.cutoff(pm)

            delta = 0.05 / lowbound
            lowbound -= delta / 100
        else:
            lowbound = delta = self._MCD

        if isnan(delta):
            print(lowbound)

        match mode.family:
            case ModeFamily.LP:
                function = self._lpcoeq
            case ModeFamily.TE:
                function = self._tecoeq
            case ModeFamily.TM:
                function = self._tmcoeq
            case ModeFamily.HE:
                function = self._hecoeq
            case ModeFamily.EH:
                function = self._ehcoeq

        return self.find_function_first_root(
            function=function,
            function_args=(mode.nu,),
            lowbound=lowbound,
            delta=delta,
            maxiter=int(250 / delta)
        )

    def get_parameters(self, V0: float) -> tuple:
        """
        { function_description }

        :param      V0:   The V0 parameter
        :type       V0:   float

        :returns:   { description_of_the_return_value }
        :rtype:     tuple
        """
        r1, r2 = self.fiber.layer_radius

        wavelength = self.fiber.V0_to_wavelength(V0=V0)
        if isinf(wavelength):
            wavelength = Wavelength(k0=1)  # because it causes troubles if 0

        layers_minimum_index_squared = [
            layer.get_minimum_index(wavelength=wavelength)**2 for layer in self.fiber
        ]

        n1sq, n2sq, n3sq = layers_minimum_index_squared

        if wavelength == 0:
            # Avoid floating point error. But there should be a
            # way to do it better.
            Usq = [numpy.inf, numpy.inf, numpy.inf]
        else:
            Usq = [wavelength.k0**2 * (nsq - n3sq) for nsq in layers_minimum_index_squared]
        s1, s2, s3 = numpy.sign(Usq)
        u1, u2, u3 = numpy.sqrt(numpy.abs(Usq))

        return u1 * r1, u2 * r1, u2 * r2, s1, s2, n1sq, n2sq, n3sq

    def __delta(self, nu, u1r1, u2r1, s1, s2, s3, n1sq, n2sq, n3sq):
        """s3 is sign of Delta"""
        if s1 < 0:
            f = ivp(nu, u1r1) / (iv(nu, u1r1) * u1r1)  # c
        else:
            jnnuu1r1 = jn(nu, u1r1)
            if jnnuu1r1 == 0:  # Avoid zero division error
                return float("inf")
            f = jvp(nu, u1r1) / (jnnuu1r1 * u1r1)  # a b d

        if s1 == s2:
            # b d
            kappa1 = -(n1sq + n2sq) * f / n2sq
            kappa2 = n1sq * f * f / n2sq - nu**2 * n3sq / n2sq * (1 / u1r1**2 - 1 / u2r1**2)**2
        else:
            # a c
            kappa1 = (n1sq + n2sq) * f / n2sq
            kappa2 = n1sq * f * f / n2sq - nu**2 * n3sq / n2sq * (1 / u1r1**2 + 1 / u2r1**2)**2

        d = kappa1**2 - 4 * kappa2
        if d < 0:
            return numpy.nan
        return u2r1 * (nu / u2r1**2 + (kappa1 + s3 * sqrt(d)) * 0.5)

    def _lpcoeq(self, v0: float, nu: int) -> float:
        u1r1, u2r1, u2r2, s1, s2, n1sq, n2sq, n3sq = self.get_parameters(v0)

        if s1 == 0:  # e
            return jn(nu + 1, u2r1) * yn(nu - 1, u2r2) - yn(nu + 1, u2r1) * jn(nu - 1, u2r2)

        (f11a, f11b) = ((jn(nu - 1, u1r1), jn(nu, u1r1)) if s1 > 0 else
                        (iv(nu - 1, u1r1), iv(nu, u1r1)))
        if s2 > 0:
            f22a, f22b = jn(nu - 1, u2r2), yn(nu - 1, u2r2)
            f2a = jn(nu, u2r1) * f22b - yn(nu, u2r1) * f22a
            f2b = jn(nu - 1, u2r1) * f22b - yn(nu - 1, u2r1) * f22a
        else:  # a
            f22a, f22b = iv(nu - 1, u2r2), kn(nu - 1, u2r2)
            f2a = iv(nu, u2r1) * f22b + kn(nu, u2r1) * f22a
            f2b = iv(nu - 1, u2r1) * f22b - kn(nu - 1, u2r1) * f22a
        return f11a * f2a * u1r1 - f11b * f2b * u2r1

    def _tecoeq(self, v0: float, nu: int) -> float:
        u1r1, u2r1, u2r2, s1, s2, n1sq, n2sq, n3sq = self.get_parameters(v0)
        (f11a, f11b) = ((j0(u1r1), jn(2, u1r1)) if s1 > 0 else
                        (i0(u1r1), -iv(2, u1r1)))
        if s2 > 0:
            f22a, f22b = j0(u2r2), y0(u2r2)
            f2a = jn(2, u2r1) * f22b - yn(2, u2r1) * f22a
            f2b = j0(u2r1) * f22b - y0(u2r1) * f22a
        else:  # a
            f22a, f22b = i0(u2r2), k0(u2r2)
            f2a = kn(2, u2r1) * f22a - iv(2, u2r1) * f22b
            f2b = i0(u2r1) * f22b - k0(u2r1) * f22a
        return f11a * f2a - f11b * f2b

    def _tmcoeq(self, v0: float, nu: int) -> float:
        u1r1, u2r1, u2r2, s1, s2, n1sq, n2sq, n3sq = self.get_parameters(v0)
        if s1 == 0:  # e
            f11a, f11b = 2, 1
        elif s1 > 0:  # a, b, d
            f11a, f11b = j0(u1r1) * u1r1, j1(u1r1)
        else:  # c
            f11a, f11b = i0(u1r1) * u1r1, i1(u1r1)
        if s2 > 0:
            f22a, f22b = j0(u2r2), y0(u2r2)
            f2a = j1(u2r1) * f22b - y1(u2r1) * f22a
            f2b = j0(u2r1) * f22b - y0(u2r1) * f22a
        else:  # a
            f22a, f22b = i0(u2r2), k0(u2r2)
            f2a = i1(u2r1) * f22b + k1(u2r1) * f22a
            f2b = i0(u2r1) * f22b - k0(u2r1) * f22a
        return f11a * n2sq * f2a - f11b * n1sq * f2b * u2r1

    def _ehcoeq(self, v0: float, nu: int) -> float:
        u1r1, u2r1, u2r2, s1, s2, n1sq, n2sq, n3sq = self.get_parameters(v0)
        if s1 == 0:
            return self.__fct3(nu, u2r1, u2r2, 2, n2sq, n3sq)
        else:
            s3 = 1 if s1 == s2 else -1

            # if n1sq > n2sq > n3sq:
            #     # s3 = 1 if nu == 1 else -1
            #     return self.__fct2(nu, u1r1, u2r1, u2r2,
            #                        s1, s2, s3,
            #                        n1sq, n2sq, n3sq)
            # else:
            return self.__fct1(nu, u1r1, u2r1, u2r2,
                               s1, s2, s3,
                               n1sq, n2sq, n3sq)

    def _hecoeq(self, v0: float, nu: int):
        u1r1, u2r1, u2r2, s1, s2, n1sq, n2sq, n3sq = self.get_parameters(v0)
        if s1 == 0:
            return self.__fct3(nu, u2r1, u2r2, -2, n2sq, n3sq)
        else:
            s3 = -1 if s1 == s2 else 1
            if n1sq > n2sq > n3sq:
                s3 = -1 if nu == 1 else 1
            #     return self.__fct1(nu, u1r1, u2r1, u2r2,
            #                        s1, s2, s3,
            #                        n1sq, n2sq, n3sq)
            # else:
            return self.__fct2(nu, u1r1, u2r1, u2r2,
                               s1, s2, s3,
                               n1sq, n2sq, n3sq)

    def __fct1(self, nu: int, u1r1: float, u2r1: float, u2r2: float, s1: float, s2: float, s3: float, n1sq: float, n2sq: float, n3sq: float) -> float:
        if s2 < 0:  # a
            b11 = iv(nu, u2r1)
            b12 = kn(nu, u2r1)
            b21 = iv(nu, u2r2)
            b22 = kn(nu, u2r2)
            b31 = iv(nu + 1, u2r1)
            b32 = kn(nu + 1, u2r1)
            f1 = b31 * b22 + b32 * b21
            f2 = b11 * b22 - b12 * b21
        else:
            b11 = jn(nu, u2r1)
            b12 = yn(nu, u2r1)
            b21 = jn(nu, u2r2)
            b22 = yn(nu, u2r2)
            if s1 == 0:
                f1 = 0
            else:
                b31 = jn(nu + 1, u2r1)
                b32 = yn(nu + 1, u2r1)
                f1 = b31 * b22 - b32 * b21
            f2 = b12 * b21 - b11 * b22
        if s1 == 0:
            delta = 1
        else:
            delta = self.__delta(nu, u1r1, u2r1, s1, s2, s3, n1sq, n2sq, n3sq)
        return f1 + f2 * delta

    def __fct2(self, nu: int, u1r1: float, u2r1: float, u2r2: float, s1: float, s2: float, s3: float, n1sq: float, n2sq: float, n3sq: float) -> float:
        with numpy.errstate(invalid='ignore'):
            delta = self.__delta(nu, u1r1, u2r1, s1, s2, s3, n1sq, n2sq, n3sq)
            n0sq = (n3sq - n2sq) / (n2sq + n3sq)

            if s2 < 0:  # a
                b11 = iv(nu, u2r1)
                b12 = kn(nu, u2r1)
                b21 = iv(nu, u2r2)
                b22 = kn(nu, u2r2)
                b31 = iv(nu + 1, u2r1)
                b32 = kn(nu + 1, u2r1)
                b41 = iv(nu - 2, u2r2)
                b42 = kn(nu - 2, u2r2)
                g1 = b11 * delta + b31
                g2 = b12 * delta - b32
                f1 = b41 * g2 - b42 * g1
                f2 = b21 * g2 - b22 * g1
            else:
                b11 = jn(nu, u2r1)
                b12 = yn(nu, u2r1)
                b21 = jn(nu, u2r2)
                b22 = yn(nu, u2r2)
                b31 = jn(nu + 1, u2r1)
                b32 = yn(nu + 1, u2r1)
                b41 = jn(nu - 2, u2r2)
                b42 = yn(nu - 2, u2r2)
                g1 = b11 * delta - b31
                g2 = b12 * delta - b32
                f1 = b41 * g2 - b42 * g1
                f2 = b22 * g1 - b21 * g2
            return f1 + n0sq * f2

    def __fct3(self, nu: int, u2r1: float, u2r2: float, dn: int, n2sq: float, n3sq: float) -> float:
        n0sq = (n3sq - n2sq) / (n2sq + n3sq)
        b11 = jn(nu, u2r1)
        b12 = yn(nu, u2r1)
        b21 = jn(nu, u2r2)
        b22 = yn(nu, u2r2)

        if dn > 0:
            b31 = jn(nu + dn, u2r1)
            b32 = yn(nu + dn, u2r1)
            f1 = b31 * b22 - b32 * b21
            f2 = b11 * b22 - b12 * b21
        else:
            b31 = jn(nu + dn, u2r2)
            b32 = yn(nu + dn, u2r2)
            f1 = b31 * b12 - b32 * b11
            f2 = b12 * b21 - b11 * b22

        return f1 - n0sq * f2
