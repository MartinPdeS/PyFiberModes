#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging

from PyFiberModes.solver.solver import FiberSolver
from PyFiberModes import Mode, ModeFamily


from scipy.special import jn, jn_zeros, kn, j0, j1, k0, k1, jvp, kvp
from PyFiberModes.constants import Y0

"""
Solver for standard layer step-index solver: SSIF
"""


class CutoffSolver(FiberSolver):
    """
    Cutoff solver for standard step-index fiber.
    """
    logger = logging.getLogger(__name__)

    def solve(self, mode: Mode) -> float:
        nu = mode.nu
        m = mode.m

        if mode.family is ModeFamily.LP:
            if nu == 0:
                nu, m = 1, -1

            else:
                nu -= 1
        elif mode.family is ModeFamily.HE:
            if nu == 1:
                m -= 1
            else:
                return self.find_HE_mode_cutoff(mode)

        return jn_zeros(nu, m)[m - 1]

    def get_cutoff_HE(self, V0: float, nu: int) -> float:
        wavelength = self.fiber.V0_to_wavelength(V0=V0)

        core, clad = self.fiber.layers

        max_index_core = core.get_maximum_index(wavelength=wavelength)

        min_index_clad = clad.get_minimum_index(wavelength=wavelength)

        ratio = max_index_core**2 / min_index_clad**2

        return (1 + ratio) * jn(nu - 2, V0) - (1 - ratio) * jn(nu, V0)

    def find_HE_mode_cutoff(self, mode: Mode) -> float:
        if mode.m > 1:

            lower_boundary_mode = Mode(
                family=mode.family,
                nu=mode.nu,
                m=mode.m - 1
            )

            lower_neff_boundary = self.fiber.get_cutoff(mode=lower_boundary_mode)

            if numpy.isnan(lower_neff_boundary) or numpy.isinf(lower_neff_boundary):
                raise AssertionError(f"find_HE_mode_cutoff: no previous cutoff for {mode} mode")

            delta = 1 / lower_neff_boundary if lower_neff_boundary else self._MCD

            lower_neff_boundary += delta
        else:
            lower_neff_boundary = delta = self._MCD

        ipoints = numpy.concatenate(
            (jn_zeros(mode.nu, mode.m), jn_zeros(mode.nu - 2, mode.m))
        )

        ipoints.sort()
        ipoints = list(ipoints[ipoints > lower_neff_boundary])

        cutoff = self.find_function_first_root(
            function=self.get_cutoff_HE,
            function_args=(mode.nu,),
            lowbound=lower_neff_boundary,
            ipoints=ipoints,
            delta=delta
        )

        if numpy.isnan(cutoff):
            self.logger.error(f"find_HE_mode_cutoff: no cutoff found for {mode} mode")
            return 0

        return cutoff


class NeffSolver(FiberSolver):
    """
    Effective index solver for standard step-index fiber
    """

    def solve(self, wavelength: float, mode: Mode, delta_neff: float, lower_neff_boundary: float):
        core, clad = self.fiber.layers

        epsilon = 1e-12

        cutoff = self.fiber.get_cutoff(mode=mode)

        if self.fiber.get_V0(wavelength=wavelength) < cutoff:
            return float("nan")

        max_core_index = core.get_maximum_index(wavelength=wavelength)

        higher_neff_boundary = numpy.sqrt(max_core_index**2 - (cutoff / (core.radius_out * wavelength.k0))**2) - epsilon

        match mode.family:
            case ModeFamily.LP:
                upper_boundary_mode = Mode(ModeFamily.LP, mode.nu + 1, mode.m)
            case ModeFamily.HE:
                upper_boundary_mode = Mode(ModeFamily.LP, mode.nu, mode.m)
            case ModeFamily.EH:
                upper_boundary_mode = Mode(ModeFamily.LP, mode.nu + 2, mode.m)
            case _:
                upper_boundary_mode = Mode(ModeFamily.LP, 1, mode.m + 1)

        cutoff = self.fiber.get_cutoff(mode=upper_boundary_mode)

        try:
            value_0 = numpy.sqrt(max_core_index**2 - (cutoff / (core.radius_out * wavelength.k0))**2) + epsilon

            value_1 = clad.get_minimum_index(wavelength=wavelength) + epsilon

            lower_neff_boundary = max(value_0, value_1)

        except ValueError:
            lower_neff_boundary = max_core_index

        match mode.family:
            case ModeFamily.LP:
                function = self._lpceq
            case ModeFamily.TE:
                function = self._teceq
            case ModeFamily.TM:
                function = self._tmceq
            case ModeFamily.EH:
                function = self._ehceq
            case ModeFamily.HE:
                function = self._heceq

        result = self._findBetween(
            function=function,
            lowbound=lower_neff_boundary,
            highbound=higher_neff_boundary,
            function_args=(wavelength, mode.nu)
        )

        return result

    def get_LP_field(self, wavelength: float, nu: int, neff: float, radius: float) -> tuple:
        """
        Gets the LP field in the form of a tuple containing two numpy arrays.
        Tuple structure is [ex, 0, 0], [0, hy, 0]

        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The neff
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        max_index_core = core.get_maximum_index(wavelength=wavelength)

        min_index_clad = clad.get_minimum_index(wavelength=wavelength)

        u = core.radius_out * wavelength.k0 * numpy.sqrt(max_index_core**2 - neff**2)
        w = core.radius_out * wavelength.k0 * numpy.sqrt(neff**2 - min_index_clad**2)

        if radius < core.radius_out:
            ex = j0(u * radius / core.radius_out) / j0(u)
        else:
            ex = k0(w * radius / core.radius_out) / k0(w)
        hy = neff * Y0 * ex  # Snyder & Love uses nco, but Bures uses neff

        return numpy.array((ex, 0, 0)), numpy.array((0, hy, 0))

    def get_TE_field(self, wavelength: float, nu, neff: float, radius: float) -> numpy.ndarray:
        r"""
        Gets the TE field in the form of a tuple containing two numpy arrays.
        Tuple structure is [0, E$_{\phi}$, 0], [H$_{r}$, 0, H$_{z}$]

        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The neff
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        max_index_core = core.get_maximum_index(wavelength=wavelength)

        min_index_clad = clad.get_minimum_index(wavelength=wavelength)

        u = core.radius_out * wavelength.k0 * numpy.sqrt(max_index_core**2 - neff**2)
        w = core.radius_out * wavelength.k0 * numpy.sqrt(neff**2 - min_index_clad**2)

        term_0 = wavelength.k0 * core.radius_out
        ratio = radius / core.radius_out

        if radius < core.radius_out:
            hz = -Y0 * u / term_0 * j0(u * ratio) / j1(u)
            ephi = -j1(u * ratio) / j1(u)
        else:
            hz = Y0 * w / term_0 * k0(w * ratio) / k1(w)
            ephi = -k1(w * ratio) / k1(w)

        hr = -neff * Y0 * ephi

        return numpy.array((0, ephi, 0)), numpy.array((hr, 0, hz))

    def get_TM_field(self, wavelength: float, nu, neff: float, radius: float) -> tuple:
        r"""
        Gets the TM field in the form of a tuple containing two numpy arrays.
        Tuple structure is [E$_{r}$, 0, E$_{z}$], [0, H$_{\phi}$, 0]

        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The neff
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        rho = core.radius_out

        k = wavelength.k0

        max_index_core = core.get_maximum_index(wavelength=wavelength)

        min_index_clad = clad.get_minimum_index(wavelength=wavelength)

        u = rho * wavelength.k0 * numpy.sqrt(max_index_core**2 - neff**2)
        w = rho * wavelength.k0 * numpy.sqrt(neff**2 - min_index_clad**2)

        radius_ratio = radius / rho
        index_ratio = max_index_core / min_index_clad

        if radius < rho:
            ez = -u / (k * neff * rho) * j0(u * radius_ratio) / j1(u)
            er = j1(u * radius_ratio) / j1(u)
            hphi = Y0 * max_index_core / neff * er
        else:
            ez = index_ratio * w / (k * neff * rho) * k0(w * radius_ratio) / k1(w)
            er = index_ratio * k1(w * radius_ratio) / k1(w)
            hphi = Y0 * index_ratio * k1(w * radius_ratio) / k1(w)

        return numpy.array((er, 0, ez)), numpy.array((0, hphi, 0))

    def get_HE_field(self, wavelength: float, nu: float, neff: float, radius: float) -> tuple:
        r"""
        Gets the HE field in the form of a tuple containing two numpy arrays.
        Tuple structure is [E$_{r}$, E$_{\phi}$, E$_{z}$], [H$_{r}$, H$_{\phi}$, H$_{z}$]

        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The neff
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        rho = core.radius_out

        k = wavelength.k0

        max_index_core_squared = core.get_maximum_index(wavelength=wavelength)**2

        max_index_clad_squared = clad.get_minimum_index(wavelength=wavelength)**2

        u = rho * k * numpy.sqrt(max_index_core_squared - neff**2)
        w = rho * k * numpy.sqrt(neff**2 - max_index_clad_squared)
        v = rho * k * numpy.sqrt(max_index_core_squared - max_index_clad_squared)

        jnu = jn(nu, u)
        knw = kn(nu, w)

        Delta = (1 - max_index_clad_squared / max_index_core_squared) / 2
        b1 = jvp(nu, u) / (u * jnu)
        b2 = kvp(nu, w) / (w * knw)
        F1 = (u * w / v)**2 * (b1 + (1 - 2 * Delta) * b2) / nu
        F2 = (v / (u * w))**2 * nu / (b1 + b2)
        a1 = (F2 - 1) / 2
        a2 = (F2 + 1) / 2
        a3 = (F1 - 1) / 2
        a4 = (F1 + 1) / 2
        a5 = (F1 - 1 + 2 * Delta) / 2
        a6 = (F1 + 1 - 2 * Delta) / 2

        if radius < rho:
            term_0 = u * radius / rho

            jmur = jn(nu - 1, term_0)
            jpur = jn(nu + 1, term_0)
            jnur = jn(nu, term_0)

            er = -(a1 * jmur + a2 * jpur) / jnu
            ephi = -(a1 * jmur - a2 * jpur) / jnu
            ez = u / (k * neff * rho) * jnur / jnu
            hr = Y0 * max_index_core_squared / neff * (a3 * jmur - a4 * jpur) / jnu
            hphi = -Y0 * max_index_core_squared / neff * (a3 * jmur + a4 * jpur) / jnu
            hz = Y0 * u * F2 / (k * rho) * jnur / jnu
        else:
            term_1 = w * radius / rho

            kmur = kn(nu - 1, term_1)
            kpur = kn(nu + 1, term_1)
            knur = kn(nu, term_1)

            er = -u / w * (a1 * kmur - a2 * kpur) / knw
            ephi = -u / w * (a1 * kmur + a2 * kpur) / knw
            ez = u / (k * neff * rho) * knur / knw
            hr = Y0 * max_index_core_squared / neff * u / w * (a5 * kmur + a6 * kpur) / knw
            hphi = -Y0 * max_index_core_squared / neff * u / w * (a5 * kmur - a6 * kpur) / knw
            hz = Y0 * u * F2 / (k * rho) * knur / knw

        return numpy.array((er, ephi, ez)), numpy.array((hr, hphi, hz))

    def get_EH_field(self, *args, **kwargs) -> tuple:
        r"""
        Gets the EH field in the form of a tuple containing two numpy arrays.
        Tuple structure is [E$_{r}$, E$_{\phi}$, E$_{z}$], [H$_{r}$, H$_{\phi}$, H$_{z}$]

        :param      wavelength:  The wavelength
        :type       wavelength:  float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The neff
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        return self.get_HE_field(*args, **kwargs)

    def get_parameter_uw(self, wavelength: float, neff: float) -> tuple:
        core, clad = self.fiber.layers

        max_index = core.get_maximum_index(wavelength=wavelength)
        min_index = clad.get_minimum_index(wavelength=wavelength)

        term_0 = core.radius_out * wavelength.k0 * numpy.sqrt(max_index**2 - neff**2)
        term_1 = core.radius_out * wavelength.k0 * numpy.sqrt(neff**2 - min_index**2)

        return term_0, term_1

    def _lpceq(self, neff: float, wavelength: float, nu: float) -> float:
        u, w = self.get_parameter_uw(
            wavelength=wavelength,
            neff=neff
        )

        return (u * jn(nu - 1, u) * kn(nu, w) + w * jn(nu, u) * kn(nu - 1, w))

    def _teceq(self, neff: float, wavelength: float, nu: float) -> float:
        u, w = self.get_parameter_uw(
            wavelength=wavelength,
            neff=neff
        )

        return u * j0(u) * k1(w) + w * j1(u) * k0(w)

    def _tmceq(self, neff: float, wavelength: float, nu: float) -> float:
        u, w = self.get_parameter_uw(
            wavelength=wavelength,
            neff=neff
        )

        max_index_core = self.fiber.get_maximum_index(
            layer_idx=0,
            wavelength=wavelength
        )

        min_index_clad = self.fiber.get_minimum_index(
            layer_idx=1,
            wavelength=wavelength
        )

        return (u * j0(u) * k1(w) * min_index_clad**2 + w * j1(u) * k0(w) * max_index_core**2)

    def _heceq(self, neff: float, wavelength: float, nu: int):
        u, w = self.get_parameter_uw(
            wavelength=wavelength,
            neff=neff
        )

        v2 = u**2 + w**2

        core = self.fiber.layers[0]
        clad = self.fiber.layers[1]

        nco = core.get_maximum_index(wavelength=wavelength)

        ncl = clad.get_minimum_index(wavelength=wavelength)

        delta = (1 - ncl**2 / nco**2) / 2
        jnu = jn(nu, u)
        knu = kn(nu, w)
        kp = kvp(nu, w)

        term_0 = jvp(nu, u) * w * knu + kp * u * jnu * (1 - delta)
        term_1 = (nu * neff * v2 * knu)
        term_2 = nco * u * w
        term_3 = u * kp * delta

        return term_0 + jnu * numpy.sqrt(term_3**2 + (term_1 / term_2)**2)

    def _ehceq(self, neff: float, wavelength: float, nu: int):
        u, w = self.get_parameter_uw(
            wavelength=wavelength,
            neff=neff
        )

        v2 = u**2 + w**2

        nco = self.fiber.get_maximum_index(
            layer_idx=0,
            wavelength=wavelength
        )

        ncl = self.fiber.get_minimum_index(
            layer_idx=1,
            wavelength=wavelength
        )

        delta = (1 - ncl**2 / nco**2) / 2
        jnu = jn(nu, u)
        knu = kn(nu, w)
        kp = kvp(nu, w)

        return (jvp(nu, u) * w * knu + kp * u * jnu * (1 - delta) - jnu * numpy.sqrt((u * kp * delta)**2 + ((nu * neff * v2 * knu) / (nco * u * w))**2))
