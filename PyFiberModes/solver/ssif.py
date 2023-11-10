#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging

from PyFiberModes.solver.solver import FiberSolver
from PyFiberModes import Mode, ModeFamily


from scipy.special import jn, jn_zeros, kn, j0, j1, k0, k1, jvp, kvp
from scipy.constants import mu_0, epsilon_0, physical_constants
eta0 = physical_constants['characteristic impedance of vacuum'][0]

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
        core, clad = self.fiber.layers

        cutoff_wavelength = self.fiber.V0_to_wavelength(V0=V0)

        normal_wavelength = core.wavelength

        self.fiber.update_wavelength(cutoff_wavelength)

        n_core = core.get_maximum_index()

        n_clad = clad.get_minimum_index()

        ratio = n_core**2 / n_clad**2

        self.fiber.update_wavelength(normal_wavelength)

        return (1 + ratio) * jn(nu - 2, V0) - (1 - ratio) * jn(nu, V0)

    def find_HE_mode_cutoff(self, mode: Mode) -> float:
        if mode.m > 1:

            lower_boundary_mode = Mode(
                family=mode.family,
                nu=mode.nu,
                m=mode.m - 1
            )

            lower_neff_boundary = self.fiber.get_cutoff_V0(mode=lower_boundary_mode)

            if numpy.isnan(lower_neff_boundary) or numpy.isinf(lower_neff_boundary):
                raise AssertionError(f"find_HE_mode_cutoff: no previous cutoff for {mode} mode")

            delta = 1 / lower_neff_boundary if lower_neff_boundary else self._MCD

            lower_neff_boundary += delta
        else:
            lower_neff_boundary = delta = self._MCD

        ipoints = numpy.concatenate(
            [jn_zeros(mode.nu, mode.m), jn_zeros(mode.nu - 2, mode.m)]
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
    def solve(self, mode: Mode, delta_neff: float, lower_neff_boundary: float):
        core, clad = self.fiber.layers

        epsilon = 1e-12

        cutoff = self.fiber.get_cutoff_V0(mode=mode)

        if self.fiber.get_V0() < cutoff:
            return float("nan")

        n_core = core.get_maximum_index()

        higher_neff_boundary = numpy.sqrt(n_core**2 - (cutoff / (core.radius_out * self.wavelength.k0))**2) - epsilon

        match mode.family:
            case ModeFamily.LP:
                upper_boundary_mode = Mode(ModeFamily.LP, mode.nu + 1, mode.m)
            case ModeFamily.HE:
                upper_boundary_mode = Mode(ModeFamily.LP, mode.nu, mode.m)
            case ModeFamily.EH:
                upper_boundary_mode = Mode(ModeFamily.LP, mode.nu + 2, mode.m)
            case _:
                upper_boundary_mode = Mode(ModeFamily.LP, 1, mode.m + 1)

        cutoff = self.fiber.get_cutoff_V0(mode=upper_boundary_mode)

        try:
            value_0 = numpy.sqrt(n_core**2 - (cutoff / (core.radius_out * self.wavelength.k0))**2) + epsilon

            value_1 = clad.get_minimum_index() + epsilon

            lower_neff_boundary = max(value_0, value_1)

        except ValueError:
            lower_neff_boundary = n_core

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

        result = self.find_root_within_range(
            function=function,
            lowbound=lower_neff_boundary,
            highbound=higher_neff_boundary,
            function_kwargs=dict(nu=mode.nu)
        )

        return result

    def get_LP_field(self, nu: int, neff: float, radius: float) -> tuple:
        r"""
        Gets the LP field in the form of a tuple containing two numpy arrays.
        Tuple structure is [:math:`E_{x}`, 0, 0], [0, :math:`H_{y}`, 0].

        The field are computed with as:

        In the core:

        .. math:
            E_x &= j_0\left( U * r \ r_{core} \right) / j_0(U) \\
            H_y &= n_{eff} * \sqrt{epsilon\_0 / mu\_0} * E_x \\

        In the clad:

        .. math:
            E_x &= k_0\left( W * r \ r_{core} \right) / k_0(W) \\
            H_y &= n_{eff} * \sqrt{epsilon\_0 / mu\_0} * E_x \\

        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The effective index
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        u, w, v = self.get_U_W_V_parameter(neff=neff)

        if radius < core.radius_out:
            ex = j0(u * radius / core.radius_out) / j0(u)
        else:
            ex = k0(w * radius / core.radius_out) / k0(w)

        hy = neff * numpy.sqrt(epsilon_0 / mu_0) * ex  # Snyder & Love uses nco, but Bures uses neff

        E_field = numpy.array((ex, 0, 0))
        H_field = numpy.array((0, hy, 0))

        return E_field, H_field

    def get_TE_field(self, nu, neff: float, radius: float) -> numpy.ndarray:
        r"""
        Gets the TE field in the form of a tuple containing two numpy arrays.
        Tuple structure is [0, :math:`E_{\phi}`, 0], [:math:`H_{r}`, 0, :math:`H_{z}`]

        The field are computed within the core and radius:

        In the core

        .. math::

            H_z    &= \frac{\sqrt{epsilon\_0 / mu\_0} * U}{k_0 r_{core}} * j_0(U * r/r_{core}) / j_1(U) \\
            E_\phi &= -j_1(U * r/r_{core}) / j_1(U) \\
            H_r &= n_{eff} * \sqrt{epsilon\_0 / mu\_0} * E_\phi \\


        In the clad

        .. math::

            H_z    &= \frac{\sqrt{epsilon\_0 / mu\_0} * W}{k_0 r_{core}} * k_0(W * r/r_{core}) / k_1(U) \\
            E_\phi &= -k_1(W * r/r_{core}) / k_1(W) \\
            H_r &= n_{eff} * \sqrt{epsilon\_0 / mu\_0} * E_\phi \\

        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The effective index
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The TE field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        u, w, _ = self.get_U_W_V_parameter(neff=neff)

        term_0 = self.wavelength.k0 * core.radius_out
        ratio = radius / core.radius_out

        if radius < core.radius_out:
            hz = -numpy.sqrt(epsilon_0 / mu_0) * u / term_0 * j0(u * ratio) / j1(u)
            ephi = -j1(u * ratio) / j1(u)
        else:
            hz = numpy.sqrt(epsilon_0 / mu_0) * w / term_0 * k0(w * ratio) / k1(w)
            ephi = -k1(w * ratio) / k1(w)

        hr = -neff * numpy.sqrt(epsilon_0 / mu_0) * ephi

        E_field = numpy.array((0, ephi, 0))
        H_field = numpy.array((hr, 0, hz))

        return E_field, H_field

    def get_TM_field(self, nu, neff: float, radius: float) -> tuple:
        r"""
        Gets the TM field in the form of a tuple containing two numpy arrays.
        Tuple structure is [:math:`E_{r}`, 0, :math:`E_{z}`], [0, :math:`H_{\phi}`, 0]


        The field are computed within the core and radius:

        In the core

        .. math::

            E_z &= \frac{-U}{k_0 * n_{eff} * r_{core}} * \frac{j_0(U * r / r_{core})}{j_1(U)} \\
            E_r &= j_1(U * r/r_{core}) / j_1(U) \\
            H_\phi &= \sqrt{epsilon\_0 / mu\_0} * n_{core} / n_{eff} * E_r \\


        In the clad

        .. math::

            E_z &= \frac{n_{core}}{n_{clad}} \frac{W}{k_0 * n_{eff} * r_{core}} * \frac{k_0(W * r / r_{core})}{k_1(W)} \\
            E_r &= \frac{n_{core}}{n_{clad}} k_1(W * r/r_{core}) / k_1(W)\\
            H_\phi &= \sqrt{epsilon\_0 / mu\_0} * \frac{n_{core}}{n_{clad}} * k_1(W * r/r_{core}) / k_1(W) \\

        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The effective index
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        rho = core.radius_out

        k = self.wavelength.k0

        n_core = core.get_maximum_index(wavelength=self.wavelength)

        n_clad = clad.get_minimum_index(wavelength=self.wavelength)

        u, w, _ = self.get_U_W_V_parameter(neff=neff)

        radius_ratio = radius / rho
        index_ratio = n_core / n_clad

        if radius < rho:
            ez = -u / (k * neff * rho) * j0(u * radius_ratio) / j1(u)
            er = j1(u * radius_ratio) / j1(u)
            hphi = numpy.sqrt(epsilon_0 / mu_0) * n_core / neff * er
        else:
            ez = index_ratio * w / (k * neff * rho) * k0(w * radius_ratio) / k1(w)
            er = index_ratio * k1(w * radius_ratio) / k1(w)
            hphi = numpy.sqrt(epsilon_0 / mu_0) * index_ratio * k1(w * radius_ratio) / k1(w)

        return numpy.array((er, 0, ez)), numpy.array((0, hphi, 0))

    def get_HE_field(self, nu: float, neff: float, radius: float) -> tuple:
        r"""
        Gets the HE field in the form of a tuple containing two numpy arrays.
        Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The effective index
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The HE field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        rho = core.radius_out

        k = self.wavelength.k0

        n_core_square = core.get_maximum_index()**2

        n_clad_square = clad.get_minimum_index()**2

        u, w, v = self.get_U_W_V_parameter(neff=neff)

        jnu = jn(nu, u)
        knw = kn(nu, w)

        delta = (1 - n_clad_square / n_core_square) / 2
        b1 = jvp(nu, u) / (u * jnu)
        b2 = kvp(nu, w) / (w * knw)
        F1 = (u * w / v)**2 * (b1 + (1 - 2 * delta) * b2) / nu
        F2 = (v / (u * w))**2 * nu / (b1 + b2)
        a1 = (F2 - 1) / 2
        a2 = (F2 + 1) / 2
        a3 = (F1 - 1) / 2
        a4 = (F1 + 1) / 2
        a5 = (F1 - 1 + 2 * delta) / 2
        a6 = (F1 + 1 - 2 * delta) / 2

        if radius < rho:
            term_0 = u * radius / rho

            jmur = jn(nu - 1, term_0)
            jpur = jn(nu + 1, term_0)
            jnur = jn(nu, term_0)

            er = -(a1 * jmur + a2 * jpur) / jnu
            ephi = -(a1 * jmur - a2 * jpur) / jnu
            ez = u / (k * neff * rho) * jnur / jnu
            hr = numpy.sqrt(epsilon_0 / mu_0) * n_core_square / neff * (a3 * jmur - a4 * jpur) / jnu
            hphi = -numpy.sqrt(epsilon_0 / mu_0) * n_core_square / neff * (a3 * jmur + a4 * jpur) / jnu
            hz = numpy.sqrt(epsilon_0 / mu_0) * u * F2 / (k * rho) * jnur / jnu
        else:
            term_1 = w * radius / rho

            kmur = kn(nu - 1, term_1)
            kpur = kn(nu + 1, term_1)
            knur = kn(nu, term_1)

            er = -u / w * (a1 * kmur - a2 * kpur) / knw
            ephi = -u / w * (a1 * kmur + a2 * kpur) / knw
            ez = u / (k * neff * rho) * knur / knw
            hr = numpy.sqrt(epsilon_0 / mu_0) * n_core_square / neff * u / w * (a5 * kmur + a6 * kpur) / knw
            hphi = -numpy.sqrt(epsilon_0 / mu_0) * n_core_square / neff * u / w * (a5 * kmur - a6 * kpur) / knw
            hz = numpy.sqrt(epsilon_0 / mu_0) * u * F2 / (k * rho) * knur / knw

        E_field = numpy.array((er, ephi, ez))
        H_field = numpy.array((hr, hphi, hz))

        return E_field, H_field

    def get_EH_field(self, *args, **kwargs) -> tuple:
        r"""
        Gets the EH field in the form of a tuple containing two numpy arrays.
        Tuple structure is [:math:`E_{r}`, :math:`E_{\phi}`, :math:`E_{z}`], [:math:`H_{r}`, :math:`H_{\phi}`, :math:`H_{z}`]

        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The effective index
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The lp field.
        :rtype:     tuple
        """
        return self.get_HE_field(*args, **kwargs)

    def get_U_W_V_parameter(self, neff: float) -> tuple:
        r"""
        Gets the U, W parameter of the fiber. Those are computed as:

        .. math:

            U &= r_{core} * k_0 * \sqrt{n_{core}^2 - n_{eff}^2} \\
            W &= r_{core} * k_0 * \sqrt{n_{eff}^2 - n_{core}^2} \\
            V &= \sqrt{U^2 + W^2} \\

        :param      neff:        The effective index
        :type       neff:        float

        :returns:   The U and W parameter.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        n_core = core.get_maximum_index()

        n_clad = clad.get_minimum_index()

        U = core.radius_out * self.wavelength.k0 * numpy.sqrt(n_core**2 - neff**2)
        W = core.radius_out * self.wavelength.k0 * numpy.sqrt(neff**2 - n_clad**2)
        V = numpy.sqrt(U**2 + W**2)

        return U, W, V

    def _lpceq(self, neff: float, nu: int) -> float:
        """
        I don't know what it returns.

        .. math::
            U * j_{\nu -1}(U) * k_{\nu}(W) + W * j_{\nu}(U) * k_{\nu - 1}(W)

        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int

        :returns:   Dont know
        :rtype:     float
        """
        u, w, _ = self.get_U_W_V_parameter(neff=neff)

        return u * jn(nu - 1, u) * kn(nu, w) + w * jn(nu, u) * kn(nu - 1, w)

    def _teceq(self, neff: float, nu: int) -> float:
        """
        I don't know what it returns.

        .. math::
            U * j_0(U) * k_1(W) + W * j_1(U) * k_0(W)

        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int

        :returns:   Dont know
        :rtype:     float
        """
        U, W, _ = self.get_U_W_V_parameter(neff=neff)

        return U * j0(U) * k1(W) + W * j1(U) * k0(W)

    def _tmceq(self, neff: float, nu: int) -> float:
        """
        I don't know what it returns.

        .. math::
            U * j_0(U) * k_1(W) * n_{clad}^2 + W * j_1(U) * k_0(W) * n_{core}^2

        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int

        :returns:   Dont know
        :rtype:     float
        """
        core, clad = self.fiber.layers

        U, W, V = self.get_U_W_V_parameter(neff=neff)

        n_core = core.get_maximum_index(wavelength=self.wavelength)

        n_clad = clad.get_minimum_index(wavelength=self.wavelength)

        return U * j0(U) * k1(W) * n_clad**2 + W * j1(U) * k0(W) * n_core**2

    def get_HE_EH_terms(self, neff, nu: int) -> float:
        """
        I don't know what it returns.

        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int

        :returns:   Dont know
        :rtype:     float
        """
        core, clad = self.fiber.layers

        u, w, v = self.get_U_W_V_parameter(neff=neff)

        n_core = core.get_maximum_index()

        n_clad = clad.get_minimum_index()

        delta = (1 - n_clad**2 / n_core**2) / 2
        jnu = jn(nu, u)
        knu = kn(nu, w)
        kp = kvp(nu, w)

        term_0 = jvp(nu, u) * w * knu + kp * u * jnu * (1 - delta)
        term_1 = (nu * neff * v**2 * knu)
        term_2 = n_core * u * w
        term_3 = u * kp * delta
        term_4 = jnu * numpy.sqrt(term_3**2 + (term_1 / term_2)**2)

        return term_0, term_4

    def _heceq(self, neff: float, nu: int) -> float:
        """
        I don't know what it returns.

        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int

        :returns:   Dont know
        :rtype:     float
        """
        term_0, term_1 = self.get_HE_EH_terms(
            neff=neff,
            nu=nu
        )

        return term_0 + term_1

    def _ehceq(self, neff: float, nu: int) -> float:
        """
        I don't know what it returns.

        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int

        :returns:   Dont know
        :rtype:     float
        """
        term_0, term_1 = self.get_HE_EH_terms(
            neff=neff,
            nu=nu
        )

        return term_0 - term_1


# -
