#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging

from PyFiberModes.solver.base_solver import BaseSolver
from PyFiberModes.mode import Mode


from scipy.special import jn, kn, j0, j1, k0, k1, jvp, kvp
from scipy.constants import mu_0, epsilon_0, physical_constants
eta0 = physical_constants['characteristic impedance of vacuum'][0]
Y0 = numpy.sqrt(epsilon_0 / mu_0)

"""
Solver for standard layer step-index solver: SSIF
"""


class NeffSolver(BaseSolver):
    """
    Effective index solver for standard step-index fiber
    """

    def get_fiber_delta(self) -> float:
        r"""
        Gets the fiber delta parameter defined parge 77 of Jacques Bures

        .. math::
            2 \Delta = 1 -\frac{n_clad^2}{n_core^2}

        :returns:   The fiber delta parameter.
        :rtype:     float
        """
        core, clad = self.fiber.layers
        n_clad = clad.refractive_index
        n_core = core.refractive_index

        return 0.5 * (1 - n_clad**2 / n_core**2)

    def get_mode_with_lower_neff(self, mode: Mode) -> Mode:
        match mode.family:
            case 'LP':
                lower_neff_mode = Mode('LP', mode.nu + 1, mode.m)
            case 'HE':
                lower_neff_mode = Mode('LP', mode.nu, mode.m)
            case 'EH':
                lower_neff_mode = Mode('LP', mode.nu + 2, mode.m)
            case _:
                lower_neff_mode = Mode('LP', 1, mode.m + 1)

        return lower_neff_mode

    def get_clad_index_from_V0(self, V0: float) -> float:
        """
        Gets a clad index value associated to a certain V0 parameter for the same exact fiber.

        :param      V0:   The V number
        :type       V0:   float

        :returns:   The clad index
        :rtype:     float
        """
        core = self.fiber.first_layer

        n_core = core.refractive_index

        NA = V0 / (self.wavelength.k0 * core.radius_out)

        n_clad_2 = n_core**2 - NA**2

        if n_clad_2 < 0:
            return numpy.nan

        return numpy.sqrt(n_clad_2)

    def get_low_neff_boundary(self, mode: Mode) -> float:
        """
        Gets the lower neff boundary using a lower neff mode.

        :param      mode:  The current mode
        :type       mode:  Mode

        :returns:   The low neff boundary.
        :rtype:     float
        """
        core, clad = self.fiber.layers

        lower_neff_mode = self.get_mode_with_lower_neff(mode=mode)

        lower_neff_cutoff_V0 = self.fiber.get_cutoff_v0(mode=lower_neff_mode)

        lower_neff_clad_index = self.get_clad_index_from_V0(V0=lower_neff_cutoff_V0)

        lower_neff_boundary = max(lower_neff_clad_index, clad.refractive_index)

        if numpy.isnan(lower_neff_boundary):
            lower_neff_boundary = clad.refractive_index

        return lower_neff_boundary

    def get_ceq_function(self, mode: Mode) -> object:
        """
        Gets the adequat phase matching mode equation function.

        :param      mode:  The mode
        :type       mode:  Mode

        :returns:   The ceq function.
        :rtype:     object
        """
        match mode.family:
            case 'LP':
                return self.get_LP_equation
            case 'TE':
                return self.get_TE_equation
            case 'TM':
                return self.get_TM_equation
            case 'EH':
                return self.get_EH_equation
            case 'HE':
                return self.get_HE_equation

    def solve(self,
            mode: Mode,
            delta_neff: float,
            max_iteration: int = 16,
            epsilon: float = 1e-12) -> float:
        """
        Solve and return the effective index (neff) for a given mode.

        :param      mode:                 The mode
        :type       mode:                 Mode
        :param      delta_neff:           The delta neff
        :type       delta_neff:           float
        :param      max_iteration:        The maximum iteration, compute time heavily depend on it.
        :type       max_iteration:        int
        :param      epsilon:              The epsilon
        :type       epsilon:              float

        :returns:   The effective index of the mode
        :rtype:     float
        """
        mode_cutoff_V0 = self.fiber.get_cutoff_v0(mode=mode)

        cutoff_wavelength = self.fiber.get_cutoff_wavelength(mode=mode)

        if cutoff_wavelength < self.wavelength:
            logging.warning(f"Mode: {mode} cutoff wavelenght: {cutoff_wavelength} is below actual wavelength: {self.fiber.wavelength}")
            return numpy.nan

        fiber_cutoff = self.fiber.get_V0()

        if fiber_cutoff < mode_cutoff_V0:
            return numpy.nan

        n_clad_equivalent = self.get_clad_index_from_V0(V0=mode_cutoff_V0)  # High neff boundary

        lower_neff_boundary = self.get_low_neff_boundary(mode=mode)

        if n_clad_equivalent < lower_neff_boundary:
            raise ValueError(f"Error in computation, most probably the given mode: {mode} does not exist in that configuration")

        function = self.get_ceq_function(mode=mode)

        result = self.find_root_within_range(
            function=function,
            x_low=lower_neff_boundary + epsilon,
            x_high=n_clad_equivalent - epsilon,
            function_args=(mode.nu, ),
            max_iteration=max_iteration
        )

        return result

    def get_LP_field(self, nu: int, neff: float, radius: float) -> tuple:
        r"""
        Gets the LP field in the form of a tuple containing two numpy arrays.
        Tuple structure is [:math:`E_{x}`, 0, 0], [0, :math:`H_{y}`, 0].

        The field are computed with as:

        In the core:

        .. math::
            E_x &= j_0\left( U * r \ r_{core} \right) / j_0(U) \\[10pt]
            H_y &= n_{eff} * \sqrt{\epsilon_0 / \mu_0} * E_x \\[10pt]

        In the clad:

        .. math::
            E_x &= k_0\left( W * r \ r_{core} \right) / k_0(W) \\[10pt]
            H_y &= n_{eff} * \sqrt{\epsilon_0 / \mu_0} * E_x \\[10pt]

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

        e_field = numpy.array((ex, 0, 0))
        h_field = numpy.array((0, hy, 0))

        return e_field, h_field

    def get_TE_field(self, nu: int, neff: float, radius: float) -> numpy.ndarray:
        r"""
        Gets the TE field in the form of a tuple containing two numpy arrays.
        Tuple structure is [0, :math:`E_{\phi}`, 0], [:math:`H_{r}`, 0, :math:`H_{z}`]

        The field are computed within the core and radius:

        In the core

        .. math::
            H_z &= \frac{\sqrt{\epsilon_0 / \mu_0} * U}{k_0 r_{core}} * \frac{j_0(U * r/r_{core})}{j_1(U)} \\[10pt]
            E_\phi &= -j_1(U * r/r_{core}) / j_1(U) \\[10pt]
            H_r &= n_{eff} * \sqrt{\epsilon_0 / \mu_0} * E_\phi \\[10pt]


        In the clad

        .. math::
            H_z &= \frac{\sqrt{\epsilon_0 / \mu_0} * W}{k_0 r_{core}} * \frac{k_0(W * r/r_{core})}{k_1(U)} \\[10pt]
            E_\phi &= -k_1(W * r/r_{core}) / k_1(W) \\[10pt]
            H_r &= n_{eff} * \sqrt{\epsilon_0 / \mu_0} * E_\phi \\[10pt]

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

        e_field = numpy.array((0, ephi, 0))
        h_field = numpy.array((hr, 0, hz))

        return e_field, h_field

    def get_TM_field(self, nu: int, neff: float, radius: float) -> tuple:
        r"""
        Gets the TM field in the form of a tuple containing two numpy arrays.
        Tuple structure is [:math:`E_{r}`, 0, :math:`E_{z}`], [0, :math:`H_{\phi}`, 0]


        The field are computed within the core and radius:

        In the core

        .. math::
            E_z &= \frac{-U}{k_0 * n_{eff} * r_{core}} * \frac{j_0(U * r / r_{core})}{j_1(U)} \\[10pt]
            E_r &= j_1(U * r/r_{core}) / j_1(U) \\[10pt]
            H_\phi &= \sqrt{\epsilon_0 / \mu_0} * n_{core} / n_{eff} * E_r \\[10pt]


        In the clad

        .. math::
            E_z &= \frac{n_{core}}{n_{clad}} \frac{W}{k_0 * n_{eff} * r_{core}} * \frac{k_0(W * r / r_{core})}{k_1(W)} \\[10pt]
            E_r &= \frac{n_{core}}{n_{clad}} k_1(W * r/r_{core}) / k_1(W)\\[10pt]
            H_\phi &= \sqrt{\epsilon_0 / \mu_0} * \frac{n_{core}}{n_{clad}} * k_1(W * r/r_{core}) / k_1(W) \\[10pt]

        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The effective index
        :type       neff:        float
        :param      radius:      The radius
        :type       radius:      float

        :returns:   The LP field.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        rho = core.radius_out

        k = self.wavelength.k0

        n_core = core.refractive_index

        n_clad = clad.refractive_index

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

        e_field = numpy.array((er, 0, ez))
        h_field = numpy.array((0, hphi, 0))

        return e_field, h_field

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

        n_core_square = core.refractive_index**2

        n_clad_square = clad.refractive_index**2

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
            hr = Y0 * n_core_square / neff * (a3 * jmur - a4 * jpur) / jnu
            hphi = -Y0 * n_core_square / neff * (a3 * jmur + a4 * jpur) / jnu
            hz = Y0 * u * F2 / (k * rho) * jnur / jnu
        else:
            term_1 = w * radius / rho

            kmur = kn(nu - 1, term_1)
            kpur = kn(nu + 1, term_1)
            knur = kn(nu, term_1)

            er = -u / w * (a1 * kmur - a2 * kpur) / knw
            ephi = -u / w * (a1 * kmur + a2 * kpur) / knw
            ez = u / (k * neff * rho) * knur / knw
            hr = Y0 * n_core_square / neff * u / w * (a5 * kmur + a6 * kpur) / knw
            hphi = -Y0 * n_core_square / neff * u / w * (a5 * kmur - a6 * kpur) / knw
            hz = Y0 * u * F2 / (k * rho) * knur / knw

        e_field = numpy.array((er, ephi, ez))
        h_field = numpy.array((hr, hphi, hz))

        return e_field, h_field

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

        :returns:   The LP field.
        :rtype:     tuple
        """
        return self.get_HE_field(*args, **kwargs)

    def get_U_W_V_parameter(self, neff: float) -> tuple:
        r"""
        Gets the U, W parameter of the fiber. Those are computed as:

        .. math:

            U &= r_{core} * k_0 * \sqrt{n_{core}^2 - n_{eff}^2} \\[10pt]
            W &= r_{core} * k_0 * \sqrt{n_{eff}^2 - n_{core}^2} \\[10pt]
            V &= \sqrt{U^2 + W^2} \\[10pt]

        :param      neff:        The effective index
        :type       neff:        float

        :returns:   The U and W parameter.
        :rtype:     tuple
        """
        core, clad = self.fiber.layers

        n_core = core.refractive_index

        n_clad = clad.refractive_index

        U = core.radius_out * self.wavelength.k0 * numpy.sqrt(n_core**2 - neff**2)
        W = core.radius_out * self.wavelength.k0 * numpy.sqrt(neff**2 - n_clad**2)
        V = numpy.sqrt(U**2 + W**2)

        return U, W, V

    def get_LP_equation(self, neff: float, nu: int) -> float:
        """
        Return the value of the phase matching equation for LP mode.

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

        value = u * jn(nu - 1, u) * kn(nu, w) + w * jn(nu, u) * kn(nu - 1, w)

        return value

    def get_TE_equation(self, neff: float, nu: int) -> float:
        """
        Return the value of the phase matching equation for TE mode.

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

    def get_TM_equation(self, neff: float, nu: int) -> float:
        """
        Return the value of the phase matching equation for TM mode.

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

        n_core = core.refractive_index

        n_clad = clad.refractive_index

        return U * j0(U) * k1(W) * n_clad**2 + W * j1(U) * k0(W) * n_core**2

    def get_HE_EH_terms(self, neff, nu: int) -> float:
        """
        Return the value of the terms for the equation for HE or EH mode.

        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int

        :returns:   Dont know
        :rtype:     float
        """
        core, clad = self.fiber.layers

        U, W, V = self.get_U_W_V_parameter(neff=neff)

        n_core = core.refractive_index

        n_clad = clad.refractive_index

        delta = (1 - n_clad**2 / n_core**2) / 2
        jnu = jn(nu, U)
        knu = kn(nu, W)
        kp = kvp(nu, W)

        term_0 = jvp(nu, U) * W * knu + kp * U * jnu * (1 - delta)
        term_1 = (nu * neff * V**2 * knu)
        term_2 = n_core * U * W
        term_3 = U * kp * delta
        term_4 = jnu * numpy.sqrt(term_3**2 + (term_1 / term_2)**2)

        return term_0, term_4

    def get_HE_equation(self, neff: float, nu: int) -> float:
        """
        Return the value of the phase matching equation for HE mode.

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

        value = term_0 + term_1

        return value

    def get_EH_equation(self, neff: float, nu: int) -> float:
        """
        Return the value of the phase matching equation for EH mode.

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
