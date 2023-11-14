#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from dataclasses import dataclass

from scipy.special import jn, yn, iv, kn
from scipy.special import j0, y0, i0, k0
from scipy.special import j1, y1, i1, k1
from scipy.special import jvp, yvp, ivp, kvp
from scipy.constants import mu_0, epsilon_0, physical_constants

eta0 = physical_constants['characteristic impedance of vacuum'][0]


@dataclass
class Geometry(object):
    radius_in: float
    """ Minimum radius of the structure """
    radius_out: float
    """ Maximum radius of the structure """
    index_list: list
    """ Refractive index of the structure """

    def __hash__(self):
        return hash((self.radius_in, self.radius_out, self.index_list))

    def __post_init__(self):
        self.index_list = self.index_list[0]
        self.refractive_index = self.index_list
        self.thickness = self.radius_out - self.radius_in


class StepIndex(Geometry):
    DEFAULT_PARAMS = []

    def get_index_at_radius(self, radius: float) -> float:
        """
        Gets the index of the local layer at a given radius.

        :param      radius:      The radius for evaluation
        :type       radius:      float

        :returns:   The index at given radius.
        :rtype:     float
        """
        if self.radius_in <= abs(radius) <= self.radius_out:
            return self.refractive_index
        else:
            return None

    def get_U_W_parameter(self, radius: float, neff: float) -> float:
        r"""
        Gets the u parameter as defined as:

        .. math::
            U &= k \rho \sqrt{ n_{core}^2 - n_{eff}^2 } \, \text{[in the core]} \\

            W &= k \rho \sqrt{ n_{neff}^2 - n_{clad}^2 } \, \text{[in the clad]} \\

        Reference: Eq: 3.2 of Jacques Bures "Optique Guidée : Fibres Optiques et Composants Passifs Tout-Fibre"

        :param      radius:      The radius
        :type       radius:      float
        :param      neff:        The neff
        :type       neff:        float

        :returns:   The u parameter at given radius.
        :rtype:     float
        """
        index = self.get_index_at_radius(radius=radius)

        U = self.wavelength.k0 * radius * numpy.sqrt(abs(index**2 - neff**2))

        return U

    def get_psi(self, radius: float, neff: float, nu: int, C: list) -> tuple:
        r"""
        Return the :math:`\psi` function

        Bessel function definition used here (from https://docs.scipy.org/doc/scipy/reference/special.html)
        | jn:  Bessel function of first kind.
        | jvp: Derivative of Bessel function of first kind.
        | iv:  Modified Bessel function of the first kind of real order.
        | ivp: Derivatives of modified Bessel functions of the first kind.

        Reference: Eq: 3.67 of Jacques Bures "Optique Guidée : Fibres Optiques et Composants Passifs Tout-Fibre"

        :param      radius:      The radius at which the field is evaluated
        :type       radius:      float
        :param      neff:        The effective index
        :type       neff:        float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      C:           I don't know.
        :type       C:           list

        :returns:   A tuple with psi (:math:`\psi`) and the derivative of psi (:math:`\dot{\psi}`)
        :rtype:     tuple
        """
        u = self.get_U_W_parameter(
            radius=radius,
            neff=neff,
        )

        layer_max_index = self.refractive_index

        if neff < layer_max_index:
            if C[1]:
                psi = C[0] * jn(nu, u) + C[1] * yn(nu, u)
                psip = u * C[0] * jvp(nu, u) + C[1] * yvp(nu, u)
            else:
                psi = C[0] * jn(nu, u)
                psip = u * C[0] * jvp(nu, u)

        else:
            if C[1]:
                psi = C[0] * iv(nu, u) + C[1] * kn(nu, u)
                psip = u * C[0] * ivp(nu, u) + C[1] * kvp(nu, u)
            else:
                psi = C[0] * iv(nu, u)
                psip = u * C[0] * ivp(nu, u)

        return psi, psip

    def get_LP_constants(self,
            radius: float,
            neff: float,
            nu: int,
            A: list) -> tuple:

        u = self.get_U_W_parameter(
            radius=radius,
            neff=neff,
        )

        if neff < self.refractive_index:
            term_0 = numpy.pi / 2 * (u * yvp(nu, u) * A[0] - yn(nu, u) * A[1])
            term_1 = numpy.pi / 2 * (jn(nu, u) * A[1] - u * jvp(nu, u) * A[0])

        else:
            term_0 = u * kvp(nu, u) * A[0] - kn(nu, u) * A[1]
            term_1 = iv(nu, u) * A[1] - u * ivp(nu, u) * A[0]

        return term_0, term_1

    def EH_fields(self,
            radius_in: float,
            radius_out: float,
            nu: int,
            neff: float,
            EH: list,
            TM: bool = True) -> list:
        """
        Returns EH field component

        :param      radius_in:   The radius in
        :type       radius_in:   float
        :param      radius_out:  The radius out
        :type       radius_out:  float
        :param      nu:          The nu parameter of the mode
        :type       nu:          int
        :param      neff:        The neff
        :type       neff:        float
        :param      EH:          { parameter_description }
        :type       EH:          list
        :param      TM:          { parameter_description }
        :type       TM:          bool

        :returns:   The EH field
        :rtype:     list
        """
        u = self.get_U_W_parameter(radius=radius_out, neff=neff)

        if radius_in == 0:
            if nu == 0:
                if TM:
                    self.C = numpy.array([1., 0., 0., 0.])
                else:
                    self.C = numpy.array([0., 0., 1., 0.])
            else:
                self.C = numpy.zeros((4, 2))
                self.C[0, 0] = 1  # Ez = 1
                self.C[2, 1] = 1  # Hz = alpha

        elif nu == 0:
            self.C = numpy.zeros(4)
            if TM:
                c = numpy.sqrt(epsilon_0 / mu_0) * self.refractive_index**2
                idx = (0, 3)

                self.C[:2] = self.get_TE_TM_constants(
                    radius_in=radius_in,
                    radius_out=radius_out,
                    neff=neff,
                    EH=EH,
                    c=c,
                    idx=idx
                )
            else:
                c = -eta0
                idx = (1, 2)

                self.C[2:] = self.get_TE_TM_constants(
                    radius_in=radius_in,
                    radius_out=radius_out,
                    neff=neff,
                    EH=EH,
                    c=c,
                    idx=idx
                )
        else:
            self.C = self.get_V_constant(
                radius_in=radius_in,
                radius_out=radius_out,
                neff=neff,
                nu=nu,
                EH=EH
            )

        # Compute EH fields
        if neff < self.refractive_index:
            c1 = self.wavelength.k0 * radius_out / u
            F3 = jvp(nu, u) / jn(nu, u)
            F4 = yvp(nu, u) / yn(nu, u)
        else:
            c1 = -self.wavelength.k0 * radius_out / u
            F3 = ivp(nu, u) / iv(nu, u)
            F4 = kvp(nu, u) / kn(nu, u)

        c2 = neff * nu / u * c1
        c3 = eta0 * c1
        c4 = numpy.sqrt(epsilon_0 / mu_0) * self.refractive_index**2 * c1

        EH[0] = self.C[0] + self.C[1]
        EH[1] = self.C[2] + self.C[3]
        EH[2] = (c2 * (self.C[0] + self.C[1]) - c3 * (F3 * self.C[2] + F4 * self.C[3]))
        EH[3] = (c4 * (F3 * self.C[0] + F4 * self.C[1]) - c2 * (self.C[2] + self.C[3]))

        return EH

    def get_V_constant(self,
            radius_in: float,
            radius_out: float,
            neff: float,
            nu,
            EH) -> float:

        a = numpy.zeros((4, 4))

        u = self.get_U_W_parameter(radius=radius_out, neff=neff)

        urp = self.get_U_W_parameter(radius=radius_in, neff=neff)

        if neff < self.refractive_index:
            B1 = jn(nu, u)
            B2 = yn(nu, u)
            F1 = jn(nu, urp) / B1
            F2 = yn(nu, urp) / B2
            F3 = jvp(nu, urp) / B1
            F4 = yvp(nu, urp) / B2
            c1 = self.wavelength.k0 * radius_out / u
        else:
            B1 = iv(nu, u)
            B2 = kn(nu, u)
            F1 = iv(nu, urp) / B1 if u else 1
            F2 = kn(nu, urp) / B2
            F3 = ivp(nu, urp) / B1 if u else 1
            F4 = kvp(nu, urp) / B2
            c1 = -self.wavelength.k0 * radius_out / u

        c2 = neff * nu / urp * c1
        c3 = eta0 * c1
        c4 = numpy.sqrt(epsilon_0 / mu_0) * self.refractive_index**2 * c1

        a[0, 0] = F1
        a[0, 1] = F2
        a[1, 2] = F1
        a[1, 3] = F2
        a[2, 0] = F1 * c2
        a[2, 1] = F2 * c2
        a[2, 2] = -F3 * c3
        a[2, 3] = -F4 * c3
        a[3, 0] = F3 * c4
        a[3, 1] = F4 * c4
        a[3, 2] = -F1 * c2
        a[3, 3] = -F2 * c2

        return numpy.linalg.solve(a, EH)

    def get_TE_TM_constants(self,
            radius_in: float,
            radius_out: float,
            neff: float,
            EH,
            c,
            idx) -> float:

        a = numpy.empty((2, 2))

        u = self.get_U_W_parameter(
            radius=radius_out,
            neff=neff,
        )

        urp = self.get_U_W_parameter(
            radius=radius_in,
            neff=neff,
        )

        if neff < self.refractive_index:
            B1 = j0(u)
            B2 = y0(u)
            F1 = j0(urp) / B1
            F2 = y0(urp) / B2
            F3 = -j1(urp) / B1
            F4 = -y1(urp) / B2
            c1 = self.wavelength.k0 * radius_out / u
        else:
            B1 = i0(u)
            B2 = k0(u)
            F1 = i0(urp) / B1
            F2 = k0(urp) / B2
            F3 = i1(urp) / B1
            F4 = -k1(urp) / B2
            c1 = -self.wavelength.k0 * radius_out / u

        c3 = c * c1

        a[0, 0] = F1
        a[0, 1] = F2
        a[1, 0] = F3 * c3
        a[1, 1] = F4 * c3

        return numpy.linalg.solve(a, EH.take(idx))
