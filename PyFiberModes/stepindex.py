#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import numpy as np
from dataclasses import dataclass

from scipy.special import jn, yn, iv, kn
from scipy.special import j0, y0, i0, k0
from scipy.special import j1, y1, i1, k1
from scipy.special import jvp, yvp, ivp, kvp
from scipy.constants import mu_0, epsilon_0, physical_constants

# Physical constants
ETA_0 = physical_constants['characteristic impedance of vacuum'][0]


@dataclass
class Geometry:
    """
    Represents a geometric structure with refractive indices and radii.

    Parameters
    ----------
    radius_in : float
        Minimum radius of the structure.
    radius_out : float
        Maximum radius of the structure.
    index_list : list
        List of refractive indices for the structure.
    """

    radius_in: float
    radius_out: float
    index_list: list

    def __post_init__(self):
        """
        Post-initialization of the class.
        Computes additional attributes such as refractive index and thickness.
        """
        self.refractive_index = self.index_list[0]
        self.thickness = self.radius_out - self.radius_in

    def __hash__(self):
        """
        Computes a unique hash for the object based on its attributes.

        Returns
        -------
        int
            Hash value for the object.
        """
        return hash((self.radius_in, self.radius_out, tuple(self.index_list)))


class StepIndex(Geometry):
    """
    Step-index structure for optical fibers and waveguides.

    Inherits
    --------
    Geometry
    """

    def get_index_at_radius(self, radius: float) -> float:
        """
        Get the refractive index at a specific radius.

        Parameters
        ----------
        radius : float
            Radius for evaluation.

        Returns
        -------
        float or None
            Refractive index at the given radius if within bounds, else None.
        """
        return self.refractive_index if self.radius_in <= abs(radius) <= self.radius_out else None

    def get_U_W_parameter(self, radius: float, neff: float) -> float:
        r"""
        Calculate the U or W parameter for waveguides.

        The U and W parameters are defined as:

        .. math::
            U = k_0 \rho \sqrt{n_{\text{core}}^2 - n_{\text{eff}}^2}

            W = k_0 \rho \sqrt{n_{\text{eff}}^2 - n_{\text{clad}}^2}

        where:
        - \(k_0\) is the free-space wave number (\(2 \pi / \lambda\)).
        - \(n_{\text{core}}\), \(n_{\text{clad}}\) are the refractive indices.
        - \(n_{\text{eff}}\) is the effective refractive index.

        Parameters
        ----------
        radius : float
            Radius for evaluation.
        neff : float
            Effective refractive index.

        Returns
        -------
        float
            U or W parameter value.
        """
        index = self.get_index_at_radius(radius)
        if index is None:
            return 0

        return (2 * numpy.pi / self.wavelength) * radius * np.sqrt(abs(index**2 - neff**2))

    def get_psi(self, radius: float, neff: float, nu: int, C: list) -> tuple:
        r"""
        Compute the :math:`\psi` function and its derivative.

        The function :math:`\psi` is computed using Bessel functions for core layers
        and modified Bessel functions for cladding layers. The equations used are:

        Core Layer:

        .. math::
            \psi = C_0 J_\nu(U) + C_1 Y_\nu(U)

            \dot{\psi} = C_0 U J'_\nu(U) + C_1 U Y'_\nu(U)

        Cladding Layer:

        .. math::
            \psi = C_0 I_\nu(U) + C_1 K_\nu(U)

            \dot{\psi} = C_0 U I'_\nu(U) + C_1 U K'_\nu(U)

        where:
            - J$_\nu$, Y$_\nu$: Bessel functions of the first and second kind.
            - I$_\nu$, K$_\nu$: Modified Bessel functions of the first and second kind.
            - C$_0$, C$_1$: Coefficients.
            - U: Parameter related to radius and refractive indices.

        Parameters
        ----------
        radius : float
            The radius at which the field is evaluated.
        neff : float
            The effective refractive index.
        nu : int
            The order of the Bessel function (mode number).
        C : list
            Coefficients for the Bessel functions. Must have at least one value,
            with the second value optional (defaults to 0).

        Returns
        -------
        tuple
            A tuple containing :math:`\psi`: The field value and :math:`\dot{\psi}`: The derivative of the field value.

        """
        u = self.get_U_W_parameter(radius=radius, neff=neff)
        layer_max_index = self.refractive_index

        # Determine the Bessel function type based on the refractive index
        if neff < layer_max_index:
            bessel_func = (jn, yn, jvp, yvp)
        else:
            bessel_func = (iv, kn, ivp, kvp)

        b1, b2, db1, db2 = bessel_func

        # Compute psi and its derivative
        C1 = C[1] if len(C) > 1 else 0  # Handle optional second coefficient
        psi = C[0] * b1(nu, u) + C1 * b2(nu, u)
        psip = C[0] * u * db1(nu, u) + C1 * db2(nu, u)

        return psi, psip

    def get_LP_constants(self, radius: float, neff: float, nu: int, A: list) -> tuple:
        r"""
        Calculate the LP mode constants.

        The constants are derived from Bessel or modified Bessel functions
        for specific mode conditions.

        Parameters
        ----------
        radius : float
            Radius for evaluation.
        neff : float
            Effective refractive index.
        nu : int
            Mode number.
        A : list
            Coefficients for Bessel functions.

        Returns
        -------
        tuple
            LP mode constants.
        """
        u = self.get_U_W_parameter(radius, neff)
        if neff < self.refractive_index:
            term_0 = np.pi / 2 * (u * yvp(nu, u) * A[0] - yn(nu, u) * A[1])
            term_1 = np.pi / 2 * (jn(nu, u) * A[1] - u * jvp(nu, u) * A[0])
        else:
            term_0 = u * kvp(nu, u) * A[0] - kn(nu, u) * A[1]
            term_1 = iv(nu, u) * A[1] - u * ivp(nu, u) * A[0]

        return term_0, term_1

    def EH_fields(self, radius_in: float, radius_out: float, nu: int, neff: float, EH: list, TM: bool = True) -> list:
        r"""
        Compute the EH field components.

        The EH field components are derived based on the mode parameters and
        the coefficients \(C_i\).

        Parameters
        ----------
        radius_in : float
            Inner radius of the structure.
        radius_out : float
            Outer radius of the structure.
        nu : int
            Mode number.
        neff : float
            Effective refractive index.
        EH : list
            List of field components.
        TM : bool, optional
            Indicates if TM mode (default is True).

        Returns
        -------
        list
            Updated EH field components.
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
                c = -ETA_0
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
            c1 = (2 * numpy.pi / self.wavelength) * radius_out / u
            F3 = jvp(nu, u) / jn(nu, u)
            F4 = yvp(nu, u) / yn(nu, u)
        else:
            c1 = -(2 * numpy.pi / self.wavelength) * radius_out / u
            F3 = ivp(nu, u) / iv(nu, u)
            F4 = kvp(nu, u) / kn(nu, u)

        c2 = neff * nu / u * c1
        c3 = ETA_0 * c1
        c4 = numpy.sqrt(epsilon_0 / mu_0) * self.refractive_index**2 * c1

        EH[0] = self.C[0] + self.C[1]
        EH[1] = self.C[2] + self.C[3]
        EH[2] = (c2 * (self.C[0] + self.C[1]) - c3 * (F3 * self.C[2] + F4 * self.C[3]))
        EH[3] = (c4 * (F3 * self.C[0] + F4 * self.C[1]) - c2 * (self.C[2] + self.C[3]))

        return EH

    def get_V_constant(self, radius_in: float, radius_out: float, neff: float, nu: int, EH: np.ndarray) -> np.ndarray:
        """
        Compute the V constants for the mode.

        This method calculates the constants using Bessel or modified Bessel
        functions based on the refractive index regime.

        Parameters
        ----------
        radius_in : float
            Inner radius of the structure.
        radius_out : float
            Outer radius of the structure.
        neff : float
            Effective refractive index.
        nu : int
            Mode number.
        EH : np.ndarray
            Electric and magnetic field components.

        Returns
        -------
        np.ndarray
            Solution to the linear system representing the V constants.
        """
        # Initialize the coefficient matrix
        a = np.zeros((4, 4))

        # Calculate U parameters for inner and outer radii
        u_out = self.get_U_W_parameter(radius=radius_out, neff=neff)
        u_in = self.get_U_W_parameter(radius=radius_in, neff=neff)

        # Determine the Bessel function types based on refractive index
        if neff < self.refractive_index:
            # Core layer: Standard Bessel functions
            bessel_func = (jn, yn, jvp, yvp)
            sign = 1
        else:
            # Cladding layer: Modified Bessel functions
            bessel_func = (iv, kn, ivp, kvp)
            sign = -1

        B1 = bessel_func[0](nu, u_out)
        B2 = bessel_func[1](nu, u_out)
        F1 = bessel_func[0](nu, u_in) / B1 if B1 != 0 else 1
        F2 = bessel_func[1](nu, u_in) / B2 if B2 != 0 else 1
        F3 = bessel_func[2](nu, u_in) / B1 if B1 != 0 else 1
        F4 = bessel_func[3](nu, u_in) / B2 if B2 != 0 else 1

        # Compute scaling coefficients
        c1 = sign * (2 * numpy.pi / self.wavelength) * radius_out / u_out if u_out != 0 else 1
        c2 = neff * nu / u_in * c1 if u_in != 0 else 1
        c3 = ETA_0 * c1
        c4 = np.sqrt(epsilon_0 / mu_0) * self.refractive_index**2 * c1

        # Fill the coefficient matrix
        a[0, 0], a[0, 1] = F1, F2
        a[1, 2], a[1, 3] = F1, F2
        a[2, 0], a[2, 1] = F1 * c2, F2 * c2
        a[2, 2], a[2, 3] = -F3 * c3, -F4 * c3
        a[3, 0], a[3, 1] = F3 * c4, F4 * c4
        a[3, 2], a[3, 3] = -F1 * c2, -F2 * c2

        # Solve the linear system for V constants
        return np.linalg.solve(a, EH)

    def get_TE_TM_constants(
            self,
            radius_in: float,
            radius_out: float,
            neff: float,
            EH: np.ndarray,
            c: float,
            idx: tuple
            ) -> np.ndarray:
        """
        Compute the TE/TM constants for the waveguide.

        This method calculates the constants using standard or modified Bessel
        functions depending on the refractive index regime.

        Parameters
        ----------
        radius_in : float
            Inner radius of the structure.
        radius_out : float
            Outer radius of the structure.
        neff : float
            Effective refractive index.
        EH : np.ndarray
            Electric and magnetic field components.
        c : float
            Scaling constant for the calculation.
        idx : tuple
            Indices for selecting specific components of the EH array.

        Returns
        -------
        np.ndarray
            Solution to the linear system representing the TE/TM constants.
        """
        # Initialize the coefficient matrix
        a = np.empty((2, 2))

        # Calculate U parameters for inner and outer radii
        u_out = self.get_U_W_parameter(radius=radius_out, neff=neff)
        u_in = self.get_U_W_parameter(radius=radius_in, neff=neff)

        # Determine the Bessel function types based on refractive index
        if neff < self.refractive_index:
            # Core layer: Standard Bessel functions
            bessel_func = (j0, y0, j1, y1)
            sign = 1
        else:
            # Cladding layer: Modified Bessel functions
            bessel_func = (i0, k0, i1, k1)
            sign = -1

        B1 = bessel_func[0](u_out)
        B2 = bessel_func[1](u_out)

        # Prevent division by zero for B1 and B2
        F1 = bessel_func[0](u_in) / B1 if B1 != 0 else 1
        F2 = bessel_func[1](u_in) / B2 if B2 != 0 else 1
        F3 = bessel_func[2](u_in) / B1 if B1 != 0 else 1
        F4 = bessel_func[3](u_in) / B2 if B2 != 0 else 1

        # Compute scaling coefficient
        c1 = sign * (2 * numpy.pi / self.wavelength) * radius_out / u_out if u_out != 0 else 1
        c3 = c * c1

        # Fill the coefficient matrix
        a[0, 0], a[0, 1] = F1, F2
        a[1, 0], a[1, 1] = F3 * c3, F4 * c3

        # Solve the linear system for TE/TM constants
        return np.linalg.solve(a, EH.take(idx))
