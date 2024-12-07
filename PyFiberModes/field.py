#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from typing import Callable
import scipy
import numpy as np
from PyFiberModes.mode_instances import HE11
from PyFiberModes.mode import Mode
from dataclasses import dataclass
from MPSPlots.styles import mps
from MPSPlots.colormaps import blue_black_red
import matplotlib.pyplot as plt
from PyFiberModes.coordinates import CartesianCoordinates


@dataclass
class CylindricalCoordinates:
    rho: numpy.ndarray
    phi: numpy.ndarray
    z: numpy.ndarray

    def to_cartesian(self) -> object:
        x = self.rho * numpy.cos(self.phi)
        y = self.rho * numpy.sin(self.phi)
        z = self.z

        cartesian_coordinate = CartesianCoordinates(x=x, y=y, z=z)

        return cartesian_coordinate

    def to_cylindrical(self):
        return self


@dataclass
class Field:
    """
    A class representing the electric and magnetic field components in a fiber.

    Parameters
    ----------
    fiber : Fiber
        The fiber associated with the mode.
    mode : Mode
        The mode to evaluate.
    limit : float
        The radius of the field to compute.
    n_point : int, optional
        The number of points for the computation grid (default is 101).

    Attributes
    ----------
    cartesian_coordinates : CartesianCoordinates
        Cartesian coordinates generated for the field computation.
    cylindrical_coordinates : CylindricalCoordinates
        Cylindrical coordinates corresponding to the cartesian grid.
    """
    fiber: object
    mode: Mode
    limit: float
    n_point: int = 101

    def __post_init__(self) -> None:
        """
        Generate the mesh coordinates that are used for field computation.
        """
        self.cartesian_coordinates = CartesianCoordinates.generate_from_square(length=2 * self.limit, n_points=self.n_point)

        self.cylindrical_coordinates = self.cartesian_coordinates.to_cylindrical()

    @staticmethod
    def initialize_array(shape: tuple) -> np.ndarray:
        """Initialize an array of zeros."""
        return np.zeros(shape)

    def get_azimuthal_dependency(self, phi: float, dependency_type: str) -> np.ndarray:
        r"""
        Compute the azimuthal dependency for the field based on the given type.

        The azimuthal dependencies are defined as:

        .. math::
            f_\nu(\phi) = \cos(\nu \phi + \phi_0)

        .. math::
            g_\nu(\phi) = -\sin(\nu \phi + \phi_0)

        where :math:`\nu` is the azimuthal mode number and :math:`\phi_0` is the
        phase shift.

        Parameters
        ----------
        phi : float
            Phase shift of the field in radians.
        dependency_type : {'f', 'g'}
            The type of azimuthal dependency to compute:
            - 'f' for cosine dependency.
            - 'g' for sine dependency.

        Returns
        -------
        numpy.ndarray
            The azimuthal dependency values, computed over the cylindrical coordinates.

        Raises
        ------
        ValueError
            If `dependency_type` is not 'f' or 'g'.

        Notes
        -----
        This method is referenced in Eq. 3.71 of Jaques Bures for the calculation
        of field azimuthal behavior in optical fibers.
        """
        angle = self.mode.nu * self.cylindrical_coordinates.phi + phi
        if dependency_type == "f":
            return np.cos(angle)
        elif dependency_type == "g":
            return -np.sin(angle)
        else:
            raise ValueError("Invalid dependency type. Use 'f' for cosine or 'g' for sine.")

    def get_index_iterator(self, array: np.ndarray):
        """
        Generate an iterator for multi-dimensional array indices.

        This method provides a way to iterate through all indices of a
        multi-dimensional NumPy array in a memory-efficient manner.

        Parameters
        ----------
        array : numpy.ndarray
            The array for which to generate the index iterator.

        Yields
        ------
        tuple
            Multi-dimensional index corresponding to the current element.

        """
        iterator = np.nditer(array, flags=['multi_index'])
        for _ in iterator:
            yield iterator.multi_index

    def wrapper_get_field(function) -> Callable:
        r"""
        Decorator for wrapping field computation methods.

        This decorator can be used to add pre- or post-processing steps
        to field computation methods in the `Field` class.

        Parameters
        ----------
        function : Callable
            The function to be wrapped.

        Returns
        -------
        callable
            The wrapped function.

        """
        def wrapper(self, *args, **kwargs):
            return function(self, *args, **kwargs)
        return wrapper

    def Ex(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the x-component of the electric field, :math:`E_x`.

        The x-component of the electric field is calculated based on the field's mode family.
        For LP (linearly polarized) modes, it uses the radial electric field and the azimuthal
        dependency :math:`f_\nu(\phi)`. For other modes, the x-component is derived using the
        transverse electric field and its polarization angle.

        .. math::
            E_x =
            \begin{cases}
            E_\rho \cdot f_\nu(\phi) & \text{if mode is LP} \\
            E_T \cdot \cos(\theta_{\text{pol}}) & \text{otherwise}
            \end{cases}

        where:
            - :math:`E_\rho` is the radial electric field.
            - :math:`f_\nu(\phi) = \cos(\nu \phi + \phi_0)` is the azimuthal dependency.
            - :math:`E_T` is the transverse electric field.
            - :math:`\theta_{\text{pol}}` is the polarization angle.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            The orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The x-component of the electric field, :math:`E_x`, computed over the Cartesian grid.

        Raises
        ------
        AttributeError
            If `self.fiber` does not implement `get_radial_field`.

        References
        ----------
        Jaques Bures, Optical Fiber Theory, Eq. 3.71.

        """
        if self.mode.family == 'LP':
            # Initialize array for the electric field component
            array = self.initialize_array(self.cartesian_coordinates.x.shape)
            azimuthal_dependency = self.get_azimuthal_dependency(phi=phi, dependency_type='f')

            # Iterate over grid points and calculate E_x
            for index in self.get_index_iterator(array):
                e_field, _ = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )
                array[index] = e_field.rho * azimuthal_dependency[index]
        else:
            # For non-LP modes, calculate using transverse field and polarization
            polarization = self.Epol(phi, theta)
            array = self.Et(phi, theta) * np.cos(polarization)

        return array

    def Ey(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the y-component of the electric field, :math:`E_y`.

        The y-component of the electric field is calculated based on the field's mode family.
        For LP (linearly polarized) modes, it uses the azimuthal dependency :math:`f_\nu(\phi)`
        and the azimuthal electric field. For other modes, it is derived using the transverse
        electric field and its polarization angle.

        .. math::
            E_y =
            \begin{cases}
            E_\phi \cdot f_\nu(\phi) & \text{if mode is LP} \\
            E_T \cdot \sin(\theta_{\text{pol}}) & \text{otherwise}
            \end{cases}

        where:
            - :math:`E_\phi` is the azimuthal electric field.
            - :math:`f_\nu(\phi) = \cos(\nu \phi + \phi_0)` is the azimuthal dependency.
            - :math:`E_T` is the transverse electric field.
            - :math:`\theta_{\text{pol}}` is the polarization angle.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            The orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The y-component of the electric field, :math:`E_y`, computed over the Cartesian grid.

        Raises
        ------
        AttributeError
            If `self.fiber` does not implement `get_radial_field`.

        References
        ----------
        Jaques Bures, Optical Fiber Theory, Eq. 3.71.
        """
        if self.mode.family == 'LP':
            # Initialize array for the electric field component
            array = self.initialize_array(self.cartesian_coordinates.x.shape)
            azimuthal_dependency = self.get_azimuthal_dependency(phi=phi, dependency_type='f')

            # Iterate over grid points and calculate E_y
            for index in self.get_index_iterator(array):
                e_field, _ = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )
                array[index] = e_field.phi * azimuthal_dependency[index]
            return array
        else:
            # For non-LP modes, calculate using transverse field and polarization
            polarization = self.Epol(phi, theta)
            array = self.Et(phi, theta) * np.sin(polarization)

        return array

    def Ez(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the z-component of the electric field, :math:`E_z`.

        The z-component of the electric field is calculated using the azimuthal dependency
        :math:`f_\nu(\phi)` and the longitudinal electric field component :math:`E_z`.

        .. math::
            E_z = E_z \cdot f_\nu(\phi)

        where:
            - :math:`E_z` is the longitudinal electric field.
            - :math:`f_\nu(\phi) = \cos(\nu \phi + \phi_0)` is the azimuthal dependency.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (not used in this calculation, default is 0).

        Returns
        -------
        numpy.ndarray
            The z-component of the electric field, :math:`E_z`, computed over the Cartesian grid.

        Raises
        ------
        AttributeError
            If `self.fiber` does not implement `get_radial_field`.

        References
        ----------
        Jaques Bures, Optical Fiber Theory, Eq. 3.71.
        """
        # Initialize the array for the z-component
        array = self.initialize_array(self.cartesian_coordinates.x.shape)

        # Compute azimuthal dependency
        azimuthal_dependency = self.get_azimuthal_dependency(phi=phi, dependency_type='f')

        # Iterate over grid points and calculate E_z
        for index in self.get_index_iterator(array):
            e_field, _ = self.fiber.get_radial_field(
                mode=self.mode,
                radius=self.cylindrical_coordinates.rho[index]
            )
            array[index] = e_field.z * azimuthal_dependency[index]

        return array

    def Er(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the radial component of the electric field, :math:`E_r`.

        The radial component of the electric field is calculated differently based on the field's
        mode family. For LP (linearly polarized) modes, it uses the transverse electric field and
        its polarization. For other modes, it incorporates the azimuthal dependency
        :math:`f_\nu(\phi)` and the radial electric field component :math:`E_\rho`.

        .. math::
            E_r =
            \begin{cases}
            E_T \cdot \cos(\theta_{\text{pol}} - \phi) & \text{if mode is LP} \\
            E_\rho \cdot f_\nu(\phi) & \text{otherwise}
            \end{cases}

        where:
            - :math:`E_T` is the transverse electric field.
            - :math:`\theta_{\text{pol}}` is the polarization angle.
            - :math:`E_\rho` is the radial electric field.
            - :math:`f_\nu(\phi) = \cos(\nu \phi + \phi_0)` is the azimuthal dependency.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            The orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The radial component of the electric field, :math:`E_r`, computed over the Cartesian grid.

        Raises
        ------
        AttributeError
            If `self.fiber` does not implement `get_radial_field`.

        References
        ----------
        Jaques Bures, Optical Fiber Theory, Eq. 3.71.
        """
        if self.mode.family == 'LP':
            # Compute for LP mode using transverse field and polarization
            polarization = self.Epol(phi, theta) - self.cylindrical_coordinates.phi
            array = self.Et(phi, theta) * np.cos(polarization)
        else:
            # Compute for non-LP mode using azimuthal dependency and radial field
            array = self.initialize_array(self.cartesian_coordinates.x.shape)
            azimuthal_dependency = self.get_azimuthal_dependency(phi=phi, dependency_type='f')

            for index in self.get_index_iterator(array):
                e_field, _ = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )
                array[index] = e_field.rho * azimuthal_dependency[index]

        return array

    def Ephi(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the azimuthal component of the electric field, :math:`E_\phi`.

        The azimuthal component of the electric field is calculated based on the field's mode family.
        For LP (linearly polarized) modes, it uses the transverse electric field and its polarization.
        For other modes, it incorporates the azimuthal dependency :math:`g_\nu(\phi)` and the azimuthal
        electric field component :math:`E_\phi`.

        .. math::
            E_\phi =
            \begin{cases}
            E_T \cdot \sin(\theta_{\text{pol}} - \phi) & \text{if mode is LP} \\
            E_\phi \cdot g_\nu(\phi) & \text{otherwise}
            \end{cases}

        where:
            - :math:`E_T` is the transverse electric field.
            - :math:`\theta_{\text{pol}}` is the polarization angle.
            - :math:`E_\phi` is the azimuthal electric field.
            - :math:`g_\nu(\phi) = -\sin(\nu \phi + \phi_0)` is the azimuthal dependency.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            The orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The azimuthal component of the electric field, :math:`E_\phi`, computed over the Cartesian grid.

        Raises
        ------
        AttributeError
            If `self.fiber` does not implement `get_radial_field`.

        References
        ----------
        Jaques Bures, Optical Fiber Theory, Eq. 3.71.
        """
        if self.mode.family == 'LP':
            # Compute for LP mode using transverse field and polarization
            polarization = self.Epol(phi, theta) - self.cylindrical_coordinates.phi
            array = self.Et(phi, theta) * np.sin(polarization)
        else:
            # Compute for non-LP mode using azimuthal dependency and azimuthal field
            array = self.initialize_array(self.cartesian_coordinates.x.shape)
            azimuthal_dependency = self.get_azimuthal_dependency(phi=phi, dependency_type='g')

            for index in self.get_index_iterator(array):
                e_field, _ = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )
                array[index] = e_field.phi * azimuthal_dependency[index]

        return array

    def Et(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the transverse component of the electric field, :math:`E_T`.

        The transverse electric field is computed as the magnitude of the
        components perpendicular to the z-axis. For LP (linearly polarized)
        modes, it is derived from the x and y components. For other modes, it
        is computed using the radial and azimuthal components.

        .. math::
            E_T =
            \begin{cases}
            \sqrt{E_x^2 + E_y^2} & \text{if mode is LP} \\
            \sqrt{E_r^2 + E_\phi^2} & \text{otherwise}
            \end{cases}

        where:
            - :math:`E_x, E_y` are the Cartesian components of the electric field.
            - :math:`E_r, E_\phi` are the cylindrical components of the electric field.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            The orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The transverse electric field, :math:`E_T`, computed over the Cartesian grid.
        """
        if self.mode.family == 'LP':
            e_x = self.Ex(phi, theta)
            e_y = self.Ey(phi, theta)
            e_transverse = np.sqrt(
                np.square(e_x) + np.square(e_y)
            )
        else:
            e_r = self.Er(phi, theta)
            e_phi = self.Ephi(phi, theta)
            e_transverse = np.sqrt(
                np.square(e_r) + np.square(e_phi)
            )

        return e_transverse

    def Epol(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the polarization angle of the transverse electric field, :math:`\theta_{\text{pol}}`.

        The polarization angle is the direction of the transverse electric field vector
        in the plane perpendicular to the z-axis. For LP (linearly polarized) modes, it is
        calculated using the x and y components. For other modes, it is calculated using
        the radial and azimuthal components.

        .. math::
            \theta_{\text{pol}} =
            \begin{cases}
            \arctan2(E_y, E_x) & \text{if mode is LP} \\
            \arctan2(E_\phi, E_r) + \phi & \text{otherwise}
            \end{cases}

        where:
            - :math:`E_x, E_y` are the Cartesian components of the electric field.
            - :math:`E_r, E_\phi` are the cylindrical components of the electric field.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            The orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The polarization angle of the transverse electric field in radians.
        """
        if self.mode.family == 'LP':
            e_y = self.Ey(phi, theta)
            e_x = self.Ex(phi, theta)
            e_polarization = np.arctan2(e_y, e_x)
        else:
            e_phi = self.Ephi(phi, theta)
            e_r = self.Er(phi, theta)
            e_polarization = np.arctan2(e_phi, e_r) + self.cylindrical_coordinates.phi

        return e_polarization

    def Emod(self, phi: float = 0, theta: float = 0) -> np.ndarray:
        r"""
        Compute the modulus (magnitude) of the electric field vector, :math:`|\vec{E}|`.

        The electric field modulus is the magnitude of the field vector, which
        includes all components (transverse and longitudinal). For LP (linearly
        polarized) modes, it uses the Cartesian components. For other modes, it
        uses the cylindrical components.

        .. math::
            |\vec{E}| =
            \begin{cases}
            \sqrt{E_x^2 + E_y^2 + E_z^2} & \text{if mode is LP} \\
            \sqrt{E_r^2 + E_\phi^2 + E_z^2} & \text{otherwise}
            \end{cases}

        where:
            - :math:`E_x, E_y, E_z` are the Cartesian components of the electric field.
            - :math:`E_r, E_\phi, E_z` are the cylindrical components of the electric field.

        Parameters
        ----------
        phi : float, optional
            The phase of the field in radians (default is 0).
        theta : float, optional
            The orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The modulus of the electric field, :math:`|\vec{E}|`, computed over the Cartesian grid.
        """
        if self.mode.family == 'LP':
            e_x = self.Ex(phi, theta)
            e_y = self.Ey(phi, theta)
            e_z = self.Ez(phi, theta)
            e_modulus = np.sqrt(
                np.square(e_x) + np.square(e_y) + np.square(e_z)
            )
        else:
            e_r = self.Er(phi, theta)
            e_phi = self.Ephi(phi, theta)
            e_z = self.Ez(phi, theta)
            e_modulus = np.sqrt(
                np.square(e_r) + np.square(e_phi) + np.square(e_z)
            )

        return e_modulus

    def Hx(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the x-component of the magnetic field :math:`H_x`.

        For LP (linearly polarized) modes, this is computed using the radial
        component :math:`H_r` and the azimuthal dependency :math:`f_\nu(\phi)`:

        .. math::
            H_x = H_r \cdot f_\nu(\phi)

        For other modes, the transverse magnetic field and polarization are used:

        .. math::
            H_x = H_T \cdot \cos(\theta_{\text{pol}})

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        np.ndarray
            The magnetic field in the x-direction over the Cartesian grid.

        Raises
        ------
        AttributeError
            If the fiber object does not implement `get_radial_field`.

        """
        if self.mode.family == 'LP':
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency(phi=phi, dependency_type='f')

            for index in self.get_index_iterator(array):

                _, h_field = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = h_field.rho * azimuthal_dependency_f[index]

        else:
            polarisation = self.Hpol(phi, theta)
            array = self.Ht(phi, theta) * numpy.cos(polarisation)

        return array

    def Hy(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the y-component of the magnetic field :math:`H_y`.

        For LP modes, this is computed using the azimuthal component :math:`H_\phi`
        and the azimuthal dependency :math:`f_\nu(\phi)`:

        .. math::
            H_y = H_\phi \cdot f_\nu(\phi)

        For other modes, the transverse magnetic field and polarization are used:

        .. math::
            H_y = H_T \cdot \sin(\theta_{\text{pol}})

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        np.ndarray
            The magnetic field in the y-direction over the Cartesian grid.

        Raises
        ------
        AttributeError
            If the fiber object does not implement `get_radial_field`.

        """
        if self.mode.family == 'LP':
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency(phi=phi, dependency_type='f')
            for index in self.get_index_iterator(array):

                _, h_field = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = h_field.phi * azimuthal_dependency_f[index]

        else:
            polarisation = self.Hpol(phi, theta)
            array = self.Ht(phi, theta) * numpy.sin(polarisation)

        return array

    def Hz(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the z-component of the magnetic field :math:`H_z`.

        The z-component is calculated for all modes using the longitudinal field
        :math:`H_z` and the azimuthal dependency :math:`f_\nu(\phi)`:

        .. math::
            H_z = H_z \cdot f_\nu(\phi)

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        np.ndarray
            The magnetic field in the z-direction over the Cartesian grid.

        Raises
        ------
        AttributeError
            If the fiber object does not implement `get_radial_field`.

        """
        array = numpy.zeros(self.cartesian_coordinates.x.shape)
        azimuthal_dependency_f = self.get_azimuthal_dependency(phi=phi, dependency_type='f')
        for index in self.get_index_iterator(array):

            _, h_field = self.fiber.get_radial_field(
                mode=self.mode,
                radius=self.cylindrical_coordinates.rho[index]
            )

            array[index] = h_field.z * azimuthal_dependency_f[index]
        return array

    def Hr(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the radial component of the magnetic field :math:`H_r`.

        For LP modes, this is computed using the transverse field :math:`H_T` and
        the polarization angle :math:`\theta_{\text{pol}}`:

        .. math::
            H_r = H_T \cdot \cos(\theta_{\text{pol}} - \phi)

        For other modes, the radial field and azimuthal dependency :math:`f_\nu(\phi)` are used:

        .. math::
            H_r = H_r \cdot f_\nu(\phi)

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        np.ndarray
            The magnetic field in the radial direction over the Cartesian grid.

        """
        if self.mode.family == 'LP':
            radial = self.Ht(phi, theta)
            polarisation = self.Hpol(phi, theta) - self.cylindrical_coordinates.phi
            azimuthal = numpy.cos(polarisation)
            array = radial * azimuthal

        else:
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency(phi=phi, dependency_type='f')

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = hr[0] * azimuthal_dependency_f[index]

        return array

    def Hphi(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the azimuthal component of the magnetic field :math:`H_\phi`.

        For LP modes, this is computed using the transverse field :math:`H_T` and
        the polarization angle :math:`\theta_{\text{pol}}`:

        .. math::
            H_\phi = H_T \cdot \sin(\theta_{\text{pol}} - \phi)

        For other modes, the azimuthal field and azimuthal dependency :math:`g_\nu(\phi)` are used:

        .. math::
            H_\phi = H_\phi \cdot g_\nu(\phi)

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        np.ndarray
            The magnetic field in the azimuthal direction over the Cartesian grid.

        """
        if self.mode.family == 'LP':
            polarisation = self.Hpol(phi, theta) - self.cylindrical_coordinates.phi
            array = self.Ht(phi, theta) * numpy.sin(polarisation)
        else:
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_g = self.get_azimuthal_dependency(phi=phi, dependency_type='g')

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = hr[1] * azimuthal_dependency_g[index]

        return array

    def Ht(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the transverse component of the magnetic field, :math:`H_T`.

        The transverse magnetic field is calculated as the magnitude of the
        components perpendicular to the z-axis. For LP (linearly polarized) modes,
        it is derived from the x and y components. For other modes, it is computed
        using the radial and azimuthal components.

        .. math::
            H_T =
            \begin{cases}
            \sqrt{H_x^2 + H_y^2} & \text{if mode is LP} \\
            \sqrt{H_r^2 + H_\phi^2} & \text{otherwise}
            \end{cases}

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The transverse magnetic field, :math:`H_T`, computed over the Cartesian grid.

        """
        if self.mode.family == 'LP':
            h_x = self.Hx(phi, theta)
            h_y = self.Hy(phi, theta)
            return numpy.sqrt(numpy.square(h_x) + numpy.square(h_y))
        else:
            h_r = self.Hr(phi, theta)
            h_phi = self.Hphi(phi, theta)
            return numpy.sqrt(numpy.square(h_r) + numpy.square(h_phi))

    def Hpol(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the polarization angle of the transverse magnetic field, :math:`\theta_{\text{pol}}`.

        The polarization angle represents the direction of the transverse magnetic
        field vector in the plane perpendicular to the z-axis. For LP modes, it is
        calculated using the x and y components. For other modes, it is calculated
        using the radial and azimuthal components.

        .. math::
            \theta_{\text{pol}} =
            \begin{cases}
            \arctan2(H_y, H_x) & \text{if mode is LP} \\
            \arctan2(H_\phi, H_r) + \phi & \text{otherwise}
            \end{cases}

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The polarization angle of the transverse magnetic field in radians.

        """
        if self.mode.family == 'LP':
            h_polarization = numpy.arctan2(
                self.Hy(phi, theta),
                self.Hx(phi, theta)
            )

        else:
            h_polarization = numpy.arctan2(
                self.Hphi(phi, theta),
                self.Hr(phi, theta)
            )
            h_polarization += self.cylindrical_coordinates.phi

        return h_polarization

    def Hmod(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Compute the modulus (magnitude) of the magnetic field vector, :math:`|\vec{H}|`.

        The magnetic field modulus is the magnitude of the field vector, which
        includes all components (transverse and longitudinal). For LP (linearly
        polarized) modes, it uses the Cartesian components. For other modes, it
        uses the cylindrical components.

        .. math::
            |\vec{H}| =
            \begin{cases}
            \sqrt{H_x^2 + H_y^2 + H_z^2} & \text{if mode is LP} \\
            \sqrt{H_r^2 + H_\phi^2 + H_z^2} & \text{otherwise}
            \end{cases}

        where:
            - :math:`H_x, H_y, H_z` are the Cartesian components of the magnetic field.
            - :math:`H_r, H_\phi, H_z` are the cylindrical components of the magnetic field.

        Parameters
        ----------
        phi : float, optional
            Phase of the field in radians (default is 0).
        theta : float, optional
            Orientation of the field in radians (default is 0).

        Returns
        -------
        numpy.ndarray
            The modulus of the magnetic field, :math:`|\vec{H}|`, computed over the Cartesian grid.

        """
        if self.mode.family == 'LP':
            h_x = self.Hx(phi, theta)
            h_y = self.Hy(phi, theta)
            h_z = self.Hz(phi, theta)
            h_modulus = numpy.sqrt(
                numpy.square(h_x) + numpy.square(h_y) + numpy.square(h_z)
            )

        else:
            h_r = self.Hr(phi, theta)
            h_phi = self.Hphi(phi, theta)
            h_z = self.Hz(phi, theta)
            h_modulus = numpy.sqrt(
                numpy.square(h_r) + numpy.square(h_phi) + numpy.square(h_z))

        return h_modulus

    def get_effective_area(self) -> float:
        r"""
        Compute the effective area of the mode, :math:`A_{\text{eff}}`.

        The effective area is defined as:

        .. math::
            A_{\text{eff}} = \frac{\left( \int |E|^2 \, dx \, dy \right)^2}{\int |E|^4 \, dx \, dy}

        where:
            - :math:`|E|` is the modulus of the electric field.
            - The integrals are evaluated over the cross-sectional area.

        Returns
        -------
        float
            The effective area, :math:`A_{\text{eff}}`.

        """
        field_array_norm = self.Emod()

        integral = self.get_integrale_square(array=self.Emod())

        term_0 = numpy.square(integral)

        term_1 = numpy.sum(numpy.power(field_array_norm, 4))

        return (term_0 / term_1)

    def get_integrale_square(self, array: numpy.ndarray) -> float:
        r"""
        Compute the square of the integral of the given array.

        .. math::
            \left( \int |F(x, y)|^2 \, dx \, dy \right)

        Parameters
        ----------
        array : numpy.ndarray
            Array representing the field values over the grid.

        Returns
        -------
        float
            The squared integral of the array.

        """
        square_field = numpy.square(array)

        sum_square_field = numpy.sum(square_field)

        integral = sum_square_field * self.cartesian_coordinates.dx * self.cartesian_coordinates.dy

        return integral

    def get_intensity(self) -> float:
        r"""
        Compute the intensity of the mode, :math:`I`.

        The intensity is defined as:

        .. math::
            I = \frac{n_{\text{eff}}}{n_{\text{eff}}^{\text{HE11}}} \cdot \int |E_T|^2 \, dx \, dy

        where:
            - :math:`n_{\text{eff}}` is the effective index of the mode.
            - :math:`n_{\text{eff}}^{\text{HE11}}` is the effective index of the fundamental HE11 mode.
            - :math:`|E_T|^2` is the squared transverse electric field.

        Returns
        -------
        float
            The intensity of the mode.

        """
        HE11_n_eff = self.fiber.get_effective_index(
            mode=HE11,
            wavelength=self.fiber.wavelength
        )

        n_eff = self.fiber.get_effective_index(
            mode=self.mode,
            wavelength=self.fiber.wavelength
        )

        norm_squared = self.get_integrale_square(array=self.Et())

        return n_eff / HE11_n_eff * norm_squared

    def get_normalization_constant(self) -> float:
        r"""
        Compute the normalization constant, :math:`N`.

        The normalization constant is defined as:

        .. math::
            N = \frac{I}{2} \cdot \epsilon_0 \cdot n_{\text{eff}}^{\text{HE11}} \cdot c

        where:
            - :math:`I` is the intensity of the mode.
            - :math:`\epsilon_0` is the permittivity of free space.
            - :math:`n_{\text{eff}}^{\text{HE11}}` is the effective index of the fundamental HE11 mode.
            - :math:`c` is the speed of light in vacuum.

        Returns
        -------
        float
            The normalization constant, :math:`N`.

        """
        neff = self.fiber.neff(mode=HE11, wavelength=self.fiber.wavelength)

        intensity = self.get_intensity()

        return 0.5 * scipy.constants.epsilon_0 * neff * scipy.constants.c * intensity

    def get_poynting_vector(self):
        """
        Gets the Poynting vector but is not implemented yet.
        """
        raise NotImplementedError('Not yet implemented')

    def plot(self, plot_type: list = (), show: bool = True, save_filename: str = None) -> plt.Figure:
        """
        Plotting function.

        Parameters
        ----------
        ax : plt.Axes, optional
            A matplotlib Axes object to draw the plot on. If None, a new figure and axes are created.
            Default is None.
        show : bool, optional
            Whether to display the plot. If False, the plot will not be shown but can still be saved
            or returned. Default is True.
        save_filename : str, optional
            A file path to save the figure. If None, the figure will not be saved. Default is None.

        Returns
        -------
        plt.Figure
            The matplotlib Figure object created or used for the plot.
        """
        with plt.style.context(mps):
            figure, axes = plt.subplots(1, len(plot_type))

        for ax, field_string in zip(axes, plot_type):
            ax.set_aspect('equal')
            field = getattr(self, field_string)()
            max_abs = abs(max(field.max(), field.min()))
            ax.pcolormesh(field, vmin=-max_abs, vmax=max_abs, cmap=blue_black_red)

        figure.tight_layout()

        if save_filename:
            figure.savefig(save_filename)

        if show:
            plt.show()

        return figure

# -
