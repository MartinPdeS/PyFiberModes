#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import scipy
from PyFiberModes.mode_instances import HE11
from PyFiberModes.mode import Mode
from dataclasses import dataclass
from MPSPlots.render2D import SceneList, Axis
from MPSTools.tools.coordinates import CartesianCoordinates


@dataclass
class Field(object):
    fiber: 'Fiber'
    """ Fiber associated to the mode """
    mode: Mode
    """ Mode to evaluate """
    limit: float
    """ Radius of the field to compute. """
    n_point: int = 101
    """ Number of points (field will be n_point x n_point) """

    def __post_init__(self) -> None:
        """
        Generate the mesh coordinates that are used for field computation.
        """
        self.cartesian_coordinates = CartesianCoordinates.generate_from_square(size=2 * self.limit, n_points=self.n_point)

        self.cylindrical_coordinates = self.cartesian_coordinates.to_cylindrical()

    def get_azimuthal_dependency_f(self, phi: float) -> numpy.ndarray:
        r"""
        Gets the azimuthal dependency :math:`f_\nu`.
        Reference: Eq. 3.71 of Jaques Bures.

        :param      phi:  Phase (rotation) of the field.
        :type       phi:  float

        :returns:   The azimuthal dependency g values in [-1, 1].
        :rtype:     numpy.ndarray
        """
        return numpy.cos(self.mode.nu * self.cylindrical_coordinates.phi + phi)

    def get_azimuthal_dependency_g(self, phi: float) -> numpy.ndarray:
        r"""
        Gets the azimuthal dependency :math:`g_\nu`.
        Reference: Eq. 3.71 of Jaques Bures.

        :param      phi:  Phase (rotation) of the field.
        :type       phi:  float

        :returns:   The azimuthal dependency g values in [-1, 1].
        :rtype:     numpy.ndarray
        """
        return -numpy.sin(self.mode.nu * self.cylindrical_coordinates.phi + phi)

    def get_index_iterator(self, array: numpy.ndarray) -> tuple:
        iterator = numpy.nditer(array, flags=['multi_index'])
        for _ in iterator:
            yield iterator.multi_index

    def wrapper_get_field(function):
        def wrapper(self, *args, **kwargs):
            return function(self, *args, **kwargs)

        return wrapper

    def Ex(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        X component of the electric field :math:`E_{x}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the x-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):
                e_field, _ = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = e_field.rho * azimuthal_dependency[index]

        else:
            polarisation = self.Epol(phi, theta)
            array = self.Et(phi, theta) * numpy.cos(polarisation)

        return array

    def Ey(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Y component of the electric field :math:`E_{y}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the y-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):
                e_field, _ = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = e_field.phi * azimuthal_dependency[index]
            return array
        else:
            polarisation = self.Epol(phi, theta)
            array = self.Et(phi, theta) * numpy.sin(polarisation)

        return array

    def Ez(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Z component of the electric field :math:`E_{z}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the z-direction
        :rtype:     numpy.ndarray
        """
        array = numpy.zeros(self.cartesian_coordinates.x.shape)

        azimuthal_dependency = self.get_azimuthal_dependency_f(phi=phi)

        for index in self.get_index_iterator(array):
            e_field, _ = self.fiber.get_radial_field(
                mode=self.mode,
                radius=self.cylindrical_coordinates.rho[index]
            )

            array[index] = e_field.z * azimuthal_dependency[index]

        return array

    def Er(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Radial component of the electric field :math:`E_{r}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the r-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            polarisation = self.Epol(phi, theta) - self.cylindrical_coordinates.phi
            array = self.Et(phi, theta) * numpy.cos(polarisation)

        else:
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):
                e_field, _ = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )
                array[index] = e_field.rho * azimuthal_dependency_f[index]

        return array

    def Ephi(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Phi component of the electric field :math:`E_{\phi}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the phi-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            polarisation = self.Epol(phi, theta) - self.cylindrical_coordinates.phi
            array = self.Et(phi, theta) * numpy.sin(polarisation)

        else:
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_g = self.get_azimuthal_dependency_g(phi=phi)

            for index in self.get_index_iterator(array):
                e_field, h_field = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = e_field.phi * azimuthal_dependency_g[index]

        return array

    def Et(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Transverse component of the electric field :math:`E_{T}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the transverse-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            e_x = self.Ex(phi, theta)
            e_y = self.Ey(phi, theta)
            e_transverse = numpy.sqrt(
                numpy.square(e_x) + numpy.square(e_y)
            )
        else:
            e_r = self.Er(phi, theta)
            e_phi = self.Ephi(phi, theta)
            e_transverse = numpy.sqrt(
                numpy.square(e_r) + numpy.square(e_phi)
            )

        return e_transverse

    def Epol(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Polarization of the transverse electric field (in radians).

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The polarisation of the transverse field
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            e_y = self.Ey(phi, theta)
            e_x = self.Ex(phi, theta)
            e_polarization = numpy.arctan2(e_y, e_x)
        else:
            e_phi = self.Ephi(phi, theta)
            e_r = self.Er(phi, theta)
            e_polarization = numpy.arctan2(e_phi, e_r) + self.cylindrical_coordinates.phi

        return e_polarization

    def Emod(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Norm of the E vector. :math:`\vec{E}`.

        :param      phi:                  The phi
        :type       phi:                  float
        :param      theta:                The theta
        :type       theta:                float

        :returns:   Norm of the H vector
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            e_x = self.Ex(phi, theta)
            e_y = self.Ey(phi, theta)
            e_z = self.Ez(phi, theta)
            e_modulus = numpy.sqrt(
                numpy.square(e_x) + numpy.square(e_y) + numpy.square(e_z)
            )
        else:
            e_r = self.Er(phi, theta)
            e_phi = self.Ephi(phi, theta)
            e_z = self.Ez(phi, theta)
            e_modulus = numpy.sqrt(
                numpy.square(e_r) + numpy.square(e_phi) + numpy.square(e_z)
            )

        return e_modulus

    def Hx(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        X component of the magnetic field :math:`H_{x}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the x-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)

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
        Y component of the magnetic field :math:`H_{y}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the y-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)
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
        Z component of the magnetic field :math:`H_{z}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the z-direction
        :rtype:     numpy.ndarray
        """
        array = numpy.zeros(self.cartesian_coordinates.x.shape)
        azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)
        for index in self.get_index_iterator(array):

            _, h_field = self.fiber.get_radial_field(
                mode=self.mode,
                radius=self.cylindrical_coordinates.rho[index]
            )

            array[index] = h_field.z * azimuthal_dependency_f[index]
        return array

    def Hr(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Radial component of the magnetic field :math:`H_{r}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the radial-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            radial = self.Ht(phi, theta)
            polarisation = self.Hpol(phi, theta) - self.cylindrical_coordinates.phi
            azimuthal = numpy.cos(polarisation)
            array = radial * azimuthal

        else:
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = hr[0] * azimuthal_dependency_f[index]

        return array

    def Hphi(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Azimuthal component of the magnetic field :math:`H_{\phi}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the phi-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family == 'LP':
            polarisation = self.Hpol(phi, theta) - self.cylindrical_coordinates.phi
            array = self.Ht(phi, theta) * numpy.sin(polarisation)
        else:
            array = numpy.zeros(self.cartesian_coordinates.x.shape)
            azimuthal_dependency_g = self.get_azimuthal_dependency_g(phi=phi)

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    radius=self.cylindrical_coordinates.rho[index]
                )

                array[index] = hr[1] * azimuthal_dependency_g[index]

        return array

    def Ht(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        r"""
        Transverse component of the magnetic field :math:`H_{T}`.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the transverse-direction
        :rtype:     numpy.ndarray
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
        """
        Polarization of the transverse magnetic field (in radians).

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The polarisation of the transverse magnetic field
        :rtype:     numpy.ndarray
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
        Norm of the H vector :math:`\vec{H}`.

        :param      phi:                  The phi
        :type       phi:                  float
        :param      theta:                The theta
        :type       theta:                float

        :returns:   Norm of the H vector
        :rtype:     numpy.ndarray
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
        """
        Gets of mode effective area. Suppose than r is large enough, such as F(r, r) = 0.

        :returns:   The effective area.
        :rtype:     float
        """
        field_array_norm = self.Emod()

        integral = self.get_integrale_square(array=self.Emod())

        term_0 = numpy.square(integral)

        term_1 = numpy.sum(numpy.power(field_array_norm, 4))

        return (term_0 / term_1)

    def get_integrale_square(self, array: numpy.ndarray) -> float:
        r"""
        Gets the integrale of the array squared.

        :param      array: The array
        :type       array: numpy.ndarray

        :returns:   The integrale square.
        :rtype:     float
        """
        square_field = numpy.square(array)

        sum_square_field = numpy.sum(square_field)

        integral = sum_square_field * self.cartesian_coordinates.dx * self.cartesian_coordinates.dy

        return integral

    def get_intensity(self) -> float:
        r"""
        Gets the intensity value of the mode.

        .. math::
            I = \frac{n_{eff}}{n_{eff}^{HE11}} * \int E_T^2 * dx * dy

        :returns:   The intensity.
        :rtype:     float
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
        Gets the normalization constant :math:`N`.

        .. math::
            N = \frac{I}{2} * \epsilon_0 * n_{eff}^{HE11} * c

        :returns:   The normalization constant.
        :rtype:     float
        """
        neff = self.fiber.neff(
            mode=HE11,
            wavelength=self.fiber.wavelength
        )

        intensity = self.get_intensity()

        return 0.5 * scipy.constants.epsilon_0 * neff * scipy.constants.c * intensity

    def get_poynting_vector(self):
        """
        Gets the Poynting vector but is not implemented yet.

        :returns:   The Poynting vector modulus.
        :rtype:     float
        """
        raise NotImplementedError('Not yet implemented')

    def add_to_ax(self, field_string: str, ax: Axis):
        field = getattr(self, field_string)()

        image = ax.add_mesh(scalar=field)

        ax.add_colorbar(symmetric=True, width='10%', artist=image)

    def plot(self, plot_type: list = []) -> SceneList:
        figure = SceneList(unit_size=(4, 4), ax_orientation='horizontal', tight_layout=True)

        for field_string in plot_type:
            ax = figure.append_ax(title=field_string, equal=True)

            self.add_to_ax(field_string=field_string, ax=ax)

        return figure

# -
