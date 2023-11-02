import numpy
import scipy
from PyFiberModes import Wavelength, ModeFamily, HE11
from PyFiberModes.mode import Mode
from dataclasses import dataclass


@dataclass
class Field(object):
    fiber: 'Fiber'
    """ Fiber associated to the mode """
    mode: Mode
    """ Mode to evaluate """
    wavelength: float
    """ Wavelength for the simulation """
    limit: float
    """ Radius of the field to compute. """
    n_point: int = 101
    """ Number of points (field will be np x np) """

    FTYPES = ('Ex', 'Ey', 'Ez', 'Er', 'Ephi', 'Et', 'Epol', 'Emod',
              'Hx', 'Hy', 'Hz', 'Hr', 'Hphi', 'Ht', 'Hpol', 'Hmod')

    def __post_init__(self):
        self.wavelength = Wavelength(self.wavelength)
        self.dx = 2 * self.limit / (self.n_point - 1)
        self.dy = 2 * self.limit / (self.n_point - 1)

        self.x_mesh, self.y_mesh = numpy.mgrid[
            -self.limit: self.limit: complex(self.n_point),
            -self.limit: self.limit: complex(self.n_point)
        ]

        self.r_mesh = numpy.sqrt(numpy.square(self.x_mesh) + numpy.square(self.y_mesh))

        self.phi_mesh = numpy.arctan2(self.y_mesh, self.x_mesh)

    def get_azimuthal_dependency_f(self, phi: float) -> numpy.ndarray:
        """
        Gets the azimuthal dependency f.

        :param      phi:  Phase (rotation) of the field.
        :type       phi:  float

        :returns:   The azimuthal dependency g values in [-1, 1].
        :rtype:     { return_type_description }
        """
        return numpy.cos(self.mode.nu * self.phi_mesh + phi)

    def get_azimuthal_dependency_g(self, phi: float) -> numpy.ndarray:
        """
        Gets the azimuthal dependency g.

        :param      phi:  Phase (rotation) of the field.
        :type       phi:  float

        :returns:   The azimuthal dependency g values in [-1, 1].
        :rtype:     { return_type_description }
        """
        return -numpy.sin(self.mode.nu * self.phi_mesh + phi)

    def get_index_iterator(self, array: numpy.ndarray) -> tuple:
        iterator = numpy.nditer(array, flags=['multi_index'])
        for _ in iterator:
            yield iterator.multi_index

    def wrapper_get_field(function):
        def wrapper(self, *args, **kwargs):
            return function(self, *args, **kwargs)

        return wrapper

    def Ex(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        X component of the electric field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the x-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):
                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = er[0] * azimuthal_dependency[index]

        else:
            polarisation = self.Epol(phi, theta)
            array = self.Et(phi, theta) * numpy.cos(polarisation)

        return array

    def Ey(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Y component of the electric field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the y-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):
                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = er[1] * azimuthal_dependency[index]
            return array
        else:
            polarisation = self.Epol(phi, theta)
            array = self.Et(phi, theta) * numpy.sin(polarisation)

        return array

    def Ez(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Z component of the electric field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the z-direction
        :rtype:     numpy.ndarray
        """
        array = numpy.zeros(self.x_mesh.shape)

        azimuthal_dependency = self.get_azimuthal_dependency_f(phi=phi)

        for index in self.get_index_iterator(array):
            er, hr = self.fiber.get_radial_field(
                mode=self.mode,
                wavelength=self.wavelength,
                radius=self.r_mesh[index]
            )

            array[index] = er[2] * azimuthal_dependency[index]

        return array

    def Er(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Radial component of the electric field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the r-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            polarisation = self.Epol(phi, theta) - self.phi_mesh
            array = self.Et(phi, theta) * numpy.cos(polarisation)

        else:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):
                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = er[0] * azimuthal_dependency_f[index]

        return array

    def Ephi(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Phi component of the electric field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the phi-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            polarisation = self.Epol(phi, theta) - self.phi_mesh
            array = self.Et(phi, theta) * numpy.sin(polarisation)

        else:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_g = self.get_azimuthal_dependency_g(phi=phi)

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = er[1] * azimuthal_dependency_g[index]

        return array

    def Et(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Transverse component of the electric field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the transverse-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
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
        """
        Polarization of the transverse electric field (in radians).

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The polarisation of the transverse field
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            e_y = self.Ey(phi, theta)
            e_x = self.Ex(phi, theta)
            e_polarization = numpy.arctan2(e_y, e_x)
        else:
            e_phi = self.Ephi(phi, theta)
            e_r = self.Er(phi, theta)
            e_polarization = numpy.arctan2(e_phi, e_r) + self.phi_mesh

        return e_polarization

    def Emod(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """modulus of the E field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
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
        """
        X component of the magnetic field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the x-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = hr[0] * azimuthal_dependency_f[index]

        else:
            polarisation = self.Hpol(phi, theta)
            array = self.Ht(phi, theta) * numpy.cos(polarisation)

        return array

    def Hy(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Y component of the magnetic field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the y-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)
            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = hr[1] * azimuthal_dependency_f[index]

        else:
            polarisation = self.Hpol(phi, theta)
            array = self.Ht(phi, theta) * numpy.sin(polarisation)

        return array

    def Hz(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Z component of the magnetic field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the z-direction
        :rtype:     numpy.ndarray
        """
        array = numpy.zeros(self.x_mesh.shape)
        azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)
        for index in self.get_index_iterator(array):

            er, hr = self.fiber.get_radial_field(
                mode=self.mode,
                wavelength=self.wavelength,
                radius=self.r_mesh[index]
            )

            array[index] = hr[2] * azimuthal_dependency_f[index]
        return array

    def Hr(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Radial component of the magnetic field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the radial-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            radial = self.Ht(phi, theta)
            polarisation = self.Hpol(phi, theta) - self.phi_mesh
            azimuthal = numpy.cos(polarisation)
            array = radial * azimuthal

        else:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = hr[0] * azimuthal_dependency_f[index]

        return array

    def Hphi(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Azimuthal component of the magnetic field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the phi-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            polarisation = self.Hpol(phi, theta) - self.phi_mesh
            array = self.Ht(phi, theta) * numpy.sin(polarisation)
        else:
            array = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_g = self.get_azimuthal_dependency_g(phi=phi)

            for index in self.get_index_iterator(array):

                er, hr = self.fiber.get_radial_field(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                array[index] = hr[1] * azimuthal_dependency_g[index]

        return array

    def Ht(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Transverse component of the magnetic field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The magnetic field in the transverse-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
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
        if self.mode.family is ModeFamily.LP:
            h_polarization = numpy.arctan2(
                self.Hy(phi, theta),
                self.Hx(phi, theta)
            )

        else:
            h_polarization = numpy.arctan2(
                self.Hphi(phi, theta),
                self.Hr(phi, theta)
            )
            h_polarization += self.phi_mesh

        return h_polarization

    def Hmod(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """modulus of the H field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
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
        Estimation of mode effective area.

        Suppose than r is large enough, such as F(r, r) = 0.

        """
        field_array_norm = self.Emod()

        term_0 = numpy.square(numpy.sum(numpy.square(field_array_norm)))

        term_1 = numpy.sum(numpy.power(field_array_norm, 4))

        return (term_0 / term_1) * self.dx * self.dy

    def get_intensity(self) -> float:
        """
        Gets the intensity.

        :returns:   The intensity.
        :rtype:     float
        """
        HE11_n_eff = self.fiber.neff(
            mode=HE11,
            wavelength=self.wavelength
        )

        n_eff = self.fiber.neff(
            mode=self.mode,
            wavelength=self.wavelength
        )

        norm_squared = numpy.sum(numpy.square(self.Et()))

        return n_eff / HE11_n_eff * norm_squared * self.dx * self.dy

    def get_normalization_constant(self) -> float:
        """
        Gets the normalization constant.

        :returns:   The normalization constant.
        :rtype:     float
        """
        neff = self.fiber.neff(
            mode=HE11,
            wavelength=self.wavelength
        )

        return 0.5 * scipy.constants.epsilon_0 * neff * scipy.constants.c * self.get_intensity()

    def S(self):
        """
        Gets the Poyinting vector.

        :returns:   The Poyinting vector modulus.
        :rtype:     float
        """
        raise NotImplementedError('Not yet implemented')

# -
