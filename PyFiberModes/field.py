import numpy
import scipy
from itertools import product
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
        Azimuthal dependency function.

        Args:
            phi0(float): Phase (rotation) of the field.

        Returns:
            2D array of values (ndarray). Values are between -1 and 1.

        """
        return numpy.cos(self.mode.nu * self.phi_mesh + phi)

    def get_azimuthal_dependency_g(self, phi: float) -> numpy.ndarray:
        """
        Azimuthal dependency function.

        Args:
            phi0(float): Phase (rotation) of the field.

        Returns:
            2D array of values (ndarray). Values are between -1 and 1.

        """
        return -numpy.sin(self.mode.nu * self.phi_mesh + phi)

    def wrapper_get_field(function):
        def wrapper(self, *args, **kwargs):
            return function(self, *args, **kwargs)

        return wrapper

    # @wrapper_get_field

    def _get_field_for_LP_mode(self, phi: float, theta: float) -> numpy.ndarray:
        if self.mode.family is ModeFamily.LP:
            data = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)
            iterator = numpy.nditer(data, flags=['multi_index'])

            for _ in iterator:
                index = iterator.multi_index

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                data[index] = er[0] * azimuthal_dependency_f[index]

            self._Ex = data
            return data
        else:
            polarisation = self.Epol(phi, theta)
            return self.Et(phi, theta) * numpy.cos(polarisation)

    def Ex(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        X component of the E field.y component of the E field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the x-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            data = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)
            iterator = numpy.nditer(data, flags=['multi_index'])

            for _ in iterator:
                index = iterator.multi_index

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                data[index] = er[0] * azimuthal_dependency_f[index]

            self._Ex = data
            return data
        else:
            polarisation = self.Epol(phi, theta)
            return self.Et(phi, theta) * numpy.cos(polarisation)

    def Ey(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """
        Y component of the E field.y component of the E field.

        :param      phi:    The phase in radian
        :type       phi:    float
        :param      theta:  The orientation in radian
        :type       theta:  float

        :returns:   The field in the y-direction
        :rtype:     numpy.ndarray
        """
        if self.mode.family is ModeFamily.LP:
            data = numpy.zeros(self.x_mesh.shape)
            azimuthal_dependency_f = self.get_azimuthal_dependency_f(phi=phi)
            iterator = numpy.nditer(data, flags=['multi_index'])

            for _ in iterator:
                index = iterator.multi_index

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[index]
                )

                data[index] = er[1] * azimuthal_dependency_f[index]
            self._Ey = data
            return data
        else:
            polarisation = self.Epol(phi, theta)
            return self.Et(phi, theta) * numpy.sin(polarisation)

    def Ez(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """z component of the E field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        self._Ez = numpy.zeros(self.x_mesh.shape)
        f = self.get_azimuthal_dependency_f(phi=phi)

        iterator = numpy.nditer(self._Ez, flags=['multi_index'])
        for _ in iterator:
            index = iterator.multi_index

            er, hr = self.fiber._rfield(
                mode=self.mode,
                wavelength=self.wavelength,
                radius=self.r_mesh[index]
            )
            self._Ez[index] = er[2] * f[index]

        return self._Ez

    def Er(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """r component of the E field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
            return self.Et(phi, theta) * numpy.cos(self.Epol(phi, theta) - self.phi_mesh)
        else:
            self._Er = numpy.zeros(self.x_mesh.shape)
            f = self.get_azimuthal_dependency_f(phi=phi)
            for i, j in product(range(self.n_point), repeat=2):
                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[j, i]
                )

                self._Er[j, i] = er[0] * f[j, i]
            return self._Er

    def Ephi(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """phi component of the E field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
            return self.Et(phi, theta) * numpy.sin(self.Epol(phi, theta) - self.phi_mesh)
        else:
            self._Ephi = numpy.zeros(self.x_mesh.shape)
            g = self.get_azimuthal_dependency_g(phi=phi)
            for i, j in product(range(self.n_point), repeat=2):

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[j, i]
                )

                self._Ephi[j, i] = er[1] * g[j, i]
            return self._Ephi

    def Et(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """transverse component of the E field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

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
        """polarization of the transverse E field (in radians).

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

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
        """x component of the H field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
            self._Hx = numpy.zeros(self.x_mesh.shape)
            f = self.get_azimuthal_dependency_f(phi=phi)
            for i, j in product(range(self.n_point), repeat=2):

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[j, i]
                )

                self._Hx[j, i] = hr[0] * f[j, i]
            return self._Hx
        else:
            return self.Ht(phi, theta) * numpy.cos(self.Hpol(phi, theta))

    def Hy(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """y component of the H field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
            self._Hy = numpy.zeros(self.x_mesh.shape)
            f = self.get_azimuthal_dependency_f(phi=phi)
            for i, j in product(range(self.n_point), repeat=2):

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[j, i]
                )

                self._Hy[j, i] = hr[1] * f[j, i]
            return self._Hy
        else:
            return self.Ht(phi, theta) * numpy.sin(self.Hpol(phi, theta))

    def Hz(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """z component of the H field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        self._Hz = numpy.zeros(self.x_mesh.shape)
        f = self.get_azimuthal_dependency_f(phi=phi)
        for i, j in product(range(self.n_point), repeat=2):

            er, hr = self.fiber._rfield(
                mode=self.mode,
                wavelength=self.wavelength,
                radius=self.r_mesh[j, i]
            )

            self._Hz[j, i] = hr[2] * f[j, i]
        return self._Hz

    def Hr(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """r component of the H field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
            radial = self.Ht(phi, theta)
            azimuthal = numpy.cos(self.Hpol(phi, theta) - self.phi_mesh)
            return radial * azimuthal

        else:
            self._Hr = numpy.zeros(self.x_mesh.shape)
            f = self.get_azimuthal_dependency_f(phi=phi)
            for i, j in product(range(self.n_point), repeat=2):

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[j, i]
                )

                self._Hr[j, i] = hr[0] * f[j, i]
            return self._Hr

    def Hphi(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """phi component of the H field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

        """
        if self.mode.family is ModeFamily.LP:
            return self.Ht(phi, theta) * numpy.sin(self.Hpol(phi, theta) - self.phi_mesh)
        else:
            self._Hphi = numpy.zeros(self.x_mesh.shape)
            g = self.get_azimuthal_dependency_g(phi=phi)
            for i, j in product(range(self.n_point), repeat=2):

                er, hr = self.fiber._rfield(
                    mode=self.mode,
                    wavelength=self.wavelength,
                    radius=self.r_mesh[j, i]
                )

                self._Hphi[j, i] = hr[1] * g[j, i]

            return self._Hphi

    def Ht(self, phi: float = 0, theta: float = 0) -> numpy.ndarray:
        """transverse component of the H field.

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

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
        """polarization of the transverse H field (in radians).

        Args:
            phi: phase (in radians)
            theta: orientation (in radians)

        Return:
            (np x np) numpy array

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

        Suppose than r is large enough, such as \|F(r, r)\| = 0.

        """
        modF = self.Emod()

        return (numpy.square(numpy.sum(numpy.square(modF))) /
                numpy.sum(numpy.power(modF, 4))) * self.dx * self.dy

    def I(self) -> float:
        neff = self.fiber.neff(
            mode=HE11,
            wavelength=self.wavelength
        )

        nm = self.fiber.neff(
            mode=self.mode,
            wavelength=self.wavelength
        )

        norm_squared = numpy.sum(numpy.square(self.Et()))

        return nm / neff * norm_squared * self.dx * self.dy

    def N(self) -> float:
        """
        Normalization constant.
        """
        neff = self.fiber.neff(
            mode=HE11,
            wavelength=self.wavelength
        )

        return 0.5 * scipy.constants.epsilon_0 * neff * scipy.constants.c * self.I()

    def S(self):
        """
        Pointing vector
        """
        raise NotImplementedError('Not yet implemented')
