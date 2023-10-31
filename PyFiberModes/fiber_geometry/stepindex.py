from PyFiberModes.fiber_geometry.geometry import Geometry
from PyFiberModes import constants
from math import sqrt
import numpy
from scipy.special import jn, yn, iv, kn
from scipy.special import j0, y0, i0, k0
from scipy.special import j1, y1, i1, k1
from scipy.special import jvp, yvp, ivp, kvp


class StepIndex(Geometry):
    DEFAULT_PARAMS = []

    def index(self, radius: float, wavelength: float) -> float:
        if self.ri <= abs(radius) <= self.ro:
            return self._m.n(wavelength, *self._mp)
        else:
            return None

    def get_minimum_index(self, wavelength: float) -> float:
        return self._m.n(wavelength, *self._mp)

    def get_maximum_index(self, wavelength: float) -> float:
        return self._m.n(wavelength, *self._mp)

    def u(self, radius: float, neff, wavelength: float) -> float:
        return wavelength.k0 * radius * sqrt(abs(self.index(radius, wavelength)**2 - neff**2))

    def Psi(self, radius: float, neff, wavelength: float, nu, C):
        u = self.u(radius, neff, wavelength)
        if neff < self.get_maximum_index(wavelength):
            psi = (C[0] * jn(nu, u) + C[1] * yn(nu, u) if C[1] else
                   C[0] * jn(nu, u))
            psip = u * (C[0] * jvp(nu, u) + C[1] * yvp(nu, u) if C[1] else
                        C[0] * jvp(nu, u))
        else:
            psi = (C[0] * iv(nu, u) + C[1] * kn(nu, u) if C[1] else
                   C[0] * iv(nu, u))
            psip = u * (C[0] * ivp(nu, u) + C[1] * kvp(nu, u) if C[1] else
                        C[0] * ivp(nu, u))
        # if numpy.isnan(psi):
        #     print(neff, self.get_maximum_index(wavelength), C, r)
        return psi, psip

    def lpConstants(self, radius: float, neff, wavelength: float, nu, A):
        u = self.u(
            radius=radius,
            neff=neff,
            wavelength=wavelength
        )

        if neff < self.get_maximum_index(wavelength):
            W = constants.pi / 2
            return (W * (u * yvp(nu, u) * A[0] - yn(nu, u) * A[1]),
                    W * (jn(nu, u) * A[1] - u * jvp(nu, u) * A[0]))
        else:
            return ((u * kvp(nu, u) * A[0] - kn(nu, u) * A[1]),
                    (iv(nu, u) * A[1] - u * ivp(nu, u) * A[0]))

    def EH_fields(self, ri, ro, nu, neff, wavelength: float, EH, tm=True):
        """

        modify EH in-place (for speed)

        """
        maximum_index = self.get_maximum_index(wavelength=wavelength)

        u = self.u(
            radius=ro,
            neff=neff,
            wavelength=wavelength
        )

        if ri == 0:
            if nu == 0:
                if tm:
                    self.C = numpy.array([1., 0., 0., 0.])
                else:
                    self.C = numpy.array([0., 0., 1., 0.])
            else:
                self.C = numpy.zeros((4, 2))
                self.C[0, 0] = 1  # Ez = 1
                self.C[2, 1] = 1  # Hz = alpha
        elif nu == 0:
            self.C = numpy.zeros(4)
            if tm:
                c = constants.Y0 * maximum_index**2
                idx = (0, 3)
                self.C[:2] = self.tetmConstants(ri, ro, neff, wavelength, EH, c, idx)
            else:
                c = -constants.eta0
                idx = (1, 2)
                self.C[2:] = self.tetmConstants(ri, ro, neff, wavelength, EH, c, idx)
        else:
            self.C = self.vConstants(ri, ro, neff, wavelength, nu, EH)

        # Compute EH fields
        if neff < maximum_index:
            c1 = wavelength.k0 * ro / u
            F3 = jvp(nu, u) / jn(nu, u)
            F4 = yvp(nu, u) / yn(nu, u)
        else:
            c1 = -wavelength.k0 * ro / u
            F3 = ivp(nu, u) / iv(nu, u)
            F4 = kvp(nu, u) / kn(nu, u)

        c2 = neff * nu / u * c1
        c3 = constants.eta0 * c1
        c4 = constants.Y0 * n * n * c1

        EH[0] = self.C[0] + self.C[1]
        EH[1] = self.C[2] + self.C[3]
        EH[2] = (c2 * (self.C[0] + self.C[1]) -
                 c3 * (F3 * self.C[2] + F4 * self.C[3]))
        EH[3] = (c4 * (F3 * self.C[0] + F4 * self.C[1]) -
                 c2 * (self.C[2] + self.C[3]))

        return EH

    def vConstants(self, ri: float, ro: float, neff, wavelength: float, nu, EH):
        a = numpy.zeros((4, 4))

        maximum_index = self.get_maximum_index(wavelength)

        u = self.u(
            radius=ro,
            neff=neff,
            wavelength=wavelength
        )

        urp = self.u(
            radius=ri,
            neff=neff,
            wavelength=wavelength
        )

        if neff < maximum_index:
            B1 = jn(nu, u)
            B2 = yn(nu, u)
            F1 = jn(nu, urp) / B1
            F2 = yn(nu, urp) / B2
            F3 = jvp(nu, urp) / B1
            F4 = yvp(nu, urp) / B2
            c1 = wavelength.k0 * ro / u
        else:
            B1 = iv(nu, u)
            B2 = kn(nu, u)
            F1 = iv(nu, urp) / B1 if u else 1
            F2 = kn(nu, urp) / B2
            F3 = ivp(nu, urp) / B1 if u else 1
            F4 = kvp(nu, urp) / B2
            c1 = -wavelength.k0 * ro / u
        c2 = neff * nu / urp * c1
        c3 = constants.eta0 * c1
        c4 = constants.Y0 * n * n * c1

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

    def tetmConstants(self, ri: float, ro: float, neff, wavelength: float, EH, c, idx):
        a = numpy.empty((2, 2))

        maximum_index = self.get_maximum_index(wavelength)

        u = self.u(
            radius=ro,
            neff=neff,
            wavelength=wavelength
        )

        urp = self.u(
            radius=ri,
            neff=neff,
            wavelength=wavelength
        )

        if neff < maximum_index:
            B1 = j0(u)
            B2 = y0(u)
            F1 = j0(urp) / B1
            F2 = y0(urp) / B2
            F3 = -j1(urp) / B1
            F4 = -y1(urp) / B2
            c1 = wavelength.k0 * ro / u
        else:
            B1 = i0(u)
            B2 = k0(u)
            F1 = i0(urp) / B1
            F2 = k0(urp) / B2
            F3 = i1(urp) / B1
            F4 = -k1(urp) / B2
            c1 = -wavelength.k0 * ro / u
        c3 = c * c1

        a[0, 0] = F1
        a[0, 1] = F2
        a[1, 0] = F3 * c3
        a[1, 1] = F4 * c3

        return numpy.linalg.solve(a, EH.take(idx))
