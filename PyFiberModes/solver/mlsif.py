from .solver import FiberSolver
from PyFiberModes import Wavelength, Mode, ModeFamily
from PyFiberModes import constants
from math import isnan
import numpy
from scipy.special import kn, kvp, k0, k1, jn, jvp, yn, yvp, iv, ivp


class NeffSolver(FiberSolver):

    def get_neff_lower_boundary(self,
            wavelength: Wavelength,
            mode: Mode,
            lower_neff_boundary: float,
            delta_neff: float) -> float:
        """
        Gets the lower boundary for neff value.

        :param      wavelength:           The wavelength
        :type       wavelength:           Wavelength
        :param      mode:                 The mode to evaluate
        :type       mode:                 Mode
        :param      lower_neff_boundary:  The lower neff boundary
        :type       lower_neff_boundary:  float
        :param      delta_neff:           The delta neff
        :type       delta_neff:           float

        :returns:   The neff lower boundary.
        :rtype:     float
        """
        if lower_neff_boundary is not None:
            return lower_neff_boundary

        lower_order_mode = None

        if mode.family is ModeFamily.HE:
            if mode.m > 1:
                lower_order_mode = Mode(ModeFamily.EH, mode.nu, mode.m - 1)

        elif mode.family is ModeFamily.EH:
            lower_order_mode = Mode(ModeFamily.HE, mode.nu, mode.m)

        elif mode.m > 1:
            lower_order_mode = Mode(mode.family, mode.nu, mode.m - 1)

        if lower_order_mode is None:
            lower_neff_boundary = self.fiber.get_maximum_index(wavelength=wavelength)

        else:
            lower_neff_boundary = self.fiber.get_effective_index(
                mode=lower_order_mode,
                wavelength=wavelength,
                delta_neff=delta_neff
            )

            if isnan(lower_neff_boundary):
                return lower_neff_boundary

        if mode.family is ModeFamily.LP and mode.nu > 0:
            pm = Mode(mode.family, mode.nu - 1, mode.m)

            lb = self.fiber.get_effective_index(
                mode=pm,
                wavelength=wavelength,
                delta_neff=delta_neff
            )

            if isnan(lb):
                return lb

            lower_neff_boundary = min(lower_neff_boundary, lb)

        return lower_neff_boundary

    def solve(self, wavelength: float, mode: Mode, delta_neff: float, lower_neff_boundary: float):
        wavelength = Wavelength(wavelength)

        lower_neff_boundary = self.get_neff_lower_boundary(
            wavelength=wavelength,
            mode=mode,
            lower_neff_boundary=lower_neff_boundary,
            delta_neff=delta_neff
        )

        last_layer = self.fiber.layers[-1]
        higher_neff_boundary = last_layer.get_minimum_index(wavelength=wavelength)

        match mode.family:
            case ModeFamily.LP:
                function = self._lpceq
            case ModeFamily.TE:
                function = self._teceq
            case ModeFamily.TM:
                function = self._tmceq
            case ModeFamily.HE:
                function = self._heceq
            case ModeFamily.EH:
                function = self._heceq

        if lower_neff_boundary <= higher_neff_boundary:
            print("impossible bound")
            return numpy.nan

        if (lower_neff_boundary - higher_neff_boundary) < 10 * delta_neff:
            delta_neff = (lower_neff_boundary - higher_neff_boundary) / 10

        return self.find_function_first_root(
            function=function,
            function_args=(wavelength, mode.nu),
            lowbound=lower_neff_boundary - 1e-15,
            highbound=higher_neff_boundary + 1e-15,
            delta=-delta_neff
        )

    def get_LP_field(self, wavelength: float, nu: int, neff: float, radius: float) -> tuple:
        N = len(self.fiber.layers)
        C = numpy.array((1, 0))

        for i in range(1, N):
            rho = self.fiber.get_inner_radius(i)
            if radius < rho:
                eval_layer = self.fiber.layers[i - 1]
                break

            current_layer = self.fiber.layers[i]
            previous_layer = self.fiber.layers[i - 1]

            A = previous_layer.Psi(rho, neff, wavelength, nu, C)
            C = current_layer.lpConstants(rho, neff, wavelength, nu, A)

        else:
            eval_layer = self.fiber.layers[-1]
            u = eval_layer.get_U_W_parameter(radius=rho, neff=neff, wavelength=wavelength)
            C = (0, A[0] / kn(nu, u))

        ex, _ = eval_layer.Psi(radius, neff, wavelength, nu, C)
        hy = neff * constants.Y0 * ex

        return numpy.array((ex, 0, 0)), numpy.array((0, hy, 0))

    def get_TE_field(self, wavelength: float, nu, neff, radius: float) -> tuple:
        pass

    def get_TM_field(self, wavelength: float, nu, neff, radius: float) -> tuple:
        n_layer = len(self.fiber.layers)
        C = numpy.array((1, 0))
        EH = numpy.zeros(4)
        radius_in = 0

        for i in range(n_layer - 1):
            radius_out = self.fiber.get_outer_radius(layer_idx=i)
            layer = self.fiber.layers[i]

            n = layer.get_maximum_index(wavelength=wavelength)

            u = layer.get_U_W_parameter(radius=radius_out, neff=neff, wavelength=wavelength)

            if i > 0:
                C = layer.tetmConstants(
                    radius_in=radius_in,
                    radius_out=radius_out,
                    neff=neff,
                    wavelength=wavelength,
                    EH=EH,
                    c=constants.Y0 * n**2,
                    idx=(0, 3)
                )

            if r < radius_out:
                break

            if neff < n:
                c1 = wavelength.k0 * radius_out / u
                F3 = jvp(nu, u) / jn(nu, u)
                F4 = yvp(nu, u) / yn(nu, u)
            else:
                c1 = -wavelength.k0 * radius_out / u
                F3 = ivp(nu, u) / iv(nu, u)
                F4 = kvp(nu, u) / kn(nu, u)

            c4 = constants.Y0 * n * n * c1

            EH[0] = C[0] + C[1]
            EH[3] = c4 * (F3 * C[0] + F4 * C[1])

            radius_in = radius_out
        else:
            last_layer = self.fiber.layers[-1]
            u = last_layer.get_U_W_parameter(radius=radius_out, neff=neff, wavelength=wavelength)

        return numpy.array((0, ephi, 0)), numpy.array((hr, 0, hz))

    def get_EH_field(self, wavelength: float, nu, neff, radius: float) -> tuple:
        return self.get_HE_field(
            wavelength=wavelength,
            nu=nu,
            neff=neff,
            radius=radius
        )

    def get_HE_field(self, wavelength: float, nu, neff, radius: float) -> tuple:
        self._heceq(neff=neff, wavelength=wavelength, nu=nu)

        for current_layer_idx, rho in enumerate(self.fiber.layer_radius):
            if radius < rho:
                break
        else:
            current_layer_idx += 1

        layer = self.fiber.layers[current_layer_idx]
        n = layer.get_maximum_index(wavelength=wavelength)
        u = layer.get_U_W_parameter(radius=rho, neff=neff, wavelength=wavelength)
        urp = u * radius / rho

        c1 = rho / u
        c2 = wavelength.k0 * c1
        c3 = nu * c1 / radius if radius else 0  # To avoid div by 0
        c6 = constants.Y0 * n * n

        if neff < n:
            B1 = jn(nu, u)
            B2 = yn(nu, u)
            F1 = jn(nu, urp) / B1
            F2 = yn(nu, urp) / B2 if current_layer_idx > 0 else 0
            F3 = jvp(nu, urp) / B1
            F4 = yvp(nu, urp) / B2 if current_layer_idx > 0 else 0
        else:
            c2 = -c2
            B1 = iv(nu, u)
            B2 = kn(nu, u)
            F1 = iv(nu, urp) / B1
            F2 = kn(nu, urp) / B2 if current_layer_idx > 0 else 0
            F3 = ivp(nu, urp) / B1
            F4 = kvp(nu, urp) / B2 if current_layer_idx > 0 else 0

        A, B, Ap, Bp = layer.C[:, 0] + layer.C[:, 1] * self.alpha

        Ez = A * F1 + B * F2
        Ezp = A * F3 + B * F4
        Hz = Ap * F1 + Bp * F2
        Hzp = Ap * F3 + Bp * F4

        if radius == 0 and nu == 1:
            # Asymptotic expansion of Ez (or Hz):
            # J1(ur/p)/r (r->0) = u/(2p)
            if neff < n:
                f = 1 / (2 * jn(nu, u))
            else:
                f = 1 / (2 * iv(nu, u))
            c3ez = A * f
            c3hz = Ap * f
        else:
            c3ez = c3 * Ez
            c3hz = c3 * Hz

        Er = c2 * (neff * Ezp - constants.eta0 * c3hz)
        Ep = c2 * (neff * c3ez - constants.eta0 * Hzp)

        Hr = c2 * (neff * Hzp - c6 * c3ez)
        Hp = c2 * (-neff * c3hz + c6 * Ezp)

        return numpy.array((Er, Ep, Ez)), numpy.array((Hr, Hp, Hz))

    def _lpceq(self, neff: float, wavelength: float, nu: int) -> float:
        n_layer = len(self.fiber.layers)
        C = numpy.zeros((n_layer - 1, 2))
        C[0, 0] = 1

        for idx in range(1, n_layer - 1):
            current_layer = self.fiber.layers[idx]
            previous_layer = self.fiber.layers[idx - 1]

            A = previous_layer.Psi(
                radius=current_layer.radius_in,
                neff=neff,
                wavelength=wavelength,
                nu=nu,
                C=C[idx - 1, :]
            )

            C[idx, :] = current_layer.lpConstants(
                radius=current_layer.radius_in,
                neff=neff,
                wavelength=wavelength,
                nu=nu,
                A=A
            )

        last_layer = self.fiber.layers[-1]
        penultimate_layer = self.fiber.layers[-2]

        A = penultimate_layer.Psi(
            radius=last_layer.radius_in,
            neff=neff,
            wavelength=wavelength,
            nu=nu,
            C=C[-1, :]
        )

        u = last_layer.get_U_W_parameter(
            radius=last_layer.radius_in,
            neff=neff,
            wavelength=wavelength
        )

        return u * kvp(nu, u) * A[0] - kn(nu, u) * A[1]

    def _teceq(self, neff: float, wavelength: float, nu: int) -> float:
        EH = numpy.empty(4)

        for layer in self.fiber.layers[:-1]:
            layer.EH_fields(
                radius_in=layer.radius_in,
                radius_out=layer.radius_out,
                nu=nu,
                neff=neff,
                wavelength=wavelength,
                EH=EH,
                TM=False
            )

        # Last layer
        _, Hz, Ep, _ = EH
        last_layer = self.fiber.layers[-1]
        u = last_layer.get_U_W_parameter(
            radius=last_layer.radius_in,
            neff=neff,
            wavelength=wavelength
        )

        F4 = k1(u) / k0(u)
        return Ep + wavelength.k0 * last_layer.radius_in / u * constants.eta0 * Hz * F4

    def _tmceq(self, neff: float, wavelength: float, nu) -> float:
        EH = numpy.empty(4)

        for layer in self.fiber.layers[:-1]:
            layer.EH_fields(
                radius_in=layer.radius_in,
                radius_out=layer.radius_out,
                nu=nu,
                neff=neff,
                wavelength=wavelength,
                EH=EH,
                TM=True
            )

        # Last layer
        Ez, _, _, Hp = EH
        last_layer = self.fiber.layers[-1]

        u = last_layer.get_U_W_parameter(
            radius=last_layer.radius_in,
            neff=neff,
            wavelength=wavelength
        )

        last_layer_index = last_layer.get_maximum_index(wavelength=wavelength)

        F4 = k1(u) / k0(u)

        return Hp - wavelength.k0 * last_layer.radius_in / u * constants.Y0 * last_layer_index**2 * Ez * F4

    def _heceq(self, neff: float, wavelength: float, nu: int) -> float:
        EH = numpy.empty((4, 2))

        for layer in self.fiber.layers[:-1]:
            try:
                layer.EH_fields(
                    radius_in=layer.radius_in,
                    radius_out=layer.radius_out,
                    nu=nu,
                    neff=neff,
                    wavelength=wavelength,
                    EH=EH
                )
            except ZeroDivisionError:
                return numpy.inf

        # Last layer
        C = numpy.zeros((4, 2))
        C[1, :] = EH[0, :]
        C[3, :] = EH[1, :]

        last_layer = self.fiber.layers[-1]

        last_layer.C = C

        u = last_layer.get_U_W_parameter(
            radius=last_layer.radius_in,
            neff=neff,
            wavelength=wavelength
        )

        last_layer_index = last_layer.get_maximum_index(wavelength=wavelength)

        F4 = kvp(nu, u) / kn(nu, u)
        c1 = -wavelength.k0 * last_layer.radius_in / u
        c2 = neff * nu / u * c1
        c3 = constants.eta0 * c1
        c4 = constants.Y0 * last_layer_index**2 * c1

        E = EH[2, :] - (c2 * EH[0, :] - c3 * F4 * EH[1, :])
        H = EH[3, :] - (c4 * F4 * EH[0, :] - c2 * EH[1, :])

        if E[1] != 0:
            self.alpha = -E[0] / E[1]
        else:
            self.alpha = -H[0] / H[1]

        return E[0] * H[1] - E[1] * H[0]

    def _ehceq(self, neff: float, wavelength: float, nu: int) -> float:
        return self._heceq(
            neff=neff,
            wavelength=wavelength,
            nu=nu
        )

# -
