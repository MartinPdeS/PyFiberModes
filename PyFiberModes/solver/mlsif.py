#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFiberModes.solver.solver import FiberSolver
from PyFiberModes import Mode, ModeFamily
import numpy
from scipy.special import kn, kvp, k0, k1, jn, jvp, yn, yvp, iv, ivp
from scipy.constants import mu_0, epsilon_0, physical_constants
eta0 = physical_constants['characteristic impedance of vacuum'][0]


class NameSpace():
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


class NeffSolver(FiberSolver):

    def get_neff_lower_boundary(self,
            mode: Mode,
            lower_neff_boundary: float,
            delta_neff: float) -> float:
        """
        Gets the lower boundary for neff value.

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
            lower_neff_boundary = self.fiber.get_maximum_index()



        else:
            lower_neff_boundary = self.fiber.get_effective_index(
                mode=lower_order_mode,
                delta_neff=delta_neff
            )

            if numpy.isnan(lower_neff_boundary):
                return lower_neff_boundary

        if mode.family is ModeFamily.LP and mode.nu > 0:
            pm = Mode(mode.family, mode.nu - 1, mode.m)

            lb = self.fiber.get_effective_index(
                mode=pm,
                delta_neff=delta_neff
            )

            if numpy.isnan(lb):
                return lb

            lower_neff_boundary = min(lower_neff_boundary, lb)

        return lower_neff_boundary

    def solve(self, mode: Mode, delta_neff: float, lower_neff_boundary: float) -> float:
        lower_neff_boundary = self.get_neff_lower_boundary(
            mode=mode,
            lower_neff_boundary=lower_neff_boundary,
            delta_neff=delta_neff
        )

        last_layer = self.fiber.layers[-1]
        higher_neff_boundary = last_layer.refractive_index

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
            print("Impossible bound")
            return numpy.nan

        if (lower_neff_boundary - higher_neff_boundary) < 10 * delta_neff:
            delta_neff = (lower_neff_boundary - higher_neff_boundary) / 10

        return self.find_function_first_root(
            function=function,
            function_args=(mode.nu,),
            lowbound=lower_neff_boundary - 1e-15,
            highbound=higher_neff_boundary + 1e-15,
            delta=-delta_neff
        )

    def get_LP_field_for_future(self, nu: int, neff: float, radius: float) -> tuple[float, float]:
        """
        Gets the :math:`LP_{\nu, m}` mode field.

        :param      nu:      The nu parameter of the LP mode
        :type       nu:      int
        :param      neff:    The effective index
        :type       neff:    float
        :param      radius:  The radius for evaluation
        :type       radius:  float

        :returns:   The LP electric and magnetic field in a tuple.
        :rtype:     tuple
        """
        n_layers = len(self.fiber.layers)
        C = numpy.array((1, 0))

        test = NameSpace(layer=[], index=[], C=[])

        for i in range(1, n_layers):
            layer_out = self.fiber.layers[i]
            layer_in = self.fiber.layers[i - 1]

            test.layer.append(layer_in)
            test.index.append(radius < layer_in.radius_out)
            test.C.append(C)

            # if radius < layer_in.radius_out:
            #     eval_layer = self.fiber.layers[i - 1]
            #     break

            A = layer_in.get_psi(
                radius=layer_in.radius_out,
                neff=neff,
                nu=nu,
                C=C
            )

            C = layer_out.get_LP_constants(
                radius=layer_in.radius_out,
                neff=neff,
                nu=nu,
                A=A
            )

        else:
            eval_layer = self.fiber.layers[-1]

            u = eval_layer.get_U_W_parameter(
                radius=layer_in.radius_out,
                neff=neff
            )

            C = (0, A[0] / kn(nu, u))

            test.layer.append(eval_layer)
            test.index.append(None)
            test.C.append(C)

        E_x = numpy.zeros(radius.shape)
        for layer, index, C in zip(test.layer, test.index, test.C):
            pass

        ex, _ = eval_layer.get_psi(
            radius=radius,
            neff=neff,
            nu=nu,
            C=C
        )

        hy = neff * numpy.sqrt(epsilon_0 / mu_0) * ex

        e_field = numpy.array((ex, 0, 0))
        h_field = numpy.array((0, hy, 0))

        return e_field, h_field

    def get_LP_field(self, nu: int, neff: float, radius: float) -> tuple[float, float]:
        """
        Gets the :math:`LP_{\nu, m}` mode field.

        :param      nu:      The nu parameter of the LP mode
        :type       nu:      int
        :param      neff:    The effective index
        :type       neff:    float
        :param      radius:  The radius for evaluation
        :type       radius:  float

        :returns:   The LP electric and magnetic field in a tuple.
        :rtype:     tuple[float, float]
        """
        C = numpy.array((1, 0))

        for layer_in, layer_out in self.fiber.iterate_interfaces():
            if radius < layer_in.radius_out:
                eval_layer = layer_in
                break

            A = layer_in.get_psi(
                radius=layer_in.radius_out,
                neff=neff,
                nu=nu,
                C=C
            )

            C = layer_out.get_LP_constants(
                radius=layer_in.radius_out,
                neff=neff,
                nu=nu,
                A=A
            )

        else:
            eval_layer = self.fiber.last_layer

            u = eval_layer.get_U_W_parameter(
                radius=eval_layer.radius_in,
                neff=neff
            )

            C = (0, A[0] / kn(nu, u))

        return self.get_LP_field_from_parameters(
            eval_layer=eval_layer,
            radius=radius,
            neff=neff,
            nu=nu,
            C=C
        )

    def get_LP_field_from_parameters(self,
            eval_layer: object,
            radius: float,
            neff: float,
            nu: int,
            C: tuple) -> tuple[float, float]:
        """
        Gets the LP field evaluation from parameters.

        :param      eval_layer:  The layer at which the field is evaluated
        :type       eval_layer:  object
        :param      nu:          The nu parameter of the LP mode
        :type       nu:          int
        :param      neff:        The effective index
        :type       neff:        float
        :param      radius:      The radius for evaluation
        :type       radius:      float
        :param      C:           Constants
        :type       C:           tuple

        :returns:   The LP field.
        :rtype:     tuple[float, float]
        """
        ex, _ = eval_layer.get_psi(
            radius=radius,
            neff=neff,
            nu=nu,
            C=C
        )

        hy = neff * numpy.sqrt(epsilon_0 / mu_0) * ex

        e_field = numpy.array((ex, 0, 0))
        h_field = numpy.array((0, hy, 0))

        return e_field, h_field

    def get_TE_field(self, wavelength: float, nu: int, neff: float, radius: float) -> tuple[float, float]:
        pass

    def get_TM_field(self, wavelength: float, nu: int, neff: float, radius: float) -> tuple[float, float]:
        n_layer = len(self.fiber.layers)
        C = numpy.array((1, 0))
        EH = numpy.zeros(4)
        radius_in = 0

        for i in range(n_layer - 1):
            radius_out = self.fiber.get_outer_radius(layer_idx=i)
            layer = self.fiber.layers[i]

            n = layer.refractive_index

            u = layer.get_U_W_parameter(radius=radius_out, neff=neff)

            if i > 0:
                C = layer.get_TE_TM_constants(
                    radius_in=radius_in,
                    radius_out=radius_out,
                    neff=neff,
                    EH=EH,
                    c=numpy.sqrt(epsilon_0 / mu_0) * n**2,
                    idx=(0, 3)
                )

            if radius < radius_out:
                break

            if neff < n:
                c1 = wavelength.k0 * radius_out / u
                F3 = jvp(nu, u) / jn(nu, u)
                F4 = yvp(nu, u) / yn(nu, u)
            else:
                c1 = -wavelength.k0 * radius_out / u
                F3 = ivp(nu, u) / iv(nu, u)
                F4 = kvp(nu, u) / kn(nu, u)

            c4 = numpy.sqrt(epsilon_0 / mu_0) * n * n * c1

            EH[0] = C[0] + C[1]
            EH[3] = c4 * (F3 * C[0] + F4 * C[1])

            radius_in = radius_out
        else:
            last_layer = self.fiber.layers[-1]
            u = last_layer.get_U_W_parameter(radius=radius_out, neff=neff)

        return numpy.array((0, ephi, 0)), numpy.array((hr, 0, hz))

    def get_EH_field(self, nu, neff, radius: float) -> tuple[float, float]:
        return self.get_HE_field(
            nu=nu,
            neff=neff,
            radius=radius
        )

    def get_HE_field(self, nu, neff, radius: float) -> tuple[float, float]:
        self._heceq(neff=neff, nu=nu)

        layer = self.fiber.get_layer_at_radius(radius)

        rho = layer.radius_out if not layer.is_last_layer else self.fiber.penultimate_layer.radius_out

        u = layer.get_U_W_parameter(radius=rho, neff=neff)

        urp = u * radius / rho

        c1 = rho / u
        c2 = self.wavelength.k0 * c1
        c3 = nu * c1 / radius if radius else 0  # To avoid div by 0
        c6 = numpy.sqrt(epsilon_0 / mu_0) * layer.refractive_index**2

        if neff < layer.refractive_index:
            B1 = jn(nu, u)
            B2 = yn(nu, u)
            F1 = jn(nu, urp) / B1
            F2 = yn(nu, urp) / B2 if layer.radius_in > 0 else 0
            F3 = jvp(nu, urp) / B1
            F4 = yvp(nu, urp) / B2 if layer.radius_in > 0 else 0
        else:
            c2 = -c2
            B1 = iv(nu, u)
            B2 = kn(nu, u)
            F1 = iv(nu, urp) / B1
            F2 = kn(nu, urp) / B2 if layer.radius_in > 0 else 0
            F3 = ivp(nu, urp) / B1
            F4 = kvp(nu, urp) / B2 if layer.radius_in > 0 else 0

        A, B, Ap, Bp = layer.C[:, 0] + layer.C[:, 1] * self.alpha

        Ez = A * F1 + B * F2
        Ezp = A * F3 + B * F4
        Hz = Ap * F1 + Bp * F2
        Hzp = Ap * F3 + Bp * F4

        if radius == 0 and nu == 1:
            # Asymptotic expansion of Ez (or Hz):
            # J1(ur/p)/r (r->0) = u/(2p)
            if neff < layer.refractive_index:
                f = 1 / (2 * jn(nu, u))
            else:
                f = 1 / (2 * iv(nu, u))
            c3ez = A * f
            c3hz = Ap * f
        else:
            c3ez = c3 * Ez
            c3hz = c3 * Hz

        Er = c2 * (neff * Ezp - eta0 * c3hz)
        Ep = c2 * (neff * c3ez - eta0 * Hzp)

        Hr = c2 * (neff * Hzp - c6 * c3ez)
        Hp = c2 * (-neff * c3hz + c6 * Ezp)

        return numpy.array((Er, Ep, Ez)), numpy.array((Hr, Hp, Hz))

    def _lpceq(self, neff: float, nu: int) -> tuple[float, float]:
        n_layer = len(self.fiber.layers)
        C = numpy.zeros((n_layer - 1, 2))
        C[0, 0] = 1

        for idx in range(1, n_layer - 1):
            current_layer = self.fiber.layers[idx]
            previous_layer = self.fiber.layers[idx - 1]

            A = previous_layer.get_psi(
                radius=current_layer.radius_in,
                neff=neff,
                nu=nu,
                C=C[idx - 1, :]
            )

            C[idx, :] = current_layer.get_LP_constants(
                radius=current_layer.radius_in,
                neff=neff,
                nu=nu,
                A=A
            )

        last_layer = self.fiber.layers[-1]
        penultimate_layer = self.fiber.layers[-2]

        A = penultimate_layer.get_psi(
            radius=last_layer.radius_in,
            neff=neff,
            nu=nu,
            C=C[-1, :]
        )

        u = last_layer.get_U_W_parameter(
            radius=last_layer.radius_in,
            neff=neff,
        )

        return u * kvp(nu, u) * A[0] - kn(nu, u) * A[1]

    def _teceq(self, neff: float, nu: int) -> tuple[float, float]:
        EH = numpy.empty(4)

        for layer in self.fiber.layers[:-1]:
            layer.EH_fields(
                radius_in=layer.radius_in,
                radius_out=layer.radius_out,
                nu=nu,
                neff=neff,
                EH=EH,
                TM=False
            )

        # Last layer
        _, Hz, Ep, _ = EH
        u = self.fiber.last_layer.get_U_W_parameter(
            radius=self.fiber.last_layer.radius_in,
            neff=neff,
        )

        F4 = k1(u) / k0(u)
        return Ep + self.wavelength.k0 * self.fiber.last_layer.radius_in / u * eta0 * Hz * F4

    def _tmceq(self, neff: float, nu) -> tuple[float, float]:
        EH = numpy.empty(4)

        for layer in self.fiber.layers[:-1]:
            layer.EH_fields(
                radius_in=layer.radius_in,
                radius_out=layer.radius_out,
                nu=nu,
                neff=neff,
                EH=EH,
                TM=True
            )

        # Last layer
        Ez, _, _, Hp = EH
        last_layer = self.fiber.layers[-1]

        u = last_layer.get_U_W_parameter(
            radius=last_layer.radius_in,
            neff=neff,
        )

        last_layer_index = last_layer.refractive_index

        F4 = k1(u) / k0(u)

        return Hp - self.wavelength.k0 * last_layer.radius_in / u * numpy.sqrt(epsilon_0 / mu_0) * last_layer_index**2 * Ez * F4

    def _heceq(self, neff: float, nu: int) -> float:
        EH = numpy.empty((4, 2))

        for layer in self.fiber.layers[:-1]:
            layer.EH_fields(
                radius_in=layer.radius_in,
                radius_out=layer.radius_out,
                nu=nu,
                neff=neff,
                EH=EH
            )

        # Last layer
        C = numpy.zeros((4, 2))
        C[1, :] = EH[0, :]
        C[3, :] = EH[1, :]

        last_layer = self.fiber.layers[-1]

        last_layer.C = C

        u = last_layer.get_U_W_parameter(
            radius=last_layer.radius_in,
            neff=neff,
        )

        last_layer_index = last_layer.refractive_index

        F4 = kvp(nu, u) / kn(nu, u)
        c1 = -self.wavelength.k0 * last_layer.radius_in / u
        c2 = neff * nu / u * c1
        c3 = eta0 * c1
        c4 = numpy.sqrt(epsilon_0 / mu_0) * last_layer_index**2 * c1

        E = EH[2, :] - (c2 * EH[0, :] - c3 * F4 * EH[1, :])
        H = EH[3, :] - (c4 * F4 * EH[0, :] - c2 * EH[1, :])

        if E[1] != 0:
            self.alpha = -E[0] / E[1]
        else:
            self.alpha = -H[0] / H[1]

        return E[0] * H[1] - E[1] * H[0]

    def _ehceq(self, neff: float, nu: int) -> float:
        return self._heceq(
            neff=neff,
            nu=nu
        )

# -