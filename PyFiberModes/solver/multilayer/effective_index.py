"""Effective-index solver for multilayer step-index fibers."""

import numpy
from scipy.special import kn, kvp, k0, k1, jn, jvp, yn, yvp, iv, ivp
from scipy.constants import mu_0, epsilon_0, physical_constants

from PyFiberModes.solver.base_solver import BaseSolver
from PyFiberModes.mode import Mode

eta0 = physical_constants['characteristic impedance of vacuum'][0]


class EffectiveIndexSolver(BaseSolver):
    """Solve effective indices for arbitrary multilayer step-index fibers.

    Parameters
    ----------
    fiber : Fiber
        Multilayer fiber to solve.
    wavelength : float
        Vacuum wavelength in meters.
    """
    def get_neff_lower_boundary(self, mode: Mode, delta_neff: float = 1e-6) -> float:
        """Gets the lower boundary for neff value.

        Parameters
        ----------
        mode : Mode
            The mode to evaluate
        delta_neff : float
            The delta neff

        Returns
        -------
        float
            The neff lower boundary.
        """
        lower_order_mode = None

        if mode.family == 'HE':
            if mode.m > 1:
                lower_order_mode = Mode('EH', mode.nu, mode.m - 1)

        elif mode.family == 'EH':
            lower_order_mode = Mode('HE', mode.nu, mode.m)

        elif mode.m > 1:
            lower_order_mode = Mode(mode.family, mode.nu, mode.m - 1)

        if lower_order_mode is None:
            lower_neff_boundary = self.fiber.maximum_index

        else:
            lower_neff_boundary = self.fiber.get_effective_index(mode=lower_order_mode)

            if numpy.isnan(lower_neff_boundary):
                return lower_neff_boundary

        if mode.family == 'LP' and mode.nu > 0:
            pm = Mode(mode.family, mode.nu - 1, mode.m)

            lb = self.fiber.get_effective_index(mode=pm)

            if numpy.isnan(lb):
                return lb

            lower_neff_boundary = min(lower_neff_boundary, lb)

        return lower_neff_boundary

    def solve(self, mode: Mode, delta_neff: float) -> float:
        """Find an effective-index root for a mode.

        Parameters
        ----------
        mode : Mode
            Mode whose effective index is requested.
        delta_neff : float
            Maximum effective-index search step.

        Returns
        -------
        float
            Effective index, or ``numpy.nan`` when no valid interval exists.
        """
        higher_neff_boundary = self.get_neff_lower_boundary(mode=mode)

        lower_neff_boundary = self.fiber.last_layer.refractive_index

        match mode.family:
            case 'LP':
                function = self.get_LP_equation
            case 'TE':
                function = self.get_TE_equation
            case 'TM':
                function = self.get_TM_equation
            case 'HE':
                function = self.get_HE_equation
            case 'EH':
                function = self.get_HE_equation

        if higher_neff_boundary <= lower_neff_boundary:
            print("Impossible bound")
            return numpy.nan

        delta_boundary = (higher_neff_boundary - lower_neff_boundary) / 100
        delta_neff = - min(delta_neff, delta_boundary)

        extra = 1e-13

        try:
            value = self.find_function_first_root(
                function=function,
                function_args=(mode.nu,),
                lowbound=higher_neff_boundary - extra,
                highbound=lower_neff_boundary + extra,
                delta=delta_neff
            )

        except ValueError:
            value = numpy.nan

        return value

    def get_LP_field(self, nu: int, neff: float, radius: float) -> tuple[float, float]:
        """        Gets the :math:`LP_{
        u, m}` mode field.

        Parameters
        ----------
        nu : int
            The nu parameter of the LP mode
        neff : float
            The effective index
        radius : float
            The radius for evaluation

        Returns
        -------
        tuple[float, float]
            The LP electric and magnetic field in a tuple.
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

    def get_LP_field_from_parameters(
            self,
            eval_layer: object,
            radius: float,
            neff: float,
            nu: int,
            C: tuple) -> tuple[float, float]:
        """Gets the LP field evaluation from parameters.

        Parameters
        ----------
        eval_layer : object
            The layer at which the field is evaluated
        radius : float
            The radius for evaluation
        neff : float
            The effective index
        nu : int
            The nu parameter of the LP mode
        C : tuple
            Constants

        Returns
        -------
        tuple[float, float]
            The LP field.
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
        """Gets the transverse electric TE field.

        Parameters
        ----------
        wavelength : float
            Vacuum wavelength in meters.
        nu : int
            The radial parameter of the mode
        neff : float
            The effective index of the mode
        radius : float
            The radius at which field is evaluated

        Returns
        -------
        tuple[float, float]
            The TE field.

        Raises
        ------
        NotImplementedError
            Method not yet implemented
        """
        raise NotImplementedError()

    def get_TM_field(self, wavelength: float, nu: int, neff: float, radius: float) -> tuple[float, float]:
        """Gets the transverse magnetic TM field.

        Parameters
        ----------
        wavelength : float
            Vacuum wavelength in meters.
        nu : int
            The radial parameter of the mode
        neff : float
            The effective index of the mode
        radius : float
            The radius at which field is evaluated

        Returns
        -------
        tuple[float, float]
            The TM field.
        """
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
                c1 = (2 * numpy.pi / self.wavelength) * radius_out / u
                F3 = jvp(nu, u) / jn(nu, u)
                F4 = yvp(nu, u) / yn(nu, u)
            else:
                c1 = -(2 * numpy.pi / self.wavelength) * radius_out / u
                F3 = ivp(nu, u) / iv(nu, u)
                F4 = kvp(nu, u) / kn(nu, u)

            c4 = numpy.sqrt(epsilon_0 / mu_0) * n * n * c1

            EH[0] = C[0] + C[1]
            EH[3] = c4 * (F3 * C[0] + F4 * C[1])

            radius_in = radius_out
        else:
            u = self.fiber.last_layer.get_U_W_parameter(radius=radius_out, neff=neff)

        return numpy.array((0, ephi, 0)), numpy.array((hr, 0, hz))

    def get_EH_field(self, nu: int, neff: float, radius: float) -> tuple[float, float]:
        """Gets the hybrid EH field.

        Parameters
        ----------
        nu : int
            The radial parameter of the mode
        neff : float
            The effective index of the mode
        radius : float
            The radius at which field is evaluated

        Returns
        -------
        tuple[float, float]
            The EH field.
        """
        return self.get_HE_field(
            nu=nu,
            neff=neff,
            radius=radius
        )

    def get_HE_field(self, nu: int, neff: float, radius: float) -> tuple[float, float]:
        """Gets the hybrid HE field.

        Parameters
        ----------
        nu : int
            The radial parameter of the mode
        neff : float
            The effective index of the mode
        radius : float
            The radius at which field is evaluated

        Returns
        -------
        tuple[float, float]
            The HE field.
        """
        self.get_HE_equation(neff=neff, nu=nu)

        layer = self.fiber.get_layer_at_radius(radius)

        rho = layer.radius_out if not layer.is_last_layer else self.fiber.penultimate_layer.radius_out

        u = layer.get_U_W_parameter(radius=rho, neff=neff)

        urp = u * radius / rho

        c1 = rho / u
        c2 = (2 * numpy.pi / self.wavelength) * c1
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

    def get_LP_equation(self, neff: float, nu: int) -> tuple[float, float]:
        """Evaluate the multilayer LP characteristic equation.

        Parameters
        ----------
        neff : float
            Trial effective index.
        nu : int
            Azimuthal order.

        Returns
        -------
        float
            Characteristic-equation residual.
        """
        C = numpy.zeros((self.fiber.n_interface, 2))
        C[0, 0] = 1

        for layer_in, layer_out in self.fiber.iterate_interfaces():
            if layer_out.is_last_layer:
                continue

            A = layer_in.get_psi(
                radius=layer_out.radius_in,
                neff=neff,
                nu=nu,
                C=C[layer_in.position, :]
            )

            C[layer_out.position, :] = layer_out.get_LP_constants(
                radius=layer_out.radius_in,
                neff=neff,
                nu=nu,
                A=A
            )

        A = self.fiber.penultimate_layer.get_psi(
            radius=self.fiber.last_layer.radius_in,
            neff=neff,
            nu=nu,
            C=C[-1, :]
        )

        u = self.fiber.last_layer.get_U_W_parameter(
            radius=self.fiber.last_layer.radius_in,
            neff=neff,
        )

        return u * kvp(nu, u) * A[0] - kn(nu, u) * A[1]

    def get_TE_equation(self, neff: float, nu: int) -> tuple[float, float]:
        """Evaluate the multilayer TE characteristic equation.

        Parameters
        ----------
        neff : float
            Trial effective index.
        nu : int
            Azimuthal order.

        Returns
        -------
        float
            Characteristic-equation residual.
        """
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
        return Ep + (2 * numpy.pi / self.wavelength) * self.fiber.last_layer.radius_in / u * eta0 * Hz * F4

    def get_TM_equation(self, neff: float, nu: int) -> tuple[float, float]:
        """Evaluate the multilayer TM characteristic equation.

        Parameters
        ----------
        neff : float
            Trial effective index.
        nu : int
            Azimuthal order.

        Returns
        -------
        float
            Characteristic-equation residual.
        """
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

        # At last layer the equation are differents.
        Ez, _, _, Hp = EH

        u = self.fiber.last_layer.get_U_W_parameter(
            radius=self.fiber.last_layer.radius_in,
            neff=neff,
        )

        F4 = k1(u) / k0(u)

        return Hp - (2 * numpy.pi / self.wavelength) * self.fiber.last_layer.radius_in / u * numpy.sqrt(epsilon_0 / mu_0) * self.fiber.last_layer.refractive_index**2 * Ez * F4

    def get_HE_equation(self, neff: float, nu: int) -> float:
        """Evaluate the multilayer hybrid-mode characteristic equation.

        Parameters
        ----------
        neff : float
            Trial effective index.
        nu : int
            Azimuthal order.

        Returns
        -------
        float
            Determinant residual for the hybrid boundary conditions.
        """
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

        F4 = kvp(nu, u) / kn(nu, u)
        c1 = -(2 * numpy.pi / self.wavelength) * last_layer.radius_in / u
        c2 = neff * nu / u * c1
        c3 = eta0 * c1
        c4 = numpy.sqrt(epsilon_0 / mu_0) * last_layer.refractive_index**2 * c1

        E = EH[2, :] - (c2 * EH[0, :] - c3 * F4 * EH[1, :])
        H = EH[3, :] - (c4 * F4 * EH[0, :] - c2 * EH[1, :])

        if E[1] != 0:
            self.alpha = -E[0] / E[1]
        else:
            self.alpha = -H[0] / H[1]

        return E[0] * H[1] - E[1] * H[0]

    def get_EH_equation(self, neff: float, nu: int) -> float:
        """Evaluate the EH equation through the shared hybrid determinant.

        Parameters
        ----------
        neff : float
            Trial effective index.
        nu : int
            Azimuthal order.

        Returns
        -------
        float
            Hybrid characteristic-equation residual.
        """
        return self.get_HE_equation(
            neff=neff,
            nu=nu
        )

# -
