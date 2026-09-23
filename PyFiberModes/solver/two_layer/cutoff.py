"""Modal cutoff solver for two-layer step-index fibers."""

import numpy
import logging

from PyFiberModes.solver.base_solver import BaseSolver
from PyFiberModes.mode import Mode


from scipy.special import jn, jn_zeros
from scipy.constants import mu_0, epsilon_0, physical_constants
eta0 = physical_constants['characteristic impedance of vacuum'][0]
Y0 = numpy.sqrt(epsilon_0 / mu_0)

class CutoffSolver(BaseSolver):
    """Solve normalized modal cutoffs for a two-layer step-index fiber.

    Parameters
    ----------
    fiber : Fiber
        Standard step-index fiber.
    wavelength : float
        Reference vacuum wavelength in meters.
    """
    logger = logging.getLogger(__name__)

    def solve(self, mode: Mode) -> float:
        """Calculate the normalized cutoff frequency of a mode.

        Parameters
        ----------
        mode : Mode
            Mode whose cutoff is requested.

        Returns
        -------
        float
            Cutoff V-number.
        """
        nu = mode.nu
        m = mode.m

        if mode.family == 'LP':
            if nu == 0:
                nu = 1
                m -= 1

            else:
                nu -= 1

        elif mode.family == 'HE':
            if nu == 1:
                m -= 1
            else:
                return self.find_HE_mode_cutoff(mode)

        return jn_zeros(nu, m)[m - 1]

    def _get_mode_cutoff_HE(self, V0: float, nu: int, mode: Mode) -> float:
        """Evaluate an HE cutoff while accounting for material dispersion.

        Parameters
        ----------
        V0 : float
            Trial normalized frequency.
        nu : int
            Azimuthal order.
        mode : Mode
            Hybrid electric mode used to obtain cutoff wavelength.

        Returns
        -------
        float
            Dispersion-adjusted cutoff residual.
        """
        core, clad = self.fiber.layers

        cutoff_wavelength = self.fiber.get_mode_cutoff_wavelength(mode=mode)

        normal_wavelength = core.wavelength

        self.fiber.update_wavelength(cutoff_wavelength)

        n_core = core.refractive_index

        n_clad = clad.refractive_index

        ratio = n_core**2 / n_clad**2

        self.fiber.update_wavelength(normal_wavelength)

        return (1 + ratio) * jn(nu - 2, V0) - (1 - ratio) * jn(nu, V0)

    def get_mode_cutoff_HE(self, V0, nu):
        """Evaluate the hybrid-mode cutoff characteristic equation.

        Parameters
        ----------
        V0 : float
            Trial normalized frequency.
        nu : int
            Azimuthal order.

        Returns
        -------
        float
            Characteristic-equation residual.
        """
        core, clad = self.fiber.layers

        n_ratio = core.refractive_index**2 / clad.refractive_index**2

        return (1 + n_ratio) * jn(nu - 2, V0) - (1 - n_ratio) * jn(nu, V0)

    def find_HE_mode_cutoff(self, mode: Mode) -> float:
        """Find a higher-order HE-mode cutoff by bracketing Bessel roots.

        Parameters
        ----------
        mode : Mode
            Hybrid electric mode.

        Returns
        -------
        float
            Cutoff V-number, or zero if no root is found.
        """
        if mode.m > 1:
            lower_neff_mode = Mode(
                family=mode.family,
                nu=mode.nu,
                m=mode.m - 1
            )

            lower_neff_boundary = self.fiber.get_mode_cutoff_v0(mode=lower_neff_mode)

            if numpy.isnan(lower_neff_boundary) or numpy.isinf(lower_neff_boundary):
                raise AssertionError(f"find_HE_mode_cutoff: no previous cutoff for {mode} mode")

            delta = 1 / lower_neff_boundary if lower_neff_boundary else self._MCD

            lower_neff_boundary += delta
        else:
            lower_neff_boundary = delta = self._MCD

        ipoints = numpy.concatenate(
            [jn_zeros(mode.nu, mode.m), jn_zeros(mode.nu - 2, mode.m)]
        )

        ipoints.sort()
        ipoints = list(ipoints[ipoints > lower_neff_boundary])

        cutoff = self.find_function_first_root(
            function=self.get_mode_cutoff_HE,
            function_args=(mode.nu,),
            lowbound=lower_neff_boundary,
            ipoints=ipoints,
            delta=delta
        )

        if numpy.isnan(cutoff):
            self.logger.error(f"find_HE_mode_cutoff: no cutoff found for {mode} mode")
            return 0

        return cutoff

# -
