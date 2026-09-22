"""Shared bracketing and root-finding facilities for modal solvers."""

import logging
import numpy
from abc import ABC, abstractmethod

from scipy.optimize import brentq, root_scalar
from PyFiberModes.solver.results import SolverResult


class BaseSolver(ABC):
    """Provide shared scalar root-finding operations for fiber solvers.

    Parameters
    ----------
    fiber : Fiber
        Fiber whose characteristic equations are solved.
    wavelength : float
        Vacuum wavelength in meters.
    """

    logger = logging.getLogger(__name__)
    _MCD = 0.1

    def __init__(self, fiber, wavelength):
        """Store the fiber and wavelength used by a solver.

        Parameters
        ----------
        fiber : Fiber
            Fiber model to solve.
        wavelength : float
            Vacuum wavelength in meters.
        """
        self.fiber = fiber
        self.wavelength = wavelength

    @abstractmethod
    def solve(self, *args, **kwargs):
        """Solve a modal equation in a concrete subclass.

        Parameters
        ----------
        *args
            Positional arguments accepted by a concrete solver.
        **kwargs
            Keyword arguments accepted by a concrete solver.

        Raises
        ------
        NotImplementedError
            Always raised by the abstract base implementation.
        """
        raise NotImplementedError()

    @staticmethod
    def evaluate_many(function, values, function_args: tuple = ()) -> numpy.ndarray:
        """Evaluate a characteristic function over an array when possible.

        Parameters
        ----------
        function : callable
            Characteristic function.
        values : array-like
            Evaluation positions.
        function_args : tuple, optional
            Extra function arguments.

        Returns
        -------
        numpy.ndarray
            Function values. Scalar-only functions transparently fall back to
            a Python iteration.
        """
        values = numpy.asarray(values, dtype=float)
        try:
            result = numpy.asarray(function(values, *function_args))
            if result.shape == values.shape:
                return result
        except (TypeError, ValueError):
            pass
        return numpy.fromiter(
            (function(value, *function_args) for value in values),
            dtype=float,
            count=values.size,
        )

    def find_function_first_root(
            self,
            function,
            function_args: tuple = (),
            lowbound: float = 0,
            highbound: float = None,
            ipoints: list = [],
            delta: float = 0.25,
            maxiter: int = numpy.inf) -> float:
        """Find the first continuous sign-changing root in an interval.

        Parameters
        ----------
        function : callable
            Scalar function whose first root is requested.
        function_args : tuple, optional
            Extra positional arguments passed to ``function``.
        lowbound : float, optional
            Initial search position.
        highbound : float, optional
            Optional terminal search position.
        ipoints : list, optional
            Explicit successive search positions.
        delta : float, optional
            Step between implicit search positions.
        maxiter : int, optional
            Maximum number of search steps.

        Returns
        -------
        float
            First accepted root, or ``numpy.nan`` when none is found.

        Notes
        -----
        Candidate discontinuities are rejected by comparing the residual at
        the Brent root with residuals at both bracket endpoints.
        """

        while True:
            if ipoints:
                maxiter = len(ipoints)
            elif highbound:
                maxiter = int((highbound - lowbound) / delta)

            a = lowbound
            fa = function(a, *function_args)
            if fa == 0:
                return a

            for i in range(1, maxiter + 1):
                b = ipoints.pop(0) if ipoints else a + delta
                if highbound:
                    if (b > highbound > lowbound) or (b < highbound < lowbound):
                        self.logger.info("find_function_first_root: no root found within allowed range")
                        return numpy.nan

                fb = function(b, *function_args)

                if fb == 0:
                    return b

                if (fa > 0 and fb < 0) or (fa < 0 and fb > 0):
                    z = brentq(function, a, b, args=function_args, xtol=1e-20)

                    fz = function(z, *function_args)
                    if abs(fa) > abs(fz) < abs(fb):  # Skip discontinuities
                        self.logger.debug(f"skipped ({fa}, {fz}, {fb})")
                        return z

                a, fa = b, fb

            if highbound and maxiter < 100:
                delta /= 10
            else:
                break

        self.logger.info(f"maxiter reached ({maxiter}, {lowbound}, {highbound})")
        return numpy.nan

    def get_new_x_low_x_high(
            self,
            function,
            function_args,
            x_low: float,
            x_high: float,
            n_slice: int = 100) -> tuple:
        """Gets the new x boundaries.
        Returns numpy.nan if no sign inversion found.

        Parameters
        ----------
        function : { type_description }
            The function
        function_args : { type_description }
            The function arguments
        x_low : float
            The x low
        x_high : float
            The x high
        n_slice : int
            The n iteration

        Returns
        -------
        tuple
            The new x low x high.
        """
        x_list = [x_low, x_high]
        x_list.sort()
        x_list = numpy.linspace(*x_list, n_slice)
        y_list = self.evaluate_many(function, x_list, function_args)

        non_nan_idx = ~numpy.isnan(y_list)

        x_list = x_list[non_nan_idx]
        y_list = y_list[non_nan_idx]

        if len(y_list) < 2:
            return numpy.nan

        sign_change = (numpy.diff(numpy.sign(y_list)) != 0) * 1

        if (sign_change == 0).all():
            return numpy.nan

        sign_change_idx = numpy.where(sign_change == 1)[0]

        sign_change_idx = sign_change_idx[0]

        x_low, x_high = x_list[sign_change_idx], x_list[sign_change_idx + 1]

        y_low, y_high = y_list[sign_change_idx], y_list[sign_change_idx + 1]

        return x_low, x_high, y_low, y_high

    def find_root_within_range(
            self,
            function,
            x_low: float,
            x_high: float,
            function_args: tuple = (),
            max_iteration: int = 100,
            tolerance: float = 1e-8) -> float:
        """Finds and return the root of a given function within range.

        Parameters
        ----------
        function : object
            The function to evaluate
        x_low : float
            The lower boundary
        x_high : float
            The higher boundary
        function_args : tuple
            The function arguments
        max_iteration : int
            The maximum iteration
        tolerance : float
            Absolute root-finding tolerance.

        Returns
        -------
        float
            The root of the function
        """
        result = self.find_root_result_within_range(
            function=function,
            x_low=x_low,
            x_high=x_high,
            function_args=function_args,
            max_iteration=max_iteration,
            tolerance=tolerance,
        )
        return result.value if result.converged else numpy.nan

    def find_root_result_within_range(
            self,
            function,
            x_low: float,
            x_high: float,
            function_args: tuple = (),
            max_iteration: int = 100,
            tolerance: float = 1e-8) -> SolverResult[float]:
        """Find a root and retain convergence diagnostics.

        Parameters
        ----------
        function : callable
            Scalar characteristic equation.
        x_low, x_high : float
            Search interval.
        function_args : tuple, optional
            Additional characteristic-function arguments.
        max_iteration : int, optional
            Maximum Brent iterations.
        tolerance : float, optional
            Absolute convergence tolerance.

        Returns
        -------
        SolverResult[float]
            Root, residual, iteration count, bracket, and status.
        """

        boundaries = self.get_new_x_low_x_high(
            function=function,
            function_args=function_args,
            x_low=x_low,
            x_high=x_high,
            n_slice=100,
        )

        if numpy.isscalar(boundaries) and numpy.isnan(boundaries):
            message = f"could not bracket a root in [{x_low}, {x_high}]"
            logging.warning(message)
            return SolverResult(value=None, converged=False, message=message)

        x_low, x_high, y_low, y_high = boundaries

        x_root = root_scalar(
            method='brentq',
            x0=(x_low + x_high) / 2,
            f=function,
            bracket=[x_low, x_high],
            args=function_args,
            maxiter=max_iteration,
            xtol=tolerance,
            options=dict(disp=True)
        )

        residual = abs(float(function(x_root.root, *function_args)))
        return SolverResult(
            value=float(x_root.root) if x_root.converged else None,
            converged=bool(x_root.converged),
            residual=residual,
            iterations=x_root.iterations,
            bracket=(float(x_low), float(x_high)),
            message=x_root.flag,
        )

# -
