from itertools import count
from scipy.optimize import brentq
import logging
import numpy


class FiberSolver(object):
    """
    Generic abstract class for callable objects used as fiber solvers.
    """

    logger = logging.getLogger(__name__)
    _MCD = 0.1

    def __init__(self, fiber, wavelength):
        self.fiber = fiber
        self.wavelength = wavelength

    def solver(self, *args, **kwargs):
        raise NotImplementedError()

    def find_function_first_root(self,
            function,
            function_args: tuple = (),
            lowbound: float = 0,
            highbound: float = None,
            ipoints: list = [],
            delta: float = 0.25,
            maxiter: int = numpy.inf):

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

    def find_root_within_range(self,
            function,
            lowbound: float,
            highbound: float,
            function_kwargs: tuple = (),
            max_iteration: int = 15):

        x_list = [lowbound, highbound]

        y_list = [
            function(lowbound, **function_kwargs),
            function(highbound, **function_kwargs)
        ]

        for j in count():  # probably not needed...
            if j == max_iteration:
                self.logger.warning("Couldn't converge to value as max iteration is reached")
                return numpy.nan

            for i in range(len(y_list) - 1):
                x0, x1 = x_list[i], x_list[i + 1]
                y0, y1 = y_list[i], y_list[i + 1]

                if (y0 > 0 and y1 < 0) or (y0 < 0 and y1 > 0):
                    args = tuple(function_kwargs.values())

                    z = brentq(
                        f=function,
                        a=x0,
                        b=x1,
                        args=args
                    )

                    y2 = function(z, **function_kwargs)

                    if abs(y0) > abs(y2) < abs(y1):  # Skip discontinuities
                        return z

            ls = len(y_list)
            for idx in range(ls - 1):
                a, b = x_list[2 * idx], x_list[2 * idx + 1]
                c = (a + b) / 2
                x_list.insert(2 * idx + 1, c)
                y_list.insert(2 * idx + 1, function(c, **function_kwargs))

# -
