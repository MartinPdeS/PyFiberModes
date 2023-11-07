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

    def __init__(self, fiber):
        self.fiber = fiber
        self.log = []
        self._logging = False

    def __call__(self, *args, **kwargs):
        raise NotImplementedError()

    def start_log(self):
        self.log = []
        self._logging = True

    def stop_log(self):
        self._logging = False

    def __record(self, fct):
        def wrapper(z, *args):
            r = fct(z, *args)
            if self._logging:
                self.log.append((z, r))
            return r
        return wrapper

    def find_function_first_root(self,
            function,
            function_args: tuple = (),
            lowbound: float = 0,
            highbound: float = None,
            ipoints: list = [],
            delta: float = 0.25,
            maxiter: int = numpy.inf):

        function = self.__record(function)  # For debug purpose.
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
                        return float("nan")
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
        return float("nan")

    def _findBetween(self,
            function,
            lowbound: float,
            highbound: float,
            function_args: tuple = (),
            max_iteration: int = 15):

        function = self.__record(function)  # For debug purpose.
        v = [lowbound, highbound]

        s = [
            function(lowbound, *function_args),
            function(highbound, *function_args)
        ]

        for j in count():  # probably not needed...
            if j == max_iteration:
                self.logger.warning("_findBetween: max iter reached")
                return float("nan")
            for i in range(len(s) - 1):
                a, b = v[i], v[i + 1]
                fa, fb = s[i], s[i + 1]

                if (fa > 0 and fb < 0) or (fa < 0 and fb > 0):
                    z = brentq(function, a, b, args=function_args)
                    fz = function(z, *function_args)
                    if abs(fa) > abs(fz) < abs(fb):  # Skip discontinuities
                        return z

            ls = len(s)
            for i in range(ls - 1):
                a, b = v[2 * i], v[2 * i + 1]
                c = (a + b) / 2
                v.insert(2 * i + 1, c)
                s.insert(2 * i + 1, function(c, *function_args))
