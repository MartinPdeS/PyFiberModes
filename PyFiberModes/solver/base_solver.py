# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import logging
import numpy

from scipy.optimize import root, brentq, bisect, brenth


class BaseSolver(object):
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

    def _find_root_within_range(self,
            function,
            x_low: float,
            x_high: float,
            function_args: tuple = (),
            max_iteration: int = 20) -> float:
        """
        Finds a root within range.

        :param      function:       The function
        :type       function:       object
        :param      x_low:          The x low
        :type       x_low:          float
        :param      x_high:         The x high
        :type       x_high:         float
        :param      function_args:  The function arguments
        :type       function_args:  tuple
        :param      max_iteration:  The maximum iteration
        :type       max_iteration:  int

        :returns:   The root value: x such as f(x) = 0
        :rtype:     float
        """
        return brenth(f=function, a=x_low, b=x_high, args=function_args)
        try:
            return brenth(f=function, a=x_low, b=x_high, args=function_args)
        except ValueError:

            delta_x = abs(x_high - x_low) / 2
            print('DSADS')
            return self.find_root_within_range(
                function=function,
                x_low=x_low + delta_x,
                x_high=x_high,
                function_args=function_args
            )

        opt = root(fun=function, x0=x_low, args=function_args)

        if not (x_low < opt.x[0] < x_high):
            self.logger.warning(f"Root found : {opt.x[0]} but outside of range: [{x_low:.16f}, {x_high:.16f}]")

        return opt.x[0]

    def get_new_x_low_x_high(self,
            function,
            function_args,
            x_low: float,
            x_high: float,
            n_iteration: int = 10):

        x_list = [x_low, x_high]
        x_list.sort()
        x_list = numpy.linspace(*x_list, n_iteration)
        y_list = [function(x, *function_args) for x in x_list]
        y_list = numpy.asarray(y_list)

        non_nan_idx = ~numpy.isnan(y_list)

        x_list = x_list[non_nan_idx]
        y_list = y_list[non_nan_idx]

        if len(y_list) < 2:
            return self.get_new_x_low_x_high(function, function_args, x_low, x_high, n_iteration=2 * n_iteration)

        sign_change = (numpy.diff(numpy.sign(y_list)) != 0) * 1

        if (sign_change == 0).all():
            return self.get_new_x_low_x_high(function, function_args, x_low, x_high, n_iteration=2 * n_iteration)

        sign_change_idx = numpy.where(sign_change == 1)[0]

        sign_change_idx = sign_change_idx[0]

        x_low, x_high = x_list[sign_change_idx], x_list[sign_change_idx + 1]

        y_low, y_high = y_list[sign_change_idx], y_list[sign_change_idx + 1]

        return x_low, x_high, y_low, y_high

    def find_root_within_range(self,
            function,
            x_low: float,
            x_high: float,
            function_args: tuple = (),
            max_iteration: int = 20):

        y_low, y_high = function(x_low, *function_args), function(x_high, *function_args)

        x_low, x_high, y_low, y_high = self.get_new_x_low_x_high(function, function_args, x_low, x_high, 100)

        try:
            x_root = brentq(f=function, a=x_low, b=x_high, args=function_args, maxiter=max_iteration, disp=True)  # Get x such as f(x) = 0
        except RuntimeError:
            self.logger.warning("Couldn't converge to value as max iteration is reached")
            return numpy.nan

        y_root = function(x_root, *function_args)  # f(x)

        if abs(y_low) > abs(y_root) < abs(y_high):  # Skip discontinuities
            return x_root

    def update_root_range(self, function, function_args: tuple, x_low: float, x_high: float, y_low: float, y_high: float) -> tuple:
        """
        Calculate the new x-range such that it updates the new x-y minimum or new x-y maximum.
        The new evalution is midway between x_low and x_high.

        :param      function:       The function
        :type       function:       object
        :param      function_args:  The function arguments
        :type       function_args:  tuple
        :param      x_low:          The x lower boundary
        :type       x_low:          float
        :param      x_high:         The x upper boundary
        :type       x_high:         float
        :param      y_low:          The y lower boundary
        :type       y_low:          float
        :param      y_high:         The y upper boundary
        :type       y_high:         float

        :returns:   The new x_low, x_high, y_low, y_high
        :rtype:     tuple
        """
        x_mid = (x_low + x_high) / 2
        y_mid = function(x_mid, *function_args)

        if y_mid > 0:
            y_high = y_mid
            x_high = x_mid

        else:
            y_low = y_mid
            x_low = x_mid

        return x_low, x_high, y_low, y_high


# -
