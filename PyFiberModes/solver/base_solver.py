# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import logging
import numpy

from scipy.optimize import brentq, root_scalar


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
            maxiter: int = numpy.inf) -> float:

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

    def get_new_x_low_x_high(self,
            function,
            function_args,
            x_low: float,
            x_high: float,
            n_slice: int = 100) -> tuple:
        """
        Gets the new x boundaries.
        Returns numpy.nan if no sign inversion found.

        :param      function:       The function
        :type       function:       { type_description }
        :param      function_args:  The function arguments
        :type       function_args:  { type_description }
        :param      x_low:          The x low
        :type       x_low:          float
        :param      x_high:         The x high
        :type       x_high:         float
        :param      n_slice:        The n iteration
        :type       n_slice:        int

        :returns:   The new x low x high.
        :rtype:     tuple
        """
        x_list = [x_low, x_high]
        x_list.sort()
        x_list = numpy.linspace(*x_list, n_slice)
        y_list = [function(x, *function_args) for x in x_list]
        y_list = numpy.asarray(y_list)

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

    def find_root_within_range(self,
            function,
            x_low: float,
            x_high: float,
            function_args: tuple = (),
            max_iteration: int = 100,
            tolerance: float = 1e-8) -> float:
        """
        Finds and return the root of a given function within range.

        :param      function:       The function to evaluate
        :type       function:       object
        :param      x_low:          The lower boundary
        :type       x_low:          float
        :param      x_high:         The higher boundary
        :type       x_high:         float
        :param      function_args:  The function arguments
        :type       function_args:  tuple
        :param      max_iteration:  The maximum iteration
        :type       max_iteration:  int

        :returns:   The root of the function
        :rtype:     float
        """
        y_low, y_high = function(x_low, *function_args), function(x_high, *function_args)

        boundaries = self.get_new_x_low_x_high(
            function=function,
            function_args=function_args,
            x_low=x_low,
            x_high=x_high,
            n_slice=100,
        )

        if numpy.isscalar(boundaries) and numpy.isnan(boundaries):
            logging.warning(f"Couldn't find neff root in range:[{x_low}, {x_high}]  for mode")
            return numpy.nan

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

        return x_root.root

# -
