#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFinitDiff.coefficients import FiniteCoefficients


def get_derivative(function, x: float, order: int, delta: float, function_kwargs: dict, accuracy: int = 4) -> float:
    """
    Returns the derivative a the given function

    :param      function:       The function
    :type       function:       function
    :param      x:              parameter to derive
    :type       x:              float
    :param      order:          Differentiation order (1 to 5)
    :type       order:          int
    :param      n_point:        Number of points (3 to 6)
    :type       n_point:        int
    :param      delta:          Distance between points
    :type       delta:          float
    :param      args:           The arguments
    :type       args:           list

    :returns:   The value of the derivative
    :rtype:     float
    """
    coefficients = FiniteCoefficients(
        derivative=order,
        accuracy=4,
        coefficient_type='central'
    )

    summation = 0
    for index, value in coefficients:
        x_eval = x + index * delta

        y = function(x_eval, **function_kwargs)

        summation += value * y

    return summation / delta**order


# -
