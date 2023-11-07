from math import sqrt

from scipy.constants import c, h, mu_0, epsilon_0, eV, physical_constants, pi

tpi = 2 * pi
eta0 = physical_constants['characteristic impedance of vacuum'][0]
Y0 = sqrt(epsilon_0 / mu_0)  #: Admitance of free-space.
