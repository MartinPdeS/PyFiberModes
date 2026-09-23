"""
Confinement by radius
=====================

Measure how the sampled electric-field energy accumulates away from the fiber
axis.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes import HE11
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1550e-9)
field = fiber.get_mode_field(HE11, limit=15e-6, n_point=151)
radii = np.linspace(0.5e-6, 14e-6, 80)
confinement = [field.get_confinement_factor(radius) for radius in radii]

figure, axis = plt.subplots()
axis.plot(radii * 1e6, confinement)
axis.set(
    xlabel="Integration radius [µm]",
    ylabel="Electric-energy fraction",
    ylim=(0, 1.02),
)
axis.grid(alpha=0.25)
plt.show()
