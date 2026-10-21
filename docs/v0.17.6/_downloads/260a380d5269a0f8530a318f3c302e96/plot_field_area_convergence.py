"""
Effective-area convergence
==========================

Repeat a derived field calculation on increasingly fine grids to check that
the reported effective area has stabilized.
"""

import matplotlib.pyplot as plt

from PyFiberModes import HE11
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1550e-9)
grid_sizes = (41, 61, 81, 101, 121)
areas = []

for grid_size in grid_sizes:
    field = fiber.get_mode_field(HE11, limit=15e-6, n_point=grid_size)
    areas.append(field.get_effective_area() * 1e12)

figure, axis = plt.subplots()
axis.plot(grid_sizes, areas, marker="o")
axis.set(xlabel="Grid width [points]", ylabel="Effective area [µm²]")
axis.grid(alpha=0.25)
plt.show()
