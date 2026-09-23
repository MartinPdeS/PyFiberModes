"""
Field-grid scaling
==================

Measure the cost of evaluating three vector components on increasing grids.
"""

from time import perf_counter

import matplotlib.pyplot as plt

from PyFiberModes import HE11
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1550e-9)
sizes = (51, 101, 151, 201)
elapsed = []

for size in sizes:
    field = fiber.get_mode_field(HE11, limit=12e-6, n_point=size)
    start = perf_counter()
    field.get_components(("Ex", "Ey", "Ez"))
    elapsed.append(perf_counter() - start)

figure, axis = plt.subplots()
axis.plot(sizes, elapsed, marker="o")
axis.set(xlabel="Grid width [points]", ylabel="Elapsed time [s]", title="Field evaluation")
plt.show()
