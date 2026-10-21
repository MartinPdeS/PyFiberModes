"""
Modal-cache speedup
===================

Compare the first effective-index solution with a repeated lookup served by
the geometry-aware modal cache.
"""

from time import perf_counter

import matplotlib.pyplot as plt

from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1550e-9)
start = perf_counter()
fiber.get_effective_index(LP01)
uncached_time = perf_counter() - start

start = perf_counter()
fiber.get_effective_index(LP01)
cached_time = perf_counter() - start

figure, axis = plt.subplots()
axis.bar(("First solution", "Cached lookup"), (uncached_time, cached_time))
axis.set(ylabel="Elapsed time [s]", title="Geometry-aware modal cache")
plt.show()
