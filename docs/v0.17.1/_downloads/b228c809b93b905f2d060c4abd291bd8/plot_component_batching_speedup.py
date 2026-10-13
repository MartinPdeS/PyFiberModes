"""
Component batching
==================

Compare separate field-component requests with the batched API that shares
radial and angular intermediates.
"""

from time import perf_counter

import matplotlib.pyplot as plt

from PyFiberModes import HE11
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1550e-9)
components = ("Ex", "Ey", "Ez")

separate_field = fiber.get_mode_field(HE11, limit=12e-6, n_point=101)
start = perf_counter()
for component in components:
    getattr(separate_field, component)()
separate_time = perf_counter() - start

batched_field = fiber.get_mode_field(HE11, limit=12e-6, n_point=101)
start = perf_counter()
batched_field.get_components(components)
batched_time = perf_counter() - start

figure, axis = plt.subplots()
axis.bar(("Separate calls", "Batched request"), (separate_time, batched_time))
axis.set(ylabel="Elapsed time [s]", title="Three electric-field components")
plt.show()
