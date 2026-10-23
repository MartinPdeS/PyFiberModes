"""
Specialty-fiber mode indices
============================

Compare the effective indices of low-order modes in a multilayer fiber.
"""

import matplotlib.pyplot as plt

from PyFiberModes import LP01, LP11, LP21
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("DCF1300S_20", wavelength=1550e-9)
modes = (LP01, LP11, LP21)
effective_indices = [fiber.get_effective_index(mode) for mode in modes]

figure, axis = plt.subplots()
axis.bar([str(mode) for mode in modes], effective_indices)
axis.set(
    xlabel="Mode",
    ylabel="Effective index",
    title="Low-order modes of DCF1300S-20",
)
axis.ticklabel_format(axis="y", useOffset=False)
plt.show()
