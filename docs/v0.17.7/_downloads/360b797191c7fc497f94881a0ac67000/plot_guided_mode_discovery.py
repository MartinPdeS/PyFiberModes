"""
Guided-mode discovery
=====================

Discover the supported linearly polarized modes and compare their effective
indices at one wavelength.
"""

import matplotlib.pyplot as plt

from PyFiberModes.analysis import find_modes
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1310e-9)
modes = find_modes(fiber, families=("LP",), max_nu=3, max_m=3)
effective_indices = [fiber.get_effective_index(mode) for mode in modes]

figure, axis = plt.subplots()
axis.bar([str(mode) for mode in modes], effective_indices)
axis.set(xlabel="Guided mode", ylabel="Effective index", title="SMF-28 at 1310 nm")
axis.ticklabel_format(axis="y", useOffset=False)
plt.show()
