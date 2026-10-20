"""
Normalized effective index
==========================

Verify that the normalized propagation parameter remains between zero and one
for a guided fundamental mode.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber

wavelengths = np.linspace(1260e-9, 1650e-9, 60)
fiber = load_fiber("SMF28", wavelength=wavelengths[0])
result = fiber.sweep(wavelengths, modes=(LP01,), metrics=("effective_index",))
effective_index = result["effective_index"][:, 0]
cladding_index = fiber.last_layer.refractive_index
core_index = fiber.maximum_index
normalized_index = (
    (effective_index**2 - cladding_index**2)
    / (core_index**2 - cladding_index**2)
)

figure, axis = plt.subplots()
axis.plot(wavelengths * 1e9, normalized_index)
axis.axhspan(0, 1, alpha=0.12, label="guided interval")
axis.set(xlabel="Wavelength [nm]", ylabel="Normalized effective index", ylim=(-0.05, 1.05))
axis.legend()
plt.show()
