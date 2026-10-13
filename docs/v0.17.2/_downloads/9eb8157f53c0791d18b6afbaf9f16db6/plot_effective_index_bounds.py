"""
Effective-index bounds
======================

A guided mode must remain between the cladding and maximum material indices.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber

wavelengths = np.linspace(1260e-9, 1650e-9, 80)
fiber = load_fiber("SMF28", wavelength=wavelengths[0])
result = fiber.sweep(wavelengths, modes=(LP01,), metrics=("effective_index",))
neff = result["effective_index"][:, 0]

figure, axis = plt.subplots()
axis.plot(wavelengths * 1e9, neff, label=r"$n_\mathrm{eff}$")
axis.axhline(fiber.minimum_index, linestyle="--", label="minimum material index")
axis.axhline(fiber.maximum_index, linestyle="--", label="maximum material index")
axis.set(xlabel="Wavelength [nm]", ylabel="Refractive index")
axis.legend()
plt.show()
