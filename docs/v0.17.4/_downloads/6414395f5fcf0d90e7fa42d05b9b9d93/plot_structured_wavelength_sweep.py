"""
Structured wavelength sweep
===========================

Sweep a mode without mutating the original fiber state.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1310e-9)
wavelengths = np.linspace(1260e-9, 1650e-9, 80)
result = fiber.sweep(
    wavelengths,
    modes=(LP01,),
    metrics=("effective_index", "group_index"),
)

figure, axis = plt.subplots()
axis.plot(wavelengths * 1e9, result["effective_index"][:, 0], label="effective index")
axis.plot(wavelengths * 1e9, result["group_index"][:, 0], label="group index")
axis.set(xlabel="Wavelength [nm]", ylabel="Index", title="SMF-28 LP01 sweep")
axis.legend()
plt.show()
