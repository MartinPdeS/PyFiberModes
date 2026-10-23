"""
Structured wavelength sweep
===========================

Sweep a mode without mutating the original fiber state.
"""
from PyFiberModes.units import ureg

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1310e-9 * ureg.meter)
wavelengths = np.linspace(1260e-9, 1650e-9, 80) * ureg.meter
result = fiber.sweep(
    wavelengths,
    modes=(LP01,),
    metrics=("effective_index", "group_index"),
)

figure, axis = plt.subplots()
axis.plot(wavelengths.to("nanometer").magnitude, result["effective_index"][:, 0], label="effective index")
axis.plot(wavelengths.to("nanometer").magnitude, result["group_index"][:, 0], label="group index")
axis.set(xlabel="Wavelength [nm]", ylabel="Index", title="SMF-28 LP01 sweep")
axis.legend()
plt.show()
