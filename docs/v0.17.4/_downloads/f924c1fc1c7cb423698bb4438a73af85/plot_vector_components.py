"""
Vector field components
=======================

Evaluate several components together so they share one radial solver pass.
"""

import matplotlib.pyplot as plt

from PyFiberModes import HE11
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1550e-9)
field = fiber.get_mode_field(HE11, limit=12e-6, n_point=151)
components = field.get_components(("Ex", "Ey", "Ez"))

figure, axes = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
for axis, (name, values) in zip(axes, components.items()):
    image = axis.imshow(values.real, origin="lower", cmap="RdBu_r")
    axis.set_title(name)
    axis.set_axis_off()
    figure.colorbar(image, ax=axis, shrink=0.75)

plt.show()
