"""
Polarization rotation
=====================

Rotate the fundamental vector mode without recomputing its shared radial
solution.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes import HE11
from PyFiberModes.fiber import load_fiber

fiber = load_fiber("SMF28", wavelength=1550e-9)
field = fiber.get_mode_field(HE11, limit=12e-6, n_point=121)
angles = (0, np.pi / 4, np.pi / 2)

figure, axes = plt.subplots(1, 3, figsize=(9, 3), constrained_layout=True)
for axis, angle in zip(axes, angles):
    component = np.real(field.Ex(theta=angle))
    limit = np.max(np.abs(component))
    image = axis.imshow(component, origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit)
    axis.set_title(f"θ = {angle / np.pi:.2g}π")
    axis.set_axis_off()
figure.colorbar(image, ax=axes, shrink=0.75, label="Electric-field amplitude")
plt.show()
