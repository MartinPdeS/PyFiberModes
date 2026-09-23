"""
Refractive-index profile
========================

Inspect the concentric refractive-index steps of a specialty fiber before
solving its modes.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes.fiber import load_fiber

fiber = load_fiber("DCF1300S_20", wavelength=1550e-9)
radii = np.linspace(0, 70e-6, 600)
indices = fiber.get_index_profile(radii)

figure, axis = plt.subplots()
axis.step(radii * 1e6, indices, where="post")
axis.set(
    xlabel="Radius [µm]",
    ylabel="Refractive index",
    title="DCF1300S-20 radial material profile",
)
axis.grid(alpha=0.25)
plt.show()
