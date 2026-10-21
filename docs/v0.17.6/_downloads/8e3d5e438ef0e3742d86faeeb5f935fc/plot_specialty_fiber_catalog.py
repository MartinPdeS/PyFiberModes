"""
Specialty-fiber catalog
=======================

Compare the layer indices of several bundled dispersion-compensating fibers.
"""

import matplotlib.pyplot as plt
import numpy as np

from PyFiberModes.fiber import load_fiber

fiber_names = ("DCF1300S_20", "DCF1300S_26", "DCF1300S_33", "DCF1300S_42")
fibers = [load_fiber(name, wavelength=1550e-9) for name in fiber_names]
layer_indices = np.asarray(
    [[layer.refractive_index for layer in fiber.layers] for fiber in fibers]
)

figure, axis = plt.subplots()
positions = np.arange(len(fiber_names))
width = 0.24
for layer_number in range(layer_indices.shape[1]):
    axis.bar(
        positions + (layer_number - 1) * width,
        layer_indices[:, layer_number],
        width=width,
        label=f"Layer {layer_number + 1}",
    )
axis.set_xticks(positions, fiber_names, rotation=20)
axis.set(ylabel="Refractive index", title="Bundled specialty-fiber designs")
axis.ticklabel_format(axis="y", useOffset=False)
axis.legend()
plt.show()
