"""
Comparing two and three layers
==============================
"""


# %%
# Imports
# ~~~~~~~
import numpy

from PyFiberModes.fiber import load_fiber
from PyFiberModes import LP01
import matplotlib.pyplot as plt

fiber_2L = load_fiber(
    fiber_name='SMF28',
    wavelength=1550e-9,
    add_air_layer=False
)

fiber_3L = load_fiber(
    fiber_name='SMF28',
    wavelength=1550e-9,
    add_air_layer=True
)


data_2L = []
data_3L = []

itr_list = numpy.linspace(1.0, 0.3, 50)

for j, itr in enumerate(itr_list):
    _fiber_2L = fiber_2L.scale(factor=itr)
    _fiber_3L = fiber_3L.scale(factor=itr)

    neff_2L = _fiber_2L.get_effective_index(mode=LP01)
    neff_3L = _fiber_3L.get_effective_index(mode=LP01)

    data_2L.append(neff_2L)
    data_3L.append(neff_3L)


figure, ax = plt.subplots(1, 1, figsize=(10, 6))

ax.set_xlabel('Inverse taper ratio [ITR]')
ax.set_ylabel('LP01 effective index')
ax.plot(itr_list, data_2L, label='2 layer', linewidth=3)
ax.plot(itr_list, data_3L, label='3 layer', marker='o', markersize=5)

plt.show()

# -
