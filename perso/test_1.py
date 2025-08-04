"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
import numpy
from PyFiberModes import LP01, LP11, LP21
from PyFiberModes.fiber import load_fiber
import matplotlib.pyplot as plt
from MPSPlots.styles import mps

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
wavelength_list = numpy.linspace(1310 - 350, 1550 + 50, 40) * 1e-9
with plt.style.context(mps):
    figure, axes = plt.subplots(2, 1)

ax = axes[0]
for mode in [LP01, LP11, LP21]:
    data = []
    for idx, wavelegnth in enumerate(wavelength_list):
        print(idx)
        smf28 = load_fiber(fiber_name='SMF28', wavelength=wavelegnth)
        smf28 = smf28.scale(1.2)
        dispersion = smf28.get_effective_index(mode=mode)
        data.append(dispersion)

    ax.plot(wavelength_list, data, label=mode)

ax.set(
    title='Mode effective index vs waveleght for FMF',
    xlabel=r'Wavelength [$m$]',
    ylabel='Effective index',
)

ax.legend()

ax = axes[1]
for mode in [LP01, LP11, LP21]:
    data = []
    for idx, wavelegnth in enumerate(wavelength_list):
        print(idx)
        smf28 = load_fiber(fiber_name='SMF28', wavelength=wavelegnth)
        smf28 = smf28.scale(1.2)
        dispersion = smf28.get_group_index(mode=mode)
        data.append(dispersion)

    ax.plot(wavelength_list, data, label=mode)

ax.set(
    title='Mode effective group index vs waveleght for FMF',
    xlabel=r'Wavelength [$m$]',
    ylabel='Group index',
)

ax.legend()

plt.show()
# -
