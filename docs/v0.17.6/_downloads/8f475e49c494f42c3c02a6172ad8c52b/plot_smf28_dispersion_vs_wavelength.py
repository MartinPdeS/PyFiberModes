"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
import numpy
from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber
import matplotlib.pyplot as plt

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
wavelength_list = numpy.linspace(1000e-9, 1800e-9, 100)
data = []
for wavelegnth in wavelength_list:
    smf28 = load_fiber(fiber_name='SMF28', wavelength=wavelegnth)
    dispersion = smf28.get_effective_index(mode=LP01)
    data.append(dispersion)


figure, ax = plt.subplots(1, 1)

ax.set(
    title='Mode fields for vectorial mode if x-direction',
    xlabel=r'Wavelength [$\mu m$]',
    ylabel='Dispersion',
)


ax.plot(wavelength_list, data)

plt.show()
# -
