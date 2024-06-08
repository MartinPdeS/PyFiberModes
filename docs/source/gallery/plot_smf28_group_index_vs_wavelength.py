"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
import numpy
from PyFiberModes import HE11, HE22, LP01, LP11, LP02, LP21, LP12, TE01, LP22
from PyFiberModes.fiber import load_fiber
from MPSPlots.render2D import SceneList

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
wavelength_list = numpy.linspace(500e-9, 3000e-9, 100)
data = []
for wavelegnth in wavelength_list:
    smf28 = load_fiber(fiber_name='SMF28', wavelength=wavelegnth)
    dispersion = smf28.get_effective_index(mode=LP01)
    data.append(dispersion)


figure = SceneList()

ax = figure.append_ax(
    x_label=r'Wavelength [$\mu m$]',
    y_label='Group index',
    x_scale_factor=1e6
)

ax.add_line(x=wavelength_list, y=data)

_ = figure.show()
# -
