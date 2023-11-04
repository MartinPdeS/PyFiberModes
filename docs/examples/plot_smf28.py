"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import FiberFactory, HE11, HE12, LP01
from PyFiberModes.fiber import load_fiber
from MPSPlots.render2D import SceneList
import numpy

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
wavelength_list = numpy.linspace(1310e-9, 1550e-9, 10)



core_indexes = numpy.linspace(1.464, 1.494, 100)
factory = FiberFactory()
factory.add_layer(name="core", radius=4e-6, index=core_indexes)
factory.add_layer(name="cladding lol", index=1.4464, radius=20e-6)
factory.add_layer(name="yo", index=1.4444, radius=32.5e-6)
fiber_test = factory[0]

print(fiber_test)



data = []

smf28 = load_fiber(fiber_name='SMF28', wavelength=1550e-9)

print('\n\n')
print(smf28)

effective_index = smf28.get_effective_index(
    mode=LP01,
    wavelength=1550e-9
)

print(effective_index)

# for wavelength in wavelength_list:
#     smf28 = load_fiber(fiber_name='SMF28', wavelength=wavelength)

#     effective_index = smf28.get_effective_index(
#         mode=LP01,
#         wavelength=wavelength
#     )

#     data.append(effective_index)

# %%
# Preparing the figure
figure = SceneList(title='Modal dispersion vs core index')

ax = figure.append_ax(show_legend=True)
# ax.add_line(
#     # x=core_indexes,
#     y=data,
#     # label=mode,
#     line_width=2
# )

# figure.show()

# -
