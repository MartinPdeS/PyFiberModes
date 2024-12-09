"""
Mode fields
===========
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import FiberFactory, HE11
from PyFiberModes.field import Field

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
core_indexes = 1.54
factory = FiberFactory(wavelength=1550e-9)
factory.add_layer(name="core", radius=4e-6, index=core_indexes)
factory.add_layer(name="cladding", index=1.4444)

mode = HE11

fiber = factory[0]

field = Field(
    fiber=fiber,
    mode=mode,
    limit=10e-6,
    n_point=201
)


figure = field.plot(plot_type=['Ex', 'Ey', 'Ez', 'Er', 'Ephi'])


# -
