"""
Mode fields
===========
"""
from PyFiberModes.units import ureg


# %%
# Imports
# ~~~~~~~
from PyFiberModes import FiberFactory, LP02
from PyFiberModes.field import Field

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
core_indexes = 1.54
factory = FiberFactory(wavelength=1550e-9 * ureg.meter)
factory.add_layer(name="core", radius=4e-6 * ureg.meter, index=(core_indexes) * ureg.dimensionless)
factory.add_layer(name="cladding", index=1.4444 * ureg.dimensionless)

mode = LP02

fiber = factory[0]

field = Field(
    fiber=fiber,
    mode=mode,
    limit=10e-6 * ureg.meter,
    n_point=201
)

figure = field.plot(plot_type=['Ex', 'Ey'])
# figure = field.plot(plot_type=['Ex', 'Ey', 'Ez', 'Er', 'Ephi'])


# -
