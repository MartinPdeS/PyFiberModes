"""
Mode fields
===========
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import HE11
from PyFiberModes.fiber import load_fiber
from PyFiberModes.field import Field

# %%
# Loading the double clad fiber [DCF]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we load a fiber from MPSTools library and define the wavelength
fiber = load_fiber(fiber_name='SMF28', wavelength=1310e-9)

# %%
# Preparing the figure

field = Field(
    fiber=fiber,
    mode=HE11,
    limit=10e-6,
    n_point=50
)


figure = field.plot(plot_type=['Ex', 'Ey', 'Ez', 'Er', 'Ephi'])

figure.show()

# -
