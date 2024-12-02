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
from MPSPlots.render2D import SceneList

# %%
# Loading the double clad fiber [DCF]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we load a fiber from MPSTools library and define the wavelength
fiber = load_fiber(fiber_name='DCF1300S_20', wavelength=1310e-9)

# %%
# Preparing the figure
figure = SceneList(
    title='Mode fields for vectorial mode if x-direction',
    unit_size=(4, 4),
    ax_orientation='horizontal'
)

field = Field(
    fiber=fiber,
    mode=HE11,
    limit=10e-6,
    n_point=30
)


figure = field.plot(plot_type=['Ex', 'Ey', 'Ez', 'Er', 'Ephi'])

figure.show()

# -
