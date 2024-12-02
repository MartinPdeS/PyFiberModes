"""
Mode fields
===========
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import HE11, LP01, LP11
from PyFiberModes.field import Field
from PyFiberModes.fiber import load_fiber
from MPSPlots.render2D import SceneList

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
fiber = load_fiber('SMF28', wavelength=1310e-9)

# %%
# Preparing the figure
figure = SceneList(
    title='Mode fields for vectorial mode if x-direction',
    unit_size=(4, 4),
    ax_orientation='horizontal',
)

for mode in [HE11, LP01, LP11]:

    field = Field(
        fiber=fiber,
        mode=mode,
        limit=10e-6,
        n_point=101
    )

    ax = figure.append_ax(title=mode)

    field.add_to_ax(field_string='Ex', ax=ax)

figure.show()

# -
