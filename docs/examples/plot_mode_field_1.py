"""
Mode fields
===========
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import FiberFactory, HE11, HE12, HE22, LP01
from PyFiberModes.field import Field
from MPSPlots.render2D import SceneList

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
core_indexes = 1.54
factory = FiberFactory()
factory.addLayer(name="core", radius=4e-6, index=core_indexes)
factory.addLayer(name="cladding", index=1.4444)


# %%
# Preparing the figure
figure = SceneList(
    title='Mode fields for vectorial mode if x-direction',
    unit_size=(4, 4),
    ax_orientation='horizontal'
)

mode = HE11

fiber = factory[0]

field = Field(
    fiber=fiber,
    mode=mode,
    wavelength=1550e-9,
    limit=10e-6,
    n_point=201
)


for field_name in ['Ex', 'Ey', 'Ez', 'Er', 'Ephi']:
    field_array = getattr(field, field_name)()

    ax = figure.append_ax(title=field_name, show_colorbar=True)

    ax.add_mesh(scalar=field_array)


figure.show()

# -
