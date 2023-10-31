"""
Mode fields
===========
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import FiberFactory, HE11, HE12, HE22
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
    title='Mode fields for vectorial modse',
    unit_size=(4, 4),
    ax_orientation='horizontal'
)

for mode in [HE11, HE12, HE22]:
    fiber = factory[0]

    field = Field(
        fiber=fiber,
        mode=mode,
        wavelength=1550e-9,
        limit=30e-6,
        n_point=101
    )

    e_x = field.Ex()
    print(e_x)

    ax = figure.append_ax(title=mode)

    ax.add_mesh(
        scalar=e_x,
    )

figure.show()

# -
