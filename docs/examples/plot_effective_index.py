"""
Effective index VS core index
=============================
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import FiberFactory, HE11, HE12, HE22
from MPSPlots.render2D import SceneList
import numpy

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
core_indexes = numpy.linspace(1.464, 1.494)
factory = FiberFactory()
factory.addLayer(name="core", radius=4e-6, index=core_indexes)
factory.addLayer(name="cladding", index=1.4444)


# %%
# Preparing the figure
figure = SceneList(title='Effective index vs core index')

ax = figure.append_ax(show_legend=True)

for mode in [HE11, HE12, HE22]:
    neff = []
    for fiber in factory:
        effective_index = fiber.neff(mode, 1550e-9)
        neff.append(effective_index)

    ax.add_line(
        x=core_indexes,
        y=neff,
        label=mode,
        line_width=2
    )

figure.show()

# -
