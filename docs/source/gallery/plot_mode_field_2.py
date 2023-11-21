"""
Mode fields
===========
"""


# %%
# Imports
# ~~~~~~~
import numpy

from PyFiberModes import LP01, LP02, LP03
from PyFiberModes.fiber import load_fiber
from MPSPlots.render2D import SceneList

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
fiber = load_fiber(fiber_name='SMF28', wavelength=1550e-9)
fiber = fiber.scale(3.5)

# %%
# Preparing the figure
radius_list = numpy.linspace(0, 2 * fiber.radius, 200)


figure = SceneList()

ax = figure.append_ax(show_legend=True, line_width=2)

for mode in [LP01, LP02, LP03]:
    data = []
    for radius in radius_list:
        e_abs = fiber.get_radial_field_norm(mode=mode, radius=radius)

        data.append(e_abs)

    ax.add_line(x=radius_list, y=data, label=mode)


figure.show()

# -
