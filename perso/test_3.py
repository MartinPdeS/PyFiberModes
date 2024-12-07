"""
Mode fields
===========
"""


# %%
# Imports
# ~~~~~~~
import numpy

from scipy.special import jn

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
    amplitudes = []
    for radius in radius_list:
        e_field, _ = fiber.get_radial_field(mode=mode, radius=radius)
        amplitudes.append(e_field.rho)

    amplitudes = numpy.asarray(amplitudes)
    amplitudes /= numpy.sign(amplitudes[0])
    ax.add_line(x=radius_list, y=amplitudes, label=mode)

    norm = fiber.get_mode_norm(mode=mode)


figure.show()

# -
