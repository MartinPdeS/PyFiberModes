"""
Figure 3.13 of Jacques Bures
============================
"""


# %%
# Imports
# ~~~~~~~
import numpy

from PyFiberModes.fiber import get_fiber_from_delta_and_V0
from PyFiberModes.mode_instances import *
from PyFiberModes.fundamentals import get_U_parameter, get_cutoff_v0
from MPSPlots.render2D import SceneList


figure = SceneList(unit_size=(7, 5))
ax = figure.append_ax(
    x_label='V number',
    y_label='U number',
    show_legend=True
)

V0_list = numpy.linspace(0.5, 12, 50)
for V0 in V0_list:
    fiber = get_fiber_from_delta_and_V0(
        delta=0.3,
        V0=V0,
        wavelength=1310e-9
    )

    mode = HE11
    data = get_U_parameter(
        fiber=fiber,
        mode=mode,
        wavelength=fiber.wavelength
    )

    print('data', mode, data)