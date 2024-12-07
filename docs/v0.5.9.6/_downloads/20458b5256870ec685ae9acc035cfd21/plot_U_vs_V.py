"""
Figure 3.13 of Jacques Bures
============================
"""


# %%
# Imports
# ~~~~~~~
import numpy

from PyFiberModes.fiber import get_fiber_from_delta_and_V0
from PyFiberModes import HE11, TE01, TM01, HE21, EH11, HE31, HE12, HE22, HE32
from PyFiberModes.fundamentals import get_U_parameter
import matplotlib.pyplot as plt

figure, ax = plt.subplots(1, 1)

ax.set(
    title='Mode fields for vectorial mode if x-direction',
    xlabel='V number',
    ylabel='U number',
)

V0_list = numpy.linspace(0.5, 12, 150)

for mode in [HE11, TE01, TM01, HE21, EH11, HE31, HE12, HE22, HE32]:
    data_list = []
    for V0 in V0_list:
        fiber = get_fiber_from_delta_and_V0(
            delta=0.3,
            V0=V0,
            wavelength=1310e-9
        )

        data = get_U_parameter(
            fiber=fiber,
            mode=mode,
            wavelength=fiber.wavelength
        )

        data_list.append(data)

    ax.plot(V0_list, data_list, label=str(mode), linewidth=1.5)

ax.plot(V0_list, V0_list, linewidth=2, label='U = V')

_ = figure.show()
