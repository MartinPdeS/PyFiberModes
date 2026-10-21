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
import matplotlib.pyplot as plt

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
fiber = load_fiber('SMF28', wavelength=1310e-9)

# %%
# Preparing the figure
figure, ax = plt.subplots(1, 1)

ax.set(
    title='Mode fields for vectorial mode if x-direction',
    xlabel='Radius',
    ylabel='Amplitude'
)

for mode in [HE11, LP01, LP11]:

    field = Field(
        fiber=fiber,
        mode=mode,
        limit=10e-6,
        n_point=101
    )

    Ex = field.Ex()
    max_abs = abs(max(Ex.max(), Ex.min()))
    ax.pcolormesh(Ex, vmin=-max_abs, vmax=max_abs)


plt.show()

# -
