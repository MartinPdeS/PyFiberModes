"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import FiberFactory, HE11, HE12, LP01, LP11, LP02
from PyFiberModes.fiber import load_fiber
from MPSPlots.render2D import SceneList
import numpy

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
wavelength_list = numpy.linspace(1310e-9, 1550e-9, 10)

smf28 = load_fiber(fiber_name='DCF1300S_20', wavelength=1550e-9)

effective_index = smf28.get_effective_index(
    mode=LP11,
    wavelength=1550e-9
)

print(effective_index)
# -
