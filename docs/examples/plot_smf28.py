"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import HE11, HE12, LP01, LP11, LP02, LP21
from PyFiberModes.fiber import load_fiber
import numpy

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
wavelength_list = numpy.linspace(1310e-9, 1550e-9, 10)

smf28 = load_fiber(fiber_name='DCF1300S_20', wavelength=1550e-9)

mode_field = smf28.get_mode_field(mode=HE12, n_point=20)

mode_field.plot(['Ex']).show()
# -
