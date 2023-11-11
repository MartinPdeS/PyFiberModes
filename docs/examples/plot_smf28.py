"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import HE11, HE12, LP01, LP11, LP02, LP21, LP12
from PyFiberModes.fiber import load_fiber

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
smf28 = load_fiber(fiber_name='SMF28', wavelength=1310e-9)


smf28.print_data(data_type_list=['cutoff_V0', 'cutoff_wavelength'], mode_list=[LP01, LP11, LP21, LP12])


# -
