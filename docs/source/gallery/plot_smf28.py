"""
Modal dispersion VS core index
==============================
"""


# %%
# Imports
# ~~~~~~~
from PyFiberModes import HE11, HE22, LP01, LP11, LP02, LP21, LP12, TE01, LP31, LP22, LP41
from PyFiberModes.fiber import load_fiber

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
smf28 = load_fiber(fiber_name='DCF1300S_20', wavelength=1310e-9)

a = smf28.does_mode_exist(LP01, LP11, )
# print(a)
smf28.print_data(data_type_list=['cutoff_wavelength'], mode_list=[LP01, HE11, HE22, TE01])


# -
