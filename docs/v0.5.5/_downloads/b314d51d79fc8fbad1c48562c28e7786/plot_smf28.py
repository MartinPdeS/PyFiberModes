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
smf28 = load_fiber(fiber_name='SMF28', wavelength=1310e-9)


smf28.print_data(
    data_type_list=['mode_cutoff_wavelength', 'effective_index', 'dispersion'], 
    mode_list=[HE11, LP01, TE01, LP11, LP02, LP21, LP12, HE22, LP31, LP22, LP41]
)


# -
