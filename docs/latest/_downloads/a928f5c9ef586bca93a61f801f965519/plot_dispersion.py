"""
Dispersion VS core index
========================
"""
from PyFiberModes.units import ureg


# %%
# Imports
# ~~~~~~~
from PyFiberModes import ConvergenceError, FiberFactory, HE11, HE12, HE22
import matplotlib.pyplot as plt
import numpy

# %%
# Generating the fiber structures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we create the different fiber design that we want to explore
core_indexes = numpy.linspace(1.464, 1.494, 100)
factory = FiberFactory(wavelength=1550e-9 * ureg.meter)
factory.add_layer(name="core", radius=4e-6 * ureg.meter, index=(core_indexes) * ureg.dimensionless)
factory.add_layer(name="cladding", index=1.4444 * ureg.dimensionless)


# %%
# Preparing the figure
figure, ax = plt.subplots(1, 1)

ax.set(
    title='Dispersion vs core index',
    xlabel='Core refractive index',
    ylabel='Modal dispersion [ps/(nm km)]'
)

for mode in [HE11, HE12, HE22]:
    data = []
    for fiber in factory:
        try:
            effective_index = fiber.get_dispersion(mode=mode).to("ps / nm / km").magnitude
        except ConvergenceError:
            # Higher-order modes can be unguided below cutoff. Use NaN so the
            # plot shows a gap without obscuring other guided solutions.
            effective_index = numpy.nan
        data.append(effective_index)

    ax.plot(core_indexes, data, label=str(mode))

plt.legend()
plt.show()

# -
