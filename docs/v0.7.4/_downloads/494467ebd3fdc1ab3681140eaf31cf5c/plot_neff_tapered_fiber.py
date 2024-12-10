"""
Effective index for tapered multi-mode fiber
============================================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from PyFiberModes.fiber import Fiber
from PyFiberModes import LP01, LP11, LP21, LP02
import matplotlib.pyplot as plt


# %%
# Computing the analytical values using FiberModes solver.
def get_mode_beta(fiber, mode_list: list, itr_list: list) -> dict:
    data_dict = {}
    for mode in mode_list:
        data_list = []
        for j, itr in enumerate(itr_list):
            _fiber = fiber.scale(factor=itr)
            data = _fiber.get_effective_index(mode=mode)
            data_list.append(data)

        data_dict[mode.__repr__()] = numpy.asarray(data_list)

    return data_dict


# %%
# Constructing fiber
fiber = Fiber(wavelength=1550e-9)
fiber.add_layer(name='core', radius=20e-6, index=1.4494)
fiber.add_layer(name='clad', radius=62.5e-6, index=1.4444)
fiber.initialize_layers()

itr_list = numpy.linspace(1, 0.1, 50)

data_dict = get_mode_beta(
    fiber=fiber,
    mode_list=[LP01, LP11, LP21, LP02],
    itr_list=itr_list
)

# %%
# Preparing the figure
figure, ax = plt.subplots(1, 1)

ax.set(
    title='Mode fields for vectorial mode if x-direction',
    xlabel='Inverse taper ratio',
    ylabel='Effective index'
)

for mode, data in data_dict.items():
    ax.plot(
        itr_list, data,
        label=mode,
        linestyle='-',
        linewidth=2
    )


plt.show()

# -
