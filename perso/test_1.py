"""
Propagation constant: SMF28
===========================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from MPSPlots.render2D import SceneList
from PyFiberModes.fiber import load_fiber
from PyFiberModes import LP01, LP02, LP11, LP21


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


smf28 = load_fiber(fiber_name='SMF28', wavelength=1550e-9, add_air_layer=False)

itr_list = numpy.linspace(1.0, 0.05, 300)

data_dict = get_mode_beta(
    fiber=smf28,
    mode_list=[LP01, LP02, LP11, LP21],
    itr_list=itr_list
)

# %%
# Preparing the figure
figure = SceneList(unit_size=(12, 4))

ax = figure.append_ax(
    x_label='Inverse taper ratio',
    y_label='Effective index',
    show_legend=True,
    font_size=18,
    tick_size=15,
    legend_font_size=18
)


for mode, data in data_dict.items():
    ax.add_line(
        x=itr_list,
        y=data,
        label=mode,
        line_style='-',
        line_width=2,
        layer_position=1
    )


_ = figure.show()

# -
