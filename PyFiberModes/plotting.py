"""Optional Matplotlib visualizations for fiber-mode results."""
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from PyFiberModes.field import Field


def plot_field(
    field: "Field",
    plot_type: list[str] | tuple[str, ...] = (),
    show: bool = True,
    save_filename: str | None = None,
) -> "Figure":
    """Plot selected components of a computed modal field.

    Parameters
    ----------
    field : Field
        Modal field whose components will be rendered.
    plot_type : list of str or tuple of str, optional
        Names of the field components to plot.
    show : bool, optional
        Display the figure interactively after rendering.
    save_filename : str, optional
        Destination at which to save the rendered figure.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing one axis per requested component.

    Raises
    ------
    ImportError
        If Matplotlib is not installed.
    ValueError
        If no field components are requested.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError as error:
        raise ImportError(
            "Field plotting requires the optional plotting dependencies. "
            "Install them with `python -m pip install 'PyFiberModes[plotting]'`."
        ) from error

    if not plot_type:
        raise ValueError("plot_type must contain at least one field component")

    color_map = LinearSegmentedColormap.from_list(
        "blue_black_red",
        ("#2455a4", "#000000", "#d73027"),
    )
    with plt.style.context("default"):
        figure, axes = plt.subplots(1, len(plot_type))

    components = field.get_components(components=plot_type)
    for axis, component_name in zip(np.atleast_1d(axes), plot_type):
        axis.set_aspect("equal")
        values = np.real_if_close(components[component_name])
        if np.iscomplexobj(values):
            values = np.abs(values)
        maximum = float(np.max(np.abs(values)))
        axis.pcolormesh(values, vmin=-maximum, vmax=maximum, cmap=color_map)

    figure.tight_layout()
    if save_filename:
        figure.savefig(save_filename)
    if show:
        plt.show()
    return figure
