"""Shared test helpers that are intentionally not part of the public package."""

import numpy


def get_mode_beta(fiber, mode_list: list, itr_list: list) -> dict[str, numpy.ndarray]:
    """Evaluate effective indices while scaling a fiber through ``itr_list``."""
    return {
        repr(mode): numpy.asarray(
            [fiber.scale(factor=itr).get_effective_index(mode=mode) for itr in itr_list]
        )
        for mode in mode_list
    }
