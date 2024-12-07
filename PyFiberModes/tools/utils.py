#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


def get_mode_beta(fiber, mode_list: list, itr_list: list):
    data_dict = {}
    for mode in mode_list:
        data_list = []
        for j, itr in enumerate(itr_list):
            _fiber = fiber.scale(factor=itr)
            neff = _fiber.get_effective_index(mode=mode)
            data_list.append(neff)

        data_dict[mode.__repr__()] = numpy.asarray(data_list)

    return data_dict

# -
