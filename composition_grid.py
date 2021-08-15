#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 2 12:36:05 2021

@author: hila
"""

import os
import numpy as np
import pandas as pd


def set_concentration_factor(x, the_wanted_lipid, new_con):
    # finds the new concentration factor for one wanted lipid
    if x['lipid'] == the_wanted_lipid:
        total_sum = np.sum(x["total_concentration"])
        factor = (total_sum-new_con*x["total_concentration"])/(total_sum-x["total_concentration"])
        return factor


def set_concentration_value(x, changed_lipid, factor, new_con):
    # calculates and returns the new concentration for the lipid
    if x["lipid"] == changed_lipid:
        new_concentration = new_con*x["total_concentration"]
    else:
        new_concentration = x["total_concentration"]*factor
    return new_concentration


old_grid = pd.read_csv("")
new_grid = "temp_until_i_find_the_function"
for i in [0, 1, 2, 3, 4]:
    wanted_lipid = ""
    f = set_concentration_factor(old_grid, wanted_lipid, i)
    new_grid["updated_concentration"] = old_grid.apply(set_concentration_value(old_grid,wanted_lipid, f, i),
                                                       axis="columns")
    new_grid.to_csv("???/new_composition_grid.csv")