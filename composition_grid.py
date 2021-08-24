#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 2 12:36:05 2021

@author: hila
"""

import numpy as np
import pandas as pd


def set_concentration_factor(x, changed_lipid, new_con):
    # finds the new concentration factor for one wanted lipid
    old_con = x.loc[x['ID'] == changed_lipid, ['total_concentration']]
    total_sum = np.sum(x["total_concentration"])
    factor = (total_sum-new_con*old_con)/(total_sum-old_con)
    return factor.at[87, 'total_concentration']


def set_concentration_value(x, changed_lipid, new_con, membrane_part, factor):
    # calculates and returns the new concentration for the lipid
    if str(x['ID']) == changed_lipid:
        new_concentration = new_con*x[membrane_part]
    else:
        new_concentration = x[membrane_part] * factor
    return new_concentration


# set updated composition grid file
old_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
new_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
membrane_parts = input("Choose membrane layer you want to change. \nFor both layers type 'total_concentration'."
                       "\nFor outer layer type 'outside'. \nFor inner layer type 'inside':")
wanted_lipid = input("Insert lipid ID:")
for i in [2, 3, 4]:
    f = set_concentration_factor(old_grid, wanted_lipid, i)
    new_grid[membrane_parts] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipid, new_con=i,
                                              membrane_part=membrane_parts, factor=f, axis="columns")
    new_grid.to_csv("C:\\Users\\hfComp\\Desktop\\WS\\new_composition_grid" + str(i) + ".csv")
