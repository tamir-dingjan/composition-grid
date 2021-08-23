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
    print(factor.at[87, 'total_concentration'])
    return factor.at[87, 'total_concentration']


def set_concentration_value(x, changed_lipid, new_con, factor):
    # calculates and returns the new concentration for the lipid
    if str(x['ID']) == changed_lipid:
        new_concentration = new_con*x["total_concentration"]
    else:
        new_concentration = x['total_concentration'] * factor
    return new_concentration


# set updated composition grid file
old_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
new_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
for i in [6,3,2]:
    wanted_lipid = "Cer4412"
    f = set_concentration_factor(old_grid, wanted_lipid, i)
    new_grid["updated_concentration"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipid, new_con=i,
                                                       factor=f, axis="columns")
    new_grid.to_csv("C:\\Users\\hfComp\\Desktop\\WS\\new_composition_grid" + str(i) + ".csv")
