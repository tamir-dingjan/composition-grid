#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 2 12:36:05 2021

@author: hila
"""

import numpy as np
import pandas as pd
# import json


def set_concentration_factor(x, changed_lipid, membrane_part, new_con):
    # finds the new concentration factor for one wanted lipid
    old_con = x.loc[x['ID'] == changed_lipid, [membrane_part]]
    total_sum = np.sum(x[membrane_part])
    factor = (total_sum-((total_sum/100)*new_con))/(total_sum-old_con)
    return factor.at[87, membrane_part]


def set_concentration_value(x, changed_lipid, new_con, membrane_part, total_sum, factor):
    # calculates and returns the new concentration for the lipid
    if str(x['ID']) == changed_lipid:
        new_concentration = (total_sum/100)*new_con
    else:
        new_concentration = x[membrane_part]*factor
    return new_concentration


# set updated composition grid file
old_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
new_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
membrane_parts = input("Choose membrane layer you want to change. \nFor both layers type 'total_concentration'."
                       "\nFor outer layer type 'outside'. \nFor inner layer type 'inside':")
wanted_lipid = "Cer4412"
for i in [3, 4, 5]:
    new_grid[membrane_parts] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipid, new_con=i,
                                              membrane_part=membrane_parts,
                                              total_sum=np.sum(old_grid[membrane_parts]),
                                              factor=set_concentration_factor(old_grid, wanted_lipid, membrane_parts,
                                                                              i), axis="columns")
    if membrane_parts == "total_concentration":
        new_grid["outside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipid, new_con=i,
                                             membrane_part="outside",
                                             total_sum=np.sum(old_grid["outside"]),
                                             factor=set_concentration_factor(old_grid, wanted_lipid,
                                                                             "outside",
                                                                             i), axis="columns")
        new_grid["inside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipid, new_con=i,
                                            membrane_part="inside",
                                            total_sum=np.sum(old_grid["inside"]),
                                            factor=set_concentration_factor(old_grid, wanted_lipid,
                                                                            "inside",
                                                                            i), axis="columns")
    else:
        new_grid["total_concentration"] = new_grid["inside"].fillna(value=0) + new_grid["outside"].fillna(value=0)
    new_grid.to_csv("C:\\Users\\hfComp\\Desktop\\WS\\new_composition_grid" + str(i) + ".csv")
