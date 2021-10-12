import numpy as np
import pandas as pd
import itertools
# import json


def set_concentration_factor(x, membrane_part, new_con, old_con):
    # finds the new concentration factor for one wanted lipid
    total_sum = np.sum(x[membrane_part])
    factor = (total_sum-((total_sum/100)*new_con))/(total_sum-old_con)
    return factor


def calc_old_con(x, changed_lipids, membrane_part):
    # calculates the old concentration for all wanted lipids in layer and returns their sum
    con = 0
    for itemb in changed_lipids:
        if (itemb[2] == membrane_part) or (itemb[2] == "total_concentration"):
            value = x.loc[x['ID'] == itemb[0], [membrane_part]]
            con = con + np.sum(value[membrane_part])
    return con


def set_concentration_value(x, changed_lipid, new_con, membrane_part, total_sum, factor):
    # calculates and returns the new concentration for the lipid
    if str(x['ID']) == changed_lipid[0][0] and (changed_lipid[0][2] == membrane_part or
                                                changed_lipid[0][2] == "total_concentration"):
        new_concentration = (total_sum/100)*changed_lipid[0][1][new_con[0]]
    elif str(x['ID']) == changed_lipid[1][0] and (changed_lipid[1][2] == membrane_part or
                                                  changed_lipid[1][2] == "total_concentration"):
        new_concentration = (total_sum / 100) * changed_lipid[1][1][new_con[1]]
    else:
        new_concentration = x[membrane_part]*factor
    return new_concentration


# set updated composition grid file
old_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
new_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
# membrane_parts = input("Choose membrane layer you want to change. \nFor both layers type 'total_concentration'."
#                       "\nFor outer layer type 'outside'. \nFor inner layer type 'inside':")
wanted_lipids = [("Cer4412", (3, 6, 1), "total_concentration"), ("Cer4202", (2, 3, 0.5), "inside")]
old_con_out = calc_old_con(old_grid, wanted_lipids, "outside")
old_con_in = calc_old_con(old_grid, wanted_lipids, "inside")
new_out = []
new_in = []
array_in = []
array_out = []
for lipid in wanted_lipids:
    if lipid[2] != "inside":
        new_out.append(np.arange(lipid[1][0], lipid[1][1], lipid[1][2]))
    if lipid[2] != "outside":
        new_in.append(np.arange(lipid[1][0], lipid[1][1], lipid[1][2]))
first_in, *rest_in = new_in
first_out, *rest_out = new_out
for tup in itertools.product(first_in, *rest_in):
    array_in.append(tup)
for tup in itertools.product(first_out, *rest_out):
    array_out.append(tup)
print(array_out, "-", array_in)
for lipid in range(len(wanted_lipids[:-1])):
    item = wanted_lipids[lipid]
    item_two = wanted_lipids[lipid+1]
    wanted_lipid = item[0]
    for i in range(len(item[1])):
        for z in range(len(item_two[1])):
            new_out = 0
            new_in = 0
            if item_two[2] != "inside":
                new_out += item_two[1][z]
            if item[2] != "inside":
                new_out += item[1][i]
            if item_two[2] != "outside":
                new_in += item_two[1][z]
            if item[2] != "outside":
                new_in += item[1][i]
            new_grid["outside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipids, new_con=[i, z],
                                                 membrane_part="outside",
                                                 total_sum=np.sum(old_grid["outside"]),
                                                 factor=set_concentration_factor(old_grid,
                                                                                 "outside",
                                                                                 new_out, old_con_out), axis="columns")
            new_grid["inside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipids, new_con=[i, z],
                                                membrane_part="inside",
                                                total_sum=np.sum(old_grid["inside"]),
                                                factor=set_concentration_factor(old_grid,
                                                                                "inside",
                                                                                new_in, old_con_in), axis="columns")
            new_grid["total_concentration"] = new_grid["inside"].fillna(value=0) + new_grid["outside"].fillna(value=0)
            new_grid.to_csv("C:\\Users\\hfComp\\Desktop\\WS\\new_composition_grid" + wanted_lipid + "_" +
                            str(item[1][i]) +
                            "_" + str(item_two[1][z]) + ".csv")
