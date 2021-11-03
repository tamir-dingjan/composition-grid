import numpy as np
import pandas as pd
import itertools
import argparse


def set_concentration_factor(x, membrane_part, new_con, old_con):
    # finds the new concentration factor for one wanted lipid
    total_sum = np.sum(x[membrane_part])
    factor = (total_sum-((total_sum/100)*new_con))/(total_sum-old_con)
    return factor


def calc_old_con(x, changed_lipids, membrane_part):
    # calculates the old concentration for all wanted lipids in layer and returns their sum
    con, index_num = 0, 0
    if membrane_part == "inside":
        index_num = 1
    elif membrane_part == "outside":
        index_num = 2
    for itemb in changed_lipids:
        if itemb[index_num] != (0, 0, 1):
            value = x.loc[x['ID'] == itemb[0], [membrane_part]]
            con = con + np.sum(value[membrane_part])
    return con


def calc_new_con_grid(membrane_part):
    global wanted_lipids
    index_num = 0
    if membrane_part == "inside":
        index_num = 1
    elif membrane_part == "outside":
        index_num = 2
    new = []
    for lipid in wanted_lipids:
        lipid[index_num] = np.arange(lipid[index_num][0], lipid[index_num][1]+lipid[index_num][2], lipid[index_num][2])
        new.append(lipid[index_num])
    array = []
    first, *rest = new
    for tup in itertools.product(first, *rest):
        array.append(tup)
    return array


def process_config_file(wanted_lipids_file):
    wanted_lipids_list = wanted_lipids_file.values.tolist()
    working_lipids = []
    for item in wanted_lipids_list:
        item[1] = [float(s) for s in item[1].split()]
        item[2] = [float(s) for s in item[2].split()]
        if len(item[0].split()) > 1:
            for new_item in item[0].split():
                working_lipids.append([new_item, item[1], item[2]])
            wanted_lipids_list.remove(item)
    for item in working_lipids:
        wanted_lipids_list.append(item)
    return wanted_lipids_list


def set_concentration_value(x, changed_lipid, new_con, membrane_part, total_sum, factor):
    # calculates and returns the new concentration for the lipid
    if factor == 0:
        factor = 1
    index_layer = 0
    if membrane_part == "inside":
        index_layer = 1
    elif membrane_part == "outside":
        index_layer = 2
    ids = []
    for lipid in changed_lipid:
        ids.append(lipid[0])
    if str(x['ID']) in ids:
        index_num = ids.index(x['ID'])
        if new_con[index_num] != 0:
            new_concentration = (total_sum/100)*new_con[index_num]
        else:
            new_concentration = x[membrane_part] * factor
    else:
        new_concentration = x[membrane_part]*factor
    return new_concentration


# set updated composition grid file
old_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
wanted_lipid_file = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\config_file.csv")
# membrane_parts = input("Choose membrane layer you want to change. \nFor both layers type 'total_concentration'."
#                       "\nFor outer layer type 'outside'. \nFor inner layer type 'inside':")
# wanted_lipids = [[lipid_name, inside_tuple, outside_tuple], [...], [...]]
wanted_lipids = process_config_file(wanted_lipid_file)
# wanted_lipids = [["Cer4412", (3, 5, 1), (3, 5, 1)], ["Cer4202", (2, 2.5, 0.5), (0, 0, 1)],
#                 ["Cer4222", (1, 2, 1), (2, 4, 2)]]
old_con_out = calc_old_con(old_grid, wanted_lipids, "outside")
old_con_in = calc_old_con(old_grid, wanted_lipids, "inside")
array_out = calc_new_con_grid("outside")
array_in = calc_new_con_grid("inside")
print(array_out, "-", array_in)
print(wanted_lipids[1])
count = 0
for outer_con in array_out:
    trash = [array_in[0]]
    while array_in[0][1] == trash[-1][1]:
        count += 1
        inner_con = array_in[0]
        new_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
        new_grid["outside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipids,
                                             new_con=outer_con,
                                             membrane_part="outside",
                                             total_sum=np.sum(old_grid["outside"]),
                                             factor=set_concentration_factor(
                                                 old_grid, "outside", sum(outer_con), old_con_out),
                                             axis="columns")
        new_grid["inside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipids,
                                            new_con=inner_con,
                                            membrane_part="inside",
                                            total_sum=np.sum(old_grid["inside"]),
                                            factor=set_concentration_factor(old_grid,
                                                                            "inside",
                                                                            sum(inner_con), old_con_in),
                                            axis="columns")
        new_grid["total_concentration"] = new_grid["inside"].fillna(value=0) + new_grid["outside"].fillna(value=0)
        new_grid.to_csv("C:\\Users\\hfComp\\Desktop\\WS\\new_composition_grid" + wanted_lipids[0][0] + "_in" +
                        str(sum(inner_con)) +
                        "_out" + str(sum(outer_con)) + "_" + str(count) + ".csv")
        trash.append(array_in.pop(0))
        print(trash[-1])
        if len(array_in) == 0:
            break



"""
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
"""