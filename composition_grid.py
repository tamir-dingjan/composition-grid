import numpy as np
import pandas as pd
import itertools
import argparse


def set_concentration_factor(x, membrane_part, new_con, old_con):
    # finds the new concentration factor for one wanted lipid
    total_sum = np.sum(x[membrane_part])
    factor = (total_sum-((total_sum/100)*new_con))/(total_sum-old_con)
    print(factor, membrane_part)
    return factor


def calc_old_con(x, changed_lipids, membrane_part):
    # calculates the old concentration for all wanted lipids in layer and returns their sum
    con, index_num = 0, 0
    if membrane_part == "inside":
        index_num = 1
    elif membrane_part == "outside":
        index_num = 2
    zero_index_list = []
    for item_lipid in changed_lipids:
        for item_name in item_lipid[0]:
            if item_lipid[index_num] != [0.0, 0.0, 1.]:
                value = x.loc[x['ID'] == item_name, [membrane_part]]
                con = con + np.sum(value[membrane_part])
                zero_index_list.append(changed_lipids.index(item_lipid))
    return con, zero_index_list


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
    global old_grid
    wanted_lipids_list = wanted_lipids_file.values.tolist()
    working_lipids = []
    for item in wanted_lipids_list:
        item[1] = [float(s) for s in item[1].split()]
        item[0] = item[0].split()
        if item[3] == "manual":
            item[2] = [float(s) for s in item[2].split()]
            if item[2] == [] or item[1] == []:
                raise Exception("Please specify wanted concentration variables for both leaflets :)")
        elif item[3] == "fixed":
            print(pd.isna(item[2]), pd.isna(item[1]))
            if (not all(pd.isna(item[1]))) and (not all(pd.isna(item[2]))):
                raise Exception("Please specify wanted concentration variables only for the chosen leaflet :)")
            out_grid = 0.0
            in_grid = 0.0
            for lipid in item[0]:
                out_grid = out_grid + ((old_grid.loc[old_grid["ID"] == str(lipid),
                                                     ["outside"]].values[0][0]) / np.sum(old_grid["outside"])) * 100
                print(out_grid, type(out_grid))
                in_grid = in_grid + ((old_grid.loc[old_grid["ID"] == str(lipid),
                                                   ["inside"]].values[0][0]) / np.sum(old_grid["inside"])) * 100
                print(in_grid)
            if (out_grid > 0) and (out_grid is not np.nan) and (in_grid is not np.nan)\
                    and (in_grid > 0):
                ratio = out_grid/in_grid
                item[2] = [float(s)*ratio for s in item[1]]
            elif (out_grid == 0) or (out_grid is np.nan):
                item[2] = [0.0, 0.0, 1.0]
            elif (in_grid is np.nan) or (in_grid == 0):
                item[2] = [float(s) for s in item[1]]
                item[1] = [0.0, 0.0, 1.0]
        if item[1] == [0.0, 0.0, 1.0] or item[2] == [0.0, 0.0, 1.0]:
            working_lipids.insert(0, [item[0], item[1], item[2]])
        else:
            working_lipids.append([item[0], item[1], item[2]])
    return working_lipids


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
        if (len(lipid[0]) > 1) and (str(x['ID']) in lipid[0]):
            ids.append(lipid[0])
            index_num = ids.index(lipid[0])
            if new_con[index_num] != 0:
                new_concentration = calc_fam_factor(lipid[0], new_con[index_num], membrane_part, total_sum,
                                                    str(x['ID']))
            else:
                new_concentration = x[membrane_part] * factor
            return new_concentration
        else:
            ids.append(lipid[0][0])
    if str(x['ID']) in ids:
        index_num = ids.index(x['ID'])
        if new_con[index_num] != 0:
            new_concentration = (total_sum/100)*new_con[index_num]
        else:
            new_concentration = x[membrane_part]*factor
    else:
        new_concentration = x[membrane_part]*factor
    return new_concentration


def calc_fam_factor(lipids, concentration, membrane_part, total_sum, chosen_lipid):
    # calculates and returns the family concentration factor
    fam_sum = 0
    for lipid in lipids:
        fam_sum = fam_sum + old_grid.loc[old_grid["ID"] == str(lipid), [membrane_part]].values[0][0]
    chosen_con = old_grid.loc[old_grid["ID"] == chosen_lipid, [membrane_part]].values[0][0]
    fam_factor = ((total_sum*(concentration/100))/fam_sum) * chosen_con
    return fam_factor


# set updated composition grid file
old_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
old_grid = old_grid.fillna(0)
wanted_lipid_file = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\config_file.csv")
# membrane_parts = input("Choose membrane layer you want to change. \nFor both layers type 'total_concentration'."
#                       "\nFor outer layer type 'outside'. \nFor inner layer type 'inside':")
# wanted_lipids = [[lipid_name, inside_tuple, outside_tuple], [...], [...]]
wanted_lipids = process_config_file(wanted_lipid_file)
print(wanted_lipids)
# wanted_lipids = [["Cer4412", (3, 5, 1), (3, 5, 1)], ["Cer4202", (2, 2.5, 0.5), (0, 0, 1)],
#                 ["Cer4222", (1, 2, 1), (2, 4, 2)]]
old_con_out, index_out = calc_old_con(old_grid, wanted_lipids, "outside")
old_con_in, index_in = calc_old_con(old_grid, wanted_lipids, "inside")
array_out = calc_new_con_grid("outside")
array_in = calc_new_con_grid("inside")
print(array_out, "-", array_in)
print(wanted_lipids[1], index_in, index_out)
count = 0
if len(array_in) >= len(array_out):
    big = array_in
    small = array_out
    chosen = "in"
    zero_index_small = index_out
    zero_index_big = index_in
else:
    big = array_out
    small = array_in
    chosen = "out"
    zero_index_big = index_out
    zero_index_small = index_in
mod_factor = len(big)/len(small)
print(len(big), "  ", len(small))
for small_con in small:
    trash = [big[0]]
    current_mod_factor = mod_factor
    reg = 0
    while current_mod_factor > 0:
        bigger_con = big[reg]
        print("reg: " + str(reg) + " last: " + str(bigger_con) + "mod: " + str(current_mod_factor), "small: ",
              str(small_con))
        val_bool = []
        for index_big in zero_index_small:
            if bigger_con[index_big] == trash[-1][index_big]:
                val_bool.append(True)
            else:
                val_bool.append(False)
        if all(val_bool):
            count += 1
            current_mod_factor -= 1
            if chosen == "in":
                inner_con = bigger_con
                outer_con = small_con
            else:
                inner_con = small_con
                outer_con = bigger_con
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
            new_grid.to_csv("C:\\Users\\hfComp\\Desktop\\WS\\new_composition_grid" + wanted_lipids[0][0][0] + "_in" +
                            str(sum(inner_con)) +
                            "_out" + str(sum(outer_con)) + "_" + str(count) + ".csv")
            trash.append(big.pop(reg))
            print(trash[-1])
            reg -= 1
        reg += 1
        if len(big) == 0:
            break
