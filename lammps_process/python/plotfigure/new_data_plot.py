# read lammps output
# transform data to array
# save array
# load saved array
# calculate new defined variable
# save new variable as array
# plot from array
#===========================================
# import site package
import os
import sys
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datapath as dp
# import my package
import read_setting as rr

# define dictionary store file name and variable name and new variable name by reading log file

# list of tuples (id, outputfile, header)
# near inner wall coor
    
# coodin -- filepath header
# valuefile -- id filepath header
# realname -- names

inwall_coord_filepath_header = (
        "output/wall/chunk/coord_inwall",
        "Chunk Coord1 Coord2",
    )
outwall_coord_filepath_header = (
        "output/wall/chunk/coord_outwall",
        "Chunk Coord1 Coord2",
    )
zbottom_coord_filepath_header = (
        "output/wall/chunk/coord_zbottom",
        "Chunk Coord1 Coord2",
    )

coord_aross_yz_filepath_header = (
        "output/wall/chunk/coord_chunk_2_3",
        "Chunk Coord1 Coord2",
    )

# coord filepath header save_npy_folder_pathrela_lmp n_sample_string for all output
chunk_contact = (
        coord_aross_yz_filepath_header,
        "output/contact/fix.contact_number.all",
        "Chunk n_contact n_contact_sq",
        "postprocess/npy/chunk_contact/",
        "repeat_ave_chunk_momentum_mass_field"
    )

chunk_mv_Ek_mass = (
        coord_aross_yz_filepath_header,
        "output/momentum_mass_field/fix.momentum_mass_field.all",
        "Chunk Ncount v_mv2 v_mv2_sq v_mv1 v_mv1_sq v_mv3 v_mv3_sq v_Ek1 v_Ek2 v_Ek3 c_m1 c_m1_sq",
        "postprocess/npy/chunk_mv_Ek_mass/",
        "repeat_ave_chunk_momentum_mass_field"
    )
chunk_omega = (
        coord_aross_yz_filepath_header,
        "output/momentum_mass_field/fix.omega.all",
        "Chunk Ncount c_omega[1] c_omega[1]_sq c_omega[2] c_omega[2]_sq c_omega[3] c_omega[3]_sq",
        "postprocess/npy/chunk_omega/",
        "repeat_ave_chunk_momentum_mass_field"
    )
chunk_stress = (
        coord_aross_yz_filepath_header,
        "output/stress/fix.stress.all",
        "Chunk Ncount c_stress[1] c_stress[1]_sq c_stress[2] c_stress[2]_sq c_stress[3] c_stress[3]_sq c_stress[4] c_stress[4]_sq c_stress[5] c_stress[5]_sq c_stress[6] c_stress[6]_sq",
        "postprocess/npy/chunk_stress/",
        "repeat_ave_chunk_momentum_mass_field",
    )

chunk_inwall_force = (
        inwall_coord_filepath_header,
        "output/wall/chunk/inwallforcefile",
        "Chunk Ncount v_inwall_per_atom_1 v_inwall_per_atom_1_sq v_inwall_per_atom_2 v_inwall_per_atom_2_sq v_inwall_per_atom_3 v_inwall_per_atom_3_sq",
        "postprocess/npy/chunk_inwall_force/",
        "repeat_ave_chunk_wallforce",
    )
chunk_outwall_force = (
        outwall_coord_filepath_header,
        "output/wall/chunk/outwallforcefile",
        "Chunk Ncount v_outwall_per_atom_1 v_outwall_per_atom_1_sq v_outwall_per_atom_2 v_outwall_per_atom_2_sq v_outwall_per_atom_3 v_outwall_per_atom_3_sq",
        "postprocess/npy/chunk_outwall_force/",
        "repeat_ave_chunk_wallforce",
    )
chunk_zbottom_force = (
        zbottom_coord_filepath_header,
        "output/wall/chunk/zbottomforcefile",
        "Chunk Ncount v_zbottom_per_atom_1 v_zbottom_per_atom_1_sq v_zbottom_per_atom_2 v_zbottom_per_atom_2_sq v_zbottom_per_atom_3 v_zbottom_per_atom_3_sq",
        "postprocess/npy/chunk_zbottom_force/",
        "repeat_ave_chunk_wallforce",
    )

# filepath, char of coord, timestep path, save_folder_path
map_chunkfile_char_save_folderpath = {
    "contact":{
        "filepath": "output/contact/fix.contact_number.all",
        "char_of_coord": 'original_2D',
        "savepath": "postprocess/npy/chunk_contact/",
    },
    "mv_Ek_mass": {
        "filepath": "output/momentum_mass_field/fix.momentum_mass_field.all",
        "char_of_coord": 'original_2D',
        "savepath": "postprocess/npy/chunk_mv_Ek_mass/",
    },
    "omega": {
        "filepath": "output/momentum_mass_field/fix.omega.all",
        "char_of_coord": 'original_2D',
        "savepath": "postprocess/npy/chunk_omega/",
    },
    "stress":{
        "filepath": "output/stress/fix.stress.all",
        "char_of_coord": 'original_2D',
        "savepath": "postprocess/npy/chunk_stress/",
    },
    "inwall":{
        "filepath": "output/wall/chunk/inwallforcefile",
        "char_of_coord": 'inwall_1D',
        "savepath": "postprocess/npy/chunk_inwall_force/",
    },
    "outwall":{
        "filepath": "output/wall/chunk/outwallforcefile",
        "char_of_coord": 'outwall_1D',
        "savepath": "postprocess/npy/chunk_outwall_force/",
    },
    "zbottom":{
        "filepath": "output/wall/chunk/zbottomforcefile",
        "char_of_coord": 'zbottom_1D',
        "savepath": "postprocess/npy/chunk_zbottom_force/",
    },
}

## add Ncount path
#for key in map_chunkfile_char_save_folderpath.keys():
#    map_chunkfile_char_save_folderpath[key]["Ncount_npy_path"] = map_chunkfile_char_save_folderpath[key]["savepath"] + "Ncount.npy"

# update timestep path
for key in map_chunkfile_char_save_folderpath.keys():
    map_chunkfile_char_save_folderpath[key]["timestep_path"] = map_chunkfile_char_save_folderpath[key]["savepath"] + "timestep.npy"

coordinate_characteristic_to_npyfilepath = {
    "original_2D": {
        "Coord1": {"headername": "Coord1", 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/Coord1.npy", "map_chunkfile": "mv_Ek_mass"},
        "Coord2": {"headername": "Coord2", 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/Coord2.npy", "map_chunkfile": "mv_Ek_mass"},
    },
    "inwall_1D": {
        "Coord2": {"headername": "Coord2", 'npyfilepath': "postprocess/npy/chunk_inwall_force/Coord2.npy", "map_chunkfile": "inwall"},
    },
    "outwall_1D": {
        "Coord2": {"headername": "Coord2", 'npyfilepath': "postprocess/npy/chunk_outwall_force/Coord2.npy", "map_chunkfile": "outwall"},
    },
    "zbottom_1D": {
        "Coord1": {"headername": "Coord1", 'npyfilepath': "postprocess/npy/chunk_zbottom_force/Coord1.npy", "map_chunkfile": "zbottom"},
    },
}
for key1 in coordinate_characteristic_to_npyfilepath.keys():
    for key2 in coordinate_characteristic_to_npyfilepath[key1].keys():
        key3 = coordinate_characteristic_to_npyfilepath[key1][key2]["map_chunkfile"]
        for key4 in map_chunkfile_char_save_folderpath[key3].keys():
            coordinate_characteristic_to_npyfilepath[key1][key2][key4] = map_chunkfile_char_save_folderpath[key3][key4]


calculation_coordinate_characteristic_to_npyfilepath = {
    "middle_2_2D": {
        "Coord1": {'npyfilepath': "postprocess/npy/calculate/coord_middle_2/Coord1.npy"},
        "Coord2": {'npyfilepath': "postprocess/npy/calculate/coord_middle_2/Coord2.npy"},
    },
    "middle_3_2D": {
        "Coord1": {'npyfilepath': "postprocess/npy/calculate/coord_middle_3/Coord1.npy"},
        "Coord2": {'npyfilepath': "postprocess/npy/calculate/coord_middle_3/Coord2.npy"},
    },
    "middle_23_2D": {
        "Coord1": {'npyfilepath': "postprocess/npy/calculate/coord_middle_23_coord/Coord1.npy"},
        "Coord2": {'npyfilepath': "postprocess/npy/calculate/coord_middle_23_coord/Coord2.npy"},
    },
}
# readfromtext variable_name_to_npyfilepath_and_coordinate_characteristic
r_npyfilepath_coord_char = {
    "Ncount": {'headername': 'Ncount', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/Ncount.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "n_contact": {'headername': 'n_contact', 'npyfilepath': "postprocess/npy/chunk_contact/n_contact.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "contact"},
    "n_contact_std": {'headername': 'n_contact_std', 'npyfilepath': "postprocess/npy/chunk_contact/n_contact_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "contact"},
    "mass": {'headername': 'c_m1', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/c_m1.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mass_std": {'headername': 'c_m1_std', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/c_m1_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mv_1": {'headername': 'v_mv1', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_mv1.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mv_1_std": {'headername': 'v_mv1_std', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_mv1_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mv_2": {'headername': 'v_mv2', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_mv2.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mv_2_std": {'headername': 'v_mv2_std', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_mv2_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mv_3": {'headername': 'v_mv3', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_mv3.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mv_3_std": {'headername': 'v_mv3_std', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_mv3_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "Ek_1": {'headername': 'v_Ek1', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_Ek1.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "Ek_1_std": {'headername': 'v_Ek1_std', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_Ek1_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "Ek_2": {'headername': 'v_Ek2', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_Ek2.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "Ek_2_std": {'headername': 'v_Ek2_std', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_Ek2_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "Ek_3": {'headername': 'v_Ek3', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_Ek3.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "Ek_3_std": {'headername': 'v_Ek3_std', 'npyfilepath': "postprocess/npy/chunk_mv_Ek_mass/v_Ek3_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "omega_1": {'headername': 'c_omega[1]', 'npyfilepath': "postprocess/npy/chunk_omega/c_omega[1].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "omega"},
    "omega_1_std": {'headername': 'c_omega[1]_std', 'npyfilepath': "postprocess/npy/chunk_omega/c_omega[1]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "omega"},
    "omega_2": {'headername': 'c_omega[2]', 'npyfilepath': "postprocess/npy/chunk_omega/c_omega[2].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "omega"},
    "omega_2_std": {'headername': 'c_omega[2]_std', 'npyfilepath': "postprocess/npy/chunk_omega/c_omega[2]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "omega"},
    "omega_3": {'headername': 'c_omega[3]', 'npyfilepath': "postprocess/npy/chunk_omega/c_omega[3].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "omega"},
    "omega_3_std": {'headername': 'c_omega[3]_std', 'npyfilepath': "postprocess/npy/chunk_omega/c_omega[3]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "omega"},
    "stress_multiply_binvolume_11": {'headername': 'c_stress[1]', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[1].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_11_std": {'headername': 'c_stress[1]_std', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[1]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_22": {'headername': 'c_stress[2]', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[2].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_22_std": {'headername': 'c_stress[2]_std', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[2]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_33": {'headername': 'c_stress[3]', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[3].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_33_std": {'headername': 'c_stress[3]_std', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[3]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_12": {'headername': 'c_stress[4]', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[4].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_12_std": {'headername': 'c_stress[4]_std', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[4]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_13": {'headername': 'c_stress[5]', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[5].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_13_std": {'headername': 'c_stress[5]_std', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[5]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_23": {'headername': 'c_stress[6]', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[6].npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "stress_multiply_binvolume_23_std": {'headername': 'c_stress[6]_std', 'npyfilepath': "postprocess/npy/chunk_stress/c_stress[6]_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "stress"},
    "inwall_force_1": {'headername': 'v_inwall_per_atom_1', 'npyfilepath': "postprocess/npy/chunk_inwall_force/v_inwall_per_atom_1.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_force_1_std": {'headername': 'v_inwall_per_atom_1_std', 'npyfilepath': "postprocess/npy/chunk_inwall_force/v_inwall_per_atom_1_std.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_force_2": {'headername': 'v_inwall_per_atom_2', 'npyfilepath': "postprocess/npy/chunk_inwall_force/v_inwall_per_atom_2.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_force_2_std": {'headername': 'v_inwall_per_atom_2_std', 'npyfilepath': "postprocess/npy/chunk_inwall_force/v_inwall_per_atom_2_std.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_force_3": {'headername': 'v_inwall_per_atom_3', 'npyfilepath': "postprocess/npy/chunk_inwall_force/v_inwall_per_atom_3.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_force_3_std": {'headername': 'v_inwall_per_atom_3_std', 'npyfilepath': "postprocess/npy/chunk_inwall_force/v_inwall_per_atom_3_std.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "outwall_force_1": {'headername': 'v_outwall_per_atom_1', 'npyfilepath': "postprocess/npy/chunk_outwall_force/v_outwall_per_atom_1.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_force_1_std": {'headername': 'v_outwall_per_atom_1_std', 'npyfilepath': "postprocess/npy/chunk_outwall_force/v_outwall_per_atom_1_std.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_force_2": {'headername': 'v_outwall_per_atom_2', 'npyfilepath': "postprocess/npy/chunk_outwall_force/v_outwall_per_atom_2.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_force_2_std": {'headername': 'v_outwall_per_atom_2_std', 'npyfilepath': "postprocess/npy/chunk_outwall_force/v_outwall_per_atom_2_std.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_force_3": {'headername': 'v_outwall_per_atom_3', 'npyfilepath': "postprocess/npy/chunk_outwall_force/v_outwall_per_atom_3.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_force_3_std": {'headername': 'v_outwall_per_atom_3_std', 'npyfilepath': "postprocess/npy/chunk_outwall_force/v_outwall_per_atom_3_std.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "zbottom_force_1": {'headername': 'v_zbottom_per_atom_1', 'npyfilepath': "postprocess/npy/chunk_zbottom_force/v_zbottom_per_atom_1.npy", 'coordinate_characteristic': 'zbottom_1D', "map_chunkfile": "zbottom"},
    "zbottom_force_1_std": {'headername': 'v_zbottom_per_atom_1_std', 'npyfilepath': "postprocess/npy/chunk_zbottom_force/v_zbottom_per_atom_1_std.npy", 'coordinate_characteristic': 'zbottom_1D', "map_chunkfile": "zbottom"},
    "zbottom_force_2": {'headername': 'v_zbottom_per_atom_2', 'npyfilepath': "postprocess/npy/chunk_zbottom_force/v_zbottom_per_atom_2.npy", 'coordinate_characteristic': 'zbottom_1D', "map_chunkfile": "zbottom"},
    "zbottom_force_2_std": {'headername': 'v_zbottom_per_atom_2_std', 'npyfilepath': "postprocess/npy/chunk_zbottom_force/v_zbottom_per_atom_2_std.npy", 'coordinate_characteristic': 'zbottom_1D', "map_chunkfile": "zbottom"},
    "zbottom_force_3": {'headername': 'v_zbottom_per_atom_3', 'npyfilepath': "postprocess/npy/chunk_zbottom_force/v_zbottom_per_atom_3.npy", 'coordinate_characteristic': 'zbottom_1D', "map_chunkfile": "zbottom"},
    "zbottom_force_3_std": {'headername': 'v_zbottom_per_atom_3_std', 'npyfilepath': "postprocess/npy/chunk_zbottom_force/v_zbottom_per_atom_3_std.npy", 'coordinate_characteristic': 'zbottom_1D', "map_chunkfile": "zbottom"},
}

for key2 in r_npyfilepath_coord_char.keys():
    key3 = r_npyfilepath_coord_char[key2]["map_chunkfile"]
    for key4 in map_chunkfile_char_save_folderpath[key3].keys():
        r_npyfilepath_coord_char[key2][key4] = map_chunkfile_char_save_folderpath[key3][key4]

# map better name to npy file name. And npy folder
# calculation variable_name_to_npyfilepath_and_coordinate_characteristic
c_npyfilepath_coord_char = {
    "velocity_1": {'npyfilepath': "postprocess/npy/calculate/velocity_1.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "velocity_1_std": {'npyfilepath': "postprocess/npy/calculate/velocity_1_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "velocity_2": {'npyfilepath': "postprocess/npy/calculate/velocity_2.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "velocity_2_std": {'npyfilepath': "postprocess/npy/calculate/velocity_2_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "velocity_3": {'npyfilepath': "postprocess/npy/calculate/velocity_3.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "velocity_3_std": {'npyfilepath': "postprocess/npy/calculate/velocity_3_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_21": {'npyfilepath': "postprocess/npy/calculate/strain_rate_21.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_21_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_21_std.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"}, # middle_2_2D means middle point along axis 2
    "strain_rate_22": {'npyfilepath': "postprocess/npy/calculate/strain_rate_22.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_22_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_22_std.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_23": {'npyfilepath': "postprocess/npy/calculate/strain_rate_23.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_23_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_23_std.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_31": {'npyfilepath': "postprocess/npy/calculate/strain_rate_31.npy", 'coordinate_characteristic': 'middle_3_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_31_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_31_std.npy", 'coordinate_characteristic': 'middle_3_2D', "map_chunkfile": "mv_Ek_mass"}, # middle_3_2D means middle point along axis 3
    "strain_rate_32": {'npyfilepath': "postprocess/npy/calculate/strain_rate_32.npy", 'coordinate_characteristic': 'middle_3_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_32_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_32_std.npy", 'coordinate_characteristic': 'middle_3_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_33": {'npyfilepath': "postprocess/npy/calculate/strain_rate_33.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_33_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_33_std.npy", 'coordinate_characteristic': 'middle_2_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_21_middle": {'npyfilepath': "postprocess/npy/calculate/strain_rate_21_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_21_middle_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_21_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"}, # middle_23_2D means middle point along axis 2
    "strain_rate_22_middle": {'npyfilepath': "postprocess/npy/calculate/strain_rate_22_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_22_middle_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_22_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_23_middle": {'npyfilepath': "postprocess/npy/calculate/strain_rate_23_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_23_middle_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_23_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_31_middle": {'npyfilepath': "postprocess/npy/calculate/strain_rate_31_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_31_middle_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_31_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"}, # middle_23_2D means middle point along axis 3
    "strain_rate_32_middle": {'npyfilepath': "postprocess/npy/calculate/strain_rate_32_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_32_middle_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_32_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_33_middle": {'npyfilepath': "postprocess/npy/calculate/strain_rate_33_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "strain_rate_33_middle_std": {'npyfilepath': "postprocess/npy/calculate/strain_rate_33_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_11": {'npyfilepath': "postprocess/npy/calculate/stress_11.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_11_std": {'npyfilepath': "postprocess/npy/calculate/stress_11_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_22": {'npyfilepath': "postprocess/npy/calculate/stress_22.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_22_std": {'npyfilepath': "postprocess/npy/calculate/stress_22_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_33": {'npyfilepath': "postprocess/npy/calculate/stress_33.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_33_std": {'npyfilepath': "postprocess/npy/calculate/stress_33_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_12": {'npyfilepath': "postprocess/npy/calculate/stress_12.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_12_std": {'npyfilepath': "postprocess/npy/calculate/stress_12_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_13": {'npyfilepath': "postprocess/npy/calculate/stress_13.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_13_std": {'npyfilepath': "postprocess/npy/calculate/stress_13_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_23": {'npyfilepath': "postprocess/npy/calculate/stress_23.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "stress_23_std": {'npyfilepath': "postprocess/npy/calculate/stress_23_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "pressure": {'npyfilepath': "postprocess/npy/calculate/pressure.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "pressure_std": {'npyfilepath': "postprocess/npy/calculate/pressure_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_12": {'npyfilepath': "postprocess/npy/calculate/mu_12.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_12_std": {'npyfilepath': "postprocess/npy/calculate/mu_12_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_13": {'npyfilepath': "postprocess/npy/calculate/mu_13.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_13_std": {'npyfilepath': "postprocess/npy/calculate/mu_13_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_23": {'npyfilepath': "postprocess/npy/calculate/mu_23.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_23_std": {'npyfilepath': "postprocess/npy/calculate/mu_23_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_12_middle": {'npyfilepath': "postprocess/npy/calculate/mu_12_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_12_middle_std": {'npyfilepath': "postprocess/npy/calculate/mu_12_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_13_middle": {'npyfilepath': "postprocess/npy/calculate/mu_13_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_13_middle_std": {'npyfilepath': "postprocess/npy/calculate/mu_13_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_23_middle": {'npyfilepath': "postprocess/npy/calculate/mu_23_middle.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_23_middle_std": {'npyfilepath': "postprocess/npy/calculate/mu_23_middle_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_tensor_12": {'npyfilepath': "postprocess/npy/calculate/mu_tensor_12.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_tensor_12_std": {'npyfilepath': "postprocess/npy/calculate/mu_tensor_12_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_tensor_13": {'npyfilepath': "postprocess/npy/calculate/mu_tensor_13.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_tensor_13_std": {'npyfilepath': "postprocess/npy/calculate/mu_tensor_13_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_tensor_23": {'npyfilepath': "postprocess/npy/calculate/mu_tensor_23.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "mu_tensor_23_std": {'npyfilepath': "postprocess/npy/calculate/mu_tensor_23_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_12": {'npyfilepath': "postprocess/npy/calculate/I_12.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_12_std": {'npyfilepath': "postprocess/npy/calculate/I_12_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_13": {'npyfilepath': "postprocess/npy/calculate/I_13.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_13_std": {'npyfilepath': "postprocess/npy/calculate/I_13_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_23": {'npyfilepath': "postprocess/npy/calculate/I_23.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_23_std": {'npyfilepath': "postprocess/npy/calculate/I_23_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_tensor": {'npyfilepath': "postprocess/npy/calculate/I_tensor.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "I_tensor_std": {'npyfilepath': "postprocess/npy/calculate/I_tensor_std.npy", 'coordinate_characteristic': 'middle_23_2D', "map_chunkfile": "mv_Ek_mass"},
    "fraction": {'npyfilepath': "postprocess/npy/calculate/fraction.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "fraction_std": {'npyfilepath': "postprocess/npy/calculate/fraction_std.npy", 'coordinate_characteristic': 'original_2D', "map_chunkfile": "mv_Ek_mass"},
    "inwall_stress_1": {'npyfilepath': "postprocess/npy/calculate/inwall_stress_1.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_stress_1_std": {'npyfilepath': "postprocess/npy/calculate/inwall_stress_1_std.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_stress_2": {'npyfilepath': "postprocess/npy/calculate/inwall_stress_2.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_stress_2_std": {'npyfilepath': "postprocess/npy/calculate/inwall_stress_2_std.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_stress_3": {'npyfilepath': "postprocess/npy/calculate/inwall_stress_3.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "inwall_stress_3_std": {'npyfilepath': "postprocess/npy/calculate/inwall_stress_3_std.npy", 'coordinate_characteristic': 'inwall_1D', "map_chunkfile": "inwall"},
    "outwall_stress_1": {'npyfilepath': "postprocess/npy/calculate/outwall_stress_1.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_stress_1_std": {'npyfilepath': "postprocess/npy/calculate/outwall_stress_1_std.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_stress_2": {'npyfilepath': "postprocess/npy/calculate/outwall_stress_2.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_stress_2_std": {'npyfilepath': "postprocess/npy/calculate/outwall_stress_2_std.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_stress_3": {'npyfilepath': "postprocess/npy/calculate/outwall_stress_3.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "outwall_stress_3_std": {'npyfilepath': "postprocess/npy/calculate/outwall_stress_3_std.npy", 'coordinate_characteristic': 'outwall_1D', "map_chunkfile": "outwall"},
    "zbottom_stress_1": {'npyfilepath': "postprocess/npy/calculate/zbottom_stress_1.npy", 'coordinate_characteristic': 'zbottom_1D',"map_chunkfile": "zbottom"},
    "zbottom_stress_1_std": {'npyfilepath': "postprocess/npy/calculate/zbottom_stress_1_std.npy", 'coordinate_characteristic': 'zbottom_1D',"map_chunkfile": "zbottom"},
    "zbottom_stress_2": {'npyfilepath': "postprocess/npy/calculate/zbottom_stress_2.npy", 'coordinate_characteristic': 'zbottom_1D',"map_chunkfile": "zbottom"},
    "zbottom_stress_2_std": {'npyfilepath': "postprocess/npy/calculate/zbottom_stress_2_std.npy", 'coordinate_characteristic': 'zbottom_1D',"map_chunkfile": "zbottom"},
    "zbottom_stress_3": {'npyfilepath': "postprocess/npy/calculate/zbottom_stress_3.npy", 'coordinate_characteristic': 'zbottom_1D',"map_chunkfile": "zbottom"},
    "zbottom_stress_3_std": {'npyfilepath': "postprocess/npy/calculate/zbottom_stress_3_std.npy", 'coordinate_characteristic': 'zbottom_1D',"map_chunkfile": "zbottom"},
}

for key2 in c_npyfilepath_coord_char.keys():
    key3 = c_npyfilepath_coord_char[key2]["map_chunkfile"]
    for key4 in map_chunkfile_char_save_folderpath[key3].keys():
        c_npyfilepath_coord_char[key2][key4] = map_chunkfile_char_save_folderpath[key3][key4]

# check if two dic contain same key
for key in r_npyfilepath_coord_char.keys():
    if key in c_npyfilepath_coord_char.keys():
        sys.exit("the key " + key + "appear in both dictionaries")

# combine two dics
c_r_npyfilepath_coord_char = c_npyfilepath_coord_char.copy()
c_r_npyfilepath_coord_char.update(r_npyfilepath_coord_char)

# create folder for npyfilepath in c_r_npyfilepath_coord_char
folder_list = [
    os.path.dirname(v['npyfilepath']) for k, v in c_r_npyfilepath_coord_char.items()
]
set_folder = set(folder_list)
for folder in set_folder:
    os.makedirs(folder, exist_ok=True)

# check if two coord dic contain same key
for key in coordinate_characteristic_to_npyfilepath.keys():
    if key in calculation_coordinate_characteristic_to_npyfilepath.keys():
        sys.exit("the key " + key + "appear in both dictionaries")

# combine two coord dics
c_r_coord_npyfilepath = coordinate_characteristic_to_npyfilepath.copy()
c_r_coord_npyfilepath.update(calculation_coordinate_characteristic_to_npyfilepath)

# create folder for npyfilepath in c_r_coord_npyfilepath
for k, dic_for_specific_char_coord in c_r_coord_npyfilepath.items():
    # create folder for npyfilepath in dic_for_specific_char_coord
    folder_list = [
        os.path.dirname(v['npyfilepath']) for k, v in dic_for_specific_char_coord.items()
    ]
    set_folder = set(folder_list)
    for folder in set_folder:
        os.makedirs(folder, exist_ok=True)

def save_v_array_to_disk(lmp_path, name, value):
    path = lmp_path + c_r_npyfilepath_coord_char[name]['npyfilepath']
    np.save(path, value)

def save_coord_to_disk(lmp_path, char, name, value):
    path = lmp_path + c_r_coord_npyfilepath[char][name]['npyfilepath']
    np.save(path, value)

def load_v_array(lmp_path, name):
    path = lmp_path + c_r_npyfilepath_coord_char[name]['npyfilepath']
    return np.load(path, mmap_mode='r')

def load_coord(lmp_path, char, name):
    path = lmp_path + c_r_coord_npyfilepath[char][name]['npyfilepath']
    return np.load(path, mmap_mode='r')

# no chunk output
inwall_force = (
    "output/wall/force_y_bottom_to_particle.allstep",
    "TimeStep v_t v_force_inwall_1 v_force_inwall_2 v_force_inwall_3", #v_force_zbottom_x v_force_zbottom_y v_force_zbottom_z
    "postprocess/npy/nochunk_wall/",
)
outwall_force = (
    "output/wall/force_y_top_to_particle.allstep",
    "TimeStep v_t v_force_outwall_1 v_force_outwall_2 v_force_outwall_3",
    "postprocess/npy/nochunk_wall/",
)
zbottom_force = (
    "output/wall/force_zbottom_to_particle.allstep",
    "TimeStep v_t v_force_zbottom_1 v_force_zbottom_2 v_force_zbottom_3",
    "postprocess/npy/nochunk_wall/",
)

nochunk_npyfilepath_coord_char = {
    "v_force_y_top_1": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_top_1.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_outwall_1": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_outwall_1.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_top_x": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_top_x.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_top_2": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_top_2.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_outwall_2": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_outwall_2.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_top_y": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_top_y.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_top_3": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_top_3.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_outwall_3": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_outwall_3.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_top_z": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_top_z.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_bottom_1": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_bottom_1.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_inwall_1": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_inwall_1.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_bottom_x": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_bottom_x.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_bottom_2": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_bottom_2.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_inwall_2": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_inwall_2.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_bottom_y": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_bottom_y.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_bottom_3": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_bottom_3.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_inwall_3": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_inwall_3.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_y_bottom_z": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_y_bottom_z.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_zbottom_1": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_zbottom_1.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_zbottom_x": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_zbottom_x.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_zbottom_2": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_zbottom_2.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_zbottom_y": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_zbottom_y.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_zbottom_3": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_zbottom_3.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
    "v_force_zbottom_z": {'npyfilepath': "postprocess/npy/nochunk_wall/v_force_zbottom_z.npy", "timestep_path":"postprocess/npy/nochunk_wall/timestep.npy"},
}

# n_1 n_2 position_index_to_array_dim_index
if rr.logfile["shearwall"] == "zcylinder":
    position_index_to_array_dim_index = {
        1: 1,
        2: 0,
    }
    n_r = int(rr.logfile['N_bin_r'])
    n_z = int(rr.logfile['N_bin_z'])
    n_1 = n_z
    n_2 = n_r
    n_12 = n_1*n_2
    
    dx = 1/n_r*int(rr.logfile['width_wall_dp_unit'])
    dy = 1/n_z*float(rr.logfile['zhi_chunk_dp_unit'])
    x_array, y_array = np.meshgrid(
                                int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*int(rr.logfile['width_wall_dp_unit']),
                                (np.arange(n_2)+0.5)/n_2*float(rr.logfile['zhi_chunk_dp_unit']),
                                )
    x_array = x_array.reshape((-1))
    y_array = y_array.reshape((-1))
    vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*float(rr.logfile['dp'])**3
elif rr.logfile["shearwall"] == "yplane":
    n_y = int(rr.logfile['N_bin_y'])
    n_z = int(rr.logfile['N_bin_z'])
    if "chunk/atom 23" in rr.logfile.keys():
        if rr.logfile["chunk/atom 23"][1] == "y":
            position_index_to_array_dim_index = {
                                            1: 0,
                                            2: 1,
                                            }
            n_1 = n_y
            n_2 = n_z
        elif rr.logfile["chunk/atom 23"][1] == "z":
            position_index_to_array_dim_index = {
                                            2: 0,
                                            1: 1,
                                            }
            n_1 = n_z
            n_2 = n_y
        else:
            sys.exit("chunk_method wrong")
    else:
        position_index_to_array_dim_index = {
                                        1: 0,
                                        2: 1,
                                        }
        n_1 = n_y
        n_2 = n_z

    n_12 = n_1*n_2
    
else:
    sys.exit("chunk_method wrong")

def count_n_2(lmp_path):
    # python read line
    with open(lmp_path + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
        # line list
        lines = f.read().strip().split('\n')
    n_line_in_file = len(lines)
    first_4_lines = lines[0: 4]
    # header
    header = first_4_lines[2].split()[1:]
    # stepsarray
    n_line_in_a_step = int(first_4_lines[3].split()[1])
    n_2 = int(n_line_in_a_step/n_1)
    return n_2

n_2 = count_n_2(rr.folder_path_list_initial_to_last[-1])
dx = 1/n_1*int(rr.logfile['width_wall_dp_unit'])
dy = 1/n_2*float(rr.logfile['zhi_chunk_dp_unit'])
vol_in_chunks = float(rr.logfile['dp'])*int(rr.logfile['x_period_dp_unit'])*dx*dy*float(rr.logfile['dp'])**2
def chunk_coord_data_transfer_text_to_npy(lmp_path, filepath_rela_lmp_path, outputfolder_rela_lmp_path):
   
    # python read line
    with open(lmp_path + filepath_rela_lmp_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    n_line_in_file = len(lines)
    first_4_lines = lines[0: 4]
    # header
    header = first_4_lines[2].split()[1:]
    # stepsarray
    n_line_in_a_step = int(first_4_lines[3].split()[1])
    step_first_in_file = int(first_4_lines[3].split()[0])
    
    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[4: 4+n_line_in_a_step]
    
    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    n_line = lines.shape[0]
    n_header = len(header)
    os.makedirs(lmp_path + outputfolder_rela_lmp_path, exist_ok=True)
    
    # save all variable in the folder
    for i in range(n_header):
        lines_column = lines[:, i]
        if n_line_in_a_step == n_1*n_2:
            #breakpoint()
            lines_column = lines_column.reshape((n_1, n_2))
        else:
            pass
        outputfilepath = lmp_path + outputfolder_rela_lmp_path + header[i]
        np.save(outputfilepath, lines_column)

# function get text data from file
# save all header variable to npy
def chunk_variable_data_transfer_text_to_npy(lmp_path, filepath_rela_lmp_path, outputfolder_rela_lmp_path, n_sample_keyword):
    # python read line
    with open(lmp_path + filepath_rela_lmp_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    n_line_in_file = len(lines)
    first_4_lines = lines[0: 4]
    # header
    header = first_4_lines[2].split()[1:]
    # stepsarray
    n_line_in_a_step = int(first_4_lines[3].split()[1])
    step_first_in_file = int(first_4_lines[3].split()[0])
    step_second_in_file = int(lines[3 + n_line_in_a_step + 1].split()[0])
    step_last_in_file = int(lines[-1 - n_line_in_a_step].split()[0])
    d_step = step_second_in_file - step_first_in_file
    totalstepsarray = np.arange(step_first_in_file, step_last_in_file + d_step, d_step, dtype=int)
    n_step = len(totalstepsarray)
    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[3: ]
    for i in reversed(range(len(totalstepsarray))):
        # del from largest
        del lines[(n_line_in_a_step + 1)*i]
    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    n_line = lines.shape[0]
    if n_line % n_step == 0:
        n_line_one_step = int(n_line/n_step)
    else:
        sys.exit("n_step not divisible n_line")
    n_header = len(header)
    os.makedirs(lmp_path + outputfolder_rela_lmp_path, exist_ok=True)
    
    # if Ncount in header
    if 'Ncount' in header:
        # reset variable if Ncount*freq_repeat is smaller than 1
        index_of_Ncount_in_header = header.index('Ncount')
        Ncount_column = lines[:, index_of_Ncount_in_header]
        reset_mask = (
            Ncount_column*int(rr.logfile[n_sample_keyword]) < 0.99
        )
    # save all variable in the folder
    for i in range(n_header):
        lines_column = lines[:, i]
        # reset value to zero if Ncount*freq_repeat is smaller than 1
        if 'Ncount' in header:
            lines_column[reset_mask] = 0
        if n_line_in_a_step == n_1*n_2:
            #breakpoint()
            lines_column = lines_column.reshape((n_step, n_1, n_2))
        elif n_line_in_a_step == n_1:
            lines_column = lines_column.reshape((n_step, n_1))
        elif n_line_in_a_step == n_2:
            lines_column = lines_column.reshape((n_step, n_2))
        else:
            sys.exit("n_line_in_a_step is not n1*n2 or n1 or n2")
        
        np.save(
            lmp_path + outputfolder_rela_lmp_path + header[i],
            lines_column,
        )
    # save timestep in the folder
    np.save(lmp_path + outputfolder_rela_lmp_path + "timestep", totalstepsarray)

def data_not_chunk_transfer_text_to_npy(lmp_path, filepath_rela_lmp_path, outputfolder_rela_lmp_path):
    # python read line
    with open(lmp_path + filepath_rela_lmp_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    n_line_in_file = len(lines)
    first_4_lines = lines[0: 4]
    # header
    if first_4_lines[1].split()[0] == "#":
        header = first_4_lines[1].split()[1:]
    else:
        header = first_4_lines[1].split()[0:]
    # stepsarray
    step_first_in_file = int(first_4_lines[2].split()[0])
    step_second_in_file = int(first_4_lines[3].split()[0])
    step_last_in_file = int(lines[-1].split()[0])
    d_step = step_second_in_file - step_first_in_file
    totalstepsarray = np.arange(step_first_in_file, step_last_in_file + d_step, d_step, dtype=int)
    n_step = len(totalstepsarray)
    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[2: ]
    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    n_line = lines.shape[0]
    if n_line % n_step == 0:
        n_line_one_step = int(n_line/n_step)
    else:
        sys.exit("n_step not divisible n_line")
    n_header = len(header)
    os.makedirs(lmp_path + outputfolder_rela_lmp_path, exist_ok=True)
    # save all variable in the folder
    for i in range(n_header):
        lines_column = lines[:, i]
        outputfilepath = lmp_path + outputfolder_rela_lmp_path + header[i]
        np.save(outputfilepath, lines_column)
    # save timestep in the folder
    np.save(lmp_path + outputfolder_rela_lmp_path + "timestep", totalstepsarray)

def text_to_npy_for_chunk(lmp_path, chunk_tuple):
    # chunk_tuple = coord filepath header save_npy_folder_path
    (coord_filepath, coord_header) = chunk_tuple[0]
    value_filepath = chunk_tuple[1]
    value_header = chunk_tuple[2]

    chunk_coord_data_transfer_text_to_npy(lmp_path, coord_filepath, chunk_tuple[3])
    chunk_variable_data_transfer_text_to_npy(lmp_path, value_filepath, chunk_tuple[3], chunk_tuple[4])

# select variable list and step
# take average
def variablesteps(lmp_path, npyfilefolderpath_rela_lmp_path, npyfilename, n_ave, inputstepsarray):
    if n_ave % 2 == 1:
        pass
    else:
        sys.exit("n_ave should be odd integer")
    # combine with previous
    # pick the smallest step and max step
    # pick the smallest step and largest step in this simulation
    def get_step_first_end(lmp_path_folder):
        timestep = np.load(lmp_path_folder + npyfilefolderpath_rela_lmp_path + "timestep" + ".npy", mmap_mode='r')
        timestep_first = timestep[int((n_ave-1)/2)]
        timestep_last = timestep[-1-int((n_ave-1)/2)]
        return (timestep_first, timestep_last)
    # input step
    step_first = inputstepsarray[0]
    # get all timestep_first timestep_last from all pre simu
    # calculate how many pre simu we need
    n_include_pre_simu = 0
    for i in range(rr.n_simu_total):
        n_include_pre_simu = n_include_pre_simu + 1
        if step_first >= get_step_first_end(
            rr.folder_path_list_initial_to_last[rr.n_simu_total-1-i]
            )[0]:
            break
    timestep = np.load(lmp_path + npyfilefolderpath_rela_lmp_path + "timestep" + ".npy", mmap_mode='r')
    d_step = timestep[1] - timestep[0]
    array = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + npyfilefolderpath_rela_lmp_path + npyfilename, mmap_mode='r')
    for j in range(n_include_pre_simu-1):
        array_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + npyfilefolderpath_rela_lmp_path + npyfilename, mmap_mode='r')
        array = np.append(array, array_append, axis=0)
    indexesarray = (inputstepsarray-timestep[0])/d_step
    indexesarray = indexesarray.astype(int)
    # average
    sum_over_n_array = 0
    for k in range(-n_ave, n_ave+1):
        sum_over_n_array = sum_over_n_array + array[indexesarray + k]
    ave_array = sum_over_n_array/(2*n_ave + 1)
    return ave_array

def calculate_std_by_ave_and_sq_ave(n_sample, ave, ave_sq):
    if np.any(ave_sq<0):
        pass
        #breakpoint()
    std_sq = ave_sq - ave**2
    std_sq[np.logical_and(-10**-20 < ave_sq, ave_sq < 0)] = 0
    std_sq[(ave_sq==0)]=0
    ratio = np.divide(std_sq, ave_sq, out=np.zeros_like(std_sq), where=ave_sq!=0)
    filter_condition_to_zero = np.logical_and(-10**-6 < ratio, ratio < 0)
    if np.any(ratio[np.logical_not(filter_condition_to_zero)]<0):
        pass
        #breakpoint()
    std_sq[filter_condition_to_zero]=0
    std = (
        std_sq
    )**0.5
    if np.any(std_sq < 0):
        A = ave_sq[std_sq<0]
        print(A)
        #breakpoint()
    return std

def calculate_std_and_save(lmp_path, chunk_tuple, variable_name):
    # calculate the std by read value and value_sq
    # save it in the same npy folder that read from
    ave = np.load(lmp_path + chunk_tuple[3] + variable_name + ".npy", mmap_mode='r')
    ave_sq = np.load(lmp_path + chunk_tuple[3] + variable_name + "_sq" + ".npy", mmap_mode='r')
    
    maskto0 = np.logical_and(
        ave_sq > -10**-50, ave_sq < 0,
    )
    ave_sq_revised = np.copy(ave_sq)
    ave_sq_revised[maskto0] = 0
    if np.any(ave_sq_revised<0):
        A = ave_sq_revised[np.logical_and(-10**-20 < ave_sq_revised, ave_sq_revised < 0)]
        print(A)
        B = ave_sq_revised[np.logical_and(-10**-6 < ave_sq_revised, ave_sq_revised < 0)]
        print(B)
        #breakpoint()
    
    #breakpoint()
    n_sample_keyword = chunk_tuple[4]
    n_sample = int(rr.logfile[n_sample_keyword])
    std = calculate_std_by_ave_and_sq_ave(n_sample, ave, ave_sq_revised)
    np.save(lmp_path + chunk_tuple[3] + variable_name + "_std" + ".npy", std)

def calculate_std_and_save_input_v_list(lmp_path, chunk_tuple, variable_name_list):
    for variable_name in variable_name_list:
        calculate_std_and_save(lmp_path, chunk_tuple, variable_name)

# https://en.wikipedia.org/wiki/Propagation_of_uncertainty
# propagation_of_std
def propagation_of_std_plus_or_minus(a, std_a, b, std_b):
    value = (
        (std_a)**2 + (std_b)**2
    )**0.5
    return value
def propagation_of_std_multi(a, std_a, b, std_b):
    value = a*b*(
        (std_a/a)**2 + (std_b/b)**2
    )**0.5
    return value
def propagation_of_std_divide(a, std_a, b, std_b):
    value = a/b*(
        (std_a/a)**2 + (std_b/b)**2
    )**0.5
    return value

bin_volume = (
    float(rr.logfile['width_wall_dp_unit'])
    *float(rr.logfile['bin_y_dp_unit_approximate'])
    *float(rr.logfile['bin_z_dp_unit_approximate'])
    *float(rr.logfile['dp'])**3
)

# data
def data_auto_nochunk_wall(lmp_path):
    # no chunk data
    for nochunk_wall_tuple in [inwall_force, outwall_force, zbottom_force]:
        data_not_chunk_transfer_text_to_npy(lmp_path, nochunk_wall_tuple[0], nochunk_wall_tuple[2])
    #### rename wall data file
    rename_dic = {
        "v_force_y_top_x": "v_force_outwall_1",
        "v_force_y_top_y": "v_force_outwall_2",
        "v_force_y_top_z": "v_force_outwall_3",
        "v_force_y_bottom_x": "v_force_inwall_1",
        "v_force_y_bottom_y": "v_force_inwall_2",
        "v_force_y_bottom_z": "v_force_inwall_3",
        "v_force_zbottom_x": "v_force_zbottom_1",
        "v_force_zbottom_y": "v_force_zbottom_2",
        "v_force_zbottom_z": "v_force_zbottom_3",
    }
    def rename_fun_for_wall(name_in, name_out):
        if os.path.exists(lmp_path + nochunk_npyfilepath_coord_char[name_in]['npyfilepath']):
            os.rename(
                lmp_path + nochunk_npyfilepath_coord_char[name_in]['npyfilepath'],
                lmp_path + nochunk_npyfilepath_coord_char[name_out]['npyfilepath'],
            )
    for key in rename_dic:
        rename_fun_for_wall(key, rename_dic[key])
    #### rename wall data file end

def data_auto(lmp_path):
    
    # chunk data
    chunk_tuple_list = [
        chunk_contact,
        chunk_mv_Ek_mass,
        chunk_omega,
        chunk_stress,
        chunk_inwall_force,
        chunk_outwall_force,
        chunk_zbottom_force,
    ]
    

    # save all variable in header to npy folder. Also save coord.
    for chunk_tuple in chunk_tuple_list:
        text_to_npy_for_chunk(lmp_path, chunk_tuple)
    # calculate std and save
    for (chunk_tuple, variable_name_list) in [
        (chunk_contact, ["n_contact"]),
        (chunk_mv_Ek_mass, ["v_mv1", "v_mv2", "v_mv3", "c_m1"]),
        (chunk_omega, ["c_omega[1]", "c_omega[2]", "c_omega[3]"]),
        (chunk_stress, ["c_stress[1]", "c_stress[2]", "c_stress[3]","c_stress[4]", "c_stress[5]", "c_stress[6]"]),
        (chunk_inwall_force, ["v_inwall_per_atom_1", "v_inwall_per_atom_2", "v_inwall_per_atom_3"]),
        (chunk_outwall_force, ["v_outwall_per_atom_1", "v_outwall_per_atom_2", "v_outwall_per_atom_3"]),
        (chunk_zbottom_force, ["v_zbottom_per_atom_1", "v_zbottom_per_atom_2", "v_zbottom_per_atom_3"]),
    ]:
        calculate_std_and_save_input_v_list(lmp_path, chunk_tuple, variable_name_list)
    
    # calculate wall stress and save
    eachbin2D_on_vertical_area = float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['bin_z_dp_unit_approximate'])*float(rr.logfile['dp'])**2
    eachbin2D_on_bottom_area = float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['bin_y_dp_unit_approximate'])*float(rr.logfile['dp'])**2
    
    for index_component_i in [1,2,3]:
        for wallstr, bin2D_area_nearwall in [
            ("inwall", eachbin2D_on_vertical_area),
            ("outwall", eachbin2D_on_vertical_area),
            ("zbottom", eachbin2D_on_bottom_area),
        ]:
            v_name = wallstr + "_force_" + str(index_component_i)
            wallforce = load_v_array(lmp_path, v_name)
            wallforce_std = load_v_array(lmp_path, v_name + "_std")
            wallstress = wallforce/bin2D_area_nearwall
            wallstress_std = wallforce_std/bin2D_area_nearwall
       
            #os.makedirs(lmp_path + "postprocess/npy/calculate/", exist_ok=True)
            save_v_array_to_disk(lmp_path, wallstr + "_stress_" + str(index_component_i), wallstress)
            save_v_array_to_disk(lmp_path, wallstr + "_stress_" + str(index_component_i) + "_std", wallstress_std)

    # calculate velocity and save
    for index_component_i in [1,2,3]:
        mv_i = load_v_array(lmp_path, "mv_" + str(index_component_i))
        mv_i_std = load_v_array(lmp_path, "mv_" + str(index_component_i) + "_std")
        mass = load_v_array(lmp_path, "mass")
        mass_std = load_v_array(lmp_path, "mass_std")
        velocity_i = mv_i/mass
        velocity_i_std = propagation_of_std_divide(mv_i, mv_i_std, mass, mass_std)
        #os.makedirs(lmp_path + "postprocess/npy/calculate/", exist_ok=True)
        save_v_array_to_disk(lmp_path, "velocity_" + str(index_component_i), velocity_i)
        save_v_array_to_disk(lmp_path, "velocity_" + str(index_component_i) + "_std", velocity_i_std)
        
    # calculate shear rate and save, revise shear rate and the grid by average
    for index_component_i in [1,2,3]:
        velocity_i = load_v_array(lmp_path, "velocity_" + str(index_component_i))
        velocity_i_std = load_v_array(lmp_path, "velocity_" + str(index_component_i) + "_std")
        Coord1 = load_coord(lmp_path, 'original_2D', "Coord1")
        Coord2 = load_coord(lmp_path, 'original_2D', "Coord2")
        # i_1
        strain_rate_2i = (velocity_i[:, :-1, :] - velocity_i[:, 1:, :])/(Coord1[:-1, :] - Coord1[ 1:, :])
        strain_rate_2i_std = propagation_of_std_plus_or_minus(
            velocity_i[:, :-1, :],
            velocity_i_std[:, :-1, :],
            velocity_i[:, 1:, :],
            velocity_i_std[:, 1:, :],
            )/np.abs(
                (Coord1[:-1, :] - Coord1[ 1:, :])
            )
        Coord1_new = (Coord1[:-1, :] + Coord1[ 1:, :])/2
        Coord2_new = (Coord2[:-1, :] + Coord2[ 1:, :])/2
        save_folder_rela_lmp = "postprocess/npy/calculate/" + "strain_rate_" + "2" + str(index_component_i) + "/"
        #os.makedirs(lmp_path + save_folder_rela_lmp, exist_ok=True)
        save_v_array_to_disk(lmp_path, "strain_rate_2" + str(index_component_i), strain_rate_2i)
        save_v_array_to_disk(lmp_path, "strain_rate_2" + str(index_component_i) + "_std", strain_rate_2i_std)
        save_coord_to_disk(lmp_path, c_r_npyfilepath_coord_char["strain_rate_2" + str(index_component_i)]["coordinate_characteristic"], "Coord1", Coord1_new)
        save_coord_to_disk(lmp_path, c_r_npyfilepath_coord_char["strain_rate_2" + str(index_component_i)]["coordinate_characteristic"], "Coord2", Coord2_new)
        
        # average to middle point
        strain_rate_2i_middle = (strain_rate_2i[:, :, :-1] + strain_rate_2i[:, :, 1:])/2
        strain_rate_2i_middle_std = propagation_of_std_plus_or_minus(
            strain_rate_2i[:, :, :-1],
            strain_rate_2i_std[:, :, :-1],
            strain_rate_2i[:, :, 1:],
            strain_rate_2i_std[:, :, 1:],
            )/2
        Coord1_new_middle = (Coord1_new[:, :-1] + Coord1_new[:, 1:])/2
        Coord2_new_middle = (Coord2_new[:, :-1] + Coord2_new[:, 1:])/2
        save_folder_rela_lmp = "postprocess/npy/calculate/" + "middle_strain_rate/"
        #os.makedirs(lmp_path + save_folder_rela_lmp, exist_ok=True)
        save_v_array_to_disk(lmp_path, "strain_rate_2" + str(index_component_i) + "_middle", strain_rate_2i_middle)
        save_v_array_to_disk(lmp_path, "strain_rate_2" + str(index_component_i) + "_middle" + "_std", strain_rate_2i_middle_std)
        save_coord_to_disk(lmp_path, 
            c_r_npyfilepath_coord_char["strain_rate_2" + str(index_component_i) + "_middle"]["coordinate_characteristic"],
            "Coord1",
            Coord1_new_middle,
        )
        save_coord_to_disk(lmp_path, 
            c_r_npyfilepath_coord_char["strain_rate_2" + str(index_component_i) + "_middle"]["coordinate_characteristic"],
            "Coord2",
            Coord2_new_middle,
        )
        # i_2
        strain_rate_3i = (velocity_i[:, :, :-1] - velocity_i[:, :, 1:])/(Coord2[:, :-1] - Coord2[:, 1:])
        strain_rate_3i_std = propagation_of_std_plus_or_minus(
            velocity_i[:, :, :-1],
            velocity_i_std[:, :, :-1],
            velocity_i[:, :, 1:],
            velocity_i_std[:, :, 1:],
            )/np.abs(
                (Coord2[:, :-1] - Coord2[:, 1:])
            )
        Coord1_new = (Coord1[:, :-1] + Coord1[:, 1:])/2
        Coord2_new = (Coord2[:, :-1] + Coord2[:, 1:])/2
        save_folder_rela_lmp = "postprocess/npy/calculate/" + "strain_rate_" + "3" + str(index_component_i) + "/"
        #os.makedirs(lmp_path + save_folder_rela_lmp, exist_ok=True)
        save_v_array_to_disk(lmp_path, "strain_rate_3" + str(index_component_i), strain_rate_3i)
        save_v_array_to_disk(lmp_path, "strain_rate_3" + str(index_component_i) + "_std", strain_rate_3i_std)
        save_coord_to_disk(lmp_path, c_r_npyfilepath_coord_char["strain_rate_3" + str(index_component_i)]["coordinate_characteristic"], "Coord1", Coord1_new)
        save_coord_to_disk(lmp_path, c_r_npyfilepath_coord_char["strain_rate_3" + str(index_component_i)]["coordinate_characteristic"], "Coord2", Coord2_new)
        # average to middle point
        # average to middle point
        strain_rate_3i_middle = (strain_rate_3i[:, :-1, :] + strain_rate_3i[:, 1:, :])/2
        strain_rate_3i_middle_std = propagation_of_std_plus_or_minus(
            strain_rate_3i[:, :-1, :],
            strain_rate_3i_std[:, :-1, :],
            strain_rate_3i[:, 1:, :],
            strain_rate_3i_std[:, 1:, :],
            )/2
        save_folder_rela_lmp = "postprocess/npy/calculate/" + "middle_strain_rate/"
        #os.makedirs(lmp_path + save_folder_rela_lmp, exist_ok=True)
        save_v_array_to_disk(lmp_path, "strain_rate_3" + str(index_component_i) + "_middle", strain_rate_3i_middle)
        save_v_array_to_disk(lmp_path, "strain_rate_3" + str(index_component_i) + "_middle" + "_std", strain_rate_3i_middle_std)
        
        # check if same coord
        Coord1_new_middle_again = (Coord1_new[:-1, :] + Coord1_new[1:, :])/2
        Coord2_new_middle_again = (Coord2_new[:-1, :] + Coord2_new[1:, :])/2
        if np.any(Coord1_new_middle_again != Coord1_new_middle):
            print("not equal coord")
            breakpoint()        

    # pressure
    bin_volume = (
        float(rr.logfile['width_wall_dp_unit'])
        *float(rr.logfile['bin_y_dp_unit_approximate'])
        *float(rr.logfile['bin_z_dp_unit_approximate'])
        *float(rr.logfile['dp'])**3
    )
    stress_11 = load_v_array(lmp_path, "stress_multiply_binvolume_11")/bin_volume
    stress_22 = load_v_array(lmp_path, "stress_multiply_binvolume_22")/bin_volume
    stress_33 = load_v_array(lmp_path, "stress_multiply_binvolume_33")/bin_volume
    stress_11_std = load_v_array(lmp_path, "stress_multiply_binvolume_11" + "_std")/bin_volume
    stress_22_std = load_v_array(lmp_path, "stress_multiply_binvolume_22" + "_std")/bin_volume
    stress_33_std = load_v_array(lmp_path, "stress_multiply_binvolume_33" + "_std")/bin_volume
    save_v_array_to_disk(lmp_path, "stress_11", stress_11)
    save_v_array_to_disk(lmp_path, "stress_11" + "_std", stress_11_std)
    save_v_array_to_disk(lmp_path, "stress_22", stress_22)
    save_v_array_to_disk(lmp_path, "stress_22" + "_std", stress_22_std)
    save_v_array_to_disk(lmp_path, "stress_33", stress_33)
    save_v_array_to_disk(lmp_path, "stress_33" + "_std", stress_33_std)
    pressure = -1/3*(stress_11 + stress_22 + stress_33)
    pressure_std = 1/3*(stress_11_std**2 + stress_22_std**2 + stress_33_std**2)**0.5
    #os.makedirs(lmp_path + "postprocess/npy/calculate/", exist_ok=True)
    np.save(lmp_path + "postprocess/npy/calculate/" + "pressure" + ".npy", pressure)
    np.save(lmp_path + "postprocess/npy/calculate/" + "pressure" + "_std" + ".npy", pressure_std)
    # stress 12 13 23
    stress_12 = load_v_array(lmp_path, "stress_multiply_binvolume_12")/bin_volume
    stress_13 = load_v_array(lmp_path, "stress_multiply_binvolume_13")/bin_volume
    stress_23 = load_v_array(lmp_path, "stress_multiply_binvolume_23")/bin_volume
    stress_12_std = load_v_array(lmp_path, "stress_multiply_binvolume_12" + "_std")/bin_volume
    stress_13_std = load_v_array(lmp_path, "stress_multiply_binvolume_13" + "_std")/bin_volume
    stress_23_std = load_v_array(lmp_path, "stress_multiply_binvolume_23" + "_std")/bin_volume
    np.save(lmp_path + "postprocess/npy/calculate/" + "stress_12" + ".npy", stress_12)
    np.save(lmp_path + "postprocess/npy/calculate/" + "stress_12_std" + ".npy", stress_12_std)
    np.save(lmp_path + "postprocess/npy/calculate/" + "stress_13" + ".npy", stress_13)
    np.save(lmp_path + "postprocess/npy/calculate/" + "stress_13_std" + ".npy", stress_13_std)
    np.save(lmp_path + "postprocess/npy/calculate/" + "stress_23" + ".npy", stress_23)
    np.save(lmp_path + "postprocess/npy/calculate/" + "stress_23_std" + ".npy", stress_23_std)
    # shear rate symmetry part
    # 01 10
    strain_rate_12_middle = 0
    strain_rate_13_middle = 0
    strain_rate_12_middle_std = 0
    strain_rate_13_middle_std = 0

    strain_rate_21_middle = load_v_array(lmp_path, "strain_rate_21_middle")
    strain_rate_21_middle_std = load_v_array(lmp_path, "strain_rate_21_middle_std") 
    strain_rate_22_middle = load_v_array(lmp_path, "strain_rate_22_middle")
    strain_rate_22_middle_std = load_v_array(lmp_path, "strain_rate_22_middle_std") 
    strain_rate_23_middle = load_v_array(lmp_path, "strain_rate_23_middle")
    strain_rate_23_middle_std = load_v_array(lmp_path, "strain_rate_23_middle_std") 
    strain_rate_31_middle = load_v_array(lmp_path, "strain_rate_31_middle")
    strain_rate_31_middle_std = load_v_array(lmp_path, "strain_rate_31_middle_std") 
    strain_rate_32_middle = load_v_array(lmp_path, "strain_rate_32_middle")
    strain_rate_32_middle_std = load_v_array(lmp_path, "strain_rate_32_middle_std") 
    strain_rate_33_middle = load_v_array(lmp_path, "strain_rate_33_middle")
    strain_rate_33_middle_std = load_v_array(lmp_path, "strain_rate_33_middle_std") 

    symmetry_strain_rate_12 = strain_rate_12_middle + strain_rate_21_middle
    symmetry_strain_rate_12_std = propagation_of_std_plus_or_minus(strain_rate_12_middle, strain_rate_12_middle_std, strain_rate_21_middle, strain_rate_21_middle_std)
    symmetry_strain_rate_13 = strain_rate_13_middle + strain_rate_31_middle
    symmetry_strain_rate_13_std = propagation_of_std_plus_or_minus(strain_rate_13_middle, strain_rate_13_middle_std, strain_rate_31_middle, strain_rate_31_middle_std)
    symmetry_strain_rate_23 = strain_rate_23_middle + strain_rate_32_middle
    symmetry_strain_rate_23_std = propagation_of_std_plus_or_minus(strain_rate_23_middle, strain_rate_23_middle_std, strain_rate_32_middle, strain_rate_32_middle_std)
    
    # shear rate second invariant
    strain_rate_second_invariant = (
        symmetry_strain_rate_12**2 + symmetry_strain_rate_13**2 + symmetry_strain_rate_23**2
    )**0.5
    strain_rate_second_invariant_std = (
        0.25*(
            4*symmetry_strain_rate_12**2*symmetry_strain_rate_12_std**2
            + 4*symmetry_strain_rate_13**2*symmetry_strain_rate_13_std**2
            + 4*symmetry_strain_rate_23**2*symmetry_strain_rate_23_std**2
        )/strain_rate_second_invariant
    )**0.5

    # symmetry shear rate ratio = strain_rate ij / strain_rate second invariant
    symmetry_strain_rate_ratio_12 = symmetry_strain_rate_12/strain_rate_second_invariant
    symmetry_strain_rate_ratio_13 = symmetry_strain_rate_13/strain_rate_second_invariant
    symmetry_strain_rate_ratio_23 = symmetry_strain_rate_23/strain_rate_second_invariant
    symmetry_strain_rate_ratio_12_std = propagation_of_std_divide(
        symmetry_strain_rate_12, symmetry_strain_rate_12_std, strain_rate_second_invariant, strain_rate_second_invariant_std,
        )
    symmetry_strain_rate_ratio_13_std = propagation_of_std_divide(
        symmetry_strain_rate_13, symmetry_strain_rate_13_std, strain_rate_second_invariant, strain_rate_second_invariant_std,
        )
    symmetry_strain_rate_ratio_23_std = propagation_of_std_divide(
        symmetry_strain_rate_23, symmetry_strain_rate_23_std, strain_rate_second_invariant, strain_rate_second_invariant_std,
        )

    # stress ratio 12 13 23 21 31 32
    stress_ratio_12 = stress_12/pressure
    #breakpoint()
    stress_ratio_13 = stress_13/pressure
    stress_ratio_23 = stress_23/pressure
    stress_ratio_12_std = propagation_of_std_divide(stress_12, stress_12_std, pressure, pressure_std)
    stress_ratio_13_std = propagation_of_std_divide(stress_13, stress_13_std, pressure, pressure_std)
    stress_ratio_23_std = propagation_of_std_divide(stress_23, stress_23_std, pressure, pressure_std)

    # mu = stress ratio (/ shear rate ratio)
    [mu_12, mu_13, mu_23] = [stress_ratio_12, stress_ratio_13, stress_ratio_23]
    [mu_12_std, mu_13_std, mu_23_std] = [stress_ratio_12_std, stress_ratio_13_std, stress_ratio_23_std]
    save_v_array_to_disk(lmp_path, "mu_12", mu_12)
    save_v_array_to_disk(lmp_path, "mu_12_std", mu_12_std)
    save_v_array_to_disk(lmp_path, "mu_13", mu_13)
    save_v_array_to_disk(lmp_path, "mu_13_std", mu_13_std)
    save_v_array_to_disk(lmp_path, "mu_23", mu_23)
    save_v_array_to_disk(lmp_path, "mu_23_std", mu_23_std)
    def ave_4_grid_for_the_last_axis(value_array):
        value_array_middle = (value_array[:, :-1, :-1] + value_array[:, 1:, :-1] + value_array[:, :-1, 1:] + value_array[:, 1:, 1:])/4
        return value_array_middle
    def ave_4_grid_for_the_last_axis_std(value_array_std):
        value_array_std_sq = value_array_std**2
        value_array_middle_std = (value_array_std_sq[:, :-1, :-1] + value_array_std_sq[:, 1:, :-1] + value_array_std_sq[:, :-1, 1:] + value_array_std_sq[:, 1:, 1:])**0.5/4
        return value_array_middle_std
    
    [mu_12_middle, mu_13_middle, mu_23_middle] = [
        ave_4_grid_for_the_last_axis(mu_12),
        ave_4_grid_for_the_last_axis(mu_13),
        ave_4_grid_for_the_last_axis(mu_23),
    ]
    [mu_12_middle_std, mu_13_middle_std, mu_23_middle_std] = [
        ave_4_grid_for_the_last_axis_std(mu_12_std),
        ave_4_grid_for_the_last_axis_std(mu_13_std),
        ave_4_grid_for_the_last_axis_std(mu_23_std),
    ]
    save_folder_rela_lmp = "postprocess/npy/calculate/" + "middle_mu_tensor/"
    #os.makedirs(lmp_path + save_folder_rela_lmp, exist_ok=True)
    save_v_array_to_disk(lmp_path, "mu_12_middle", mu_12_middle)
    save_v_array_to_disk(lmp_path, "mu_12_middle_std", mu_12_middle_std)
    save_v_array_to_disk(lmp_path, "mu_13_middle", mu_13_middle)
    save_v_array_to_disk(lmp_path, "mu_13_middle_std", mu_13_middle_std)
    save_v_array_to_disk(lmp_path, "mu_23_middle", mu_23_middle)
    save_v_array_to_disk(lmp_path, "mu_23_middle_std", mu_23_middle_std)
    
    mu_tensor_12 = mu_12_middle/symmetry_strain_rate_ratio_12
    mu_tensor_13 = mu_13_middle/symmetry_strain_rate_ratio_13
    mu_tensor_23 = mu_23_middle/symmetry_strain_rate_ratio_23
    mu_tensor_12_std = propagation_of_std_divide(mu_12_middle, mu_12_middle_std, symmetry_strain_rate_ratio_12, symmetry_strain_rate_ratio_12_std)
    mu_tensor_13_std = propagation_of_std_divide(mu_13_middle, mu_13_middle_std, symmetry_strain_rate_ratio_13, symmetry_strain_rate_ratio_13_std)
    mu_tensor_23_std = propagation_of_std_divide(mu_23_middle, mu_23_middle_std, symmetry_strain_rate_ratio_23, symmetry_strain_rate_ratio_23_std)

    save_v_array_to_disk(lmp_path, "mu_tensor_12", mu_tensor_12)
    save_v_array_to_disk(lmp_path, "mu_tensor_12_std", mu_tensor_12_std)
    save_v_array_to_disk(lmp_path, "mu_tensor_13", mu_tensor_13)
    save_v_array_to_disk(lmp_path, "mu_tensor_13_std", mu_tensor_13_std)
    save_v_array_to_disk(lmp_path, "mu_tensor_23", mu_tensor_23)
    save_v_array_to_disk(lmp_path, "mu_tensor_23_std", mu_tensor_23_std)
    # I = symmetry strain_rate(second_invariant) d / sqrt(P/density of particle)
    diameter = float(rr.logfile["dp"])
    density = float(rr.logfile['den'])
    pressure_sqrt_std = (0.25*pressure_std**2/pressure)**0.5
    pressure_middle = ave_4_grid_for_the_last_axis(pressure)
    pressure_sqrt_middle_std = ave_4_grid_for_the_last_axis_std(pressure_sqrt_std)
    I_12 = symmetry_strain_rate_12*diameter/(pressure_middle/density)**0.5
    I_13 = symmetry_strain_rate_13*diameter/(pressure_middle/density)**0.5
    I_23 = symmetry_strain_rate_23*diameter/(pressure_middle/density)**0.5
    I_12_std = propagation_of_std_divide(symmetry_strain_rate_12, symmetry_strain_rate_12_std, pressure_middle**0.5, pressure_sqrt_middle_std)*diameter*density**0.5
    I_13_std = propagation_of_std_divide(symmetry_strain_rate_13, symmetry_strain_rate_13_std, pressure_middle**0.5, pressure_sqrt_middle_std)*diameter*density**0.5
    I_23_std = propagation_of_std_divide(symmetry_strain_rate_23, symmetry_strain_rate_23_std, pressure_middle**0.5, pressure_sqrt_middle_std)*diameter*density**0.5
    I_tensor = strain_rate_second_invariant*diameter/(pressure_middle/density)**0.5
    I_tensor_std = propagation_of_std_divide(strain_rate_second_invariant, strain_rate_second_invariant_std, pressure_middle**0.5, pressure_sqrt_middle_std)*diameter*density**0.5
    save_v_array_to_disk(lmp_path, "I_12", I_12)
    save_v_array_to_disk(lmp_path, "I_12_std", I_12_std)
    save_v_array_to_disk(lmp_path, "I_13", I_13)
    save_v_array_to_disk(lmp_path, "I_13_std", I_13_std)
    save_v_array_to_disk(lmp_path, "I_23", I_23)
    save_v_array_to_disk(lmp_path, "I_23_std", I_23_std)
    save_v_array_to_disk(lmp_path, "I_tensor", I_tensor)
    save_v_array_to_disk(lmp_path, "I_tensor_std", I_tensor_std)
    #breakpoint()
    # calculate fraction and save
    mass = load_v_array(lmp_path, "mass")
    mass_std = load_v_array(lmp_path, "mass_std")
    fraction = mass/float(rr.logfile['den'])/vol_in_chunks
    fraction_std = mass_std/float(rr.logfile['den'])/vol_in_chunks
    np.save(lmp_path + "postprocess/npy/calculate/" + "fraction" + ".npy", fraction)
    np.save(lmp_path + "postprocess/npy/calculate/" + "fraction_std" + ".npy", fraction_std)

# main exclusive
def main():
    data_auto_nochunk_wall(lmp_path = rr.folder_path_list_initial_to_last[-1])
    data_auto_nochunk_wall(lmp_path = rr.folder_path_list_initial_to_last[-2])
    data_auto_nochunk_wall(lmp_path = rr.folder_path_list_initial_to_last[-3])
    data_auto(lmp_path = rr.folder_path_list_initial_to_last[-1])
    data_auto(lmp_path = rr.folder_path_list_initial_to_last[-2])
    data_auto(lmp_path = rr.folder_path_list_initial_to_last[-3])
    print('\a')
# main exclusive
if __name__ == "__main__":
    
    main()