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

# setting matplotlib
# Here are some matplotlib style settings that I (JSW) find useful for publication-quality plots
# http://tobackgroup.physics.tamu.edu/useful/matplotlib-formatting.txt
# First, we set the width of a figure to be exactly the width of a single column in PRL or PRD (in
# inches), and set the aspect ratio to be 4:3.

plt.rc('figure', figsize=(3.375, 3.375*(.75)), titlesize=10)

# This line sets the resolution (number of dots per inch) for saving figures to raster formats like
# PNG.  It's mainly useful to set the size of the figures displayed in the webbrowser in the Jupyter
# or IPython notebook
plt.rc('savefig', dpi=180)
#plt.rc('axes.formatter', limits=(0,0), use_mathtext=True)
plt.rc('savefig', bbox='tight')
plt.ticklabel_format(style='sci', axis='both')

# The default for scatter plots is to put three points in each legend entry.  Usually one is
# sufficient.
#plt.rc('legend', scatterpoints=1, fontsize=10, frameon=True)

# We would like the text to be created using TeX, so we get full support for LaTeX math syntax.
#plt.rc('text', usetex=True)

# We want all the text to be Computer Modern, and 10pt.  This way, combined with the correct setting
# of the size of the figure, the text in the figure will exactly match the text in the rest of the
# manuscript.

plt.rc('font', size=8, weight='normal', family='sans-serif', serif=['Helvetica', 'Arial'])

plt.rc('axes', labelsize=10, titlesize=8, labelweight='normal')

# http://aeturrell.com/2018/01/31/publication-quality-plots-in-python/
#plt.rc('xtick', labelsize=10)
#plt.rc('ytick', labelsize=10)
plt.rc('figure', autolayout=False)
plt.rc('lines', linewidth=2)
plt.rc('lines', markersize=4)
plt.rc('mathtext', fontset="stix")



# more setting
# http://physicalmodelingwithpython.blogspot.com/2015/06/making-plots-for-publication.html

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
    "TimeStep v_t v_force_inwall_1 v_force_inwall_2 v_force_inwall_3",
)
outwall_force = (
    "output/wall/force_y_top_to_particle.allstep",
    "TimeStep v_t v_force_outwall_1 v_force_outwall_2 v_force_outwall_3",
)
zbottom_force = (
    "output/wall/force_zbottom_to_particle.allstep",
    "TimeStep v_t v_force_zbottom_1 v_force_zbottom_2 v_force_zbottom_3",
)

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
    os.makedirs(outputfolder_rela_lmp_path, exist_ok=True)
    
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
    os.makedirs(outputfolder_rela_lmp_path, exist_ok=True)
    
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
    header = first_4_lines[1].split()[1:]
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
    os.makedirs(outputfolder_rela_lmp_path, exist_ok=True)
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

def calculate_std_and_save(chunk_tuple, variable_name):
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

def calculate_std_and_save_input_v_list(chunk_tuple, variable_name_list):
    for variable_name in variable_name_list:
        calculate_std_and_save(chunk_tuple, variable_name)

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

def time_from_step_0(step):
    return step*float(rr.logfile["ts"])

def time_from_start_rotate(step):
    return time_from_step_0(step)-rr.logfile["rotate_start_time"]

def strain_from_rotate_start(step):
    shear_rate = float(rr.logfile['in_velocity'])/(float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp']))
    return time_from_start_rotate(step)*shear_rate

def transfer_time_to_str(timearray):
    if isinstance(timearray, str):
        if timearray == 'all':
            string = 'all'
    else:
        string = "-".join([
            "{:.2e}".format(a) for a in [timearray[0], timearray[-1]]
        ])
    return string

def transfer_coor_to_str(coordarray, ifave=False):
    if isinstance(coordarray, str):
        if coordarray == 'all':
            string = 'all'
        elif coordarray == 'ave':
            string = 'ave'
        elif coordarray == 'sumdividemass':
            string = 'sumdividemass'
    else:
        if ifave:
            string ="ave" + "-".join([
                "{:d}".format(a) for a in [coordarray[0], coordarray[-1]]
            ])
        else:
            string ="-".join([
                "{:d}".format(a) for a in [coordarray[0], coordarray[-1]]
            ])
    return string

def get_value(lmp_path, name, n_ave, inputstepsarray):
    if n_ave % 2 == 1:
        pass
    else:
        sys.exit("n_ave should be odd integer")
    # combine with previous
    # pick the smallest step and max step
    # pick the smallest step and largest step in this simulation
    
    def get_step_first_end(lmp_path_folder):
        timestep_local = np.load(lmp_path_folder + c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        
        timestep_first = timestep_local[int((n_ave-1)/2)]
        timestep_last = timestep_local[-1-int((n_ave-1)/2)]
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
    timestep = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
    
    array = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + c_r_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
    for j in range(n_include_pre_simu-1):
        timestep_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        array_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + c_r_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
        timestep = np.append(timestep, timestep_append, axis=0)
        array = np.append(array, array_append, axis=0)
    d_step = timestep[1] - timestep[0]
    indexesarray = (inputstepsarray-timestep[0])/d_step
    indexesarray = indexesarray.astype(int)
    # average
    sum_over_n_array = 0
    for k in range(-n_ave, n_ave+1):
        sum_over_n_array = sum_over_n_array + array[indexesarray + k]
    ave_array = sum_over_n_array/(2*n_ave + 1)
    return ave_array

# data
def data_auto(lmp_path):
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
        calculate_std_and_save_input_v_list(chunk_tuple, variable_name_list)
    
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

lmp_path = rr.folder_path_list_initial_to_last[-1]

# check if all timestep the same for chunk
for key in map_chunkfile_char_save_folderpath.keys():
    path = lmp_path + map_chunkfile_char_save_folderpath[key]["timestep_path"]
    if np.any(
        np.load(path, mmap_mode='r') != np.load(lmp_path + map_chunkfile_char_save_folderpath["mv_Ek_mass"]["timestep_path"], mmap_mode='r')
    ):
        print("timestep file not match for all chunk file")
        breakpoint()

# scale
stress_scale = float(rr.logfile['den'])*float(rr.logfile['g'])*float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp'])
stress_scale_str = r'$\rho_s g h $'
velocity_scale = float(rr.logfile['in_velocity'])
if velocity_scale < 0:
    velocity_scale = -velocity_scale
velocity_scale_str = r'$V_{inwall}$'
strain_rate_scale = velocity_scale/(float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp']))
strain_rate_scale_str = r'$V_{inwall}/Width$'
coord_scale = float(rr.logfile['dp'])
coord_scale_str = r'$d_p$'
mu_scale = -1
mu_scale_str = None #r'$-1$'
I_scale = -1
I_scale_str = r'$-1$'
mu_tensor_scale = 1
mu_tensor_scale_str = None
I_tensor_scale = 1
I_tensor_scale_str = None
coord_scale = float(rr.logfile["dp"])
coord_scale_str = r'$d_p$'
# variable_name to string in pic label or legend
v_name_to_labal_str = {
    "pressure":"P",
    "n_contact":"Z",
    "velocity_1": "V",
    "velocity_2": "V",
    "velocity_3": "V",
    "strain": r'$\gamma$',
    "inwall_force_1": r'$F_x$',
    "inwall_force_2": r'$F_y$',
    "inwall_force_3": r'$F_z$',
    "outwall_force_1": r'$F_x$',
    "outwall_force_2": r'$F_y$',
    "outwall_force_3": r'$F_z$',
    "zbottom_force_1": r'$F_x$',
    "zbottom_force_2": r'$F_y$',
    "zbottom_force_3": r'$F_z$',
    "strain_rate_21": r'$\dot{\gamma_{21}}$',
    "strain_rate_22": r'$\dot{\gamma_{22}}$',
    "strain_rate_23": r'$\dot{\gamma_{23}}$',
    "strain_rate_31": r'$\dot{\gamma_{31}}$',
    "strain_rate_32": r'$\dot{\gamma_{32}}$',
    "strain_rate_33": r'$\dot{\gamma_{33}}$',
    "strain_rate_21_middle": r'$\dot{\gamma_{21}}$',
    "strain_rate_22_middle": r'$\dot{\gamma_{22}}$',
    "strain_rate_23_middle": r'$\dot{\gamma_{23}}$',
    "strain_rate_31_middle": r'$\dot{\gamma_{31}}$',
    "strain_rate_32_middle": r'$\dot{\gamma_{32}}$',
    "strain_rate_33_middle": r'$\dot{\gamma_{33}}$',
    "stress_11": r'$\sigma_{11}$',
    "stress_22": r'$\sigma_{22}$',
    "stress_33": r'$\sigma_{33}$',
    "stress_12": r'$\sigma_{12}$',
    "stress_13": r'$\sigma_{13}$',
    "stress_23": r'$\sigma_{23}$',
    "mu_12": r'$\mu_{12}$',
    "mu_13": r'$\mu_{13}$',
    "mu_23": r'$\mu_{23}$',
    "mu_12_middle": r'$\mu_{12}$',
    "mu_13_middle": r'$\mu_{13}$',
    "mu_23_middle": r'$\mu_{23}$',
    "mu_tensor_12": r'$\mu_{12}$',
    "mu_tensor_13": r'$\mu_{13}$',
    "mu_tensor_23": r'$\mu_{23}$',
    "I_12": r'$I_{12}$',
    "I_13": r'$I_{13}$',
    "I_23": r'$I_{23}$',
    "I_tensor": r'$I$',
    "fraction": r'$\phi$',
    "inwall_stress_1": r'$\sigma_{21}$',
    "inwall_stress_2": r'$\sigma_{22}$',
    "inwall_stress_3": r'$\sigma_{23}$',
    "outwall_stress_1": r'$\sigma_{21}$',
    "outwall_stress_2": r'$\sigma_{22}$',
    "outwall_stress_3": r'$\sigma_{23}$',
    "zbottom_stress_1": r'$\sigma_{31}$',
    "zbottom_stress_2": r'$\sigma_{32}$',
    "zbottom_stress_3": r'$\sigma_{33}$',
    "Coord1": r'y',
    "Coord2": r'z',
}

def transfer_v_name_to_label_str_in_figure(v_name):
    if v_name in v_name_to_labal_str.keys():
        v_name_label_in_figure = v_name_to_labal_str[v_name]
    else:
        v_name_label_in_figure = v_name
    return v_name_label_in_figure

# sum mean select (plus minus multiply divide) operation for variable coord1 coord2 time
def sum_variable(variable):

    return variable

def plot_1D_from_chunk2D(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray, coord1_index_array, coord2_index_array,
        ave_over_coord1=False, ave_over_coord2=False,
        figure_class=None, legend_class=None, spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        x_scale = 'linear', y_scale = 'linear',
        useerrorbar = True,
        ifdivideNcount = False,
        maskstatic=None,
        masknonstatic=None,
    ):
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    # str for variable used in figure label legend
    
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')

    if coord2_index_array=="sumdividemass" and y_name=="mv_1":
        y_name_label_in_figure = "V"
    # legend_class: time, coord1, coord2
    # spaceave: coord1, coord2
    char = c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord1 = load_coord(lmp_path, char, "Coord1")
    Coord2 = load_coord(lmp_path, char, "Coord2")

    time_array = time_from_start_rotate(inputstepsarray)
    if x_name == "timestep":
        x_value = inputstepsarray
        x_value_std = 0
    elif x_name == "time":
        x_value = time_array
        x_value_std = 0
    elif x_name == "strain":
        x_value = strain_from_rotate_start(inputstepsarray)
        x_value_std = 0
    elif x_name == "Coord1":
        x_value = Coord1
        x_value_std = 0
    else:
        x_value = get_value(lmp_path, x_name, n_ave, inputstepsarray)
        x_value_std = get_value(lmp_path, x_name + "_std", n_ave, inputstepsarray)

    if y_name == "timestep":
        y_value = inputstepsarray
        y_value_std = 0
    elif y_name == "time":
        y_value = time_array
        y_value_std = 0
    else:
        y_value = get_value(lmp_path, y_name, n_ave, inputstepsarray)
        y_value_std = get_value(lmp_path, y_name + "_std", n_ave, inputstepsarray)
    if ifdivideNcount:
        Ncount_value = get_value(lmp_path, 'Ncount', n_ave, inputstepsarray)
    if coord2_index_array=="sumdividemass":
        mass_value = get_value(lmp_path, 'mass', n_ave, inputstepsarray)
    # scale factor
    x_value = x_value/x_scale_factor
    x_value_std = x_value_std/x_scale_factor
    y_value = y_value/y_scale_factor
    y_value_std = y_value_std/y_scale_factor

    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    if ifdivideNcount:
        subfoldername = y_name + "_divideN_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    if useerrorbar:
        os.makedirs(diagram_path_add_nve + subfoldername + "errorbar/", exist_ok=True)

    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    # title

    if x_scale_str is None:
        x_label_str = x_name_label_in_figure
    else:
        x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
    ax.set_xlabel(x_label_str)
    
    if y_scale_str is None:
        y_label_str = y_name_label_in_figure
    else:
        y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
    ax.set_ylabel(y_label_str)
    if ifdivideNcount and y_label_str=="n_contact":
        ax.set_ylabel("<Z>")
    
    
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    
    # legend_class: time, coord1, coord2
    if legend_class == 'tc1c2':
        for indexstep, step in enumerate(inputstepsarray):
            for coord1_index in coord1_index_array:
                for coord2_index in coord2_index_array:
                    time = time_from_start_rotate(step)
                    strain = strain_from_rotate_start(step)
                    label = ", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                        "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                        "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[indexstep]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[indexstep, coord1_index, coord2_index]
                        x_value_std_plot = x_value_std[indexstep, coord1_index, coord2_index]
                    
                    y_value_plot = y_value[indexstep, coord1_index, coord2_index]
                    y_value_std_plot = y_value_std[indexstep, coord1_index, coord2_index]
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[indexstep, coord1_index, coord2_index]
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c1c2':
        for coord1_index in coord1_index_array:
            for coord2_index in coord2_index_array:
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                ])

                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index, coord2_index]
                    x_value_std_plot = x_value_std[:, coord1_index, coord2_index]
                
                y_value_plot = y_value[:, coord1_index, coord2_index]
                y_value_std_plot = y_value_std[:, coord1_index, coord2_index]
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index, coord2_index]
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    elif legend_class == 'c1':
        if coord2_index_array=='ave':
            if ave_over_coord1:
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :]
                    x_value_plot = np.nanmean(x_value_plot, axis=(1,2))
                    x_value_std_plot = np.nanmean(x_value_std_plot, axis=(1,2))
                y_value_plot = y_value[:, coord1_index_array, :]
                y_value_std_plot = y_value_std[:, coord1_index_array, :]
                y_value_plot = np.nanmean(y_value_plot, axis=(1,2))
                y_value_std_plot = np.nanmean(y_value_std_plot, axis=(1,2))
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "({})".format(coord_scale_str),
                    "z=average",
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord1_index in coord1_index_array:
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y={:.1f} ({})".format(Coord1[coord1_index, 0]/coord_scale, coord_scale_str),
                        "z=average",
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, coord1_index, :]
                        x_value_std_plot = x_value_std[:, coord1_index, :]
                        #breakpoint()
                        x_value_plot = np.nanmean(x_value_plot, axis=1)
                        x_value_std_plot = np.nanmean(x_value_std_plot, axis=1)
                    y_value_plot = y_value[:, coord1_index, :]
                    y_value_std_plot = y_value_std[:, coord1_index, :]
                    y_value_plot = np.nanmean(y_value_plot, axis=1)
                    y_value_std_plot = np.nanmean(y_value_std_plot, axis=1)
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, coord1_index, :]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=1)
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c2':
        if coord1_index_array=='ave':
            if ave_over_coord2:
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot, axis=(1,2))
                    x_value_std_plot = np.nanmean(x_value_std_plot, axis=(1,2))
                y_value_plot = y_value[:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot, axis=(1,2))
                y_value_std_plot = np.nanmean(y_value_std_plot, axis=(1,2))
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),"y=average",
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord2_index in coord2_index_array:
                    
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y=average",
                        "z={:.1f} ({})".format(Coord2[0, coord2_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, :, coord2_index]
                        x_value_std_plot = x_value_std[:, :, coord2_index]
                        x_value_plot = np.nanmean(x_value_plot, axis=1)
                        x_value_std_plot = np.nanmean(x_value_std_plot, axis=1)
                    y_value_plot = y_value[:, :, coord2_index]
                    y_value_std_plot = y_value_std[:, :, coord2_index]
                    y_value_plot = np.nanmean(y_value_plot, axis=1)
                    y_value_std_plot = np.nanmean(y_value_std_plot, axis=1)
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, :, coord2_index]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=1)
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 't':
        if coord2_index_array=="sumdividemass":
            if x_name == "Coord1":
                for index_step, step in enumerate(inputstepsarray):
                    strain = strain_from_rotate_start(step)
                    label=", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                    ])
                    x_value_plot = x_value[coord1_index_array, 0]
                    x_value_std_plot = 0
                    y_value_plot = y_value[index_step, coord1_index_array, :]
                    y_value_std_plot = y_value_std[index_step, coord1_index_array, :]
                    mass_value_plot = mass_value[index_step, coord1_index_array, :]
                    
                    y_value_plot = np.nanmean(y_value_plot, axis=1)
                    y_value_std_plot = np.nanmean(y_value_std_plot, axis=1)
                    mass_value_plot = np.nanmean(mass_value_plot, axis=1)
                    y_value_plot = y_value_plot/mass_value_plot
                    y_value_std_plot = y_value_std_plot/mass_value_plot
                    
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=6,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=6,
                        )

    elif legend_class == None:
        if ave_over_coord1:
            if ave_over_coord2:
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot, axis=(1,2))
                    x_value_std_plot = np.nanmean(x_value_std_plot, axis=(1,2))
                y_value_plot = y_value[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot, axis=(1,2))
                y_value_std_plot = np.nanmean(y_value_std_plot, axis=(1,2))
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    if useerrorbar:
        fig.savefig(
        "_".join([
            diagram_path_add_nve + subfoldername + "errorbar/" + 'step',
            transfer_time_to_str(inputstepsarray),
            'coord1',
            transfer_coor_to_str(coord1_index_array, ave_over_coord1),
            'coord2',
            transfer_coor_to_str(coord2_index_array, ave_over_coord2),
            ".png",
        ]),
        format="png",
        )
    else:
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + 'step',
                transfer_time_to_str(inputstepsarray),
                'coord1',
                transfer_coor_to_str(coord1_index_array, ave_over_coord1),
                'coord2',
                transfer_coor_to_str(coord2_index_array, ave_over_coord2),
                ".png",
            ]),
            format="png",
        )
    # close figure after save
    plt.close('all')

def plot_1D_from_chunk2D_mask(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray, coord1_index_array, coord2_index_array,
        ave_over_coord1=False, ave_over_coord2=False,
        figure_class=None, legend_class=None, spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        x_scale = 'linear', y_scale = 'linear',
        useerrorbar = True,
        ifdivideNcount = False,
        ifmaskstatic=False,
        ifmasknonstatic=True,
    ):
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    # str for variable used in figure label legend
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
    # legend_class: time, coord1, coord2
    # spaceave: coord1, coord2
    
    time_array = time_from_start_rotate(inputstepsarray)
    if x_name == "timestep":
        x_value = inputstepsarray
        x_value_std = 0
    elif x_name == "time":
        x_value = time_array
        x_value_std = 0
    elif x_name == "strain":
        x_value = strain_from_rotate_start(inputstepsarray)
        x_value_std = 0
    else:
        x_value = get_value(lmp_path, x_name, n_ave, inputstepsarray)
        x_value_std = get_value(lmp_path, x_name + "_std", n_ave, inputstepsarray)

    if y_name == "timestep":
        y_value = inputstepsarray
        y_value_std = 0
    elif y_name == "time":
        y_value = time_array
        y_value_std = 0
    else:
        y_value = get_value(lmp_path, y_name, n_ave, inputstepsarray)
        y_value_std = get_value(lmp_path, y_name + "_std", n_ave, inputstepsarray)
    if ifdivideNcount:
        Ncount_value = get_value(lmp_path, 'Ncount', n_ave, inputstepsarray)
    # scale factor
    x_value = x_value/x_scale_factor
    x_value_std = x_value_std/x_scale_factor
    y_value = y_value/y_scale_factor
    y_value_std = y_value_std/y_scale_factor

    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    if ifdivideNcount:
        subfoldername = y_name + "_divideN_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    if useerrorbar:
        os.makedirs(diagram_path_add_nve + subfoldername + "errorbar/", exist_ok=True)

    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    # title

    if x_scale_str is None:
        x_label_str = x_name_label_in_figure
    else:
        x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
    ax.set_xlabel(x_label_str)
    
    if y_scale_str is None:
        y_label_str = y_name_label_in_figure
    else:
        y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
    ax.set_ylabel(y_label_str)
    if ifdivideNcount and y_label_str=="n_contact":
        ax.set_ylabel("<Z>")
    
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    char = c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord1 = load_coord(lmp_path, char, "Coord1")
    Coord2 = load_coord(lmp_path, char, "Coord2")
    x = Coord1/float(rr.logfile["dp"])
    y = Coord2/float(rr.logfile["dp"])
    maskstatic = (y-9) - (9-0.5)/(15.5-10.5)*(x-15.5) < 0
    masknonstatic = (y-25) - (25-0.5)/(15.5-5)*(x-5) > 0
    if ifmaskstatic:
        maskstaticornot = maskstatic
    if ifmasknonstatic:
        maskstaticornot = masknonstatic
    # legend_class: time, coord1, coord2
    if legend_class == 'tc1c2':
        for indexstep, step in enumerate(inputstepsarray):
            for coord1_index in coord1_index_array:
                for coord2_index in coord2_index_array:
                    strain = strain_from_rotate_start(step)
                    label = ", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                        "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                        "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[indexstep]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[indexstep, coord1_index, coord2_index]
                        x_value_std_plot = x_value_std[indexstep, coord1_index, coord2_index]
                    
                    y_value_plot = y_value[indexstep, coord1_index, coord2_index]
                    y_value_std_plot = y_value_std[indexstep, coord1_index, coord2_index]
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[indexstep, coord1_index, coord2_index]
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c1c2':
        for coord1_index in coord1_index_array:
            for coord2_index in coord2_index_array:
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                ])

                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index, coord2_index]
                    x_value_std_plot = x_value_std[:, coord1_index, coord2_index]
                
                y_value_plot = y_value[:, coord1_index, coord2_index]
                y_value_std_plot = y_value_std[:, coord1_index, coord2_index]
                
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index, coord2_index]
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    elif legend_class == 'c1':
        if coord2_index_array=='ave':
            if ave_over_coord1:
                maskstaticornot_plot = maskstaticornot[coord1_index_array, :]
                aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :]
                    x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_plot = y_value[:, coord1_index_array, :]
                y_value_std_plot = y_value_std[:, coord1_index_array, :]
                y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "({})".format(coord_scale_str),
                    "z=average",
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord1_index in coord1_index_array:
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y={:.1f} ({})".format(Coord1[coord1_index, 0]/coord_scale, coord_scale_str),
                        "z=average",
                    ])
                    maskstaticornot_plot = maskstaticornot[coord1_index, :]
                    aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, coord1_index, :]
                        x_value_std_plot = x_value_std[:, coord1_index, :]
                        #breakpoint()
                        x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_plot = y_value[:, coord1_index, :]
                    y_value_std_plot = y_value_std[:, coord1_index, :]
                    y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, coord1_index, :]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c2':
        if coord1_index_array=='ave':
            if ave_over_coord2:
                maskstaticornot_plot = maskstaticornot[:, coord2_index_array]
                aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_plot = y_value[:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "y=average",
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord2_index in coord2_index_array:
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y=average",
                        "z={:.1f} ({})".format(Coord2[0, coord2_index]/coord_scale, coord_scale_str),
                    ])
                    maskstaticornot_plot = maskstaticornot[:, coord2_index]
                    aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, :, coord2_index]
                        x_value_std_plot = x_value_std[:, :, coord2_index]
                        x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_plot = y_value[:, :, coord2_index]
                    y_value_std_plot = y_value_std[:, :, coord2_index]
                    y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, :, coord2_index]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == None:
        if ave_over_coord1:
            if ave_over_coord2:
                maskstaticornot_plot = maskstaticornot[coord1_index_array, :][:, coord2_index_array]
                aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_plot = y_value[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    if useerrorbar:
        fig.savefig(
        "_".join([
            diagram_path_add_nve + subfoldername + "errorbar/" + "mask_" + 'step',
            transfer_time_to_str(inputstepsarray),
            'coord1',
            transfer_coor_to_str(coord1_index_array, ave_over_coord1),
            'coord2',
            transfer_coor_to_str(coord2_index_array, ave_over_coord2),
            ".png",
        ]),
        format="png",
        )
    else:
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + "mask_" + 'step',
                transfer_time_to_str(inputstepsarray),
                'coord1',
                transfer_coor_to_str(coord1_index_array, ave_over_coord1),
                'coord2',
                transfer_coor_to_str(coord2_index_array, ave_over_coord2),
                ".png",
            ]),
            format="png",
        )
    # close figure after save
    plt.close('all')



def plot_1D_for_chunk1D_near_wall(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray, coord_index_array,
        figure_class=None, legend_class=None, spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        x_scale = 'linear', y_scale = 'linear',
        useerrorbar = True,
        ifsum = False,
    ):
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    # str for variable used in figure label legend
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
    # legend_class: time, coord
    # spaceave: coord

    time_array = time_from_start_rotate(inputstepsarray)
    if x_name == "timestep":
        x_value = inputstepsarray
        x_value_std = 0
    elif x_name == "time":
        x_value = time_array
        x_value_std = 0
    elif x_name == "strain":
        x_value = strain_from_rotate_start(inputstepsarray)
        x_value_std = 0
    else:
        x_value = get_value(lmp_path, x_name, n_ave, inputstepsarray)
        x_value_std = get_value(lmp_path, x_name + "_std", n_ave, inputstepsarray)

    if y_name == "timestep":
        y_value = inputstepsarray
        y_value_std = 0
    elif y_name == "time":
        y_value = time_array
        y_value_std = 0
    else:
        y_value = get_value(lmp_path, y_name, n_ave, inputstepsarray)
        y_value_std = get_value(lmp_path, y_name + "_std", n_ave, inputstepsarray)

    # scale factor
    x_value = x_value/x_scale_factor
    x_value_std = x_value_std/x_scale_factor
    y_value = y_value/y_scale_factor
    y_value_std = y_value_std/y_scale_factor

    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    if useerrorbar:
        os.makedirs(diagram_path_add_nve + subfoldername + "errorbar/", exist_ok=True)

    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    # title
    
    
    if x_scale_str is None:
        x_label_str = x_name_label_in_figure
    else:
        x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
    ax.set_xlabel(x_label_str)
    if y_scale_str is None:
        y_label_str = y_name_label_in_figure
    else:
        y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
    ax.set_ylabel(y_label_str)
    
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    char = c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    # check if only 1 coord
    if len(c_r_coord_npyfilepath[char].keys()) == 1:
        pass
    else:
        print("the coord near wall is not 1")
        breakpoint()
        sys.exit()
    Coord_str = list(c_r_coord_npyfilepath[char].keys())[0]
    Coord = load_coord(lmp_path, char, Coord_str)
    if Coord_str == "Coord1":
        coord_str_in_plot = "y"
    elif Coord_str == "Coord2":
        coord_str_in_plot = "z"
    else:
        sys.exit("not y not z")
    # legend_class: time, coord1, coord2
    if ifsum==False:
        if legend_class == 'tc':
            for indexstep, step in enumerate(inputstepsarray):
                for coord_index in coord_index_array:
                    time = time_from_start_rotate(step)
                    strain = strain_from_rotate_start(time)
                    label = ", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                        coord_str_in_plot + "={:.1f} ({})".format(Coord[coord_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[indexstep]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[indexstep, coord_index]
                        x_value_std_plot = x_value_std[indexstep, coord_index]
                    
                    y_value_plot = y_value[indexstep, coord_index]
                    y_value_std_plot = y_value_std[indexstep, coord_index]
                    
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
        elif legend_class == 'c':
            for coord_index in coord_index_array:                    
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    coord_str_in_plot + "={:.1f} ({})".format(Coord[coord_index]/coord_scale, coord_scale_str),
                ])
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord_index]
                    x_value_std_plot = x_value_std[:, coord_index]
                
                y_value_plot = y_value[:, coord_index]
                y_value_std_plot = y_value_std[:, coord_index]
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    else:
        strain1 = strain_from_rotate_start(inputstepsarray[0])
        strain2 = strain_from_rotate_start(inputstepsarray[-1])
        label=", ".join([
            r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
        ])
        if x_name == "timestep" or x_name == "time" or x_name == "strain":
            x_value_plot = x_value[:]
            x_value_std_plot = 0
        else:
            x_value_plot = x_value[:, :]
            x_value_std_plot = x_value_std[:, :]
            x_value_plot = np.average(x_value_plot, axis=1)
        
        y_value_plot = y_value[:, :]
        y_value_std_plot = y_value_std[:, :]
        
        y_value_plot = np.average(y_value_plot, axis=1)
        ax.plot(
            x_value_plot, y_value_plot,
            label=label,
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    if ifsum==False:
        if useerrorbar:
            fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + "errorbar/step",
                transfer_time_to_str(inputstepsarray),
                Coord_str,
                transfer_coor_to_str(coord_index_array),
                ".png",
            ]),
            format="png",
            )
        else:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + 'step',
                    transfer_time_to_str(inputstepsarray),
                    Coord_str,
                    transfer_coor_to_str(coord_index_array),
                    ".png",
                ]),
                format="png",
                )
    else:
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + 'step',
                transfer_time_to_str(inputstepsarray),
                ".png",
            ]),
            format="png",
        )
    # close figure after save
    plt.close('all')

def plot_quiver_from_chunk2D(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray,
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=0.1, label_scale=0.2,
        ifloglength=False,
        logvaluemin=10**-4,
        valueminus=0,
        ifplotseparateupdown=False,
        quiver_scale_up=0.1, label_scale_up=0.2,
        quiver_scale_down=0.1, label_scale_down=0.2,
        ifstreamplot=True,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=-10, vmax=10,
    ):
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    if ifloglength:
        quiver_scale = np.log10(quiver_scale/logvaluemin)
        label_scale = np.log10(label_scale/logvaluemin)
        quiver_scale_up = np.log10(quiver_scale_up/logvaluemin)
        label_scale_up = np.log10(label_scale_up/logvaluemin)
        quiver_scale_down = np.log10(quiver_scale_down/logvaluemin)
        label_scale_down = np.log10(label_scale_down/logvaluemin)
    # str for variable used in figure label legend
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
    # spaceave: coord1, coord2
    time_array = time_from_start_rotate(inputstepsarray)
    
    x_value = get_value(lmp_path, x_name, n_ave, inputstepsarray)
    y_value = get_value(lmp_path, y_name, n_ave, inputstepsarray)
    mask = x_value>0.1
    x_value[mask] = x_value[mask] - valueminus
    mask = y_value>0.1
    y_value[mask] = y_value[mask] - valueminus
    # scale factor
    x_value = x_value/x_scale_factor
    y_value = y_value/y_scale_factor

    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    os.makedirs(diagram_path_add_nve + subfoldername + "streamplot/", exist_ok=True)
    if ifplotseparateupdown:
        os.makedirs(diagram_path_add_nve + subfoldername + "up/", exist_ok=True)
        os.makedirs(diagram_path_add_nve + subfoldername + "down/", exist_ok=True)
    if ifloglength:
        os.makedirs(diagram_path_add_nve + subfoldername + "log/", exist_ok=True)
        if ifplotseparateupdown:
            os.makedirs(diagram_path_add_nve + subfoldername + "log/" + "up/", exist_ok=True)
            os.makedirs(diagram_path_add_nve + subfoldername + "log/" + "down/", exist_ok=True)
    char = c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord1 = load_coord(lmp_path, char, "Coord1")
    Coord2 = load_coord(lmp_path, char, "Coord2")
    # plot quiver vector field
    for indexstep, step in enumerate(inputstepsarray):
        x_value_plot = x_value[indexstep]
        y_value_plot = y_value[indexstep]
        # plot ave_z velocity across y
        fig, ax = plt.subplots()
        
        # title

        if x_scale_str is None:
            x_label_str = x_name_label_in_figure
        else:
            x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
        ax.set_xlabel(x_label_str)
        
        if y_scale_str is None:
            y_label_str = y_name_label_in_figure
        else:
            y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
        ax.set_ylabel(y_label_str)            
        # plot
        if ifloglength:
            length = (x_value_plot**2 + y_value_plot**2)**0.5
            masktolog = (length>=logvaluemin)
            length_log = np.copy(length)
            length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
            #breakpoint()
            x_value_plot_log = length_log*x_value_plot/length
            y_value_plot_log = length_log*y_value_plot/length
            Q = ax.quiver(
                Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                cmap='hsv',
            )
        else:
            Q = ax.quiver(
                Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                cmap='hsv',
            )
            
        
        ax.set_xlabel(
            'y ({})'.format(coord_scale_str)
        )
        ax.set_ylabel(
            'z ({})'.format(coord_scale_str)
        )
        time = time_from_start_rotate(step)
        if ifloglength:
            ax.quiverkey(
                Q, 0.1, 0.95, label_scale,
                label = "logscale, : {:.2e} in 45 degree, ".format(label_scale) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time),
                labelpos='E',
                coordinates='figure', angle=45,
            )
        else:
            ax.quiverkey(
                Q, 0.1, 0.95, label_scale,
                label = " : {:.2e} in 45 degree, ".format(label_scale) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time),
                labelpos='E',
                coordinates='figure', angle=45,
            )
            
        plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
        plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

        if ifloglength:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + "log/" + 'step',
                    str(step),
                    ".png",
                ]),
                format="png",
                bbox_inches=None,
            )
        else:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + 'step',
                    str(step),
                    ".png",
                ]),
                format="png",
                bbox_inches=None,
            )

        # close figure after save
        plt.close('all')
    # plot up and down separate figure
    if ifplotseparateupdown:
        maskup = y_value>0
        maskdown = y_value<0
        x_value_up = np.zeros_like(x_value)
        y_value_up = np.zeros_like(y_value)
        x_value_down = np.zeros_like(x_value)
        y_value_down = np.zeros_like(y_value)
        x_value_up[maskup] = x_value[maskup]
        y_value_up[maskup] = y_value[maskup]
        x_value_down[maskdown] = x_value[maskdown]
        y_value_down[maskdown] = y_value[maskdown]
        for indexstep, step in enumerate(inputstepsarray):
            x_value_plot = x_value_up[indexstep]
            y_value_plot = y_value_up[indexstep]
            total_up = np.sum(y_value_plot)
            # plot ave_z velocity across y
            fig, ax = plt.subplots()
            
            # title

            if x_scale_str is None:
                x_label_str = x_name_label_in_figure
            else:
                x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
            ax.set_xlabel(x_label_str)
            
            if y_scale_str is None:
                y_label_str = y_name_label_in_figure
            else:
                y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
            ax.set_ylabel(y_label_str)            
            # plot
            if ifloglength:
                length = (x_value_plot**2 + y_value_plot**2)**0.5
                masktolog = (length>=logvaluemin)
                length_log = np.copy(length)
                length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
                #breakpoint()
                x_value_plot_log = length_log*x_value_plot/length
                y_value_plot_log = length_log*y_value_plot/length
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_up,
                    cmap='hsv',
                )
            else:
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_up,
                    cmap='hsv',
                )
                
            
            ax.set_xlabel(
                'y ({})'.format(coord_scale_str)
            )
            ax.set_ylabel(
                'z ({})'.format(coord_scale_str)
            )
            time = time_from_start_rotate(step)
            if ifloglength:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_up,
                    label = "up, logscale, : {:.2e} in 45 degree, ".format(label_scale_up) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total up is {:.2e}".format(total_up),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
            else:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_up,
                    label = " : up, {:.2e} in 45 degree, ".format(label_scale_up) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total up is {:.2e}".format(total_up),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
                
            plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
            plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

            if ifloglength:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "log/" + "up/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )
            else:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "up/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )

            # close figure after save
            plt.close('all')

        for indexstep, step in enumerate(inputstepsarray):
            x_value_plot = x_value_down[indexstep]
            y_value_plot = y_value_down[indexstep]
            total_down = np.sum(y_value_plot)
            # plot ave_z velocity across y
            fig, ax = plt.subplots()
            
            # title

            if x_scale_str is None:
                x_label_str = x_name_label_in_figure
            else:
                x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
            ax.set_xlabel(x_label_str)
            
            if y_scale_str is None:
                y_label_str = y_name_label_in_figure
            else:
                y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
            ax.set_ylabel(y_label_str)            
            # plot
            if ifloglength:
                length = (x_value_plot**2 + y_value_plot**2)**0.5
                masktolog = (length>=logvaluemin)
                length_log = np.copy(length)
                length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
                #breakpoint()
                x_value_plot_log = length_log*x_value_plot/length
                y_value_plot_log = length_log*y_value_plot/length
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_down,
                    cmap='hsv',
                )
            else:
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_down,
                    cmap='hsv',
                )
                
            
            ax.set_xlabel(
                'y ({})'.format(coord_scale_str)
            )
            ax.set_ylabel(
                'z ({})'.format(coord_scale_str)
            )
            time = time_from_start_rotate(step)
            if ifloglength:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_down,
                    label = "down, logscale, : {:.2e} in 45 degree, ".format(label_scale_down) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total down is {:.2e}".format(total_down),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
            else:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_down,
                    label = " : down, {:.2e} in 45 degree, ".format(label_scale_down) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total down is {:.2e}".format(total_down),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
                
            plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
            plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

            if ifloglength:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "log/" + "down/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )
            else:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "down/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )

            # close figure after save
            plt.close('all')
    # plot streamplot
    for indexstep, step in enumerate(inputstepsarray):
        x_value_plot = x_value[indexstep]
        y_value_plot = y_value[indexstep]
        # plot ave_z velocity across y
        fig, ax = plt.subplots()
        
        # title

        if x_scale_str is None:
            x_label_str = x_name_label_in_figure
        else:
            x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
        ax.set_xlabel(x_label_str)
        
        if y_scale_str is None:
            y_label_str = y_name_label_in_figure
        else:
            y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
        ax.set_ylabel(y_label_str)
        # check if line0
        if ifstreamplot:
            # plot
            strm = ax.streamplot(
                Coord1[:,0]/coord_scale, Coord2[0,:]/coord_scale, np.transpose(x_value_plot), np.transpose(y_value_plot),
                linewidth=1, color='k',
                density=[0.8, 0.8],
            )
        vector_length = np.sqrt(x_value_plot**2 + y_value_plot**2)
        vector_length = np.ma.masked_where(np.logical_not(vector_length > 0), vector_length)
        #breakpoint()
        d_Coord1 = Coord1[:,0][1]-Coord1[:,0][0]
        d_Coord2 = Coord2[0,:][1]-Coord2[0,:][0]
        Coord1_expand = Coord1[:,0] - d_Coord1/2
        Coord1_expand = np.append(Coord1_expand, (Coord1_expand[-1]+d_Coord1))
        Coord2_expand = Coord2[0,:] - d_Coord2/2
        Coord2_expand = np.append(Coord2_expand, (Coord2_expand[-1]+d_Coord2))
        if ifcontour:
            if contour_v_min_max == "constant":
                if contour_norm == "linear":
                    contour = ax.pcolor(
                        Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                        norm=colors.Normalize(vmin=vmin, vmax=vmax), #norm=colors.LogNorm(vmin=10**-1, vmax=10), #norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                        cmap='coolwarm',
                    )
                elif contour_norm == "log":
                    contour = ax.pcolor(
                    Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                    norm=colors.LogNorm(vmin=vmin, vmax=vmax), #norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                    cmap='coolwarm',
                )
                else:
                    sys.exit("Error: contour_norm not defined")
            elif contour_v_min_max == "min_to_max":
                if contour_norm == "linear":
                    contour = ax.pcolor(
                        Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                        norm=colors.Normalize(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                        cmap='coolwarm',
                    )
                elif contour_norm == "log":
                    contour = ax.pcolor(
                    Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                    norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                    cmap='coolwarm',
                )
                else:
                    sys.exit("Error: contour_norm not defined")
            else:
                sys.exit("Error: contour_v_min_max not defined")
            fig.colorbar(contour, ax=ax, extend='max')

        time = time_from_start_rotate(step)
        strain = strain_from_rotate_start(step)
        ax.set_title(
            "t = {:.2f} s".format(time) + ", strain is {:.2f}".format(strain)
        )
        ax.set_xlabel(
            'y ({})'.format(coord_scale_str)
        )
        ax.set_ylabel(
            'z ({})'.format(coord_scale_str)
        )
        
            
        plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
        plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + 'streamplot/' + 'step',
                str(step),
                ".png",
            ]),
            format="png",
            bbox_inches=None,
        )

        # close figure after save
        plt.close('all')

    
# main exclusive
def main():
    #data_auto(lmp_path = rr.folder_path_list_initial_to_last[-1])
    #data_auto(lmp_path = rr.folder_path_list_initial_to_last[-2])
    #data_auto(lmp_path = rr.folder_path_list_initial_to_last[-3])
    # set plot path 
    lmp_path = rr.folder_path_list_initial_to_last[-1]
    # contour for pressure and stress12 mu_12
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "strain_rate_21_middle", "strain_rate_21_middle", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*strain_rate_scale, x_scale_str=None, y_scale_factor=2**0.5*strain_rate_scale, y_scale_str=None,
        quiver_scale=strain_rate_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=5,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "mu_12", "mu_12", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=0.8,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "pressure", "pressure", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale, x_scale_str=None, y_scale_factor=2**0.5*stress_scale, y_scale_str=None,
        quiver_scale=stress_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=10,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "stress_12", "stress_12", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale, x_scale_str=None, y_scale_factor=2**0.5*stress_scale, y_scale_str=None,
        quiver_scale=stress_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=3,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "fraction", "fraction", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0.3, vmax=0.7,
    )
    

    # contour for pressure and stress12 mu_12
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "strain_rate_21_middle", "strain_rate_21_middle", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*strain_rate_scale, x_scale_str=None, y_scale_factor=2**0.5*strain_rate_scale, y_scale_str=None,
        quiver_scale=strain_rate_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=5,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "mu_12", "mu_12", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=0.8,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "pressure", "pressure", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale, x_scale_str=None, y_scale_factor=2**0.5*stress_scale, y_scale_str=None,
        quiver_scale=stress_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=10,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "stress_12", "stress_12", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale, x_scale_str=None, y_scale_factor=2**0.5*stress_scale, y_scale_str=None,
        quiver_scale=stress_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=3,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "fraction", "fraction", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0.3, vmax=0.7,
    )
    
    # 1Dplot from chunk2D
    # velocity
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.append(np.arange(4500000, 7000000, 500000), np.array([40000000])), "Coord1", "mv_1", coord_scale, coord_scale_str, velocity_scale, velocity_scale_str, 'log'),
        (np.append(np.arange(4500000, 7000000, 500000), np.array([40000000])), "Coord1", "mv_1", coord_scale, coord_scale_str, velocity_scale, velocity_scale_str, 'linear'),
    ]:
        coord1_index_array = np.arange(16)
        plot_1D_from_chunk2D(
            lmp_path,
            51, x_name, y_name, input_stepsarray, coord1_index_array, 'sumdividemass', legend_class = 't',
            x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
            x_scale=x_scale,
            useerrorbar=False,
        )
    
    # wall stress
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.arange(5000000,40000000,1000000),"strain", "inwall_stress_1", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "inwall_stress_2", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "inwall_stress_3", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "outwall_stress_1", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "outwall_stress_2", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "outwall_stress_3", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "zbottom_stress_1", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "zbottom_stress_2", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "zbottom_stress_3", 1, None, stress_scale, stress_scale_str, 'linear'),
    ]:
        plot_1D_for_chunk1D_near_wall(
            lmp_path,
            51, x_name, y_name, input_stepsarray, coord_index_array=None, legend_class = 'c',
            x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
            x_scale=x_scale,
            useerrorbar=False,
            ifsum=True
        )
    # contact number
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale, ifdivideNcount) in [
        (np.arange(1000000,40000000,1000000),"strain", "n_contact", 1, None, 1, None, 'linear', True),
        (np.arange(1000000,40000000,1000000),"strain", "fraction", 1, None, 1, None, 'linear', False),
        ]:
        for (coord1_index_array, coord2_index_array, legend_class, ave_over_coord1, ave_over_coord2) in [
            ([0, 1, 2], 'ave', 'c1', True, False),
            ([-1, -2, -3], 'ave', 'c1', True, False),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                )
                plot_1D_from_chunk2D_mask(
                    lmp_path,
                    51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                    ifmaskstatic=False,
                    ifmasknonstatic=True,
                )
    # mu I
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, y_scale) in [
        (np.arange(6000000,40000000,1000000), "mu_12_middle", "I_12", mu_scale, mu_scale_str, I_scale, None, 'log'),
        (np.array([30000000]), "mu_12_middle", "I_12", mu_scale, mu_scale_str, I_scale, None, 'log'),
        # strain
        (np.array([30000000]), "strain", "pressure", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.array([30000000]), "strain", "stress_12", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.array([30000000]), "strain", "strain_rate_21_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    ]:
        for (coord1_index_array, coord2_index_array) in [
                ([0], np.arange(0,16,2)),
                ([-1], np.arange(0,16,2)),
                (np.arange(0,15,2), [0]),
                ]:
            plot_1D_from_chunk2D_mask(
                lmp_path,
                51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                x_scale=x_scale,
                y_scale=y_scale,
                useerrorbar=useerrorbar,
                ifmaskstatic=False,
                ifmasknonstatic=True,
            )
            plot_1D_from_chunk2D(
                lmp_path,
                51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                x_scale=x_scale,
                y_scale=y_scale,
                useerrorbar=useerrorbar,
            )
    
    
    """
    # 1Dplot from chunk2D

    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.arange(1000000,40000000,1000000),"strain", "fraction", 1, None, 1, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "n_contact", 1, None, 1, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"I_12", "mu_12_middle", I_scale, None, mu_scale, mu_scale_str, 'log'),
        (np.arange(1000000,40000000,1000000),"I_tensor", "mu_tensor_12", I_tensor_scale, I_tensor_scale_str, mu_tensor_scale, mu_tensor_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "I_12", 1, None, I_scale, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "mu_12_middle", 1, None, mu_scale, mu_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "I_tensor", 1, None, I_scale, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "mu_tensor_12", 1, None, mu_scale, mu_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "velocity_1", 1, None, velocity_scale, velocity_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "velocity_2", 1, None, velocity_scale, velocity_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "velocity_3", 1, None, velocity_scale, velocity_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_21_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_31_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_22_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_32_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_23_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_33_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_11", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_22", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_33", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_12", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_13", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_23", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "mu_12_middle", 1, None, mu_scale, mu_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "I_12", 1, None, I_scale, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "mu_tensor_12", 1, None, mu_tensor_scale, mu_tensor_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "I_tensor", 1, None, I_tensor_scale, I_tensor_scale_str, 'linear'),
    ]:
        
        for (coord1_index_array, coord2_index_array) in [
            ([0], [0,1,5,7,8]),
            ([-1], [0,1,5,7,8]),
            ([0,1,2,3,4,5,6,-2,-1], [0]),
            ([0,1,2,3,4,5,6,-2,-1], [1]),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                )
                plot_1D_from_chunk2D_mask(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ifmaskstatic=False,
                    ifmasknonstatic=True,
                )
        for (coord1_index_array, coord2_index_array, legend_class) in [
            ([0,-1], 'ave', 'c1'),
            ('ave', [0], 'c2'),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                )
                #plot_1D_from_chunk2D_mask(
                #    lmp_path,
                #    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                #    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                #    x_scale=x_scale,
                #    useerrorbar=useerrorbar,
                #    ifmaskstatic=False,
                #    ifmasknonstatic=True,
                #)
    
    
    
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale, ifdivideNcount) in [
        (np.arange(1000000,40000000,1000000),"strain", "n_contact", 1, None, 1, None, 'linear', True),
        (np.arange(1000000,40000000,1000000),"strain", "fraction", 1, None, 1, None, 'linear', False),
        ]:
        for (coord1_index_array, coord2_index_array, legend_class, ave_over_coord1, ave_over_coord2) in [
            ([0, 1], 'ave', 'c1', True, False),
            ([-1, -2], 'ave', 'c1', True, False),
            ([0, 1, 2], 'ave', 'c1', True, False),
            ([1, 2, 3], 'ave', 'c1', True, False),
            ([2, 3, 4], 'ave', 'c1', True, False),
            ([3, 4, 5], 'ave', 'c1', True, False),
            ([4, 5, 6], 'ave', 'c1', True, False),
            ([5, 6, 7], 'ave', 'c1', True, False),
            ([-1, -2, -3], 'ave', 'c1', True, False),
            ('ave', [0], 'c2', False, False),
            ([0, 1], [0,1,2,3,4,5,6,7], None, True, True),
            ([-1, -2], [0,1,2,3,4,5,6,7], None, True, True),
            ([0, 1, 2], [0,1,2,3,4,5,6,7], None, True, True),
            ([1, 2, 3], [0,1,2,3,4,5,6,7], None, True, True),
            ([2, 3, 4], [0,1,2,3,4,5,6,7], None, True, True),
            ([3, 4, 5], [0,1,2,3,4,5,6,7], None, True, True),
            ([4, 5, 6], [0,1,2,3,4,5,6,7], None, True, True),
            ([5, 6, 7], [0,1,2,3,4,5,6,7], None, True, True),
            ([6, 7, 8], [0,1,2,3,4,5,6,7], None, True, True),
            ([7, 8, 9], [0,1,2,3,4,5,6,7], None, True, True),
            ([8, 9, 10], [0,1,2,3,4,5,6,7], None, True, True),
            ([9, 10, 11], [0,1,2,3,4,5,6,7], None, True, True),
            ([10, 11, 12], [0,1,2,3,4,5,6,7], None, True, True),
            ([11, 12, 13], [0,1,2,3,4,5,6,7], None, True, True),
            ([12, 13, 14], [0,1,2,3,4,5,6,7], None, True, True),
            ([-1, -2, -3], [0,1,2,3,4,5,6,7], None, True, True),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                )
                plot_1D_from_chunk2D_mask(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                    ifmaskstatic=False,
                    ifmasknonstatic=True,
                )
    
    
    # 1Dplot from chunk1D
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.arange(1000000,40000000,1000000),"strain", "inwall_stress_1", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "inwall_stress_2", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "inwall_stress_3", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "outwall_stress_1", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "outwall_stress_2", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "outwall_stress_3", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "zbottom_stress_1", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "zbottom_stress_2", 1, None, stress_scale, stress_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "zbottom_stress_3", 1, None, stress_scale, stress_scale_str, 'linear'),
        ]:
        for (coord_index_array) in [
            [0,1,2,3,4,5,-2,-1],
            ]:
            for useerrorbar in [True, False]:
                plot_1D_for_chunk1D_near_wall(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord_index_array, legend_class = 'c',
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                )
    
    

    # 2D plot (vector field, streamlines, etc)
    #plot_quiver
    for ifloglength in [True, False]:
        # mv
        
        plot_quiver_from_chunk2D(
            lmp_path,
            n_ave, "mv_2", "mv_3", np.arange(1000000,40000000,1000000),
            spaceave=None,
            x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
            quiver_scale=5*10**-9, label_scale=5*10**-9,
            ifloglength=ifloglength,
            ifplotseparateupdown=True,
            quiver_scale_up=5*10**-9, label_scale_up=5*10**-9,
            quiver_scale_down=5*10**-9, label_scale_down=5*10**-9,
        )
        plot_quiver_from_chunk2D(
            lmp_path,
            n_ave, "velocity_2", "velocity_3", np.arange(1000000,40000000,1000000),
            spaceave=None,
            x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
            quiver_scale=0.0001, label_scale=0.0001,
            ifloglength=ifloglength,
            ifplotseparateupdown=True,
            quiver_scale_up=0.0001, label_scale_up=0.0001,
            quiver_scale_down=0.0005, label_scale_down=0.0005,
        )

        plot_quiver_from_chunk2D(
            lmp_path,
            n_ave, "velocity_1", "velocity_1", np.arange(1000000,40000000,1000000),
            spaceave=None,
            x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
            quiver_scale=0.002, label_scale=0.002,
            ifloglength=ifloglength,
        )
    # fraction
    plot_quiver_from_chunk2D(
        lmp_path,
        n_ave, "fraction", "fraction", np.arange(1000000,40000000,1000000),
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=0.2, label_scale=1,
        ifloglength=False,
        valueminus=0.5,
    )
    

    # plot shear rate for compare with quiver
    """
    
    print('\a')
# main exclusive
if __name__ == "__main__":
    
    main()