import os
import sys
lmp_folder_path_dic = {
    1: os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '01_raw',
        'Link to lmp_run',
        'block_xp_50_w_16_h_30_Sa_2e-6_his_yes_xmu_5e-1_noshearfor_5e6_and_run_1e7',
        'f_5e6_run_1e7',
        'f_15e6_run_15e7',
    ),
    2: os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '01_raw',
        'Link to lmp_run',
        '20200921_nott_H_60_W_16_L_50',
        'f_5e6',
        'f_15e6',
        'f_45e6',
    ),
    3: os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '01_raw',
        'Link to lmp_run',
        '20200921_nott_H_90_W_16_L_50',
        'f_5e6',
        'f_15e6',
        'f_35e6',
    ),
    4: os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '01_raw',
        'Link to lmp_run',
        'block_xp_50_w_16_h_120_Sa_2e-6_his_yes_xmu_5e-1_noshear_2e7',
        'f_2e7_run_5e7',
    ),
    5: os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '01_raw',
        'Link to lmp_run',
        'block_xp_50_w_16_h_30_Sa_2e-6_his_yes_xmu_5e-1_noshearfor_5e6_and_run_1e7',
        'f_5e6_run_1e7_v1_20191111_Sa_2e-7',
    ),
    6: os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '01_raw',
        'Link to lmp_run',
        'block_xp_50_w_16_h_30_Sa_2e-6_his_yes_xmu_5e-1_noshearfor_5e6_and_run_1e7',
        'f_5e6_run_1e7_v1_20191111_Sa_2e-5',
    ),
}
lmp_folder_path = lmp_folder_path_dic[2]

log_output_folder_name = "log"
logpicklefilename = 'log.pickle'



def folder_name_under_subfolder_of_data(log_variable_dic_list):
    # L W H Sa rotate_start_time endtime history wallshape xmu hooke/hertz kn kt gamma_n gamma_t timestep
    variable_name = [

    ]
    variable_label_string = {

    }
    log_variable = log_variable_dic_list[-1]
    name = (
        "W_"
        + str(log_variable["width_wall_dp_unit"])
        + "_H_" + str(log_variable["z_length_create_dp_unit"])
        + "_L_" + str(log_variable["x_period_dp_unit"])
        + "_Sa_" + "{:.2e}".format(float(log_variable["Sa"]))
        + "_total_time_" + "{:.2e}".format(float(log_variable["total_time"]))
    )
    return name

def intermediate_subfolder_path(folder_name_under_subfolder_of_data):
    folder_path = os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '02_intermediate',
        folder_name_under_subfolder_of_data,
    )
    return folder_path

coord_name_list = [
    "Coord1",
    "Coord2",
]

coord_2_3_fixtimeave_id_name = 'coord1and2_chunk_2_3'

coord_chunk_id_23 = 'coord1and2_chunk_2_3'
coord_chunk_id_inwall = 'coord1and2_chunk_near_inwall'
coord_chunk_id_outwall = 'coord1and2_chunk_near_outwall'
coord_chunk_id_zbottom = 'coord1and2_chunk_near_zbottom'

coord_chunk_id_23_replace = 'avspatial'

coord_chunk_id_list = [
    coord_chunk_id_23,
    coord_chunk_id_inwall,
    coord_chunk_id_outwall,
    coord_chunk_id_zbottom,
]

map_fixtimeave_value_to_coord_by_id = {
    'avspatial_ave': coord_chunk_id_23,
    'avspatial_omega_ave': coord_chunk_id_23,
    'avspatialstress_ave': coord_chunk_id_23,
    'contact_ave': coord_chunk_id_23,
    'ave_std_inwall': coord_chunk_id_inwall,
    'ave_std_outwall': coord_chunk_id_outwall,
    'ave_std_zbottom': coord_chunk_id_zbottom,
}

map_fixtimeave_to_fixchunkave = {
    'avspatial_ave': 'avspatial',
}

no_coord_fixtimeave = [
    'timeav_inwall_force',
    'timeav_outwall_force',
    'timeav_zbottom_force',
    'timeav_y_bottom_force',
    'timeav_y_top_force',
]

chunk_output_list = list(map_fixtimeave_value_to_coord_by_id.keys()) + coord_chunk_id_list
outputlist = list(map_fixtimeave_value_to_coord_by_id.keys()) + coord_chunk_id_list + no_coord_fixtimeave

npy_output_subfolder_name_map_from_id = {
    coord_chunk_id_23: "coord_chunk_2_3",
    coord_chunk_id_inwall: "coord_chunk_inwall",
    coord_chunk_id_outwall: "coord_chunk_outwall",
    coord_chunk_id_zbottom: "coord_chunk_zbottom",
    'avspatial_ave': "chunk_spatial",
    'avspatial_omega_ave': "chunk_spatial_omega",
    'avspatialstress_ave': "chunk_spatial_stress",
    'contact_ave': "chunk_spatial_contact",
    'ave_std_inwall': "chunk_inwall",
    'ave_std_outwall': "chunk_outwall",
    'ave_std_zbottom': "chunk_zbottom",
    'timeav_inwall_force': "inwall_force",
    'timeav_outwall_force': "outwall_force",
    'timeav_zbottom_force': "zbottom_force",
    'timeav_y_bottom_force': "inwall_force",
    'timeav_y_top_force': "outwall_force",
}
for id_fix in outputlist:
    if id_fix not in npy_output_subfolder_name_map_from_id:
        sys.exit("some id not includeed in npy_output_subfolder_name_map_from_id")

output_shape_map_from_id = {
    coord_chunk_id_23:  ['n_1', 'n_2'],
    'coord1and2_chunk_near_inwall': ['n_2'],
    'coord1and2_chunk_near_outwall': ['n_2'],
    'coord1and2_chunk_near_zbottom': ['n_1'],
    'avspatial_ave': ['t', 'n_1', 'n_2'],
    'avspatial_omega_ave': ['t', 'n_1', 'n_2'],
    'avspatialstress_ave': ['t', 'n_1', 'n_2'],
    'contact_ave': ['t', 'n_1', 'n_2'],
    'ave_std_inwall': ['t', 'n_2'],
    'ave_std_outwall': ['t', 'n_2'],
    'ave_std_zbottom': ['t', 'n_1'],
}

for id_fix in chunk_output_list:
    if id_fix not in output_shape_map_from_id:
        sys.exit("some id not includeed in output_shape_map_from_id")

# npy
def npy_raw_folder_path(n, log_variable_dic_list):
    folder_path = os.path.join(
        intermediate_subfolder_path(folder_name_under_subfolder_of_data(log_variable_dic_list)),
        'npy',
        'simu_' + str(n),
    )
    return folder_path

def fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list):
    folder_path = os.path.join(
        npy_raw_folder_path(n, log_variable_dic_list),
        npy_output_subfolder_name_map_from_id[fixtimeave_id_name],
    )
    return folder_path

def fixtimeave_text_file_path(n, fixtimeave_id_name, log_variable_dic_list, folder_path_list_initial_to_last):

    folder_path = os.path.join(
        folder_path_list_initial_to_last[n],
        log_variable_dic_list[n]['fixavetime'][fixtimeave_id_name]['file'],
    )
    return folder_path

def fixchunkave_text_file_path(n, fixchunkave_id_name, log_variable_dic_list, folder_path_list_initial_to_last):

    folder_path = os.path.join(
        folder_path_list_initial_to_last[n],
        log_variable_dic_list[n]['fixavechunk'][fixchunkave_id_name]['file'],
    )
    return folder_path

def add_npy_if_not(filename):
    if filename[-4:] == ".npy":
        pass
    else:
        filename = filename + ".npy"
    return filename

def eliminate_npy_if_yes(filename):
    if filename[-4:] == ".npy":
        filename = filename[:-4]
    else:
        pass
    return filename

def fixtimeave_npy_output_file_path(n, fixtimeave_id_name, log_variable_dic_list, filename):
    filepath = os.path.join(
        fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list),
        add_npy_if_not(filename),
    )
    return filepath

def fixtimeave_npy_output_coord_file_path(n, fixtimeave_id_name_for_v, log_variable_dic_list, coord_name):
    fixtimeave_id_name_for_coord = map_fixtimeave_value_to_coord_by_id[fixtimeave_id_name_for_v]
    filepath = fixtimeave_npy_output_file_path(n, fixtimeave_id_name_for_coord, log_variable_dic_list, coord_name)
    return filepath

dic_rename_fixtimeave_npy_headername_to_use = {
    'Coord1': 'Coord1',
    'Coord2': 'Coord2',
    'TimeStep': 'timestep',
    'timestep': 'timestep',
    'v_t': 'time',
    'Ncount': 'Ncount',
    'n_contact': 'n_contact',
    'c_m1': 'mass',
    'v_mv1': 'mv_1',
    'v_mv2': 'mv_2',
    'v_mv3': 'mv_3',
    'v_mvx': 'mv_1',
    'v_mvy': 'mv_2',
    'v_mvz': 'mv_3',
    'v_Ek1': 'Ek_1',
    'v_Ek2': 'Ek_2',
    'v_Ek3': 'Ek_3',
    'c_omega[1]': 'omega_1',
    'c_omega[2]': 'omega_2',
    'c_omega[3]': 'omega_3',
    'c_stress[1]': 'stress_multiply_binvolume_11',
    'c_stress[2]': 'stress_multiply_binvolume_22',
    'c_stress[3]': 'stress_multiply_binvolume_33',
    'c_stress[4]': 'stress_multiply_binvolume_12',
    'c_stress[5]': 'stress_multiply_binvolume_13',
    'c_stress[6]': 'stress_multiply_binvolume_23',
    'v_inwall_per_atom_1': 'chunk_inwall_force_1',
    'v_inwall_per_atom_2': 'chunk_inwall_force_2',
    'v_inwall_per_atom_3': 'chunk_inwall_force_3',
    'v_outwall_per_atom_1': 'chunk_outwall_force_1',
    'v_outwall_per_atom_2': 'chunk_outwall_force_2',
    'v_outwall_per_atom_3': 'chunk_outwall_force_3',
    'v_zbottom_per_atom_1': 'chunk_zbottom_force_1',
    'v_zbottom_per_atom_2': 'chunk_zbottom_force_2',
    'v_zbottom_per_atom_3': 'chunk_zbottom_force_3',
    'v_force_y_bottom_1': 'inwall_force_1',
    'v_force_inwall_1': 'inwall_force_1',
    'v_force_y_bottom_x': 'inwall_force_1',
    'v_force_y_bottom_2': 'inwall_force_2',
    'v_force_inwall_2': 'inwall_force_2',
    'v_force_y_bottom_y': 'inwall_force_2',
    'v_force_y_bottom_3': 'inwall_force_3',
    'v_force_inwall_3': 'inwall_force_3',
    'v_force_y_bottom_z': 'inwall_force_3',
    'v_force_y_top_1': 'outwall_force_1',
    'v_force_outwall_1': 'outwall_force_1',
    'v_force_y_top_x': 'outwall_force_1',
    'v_force_y_top_2': 'outwall_force_2',
    'v_force_outwall_2': 'outwall_force_2',
    'v_force_y_top_y': 'outwall_force_2',
    'v_force_y_top_3': 'outwall_force_3',
    'v_force_outwall_3': 'outwall_force_3',
    'v_force_y_top_z': 'outwall_force_3',
    'v_force_zbottom_1': 'zbottom_force_1',
    'v_force_z_bottom_1': 'zbottom_force_1',
    'v_force_zbottom_x': 'zbottom_force_1',
    'v_force_z_bottom_x': 'zbottom_force_1',
    'v_force_zbottom_2': 'zbottom_force_2',
    'v_force_z_bottom_2': 'zbottom_force_2',
    'v_force_zbottom_y': 'zbottom_force_2',
    'v_force_z_bottom_y': 'zbottom_force_2',
    'v_force_zbottom_3': 'zbottom_force_3',
    'v_force_z_bottom_3': 'zbottom_force_3',
    'v_force_zbottom_z': 'zbottom_force_3',
    'v_force_z_bottom_z': 'zbottom_force_3',
}

v_name_after_rename_list = [
    dic_rename_fixtimeave_npy_headername_to_use[key] for key in dic_rename_fixtimeave_npy_headername_to_use
]
def map_name_to_sq_name(name):
    return name + "_sq"

def map_name_to_std_name(name):
    return name + "_std"
dic_rename_fixtimeave_npy_headername_to_use_sq = {}
for key in dic_rename_fixtimeave_npy_headername_to_use:
    dic_rename_fixtimeave_npy_headername_to_use_sq[map_name_to_sq_name(key)] = map_name_to_sq_name(
        dic_rename_fixtimeave_npy_headername_to_use[key]
    )

dic_rename_fixtimeave_npy_headername_to_use.update(dic_rename_fixtimeave_npy_headername_to_use_sq)

def sq_file_path(n, log_variable_dic_list, fixtimeave_id_name, filename):
    sq_filename = map_name_to_sq_name(eliminate_npy_if_yes(filename))
    sq_filename = add_npy_if_not(sq_filename)
    filepath = fixtimeave_npy_output_file_path(n, fixtimeave_id_name, log_variable_dic_list, sq_filename)
    return filepath


# processed folder
def processed_subfolder_path(folder_name_under_subfolder_of_data):
    folder_path = os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '03_processed',
        folder_name_under_subfolder_of_data,
    )
    return folder_path
# calculate npy folder

def npy_calculated_folder_path(n, log_variable_dic_list):
    folder_path = os.path.join(
        processed_subfolder_path(folder_name_under_subfolder_of_data(log_variable_dic_list)),
        'npy_calculated',
        'simu_' + str(n),
    )
    return folder_path

# std
def npy_calculated_std_file_path(v_name, n, log_variable_dic_list):
    file_path = os.path.join(
        npy_calculated_folder_path(n, log_variable_dic_list),
        add_npy_if_not(v_name),
    )
    return file_path

# calculate variable name to filename
calculated_coord_filename = {
    # coord
    "coord_1_middle_23": "coord_1_middle_23",
    "coord_2_middle_23": "coord_2_middle_23",
}
def npy_calculated_coord_file_path(v_name, n, log_variable_dic_list):
    file_path = os.path.join(
        npy_calculated_folder_path(n, log_variable_dic_list),
        add_npy_if_not(calculated_coord_filename[v_name]),
    )
    return file_path

# calculate variable name to filename
calculated_variable_name_map_to_filename = {
    # velocity
    "velocity_1": "velocity_1",
    "velocity_2": "velocity_2",
    "velocity_3": "velocity_3",
    # strain rate
    "strain_rate_11": "strain_rate_11",
    "strain_rate_22": "strain_rate_22",
    "strain_rate_33": "strain_rate_33",
    "strain_rate_12": "strain_rate_12",
    "strain_rate_13": "strain_rate_13",
    "strain_rate_23": "strain_rate_23",
    "strain_rate_21": "strain_rate_21",
    "strain_rate_31": "strain_rate_31",
    "strain_rate_32": "strain_rate_32",
    # middle strain_rate
    "strain_rate_11_middle": "strain_rate_11_middle",
    "strain_rate_22_middle": "strain_rate_22_middle",
    "strain_rate_33_middle": "strain_rate_33_middle",
    "strain_rate_12_middle": "strain_rate_12_middle",
    "strain_rate_13_middle": "strain_rate_13_middle",
    "strain_rate_23_middle": "strain_rate_23_middle",
    "strain_rate_21_middle": "strain_rate_21_middle",
    "strain_rate_31_middle": "strain_rate_31_middle",
    "strain_rate_32_middle": "strain_rate_32_middle",
    # fraction
    "fraction": "fraction",
    # wall stress
    "inwall_stress_1": "inwall_stress_1",
    "inwall_stress_2": "inwall_stress_2",
    "inwall_stress_3": "inwall_stress_3",
    "outwall_stress_1": "outwall_stress_1",
    "outwall_stress_2": "outwall_stress_2",
    "outwall_stress_3": "outwall_stress_3",
    "zbottom_stress_1": "zbottom_stress_1",
    "zbottom_stress_2": "zbottom_stress_2",
    "zbottom_stress_3": "zbottom_stress_3",
    # stress
    "stress_11": "stress_11",
    "stress_22": "stress_22",
    "stress_33": "stress_33",
    "stress_12": "stress_12",
    "stress_13": "stress_13",
    "stress_23": "stress_23",
    # pressure
    "pressure": "pressure",
    # mu
    "mu_11": "mu_11",
    "mu_22": "mu_22",
    "mu_33": "mu_33",
    "mu_12": "mu_12",
    "mu_13": "mu_13",
    "mu_23": "mu_23",
    "mu_21": "mu_21",
    "mu_31": "mu_31",
    "mu_32": "mu_32",
    # mu_ij_middle
    "mu_11_middle": "mu_11_middle",
    "mu_22_middle": "mu_22_middle",
    "mu_33_middle": "mu_33_middle",
    "mu_12_middle": "mu_12_middle",
    "mu_13_middle": "mu_13_middle",
    "mu_23_middle": "mu_23_middle",
    "mu_21_middle": "mu_21_middle",
    "mu_31_middle": "mu_31_middle",
    "mu_32_middle": "mu_32_middle",
    # mu_tensor
    "mu_tensor_11": "mu_tensor_11",
    "mu_tensor_22": "mu_tensor_22",
    "mu_tensor_33": "mu_tensor_33",
    "mu_tensor_12": "mu_tensor_12",
    "mu_tensor_13": "mu_tensor_13",
    "mu_tensor_23": "mu_tensor_23",
    "mu_tensor_21": "mu_tensor_21",
    "mu_tensor_31": "mu_tensor_31",
    "mu_tensor_32": "mu_tensor_32",
    # I
    "I_11": "I_11",
    "I_22": "I_22",
    "I_33": "I_33",
    "I_12": "I_12",
    "I_13": "I_13",
    "I_23": "I_23",
    "I_21": "I_21",
    "I_31": "I_31",
    "I_32": "I_32",
    # I_tensor
    "I_tensor": "I_tensor",
}

def npy_calculated_variable_file_path(v_name, n, log_variable_dic_list):
    file_path = os.path.join(
        npy_calculated_folder_path(n, log_variable_dic_list),
        add_npy_if_not(calculated_variable_name_map_to_filename[v_name]),
    )
    return file_path

def v_name_to_path(n, v_name, log_variable_dic_list, fixtimeave_id_name=None, is_std=False, is_calculated_v=False):
    if fixtimeave_id_name is not None:
        path = fixtimeave_npy_output_file_path(n, fixtimeave_id_name, log_variable_dic_list, v_name)
    elif is_std:
        path = npy_calculated_std_file_path(v_name, n, log_variable_dic_list)
    elif is_calculated_v:
        path = npy_calculated_variable_file_path(v_name, n, log_variable_dic_list)
    else:
        sys.exit("input wrong for fixtimeave_id_name or is_std or is_calculated_v")
    return path

view_folder_path = os.path.join(
    os.path.expanduser("~"),
    'simulation',
    'lammps_process',
    'data',
    '07_plots',
    'view',
)
paper_folder_path = os.path.join(
    os.path.expanduser("~"),
    'simulation',
    'lammps_process',
    'data',
    '07_plots',
    'paper',
)
def plots_for_view_folder(log_variable_dic_list):
    folder_path = os.path.join(
        view_folder_path,
        folder_name_under_subfolder_of_data(log_variable_dic_list),
    )
    return folder_path

def plots_for_paper_folder(log_variable_dic_list):
    folder_path = os.path.join(
        paper_folder_path,
        folder_name_under_subfolder_of_data(log_variable_dic_list),
    )
    return folder_path

def plots_for_view_file_path(log_variable_dic_list, figurename):
    file_path = os.path.join(
        plots_for_view_folder(log_variable_dic_list),
        figurename,
    )
    return file_path

def plots_for_paper_file_path(log_variable_dic_list, figurename):
    file_path = os.path.join(
        plots_for_paper_folder(log_variable_dic_list),
        figurename,
    )
    return file_path

def latex_folder_path(log_variable_dic_list):
    folder_path = os.path.join(
        os.path.expanduser("~"),
        'simulation',
        'lammps_process',
        'data',
        '08_latex',
        folder_name_under_subfolder_of_data(log_variable_dic_list),
    )
    return folder_path

def latex_for_view_file_path(log_variable_dic_list):
    path = os.path.join(
        latex_folder_path(log_variable_dic_list),
        'view' + folder_name_under_subfolder_of_data(log_variable_dic_list),
    )
    return path

def latex_for_paper_file_path(log_variable_dic_list):
    path = os.path.join(
        latex_folder_path(log_variable_dic_list),
        'paper' + folder_name_under_subfolder_of_data(log_variable_dic_list),
    )
    return path

###### plot


###### sync to overleaf dropbox
dropbox_thesis_figures_folder_path = '/home/ic6413/Dropbox/應用程式/Overleaf/Thesis Figures/images'
dropbox_thesis_figures_folder_path_simulation = os.path.join(
    dropbox_thesis_figures_folder_path,
    'simulation',
)
def plots_for_sync_paper_folder(log_variable_dic_list):
    folder_path = os.path.join(
        dropbox_thesis_figures_folder_path_simulation,
        folder_name_under_subfolder_of_data(log_variable_dic_list),
    )
    return folder_path