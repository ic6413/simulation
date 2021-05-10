import os
lmp_folder_path = os.path.join(
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
)

log_output_folder_name = "log"
logpicklefilename = 'log.pickle'

def folder_name_under_subfolder_of_data(log_variable_dic_list):
    # L W H Sa rotate_start_time endtime history wallshape xmu hooke/hertz kn kt gamma_n gamma_t timestep
    variable_name = [

    ]
    variable_label_string = {

    }
    name = 'test_subfoldername'
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

coord_2_3_fixtimeave_id_name = 'coord1and2_chunk_2_3'

coord_chunk_id_list = [
    'coord1and2_chunk_2_3',
    'coord1and2_chunk_near_inwall',
    'coord1and2_chunk_near_outwall',
    'coord1and2_chunk_near_zbottom',
]

map_fixtimeave_value_to_coord_by_id = {
    'avspatial_ave': 'coord1and2_chunk_2_3',
    'avspatial_omega_ave': 'coord1and2_chunk_2_3',
    'avspatialstress_ave': 'coord1and2_chunk_2_3',
    'contact_ave': 'coord1and2_chunk_2_3',
    'ave_std_inwall': 'coord1and2_chunk_near_inwall',
    'ave_std_outwall': 'coord1and2_chunk_near_outwall',
    'ave_std_zbottom': 'coord1and2_chunk_near_zbottom',
}

no_coord_fixtimeave = [
    'timeav_inwall_force',
    'timeav_outwall_force',
    'timeav_zbottom_force',
]

npy_output_subfolder_name_map_from_id = {
    'coord1and2_chunk_2_3': "coord_chunk_2_3",
    'coord1and2_chunk_near_inwall': "coord_chunk_inwall",
    'coord1and2_chunk_near_outwall': "coord_chunk_outwall",
    'coord1and2_chunk_near_zbottom': "coord_chunk_zbottom",
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
}

output_shape_map_from_id = {
    'coord1and2_chunk_2_3':  ['n_1', 'n_2'],
    'coord1and2_chunk_near_inwall': ['n_2'],
    'coord1and2_chunk_near_outwall': ['n_2'],
    'coord1and2_chunk_near_zbottom': ['n_1'],
    'avspatial_ave': ['t', 'n_1', 'n_2'],
    'avspatial_omega_ave': ['t', 'n_1', 'n_2'],
    'avspatialstress_ave': ['t', 'n_1', 'n_2'],
    'ave_std_inwall': ['t', 'n_2'],
    'ave_std_outwall': ['t', 'n_2'],
    'ave_std_zbottom': ['t', 'n_1'],
}

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

dic_rename_fixtimeave_npy_headername_to_use = {
    'TimeStep': 'timestep',
    'timestep': 'timestep',
    'Ncount': 'Ncount',
    'n_contact': 'n_contact',
    'c_m1': 'mass',
    'v_mv1': 'mv_1',
    'v_mv2': 'mv_2',
    'v_mv3': 'mv_3',
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

dic_rename_fixtimeave_npy_headername_to_use_sq = {}
for key in dic_rename_fixtimeave_npy_headername_to_use:
    dic_rename_fixtimeave_npy_headername_to_use_sq[key + "_sq"] = dic_rename_fixtimeave_npy_headername_to_use[key] + "_sq"

dic_rename_fixtimeave_npy_headername_to_use.update(dic_rename_fixtimeave_npy_headername_to_use_sq)

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
calculated_variable_filename_map_to_filename = {
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
        add_npy_if_not(calculated_variable_filename_map_to_filename[v_name]),
    )
    return file_path


# stress
# pressure
# strain rate

# mu_tensor
# mu
# I_tensor
# I



calculated_variable_name_list = [

]