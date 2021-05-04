import os
lmp_folder_path = os.path.join(
    os.path.expanduser("~"),
    'simulation',
    'lammps_process',
    'data',
    '01_raw',
    'Link to lmp_run',
    'test',
    '20200511',
    'f_5e6',
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

def post_process_folder_path(folder_name_under_subfolder_of_data):
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
def npy_folder_path(n, log_variable_dic_list):
    folder_path = os.path.join(
        post_process_folder_path(folder_name_under_subfolder_of_data(log_variable_dic_list)),
        'npy',
        'simu_' + str(n),
    )
    return folder_path

def fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list):
    folder_path = os.path.join(
        npy_folder_path(n, log_variable_dic_list),
        npy_output_subfolder_name_map_from_id[fixtimeave_id_name],
    )
    return folder_path

def fixtimeave_text_file_path(n, fixtimeave_id_name, log_variable_dic_list, folder_path_list_initial_to_last):
    folder_path = os.path.join(
        folder_path_list_initial_to_last[n],
        log_variable_dic_list[n]['fixavetime'][fixtimeave_id_name]['file'],
    )
    return folder_path