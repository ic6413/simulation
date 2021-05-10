#!/usr/bin/env python3
import os
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data as dd

lmp_folder_path = di.lmp_folder_path

# pickle variable
(log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)

folder_name_under_subfolder_of_data = di.folder_name_under_subfolder_of_data(log_variable_dic_list)

intermediate_subfolder_path = di.intermediate_subfolder_path(folder_name_under_subfolder_of_data)

os.makedirs(intermediate_subfolder_path, exist_ok=True)

log_folder_path = os.path.join(intermediate_subfolder_path, di.log_output_folder_name)
os.makedirs(log_folder_path, exist_ok=True)
logpicklepath = os.path.join(log_folder_path,  di.logpicklefilename)
dr.dump_variable_save(lmp_folder_path, logpicklepath, None, log_variable_name = di.logpicklefilename)

# check if all fixavetime included
outputlist = list(di.map_fixtimeave_value_to_coord_by_id.keys()) + list(di.output_shape_map_from_id.keys()) + di.no_coord_fixtimeave

for key in log_variable['fixavetime'].keys():
    if key not in outputlist:
        breakpoint()
        sys.exit("some fixavetime not included")

# mode scalar
for fixtimeave_id_name in di.no_coord_fixtimeave:
    dd.multi_simu_save_non_coord_to_npy_scalar(
        folder_path_list_initial_to_last,
        log_variable_dic_list, fixtimeave_id_name,
        n_row_header=2,
    )

# len_in_each_dim_coord23
coord_2_3_path = os.path.join(
    lmp_folder_path,
    log_variable['fixavetime'][di.coord_2_3_fixtimeave_id_name]['file'],
)
len_in_each_dim_coord23 = dd.len_in_each_dim_coord23(coord_2_3_path, n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2)
# coord chunk
for fixtimeave_id_name in di.coord_chunk_id_list:
    dd.multi_simu_save_coord_to_npy(
    folder_path_list_initial_to_last,
    len_in_each_dim_coord23,
    log_variable_dic_list, fixtimeave_id_name,
    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    )

# chunk data
# The number of row of header
# count from 1
for fixtimeave_id_name in di.map_fixtimeave_value_to_coord_by_id.keys():
    dd.multi_simu_save_non_coord_to_npy(
    folder_path_list_initial_to_last,
    len_in_each_dim_coord23,
    log_variable_dic_list, fixtimeave_id_name,
    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    )