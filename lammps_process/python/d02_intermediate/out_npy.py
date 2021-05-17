#!/usr/bin/env python3
import os
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data as dd
lmp_folder_path = di.lmp_folder_path
# pickle variable
(log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)

# check if all fixavetime included
outputlist = list(di.map_fixtimeave_value_to_coord_by_id.keys()) + di.coord_chunk_id_list + di.no_coord_fixtimeave

for key in log_variable['fixavetime'].keys():
    if key not in outputlist:
        sys.exit("fixavetime id " + key + " is not included in outputlist")

# mode scalar
for fixtimeave_id_name in di.no_coord_fixtimeave:
    dd.multi_simu_save_non_coord_to_npy_scalar(
        folder_path_list_initial_to_last,
        log_variable_dic_list, fixtimeave_id_name,
        n_row_header=2,
    )
if di.coord_chunk_id_23 in log_variable['fixavetime']:
    # len_in_each_dim_coord23
    coord_2_3_path = os.path.join(
        lmp_folder_path,
        log_variable['fixavetime'][di.coord_chunk_id_23]['file'],
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
else:
    # len_in_each_dim_coord23
    coord_2_3_path = os.path.join(
        lmp_folder_path,
        log_variable['fixavechunk'][di.coord_chunk_id_23_replace]['file'],
    )
    len_in_each_dim_coord23 = dd.len_in_each_dim_coord23(coord_2_3_path, n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2)
    # coord chunk
    dd.multi_simu_save_coord_23_to_npy_replace(
        folder_path_list_initial_to_last,
        len_in_each_dim_coord23,
        log_variable_dic_list,
        n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    )
    dd.save_inwall_coord_from_chunk23(log_variable_dic_list)
    dd.save_outwall_coord_from_chunk23(log_variable_dic_list)
    dd.save_zbottom_coord_from_chunk23(log_variable_dic_list)
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