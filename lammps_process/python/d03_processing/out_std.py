#!/usr/bin/env python3
import os
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.calculate_new_variable as dc

lmp_folder_path = di.lmp_folder_path

# pickle variable
(log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)

# create folder
for n in range(n_simu_total):
    os.makedirs(di.npy_calculated_folder_path(n, log_variable_dic_list), exist_ok=True)


# calculate std
dc.multi_calculate_std_and_save(log_variable_dic_list)