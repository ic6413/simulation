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
