#!/usr/bin/env python3
import os
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data as dd

def main():
    # check if there are diff path map to the same subfolder name
    folder_name_dic = {}
    for key in di.lmp_folder_path_dic:
        check_lmp_folder_path = di.lmp_folder_path_dic[key]

        # pickle variable
        (log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(check_lmp_folder_path)

        folder_name_under_subfolder_of_data = di.folder_name_under_subfolder_of_data(log_variable_dic_list)

        folder_name_dic[key] = folder_name_under_subfolder_of_data
    if len(set(folder_name_dic.values()))!=len(list(folder_name_dic.values())):
        sys.exit('There are duplicates value for folder_name_dic')
    else:
        print('No duplicates value for folder_name_dic')

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

if __name__ == "__main__":
    main()
