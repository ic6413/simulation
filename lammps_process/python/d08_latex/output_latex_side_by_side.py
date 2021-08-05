#!/usr/bin/env python3
import os
import d00_utils.input_text as di
import d00_utils.latex as dl
import d00_utils.read_log as dr

lmp_folder_path_list = [
    di.lmp_folder_path_dic[1],
    di.lmp_folder_path_dic[2],
    di.lmp_folder_path_dic[3],
]
def get_list_of_log_dic_list(lmp_folder_path_list):

    list_of_log_variable_dic_list = []
    for lmp_folder_path in lmp_folder_path_list:
        (log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)
        list_of_log_variable_dic_list.append(log_variable_dic_list)
    return list_of_log_variable_dic_list

def main():
    lmp_folder_path = di.lmp_folder_path
    # pickle variable
    (log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)
    default_folder = os.path.join(
        di.latex_folder_path_08,
        'multi',
    )
    os.makedirs(default_folder, exist_ok=True)
    default_filepath = os.path.join(
        di.latex_folder_path_08,
        'multi',
        'multi',
    )
    pic_rela_path_list = dl.plots_rela_path_list(log_variable_dic_list, di.plots_for_view_folder(log_variable_dic_list))
    list_of_log_variable_dic_list = get_list_of_log_dic_list(lmp_folder_path_list)
    dl.produce_latex_side_by_side(list_of_log_variable_dic_list, pic_rela_path_list, default_filepath, False)


if __name__ == "__main__":
    main()