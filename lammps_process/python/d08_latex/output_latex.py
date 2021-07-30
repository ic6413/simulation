#!/usr/bin/env python3

import d00_utils.input_text as di
import d00_utils.latex as dl
import d00_utils.read_log as dr

def main():
    lmp_folder_path = di.lmp_folder_path
    # pickle variable
    (log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)

    dl.produce_latex_paper(log_variable_dic_list)

    dl.produce_latex_view(log_variable_dic_list)

if __name__ == "__main__":
    main()