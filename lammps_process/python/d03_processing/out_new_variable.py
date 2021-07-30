#!/usr/bin/env python3
import os
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.calculate_new_variable as dc

def main():
    lmp_folder_path = di.lmp_folder_path

    # pickle variable
    (log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)

    # create folder
    for n in range(n_simu_total):
        os.makedirs(di.npy_calculated_folder_path(n, log_variable_dic_list), exist_ok=True)


    # calculate wall stress and save

    # calculate velocity and save
    dc.multi_save_velocity_by_mv('mv_1', 'velocity_1', log_variable_dic_list, mass_name='mass')
    dc.multi_save_velocity_by_mv('mv_2', 'velocity_2', log_variable_dic_list, mass_name='mass')
    dc.multi_save_velocity_by_mv('mv_3', 'velocity_3', log_variable_dic_list, mass_name='mass')

    # calculate fraction
    dc.multi_save_fraction_by_mass('fraction', log_variable_dic_list, mass_name='mass')

if __name__ == "__main__":
    main()