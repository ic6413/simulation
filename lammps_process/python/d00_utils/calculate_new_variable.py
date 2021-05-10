import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr

def save_velocity_by_mv(n, in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass'):
    in_file_path_1 = di.fixtimeave_npy_output_file_path(n, 'avspatial_ave', log_variable_dic_list, in_name_mv)
    mass_file_path = di.fixtimeave_npy_output_file_path(n, 'avspatial_ave', log_variable_dic_list, mass_name)
    if os.path.exists(in_file_path_1) and os.path.exists(mass_file_path):
        out_file_path = di.npy_calculated_variable_file_path(out_v_name, n, log_variable_dic_list)
        mv_i = np.load(in_file_path_1, mmap_mode='r')
        mass = np.load(mass_file_path, mmap_mode='r')
        velocity_i = mv_i/mass
        np.save(out_file_path, velocity_i)

def multi_save_velocity_by_mv(in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass'):
    for n in range(len(log_variable_dic_list)):
        save_velocity_by_mv(n, in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass')