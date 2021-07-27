import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data_for_plot as ddfp
lmp_folder_path = di.lmp_folder_path
# pickle variable
(log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)

def check_support_weight(n_ave, inputstepsarray, error_tolerence = 0.01):

    inwall_force_1 = ddfp.get_ave_value(n_ave, 'timeav_inwall_force', 'inwall_force_1', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_inwall_force', is_std=False, is_calculated_v=False,)
    outwall_force_1 = ddfp.get_ave_value(n_ave, 'timeav_outwall_force', 'outwall_force_1', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_outwall_force', is_std=False, is_calculated_v=False,)
    zbottom_force_1 = ddfp.get_ave_value(n_ave, 'timeav_zbottom_force', 'zbottom_force_1', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_zbottom_force', is_std=False, is_calculated_v=False,)

    inwall_force_2 = ddfp.get_ave_value(n_ave, 'timeav_inwall_force', 'inwall_force_2', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_inwall_force', is_std=False, is_calculated_v=False,)
    outwall_force_2 = ddfp.get_ave_value(n_ave, 'timeav_outwall_force', 'outwall_force_2', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_outwall_force', is_std=False, is_calculated_v=False,)
    zbottom_force_2 = ddfp.get_ave_value(n_ave, 'timeav_zbottom_force', 'zbottom_force_2', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_zbottom_force', is_std=False, is_calculated_v=False,)

    inwall_force_3 = ddfp.get_ave_value(n_ave, 'timeav_inwall_force', 'inwall_force_3', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_inwall_force', is_std=False, is_calculated_v=False,)
    outwall_force_3 = ddfp.get_ave_value(n_ave, 'timeav_outwall_force', 'outwall_force_3', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_outwall_force', is_std=False, is_calculated_v=False,)
    zbottom_force_3 = ddfp.get_ave_value(n_ave, 'timeav_zbottom_force', 'zbottom_force_3', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='timeav_zbottom_force', is_std=False, is_calculated_v=False,)

    total_wall_force_from_total_1 = inwall_force_1 + outwall_force_1 + zbottom_force_1
    total_wall_force_from_total_2 = inwall_force_2 + outwall_force_2 + zbottom_force_2
    total_wall_force_from_total_3 = inwall_force_3 + outwall_force_3 + zbottom_force_3
        
    chunk_inwall_force_1 = ddfp.get_ave_value(n_ave, 'ave_std_inwall', 'chunk_inwall_force_1', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_inwall', is_std=False, is_calculated_v=False,)
    chunk_outwall_force_1 = ddfp.get_ave_value(n_ave, 'ave_std_outwall', 'chunk_outwall_force_1', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_outwall', is_std=False, is_calculated_v=False,)
    chunk_zbottom_force_1 = ddfp.get_ave_value(n_ave, 'ave_std_zbottom', 'chunk_zbottom_force_1', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_zbottom', is_std=False, is_calculated_v=False,)

    chunk_inwall_force_2 = ddfp.get_ave_value(n_ave, 'ave_std_inwall', 'chunk_inwall_force_2', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_inwall', is_std=False, is_calculated_v=False,)
    chunk_outwall_force_2 = ddfp.get_ave_value(n_ave, 'ave_std_outwall', 'chunk_outwall_force_2', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_outwall', is_std=False, is_calculated_v=False,)
    chunk_zbottom_force_2 = ddfp.get_ave_value(n_ave, 'ave_std_zbottom', 'chunk_zbottom_force_2', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_zbottom', is_std=False, is_calculated_v=False,)

    chunk_inwall_force_3 = ddfp.get_ave_value(n_ave, 'ave_std_inwall', 'chunk_inwall_force_3', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_inwall', is_std=False, is_calculated_v=False,)
    chunk_outwall_force_3 = ddfp.get_ave_value(n_ave, 'ave_std_outwall', 'chunk_outwall_force_3', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_outwall', is_std=False, is_calculated_v=False,)
    chunk_zbottom_force_3 = ddfp.get_ave_value(n_ave, 'ave_std_zbottom', 'chunk_zbottom_force_3', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='ave_std_zbottom', is_std=False, is_calculated_v=False,)

    total_chunk_wall_force_from_total_1 = np.sum(chunk_inwall_force_1, axis=1) + np.sum(chunk_outwall_force_1, axis=1) + np.sum(chunk_zbottom_force_1, axis=1)
    total_chunk_wall_force_from_total_2 = np.sum(chunk_inwall_force_2, axis=1) + np.sum(chunk_outwall_force_2, axis=1) + np.sum(chunk_zbottom_force_2, axis=1)
    total_chunk_wall_force_from_total_3 = np.sum(chunk_inwall_force_3, axis=1) + np.sum(chunk_outwall_force_3, axis=1) + np.sum(chunk_zbottom_force_3, axis=1)

    mass = ddfp.get_ave_value(n_ave, 'avspatial_ave', 'mass', inputstepsarray, log_variable_dic_list, fixtimeave_id_name='avspatial_ave')
    mg = np.sum(mass, axis=(1,2))*float(log_variable['g'])
    check_wall_support_error_from_chunk_wall_force = np.abs(total_wall_force_from_total_3 - mg)/mg
    check_wall_support_error_from_total_wall_force = np.abs(total_chunk_wall_force_from_total_3 - mg)/mg


    def check_error(source_wall_force_3, string_source_wall_force):
        for n, value in enumerate(source_wall_force_3):
            if value >= error_tolerence:
                sys.exit('wall support wrong from {0}, for step {1}, error is {2}'.format(string_source_wall_force, inputstepsarray[n], value))

    check_error(check_wall_support_error_from_chunk_wall_force, 'chunk_wall_force')
    check_error(check_wall_support_error_from_total_wall_force, 'total_wall_force')