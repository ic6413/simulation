import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr

def calculate_std_by_ave_and_sq_ave(n_sample, ave, ave_sq):
    if np.any(ave_sq<0):
        pass
        #breakpoint()
    std_sq = ave_sq - ave**2
    std_sq[np.logical_and(-10**-20 < ave_sq, ave_sq < 0)] = 0
    std_sq[(ave_sq==0)]=0
    ratio = np.divide(std_sq, ave_sq, out=np.zeros_like(std_sq), where=ave_sq!=0)
    filter_condition_to_zero = np.logical_and(-10**-6 < ratio, ratio < 0)
    if np.any(ratio[np.logical_not(filter_condition_to_zero)]<0):
        pass
        #breakpoint()
    std_sq[filter_condition_to_zero]=0
    std = (
        std_sq
    )**0.5
    if np.any(std_sq < 0):
        A = ave_sq[std_sq<0]
        print(A)
        #breakpoint()
    return std

def calculate_std_and_save(n, log_variable_dic_list, filename, filename_sq, fixtimeave_id_name, out_file_name):
    # calculate the std by read value and value_sq
    # save it in the same npy folder that read from
    ave = np.load(di.fixtimeave_npy_output_file_path(n, fixtimeave_id_name, log_variable_dic_list, filename), mmap_mode='r')
    ave_sq = np.load(di.fixtimeave_npy_output_file_path(n, fixtimeave_id_name, log_variable_dic_list, filename_sq), mmap_mode='r')
    
    maskto0 = np.logical_and(
        ave_sq > -10**-50, ave_sq < 0,
    )
    ave_sq_revised = np.copy(ave_sq)
    ave_sq_revised[maskto0] = 0
    if np.any(ave_sq_revised<0):
        A = ave_sq_revised[np.logical_and(-10**-20 < ave_sq_revised, ave_sq_revised < 0)]
        print(A)
        B = ave_sq_revised[np.logical_and(-10**-6 < ave_sq_revised, ave_sq_revised < 0)]
        print(B)
        #breakpoint()
    
    #breakpoint()
    n_sample = int(log_variable_dic_list[n]['fixavetime'][fixtimeave_id_name]['n_repeat'])
    std = calculate_std_by_ave_and_sq_ave(n_sample, ave, ave_sq_revised)
    np.save(di.npy_calculated_std_file_path(out_file_name, n, log_variable_dic_list), std)

def multi_calculate_std_and_save(log_variable_dic_list):
    for n in range(len(log_variable_dic_list)):
        for fixtimeave_id_name in di.map_fixtimeave_value_to_coord_by_id:
            folder_path = di.fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list)
            for root, dirs, files in os.walk(folder_path):
                for name in files:
                    name = di.eliminate_npy_if_yes(name)
                    filepath_sq = di.sq_file_path(n, log_variable_dic_list, fixtimeave_id_name, name)
                    if os.path.exists(filepath_sq):
                        calculate_std_and_save(
                            n, log_variable_dic_list, name, di.map_name_to_sq_name(name), fixtimeave_id_name, di.map_name_to_std_name(name),
                        )


def propagation_of_std_plus_or_minus(a, std_a, b, std_b):
    value = (
        (std_a)**2 + (std_b)**2
    )**0.5
    return value
def propagation_of_std_multi(a, std_a, b, std_b):
    value = a*b*(
        (std_a/a)**2 + (std_b/b)**2
    )**0.5
    return value
def propagation_of_std_divide(a, std_a, b, std_b):
    value = a/b*(
        (std_a/a)**2 + (std_b/b)**2
    )**0.5
    return value

def save_velocity_by_mv(n, in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass'):
    in_file_path_1 = di.fixtimeave_npy_output_file_path(n, 'avspatial_ave', log_variable_dic_list, in_name_mv)
    mass_file_path = di.fixtimeave_npy_output_file_path(n, 'avspatial_ave', log_variable_dic_list, mass_name)
    if os.path.exists(in_file_path_1) and os.path.exists(mass_file_path):
        out_file_path = di.npy_calculated_variable_file_path(out_v_name, n, log_variable_dic_list)
        mv_i = np.load(in_file_path_1, mmap_mode='r')
        mass = np.load(mass_file_path, mmap_mode='r')
        velocity_i = mv_i/mass
        np.save(out_file_path, velocity_i)
    else:sys.exit("in_file_path_1 and mass_file_path not exist")

def multi_save_velocity_by_mv(in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass'):
    for n in range(len(log_variable_dic_list)):
        save_velocity_by_mv(n, in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass')