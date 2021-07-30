import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data as dd

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

def vol_in_chunks(log_variable):
    n_1 = int(log_variable['N_bin_y'])
    n_2 = int(log_variable['N_bin_z'])
    dx = 1/n_1*int(log_variable['width_wall_dp_unit'])
    dy = 1/n_2*float(log_variable['zhi_chunk_dp_unit'])
    vol_in_chunks = float(log_variable['dp'])**3*int(log_variable['x_period_dp_unit'])*dx*dy
    return vol_in_chunks

def save_fraction_by_mass(n, out_v_name, log_variable_dic_list, mass_name='mass'):
    mass_file_path = di.fixtimeave_npy_output_file_path(n, 'avspatial_ave', log_variable_dic_list, mass_name)
    if os.path.exists(mass_file_path):
        out_file_path = di.npy_calculated_variable_file_path(out_v_name, n, log_variable_dic_list)
        mass = np.load(mass_file_path, mmap_mode='r')
        fraction = mass/float(log_variable_dic_list[-1]['den'])/vol_in_chunks(log_variable_dic_list[-1])
        np.save(out_file_path, fraction)
    else:
        sys.exit(
            "mass_file_path not exist for n=" + str(n)
            + ", mass_file_path is " + mass_file_path
        )

def save_velocity_by_mv(n, in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass'):
    in_file_path_1 = di.fixtimeave_npy_output_file_path(n, 'avspatial_ave', log_variable_dic_list, in_name_mv)
    mass_file_path = di.fixtimeave_npy_output_file_path(n, 'avspatial_ave', log_variable_dic_list, mass_name)
    if os.path.exists(in_file_path_1) and os.path.exists(mass_file_path):
        out_file_path = di.npy_calculated_variable_file_path(out_v_name, n, log_variable_dic_list)
        mv_i = np.load(in_file_path_1, mmap_mode='r')
        mass = np.load(mass_file_path, mmap_mode='r')
        velocity_i = mv_i/mass
        np.save(out_file_path, velocity_i)
    else:
        sys.exit(
            "in_file_path_1 and mass_file_path not exist for n=" + str(n)
            + ", in_file_path_1 is " + in_file_path_1 + ", mass_file_path is " + mass_file_path
        )

# stress correction due to average velocity in the chunk
def stress_correction_by_ave_velocity(i, j, n, log_variable_dic_list, mass_name='mass'):
    
    in_name_mv_i = "mv_" + str(i)
    in_name_mv_j = "mv_" + str(j)
    mv_i = dd.get_variable(n, in_name_mv_i, log_variable_dic_list, fixtimeave_id_name='avspatial_ave', is_std=False, is_calculated_v=False)
    mv_j = dd.get_variable(n, in_name_mv_j, log_variable_dic_list, fixtimeave_id_name='avspatial_ave', is_std=False, is_calculated_v=False)
    mass = dd.get_variable(n, mass_name, log_variable_dic_list, fixtimeave_id_name='avspatial_ave', is_std=False, is_calculated_v=False)
    stress_correction = np.divide(mv_i*mv_j, mass, out=np.zeros_like(mv_i*mv_j), where=mass!=0)
    
    return stress_correction

def stress_fix_by_velocity_and_wall(i, j, n, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    if (i, j) not in [
        (1, 1),
        (2, 2),
        (3, 3),
        (1, 2),
        (1, 3),
        (2, 3),
    ]:
        sys.exit("i j not correct")
    in_name = 'stress_multiply_binvolume' + "_" + str(i) + str(j)
    stress_multiply_binvolume = dd.get_variable(
        n, in_name, log_variable_dic_list, fixtimeave_id_name='avspatialstress_ave',
        is_std=False, is_calculated_v=False,
    )
    
    if ifcorrect_by_ave_velocity:
        stress_multiply_binvolume = stress_multiply_binvolume + stress_correction_by_ave_velocity(i, j, n, log_variable_dic_list)
    
    coord1 = dd.get_coord_by_variable(n, "Coord1", log_variable_dic_list, 'avspatial_ave')
    coord2 = dd.get_coord_by_variable(n, "Coord2", log_variable_dic_list, 'avspatial_ave')
    d_coord1 = coord1[1,0] - coord1[0,0]
    d_coord2 = coord2[0,1] - coord2[0,0]
    log_variable = log_variable_dic_list[-1]
    if ifcorrect_by_wall:
        if i==3 and j==3:
            chunk_zbottom_force_3 = dd.get_variable(
                n, 'chunk_zbottom_force_3', log_variable_dic_list, fixtimeave_id_name='ave_std_zbottom',
                is_std=False, is_calculated_v=False,
            )
            stress_volume_bottom_contribution = -float(log_variable['dp'])*chunk_zbottom_force_3/2
            stress_multiply_binvolume[:, :, 0] += stress_volume_bottom_contribution
            
        if i==2 and j==2:
            chunk_inwall_force_2 = dd.get_variable(n, 'chunk_inwall_force_2', log_variable_dic_list, fixtimeave_id_name='ave_std_inwall', is_std=False, is_calculated_v=False,)
            stress_multiply_binvolume[:, 0, :] += -float(log_variable['dp'])*chunk_inwall_force_2/2

            chunk_outwall_force_2 = dd.get_variable(n, 'chunk_outwall_force_2', log_variable_dic_list, fixtimeave_id_name='ave_std_outwall', is_std=False, is_calculated_v=False,)
            stress_multiply_binvolume[:, -1, :] += -1*(-float(log_variable['dp']))*chunk_outwall_force_2/2
        if i==2 and j==3:
            chunk_inwall_force_3 = dd.get_variable(n, 'chunk_inwall_force_3', log_variable_dic_list, fixtimeave_id_name='ave_std_inwall', is_std=False, is_calculated_v=False,)
            stress_multiply_binvolume[:, 0, :] += -float(log_variable['dp'])*chunk_inwall_force_3/2

            chunk_outwall_force_3 = dd.get_variable(n, 'chunk_outwall_force_3', log_variable_dic_list, fixtimeave_id_name='ave_std_outwall', is_std=False, is_calculated_v=False,)
            stress_multiply_binvolume[:, -1, :] += -1*(-float(log_variable['dp']))*chunk_outwall_force_3/2
    binvolume = (
        d_coord1
        *d_coord2
        *float(log_variable_dic_list[-1]['x_period_dp_unit'])
        *float(log_variable_dic_list[-1]['dp'])
    )
    stress = stress_multiply_binvolume/binvolume
    return stress

def ave_stress_fix_by_velocity_and_wall(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, fixid_for_timestep='avspatialstress_ave', ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    def fun(i, j, n, log_variable_dic_list):
        return stress_fix_by_velocity_and_wall(i, j, n, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    return ave_n(fun, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, i, j, log_variable_dic_list, fixid_for_timestep)

def ave_4_grid_for_the_last_axis(value_array):
    value_array_middle = (value_array[:, :-1, :-1] + value_array[:, 1:, :-1] + value_array[:, :-1, 1:] + value_array[:, 1:, 1:])/4
    return value_array_middle

def stress_middle_fix_by_velocity_and_wall(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    stress = ave_stress_fix_by_velocity_and_wall(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    stress = ave_4_grid_for_the_last_axis(stress)
    return stress

def pressure_middle(
        n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list,
        ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True,
    ):
    pressure = 0
    for (i,j) in [
        (1,1),
        (2,2),
        (3,3),
    ]:
        pressure += -1/3*stress_middle_fix_by_velocity_and_wall(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    return pressure

def stress_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    stress_tensor_deviator = stress_middle_fix_by_velocity_and_wall(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    if i == j:
        stress_tensor_deviator += pressure_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    return stress_tensor_deviator

def find_index_to_read_from(inputsteps, fixid_for_timestep, log_variable_dic_list, n_ave):
    input_first_step = inputsteps[0]
    for i in range(len(log_variable_dic_list)):
        timestep = dd.get_timestep(len(log_variable_dic_list)-1-i, fixid_for_timestep, log_variable_dic_list)

        if input_first_step >= timestep[int((n_ave-1)/2)]:
            break
    index_to_read_from = len(log_variable_dic_list) - 1 -i
    return index_to_read_from

def average_over_n_ave(n_ave, array, indexesarray):
    sum_over_n_array = 0
    m = int((n_ave-1)/2)
    for k in range(-m, m+1):
        sum_over_n_array = sum_over_n_array + array[indexesarray + k]
    ave_array = sum_over_n_array/n_ave
    return ave_array

def stress_abs_sum_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=False, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    stress_abs_sum = 0
    if not ignore_diagonal:
        for (i,j) in [
            (1,1),
            (2,2),
            (3,3),
        ]:
            stress_abs_sum += 1/2*(stress_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall))**2

    for (i,j) in [
        (1,2),
        (1,3),
        (2,3),
    ]:
        stress_abs_sum += 2*1/2*(stress_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall))**2
    
    stress_abs_sum = stress_abs_sum**0.5
    return stress_abs_sum

def d_vi_d_x2(i, n, log_variable_dic_list):
    mass = dd.get_variable(n, 'mass', log_variable_dic_list, fixtimeave_id_name='avspatial_ave', is_std=False, is_calculated_v=False)
    mv = dd.get_variable(n, 'mv_' + str(i), log_variable_dic_list, fixtimeave_id_name='avspatial_ave', is_std=False, is_calculated_v=False)
    velocity = mv/mass
    Coord1 = dd.get_coord_by_variable(n, "Coord1", log_variable_dic_list, 'avspatial_ave')
    dvdx = (velocity[:, :-1, :] - velocity[:, 1:, :])/(Coord1[:-1, :] - Coord1[ 1:, :])
    dvdx = 1/2*(dvdx[:, :, :-1] + dvdx[:, :, 1:])
    return dvdx

def d_vi_d_x3(i, n, log_variable_dic_list):
    mass = dd.get_variable(n, 'mass', log_variable_dic_list, fixtimeave_id_name='avspatial_ave', is_std=False, is_calculated_v=False)
    mv = dd.get_variable(n, 'mv_' + str(i), log_variable_dic_list, fixtimeave_id_name='avspatial_ave', is_std=False, is_calculated_v=False)
    velocity = mv/mass
    Coord2 = dd.get_coord_by_variable(n, "Coord2", log_variable_dic_list, 'avspatial_ave')
    dvdx = (velocity[:, :, :-1] - velocity[:, :, 1:])/(Coord2[:, :-1] - Coord2[:, 1:])
    dvdx = 1/2*(dvdx[:, :-1, :] + dvdx[:, 1:, :])
    return dvdx

def d_vi_d_xj(i, j, n, log_variable_dic_list):
    if j == 1:
        dvdx = 0
    elif j == 2:
        dvdx = d_vi_d_x2(i, n, log_variable_dic_list)
    elif j == 3:
        dvdx = d_vi_d_x3(i, n, log_variable_dic_list)
    else:
        sys.exit('j not 1 2 3')
    return dvdx

def ave_n_only_time(fun, n_ave, inputstepsarray, i, j, log_variable_dic_list, fixid_for_timestep):
    if n_ave % 2 == 1:
        pass
    else:
        sys.exit("n_ave should be odd integer")
    index_to_read_from = find_index_to_read_from(inputstepsarray, fixid_for_timestep, log_variable_dic_list, n_ave)
    timestep = dd.get_timestep(index_to_read_from, fixid_for_timestep, log_variable_dic_list)
    
    array = fun(i, j, index_to_read_from, log_variable_dic_list)
    if isinstance(array, (list, tuple, np.ndarray)):
        if index_to_read_from + 1 < len(log_variable_dic_list):
            for k in range(index_to_read_from+1, len(log_variable_dic_list)):
                timestep_append = dd.get_timestep(k, fixid_for_timestep, log_variable_dic_list)
                array_append = fun(i, j, k, log_variable_dic_list)
                timestep = np.append(timestep, timestep_append, axis=0)
                array = np.append(array, array_append, axis=0)
        d_step = timestep[1] - timestep[0]
        indexesarray = (inputstepsarray-timestep[0])/d_step
        indexesarray = indexesarray.astype(int)
        ave_array = average_over_n_ave(n_ave, array, indexesarray)
    else:
        ave_array = array
    return ave_array

def ave_coord1_coord2(array, n_ave_coord1, n_ave_coord2, arraytype='tc1c2'):
    if isinstance(array, (list, tuple, np.ndarray)):
        if arraytype == 'tc1c2':
            if len(array.shape) != 3:
                sys.exit('array dim is should be 3 but input is {}'.format(len(array.shape)))
            output_array = 0
            for n_c1 in range(n_ave_coord1):
                for n_c2 in range(n_ave_coord2):
                    index = np.s_[:, n_c1: (n_c1 + array.shape[1] - n_ave_coord1 + 1), n_c2: (n_c2 + array.shape[2] - n_ave_coord2 + 1)]
                    output_array += array[index]
            output_array = output_array/(n_ave_coord1*n_ave_coord2)
    else:
        output_array = array
    
    return output_array

def ave_n(fun, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, i, j, log_variable_dic_list, fixid_for_timestep, arraytype='tc1c2'):
    array = ave_n_only_time(fun, n_ave, inputstepsarray, i, j, log_variable_dic_list, fixid_for_timestep)
    array = ave_coord1_coord2(array, n_ave_coord1, n_ave_coord2, arraytype=arraytype)
    return array

def ave_d_vi_d_xj(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, fixid_for_timestep='avspatial_ave'):
    return ave_n(d_vi_d_xj, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, i, j, log_variable_dic_list, fixid_for_timestep)

def strain_rate_tensor(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list):
    strain_rate_tensor = 1/2*(
        ave_d_vi_d_xj(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list)
        + ave_d_vi_d_xj(j, i, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list)
    )
    return strain_rate_tensor

def strain_rate_diagonal_ave(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list):
    strain_rate_diagonal_ave = 0
    for (i,j) in [
        (1,1),
        (2,2),
        (3,3),
    ]:
        strain_rate_diagonal_ave += 1/3*strain_rate_tensor(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list)
    return strain_rate_diagonal_ave

def strain_rate_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list):
    strain_rate_tensor_deviator = strain_rate_tensor(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list)
    if i == j:
        strain_rate_tensor_deviator -= strain_rate_diagonal_ave(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list)
    return strain_rate_tensor_deviator

def strain_rate_second_invariant(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=False):
    strain_rate_second_invariant = 0
    if not ignore_diagonal:
        for (i,j) in [
            (1,1),
            (2,2),
            (3,3),
        ]:
            strain_rate_second_invariant += 1/2*(strain_rate_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list))**2
    for (i,j) in [
        (1,2),
        (1,3),
        (2,3),
    ]:
        strain_rate_second_invariant += 2*1/2*(strain_rate_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list))**2

    strain_rate_second_invariant = strain_rate_second_invariant**0.5
    return strain_rate_second_invariant

def mu(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=False, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    stress = stress_abs_sum_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=ignore_diagonal, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    pressure = pressure_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    mu = stress/pressure
    return mu

def mu_ij(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    stress = stress_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    pressure = pressure_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    mu_ij = stress/pressure
    return mu_ij

def I(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=False, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    strain_rate = strain_rate_second_invariant(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=ignore_diagonal)
    log_variable = log_variable_dic_list[-1]
    d = float(log_variable['dp'])
    P = pressure_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    density = float(log_variable['den'])
    I = strain_rate*d/(P/density)**0.5
    return I

def I_ij(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    strain_rate = strain_rate_tensor_deviator(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list)
    log_variable = log_variable_dic_list[-1]
    d = float(log_variable['dp'])
    P = pressure_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    density = float(log_variable['den'])
    I_ij = strain_rate*d/(P/density)**0.5
    return I_ij

def trace_I(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True):
    strain_rate = strain_rate_diagonal_ave(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list)
    log_variable = log_variable_dic_list[-1]
    d = float(log_variable['dp'])
    P = pressure_middle(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=ifcorrect_by_ave_velocity, ifcorrect_by_wall=ifcorrect_by_wall)
    density = float(log_variable['den'])
    trace_I = strain_rate*d/(P/density)**0.5
    return trace_I

def multi_save_velocity_by_mv(in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass'):
    for n in range(len(log_variable_dic_list)):
        save_velocity_by_mv(n, in_name_mv, out_v_name, log_variable_dic_list, mass_name='mass')

def multi_save_fraction_by_mass(out_v_name, log_variable_dic_list, mass_name='mass'):
    for n in range(len(log_variable_dic_list)):
        save_fraction_by_mass(n, out_v_name, log_variable_dic_list, mass_name='mass')