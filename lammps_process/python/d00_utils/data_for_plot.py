import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data as dd

def get_timestep(n, fixid_for_timestep, log_variable_dic_list):
    return dd.get_timestep(n, fixid_for_timestep, log_variable_dic_list)

def get_d_step(n, fixid_for_timestep, log_variable_dic_list):
    timestep = get_timestep(n, fixid_for_timestep, log_variable_dic_list)
    d_step = timestep[1] - timestep[0]
    return d_step
    
def get_variable(n, v_name, log_variable_dic_list, fixtimeave_id_name=None, is_std=False, is_calculated_v=False):
    return dd.get_variable(n, v_name, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_std=is_std, is_calculated_v=is_calculated_v)

def get_coord_by_variable(n, coord_name, log_variable_dic_list, fixtimeave_id_name):
    return dd.get_coord_by_variable(n, coord_name, log_variable_dic_list, fixtimeave_id_name)

def find_index_to_read_from(inputsteps, fixid_for_timestep, log_variable_dic_list, n_ave):
    input_first_step = inputsteps[0]
    for i in range(len(log_variable_dic_list)):
        timestep = get_timestep(len(log_variable_dic_list)-1-i, fixid_for_timestep, log_variable_dic_list)

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

def average_over_n_ave_std(n_ave, std_array, indexesarray):
    std_array_sq = std_array**2
    sum_over_n_array_sq = 0
    for k in range(-(n_ave-1)/2, (n_ave-1)/2+1):
        sum_over_n_array_sq = sum_over_n_array_sq + std_array_sq[indexesarray + k]
    std_ave_array = (sum_over_n_array_sq/n_ave)**0.5
    return std_ave_array

def get_array_and_index(
    n_ave, fixid_for_timestep, v_name, inputstepsarray, log_variable_dic_list,
    fixtimeave_id_name=None, is_std=False, is_calculated_v=False,
    ):
    if n_ave % 2 == 1:
        pass
    else:
        sys.exit("n_ave should be odd integer")
    index_to_read_from = find_index_to_read_from(inputstepsarray, fixid_for_timestep, log_variable_dic_list, n_ave)
    timestep = get_timestep(index_to_read_from, fixid_for_timestep, log_variable_dic_list)
    array = get_variable(index_to_read_from, v_name, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_std=is_std, is_calculated_v=is_calculated_v)

    if index_to_read_from + 1 < len(log_variable_dic_list):
        for j in range(index_to_read_from+1, len(log_variable_dic_list)):
            timestep_append = get_timestep(j, fixid_for_timestep, log_variable_dic_list)
            array_append = get_variable(j, v_name, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_std=is_std, is_calculated_v=is_calculated_v)
            timestep = np.append(timestep, timestep_append, axis=0)
            array = np.append(array, array_append, axis=0)
    d_step = timestep[1] - timestep[0]
    indexesarray = (inputstepsarray-timestep[0])/d_step
    indexesarray = indexesarray.astype(int)
    return (array, indexesarray)

def get_ave_value(n_ave, fixid_for_timestep, v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=None, is_std=False, is_calculated_v=False,):
    (array, indexesarray) = get_array_and_index(
        n_ave, fixid_for_timestep, v_name, inputstepsarray, log_variable_dic_list,
        fixtimeave_id_name=fixtimeave_id_name, is_std=is_std, is_calculated_v=is_calculated_v,
    )
    ave_array = average_over_n_ave(n_ave, array, indexesarray)
    return ave_array

def get_std_ave_value(n_ave, fixid_for_timestep, std_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=None, is_std=False, is_calculated_v=False,):
    (std_array, indexesarray) = get_array_and_index(
        n_ave, fixid_for_timestep, std_name, inputstepsarray, log_variable_dic_list,
        fixtimeave_id_name=fixtimeave_id_name, is_std=is_std, is_calculated_v=is_calculated_v,
    )
    ave_array = average_over_n_ave_std(n_ave, std_array, indexesarray)
    return ave_array

def time_from_step_0(step, log_variable):
    return step*float(log_variable["ts"])

def time_from_start_rotate(step, log_variable):
    return time_from_step_0(step, log_variable)-log_variable["rotate_start_time"]

def strain_from_rotate_start(step, log_variable):
    shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
    return time_from_start_rotate(step, log_variable)*shear_rate

def api_ave_value_for_select_region_yz(
    v_name,
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    n_y_0, d_n_y,
    n_z_0, d_n_z,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name=None,
    is_std=False,
    is_calculated_v=False,
    ):
    value = get_ave_value(n_ave, fixtimeave_id_fortime, v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_calculated_v=is_calculated_v)
    value = value[:, n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z]
    value = np.nanmean(value, axis=(1,2))
    coord_1 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord1", log_variable_dic_list, fixtimeave_id_forcoord)
    coord_2 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord2", log_variable_dic_list, fixtimeave_id_forcoord)
    first_coord_1 = coord_1[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][0, 0]
    last_coord_1 = coord_1[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][-1, 0]
    first_coord_2 = coord_2[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][0, 0]
    last_coord_2 = coord_2[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][0, -1]
    return (value, first_coord_1, last_coord_1, first_coord_2, last_coord_2)

def api_ave_ratio_value_for_select_region_yz(
    v1_name, v2_name,
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    n_y_0, d_n_y,
    n_z_0, d_n_z,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name1=None, fixtimeave_id_name2=None,
    is_calculated_v1=False, is_calculated_v2=False,
    ):
    value1 = get_ave_value(n_ave, fixtimeave_id_fortime, v1_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name1, is_calculated_v=is_calculated_v1)
    value2 = get_ave_value(n_ave, fixtimeave_id_fortime, v2_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name2, is_calculated_v=is_calculated_v2)
    value1 = value1[:, n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z]
    value2 = value2[:, n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z]
    value = np.nanmean(value1, axis=(1,2))/np.nanmean(value2, axis=(1,2))
    coord_1 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord1", log_variable_dic_list, fixtimeave_id_forcoord)
    coord_2 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord2", log_variable_dic_list, fixtimeave_id_forcoord)
    first_coord_1 = coord_1[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][0, 0]
    last_coord_1 = coord_1[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][-1, 0]
    first_coord_2 = coord_2[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][0, 0]
    last_coord_2 = coord_2[n_y_0: n_y_0+d_n_y, n_z_0: n_z_0+d_n_z][0, -1]
    return (value, first_coord_1, last_coord_1, first_coord_2, last_coord_2)

def velocity(
    mv_v_name,
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    mass_v_name='mass',
    fixtimeave_id_name='avspatial_ave',
    ):
    mass = get_ave_value(n_ave, 'avspatial_ave', mass_v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name)
    mv = get_ave_value(n_ave, 'avspatial_ave', mv_v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name)
    velocity = mv/mass
    return velocity

def velocity_ave_z(
    mv_v_name,
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    mass_v_name='mass',
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ):

    mass = get_ave_value(n_ave, 'avspatial_ave', mass_v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name)
    mv = get_ave_value(n_ave, 'avspatial_ave', mv_v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name)

    # sum over axis 
    mass = np.sum(mass, axis=n_sum_over_axis)
    mv = np.sum(mv, axis=n_sum_over_axis)

    velocity = mv/mass
    coord1 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord1", log_variable_dic_list, fixtimeave_id_name)[:,0]
    
    strain_list = strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    labels_list = [r'$\gamma$' + "={:.3f}".format(strain) for strain in strain_list]
    return (coord1, velocity, labels_list)

def chunk_wall_stress_coord(
    log_variable_dic_list,
    fixtimeave_id_name,
    v_name,
    n_ave,
    inputstepsarray,
    coord_name,
    ):
    force = get_ave_value(
        n_ave, fixtimeave_id_name, v_name, inputstepsarray, log_variable_dic_list,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    Ncount = get_ave_value(
        n_ave, fixtimeave_id_name, 'Ncount', inputstepsarray, log_variable_dic_list,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    force = np.ma.masked_where(Ncount<0.5, force)
    coord = get_coord_by_variable(len(log_variable_dic_list)-1, coord_name, log_variable_dic_list, fixtimeave_id_name)
    d_coord = coord[1] - coord[0]
    grid_area = (
        d_coord
        *float(log_variable_dic_list[-1]['x_period_dp_unit'])
        *float(log_variable_dic_list[-1]['dp'])
    )
    stress = force/grid_area
    return (stress, coord)

def chunk_wall_stress_ratio_coord(
    log_variable_dic_list,
    fixtimeave_id_name,
    v_name_1,
    v_name_2,
    n_ave,
    inputstepsarray,
    coord_name,
    ):
    force_1 = get_ave_value(
        n_ave, fixtimeave_id_name, v_name_1, inputstepsarray, log_variable_dic_list,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    force_2 = get_ave_value(
        n_ave, fixtimeave_id_name, v_name_2, inputstepsarray, log_variable_dic_list,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    coord = get_coord_by_variable(len(log_variable_dic_list)-1, coord_name, log_variable_dic_list, fixtimeave_id_name)
    stress_ratio = force_1/force_2
    Ncount = get_ave_value(
        n_ave, fixtimeave_id_name, 'Ncount', inputstepsarray, log_variable_dic_list,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    stress_ratio = np.ma.masked_where(Ncount<0.5, stress_ratio)
    return (stress_ratio, coord)

def total_wall_force_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    force_v_name,
    fixtimeave_id_name='timeav_inwall_force',
    ):
    force = get_ave_value(n_ave, 'timeav_inwall_force', force_v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name)
    strain = strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    return (strain, force)

def total_wall_force_ratio_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    force_v1_name, force_v2_name,
    fixtimeave_id_name='timeav_inwall_force',
    ):
    force1 = get_ave_value(n_ave, 'timeav_inwall_force', force_v1_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name)
    force2 = get_ave_value(n_ave, 'timeav_inwall_force', force_v2_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name)
    force_ratio = force1/force2
    strain = strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    return (strain, force_ratio)

def quiver_data(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_quiver_vector_1,
    v_name_quiver_vector_2,
    v_name_quiver_axis_1,
    v_name_quiver_axis_2,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name='avspatial_ave',
    is_calculated_v=False,
    ):
    V_1 = get_ave_value(n_ave, fixtimeave_id_fortime, v_name_quiver_vector_1, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_calculated_v=is_calculated_v)
    V_2 = get_ave_value(n_ave, fixtimeave_id_fortime, v_name_quiver_vector_2, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_calculated_v=is_calculated_v)

    coord_1 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord1", log_variable_dic_list, fixtimeave_id_forcoord)
    coord_2 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord2", log_variable_dic_list, fixtimeave_id_forcoord)
    return (coord_1, coord_2, V_1, V_2)

def contour_data(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_quiver_vector,
    v_name_quiver_axis_1,
    v_name_quiver_axis_2,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name='avspatial_ave',
    is_calculated_v=False,
    ):
    V = get_ave_value(n_ave, fixtimeave_id_fortime, v_name_quiver_vector, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_calculated_v=is_calculated_v)

    coord_1 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord1", log_variable_dic_list, fixtimeave_id_forcoord)
    coord_2 = get_coord_by_variable(len(log_variable_dic_list)-1, "Coord2", log_variable_dic_list, fixtimeave_id_forcoord)
    return (coord_1, coord_2, V)