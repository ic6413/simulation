#### function to reading and writing data etc
import sys
import os
import numpy as np
import d00_utils.input_text as di

# determine if Coord1 is y or z
def if_Coord1_is_y_or_z(log_variable):
    if log_variable["shearwall"] == "yplane":
        if "chunk/atom 23" in log_variable.keys():
            if log_variable["chunk/atom 23"][1] == "y":
                Coord1_is_y_or_z = 'y'
            elif log_variable["chunk/atom 23"][1] == "z":
                Coord1_is_y_or_z = 'z'
            else:
                sys.exit("chunk_method wrong")
    return Coord1_is_y_or_z

def len_in_each_dim_coord23(coord_2_3_path, n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2):
    # find n_1
    # python read line
    with open(coord_2_3_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    first_4_lines = lines[0: 4]
    # header
    header_line = first_4_lines[n_row_header-1]
    if header_line[0] == "#":
        n_column_variable_start = 2
    else:
        n_column_variable_start = 1
    header = header_line.split()[n_column_variable_start-1:]
    # stepsarray
    n_line_in_a_step = int(first_4_lines[n_row_header].split()[n_column_of_chunk_number-1])

    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[4: 4+n_line_in_a_step]

    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    if header[1] != "Coord1":
        sys.exit("header[1] is not Coord1")
    array_coord1 = lines[:, 1]
    if array_coord1[0] != array_coord1[1]:
        sys.exit("array_coord1[0] != array_coord1[1]")
    for n, array_value in enumerate(array_coord1):
        if array_value != array_coord1[0]:
            break
    n_2 = n
    if array_coord1.shape[0]%n_2 != 0:
        sys.exit('array number is not divided by n_2')
    n_1 = int(array_coord1.shape[0]/n_2)
    len_in_each_dim_coord23 = {
        "n_1": n_1,
        "n_2": n_2,
    }
    return len_in_each_dim_coord23

# reshape method for coord_2_3 chunk_1D  [timestep, Coord1, Coord2]
def reshape_array(array, output_shape, len_in_each_dim_coord23, n_step=None):
    if len(array.shape) != 1:
        sys.exit('shape of input array is not 1')
    if len(output_shape) == 1:
        pass
    elif len(output_shape) == 2:
        if output_shape == ['n_1', 'n_2']:
            array = array.reshape((len_in_each_dim_coord23['n_1'], len_in_each_dim_coord23['n_2']))
        elif output_shape == ['t', 'n_2']:
            array = array.reshape((n_step, len_in_each_dim_coord23['n_2']))
        elif output_shape == ['t', 'n_1']:
            array = array.reshape((n_step, len_in_each_dim_coord23['n_1']))
        else:
            sys.exit('output shape wrong')
    elif len(output_shape) == 3:
        if output_shape == ['t', 'n_1', 'n_2']:
            array = array.reshape((n_step, len_in_each_dim_coord23['n_1'], len_in_each_dim_coord23['n_2']))
        else:
            sys.exit('output shape wrong')
    else:
        sys.exit('len of output_shape wrong')
    array_with_shape = {
        'array': array,
        'shape': output_shape,
    }
    return array_with_shape

# reorder method for coord_2_3 chunk_1D  [timestep, y, z] or [timestep, y] or [timestep, z]
def reorder_array(array_with_shape, Coord1_is_y_or_z):
    if Coord1_is_y_or_z == 'y':
        reordered_array = array_with_shape['array']
    elif Coord1_is_y_or_z == 'z':
        if array_with_shape['shape'] == ['n_1', 'n_2']:
            reordered_array = np.transpose(
                array_with_shape['array'],
                axes=(1, 0),
            )
        elif array_with_shape['shape'] == ['t', 'n_1', 'n_2']:
            reordered_array = np.transpose(
                array_with_shape['array'],
                axes=(0, 2, 1),
            )
    return reordered_array

# check 3 line in vector, 2 line in scalar
def check_fixtimeave_mode_consistency(mode, n_first_number_row_expected, lines):
    for n in range(n_first_number_row_expected-1):
        if lines[n].split()[0][0].isdigit():
            sys.exit('mode not consistency with number of headers')
    if not lines[n_first_number_row_expected-1].split()[0][0].isdigit():
        sys.exit('mode not consistency with number of headers')

def check_npy_file_exist(filepath):
    if filepath[-4:] == ".npy" and os.path.exists(filepath):
        print("npy file exist, path is " + filepath)
    if filepath[-4:] != ".npy" and os.path.exists(filepath + ".npy"):
        print("npy file exist, path is " + filepath)
# mode vector
def save_coord_to_npy(
    input_text_path, output_folder_path, output_shape, len_in_each_dim_coord23,
    log_variable, fixtimeave_id_name,
    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    fixtype='fixavetime',
    ):
    # python read line
    with open(input_text_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    if fixtype=='fixavetime':
        check_fixtimeave_mode_consistency(
            log_variable[fixtype][fixtimeave_id_name]['mode'],
            n_row_header+1,
            lines,
        )
    first_4_lines = lines[0: n_row_header+1]
    # header
    header_line = first_4_lines[n_row_header-1]
    if header_line[0] == "#":
        n_column_variable_start = 2
    else:
        n_column_variable_start = 1
    header = header = header_line.split()[n_column_variable_start-1:]
    # stepsarray
    n_line_in_a_step = int(first_4_lines[n_row_header].split()[n_column_of_chunk_number-1])

    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[n_row_header+1: n_row_header+1+n_line_in_a_step]

    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    n_header = len(header)
    os.makedirs(output_folder_path, exist_ok=True)

    # save all variable in the folder
    for i in range(n_header):
        if header[i] in di.coord_name_list:
            lines_column = lines[:, i]
            array_with_shape = reshape_array(lines_column, output_shape, len_in_each_dim_coord23, n_step=None)
            Coord1_is_y_or_z = if_Coord1_is_y_or_z(log_variable)
            ordered_array = reorder_array(array_with_shape, Coord1_is_y_or_z)
            outputfilepath = os.path.join(
                output_folder_path,
                header[i],
            )
            check_npy_file_exist(outputfilepath)
            np.save(outputfilepath, ordered_array)
# mode vector
def save_non_coord_to_npy(
    input_text_path, output_folder_path, output_shape, len_in_each_dim_coord23,
    log_variable, fixtimeave_id_name,
    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    fixtype='fixavetime',
    ):
    # python read line
    with open(input_text_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    if fixtype=='fixavetime':
        check_fixtimeave_mode_consistency(
            log_variable[fixtype][fixtimeave_id_name]['mode'],
            n_row_header+1,
            lines,
        )
    n_line_in_file = len(lines)
    first_4_lines = lines[0: n_row_header + 1]
    # header
    header_line = first_4_lines[n_row_header-1]
    if header_line[0] == "#":
        n_column_variable_start = 2
    else:
        n_column_variable_start = 1
    header = header_line.split()[n_column_variable_start-1:]
    # stepsarray
    n_line_in_a_step = int(first_4_lines[n_row_header].split()[1])
    step_first_in_file = int(first_4_lines[n_row_header].split()[n_column_of_step - 1])
    step_second_in_file = int(lines[n_row_header + n_line_in_a_step + 1].split()[n_column_of_step - 1])
    step_last_in_file = int(lines[-1 - n_line_in_a_step].split()[n_column_of_step - 1])
    d_step = step_second_in_file - step_first_in_file
    totalstepsarray = np.arange(step_first_in_file, step_last_in_file + d_step, d_step, dtype=int)
    n_step = len(totalstepsarray)
    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[n_row_header: ]
    for i in reversed(range(len(totalstepsarray))):
        # del from largest
        del lines[(n_line_in_a_step + 1)*i]
    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    n_line = lines.shape[0]
    if n_line % n_step == 0:
        n_line_one_step = int(n_line/n_step)
    else:
        sys.exit("n_step not divisible n_line")
    n_header = len(header)
    os.makedirs(output_folder_path, exist_ok=True)

    # if Ncount in header
    if 'Ncount' in header:
        # reset variable if Ncount*freq_repeat is smaller than 1
        index_of_Ncount_in_header = header.index('Ncount')
        Ncount_column = lines[:, index_of_Ncount_in_header]
        reset_mask = (
            Ncount_column*int(log_variable[fixtype][fixtimeave_id_name]['n_repeat']) < 0.99
        )
    # save all variable in the folder
    for i in range(n_header):
        lines_column = lines[:, i]
        # reset value to zero if Ncount*freq_repeat is smaller than 1
        if 'Ncount' in header:
            lines_column[reset_mask] = 0
        array_with_shape = reshape_array(lines_column, output_shape, len_in_each_dim_coord23, n_step=n_step)
        Coord1_is_y_or_z = if_Coord1_is_y_or_z(log_variable)
        ordered_array = reorder_array(array_with_shape, Coord1_is_y_or_z)
        
        outputfilepath = os.path.join(
            output_folder_path,
            header[i],
        )
        check_npy_file_exist(outputfilepath)
        np.save(
            outputfilepath,
            ordered_array,
        )
    # save timestep in the folder
    timestep_filepath = os.path.join(
        output_folder_path,
        "timestep",
    )
    check_npy_file_exist(timestep_filepath)
    np.save(timestep_filepath, totalstepsarray)
# mode scalar
def save_non_coord_to_npy_scalar(
    input_text_path, output_folder_path,
    log_variable, fixtimeave_id_name,
    n_row_header=2,
    ):
    # python read line
    with open(input_text_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    check_fixtimeave_mode_consistency(
        log_variable['fixavetime'][fixtimeave_id_name]['mode'],
        n_row_header+1,
        lines,
    )

    # header
    header_line = lines[n_row_header-1]
    if header_line[0] == "#":
        n_column_variable_start = 2
    else:
        n_column_variable_start = 1
    header = header_line.split()[n_column_variable_start-1:]

    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[n_row_header: ]
    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    n_header = len(header)
    os.makedirs(output_folder_path, exist_ok=True)

    # save all variable in the folder
    for i in range(n_header):
        lines_column = lines[:, i]
        
        outputfilepath = os.path.join(
            output_folder_path,
            header[i],
        )
        check_npy_file_exist(outputfilepath)
        np.save(
            outputfilepath,
            lines_column,
        )

def multi_simu_save_non_coord_to_npy_scalar(
    folder_path_list_initial_to_last,
    log_variable_dic_list, fixtimeave_id_name,
    n_row_header=2,
    ):
    for n, log_variable in enumerate(log_variable_dic_list):
        # mode scalar
        if fixtimeave_id_name in log_variable_dic_list[n]['fixavetime'].keys():
            output_folder_path = di.fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list)
            input_text_path = di.fixtimeave_text_file_path(n, fixtimeave_id_name, log_variable_dic_list, folder_path_list_initial_to_last)
            if os.path.exists(input_text_path):
                save_non_coord_to_npy_scalar(
                    input_text_path, output_folder_path,
                    log_variable, fixtimeave_id_name,
                    n_row_header=2,
                )

def multi_simu_save_coord_23_to_npy_replace(
    folder_path_list_initial_to_last,
    len_in_each_dim_coord23,
    log_variable_dic_list,
    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    ):
    for n, log_variable in enumerate(log_variable_dic_list):
        # mode scalar
        if di.coord_chunk_id_23_replace in log_variable_dic_list[n]['fixavechunk'].keys():
            output_folder_path = di.fixtimeave_npy_output_folder_path(n, di.coord_chunk_id_23, log_variable_dic_list)
            input_text_path = di.fixchunkave_text_file_path(n, di.coord_chunk_id_23_replace, log_variable_dic_list, folder_path_list_initial_to_last)
            if os.path.exists(input_text_path):
                save_coord_to_npy(
                    input_text_path, output_folder_path, di.output_shape_map_from_id[di.coord_chunk_id_23], len_in_each_dim_coord23,
                    log_variable, di.coord_chunk_id_23_replace,
                    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
                    fixtype="fixavechunk",
                )

def inwall_cood_from_chunk23_coord(coord_array_from_chunk23):
    if len(coord_array_from_chunk23.shape) != 2:
        sys.exit("len(coord_array_from_chunk23.shape) not 2")
    coord_inwall = coord_array_from_chunk23[0,:]
    return coord_inwall

def outwall_cood_from_chunk23_coord(coord_array_from_chunk23):
    if len(coord_array_from_chunk23.shape) != 2:
        sys.exit("len(coord_array_from_chunk23.shape) not 2")
    coord_outwall = coord_array_from_chunk23[-1,:]
    return coord_outwall

def zbottom_cood_from_chunk23_coord(coord_array_from_chunk23):
    if len(coord_array_from_chunk23.shape) != 2:
        sys.exit("len(coord_array_from_chunk23.shape) not 2")
    coord_zbottom = coord_array_from_chunk23[:,0]
    return coord_zbottom

def save_inwall_coord_from_chunk23(log_variable_dic_list):
    for name in di.coord_name_list:
        for n in range(len(log_variable_dic_list)):
            value = np.load(di.fixtimeave_npy_output_file_path(n, di.coord_chunk_id_23, log_variable_dic_list, name))
            os.makedirs(di.fixtimeave_npy_output_folder_path(n, di.coord_chunk_id_inwall, log_variable_dic_list), exist_ok=True)
            np.save(
                di.fixtimeave_npy_output_file_path(n, di.coord_chunk_id_inwall,
                log_variable_dic_list,name), inwall_cood_from_chunk23_coord(value),
            )

def save_outwall_coord_from_chunk23(log_variable_dic_list):
    for name in di.coord_name_list:
        for n in range(len(log_variable_dic_list)):
            value = np.load(di.fixtimeave_npy_output_file_path(n, di.coord_chunk_id_23, log_variable_dic_list, name))
            os.makedirs(di.fixtimeave_npy_output_folder_path(n, di.coord_chunk_id_outwall, log_variable_dic_list), exist_ok=True)
            np.save(
                di.fixtimeave_npy_output_file_path(n, di.coord_chunk_id_outwall,
                log_variable_dic_list,name), outwall_cood_from_chunk23_coord(value),
            )

def save_zbottom_coord_from_chunk23(log_variable_dic_list):
    for name in di.coord_name_list:
        for n in range(len(log_variable_dic_list)):
            value = np.load(di.fixtimeave_npy_output_file_path(n, di.coord_chunk_id_23, log_variable_dic_list, name))
            os.makedirs(di.fixtimeave_npy_output_folder_path(n, di.coord_chunk_id_zbottom, log_variable_dic_list), exist_ok=True)
            np.save(
                di.fixtimeave_npy_output_file_path(n, di.coord_chunk_id_zbottom,
                log_variable_dic_list,name), zbottom_cood_from_chunk23_coord(value),
            )

def multi_simu_save_non_coord_to_npy(
    folder_path_list_initial_to_last,
    len_in_each_dim_coord23,
    log_variable_dic_list, fixtimeave_id_name,
    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    ):
    for n, log_variable in enumerate(log_variable_dic_list):
        # mode scalar
        if fixtimeave_id_name in log_variable_dic_list[n]['fixavetime'].keys():
            output_folder_path = di.fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list)
            input_text_path = di.fixtimeave_text_file_path(n, fixtimeave_id_name, log_variable_dic_list, folder_path_list_initial_to_last)
            if os.path.exists(input_text_path):
                save_non_coord_to_npy(
                    input_text_path, output_folder_path, di.output_shape_map_from_id[fixtimeave_id_name], len_in_each_dim_coord23,
                    log_variable, fixtimeave_id_name,
                    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
                )
        elif fixtimeave_id_name in di.map_fixtimeave_to_fixchunkave:
            if di.map_fixtimeave_to_fixchunkave[fixtimeave_id_name] in log_variable_dic_list[n]['fixavechunk'].keys():
                output_folder_path = di.fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list)
                input_text_path_chunk = di.fixchunkave_text_file_path(n, di.map_fixtimeave_to_fixchunkave[fixtimeave_id_name], log_variable_dic_list, folder_path_list_initial_to_last)
                
                if os.path.exists(input_text_path_chunk):
                    save_non_coord_to_npy(
                        input_text_path_chunk, output_folder_path, di.output_shape_map_from_id[fixtimeave_id_name], len_in_each_dim_coord23,
                        log_variable, di.map_fixtimeave_to_fixchunkave[fixtimeave_id_name],
                        n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
                        fixtype='fixavechunk',
                    )

def multi_simu_save_coord_to_npy(
    folder_path_list_initial_to_last,
    len_in_each_dim_coord23,
    log_variable_dic_list, fixtimeave_id_name,
    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
    ):
    for n, log_variable in enumerate(log_variable_dic_list):
        # mode scalar
        if fixtimeave_id_name in log_variable_dic_list[n]['fixavetime'].keys():
            output_folder_path = di.fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list)
            input_text_path = di.fixtimeave_text_file_path(n, fixtimeave_id_name, log_variable_dic_list, folder_path_list_initial_to_last)
            if os.path.exists(input_text_path):
                save_coord_to_npy(
                    input_text_path, output_folder_path, di.output_shape_map_from_id[fixtimeave_id_name], len_in_each_dim_coord23,
                    log_variable, fixtimeave_id_name,
                    n_row_header=3, n_column_of_step=1, n_column_of_chunk_number=2,
                )

def rename_all_variable(n_simu_total, log_variable_dic_list):
    for fixtimeave_id_name in di.outputlist:
        for n in range(n_simu_total):
            folder_path = di.fixtimeave_npy_output_folder_path(n, fixtimeave_id_name, log_variable_dic_list)
            for root, dirs, files in os.walk(folder_path):
                for name in files:
                    if di.eliminate_npy_if_yes(name) in di.dic_rename_fixtimeave_npy_headername_to_use:
                        os.rename(
                            os.path.join(root, name),
                            os.path.join(root, di.add_npy_if_not(di.dic_rename_fixtimeave_npy_headername_to_use[di.eliminate_npy_if_yes(name)])),
                        )

def get_timestep(n, fixid_for_timestep, log_variable_dic_list):
    timestepfilepath = di.fixtimeave_npy_output_file_path(n, fixid_for_timestep, log_variable_dic_list, 'timestep')
    timestep = np.load(timestepfilepath, mmap_mode='r')
    return timestep
    
def get_variable(n, v_name, log_variable_dic_list, fixtimeave_id_name=None, is_std=False, is_calculated_v=False):
    v_path = di.v_name_to_path(n, v_name, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_std=is_std, is_calculated_v=is_calculated_v)
    value = np.load(v_path, mmap_mode='r')
    return value

def get_coord_by_variable(n, coord_name, log_variable_dic_list, fixtimeave_id_name):
    v_path = di.fixtimeave_npy_output_coord_file_path(n, fixtimeave_id_name, log_variable_dic_list, coord_name)
    value = np.load(v_path, mmap_mode='r')
    return value                    