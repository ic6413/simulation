import pickle
import json
import sys
import os
import pandas as pd


def get_lines_fromlog_variable(log_file_path):
    
    if os.path.isfile(log_file_path):
        with open(log_file_path, mode='r') as f:
            lines = f.read().strip().split('\n')
    else:
        sys.exit("log file not exist")

    return lines

def thermo_dataframe(lines):
    n_line_header_thermo = 0
    n_line_end_thermo = 0
    for n, line in enumerate(lines):
        if len(line) >= 2 and line.split()[0]=="Step":
            n_line_header_thermo = n
        if len(line) >= 2 and line.split()[0]=="Loop" and line.split()[1]=="time":           
            n_line_end_thermo = n - 1
    # check if n_line_end_thermo > n_line_header_thermo > 0
    if n_line_end_thermo <= n_line_header_thermo or n_line_header_thermo <= 0:
        sys.exit('thermo output wrong')
    header = lines[n_line_header_thermo].split()
    datastring = lines[n_line_header_thermo + 1: n_line_end_thermo + 1]
    splitted_datastring = [sub.split() for sub in datastring]
    df = pd.DataFrame(
        splitted_datastring,
        columns = header,
        dtype = 'float64'
    )
    return df

def begin_time(thermo_dataframe):
    begin_time = thermo_dataframe['Time'].values[0]
    return begin_time

def end_time(thermo_dataframe):
    end_time = thermo_dataframe['Time'].values[-1]
    return end_time

def begin_to_end_time(thermo_dataframe):
    begin_to_end_time = end_time(thermo_dataframe) - begin_time(thermo_dataframe)
    return begin_to_end_time

def begin_step(thermo_dataframe):
    begin_step = thermo_dataframe['Step'].values[0]
    return begin_step

def end_step(thermo_dataframe):
    end_step = thermo_dataframe['Step'].values[-1]
    return end_step

def begin_to_end_step(thermo_dataframe):
    begin_to_end_step = end_step(thermo_dataframe) - begin_step(thermo_dataframe)
    return begin_to_end_step

def get_a_log_variable_dic_from_a_logfile(log_file_folder_path, logfilename = 'log.lammps'):
    log_file_path = os.path.join(log_file_folder_path, logfilename)
    lines = get_lines_fromlog_variable(log_file_path)
    # get input parameters from lines in log files
    # create log_variable dictionary from log file
    # write values of variables to dictionary 

    lines_start_variable = [line for line in lines if line.startswith("variable")]
    variable_names = [line.split()[1] for line in lines_start_variable]
    log_variable = dict.fromkeys(variable_names)
    
    # get variable
    for variable_name in log_variable.keys():
        satisfy_lines = [line for line in lines_start_variable if line.split()[1] == variable_name]
        
        if len(satisfy_lines) != 0:
            first_line_words = satisfy_lines[0].split()
        
            if first_line_words[2] == "index" and first_line_words[3][0] != '$':
                variable_value = first_line_words[3]

            elif first_line_words[2] == "index" and first_line_words[3][0] == '$':
                second_satisfy_line = satisfy_lines[1]
                second_line_words = second_satisfy_line.split()
                variable_value = second_line_words[3]

            elif first_line_words[2] == "equal" or first_line_words[2] == "string":
                last_satisfy_line = satisfy_lines[-1]
                last_line_words = last_satisfy_line.split()
                variable_value = last_line_words[3]
            
            elif first_line_words[2] == "getenv":
                variable_value = first_line_words[3]

            else:
                pass
            log_variable[variable_name] = variable_value

        else:
            sys.exit("can not find variable {} in log file".format(variable_name))

    # write shearwall parameters to log_variable dictionary
    if "if_inwall_wall_gran" in log_variable.keys():
        if log_variable["if_inwall_wall_gran"]=="yes":
            log_variable["shearwall"] = "yplane"
        else:
            for line in lines:
                if line.startswith("fix") and line.split()[3] == "wall/gran":
                    if line.split()[1] == "inwall": 
                        log_variable["shearwall"] = line.split()[11]
                        break
                    elif line.split()[1] == "y_bottom":
                        log_variable["shearwall"] = line.split()[11]
                        break
                    else:
                        sys.exit("shearwall missed")
    else:
        for line in lines:
            if line.startswith("fix") and line.split()[3] == "wall/gran":
                if line.split()[1] == "inwall": 
                    log_variable["shearwall"] = line.split()[11]
                    break
                elif line.split()[1] == "y_bottom":
                    log_variable["shearwall"] = line.split()[11]
                    break
                else:
                    sys.exit("shearwall missed")

    # write chunk_2_3 parameters to log_variable dictionary
    # select compute line
    lines_start_compute = [line for line in lines if line.startswith("compute")]
    # get chunk 23 info
    satisfy_lines = [line for line in lines_start_compute if line.split()[3] == 'chunk/atom' and line.split()[1] == "chunk_2_3"]
    # get chunk/atom
    if len(satisfy_lines) != 0:
        log_variable["chunk/atom 23"] = [
            satisfy_lines[0].split()[1],
            satisfy_lines[0].split()[5],
            satisfy_lines[0].split()[4],
        ]
    else:
        pass
    
    if "lines" in log_variable.keys():
        sys.exit("lines in log key")
    else:
        log_variable['lines'] = lines

    if "thermo_dataframe" in log_variable.keys():
        sys.exit("thermo_dataframe in log key")
    else:
        log_variable['thermo_dataframe'] = thermo_dataframe(lines)

    if "log_file_path" in log_variable.keys():
        sys.exit("log_file_path in log key")
    else:
        log_variable['log_file_path'] = log_file_path

    if "log_file_folder_path" in log_variable.keys():
        sys.exit("log_file_folder_path in log key")
    else:
        log_variable['log_file_folder_path'] = log_file_folder_path

    return log_variable

def folder_name_under_subfolder_of_data(log_variable_dic_list):
    # L W H Sa rotate_start_time endtime history wallshape xmu hooke/hertz kn kt gamma_n gamma_t timestep
    variable_name = [

    ]
    variable_label_string = {

    }
    name = 0
    return name

def get_all_log_variable(lmp_folder_path):
    # get previous log_variable
    folder_path_list_last_to_initial = [lmp_folder_path]
    log_variable = get_a_log_variable_dic_from_a_logfile(lmp_folder_path)
    log_variable_dic_list_last_to_initial = [log_variable]
    parent_folder_path = os.path.abspath(os.path.join(lmp_folder_path, os.pardir))
    while log_variable['rst_from'] != '0' and os.path.isfile(os.path.join(parent_folder_path, 'log.lammps')):
        log_variable = get_a_log_variable_dic_from_a_logfile(parent_folder_path)
        log_variable_dic_list_last_to_initial.append(
            log_variable
        )
        folder_path_list_last_to_initial.append(
            parent_folder_path
        )
        parent_folder_path = os.path.abspath(os.path.join(parent_folder_path, os.pardir))

    log_variable_dic_list = log_variable_dic_list_last_to_initial[::-1]
    if log_variable_dic_list[0]['rst_from'] != '0':
        sys.exit('log_variable wrong')
    return log_variable_dic_list

def add_total_time(log_variable_dic_list):
    if not int(log_variable_dic_list[0]["rst_from"]) == 0:
        sys.exit('first simu rst_from is not 0')
    total_time = 0
    for log_variable in log_variable_dic_list:
        total_time += begin_to_end_time(log_variable['thermo_dataframe'])
        log_variable['total_time'] = total_time
    return log_variable_dic_list

def add_rotate_start_time(log_variable_dic_list):
    for i in range(len(log_variable_dic_list)):
        if log_variable_dic_list[i]["ifrotate"] == "yes" and float(log_variable_dic_list[i]["Sa"]) != 0:
            break
    if i == 0:
        rotate_start_time = 0
    else:
        rotate_start_time = log_variable_dic_list[i-1]['total_time']
    for log_variable in log_variable_dic_list:
        log_variable['rotate_start_time'] = rotate_start_time
    return log_variable_dic_list

def add_log_variable_from_log_variable_dic_list(log_variable_dic_list):
    log_variable_dic_list = add_total_time(log_variable_dic_list)
    log_variable_dic_list = add_rotate_start_time(log_variable_dic_list)
    return log_variable_dic_list

def intermediate_variables(log_variable_dic_list):
    for log_variable in log_variable_dic_list:
        del log_variable['thermo_dataframe']
        del log_variable['lines']
    return log_variable_dic_list

def dump_variable(lmp_folder_path):
    log_variable_dic_list = get_all_log_variable(lmp_folder_path)
    log_variable_dic_list = add_log_variable_from_log_variable_dic_list(log_variable_dic_list)
    log_variable_dic_list = intermediate_variables(log_variable_dic_list)
    n_simu_total = len(log_variable_dic_list)
    folder_path_list_initial_to_last = [log_variable['log_file_folder_path'] for log_variable in log_variable_dic_list]
    log_variable = log_variable_dic_list[-1]
    
    return [log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last]

def dump_variable_save(lmp_folder_path, outputpicklepath, outputjsonpath, log_variable_name = 'log.lammps'):
    [log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last] = dump_variable(lmp_folder_path)
    # save all global variable get from log_variable to pickle
    if outputpicklepath is not None:
        with open(outputpicklepath, 'wb') as f:
            pickle.dump([log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last], f)
    if outputjsonpath is not None:
        with open(outputjsonpath, 'w') as f:
            json.dump([log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last], f)
    return [log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last]

def get_dumped_variable(lmp_folder_path, inputpicklepath):
    if os.path.isfile(inputpicklepath):
        with open(inputpicklepath, 'rb') as f:
            [log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last] = pickle.load(f)
    else:
        [log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last] = dump_variable_save(lmp_folder_path, inputpicklepath, None)
    
    return [log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last]

