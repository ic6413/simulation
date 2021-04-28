##### set all global variable ######
import pickle
import sys
import os

# executive path 
lammps_directory = os.getcwd() + '/'
# read from pickle if exist
if os.path.isfile(lammps_directory + 'python_global.pckl'):
    f = open(lammps_directory + 'python_global.pckl', 'rb')
    [log_variable_dic_list, n_simu_total, line_log_list, log_variable, folder_path_list_initial_to_last] = pickle.load(f)
    f.close()
# if global pickle not exist
else: 
    # read line_log_list, all log_variable for all current and previous simulation  
    # and count n_simu_total, the number of simulation including current and previous 
    n_simu_total = 0
    line_log_list = []
    log_variable_dic_list = []
    
    def parent_dir(mypath):
        return os.path.abspath(os.path.join(mypath, os.pardir))

    def step_from_current_dir(dir):
        if os.path.isfile(dir):
            with open(dir, mode='r') as f:
                lines = f.read().strip().split('\n')
        else:
            sys.exit("file not exist")
        lines_start_variable = [line for line in lines if line.startswith("variable")]
        satisfy_lines = [line for line in lines_start_variable if line.split()[1] == 'rst_from']
        step = int(satisfy_lines[0].split()[3])
        return step

    def rst_from_in_parent_log(mypath):
        parent_log_path = parent_dir(mypath) + '/log.lammps'
        return step_from_current_dir(parent_log_path)

    # initialize dir as lammps_directory
    dir = lammps_directory
    while True:
        # read log file
        if os.path.isfile(dir + '/log.lammps'):
            with open(dir + '/log.lammps', mode='r') as f:
                lines = f.read().strip().split('\n')
        else:
            sys.exit("log file not exist")
        # append lines to line_log_list and count n_simu_total
        line_log_list.append(lines)
        n_simu_total += 1
        # find rst_from
        # read variable from log file
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

        # shearwall
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
        # append log_variable to log_variable_dic_list
        log_variable_dic_list.append(log_variable)
        # if rst_from = 0 stop reading otherwise continue to set dir to parent folder
        rst_from = int(
            [line for line in lines_start_variable if line.split()[1] == 'rst_from'][0].split()[3]
        )
        if rst_from == 0:
            break
        else:
            dir = os.path.abspath(os.path.join(dir, os.pardir))
    log_variable_dic_list = log_variable_dic_list[::-1]
    line_log_list = line_log_list[::-1]
    folder_path_list_initial_to_last = [lammps_directory+"../"*(n_simu_total-1-n) for n in range(n_simu_total)]
    folder_path_list_last_to_initial = folder_path_list_initial_to_last[::-1]
    # variable names in lammps log file should be contained 
    """
    variable_names_must_be_contained = [
        "rst_from",
        "runstep",
        "runstep_loop",
        "n_loop",
        "KEtr_to_jumpout",
        "contact_change",
        "contact_model",
        "n_processor_x",
        "ts",
        "eliminate_history",
        "ifrotate",
        "trace_id1",
        "trace_id2",
        "trace_id3",
        "iftrace_maxKEt",
        "iftrace_maxKEr",
        "iftrace_maxKEtr",
        "ifappend",
        "iffreeze",
        "create_atom_method",
        "ifpour",
        "pourstep",
        "pour_number",
        "ifdeposit",
        "deposit_number",
        "chunk_method",
        "ifsetvelocity",
        "ifresetdt",
        "recount_dof",
        "freq_restart_big",
        "freq_restart_small",
        "freq_dump_trace_image",
        "freq_dump_all_image",
        "freq_dump_single_trace",
        "freq_dump_single_all",
        "freq_dump_pair_trace",
        "freq_dump_pair_all",
        "freq_thermo",
        "freq_print_trace",
        "freq_fixvector",
        "every_fixvector",
        "freq_ave_chunk_momentum_mass_field",
        "repeat_ave_chunk_momentum_mass_field",
        "every_ave_chunk_momentum_mass_field",
        "freq_ave_wall",
        "repeat_ave_wall",
        "every_ave_wall",
        "freq_balance",
        "cutoff_dumpnb",
        "n_type",
        "ifairviscous",
        "gamma_air",
        "g",
        "dp",
        "rp",
        "den",
        "mp",
        "dp_big_dp_unit",
        "dp_small_dp_unit",
        "skin_dp_unit",
        "ri_wall_dp_unit",
        "width_wall_dp_unit",
        "zlo_wall_dp_unit",
        "ro_wall_dp_unit",
        "ri_create_dp_unit",
        "ro_create_dp_unit",
        "zlo_create_dp_unit",
        "zhi_create_dp_unit",
        "zhi_box_dp_unit",
        "zhi_chunk_dp_unit",
        "zlo_pour_dp_unit",
        "zhi_pour_dp_unit",
        "zhi_freeze_dp_unit",
        "ri_freeze_dp_unit",
        "ro_freeze_dp_unit",
        "kn_scale",
        "kn",
        "kt",
        "gamma_n_scale",
        "gamma_n",
        "gamma_t",
        "xmu",
        "ifdamp_tangent",
        "N_bin_x",
        "N_bin_y",
        "N_bin_z",
        "N_bin_r",
        "N_bin_theta",
        "Sa",
        "ro_wall",
        "ri_wall",
        "zlo_wall",
        "zhi_wall",
        "zhi_chunk",
        "ro_box",
        "zhi_box",
        "dp_big",
        "dp_small",
        "skin",
        "ri_create",
        "ro_create",
        "zlo_create",
        "zhi_create",
        "zlo_pour",
        "zhi_pour",
        "N_total_bin",
        "N_total_rz",
        "omega_in",
        "in_velocity",
        "zhi_freeze",
        "ri_freeze",
        "ro_freeze",
        "edge_lattice",
        "z_fraction_of_lattice",
    ]
    variable_names_must_be_contained_no_need_temporarily = [
        "ro_box_dp_unit",
        "freq_dump_all_vtk",
        "freq_dump_all_stl",
        "freq_dump_all_movie",
        "ratio_type1",
        "ratio_type2",
        "ratio_type3",
    ]

    variable_names_must_be_contained_add_periodic_block = [
        "x_period_dp_unit",
    ]
    """
    # calculate previous time for every log
    restart_time = 0
    for n in range(n_simu_total):
        if n == 0:
            log_variable_dic_list[n]["previous_time"] = restart_time
        else:
            for n_line1, line in enumerate(line_log_list[n-1]):
                if len(line) >= 2:
                    if line.split()[0]=="Step" and line.split()[1]=="Time":
                        break
            
            for n_line2, line in enumerate(line_log_list[n-1]):
                if len(line) >= 2:
                    if line.split()[0]=="Loop" and line.split()[1]=="time":           
                        break

            for n_line, line in enumerate(line_log_list[n-1]):
                if n_line>n_line1 and n_line<n_line2:
                    if line.split()[0]==str(log_variable_dic_list[n]['rst_from']):
                        restart_time += float(line.split()[1])
            log_variable_dic_list[n]["previous_time"] = restart_time
    
    # count the time no rotation
    rotate_start_time = 0
    for i in range(n_simu_total):
        if int(log_variable_dic_list[i]["rst_from"]) == 0:
            previous_time_step = 0
        else:
            previous_time_step = float(log_variable_dic_list[i-1]["ts"])
        rotate_start_time += previous_time_step*int(log_variable_dic_list[i]["rst_from"])
        if log_variable_dic_list[i]["ifrotate"] == "yes" and float(log_variable_dic_list[i]["Sa"]) != 0:
            log_variable_dic_list[i]["rotate_start_time"] = rotate_start_time
            break
        else:
            log_variable_dic_list[i]["rotate_start_time"] = None # have not rotate

    for i in range(n_simu_total):    
        log_variable_dic_list[i]["rotate_start_time"] = rotate_start_time
        
    # log_variable
    log_variable = log_variable_dic_list[-1]
    # save all global variable get from log_variable to pickle
    f = open(lammps_directory + 'python_global.pckl', 'wb')
    pickle.dump([log_variable_dic_list, n_simu_total, line_log_list, log_variable, folder_path_list_initial_to_last], f)
    f.close()