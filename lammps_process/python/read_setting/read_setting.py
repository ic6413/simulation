import sys
import os

# === current module inputvariable ===
# set lammps directory (current workspace directory or path)

lammps_directory = os.getcwd() + '/'
# attribute
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

def read_log(folderpath):
    log_path = folderpath + 'log.lammps'
    if os.path.isfile(log_path):        
        with open(log_path, mode='r') as f:
            lines = f.read().strip().split('\n')
    else:
        sys.exit("file not exist")

    lines_start_variable = [line for line in lines if line.startswith("variable")]
    variable_names_must_be_contained = [line.split()[1] for line in lines if line.startswith("variable")]
    lines_start_compute = [line for line in lines if line.startswith("compute")]

    # read variable from log file
    def read_variable():
        logfile_in_folder_in = {}
        variable_names = [line.split()[1] for line in lines_start_variable]
        
        for variable_name in variable_names:

            satisfy_lines = [line for line in lines_start_variable if line.split()[1] == variable_name]
            
            if len(satisfy_lines) != 0:
                first_line_words = satisfy_lines[0].split()
            
                if first_line_words[2] == "index":
                    variable_value = first_line_words[3]

                elif first_line_words[2] == "equal" or first_line_words[2] == "string":
                    last_satisfy_line = satisfy_lines[-1]
                    last_line_words = last_satisfy_line.split()
                    variable_value = last_line_words[3]

                else:
                    pass

                logfile_in_folder_in[variable_name] = variable_value

            else:
                sys.exit("can not find variable {} in log file".format(variable_name))

        return logfile_in_folder_in


    logfile_in_folder = read_variable()
    if "if_ybottom_wall_gran" in logfile_in_folder.keys():
        if logfile_in_folder["if_ybottom_wall_gran"]=="yes":
            logfile_in_folder["shearwall"] = "yplane"
        else:
            for line in lines:
                if line.startswith("fix") and line.split()[3] == "wall/gran":
                    if line.split()[1] == "inwall": 
                        logfile_in_folder["shearwall"] = line.split()[11]
                        break
                    elif line.split()[1] == "y_bottom":
                        logfile_in_folder["shearwall"] = line.split()[11]
                        break
                    else:
                        sys.exit("shearwall missed")
    else:
        for line in lines:
            if line.startswith("fix") and line.split()[3] == "wall/gran":
                if line.split()[1] == "inwall": 
                    logfile_in_folder["shearwall"] = line.split()[11]
                    break
                elif line.split()[1] == "y_bottom":
                    logfile_in_folder["shearwall"] = line.split()[11]
                    break
                else:
                    sys.exit("shearwall missed")
            
    
    # start to read logfile_in_folder for starting from 0
    def read_log_0():
        logfile_in_folder_0 = {}
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


        dir = folderpath
        step = int(logfile_in_folder["rst_from"])
        while step != 0:
            step = rst_from_in_parent_log(dir)
            dir = parent_dir(dir)
        
        log_0_path = dir + '/log.lammps'

        if os.path.isfile(log_0_path):
                
            with open(log_0_path, mode='r') as f:
                lines_0 = f.read().strip().split('\n')

        else:
            sys.exit("file not exist")
        lines_0_start_compute = [line for line in lines_0 if line.startswith("compute")]
        satisfy_lines = [line for line in lines_0_start_compute if line.split()[3] == 'chunk/atom']
        if len(satisfy_lines) != 0:
            logfile_in_folder_0["chunk/atom"] = [
                                    satisfy_lines[0].split()[1],
                                    satisfy_lines[0].split()[5],
                                    ]
        else:
            pass
        
        return logfile_in_folder_0
    
    for key in read_log_0().keys():
        if key not in logfile_in_folder.keys():
            logfile_in_folder[key] = read_log_0()[key]
    
    return logfile_in_folder

logfile = read_log(lammps_directory)

def log_current_plus_previous(currentfolderpath):
    logfilelist_from_lastest_to_initial = []
    logfilelist_from_lastest_to_initial.append(
        read_log(currentfolderpath)
        )

    folderpath = currentfolderpath

    while int(read_log(folderpath)["rst_from"]) > 0:
        folderpath = folderpath + "../"
        logfilelist_from_lastest_to_initial.append(
        read_log(folderpath)
        )
    return logfilelist_from_lastest_to_initial

n_loglist = len(log_current_plus_previous(lammps_directory))
logfilelist_from_initial_to_lastest = [log_current_plus_previous(lammps_directory)[n_loglist-1-n] for n in range(n_loglist)]

def log_current_plus_previousfrom_initial_to_lastest(currentfolderpath):
    logfilelist_from_lastest_to_initial=log_current_plus_previous(currentfolderpath)
    n_loglist = len(logfilelist_from_lastest_to_initial)
    logfilelist_from_initial_to_lastest = [logfilelist_from_lastest_to_initial[n_loglist-1-n] for n in range(n_loglist)]
    return logfilelist_from_initial_to_lastest