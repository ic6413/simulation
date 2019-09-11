import sys
import os

# === current module inputvariable ===
# set lammps directory (current workspace directory or path)
lammps_directory = os.getcwd() + '/'
# attribute

log_path = lammps_directory + 'log.lammps'

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
    "ifstoreforce",
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
    "ri_dp_unit",
    "width_dp_unit",
    "z_bottom_dp_unit",
    "ro_dp_unit",
    "ri_delete_dp_unit",
    "ro_insert_dp_unit",
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
    "ro",
    "ri",
    "z_bottom",
    "z_top",
    "zhi_chunk",
    "r_box",
    "zhi_box",
    "dp_big",
    "dp_small",
    "skin",
    "ri_delete",
    "ro_insert",
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
    "r_box_dp_unit",
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
    logfile = {}
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

            logfile[variable_name] = variable_value

        else:
            sys.exit("can not find variable {} in log file".format(variable_name))

    return logfile


logfile = read_variable()


satisfy_lines = [line for line in lines_start_compute if line.split()[3] == 'chunk/atom']
if len(satisfy_lines) != 0:
    logfile["chunk/atom"] = [
                             satisfy_lines[0].split()[1],
                             satisfy_lines[0].split()[5],
                            ]
else:
    pass


for line in lines:
    if line.startswith("fix") and line.split()[3] == "wall/gran":
        if line.split()[1] == "inwall": 
            logfile["shearwall"] = line.split()[11]
            break
        elif line.split()[1] == "y_bottom":
            logfile["shearwall"] = line.split()[11]
            break
        else:
            sys.exit("shearwall missed")

