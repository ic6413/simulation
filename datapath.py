# import
import os
import initialize

# === current module inputvariable ===
# set lammps directory (current workspace directory or path)
lammps_directory = os.getcwd() + '/'  #os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')
print('current working directory:' + lammps_directory)

# Lammps start step
stepstart = initialize.stepstart

def put_id_on_file(id_list, f_name_without_id):
    
    if id_list != 'all':
        f_name_add_id = (f_name_without_id[:-3] + '_id_' + '_'.join(str(i) for i in id_list)) + ".h5"
    else:
        f_name_add_id = f_name_without_id

    return f_name_add_id

def put_maxlabel_on_file(maxlabel, f_name_without_id):
    f_name_add_id = f_name_without_id[:-3] + '_maxlabel_' + maxlabel + ".h5"
    return f_name_add_id

# Lammps output file name
dumpname_custom = "nb_maxKEatom_" + str(stepstart) + ".dump"
dumpname_pair = "pair_trace_" + str(stepstart) + ".dump"
logname = "log_" + str(stepstart) + ".lammps"
# specify the path of inputvariable files
pair_path = lammps_directory + 'output/dump/' + dumpname_pair
custom_path = lammps_directory + 'output/dump/' + dumpname_custom
thermo_path = lammps_directory + logname
## folder path of output
post_process_path = lammps_directory + 'postprocess/'
# attribute
attribute_py_path = lammps_directory + 'attributes.py'
attribute_json_path = post_process_path + 'attribute/'
# hdf5
hdf5_csv_path = post_process_path + 'hdf5_csv/'
f_thermo = hdf5_csv_path + "thermo.h5"
f_custom = hdf5_csv_path + "custom.h5"
f_pair = hdf5_csv_path + "pair.h5"
f_cipcj = hdf5_csv_path + "cipcj.h5"

# debug
debug_print_path = post_process_path + 'debug/'
debug_fig_path = debug_print_path + 'fig/'
debug_fig_thermo_path = debug_fig_path + 'thermo/'
debug_fig_oneatom_path = debug_fig_path + 'oneatom/'
debug_fig_atomij_path = debug_fig_path + 'atomij/'

# interactive
interactive_path = post_process_path + 'interactive/'
# diagram
diagram_path = post_process_path + 'diagram/'
# attribute path
f_attribute = attribute_json_path + 'attribute'
