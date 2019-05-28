# import
import os

# === current module inputvariable ===
# set lammps directory (current workspace directory or path)
lammps_directory = os.getcwd() + '/'  #os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')
print('current working directory:' + lammps_directory)

# Lammps start step
stepstart = 0

def put_id_on_file(id_list, f_name_without_id):
    
    if id_list != 'all':
        f_name_add_id = (f_name_without_id[:-3] + '_id_' + '_'.join(str(i) for i in id_list))
    else:
        f_name_add_id = f_name_without_id

    return f_name_add_id

# Lammps output file name
dumpname_custom = "nb_maxKEatom_0.dump"
dumpname_pair = "pair_trace_0.dump"
logname = "log_" + str(stepstart) + ".lammps"
# specify the path of inputvariable files
pair_path = lammps_directory + 'output/dump/' + dumpname_pair
custom_path = lammps_directory + 'output/dump/' + dumpname_custom
thermo_path = lammps_directory + logname
## folder path of output
post_process_path = lammps_directory + 'postprocess/'
# attribute
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
