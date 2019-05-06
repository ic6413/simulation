# import
import os

# input directory
lammps_directory = os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')

# specify the path of input files
pair_path = lammps_directory + 'output/dump/pair_trace_0.dump'
custom_path = lammps_directory + 'output/dump/nb_maxKEatom_0.dump'
thermo_path = lammps_directory + 'log_0.lammps'
## folder path of output
post_process_path = lammps_directory + 'postprocess/'
# attribute
attribute_json_path = post_process_path + 'attribute/'
# hdf5
hdf5_csv_path = post_process_path + 'hdf5_csv/'
# debug
debug_print_path = post_process_path + 'debug/'
# interactive
interactive_path = post_process_path + 'interactive/'
# diagram
diagram_path = post_process_path + 'diagram/'


# file path of output
f_attribute = attribute_json_path + 'attribute'

# output filename, no filename extension
#=========== data arrange ========================

# select combine id_i
combine_id_i_list = 'all'
custom_id_i_list = 'all' # id_list or 'all'

def put_id_on_file(id_list, f_name_without_id):
    
    if custom_id_i_list != 'all':
        f_name_add_id = (f_name_without_id + '_' + '_'.join(str(i) for i in custom_id_i_list))
    else:
        f_name_add_id = (f_name_without_id + '_' + 'all')

    return f_name_add_id 

f_cipcj = hdf5_csv_path + put_id_on_file(combine_id_i_list, 'cipcj_id')
f_custom = hdf5_csv_path + put_id_on_file(custom_id_i_list, 'custom_id')
f_thermo = hdf5_csv_path + 'thermo'