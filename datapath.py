# import
import os
# import project module
import output_control as oc

# input directory
lammps_directory = oc.folder_path

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
f_cipcj = hdf5_csv_path + oc.f_cipcj
f_custom = hdf5_csv_path + oc.f_custom
f_thermo = hdf5_csv_path + 'thermo'