# import
import os

# import simulation module
import inputvariable

# inputvariable
lammps_directory = inputvariable.lammps_directory
# end inputvariable

# specify the path of inputvariable files
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
