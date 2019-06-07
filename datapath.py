# import
import sys
import os
import numpy as np

# === current module inputvariable ===
# set lammps directory (current workspace directory or path)
lammps_directory = os.getcwd() + '/'  #os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')
print('current working directory:' + lammps_directory)

# attribute
attribute_lammps_path = lammps_directory + 'output/simulation_setting/attribute.lammps'
attribute_py_path = lammps_directory + 'attributes.py'

# ====================================== import attribute
if os.path.isfile(attribute_lammps_path):
    with open(attribute_lammps_path, mode='r') as f:
        lines = f.read().strip().split('\n')
        
        attribute_dict = {}
        for line in lines:
            linesplit = line.split()
            first_column = linesplit[0]
            second_column = linesplit[1]
            attribute_dict[first_column] = float(second_column)
        # startstep
        startstep = attribute_dict['startstep']
        startstep = int(startstep)
        # timestep
        ts = attribute_dict['timestep']
        # atom radius
        dp0 = attribute_dict['diameter']
        density = attribute_dict['density']
        # intersection pointmu
        r_in = attribute_dict['ri']
        r_out = attribute_dict['ro']
        # gravity
        g = attribute_dict['gravitational_acceleration']
        # parameter
        mu = attribute_dict['friction_coefficient']
        kn = attribute_dict['kn'] 
        kt = attribute_dict['kt'] 
        gamma_n = attribute_dict['gamma_n']
        gamma_t = attribute_dict['gamma_t']
        n_type = int(attribute_dict['n_type'])
        if n_type == 1:
            type_radius_list = [
                [1, 1.0*dp0/2]
            ]
        elif n_type == 3:
            dp_big = attribute_dict['dp_big']
            dp_small = attribute_dict['dp_small']
            type_radius_list = [
                [1, dp_small/2],
                [2, dp0/2],
                [3, dp_big/2],
            ]
        else:
            sys.exit("n_type not 3 or 1")
        type_radius_array = np.transpose(np.asarray(type_radius_list))
        z_bottom = attribute_dict['z_bottom']
        walls_p = [
            [
                'p',
                [0,0,z_bottom],
                [0,0,1],
            ],
        ]
        walls_cy = [
            [
                'cy',
                [0,0,0],
                [0,0,1],
                r_in,
            ],
            [
                'cy',
                [0,0,0],
                [0,0,1],
                r_out,
            ]
        ]

elif os.path.isfile(attribute_py_path):
    import imp
    attributes = imp.load_source('attributes', attribute_py_path)
    print ('imported attribute by .py')
    # startstep
    startstep = attributes.startstep
    # timestep
    ts = attributes.ts
    # atom radius
    dp0 = attributes.dp0
    density = attributes.density
    # intersection pointmu
    r_in = attributes.r_in
    r_out = attributes.r_out
    # gravity
    g = attributes.g
    # parameter
    mu = attributes.mu
    kn = attributes.kn 
    kt = attributes.kt 
    gamma_n = attributes.gamma_n
    gamma_t = attributes.gamma_t

    type_radius_array = np.transpose(np.asarray(attributes.type_radius_list))
    walls_p = attributes.walls_p
    walls_cy = attributes.walls_cy

else:
    sys.exit("no attribute .lammps or .py to get attribute")

print ("finish creating attribute")

# ====================================== end import attribute


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
dumpname_custom_all = "all_" + str(startstep) + ".dump"
dumpname_custom = "nb_maxKEatom_" + str(startstep) + ".dump"
dumpname_pair = "pair_trace_" + str(startstep) + ".dump"
logname = "log_" + str(startstep) + ".lammps"
# specify the path of inputvariable files
pair_path = lammps_directory + 'output/dump/' + dumpname_pair
custom_all_path = lammps_directory + 'output/dump/' + dumpname_custom_all
custom_path = lammps_directory + 'output/dump/' + dumpname_custom
thermo_path = lammps_directory + logname
def add_label_max_tracepath(label):
    path = lammps_directory + 'output/trace_print/' + "trace_atom_" + "max" + "_" + label + ".txt"
    return path
def add_label_put_id_tracepath(label, id):
    path = lammps_directory + 'output/trace_print/' + "trace_atom_" + str(id) + "_" + label + ".txt" 
    return path
## folder path of output
post_process_path = lammps_directory + 'postprocess/'
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

