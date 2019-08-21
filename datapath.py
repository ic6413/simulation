# import
import sys
import os
import numpy as np

# === current module inputvariable ===
# set lammps directory (current workspace directory or path)
lammps_directory = os.getcwd() + '/'  #os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')
print('current working directory:' + lammps_directory)

# attribute
attribute_lammps_path = lammps_directory + 'output/setting/attribute.lammps'
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
        try:
            attribute_dict[first_column] = float(second_column)
        except:
            attribute_dict[first_column] = second_column
            
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
    try:
        om_in = attribute_dict['om_in']
    except:
        om_in = 0
    # gravity
    gravitational_acceleration = attribute_dict['gravitational_acceleration']
    gravitation_direction_x = 0
    gravitation_direction_y = 0
    gravitation_direction_z = -1
    g = gravitational_acceleration*np.array([gravitation_direction_x, gravitation_direction_y, gravitation_direction_z])
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
            None,
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [0,0,0],
        ],
    ]
    walls_cy = [
        [
            'cy',
            [0,0,0],
            [0,0,1],
            r_in,
            None,
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [0,0,om_in],
            [0,0,0],
        ],
        [
            'cy',
            [0,0,0],
            [0,0,1],
            r_out,
            None,
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [0,0,0],
        ]
    ]
    walls = walls_p + walls_cy

    walls_name = [None]*len(walls)
    for i, wall in enumerate(walls):
        
        if wall[2] == [0, 0, 1] and wall[0]=='p':
            walls_name[i] = 'z_plane'
        
        if wall[2] == [0, 0, 1] and wall[0]=='cy' and wall[3] == r_in:
            walls_name[i] = 'z_cylinder_in'

        if wall[2] == [0, 0, 1] and wall[0]=='cy' and wall[3] == r_out:
            walls_name[i] = 'z_cylinder_out'

    walls_id_name = [(-1-id, name) for id, name in enumerate(walls_name)]


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
    g = np.asarray(g)
    # parameter
    mu = attributes.mu
    kn = attributes.kn 
    kt = attributes.kt 
    gamma_n = attributes.gamma_n
    gamma_t = attributes.gamma_t
    z_bottom = attributes.z_bottom
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

# Lammps output path of inputvariable files
pair_all_path = lammps_directory + 'output/pair_all/dump.all.pair.allstep'
custom_all_path = lammps_directory + 'output/single_all/dump.all.single.allstep'
custom_near_trace_path = lammps_directory + 'output/single_trace/dump.trace.single.allstep'
thermo_path = lammps_directory + "log.lammps"
def trace_print_path(label, id):
    path = lammps_directory + 'output/trace/' + str(id) + "/" + label
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

