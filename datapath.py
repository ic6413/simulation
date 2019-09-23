# import
import sys
import os
import numpy as np
import osmanage as om
# import read_setting
import read_setting.read_setting as rr

# === current module inputvariable ===
# set lammps directory (current workspace directory or path)
lammps_directory = os.getcwd() + '/'  #os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')
print('current working directory:' + lammps_directory)

# attribute
attribute_lammps_path = lammps_directory + 'output/setting/attribute.lammps'
attribute_py_path = lammps_directory + 'attributes.py'

# ====================================== import attribute

# startstep
startstep = int(rr.logfile['rst_from'])
# timestep
ts = eval(rr.logfile['ts'])
# atom radius
dp0 = eval(rr.logfile['dp'])
density = eval(rr.logfile['den'])
# intersection pointmu
breakpoint()
if rr.logfile["shearwall"] == "zcylinder":
    r_in = eval(rr.logfile['ri_wall'])
    r_out = eval(rr.logfile['ro_wall'])

try:
    omega_in = eval(rr.logfile['omega_in'])
except:
    omega_in = 0
try:
    N_bin_r = int(rr.logfile['N_bin_r'])
except:
    pass
try:
    N_bin_z = int(rr.logfile['N_bin_z'])
except:
    pass
# gravity
gravitation_direction_x = 0
gravitation_direction_y = 0
gravitation_direction_z = -1
g = eval(rr.logfile['g'])*np.array([gravitation_direction_x, gravitation_direction_y, gravitation_direction_z])
# parameter
mu = eval(rr.logfile['xmu'])
kn = eval(rr.logfile['kn']) 
kt = eval(rr.logfile['kt'])
gamma_n = eval(rr.logfile['gamma_n'])
gamma_t = eval(rr.logfile['gamma_t'])
n_type = int(rr.logfile['n_type'])
if n_type == 1:
    type_radius_list = [
        [1, 1.0*dp0/2]
    ]
elif n_type == 3:
    dp_big = eval(rr.logfile['dp_big'])
    dp_small = eval(rr.logfile['dp_small'])
    type_radius_list = [
        [1, dp_small/2],
        [2, dp0/2],
        [3, dp_big/2],
    ]
else:
    sys.exit("n_type not 3 or 1")
type_radius_array = np.transpose(np.asarray(type_radius_list))
zlo_wall = eval(rr.logfile['zlo_wall'])
walls_p = [
    [
        'p',
        [0,0,zlo_wall],
        [0,0,1],
        None,
        [0,0,0],
        [0,0,0],
        [0,0,0],
        [0,0,0],
        [0,0,0],
    ],
]
if rr.logfile["shearwall"] == "zcylinder":
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
            [0,0,omega_in],
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

filelist = os.listdir(lammps_directory + 'output/single_all/')
try:
    for i in range(2):
        file = filelist[i]
        if not file.endswith("allstep"):
            custom_all_firststep_path = lammps_directory + 'output/single_all/' + file
        else:
            pass
except:
    pass

custom_near_trace_path = lammps_directory + 'output/single_trace/dump.trace.single.allstep'
custom_near_trace_firststep_path = lammps_directory + 'output/single_trace/dump.trace.single.' + str(startstep)
thermo_path = lammps_directory + "log.lammps"
def trace_print_path(label, id):
    path = lammps_directory + 'output/trace/' + str(id) + "/" + label
    return path
## folder path of output
post_process_path = lammps_directory + 'postprocess/'
om.create_directory(post_process_path)
# hdf5
hdf5_csv_path = post_process_path + 'hdf5_csv/'
f_thermo = hdf5_csv_path + "thermo.h5"
f_custom = hdf5_csv_path + "custom.h5"
f_pair = hdf5_csv_path + "pair.h5"
f_cipcj = hdf5_csv_path + "cipcj.h5"

# debug
debug_print_path = post_process_path + 'debug/'
om.create_directory(debug_print_path)
debug_fig_path = debug_print_path + 'fig/'
debug_fig_thermo_path = debug_fig_path + 'thermo/'
debug_fig_oneatom_path = debug_fig_path + 'oneatom/'
debug_fig_atomij_path = debug_fig_path + 'atomij/'

# interactive
interactive_path = post_process_path + 'interactive/'
# diagram
diagram_path = post_process_path + 'diagram/'
om.create_directory(diagram_path)
f_momentum_mass_field_path = diagram_path + "momentum_mass_field/"
om.create_directory(f_momentum_mass_field_path)
f_momentum_mass_field_rz_path = f_momentum_mass_field_path + "rz/"
om.create_directory(f_momentum_mass_field_rz_path)
f_momentum_mass_field_rtheta_path = f_momentum_mass_field_path + "rtheta/"
om.create_directory(f_momentum_mass_field_rtheta_path)
f_momentum_mass_field_density_path = f_momentum_mass_field_path + "density/"
om.create_directory(f_momentum_mass_field_density_path)