# import
import sys
import os
import numpy as np
import os
# import read_setting
import read_setting as rr
# setting
abs_error_tolerence = 1e-13
if_plot_velocity_field_scale_same = "yes"
quiver_scale_velocity_xaxis_shearplanenormal_yaxis_z = 0.1
quiver_scale_velocity_xaxis_shearplaneshear_yaxis_z = 1

# ====================================== import attribute

# startstep
startstep = int(rr.logfile['rst_from'])
# timestep
ts = eval(rr.logfile['ts'])
# atom radius
dp0 = eval(rr.logfile['dp'])
density = eval(rr.logfile['den'])
# intersection pointmu
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
#kn = eval(rr.logfile['kn']) 
#kt = eval(rr.logfile['kt'])
#gamma_n = eval(rr.logfile['gamma_n'])
#gamma_t = eval(rr.logfile['gamma_t'])
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
pair_all_path = rr.lammps_directory + 'output/pair_all/dump.all.pair.allstep'
custom_all_path = rr.lammps_directory + 'output/single_all/dump.all.single.allstep'

filelist = os.listdir(rr.lammps_directory + 'output/single_all/')
try:
    for i in range(2):
        file = filelist[i]
        if not file.endswith("allstep"):
            custom_all_firststep_path = rr.lammps_directory + 'output/single_all/' + file
        else:
            pass
except:
    pass

custom_near_trace_path = rr.lammps_directory + 'output/single_trace/dump.trace.single.allstep'
custom_near_trace_firststep_path = rr.lammps_directory + 'output/single_trace/dump.trace.single.' + str(startstep)
thermo_path = rr.lammps_directory + "log.lammps"
def trace_print_path(label, id):
    path = rr.lammps_directory + 'output/trace/' + str(id) + "/" + label
    return path



## folder path of output

# postprocess
post_process_path = rr.lammps_directory + 'postprocess/'
os.makedirs(post_process_path, exist_ok=True)

# processed data
post_data_path = post_process_path + "post_data/"
os.makedirs(post_data_path, exist_ok=True)
# combine previous
combine_previous_path = post_data_path + "combine_previous/"
os.makedirs(combine_previous_path, exist_ok=True)
# combine previous output
combine_previous_from_output = combine_previous_path + "output/"
os.makedirs(combine_previous_from_output, exist_ok=True)
# copy folder in output to folder for previous data 
for root, dirnames, filenames in os.walk(rr.lammps_directory + 'output/'):
    for dirname in dirnames:
        folderpath = os.path.join(root, dirname)
        createfolderpath = folderpath.replace(rr.lammps_directory + 'output/', combine_previous_from_output)
        os.makedirs(createfolderpath, exist_ok=True)

# hdf5
hdf5_csv_path = post_process_path + 'hdf5_csv/'
f_thermo = hdf5_csv_path + "thermo.h5"
f_custom = hdf5_csv_path + "custom.h5"
f_pair = hdf5_csv_path + "pair.h5"
f_cipcj = hdf5_csv_path + "cipcj.h5"

# debug
debug_print_path = post_process_path + 'debug/'
os.makedirs(debug_print_path, exist_ok=True)
debug_fig_path = debug_print_path + 'fig/'
debug_fig_thermo_path = debug_fig_path + 'thermo/'
debug_fig_oneatom_path = debug_fig_path + 'oneatom/'
debug_fig_atomij_path = debug_fig_path + 'atomij/'

# interactive
interactive_path = post_process_path + 'interactive/'
# diagram
diagram_path = post_process_path + 'diagram/'
os.makedirs(diagram_path, exist_ok=True)
# momentum
f_momentum_mass_field_path = diagram_path + "momentum_mass_field/"
os.makedirs(f_momentum_mass_field_path, exist_ok=True)
# stress
f_stress_field_path = diagram_path + "stress/"
os.makedirs(f_stress_field_path, exist_ok=True)
# mu I
f_mu_I_path = diagram_path + "mu_I/"
os.makedirs(f_mu_I_path, exist_ok=True)
# same scale for each plot
f_momentum_mass_field_samescale_path = f_momentum_mass_field_path + "same_scale/"
os.makedirs(f_momentum_mass_field_samescale_path, exist_ok=True)
f_stress_field_samescale_path = f_stress_field_path + "same_scale/"
os.makedirs(f_stress_field_samescale_path, exist_ok=True)

# stress path
def stress_folder_path(i,j):
    return f_stress_field_samescale_path + str(i) + str(j) + "/"
def stress_fix_k_folder_path(i,j,k):
    return f_stress_field_samescale_path + str(i) + str(j) + "_k_" + str(k) + "/"
def stress_time_folder_path(i,j):
    return f_stress_field_samescale_path + str(i) + str(j) + "_vs_time/"
# create stress folder
for [i, j] in [
                [0, 0],
                [1, 1],
                [2, 2],
                [0, 1],
                [0, 2],
                [1, 2],
                ]:
    os.makedirs(stress_folder_path(i,j), exist_ok=True)
# wall stress
def wall_stress_fix_k_folder_path(i,j,k):
    return f_stress_field_samescale_path + "wall/" + str(i) + str(j) + "_k_" + str(k) + "/"

if if_plot_velocity_field_scale_same == "yes":
    if rr.logfile["shearwall"] == "zcylinder":
        f_momentum_mass_field_v23x23_path = f_momentum_mass_field_samescale_path + "Vr_Vz_Xr_Xz/"
        f_momentum_mass_field_v13x23_path = f_momentum_mass_field_samescale_path + "Vr_Vtheta_Xr_Xz/"
        os.makedirs(f_momentum_mass_field_v23x23_path, exist_ok=True)
        os.makedirs(f_momentum_mass_field_v13x23_path, exist_ok=True)
    if rr.logfile["shearwall"] == "yplane":
        f_momentum_mass_field_v23x23_path = f_momentum_mass_field_samescale_path + "Vx_Vz_Xx_Xz/"
        f_momentum_mass_field_v13x23_path = f_momentum_mass_field_samescale_path + "Vx_Vtheta_Xx_Xz/"
        os.makedirs(f_momentum_mass_field_v23x23_path, exist_ok=True)
        os.makedirs(f_momentum_mass_field_v13x23_path, exist_ok=True)

    f_momentum_mass_field_volumnfraction_x23_path = f_momentum_mass_field_samescale_path + "density_Xx_Xz/"
    os.makedirs(f_momentum_mass_field_volumnfraction_x23_path, exist_ok=True)
else:
    if rr.logfile["shearwall"] == "zcylinder":
        f_momentum_mass_field_v23x23_path = f_momentum_mass_field_path + "Vr_Vz_Xr_Xz/"
        f_momentum_mass_field_v13x23_path = f_momentum_mass_field_path + "Vr_Vtheta_Xr_Xz/"
        os.makedirs(f_momentum_mass_field_v23x23_path, exist_ok=True)
        os.makedirs(f_momentum_mass_field_v13x23_path, exist_ok=True)
    if rr.logfile["shearwall"] == "yplane":
        f_momentum_mass_field_v23x23_path = f_momentum_mass_field_path + "Vx_Vz_Xx_Xz/"
        f_momentum_mass_field_v13x23_path = f_momentum_mass_field_path + "Vx_Vtheta_Xx_Xz/"
        os.makedirs(f_momentum_mass_field_v23x23_path, exist_ok=True)
        os.makedirs(f_momentum_mass_field_v13x23_path, exist_ok=True)

    f_momentum_mass_field_volumnfraction_x23_path = f_momentum_mass_field_path + "density_Xx_Xz/"
    os.makedirs(f_momentum_mass_field_volumnfraction_x23_path, exist_ok=True) 

f_max_velocity_near_wall = diagram_path + "max_velocity_near_wall/"

f_wall_force_plot_path = diagram_path + "wall_force/"

f_strain_rate_path = diagram_path + "strain_rate/"
os.makedirs(f_strain_rate_path, exist_ok=True)

f_ek_path = diagram_path + "ek/"
os.makedirs(f_ek_path, exist_ok=True)

f_velocity_path = diagram_path + "velocity/"
os.makedirs(f_velocity_path, exist_ok=True)

f_ekminusekave_path = diagram_path + "ekminusekave/"
os.makedirs(f_ekminusekave_path, exist_ok=True)

f_fraction_check_everygrid = diagram_path + "fraction_check_everygrid/"

f_velocity_i_time_ave_j_fix_k_ave_path = diagram_path + "velocity_i_time_ave_j_fix_k_ave/"

def f_path_strain_rate_i_j_ave_k(i,j,k):
    f_path = f_strain_rate_path + "strain_rate_" + str(i) + str(j) + "ave" + str(k) + "/"
    return f_path

def f_path_ekovermass_i_j_ave_k(i,j,k):
    f_path = f_ek_path + "ek" + str(i) + "_ave" + str(j) + "_x" + str(k) + "/"
    return f_path

def f_path_ekminusekaveovermass_i_j_ave_k(i,j,k):
    f_path = f_ekminusekave_path + "ekminusekav" + str(i) + "_ave" + str(k) + "_x" + str(j) + "/"
    return f_path

def f_path_strain_rate_i_j_x23(i,j):
    f_path = f_strain_rate_path + "strain_rate_" + str(i) + str(j) + "_x23" + "/"
    return f_path

def f_path_velocity_i_j_ave_k(i,j,k):
    f_path = f_velocity_path + "velocity" + str(i) + "_ave" + str(j) + "_x" + str(k) + "/"
    return f_path

def f_path_velocity_i_j_ave_k_no_top(i,j,k):
    f_path = f_velocity_path + "velocity_no_top" + str(i) + "_ave" + str(j) + "_x" + str(k) + "/"
    return f_path

def f_path_ek_i_j_ave_k_no_top(i,j,k):
    f_path = f_ek_path + "ek_no_top" + str(i) + "_ave" + str(j) + "_x" + str(k) + "/"
    return f_path

# Latex report
latex_path = post_process_path + 'latex/'
os.makedirs(latex_path, exist_ok=True)
latex_pics_path = latex_path + 'pics/'
os.makedirs(latex_pics_path, exist_ok=True)