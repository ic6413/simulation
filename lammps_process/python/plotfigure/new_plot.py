# read lammps output
# transform data to array
# save array
# load saved array
# calculate new defined variable
# save new variable as array
# plot from array
#===========================================
# import site package
import os
import sys
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datapath as dp
# import my package
import read_setting as rr
import plotfigure.new_data_plot as pn

# setting matplotlib
# Here are some matplotlib style settings that I (JSW) find useful for publication-quality plots
# http://tobackgroup.physics.tamu.edu/useful/matplotlib-formatting.txt
# First, we set the width of a figure to be exactly the width of a single column in PRL or PRD (in
# inches), and set the aspect ratio to be 4:3.

plt.rc('figure', figsize=(3.375, 3.375*(.75)), titlesize=10)

# This line sets the resolution (number of dots per inch) for saving figures to raster formats like
# PNG.  It's mainly useful to set the size of the figures displayed in the webbrowser in the Jupyter
# or IPython notebook
plt.rc('savefig', dpi=180)
#plt.rc('axes.formatter', limits=(0,0), use_mathtext=True)
plt.rc('savefig', bbox='tight')
plt.ticklabel_format(style='sci', axis='both')

# The default for scatter plots is to put three points in each legend entry.  Usually one is
# sufficient.
#plt.rc('legend', scatterpoints=1, fontsize=10, frameon=True)

# We would like the text to be created using TeX, so we get full support for LaTeX math syntax.
#plt.rc('text', usetex=True)

# We want all the text to be Computer Modern, and 10pt.  This way, combined with the correct setting
# of the size of the figure, the text in the figure will exactly match the text in the rest of the
# manuscript.

plt.rc('font', size=8, weight='normal', family='sans-serif', serif=['Helvetica', 'Arial'])

plt.rc('axes', labelsize=10, titlesize=8, labelweight='normal')

# http://aeturrell.com/2018/01/31/publication-quality-plots-in-python/
#plt.rc('xtick', labelsize=10)
#plt.rc('ytick', labelsize=10)
plt.rc('figure', autolayout=False)
plt.rc('lines', linewidth=2)
plt.rc('lines', markersize=4)
plt.rc('mathtext', fontset="stix")

# more setting
# http://physicalmodelingwithpython.blogspot.com/2015/06/making-plots-for-publication.html


def time_from_step_0(step):
    return step*float(rr.logfile["ts"])

def time_from_start_rotate(step):
    return time_from_step_0(step)-rr.logfile["rotate_start_time"]

def strain_from_rotate_start(step):
    shear_rate = float(rr.logfile['in_velocity'])/(float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp']))
    return time_from_start_rotate(step)*shear_rate


def transfer_time_to_str(timearray):
    if isinstance(timearray, str):
        if timearray == 'all':
            string = 'all'
    else:
        string = "-".join([
            "{:.2e}".format(a) for a in [timearray[0], timearray[-1]]
        ])
    return string

def transfer_coor_to_str(coordarray, ifave=False):
    if isinstance(coordarray, str):
        if coordarray == 'all':
            string = 'all'
        elif coordarray == 'ave':
            string = 'ave'
        elif coordarray == 'sumdividemass':
            string = 'sumdividemass'
    else:
        if ifave:
            string ="ave" + "-".join([
                "{:d}".format(a) for a in [coordarray[0], coordarray[-1]]
            ])
        else:
            string ="-".join([
                "{:d}".format(a) for a in [coordarray[0], coordarray[-1]]
            ])
    return string

def get_value_nochunk_include_ts_t_st(lmp_path, name, n_ave, inputstepsarray, y_name):

    if name == "timestep":
        value = inputstepsarray
    elif name == "time":
        value = time_from_start_rotate(inputstepsarray)
    elif name == "strain":
        value = strain_from_rotate_start(inputstepsarray)
    elif name == 'Coord1' or name == 'Coord2':
        char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
        value = pn.load_coord(lmp_path, char, name)
    else:
        value = get_value_nochunk(lmp_path, name, n_ave, inputstepsarray)
    return value
# get value for timestep time strain
def get_value_include_ts_t_st(lmp_path, name, n_ave, inputstepsarray, y_name):

    if name == "timestep":
        value = inputstepsarray
    elif name == "time":
        value = time_from_start_rotate(inputstepsarray)
    elif name == "strain":
        value = strain_from_rotate_start(inputstepsarray)
    elif name == 'Coord1' or name == 'Coord2':
        char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
        value = pn.load_coord(lmp_path, char, name)
    else:
        value = get_value(lmp_path, name, n_ave, inputstepsarray)
    return value
# get std value for timestep time strain
def get_std_value_include_ts_t_st(lmp_path, name, n_ave, inputstepsarray, y_name):
    if name == "timestep":
        value_std = 0
    elif name == "time":
        value_std = 0
    elif name == "strain":
        value_std = 0
    elif name == 'Coord1' or name == 'Coord2':
        value_std = 0
    else:
        value_std = get_std_value(lmp_path, name, n_ave, inputstepsarray)
    return value_std

def get_value_nochunk(lmp_path, name, n_ave, inputstepsarray):
    if n_ave % 2 == 1:
        pass
    else:
        sys.exit("n_ave should be odd integer")
    # combine with previous
    # pick the smallest step and max step
    # pick the smallest step and largest step in this simulation
    
    def get_step_first_end(lmp_path_folder):
        timestep_local = np.load(lmp_path_folder + pn.nochunk_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        
        timestep_first = timestep_local[int((n_ave-1)/2)]
        timestep_last = timestep_local[-1-int((n_ave-1)/2)]
        return (timestep_first, timestep_last)
    # input step
    step_first = inputstepsarray[0]
    # get all timestep_first timestep_last from all pre simu
    # calculate how many pre simu we need
    n_include_pre_simu = 0
    for i in range(rr.n_simu_total):
        n_include_pre_simu = n_include_pre_simu + 1
        if step_first >= get_step_first_end(
            rr.folder_path_list_initial_to_last[rr.n_simu_total-1-i]
            )[0]:
            
            break
    timestep = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + pn.nochunk_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
    
    array = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + pn.nochunk_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
    for j in range(n_include_pre_simu-1):
        timestep_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + pn.nochunk_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        array_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + pn.nochunk_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
        timestep = np.append(timestep, timestep_append, axis=0)
        array = np.append(array, array_append, axis=0)
    d_step = timestep[1] - timestep[0]
    indexesarray = (inputstepsarray-timestep[0])/d_step
    indexesarray = indexesarray.astype(int)
    # average
    sum_over_n_array = 0
    m = int((n_ave-1)/2)
    for k in range(-m, m+1):
        sum_over_n_array = sum_over_n_array + array[indexesarray + k]
    ave_array = sum_over_n_array/n_ave
    return ave_array

def get_value(lmp_path, name, n_ave, inputstepsarray):
    if n_ave % 2 == 1:
        pass
    else:
        sys.exit("n_ave should be odd integer")
    # combine with previous
    # pick the smallest step and max step
    # pick the smallest step and largest step in this simulation
    
    def get_step_first_end(lmp_path_folder):
        timestep_local = np.load(lmp_path_folder + pn.c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        
        timestep_first = timestep_local[int((n_ave-1)/2)]
        timestep_last = timestep_local[-1-int((n_ave-1)/2)]
        return (timestep_first, timestep_last)
    # input step
    step_first = inputstepsarray[0]
    # get all timestep_first timestep_last from all pre simu
    # calculate how many pre simu we need
    n_include_pre_simu = 0
    for i in range(rr.n_simu_total):
        n_include_pre_simu = n_include_pre_simu + 1
        if step_first >= get_step_first_end(
            rr.folder_path_list_initial_to_last[rr.n_simu_total-1-i]
            )[0]:
            
            break
    timestep = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + pn.c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
    
    array = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + pn.c_r_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
    for j in range(n_include_pre_simu-1):
        timestep_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + pn.c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        array_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + pn.c_r_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
        timestep = np.append(timestep, timestep_append, axis=0)
        array = np.append(array, array_append, axis=0)
    d_step = timestep[1] - timestep[0]
    indexesarray = (inputstepsarray-timestep[0])/d_step
    indexesarray = indexesarray.astype(int)
    # average
    sum_over_n_array = 0
    m = int((n_ave-1)/2)
    for k in range(-m, m+1):
        sum_over_n_array = sum_over_n_array + array[indexesarray + k]
    ave_array = sum_over_n_array/n_ave
    return ave_array

def get_std_value(lmp_path, name, n_ave, inputstepsarray):
    if n_ave % 2 == 1:
        pass
    else:
        sys.exit("n_ave should be odd integer")
    # combine with previous
    # pick the smallest step and max step
    # pick the smallest step and largest step in this simulation
    
    def get_step_first_end(lmp_path_folder):
        timestep_local = np.load(lmp_path_folder + pn.c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        
        timestep_first = timestep_local[int((n_ave-1)/2)]
        timestep_last = timestep_local[-1-int((n_ave-1)/2)]
        return (timestep_first, timestep_last)
    # input step
    step_first = inputstepsarray[0]
    # get all timestep_first timestep_last from all pre simu
    # calculate how many pre simu we need
    n_include_pre_simu = 0
    for i in range(rr.n_simu_total):
        n_include_pre_simu = n_include_pre_simu + 1
        if step_first >= get_step_first_end(
            rr.folder_path_list_initial_to_last[rr.n_simu_total-1-i]
            )[0]:
            
            break
    timestep = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + pn.c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
    std_array = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total-n_include_pre_simu] + pn.c_r_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
    for j in range(n_include_pre_simu-1):
        timestep_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + pn.c_r_npyfilepath_coord_char[name]["timestep_path"], mmap_mode='r')
        array_append = np.load(rr.folder_path_list_initial_to_last[rr.n_simu_total - n_include_pre_simu + j + 1] + pn.c_r_npyfilepath_coord_char[name]['npyfilepath'], mmap_mode='r')
        timestep = np.append(timestep, timestep_append, axis=0)
        std_array = np.append(std_array, array_append, axis=0)
    d_step = timestep[1] - timestep[0]
    indexesarray = (inputstepsarray-timestep[0])/d_step
    indexesarray = indexesarray.astype(int)
    # average
    std_array_sq = std_array**2
    sum_over_n_array_sq = 0
    for k in range(-(n_ave-1)/2, (n_ave-1)/2+1):
        sum_over_n_array_sq = sum_over_n_array_sq + std_array_sq[indexesarray + k]
    std_ave_array = (sum_over_n_array_sq/n_ave)**0.5
    return std_ave_array

lmp_path = rr.folder_path_list_initial_to_last[-1]
"""
# check if all timestep the same for chunk
for key in pn.map_chunkfile_char_save_folderpath.keys():
    path = lmp_path + pn.map_chunkfile_char_save_folderpath[key]["timestep_path"]
    if np.any(
        np.load(path, mmap_mode='r') != np.load(lmp_path + pn.map_chunkfile_char_save_folderpath["mv_Ek_mass"]["timestep_path"], mmap_mode='r')
    ):
        print("timestep file not match for all chunk file")
        breakpoint()
"""

stress_scale_width = float(rr.logfile['den'])*float(rr.logfile['g'])*float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp'])
stress_scale_width_str = r'$\rho_s g w $'
stress_scale_height = float(rr.logfile['den'])*float(rr.logfile['g'])*float(rr.logfile['z_length_create_dp_unit'])*float(rr.logfile['dp'])
stress_scale_height_str = r'$\rho_s g h $'
inwall_area = (
    float(rr.logfile['x_period_dp_unit'])
    *float(rr.logfile['z_length_create_dp_unit'])
    *float(rr.logfile['dp'])**2
)
bottom_area = (
    float(rr.logfile['width_wall_dp_unit'])
    *float(rr.logfile['x_period_dp_unit'])
    *float(rr.logfile['dp'])**2
)
velocity_scale = float(rr.logfile['in_velocity'])
if velocity_scale < 0:
    velocity_scale = -velocity_scale
velocity_scale_str = r'$V_{inwall}$'
strain_rate_scale = velocity_scale/(float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp']))
strain_rate_scale_str = r'$V_{inwall}/Width$'
coord_scale = float(rr.logfile['dp'])
coord_scale_str = r'$d_p$'
mu_scale = -1
mu_scale_str = None #r'$-1$'
I_scale = -1
I_scale_str = r'$-1$'
mu_tensor_scale = 1
mu_tensor_scale_str = None
I_tensor_scale = 1
I_tensor_scale_str = None

# scale dictionary
scale = {
    'stress': stress_scale_width,
    'velocity': velocity_scale,
    'strain_rate': strain_rate_scale,
    'coord_scale': coord_scale,
    'mu': mu_scale,
    'I': I_scale,
    'mu_tensor': mu_tensor_scale,
    'I_tensor': I_tensor_scale,
    'Coord': float(rr.logfile['dp']),
}
# variable_name to string in pic label or legend
v_name_to_labal_str = {
    "pressure":"P",
    "n_contact":"Z",
    "velocity_1": "V",
    "velocity_2": "V",
    "velocity_3": "V",
    "strain": r'$\gamma$',
    "inwall_force_1": r'$F_x$',
    "inwall_force_2": r'$F_y$',
    "inwall_force_3": r'$F_z$',
    "outwall_force_1": r'$F_x$',
    "outwall_force_2": r'$F_y$',
    "outwall_force_3": r'$F_z$',
    "zbottom_force_1": r'$F_x$',
    "zbottom_force_2": r'$F_y$',
    "zbottom_force_3": r'$F_z$',
    "strain_rate_21": r'$\dot{\gamma_{21}}$',
    "strain_rate_22": r'$\dot{\gamma_{22}}$',
    "strain_rate_23": r'$\dot{\gamma_{23}}$',
    "strain_rate_31": r'$\dot{\gamma_{31}}$',
    "strain_rate_32": r'$\dot{\gamma_{32}}$',
    "strain_rate_33": r'$\dot{\gamma_{33}}$',
    "strain_rate_21_middle": r'$\dot{\gamma_{21}}$',
    "strain_rate_22_middle": r'$\dot{\gamma_{22}}$',
    "strain_rate_23_middle": r'$\dot{\gamma_{23}}$',
    "strain_rate_31_middle": r'$\dot{\gamma_{31}}$',
    "strain_rate_32_middle": r'$\dot{\gamma_{32}}$',
    "strain_rate_33_middle": r'$\dot{\gamma_{33}}$',
    "stress_11": r'$\sigma_{11}$',
    "stress_22": r'$\sigma_{22}$',
    "stress_33": r'$\sigma_{33}$',
    "stress_12": r'$\sigma_{12}$',
    "stress_13": r'$\sigma_{13}$',
    "stress_23": r'$\sigma_{23}$',
    "mu_12": r'$\mu_{12}$',
    "mu_13": r'$\mu_{13}$',
    "mu_23": r'$\mu_{23}$',
    "mu_12_middle": r'$\mu_{12}$',
    "mu_13_middle": r'$\mu_{13}$',
    "mu_23_middle": r'$\mu_{23}$',
    "mu_tensor_12": r'$\mu_{12}$',
    "mu_tensor_13": r'$\mu_{13}$',
    "mu_tensor_23": r'$\mu_{23}$',
    "I_12": r'$I_{12}$',
    "I_13": r'$I_{13}$',
    "I_23": r'$I_{23}$',
    "I_tensor": r'$I$',
    "fraction": r'$\phi$',
    "inwall_stress_1": r'$\sigma_{21}$',
    "inwall_stress_2": r'$\sigma_{22}$',
    "inwall_stress_3": r'$\sigma_{23}$',
    "outwall_stress_1": r'$\sigma_{21}$',
    "outwall_stress_2": r'$\sigma_{22}$',
    "outwall_stress_3": r'$\sigma_{23}$',
    "zbottom_stress_1": r'$\sigma_{31}$',
    "zbottom_stress_2": r'$\sigma_{32}$',
    "zbottom_stress_3": r'$\sigma_{33}$',
    "Coord1": r'y',
    "Coord2": r'z',
}

def transfer_v_name_to_label_str_in_figure(v_name):
    if v_name in v_name_to_labal_str.keys():
        v_name_label_in_figure = v_name_to_labal_str[v_name]
    else:
        v_name_label_in_figure = v_name
    return v_name_label_in_figure

# sum mean select (plus minus multiply divide) operation for variable coord1 coord2 time
def sum_variable(variable):

    return variable

def plot_1D_from_chunk2D(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray, coord1_index_array, coord2_index_array,
        ave_over_coord1=False, ave_over_coord2=False,
        figure_class=None, legend_class=None, spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        x_scale = 'linear', y_scale = 'linear',
        useerrorbar = True,
        ifdivideNcount = False,
        maskstatic=None,
        masknonstatic=None,
    ):
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    # str for variable used in figure label legend
    
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')

    if coord2_index_array=="sumdividemass" and y_name=="mv_1":
        y_name_label_in_figure = "V"
    # legend_class: time, coord1, coord2
    # spaceave: coord1, coord2
    char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord1 = pn.load_coord(lmp_path, char, "Coord1")
    Coord2 = pn.load_coord(lmp_path, char, "Coord2")
    x_value = get_value_include_ts_t_st(lmp_path, x_name, n_ave, inputstepsarray, y_name)
    x_value_std = get_std_value_include_ts_t_st(lmp_path, x_name + "_std", n_ave, inputstepsarray, y_name)
    y_value = get_value_include_ts_t_st(lmp_path, y_name, n_ave, inputstepsarray, y_name)
    y_value_std = get_std_value_include_ts_t_st(lmp_path, y_name + "_std", n_ave, inputstepsarray, y_name)
    if ifdivideNcount:
        Ncount_value = get_value_include_ts_t_st(lmp_path, 'Ncount', n_ave, inputstepsarray, y_name)
    if coord2_index_array=="sumdividemass":
        mass_value = get_value_include_ts_t_st(lmp_path, 'mass', n_ave, inputstepsarray, y_name)
    # scale factor
    x_value = x_value/x_scale_factor
    x_value_std = x_value_std/x_scale_factor
    y_value = y_value/y_scale_factor
    y_value_std = y_value_std/y_scale_factor

    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    if ifdivideNcount:
        subfoldername = y_name + "_divideN_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    if useerrorbar:
        os.makedirs(diagram_path_add_nve + subfoldername + "errorbar/", exist_ok=True)

    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    # title

    if x_scale_str is None:
        x_label_str = x_name_label_in_figure
    else:
        x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
    ax.set_xlabel(x_label_str)
    
    if y_scale_str is None:
        y_label_str = y_name_label_in_figure
    else:
        y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
    ax.set_ylabel(y_label_str)
    if ifdivideNcount and y_label_str=="n_contact":
        ax.set_ylabel("<Z>")
    
    
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    
    # legend_class: time, coord1, coord2
    if legend_class == 'tc1c2':
        for indexstep, step in enumerate(inputstepsarray):
            for coord1_index in coord1_index_array:
                for coord2_index in coord2_index_array:
                    time = time_from_start_rotate(step)
                    strain = strain_from_rotate_start(step)
                    label = ", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                        "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                        "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[indexstep]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[indexstep, coord1_index, coord2_index]
                        x_value_std_plot = x_value_std[indexstep, coord1_index, coord2_index]
                    
                    y_value_plot = y_value[indexstep, coord1_index, coord2_index]
                    y_value_std_plot = y_value_std[indexstep, coord1_index, coord2_index]
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[indexstep, coord1_index, coord2_index]
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c1c2':
        for coord1_index in coord1_index_array:
            for coord2_index in coord2_index_array:
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                ])

                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index, coord2_index]
                    x_value_std_plot = x_value_std[:, coord1_index, coord2_index]
                
                y_value_plot = y_value[:, coord1_index, coord2_index]
                y_value_std_plot = y_value_std[:, coord1_index, coord2_index]
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index, coord2_index]
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    elif legend_class == 'c1':
        if coord2_index_array=='ave':
            if ave_over_coord1:
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :]
                    x_value_plot = np.nanmean(x_value_plot, axis=(1,2))
                    x_value_std_plot = np.nanmean(x_value_std_plot, axis=(1,2))
                y_value_plot = y_value[:, coord1_index_array, :]
                y_value_std_plot = y_value_std[:, coord1_index_array, :]
                y_value_plot = np.nanmean(y_value_plot, axis=(1,2))
                y_value_std_plot = np.nanmean(y_value_std_plot, axis=(1,2))
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "({})".format(coord_scale_str),
                    "z=average",
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord1_index in coord1_index_array:
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y={:.1f} ({})".format(Coord1[coord1_index, 0]/coord_scale, coord_scale_str),
                        "z=average",
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, coord1_index, :]
                        x_value_std_plot = x_value_std[:, coord1_index, :]
                        #breakpoint()
                        x_value_plot = np.nanmean(x_value_plot, axis=1)
                        x_value_std_plot = np.nanmean(x_value_std_plot, axis=1)
                    y_value_plot = y_value[:, coord1_index, :]
                    y_value_std_plot = y_value_std[:, coord1_index, :]
                    y_value_plot = np.nanmean(y_value_plot, axis=1)
                    y_value_std_plot = np.nanmean(y_value_std_plot, axis=1)
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, coord1_index, :]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=1)
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c2':
        if coord1_index_array=='ave':
            if ave_over_coord2:
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot, axis=(1,2))
                    x_value_std_plot = np.nanmean(x_value_std_plot, axis=(1,2))
                y_value_plot = y_value[:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot, axis=(1,2))
                y_value_std_plot = np.nanmean(y_value_std_plot, axis=(1,2))
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),"y=average",
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord2_index in coord2_index_array:
                    
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y=average",
                        "z={:.1f} ({})".format(Coord2[0, coord2_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, :, coord2_index]
                        x_value_std_plot = x_value_std[:, :, coord2_index]
                        x_value_plot = np.nanmean(x_value_plot, axis=1)
                        x_value_std_plot = np.nanmean(x_value_std_plot, axis=1)
                    y_value_plot = y_value[:, :, coord2_index]
                    y_value_std_plot = y_value_std[:, :, coord2_index]
                    y_value_plot = np.nanmean(y_value_plot, axis=1)
                    y_value_std_plot = np.nanmean(y_value_std_plot, axis=1)
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, :, coord2_index]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=1)
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 't':
        if coord2_index_array=="sumdividemass":
            if x_name == "Coord1":
                for index_step, step in enumerate(inputstepsarray):
                    strain = strain_from_rotate_start(step)
                    label=", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                    ])
                    x_value_plot = x_value[coord1_index_array, 0]
                    x_value_std_plot = 0
                    y_value_plot = y_value[index_step, coord1_index_array, :]
                    y_value_std_plot = y_value_std[index_step, coord1_index_array, :]
                    mass_value_plot = mass_value[index_step, coord1_index_array, :]
                    
                    y_value_plot = np.nanmean(y_value_plot, axis=1)
                    y_value_std_plot = np.nanmean(y_value_std_plot, axis=1)
                    mass_value_plot = np.nanmean(mass_value_plot, axis=1)
                    y_value_plot = y_value_plot/mass_value_plot
                    y_value_std_plot = y_value_std_plot/mass_value_plot
                    
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=6,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=6,
                        )

    elif legend_class == None:
        if ave_over_coord1:
            if ave_over_coord2:
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot, axis=(1,2))
                    x_value_std_plot = np.nanmean(x_value_std_plot, axis=(1,2))
                y_value_plot = y_value[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot, axis=(1,2))
                y_value_std_plot = np.nanmean(y_value_std_plot, axis=(1,2))
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    if useerrorbar:
        fig.savefig(
        "_".join([
            diagram_path_add_nve + subfoldername + "errorbar/" + 'step',
            transfer_time_to_str(inputstepsarray),
            'coord1',
            transfer_coor_to_str(coord1_index_array, ave_over_coord1),
            'coord2',
            transfer_coor_to_str(coord2_index_array, ave_over_coord2),
            ".png",
        ]),
        format="png",
        )
    else:
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + 'step',
                transfer_time_to_str(inputstepsarray),
                'coord1',
                transfer_coor_to_str(coord1_index_array, ave_over_coord1),
                'coord2',
                transfer_coor_to_str(coord2_index_array, ave_over_coord2),
                ".png",
            ]),
            format="png",
        )
    # close figure after save
    plt.close('all')

def plot_1D_from_chunk2D_mask(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray, coord1_index_array, coord2_index_array,
        ave_over_coord1=False, ave_over_coord2=False,
        figure_class=None, legend_class=None, spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        x_scale = 'linear', y_scale = 'linear',
        useerrorbar = True,
        ifdivideNcount = False,
        ifmaskstatic=False,
        ifmasknonstatic=True,
    ):
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    # str for variable used in figure label legend
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
    # legend_class: time, coord1, coord2
    # spaceave: coord1, coord2
    x_value = get_value_include_ts_t_st(lmp_path, x_name, n_ave, inputstepsarray, y_name)
    x_value_std = get_std_value_include_ts_t_st(lmp_path, x_name + "_std", n_ave, inputstepsarray, y_name)
    y_value = get_value_include_ts_t_st(lmp_path, y_name, n_ave, inputstepsarray, y_name)
    y_value_std = get_std_value_include_ts_t_st(lmp_path, y_name + "_std", n_ave, inputstepsarray, y_name)
    if ifdivideNcount:
        Ncount_value = get_value_include_ts_t_st(lmp_path, 'Ncount', n_ave, inputstepsarray, y_name)
    # scale factor
    x_value = x_value/x_scale_factor
    x_value_std = x_value_std/x_scale_factor
    y_value = y_value/y_scale_factor
    y_value_std = y_value_std/y_scale_factor

    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    if ifdivideNcount:
        subfoldername = y_name + "_divideN_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    if useerrorbar:
        os.makedirs(diagram_path_add_nve + subfoldername + "errorbar/", exist_ok=True)

    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    # title

    if x_scale_str is None:
        x_label_str = x_name_label_in_figure
    else:
        x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
    ax.set_xlabel(x_label_str)
    
    if y_scale_str is None:
        y_label_str = y_name_label_in_figure
    else:
        y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
    ax.set_ylabel(y_label_str)
    if ifdivideNcount and y_label_str=="n_contact":
        ax.set_ylabel("<Z>")
    
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord1 = pn.load_coord(lmp_path, char, "Coord1")
    Coord2 = pn.load_coord(lmp_path, char, "Coord2")
    x = Coord1/float(rr.logfile["dp"])
    y = Coord2/float(rr.logfile["dp"])
    maskstatic = (y-9) - (9-0.5)/(15.5-10.5)*(x-15.5) < 0
    masknonstatic = (y-25) - (25-0.5)/(15.5-5)*(x-5) > 0
    if ifmaskstatic:
        maskstaticornot = maskstatic
    if ifmasknonstatic:
        maskstaticornot = masknonstatic
    # legend_class: time, coord1, coord2
    if legend_class == 'tc1c2':
        for indexstep, step in enumerate(inputstepsarray):
            for coord1_index in coord1_index_array:
                for coord2_index in coord2_index_array:
                    strain = strain_from_rotate_start(step)
                    label = ", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                        "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                        "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[indexstep]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[indexstep, coord1_index, coord2_index]
                        x_value_std_plot = x_value_std[indexstep, coord1_index, coord2_index]
                    
                    y_value_plot = y_value[indexstep, coord1_index, coord2_index]
                    y_value_std_plot = y_value_std[indexstep, coord1_index, coord2_index]
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[indexstep, coord1_index, coord2_index]
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c1c2':
        for coord1_index in coord1_index_array:
            for coord2_index in coord2_index_array:
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                ])

                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index, coord2_index]
                    x_value_std_plot = x_value_std[:, coord1_index, coord2_index]
                
                y_value_plot = y_value[:, coord1_index, coord2_index]
                y_value_std_plot = y_value_std[:, coord1_index, coord2_index]
                
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index, coord2_index]
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    elif legend_class == 'c1':
        if coord2_index_array=='ave':
            if ave_over_coord1:
                maskstaticornot_plot = maskstaticornot[coord1_index_array, :]
                aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :]
                    x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_plot = y_value[:, coord1_index_array, :]
                y_value_std_plot = y_value_std[:, coord1_index_array, :]
                y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "({})".format(coord_scale_str),
                    "z=average",
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord1_index in coord1_index_array:
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y={:.1f} ({})".format(Coord1[coord1_index, 0]/coord_scale, coord_scale_str),
                        "z=average",
                    ])
                    maskstaticornot_plot = maskstaticornot[coord1_index, :]
                    aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, coord1_index, :]
                        x_value_std_plot = x_value_std[:, coord1_index, :]
                        #breakpoint()
                        x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_plot = y_value[:, coord1_index, :]
                    y_value_std_plot = y_value_std[:, coord1_index, :]
                    y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, coord1_index, :]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == 'c2':
        if coord1_index_array=='ave':
            if ave_over_coord2:
                maskstaticornot_plot = maskstaticornot[:, coord2_index_array]
                aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_plot = y_value[:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "y=average",
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
            else:
                for coord2_index in coord2_index_array:
                    strain1 = strain_from_rotate_start(inputstepsarray[0])
                    strain2 = strain_from_rotate_start(inputstepsarray[-1])
                    label=", ".join([
                        r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                        "y=average",
                        "z={:.1f} ({})".format(Coord2[0, coord2_index]/coord_scale, coord_scale_str),
                    ])
                    maskstaticornot_plot = maskstaticornot[:, coord2_index]
                    aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, :, coord2_index]
                        x_value_std_plot = x_value_std[:, :, coord2_index]
                        x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_plot = y_value[:, :, coord2_index]
                    y_value_std_plot = y_value_std[:, :, coord2_index]
                    y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                    y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=1)/aveNratio
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, :, coord2_index]
                        Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=1)/aveNratio
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
    elif legend_class == None:
        if ave_over_coord1:
            if ave_over_coord2:
                maskstaticornot_plot = maskstaticornot[coord1_index_array, :][:, coord2_index_array]
                aveNratio = np.sum(maskstaticornot_plot)/np.sum(np.zeros_like(maskstaticornot_plot)+1)
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_std_plot = x_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                    x_value_plot = np.nanmean(x_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    x_value_std_plot = np.nanmean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_plot = y_value[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_std_plot = y_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                y_value_plot = np.nanmean(y_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                y_value_std_plot = np.nanmean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                if ifdivideNcount:
                    Ncount_value_plot = Ncount_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    Ncount_value_plot = np.nanmean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))/aveNratio
                    y_value_plot = y_value_plot/Ncount_value_plot
                    y_value_std_plot = y_value_std_plot/Ncount_value_plot
                zstring = "_".join(
                    ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                )
                ystring = "_".join(
                    ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                )
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    "average y=",
                    ystring,
                    "average z=",
                    zstring,
                    "({})".format(coord_scale_str),
                ])
                
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    if useerrorbar:
        fig.savefig(
        "_".join([
            diagram_path_add_nve + subfoldername + "errorbar/" + "mask_" + 'step',
            transfer_time_to_str(inputstepsarray),
            'coord1',
            transfer_coor_to_str(coord1_index_array, ave_over_coord1),
            'coord2',
            transfer_coor_to_str(coord2_index_array, ave_over_coord2),
            ".png",
        ]),
        format="png",
        )
    else:
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + "mask_" + 'step',
                transfer_time_to_str(inputstepsarray),
                'coord1',
                transfer_coor_to_str(coord1_index_array, ave_over_coord1),
                'coord2',
                transfer_coor_to_str(coord2_index_array, ave_over_coord2),
                ".png",
            ]),
            format="png",
        )
    # close figure after save
    plt.close('all')


def plot_1D_for_chunk1D_near_wall(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray, coord_index_array,
        figure_class=None, legend_class=None, spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        x_scale = 'linear', y_scale = 'linear',
        useerrorbar = True, 
        ifsum = False,
        setymin = False,
        ylim_bottom = 0,
        showlegend=True,
    ):
    # str for variable used in figure label legend
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
    # legend_class: time, coord    
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    # spaceave: coord

    x_value = get_value_include_ts_t_st(lmp_path, x_name, n_ave, inputstepsarray, y_name)
    x_value_std = get_std_value_include_ts_t_st(lmp_path, x_name + "_std", n_ave, inputstepsarray, y_name)
    y_value = get_value_include_ts_t_st(lmp_path, y_name, n_ave, inputstepsarray, y_name)
    y_value_std = get_std_value_include_ts_t_st(lmp_path, y_name + "_std", n_ave, inputstepsarray, y_name)
    # scale factor
    x_value = x_value/x_scale_factor
    x_value_std = x_value_std/x_scale_factor
    y_value = y_value/y_scale_factor
    y_value_std = y_value_std/y_scale_factor

    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    if useerrorbar:
        os.makedirs(diagram_path_add_nve + subfoldername + "errorbar/", exist_ok=True)

    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    # title
    
    
    if x_scale_str is None:
        x_label_str = x_name_label_in_figure
    else:
        x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
    ax.set_xlabel(x_label_str)
    if y_scale_str is None:
        y_label_str = y_name_label_in_figure
    else:
        y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
    ax.set_ylabel(y_label_str)
    
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    # check if only 1 coord
    if len(pn.c_r_npyfilepath_coord_char[char].keys()) == 1:
        pass
    else:
        print("the coord near wall is not 1")
        breakpoint()
        sys.exit()
    Coord_str = list(pn.c_r_npyfilepath_coord_char[char].keys())[0]
    Coord = pn.load_coord(lmp_path, char, Coord_str)
    if Coord_str == "Coord1":
        coord_str_in_plot = "y"
    elif Coord_str == "Coord2":
        coord_str_in_plot = "z"
    else:
        sys.exit("not y not z")
    # legend_class: time, coord1, coord2
    if ifsum==False:
        if legend_class == 'tc':
            for indexstep, step in enumerate(inputstepsarray):
                for coord_index in coord_index_array:
                    time = time_from_start_rotate(step)
                    strain = strain_from_rotate_start(time)
                    label = ", ".join([
                        r'$\gamma$' + "={:.2f}".format(strain),
                        coord_str_in_plot + "={:.1f} ({})".format(Coord[coord_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time" or x_name == "strain":
                        x_value_plot = x_value[indexstep]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[indexstep, coord_index]
                        x_value_std_plot = x_value_std[indexstep, coord_index]
                    
                    y_value_plot = y_value[indexstep, coord_index]
                    y_value_std_plot = y_value_std[indexstep, coord_index]
                    
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
        elif legend_class == 'c':
            for coord_index in coord_index_array:                    
                strain1 = strain_from_rotate_start(inputstepsarray[0])
                strain2 = strain_from_rotate_start(inputstepsarray[-1])
                label=", ".join([
                    r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
                    coord_str_in_plot + "={:.1f} ({})".format(Coord[coord_index]/coord_scale, coord_scale_str),
                ])
                if x_name == "timestep" or x_name == "time" or x_name == "strain":
                    x_value_plot = x_value[:]
                    x_value_std_plot = 0
                else:
                    x_value_plot = x_value[:, coord_index]
                    x_value_std_plot = x_value_std[:, coord_index]
                
                y_value_plot = y_value[:, coord_index]
                y_value_std_plot = y_value_std[:, coord_index]
                if useerrorbar:
                    ax.errorbar(
                        x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
                else:
                    ax.plot(
                        x_value_plot, y_value_plot,
                        label=label,
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                    )
    else:
        strain1 = strain_from_rotate_start(inputstepsarray[0])
        strain2 = strain_from_rotate_start(inputstepsarray[-1])
        label=", ".join([
            r'$\gamma$' + "={:.2f} to {:.2f}".format(strain1, strain2),
        ])
        if x_name == "timestep" or x_name == "time" or x_name == "strain":
            x_value_plot = x_value[:]
            x_value_std_plot = 0
        else:
            x_value_plot = x_value[:, :]
            x_value_std_plot = x_value_std[:, :]
            x_value_plot = np.average(x_value_plot, axis=1)
        
        y_value_plot = y_value[:, :]
        y_value_std_plot = y_value_std[:, :]
        
        y_value_plot = np.average(y_value_plot, axis=1)
        ax.plot(
            x_value_plot, y_value_plot,
            label=label,
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    if showlegend:
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    if setymin:
        ax.set_ylim(bottom=ylim_bottom)
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    if ifsum==False:
        if useerrorbar:
            fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + "errorbar/step",
                transfer_time_to_str(inputstepsarray),
                Coord_str,
                transfer_coor_to_str(coord_index_array),
                ".png",
            ]),
            format="png",
            )
        else:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + 'step',
                    transfer_time_to_str(inputstepsarray),
                    Coord_str,
                    transfer_coor_to_str(coord_index_array),
                    ".png",
                ]),
                format="png",
                )
    else:
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + 'step',
                transfer_time_to_str(inputstepsarray),
                ".png",
            ]),
            format="png",
        )
    # close figure after save
    plt.close('all')

def plot_quiver_from_chunk2D(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray,
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=0.1, label_scale=0.2,
        ifloglength=False,
        logvaluemin=10**-4,
        valueminus=0,
        ifplotseparateupdown=False,
        quiver_scale_up=0.1, label_scale_up=0.2,
        quiver_scale_down=0.1, label_scale_down=0.2,
        ifstreamplot=True,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=-10, vmax=10,
    ):
    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    if ifloglength:
        quiver_scale = np.log10(quiver_scale/logvaluemin)
        label_scale = np.log10(label_scale/logvaluemin)
        quiver_scale_up = np.log10(quiver_scale_up/logvaluemin)
        label_scale_up = np.log10(label_scale_up/logvaluemin)
        quiver_scale_down = np.log10(quiver_scale_down/logvaluemin)
        label_scale_down = np.log10(label_scale_down/logvaluemin)
    # str for variable used in figure label legend
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
    # spaceave: coord1, coord2
    time_array = time_from_start_rotate(inputstepsarray)
    
    x_value = get_value_include_ts_t_st(lmp_path, x_name, n_ave, inputstepsarray, y_name)
    y_value = get_value_include_ts_t_st(lmp_path, y_name, n_ave, inputstepsarray, y_name)
    mask = x_value>0.1
    x_value[mask] = x_value[mask] - valueminus
    mask = y_value>0.1
    y_value[mask] = y_value[mask] - valueminus
    # scale factor
    x_value = x_value/x_scale_factor
    y_value = y_value/y_scale_factor
    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    os.makedirs(diagram_path_add_nve + subfoldername + "streamplot/", exist_ok=True)
    if ifplotseparateupdown:
        os.makedirs(diagram_path_add_nve + subfoldername + "up/", exist_ok=True)
        os.makedirs(diagram_path_add_nve + subfoldername + "down/", exist_ok=True)
    if ifloglength:
        os.makedirs(diagram_path_add_nve + subfoldername + "log/", exist_ok=True)
        if ifplotseparateupdown:
            os.makedirs(diagram_path_add_nve + subfoldername + "log/" + "up/", exist_ok=True)
            os.makedirs(diagram_path_add_nve + subfoldername + "log/" + "down/", exist_ok=True)
    char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord1 = pn.load_coord(lmp_path, char, "Coord1")
    Coord2 = pn.load_coord(lmp_path, char, "Coord2")
    # plot quiver vector field
    for indexstep, step in enumerate(inputstepsarray):
        x_value_plot = x_value[indexstep]
        y_value_plot = y_value[indexstep]
        # plot ave_z velocity across y
        fig, ax = plt.subplots()
        
        # title

        if x_scale_str is None:
            x_label_str = x_name_label_in_figure
        else:
            x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
        ax.set_xlabel(x_label_str)
        
        if y_scale_str is None:
            y_label_str = y_name_label_in_figure
        else:
            y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
        ax.set_ylabel(y_label_str)            
        # plot
        if ifloglength:
            length = (x_value_plot**2 + y_value_plot**2)**0.5
            masktolog = (length>=logvaluemin)
            length_log = np.copy(length)
            length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
            #breakpoint()
            x_value_plot_log = length_log*x_value_plot/length
            y_value_plot_log = length_log*y_value_plot/length
            Q = ax.quiver(
                Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                cmap='hsv',
            )
        else:
            Q = ax.quiver(
                Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                cmap='hsv',
            )
            
        
        ax.set_xlabel(
            'y ({})'.format(coord_scale_str)
        )
        ax.set_ylabel(
            'z ({})'.format(coord_scale_str)
        )
        time = time_from_start_rotate(step)
        if ifloglength:
            ax.quiverkey(
                Q, 0.1, 0.95, label_scale,
                label = "logscale, : {:.2e} in 45 degree, ".format(label_scale) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time),
                labelpos='E',
                coordinates='figure', angle=45,
            )
        else:
            ax.quiverkey(
                Q, 0.1, 0.95, label_scale,
                label = " : {:.2e} in 45 degree, ".format(label_scale) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time),
                labelpos='E',
                coordinates='figure', angle=45,
            )
            
        plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
        plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

        if ifloglength:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + "log/" + 'step',
                    str(step),
                    ".png",
                ]),
                format="png",
                bbox_inches=None,
            )
        else:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + 'step',
                    str(step),
                    ".png",
                ]),
                format="png",
                bbox_inches=None,
            )

        # close figure after save
        plt.close('all')
    # plot up and down separate figure
    if ifplotseparateupdown:
        maskup = y_value>0
        maskdown = y_value<0
        x_value_up = np.zeros_like(x_value)
        y_value_up = np.zeros_like(y_value)
        x_value_down = np.zeros_like(x_value)
        y_value_down = np.zeros_like(y_value)
        x_value_up[maskup] = x_value[maskup]
        y_value_up[maskup] = y_value[maskup]
        x_value_down[maskdown] = x_value[maskdown]
        y_value_down[maskdown] = y_value[maskdown]
        for indexstep, step in enumerate(inputstepsarray):
            x_value_plot = x_value_up[indexstep]
            y_value_plot = y_value_up[indexstep]
            total_up = np.sum(y_value_plot)
            # plot ave_z velocity across y
            fig, ax = plt.subplots()
            
            # title

            if x_scale_str is None:
                x_label_str = x_name_label_in_figure
            else:
                x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
            ax.set_xlabel(x_label_str)
            
            if y_scale_str is None:
                y_label_str = y_name_label_in_figure
            else:
                y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
            ax.set_ylabel(y_label_str)            
            # plot
            if ifloglength:
                length = (x_value_plot**2 + y_value_plot**2)**0.5
                masktolog = (length>=logvaluemin)
                length_log = np.copy(length)
                length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
                #breakpoint()
                x_value_plot_log = length_log*x_value_plot/length
                y_value_plot_log = length_log*y_value_plot/length
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_up,
                    cmap='hsv',
                )
            else:
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_up,
                    cmap='hsv',
                )
                
            
            ax.set_xlabel(
                'y ({})'.format(coord_scale_str)
            )
            ax.set_ylabel(
                'z ({})'.format(coord_scale_str)
            )
            time = time_from_start_rotate(step)
            if ifloglength:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_up,
                    label = "up, logscale, : {:.2e} in 45 degree, ".format(label_scale_up) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total up is {:.2e}".format(total_up),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
            else:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_up,
                    label = " : up, {:.2e} in 45 degree, ".format(label_scale_up) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total up is {:.2e}".format(total_up),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
                
            plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
            plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

            if ifloglength:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "log/" + "up/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )
            else:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "up/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )

            # close figure after save
            plt.close('all')

        for indexstep, step in enumerate(inputstepsarray):
            x_value_plot = x_value_down[indexstep]
            y_value_plot = y_value_down[indexstep]
            total_down = np.sum(y_value_plot)
            # plot ave_z velocity across y
            fig, ax = plt.subplots()
            
            # title

            if x_scale_str is None:
                x_label_str = x_name_label_in_figure
            else:
                x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
            ax.set_xlabel(x_label_str)
            
            if y_scale_str is None:
                y_label_str = y_name_label_in_figure
            else:
                y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
            ax.set_ylabel(y_label_str)            
            # plot
            if ifloglength:
                length = (x_value_plot**2 + y_value_plot**2)**0.5
                masktolog = (length>=logvaluemin)
                length_log = np.copy(length)
                length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
                #breakpoint()
                x_value_plot_log = length_log*x_value_plot/length
                y_value_plot_log = length_log*y_value_plot/length
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_down,
                    cmap='hsv',
                )
            else:
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_down,
                    cmap='hsv',
                )
                
            
            ax.set_xlabel(
                'y ({})'.format(coord_scale_str)
            )
            ax.set_ylabel(
                'z ({})'.format(coord_scale_str)
            )
            time = time_from_start_rotate(step)
            if ifloglength:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_down,
                    label = "down, logscale, : {:.2e} in 45 degree, ".format(label_scale_down) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total down is {:.2e}".format(total_down),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
            else:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_down,
                    label = " : down, {:.2e} in 45 degree, ".format(label_scale_down) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total down is {:.2e}".format(total_down),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
                
            plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
            plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

            if ifloglength:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "log/" + "down/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )
            else:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "down/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )

            # close figure after save
            plt.close('all')
    # plot streamplot
    for indexstep, step in enumerate(inputstepsarray):
        x_value_plot = x_value[indexstep]
        y_value_plot = y_value[indexstep]
        # plot ave_z velocity across y
        fig, ax = plt.subplots()
        
        # title

        if x_scale_str is None:
            x_label_str = x_name_label_in_figure
        else:
            x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
        ax.set_xlabel(x_label_str)
        
        if y_scale_str is None:
            y_label_str = y_name_label_in_figure
        else:
            y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
        ax.set_ylabel(y_label_str)
        # check if line0
        if ifstreamplot:
            # plot
            strm = ax.streamplot(
                Coord1[:,0]/coord_scale, Coord2[0,:]/coord_scale, np.transpose(x_value_plot), np.transpose(y_value_plot),
                linewidth=1, color='k',
                density=[0.8, 0.8],
            )
        vector_length = np.sqrt(x_value_plot**2 + y_value_plot**2)
        vector_length = np.ma.masked_where(np.logical_not(vector_length > 0), vector_length)
        #breakpoint()
        d_Coord1 = Coord1[:,0][1]-Coord1[:,0][0]
        d_Coord2 = Coord2[0,:][1]-Coord2[0,:][0]
        Coord1_expand = Coord1[:,0] - d_Coord1/2
        Coord1_expand = np.append(Coord1_expand, (Coord1_expand[-1]+d_Coord1))
        Coord2_expand = Coord2[0,:] - d_Coord2/2
        Coord2_expand = np.append(Coord2_expand, (Coord2_expand[-1]+d_Coord2))
        if ifcontour:
            if contour_v_min_max == "constant":
                if contour_norm == "linear":
                    contour = ax.pcolor(
                        Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                        norm=colors.Normalize(vmin=vmin, vmax=vmax), #norm=colors.LogNorm(vmin=10**-1, vmax=10), #norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                        cmap='coolwarm',
                    )
                elif contour_norm == "log":
                    contour = ax.pcolor(
                    Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                    norm=colors.LogNorm(vmin=vmin, vmax=vmax), #norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                    cmap='coolwarm',
                )
                else:
                    sys.exit("Error: contour_norm not defined")
            elif contour_v_min_max == "min_to_max":
                if contour_norm == "linear":
                    contour = ax.pcolor(
                        Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                        norm=colors.Normalize(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                        cmap='coolwarm',
                    )
                elif contour_norm == "log":
                    contour = ax.pcolor(
                    Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                    norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                    cmap='coolwarm',
                )
                else:
                    sys.exit("Error: contour_norm not defined")
            else:
                sys.exit("Error: contour_v_min_max not defined")
            fig.colorbar(contour, ax=ax, extend='max')

        time = time_from_start_rotate(step)
        strain = strain_from_rotate_start(step)
        ax.set_title(
            "t = {:.2f} s".format(time) + ", strain is {:.2f}".format(strain)
        )
        ax.set_xlabel(
            'y ({})'.format(coord_scale_str)
        )
        ax.set_ylabel(
            'z ({})'.format(coord_scale_str)
        )

        plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
        plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + 'streamplot/' + 'step',
                str(step),
                ".png",
            ]),
            format="png",
            bbox_inches=None,
        )

        # close figure after save
        plt.close('all')


def plot_quiver_from_chunk2D_fraction(
        lmp_path,
        n_ave, x_name, y_name, inputstepsarray,
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=0.1, label_scale=0.2,
        ifloglength=False,
        logvaluemin=10**-4,
        valueminus=0,
        ifplotseparateupdown=False,
        quiver_scale_up=0.1, label_scale_up=0.2,
        quiver_scale_down=0.1, label_scale_down=0.2,
        ifstreamplot=True,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=-10, vmax=10,
        ave_y=1, ave_z=1,
    ):

    diagram_path_add_nve = dp.diagram_path + "n_ave_" + str(n_ave) + "/"
    if ifloglength:
        quiver_scale = np.log10(quiver_scale/logvaluemin)
        label_scale = np.log10(label_scale/logvaluemin)
        quiver_scale_up = np.log10(quiver_scale_up/logvaluemin)
        label_scale_up = np.log10(label_scale_up/logvaluemin)
        quiver_scale_down = np.log10(quiver_scale_down/logvaluemin)
        label_scale_down = np.log10(label_scale_down/logvaluemin)
    # str for variable used in figure label legend
    x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
    y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
    # spaceave: coord1, coord2
    time_array = time_from_start_rotate(inputstepsarray)
    
    x_value = get_value_include_ts_t_st(lmp_path, x_name, n_ave, inputstepsarray, y_name)
    y_value = get_value_include_ts_t_st(lmp_path, y_name, n_ave, inputstepsarray, y_name)
    mask = x_value>0.1
    x_value[mask] = x_value[mask] - valueminus
    mask = y_value>0.1
    y_value[mask] = y_value[mask] - valueminus
    # scale factor
    x_value = x_value/x_scale_factor
    y_value = y_value/y_scale_factor
    # subfolder name
    subfoldername = y_name + "_" + x_name + "/"
    os.makedirs(diagram_path_add_nve + subfoldername, exist_ok=True)
    os.makedirs(diagram_path_add_nve + subfoldername + "streamplot/", exist_ok=True)
    if ifplotseparateupdown:
        os.makedirs(diagram_path_add_nve + subfoldername + "up/", exist_ok=True)
        os.makedirs(diagram_path_add_nve + subfoldername + "down/", exist_ok=True)
    if ifloglength:
        os.makedirs(diagram_path_add_nve + subfoldername + "log/", exist_ok=True)
        if ifplotseparateupdown:
            os.makedirs(diagram_path_add_nve + subfoldername + "log/" + "up/", exist_ok=True)
            os.makedirs(diagram_path_add_nve + subfoldername + "log/" + "down/", exist_ok=True)
    char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord1 = pn.load_coord(lmp_path, char, "Coord1")
    Coord2 = pn.load_coord(lmp_path, char, "Coord2")
    
    x_value_ave = 0
    y_value_ave = 0
    Coord1_ave = 0
    Coord2_ave = 0
    n_y = Coord1.shape[0]
    n_z = Coord1.shape[1] 
    for i in range(ave_y):
        for j in range(ave_z):
            x_value_ave += x_value[:, i:i+1-ave_y+n_y, j:j+1-ave_z+n_z]
            y_value_ave += y_value[:, i:i+1-ave_y+n_y, j:j+1-ave_z+n_z]
            Coord1_ave += Coord1[i:i+1-ave_y+n_y, j:j+1-ave_z+n_z]
            Coord2_ave += Coord2[i:i+1-ave_y+n_y, j:j+1-ave_z+n_z]
    x_value_ave = x_value_ave/ave_y/ave_z
    y_value_ave = y_value_ave/ave_y/ave_z
    Coord1_ave = Coord1_ave/ave_y/ave_z
    Coord2_ave = Coord2_ave/ave_y/ave_z
    x_value = x_value_ave
    y_value = y_value_ave
    Coord1 = Coord1_ave
    Coord2 = Coord2_ave
    # plot quiver vector field
    for indexstep, step in enumerate(inputstepsarray):
        x_value_plot = x_value[indexstep]
        y_value_plot = y_value[indexstep]
        # plot ave_z velocity across y
        fig, ax = plt.subplots()
        
        # title

        if x_scale_str is None:
            x_label_str = x_name_label_in_figure
        else:
            x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
        ax.set_xlabel(x_label_str)
        
        if y_scale_str is None:
            y_label_str = y_name_label_in_figure
        else:
            y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
        ax.set_ylabel(y_label_str)            
        # plot
        if ifloglength:
            length = (x_value_plot**2 + y_value_plot**2)**0.5
            masktolog = (length>=logvaluemin)
            length_log = np.copy(length)
            length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
            #breakpoint()
            x_value_plot_log = length_log*x_value_plot/length
            y_value_plot_log = length_log*y_value_plot/length
            Q = ax.quiver(
                Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                cmap='hsv',
            )
        else:
            Q = ax.quiver(
                Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                cmap='hsv',
            )
            
        
        ax.set_xlabel(
            'y ({})'.format(coord_scale_str)
        )
        ax.set_ylabel(
            'z ({})'.format(coord_scale_str)
        )
        time = time_from_start_rotate(step)
        if ifloglength:
            ax.quiverkey(
                Q, 0.1, 0.95, label_scale,
                label = "logscale, : {:.2e} in 45 degree, ".format(label_scale) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time),
                labelpos='E',
                coordinates='figure', angle=45,
            )
        else:
            ax.quiverkey(
                Q, 0.1, 0.95, label_scale,
                label = " : {:.2e} in 45 degree, ".format(label_scale) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time),
                labelpos='E',
                coordinates='figure', angle=45,
            )
            
        plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
        plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

        if ifloglength:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + "log/" + 'step',
                    str(step),
                    ".png",
                ]),
                format="png",
                bbox_inches=None,
            )
        else:
            fig.savefig(
                "_".join([
                    diagram_path_add_nve + subfoldername + 'step',
                    str(step),
                    ".png",
                ]),
                format="png",
                bbox_inches=None,
            )

        # close figure after save
        plt.close('all')
    # plot up and down separate figure
    if ifplotseparateupdown:
        maskup = y_value>0
        maskdown = y_value<0
        x_value_up = np.zeros_like(x_value)
        y_value_up = np.zeros_like(y_value)
        x_value_down = np.zeros_like(x_value)
        y_value_down = np.zeros_like(y_value)
        x_value_up[maskup] = x_value[maskup]
        y_value_up[maskup] = y_value[maskup]
        x_value_down[maskdown] = x_value[maskdown]
        y_value_down[maskdown] = y_value[maskdown]
        for indexstep, step in enumerate(inputstepsarray):
            x_value_plot = x_value_up[indexstep]
            y_value_plot = y_value_up[indexstep]
            total_up = np.sum(y_value_plot)
            # plot ave_z velocity across y
            fig, ax = plt.subplots()
            
            # title

            if x_scale_str is None:
                x_label_str = x_name_label_in_figure
            else:
                x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
            ax.set_xlabel(x_label_str)
            
            if y_scale_str is None:
                y_label_str = y_name_label_in_figure
            else:
                y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
            ax.set_ylabel(y_label_str)            
            # plot
            if ifloglength:
                length = (x_value_plot**2 + y_value_plot**2)**0.5
                masktolog = (length>=logvaluemin)
                length_log = np.copy(length)
                length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
                #breakpoint()
                x_value_plot_log = length_log*x_value_plot/length
                y_value_plot_log = length_log*y_value_plot/length
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_up,
                    cmap='hsv',
                )
            else:
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_up,
                    cmap='hsv',
                )
                
            
            ax.set_xlabel(
                'y ({})'.format(coord_scale_str)
            )
            ax.set_ylabel(
                'z ({})'.format(coord_scale_str)
            )
            time = time_from_start_rotate(step)
            if ifloglength:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_up,
                    label = "up, logscale, : {:.2e} in 45 degree, ".format(label_scale_up) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total up is {:.2e}".format(total_up),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
            else:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_up,
                    label = " : up, {:.2e} in 45 degree, ".format(label_scale_up) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total up is {:.2e}".format(total_up),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
                
            plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
            plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

            if ifloglength:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "log/" + "up/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )
            else:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "up/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )

            # close figure after save
            plt.close('all')

        for indexstep, step in enumerate(inputstepsarray):
            x_value_plot = x_value_down[indexstep]
            y_value_plot = y_value_down[indexstep]
            total_down = np.sum(y_value_plot)
            # plot ave_z velocity across y
            fig, ax = plt.subplots()
            
            # title

            if x_scale_str is None:
                x_label_str = x_name_label_in_figure
            else:
                x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
            ax.set_xlabel(x_label_str)
            
            if y_scale_str is None:
                y_label_str = y_name_label_in_figure
            else:
                y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
            ax.set_ylabel(y_label_str)            
            # plot
            if ifloglength:
                length = (x_value_plot**2 + y_value_plot**2)**0.5
                masktolog = (length>=logvaluemin)
                length_log = np.copy(length)
                length_log[masktolog] = 1+logvaluemin + np.log10(length[masktolog]/logvaluemin)
                #breakpoint()
                x_value_plot_log = length_log*x_value_plot/length
                y_value_plot_log = length_log*y_value_plot/length
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot_log, y_value_plot_log,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_down,
                    cmap='hsv',
                )
            else:
                Q = ax.quiver(
                    Coord1/coord_scale, Coord2/coord_scale, x_value_plot, y_value_plot,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale_down,
                    cmap='hsv',
                )
                
            
            ax.set_xlabel(
                'y ({})'.format(coord_scale_str)
            )
            ax.set_ylabel(
                'z ({})'.format(coord_scale_str)
            )
            time = time_from_start_rotate(step)
            if ifloglength:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_down,
                    label = "down, logscale, : {:.2e} in 45 degree, ".format(label_scale_down) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total down is {:.2e}".format(total_down),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
            else:
                ax.quiverkey(
                    Q, 0.1, 0.95, label_scale_down,
                    label = " : down, {:.2e} in 45 degree, ".format(label_scale_down) + "\n " + x_name + ", " + y_name + " at t {:.2f} s".format(time) + ", total down is {:.2e}".format(total_down),
                    labelpos='E',
                    coordinates='figure', angle=45,
                )
                
            plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
            plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))

            if ifloglength:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "log/" + "down/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )
            else:
                fig.savefig(
                    "_".join([
                        diagram_path_add_nve + subfoldername + "down/" + 'step',
                        str(step),
                        ".png",
                    ]),
                    format="png",
                    bbox_inches=None,
                )

            # close figure after save
            plt.close('all')
    # plot streamplot
    for indexstep, step in enumerate(inputstepsarray):
        x_value_plot = x_value[indexstep]
        y_value_plot = y_value[indexstep]
        # plot ave_z velocity across y
        fig, ax = plt.subplots()
        
        # title

        if x_scale_str is None:
            x_label_str = x_name_label_in_figure
        else:
            x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
        ax.set_xlabel(x_label_str)
        
        if y_scale_str is None:
            y_label_str = y_name_label_in_figure
        else:
            y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
        ax.set_ylabel(y_label_str)
        # check if line0
        if ifstreamplot:
            # plot
            strm = ax.streamplot(
                Coord1[:,0]/coord_scale, Coord2[0,:]/coord_scale, np.transpose(x_value_plot), np.transpose(y_value_plot),
                linewidth=1, color='k',
                density=[0.8, 0.8],
            )
        vector_length = np.sqrt(x_value_plot**2 + y_value_plot**2)
        vector_length = np.ma.masked_where(np.logical_not(vector_length > 0), vector_length)
        #breakpoint()
        d_Coord1 = Coord1[:,0][1]-Coord1[:,0][0]
        d_Coord2 = Coord2[0,:][1]-Coord2[0,:][0]
        Coord1_expand = Coord1[:,0] - d_Coord1/2
        Coord1_expand = np.append(Coord1_expand, (Coord1_expand[-1]+d_Coord1))
        Coord2_expand = Coord2[0,:] - d_Coord2/2
        Coord2_expand = np.append(Coord2_expand, (Coord2_expand[-1]+d_Coord2))
        if ifcontour:
            if contour_v_min_max == "constant":
                if contour_norm == "linear":
                    contour = ax.pcolor(
                        Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                        norm=colors.Normalize(vmin=vmin, vmax=vmax), #norm=colors.LogNorm(vmin=10**-1, vmax=10), #norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                        cmap='coolwarm',
                    )
                elif contour_norm == "log":
                    contour = ax.pcolor(
                    Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                    norm=colors.LogNorm(vmin=vmin, vmax=vmax), #norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                    cmap='coolwarm',
                )
                else:
                    sys.exit("Error: contour_norm not defined")
            elif contour_v_min_max == "min_to_max":
                if contour_norm == "linear":
                    contour = ax.pcolor(
                        Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                        norm=colors.Normalize(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                        cmap='coolwarm',
                    )
                elif contour_norm == "log":
                    contour = ax.pcolor(
                    Coord1_expand/coord_scale, Coord2_expand/coord_scale, np.transpose(vector_length),
                    norm=colors.LogNorm(vmin=np.transpose(vector_length).min(), vmax=np.transpose(vector_length).max()),
                    cmap='coolwarm',
                )
                else:
                    sys.exit("Error: contour_norm not defined")
            else:
                sys.exit("Error: contour_v_min_max not defined")
            fig.colorbar(contour, ax=ax, extend='max')

        time = time_from_start_rotate(step)
        strain = strain_from_rotate_start(step)
        ax.set_title(
            "t = {:.2f} s".format(time) + ", strain is {:.2f}".format(strain)
        )
        ax.set_xlabel(
            'y ({})'.format(coord_scale_str)
        )
        ax.set_ylabel(
            'z ({})'.format(coord_scale_str)
        )

        plt.xticks(np.arange(0, (max(Coord1[:,0])-min(Coord1[:,0]))/coord_scale+2, 2))
        plt.yticks(np.arange(0, (max(Coord2[0,:])-min(Coord2[0,:]))/coord_scale+6, 6))
        fig.savefig(
            "_".join([
                diagram_path_add_nve + subfoldername + 'streamplot/' + 'step',
                str(step),
                ".png",
            ]),
            format="png",
            bbox_inches=None,
        )

        # close figure after save
        plt.close('all')


def plotdata_nochunk(
        lmp_path,
        if_on_paper, n_ave, 
        inputstepsarray,
        x_name, x_scale_factor, fig_x_label,
        y_name, y_scale_factor, fig_y_label,
        useerrorbar,
        if_mv_over_m=False,
        if_include_0_y_axis=False,
    ):
    # onpaper

    # plot data from chunk2D chunk1D
    # choose lmp_path path for read data
    # n_ave for average data output
    # x_name, y_name
    # choose timestep inputstepsarray
    # choose coord1_index_array, coord2_index_array,
    # n_ave_coord1, n_ave_coord2 for number of average
    # x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        # ifdivideNcount = False,
        # maskstatic=None,
        # masknonstatic=None,
    # output x y coord1 coord2 time
    
    # useerrorbar = True
    # output std x_value_std
    

    # figure_class=None, legend_class=None,
    # x_scale = 'linear', y_scale = 'linear',
    
    

    # output figure

    # create folder
    plotdata_createfolder(if_on_paper, n_ave, x_name, y_name)

    (fig, ax) = plotdata_figure_process_nochunk(
        lmp_path,
        if_on_paper,
        n_ave,
        inputstepsarray,
        x_name, x_scale_factor, fig_x_label,
        y_name, y_scale_factor, fig_y_label,
        useerrorbar,
        if_mv_over_m=if_mv_over_m,
        if_include_0_y_axis=if_include_0_y_axis,
    )
    file_name = "_".join([
        y_name,
        x_name,
        "nave_" + str(n_ave),
        "st_" + transfer_time_to_str(inputstepsarray),
    ])
    fig.savefig(
    "".join([
        plotdata_folder_path(if_on_paper, n_ave, x_name, y_name),
        file_name,
        ".png",
    ]),
    format="png",
    )
    # close figure after save
    plt.close('all')


def plotdata(
        lmp_path,
        if_on_paper, n_ave, 
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        x_name, x_scale_factor, fig_x_label, x_array_dim_order, x_array_legend_input_array_dim_order, x_array_legend_output_array_dim_order,
        y_name, y_scale_factor, fig_y_label, y_array_dim_order, y_array_legend_input_array_dim_order, y_array_legend_output_array_dim_order,
        useerrorbar,
        if_mv_over_m=False,
    ):
    # onpaper

    # plot data from chunk2D chunk1D
    # choose lmp_path path for read data
    # n_ave for average data output
    # x_name, y_name
    # choose timestep inputstepsarray
    # choose coord1_index_array, coord2_index_array,
    # n_ave_coord1, n_ave_coord2 for number of average
    # x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        # ifdivideNcount = False,
        # maskstatic=None,
        # masknonstatic=None,
    # output x y coord1 coord2 time
    
    # useerrorbar = True
    # output std x_value_std
    

    # figure_class=None, legend_class=None,
    # x_scale = 'linear', y_scale = 'linear',
    
    

    # output figure

    # create folder
    plotdata_createfolder(if_on_paper, n_ave, x_name, y_name)

    (fig, ax) = plotdata_figure_process(
        lmp_path,
        if_on_paper,
        n_ave,
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        x_name, x_scale_factor, fig_x_label, x_array_dim_order, x_array_legend_input_array_dim_order, x_array_legend_output_array_dim_order,
        y_name, y_scale_factor, fig_y_label, y_array_dim_order, y_array_legend_input_array_dim_order, y_array_legend_output_array_dim_order,
        useerrorbar,
        if_mv_over_m=if_mv_over_m,
    )
    plotdata_figure_save(
        fig, ax,
        n_ave,
        x_name,
        y_name,
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        y_array_legend_input_array_dim_order,
        y_array_legend_output_array_dim_order,
        if_on_paper,
    )

def plotdata1(
        lmp_path,
        if_on_paper, n_ave, 
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        x_name, x_scale_factor, fig_x_label, x_array_dim_order, x_array_legend_input_array_dim_order, x_array_legend_output_array_dim_order,
        y_name, y_scale_factor, fig_y_label, y_array_dim_order, y_array_legend_input_array_dim_order, y_array_legend_output_array_dim_order,
        useerrorbar,
        if_mv_over_m=False,
    ):
    y_dict = {
        'name': y_name, 'scale_factor': y_scale_factor, 'fig_label': fig_y_label,
        'array_dim_order': y_array_dim_order,
        'legend_input_array_dim_order': y_array_legend_input_array_dim_order, 'legend_output_array_dim_order': y_array_legend_output_array_dim_order,
    }
    # create folder
    plotdata_createfolder(if_on_paper, n_ave, x_name, y_name)

    (fig, ax) = plotdata_figure_process(
        lmp_path,
        if_on_paper,
        n_ave,
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        x_name, x_scale_factor, fig_x_label, x_array_dim_order, x_array_legend_input_array_dim_order, x_array_legend_output_array_dim_order,
        y_name, y_scale_factor, fig_y_label, y_array_dim_order, y_array_legend_input_array_dim_order, y_array_legend_output_array_dim_order,
        useerrorbar,
        if_mv_over_m=if_mv_over_m,
    )
    listforcolor = [0, 5, 10, 15, 20, 25, 30]
    colorlist = ['black', 'blue', 'green', 'red', 'c', 'purple', 'chocolate', 'pink', 'navy']
    for i in range(len(listforcolor)-1):
        for n in range(listforcolor[i], listforcolor[i+1]):
            ax.properties()['children'][n].set_color(colorlist[i])
    # legend x y
    legend_dic = legend_dic_produce(y_dict, inputstepsarray, coord1_index_array, coord2_index_array)
    if not if_on_paper:
        ax.legend(title=legend_dic['title'], loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xscale('log')
    plotdata_figure_save(
        fig, ax,
        n_ave,
        x_name,
        y_name,
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        y_array_legend_input_array_dim_order,
        y_array_legend_output_array_dim_order,
        if_on_paper,
    )

def plotdata2(
        lmp_path,
        if_on_paper, n_ave, 
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        x_name, x_scale_factor, fig_x_label, x_array_dim_order, x_array_legend_input_array_dim_order, x_array_legend_output_array_dim_order,
        y_name, y_scale_factor, fig_y_label, y_array_dim_order, y_array_legend_input_array_dim_order, y_array_legend_output_array_dim_order,
        useerrorbar,
        if_mv_over_m=False,
    ):
    y_dict = {
        'name': y_name, 'scale_factor': y_scale_factor, 'fig_label': fig_y_label,
        'array_dim_order': y_array_dim_order,
        'legend_input_array_dim_order': y_array_legend_input_array_dim_order, 'legend_output_array_dim_order': y_array_legend_output_array_dim_order,
    }
    # create folder
    plotdata_createfolder(if_on_paper, n_ave, x_name, y_name)

    (fig, ax) = plotdata_figure_process(
        lmp_path,
        if_on_paper,
        n_ave,
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        x_name, x_scale_factor, fig_x_label, x_array_dim_order, x_array_legend_input_array_dim_order, x_array_legend_output_array_dim_order,
        y_name, y_scale_factor, fig_y_label, y_array_dim_order, y_array_legend_input_array_dim_order, y_array_legend_output_array_dim_order,
        useerrorbar,
        if_mv_over_m=if_mv_over_m,
    )
    listforcolor = np.arange(0,6,1)
    colorlist = ['black', 'blue', 'green', 'red', 'c', 'purple', 'chocolate', 'pink', 'navy','black', 'blue', 'green', 'red', 'c', 'purple', 'chocolate', 'pink', 'navy']
    for i in range(len(listforcolor)-1):
        for n in range(listforcolor[i], listforcolor[i+1]):
            ax.properties()['children'][n].set_color(colorlist[i])
    # legend x y
    legend_dic = legend_dic_produce(y_dict, inputstepsarray, coord1_index_array, coord2_index_array)
    if not if_on_paper:
        ax.legend(title=legend_dic['title'], loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xscale('log')
    plotdata_figure_save(
        fig, ax,
        n_ave,
        x_name,
        y_name,
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        y_array_legend_input_array_dim_order,
        y_array_legend_output_array_dim_order,
        if_on_paper,
    )

def legend_dic_produce(y_dict, inputstepsarray, coord1_index_array, coord2_index_array):
    if coord1_index_array is not None and coord2_index_array is not None:
        if y_dict["legend_output_array_dim_order"] == '(c1c2)t':
            coord_key_list = ["Coord1", "Coord2"]
            legend_input_coord_array_dim_order = 'c1c2'
            legend_output_coord_array_dim_order = '(c1c2)'
        elif y_dict["legend_output_array_dim_order"] in ['tc1c2ave']:
            coord_key_list = ["Coord1", "Coord2"]
            legend_input_coord_array_dim_order = 'c1c2'
            legend_output_coord_array_dim_order = 'c1c2ave'

    elif coord1_index_array is not None and coord2_index_array is None:
        coord_key_list = ["Coord1"]
        legend_input_coord_array_dim_order = 'c1'
        if y_dict["legend_output_array_dim_order"] in ['c1t', 'c1c2avet']:
            legend_output_coord_array_dim_order = 'c1'
        elif y_dict["legend_output_array_dim_order"] == 'c1avet':
            legend_output_coord_array_dim_order = 'c1ave'

    elif coord1_index_array is None and coord2_index_array is not None:
        coord_key_list = ["Coord2"]
        legend_input_coord_array_dim_order = 'c2'
        if y_dict["legend_output_array_dim_order"] in ['c2t', 'c2c1avet']:
            legend_output_coord_array_dim_order = 'c2'
        elif y_dict["legend_output_array_dim_order"] == 'c2avet':
            legend_output_coord_array_dim_order = 'c2ave'
    
    char = pn.c_r_npyfilepath_coord_char[y_dict['name']]["coordinate_characteristic"]
    Coord = {}
    for key in coord_key_list:
        Coord[key] = pn.load_coord(lmp_path, char, key)
        # select corrd
        Coord[key] = select_coord(Coord[key], legend_input_coord_array_dim_order, coord1_index_array, coord2_index_array)
        Coord[key] = value_array_to_plot(Coord[key], legend_input_coord_array_dim_order, legend_output_coord_array_dim_order)
    
    legend_dic = {}
    # create legend label string
    if y_dict["legend_output_array_dim_order"] != 'tc1c2ave':

        if coord1_index_array is not None and coord2_index_array is not None:
            legend_dic['label_list'] = [str(c1)+ "," +str(c2) for (c1, c2) in zip(Coord["Coord1"], Coord["Coord2"])]
            legend_dic['title'] = plot_func_input["Coord1"]["axis_label"] + "," + plot_func_input["Coord2"]["axis_label"]

        elif coord1_index_array is not None and coord2_index_array is None:
            legend_dic['label_list'] = [str(c1) for c1 in Coord["Coord1"]]
            legend_dic['title'] = plot_func_input["Coord1"]["axis_label"]

        elif coord1_index_array is None and coord2_index_array is not None:
            legend_dic['label_list'] = [str(c2) for c2 in Coord["Coord2"]]
            legend_dic['title'] = plot_func_input["Coord2"]["axis_label"]
            
    elif y_dict["legend_output_array_dim_order"] == 'tc1c2ave':
        legend_dic['label_list'] = ["{:.2f}".format(s) for s in strain_from_rotate_start(inputstepsarray)]
        legend_dic['title'] = 'strain'
    return legend_dic

def plotdata_folder_path(if_on_paper, n_ave, x_name, y_name):
    # folderpath
    if if_on_paper:
        folderpath = "".join([
            dp.diagram_path, 
            "onpaper/",
            "n_ave_",
            str(n_ave),
            "/",
            y_name,
            "_",
            x_name,
            "/",
        ])
    else:
        folderpath = "".join([
            dp.diagram_path,
            "notonpaper/",
            "n_ave_",
            str(n_ave),
            "/",
            y_name,
            "_",
            x_name,
            "/",
        ])
    return folderpath

def plotdata_createfolder(if_on_paper, n_ave, x_name, y_name):
    # create folder
    os.makedirs(
        plotdata_folder_path(if_on_paper, n_ave, x_name, y_name),
        exist_ok=True,
    )

def plotdata_file_name(
        if_on_paper, n_ave, x_name, y_name, stepsarray, coord1array, coord2array,
        y_array_legend_input_array_dim_order,
        y_array_legend_output_array_dim_order,
    ):
    # file_name
    if coord1array is not None and coord2array is not None:
        file_name = "_".join([
            y_name,
            x_name,
            "nave_" + str(n_ave),
            "st_" + transfer_time_to_str(stepsarray),
            "c1_" + transfer_coor_to_str(coord1array),
            "c2_" + transfer_coor_to_str(coord2array),
            "yli_" + y_array_legend_input_array_dim_order,
            "ylo_" + y_array_legend_output_array_dim_order,
        ])
    elif coord1array is None and coord2array is not None:
        file_name = "_".join([
            y_name,
            x_name,
            "nave_" + str(n_ave),
            "st_" + transfer_time_to_str(stepsarray),
            "c2_" + transfer_coor_to_str(coord2array),
            "yli_" + y_array_legend_input_array_dim_order,
            "ylo_" + y_array_legend_output_array_dim_order,
        ])
    elif coord1array is not None and coord2array is None:
        file_name = "_".join([
            y_name,
            x_name,
            "nave_" + str(n_ave),
            "st_" + transfer_time_to_str(stepsarray),
            "c1_" + transfer_coor_to_str(coord1array),
            "yli_" + y_array_legend_input_array_dim_order,
            "ylo_" + y_array_legend_output_array_dim_order,
        ])
    return file_name

def value_array_to_plot(array, input_array_dim_order, output_array_dim_order, motion=None):
    """
    reorganize array for easily plot
    ----------
    array : array_like
            Input array.
    input_array_dim_order : string
        'tc1c2' or other permutation for t, c1, c2.
    output_array_dim_order : string
        'tc1c2' or other permutation for t, c1, c2.
    motion : string
        all nochange separate choose1 sum ave
    Returns
    ----------
    p : ndarray
        `a` with its axes permuted.  A view is returned whenever
        possible.
    """
    if input_array_dim_order == "tc1c2" and output_array_dim_order == "c1c2avet" and motion == None:
        array = np.transpose(array, axes=[1, 2, 0])
        array = np.nanmean(array, axis=1)
    if input_array_dim_order == "tc1c2" and output_array_dim_order == "tc1c2ave" and motion == None:
        array = np.nanmean(array, axis=2)
    if input_array_dim_order == "c1c2" and output_array_dim_order == "c1c2ave" and motion == None:
        array = np.nanmean(array, axis=1)

    if input_array_dim_order == "tc1c2" and output_array_dim_order == "(c1c2)t" and motion == None:
        # order tc1c2 -> (c1c2)t
        # order tc1c2 -> c1c2t
        array = np.transpose(array, axes=[1, 2, 0])
        # order tc1c2 -> c1c2t
        # flatten all but the last dimension
        array = array.reshape(-1, array.shape[-1])
    
    if input_array_dim_order == "c1c2" and output_array_dim_order == "(c1c2)" and motion == None:
        # order c1c2 -> (c1c2)
        # flatten all
        array = array.reshape(-1)

    if input_array_dim_order == "t" and output_array_dim_order == "t" and motion == None:
        # order c1c2 -> (c1c2)
        # flatten all
        pass

    if input_array_dim_order == "tc1" and output_array_dim_order == "c1t" and motion == None:
        # order tc1 -> c1t
        array = np.transpose(array, axes=[1, 0])
    if input_array_dim_order == "c1" and output_array_dim_order == "c1" and motion == None:
        pass

    if input_array_dim_order == "c2" and output_array_dim_order == "c2" and motion == None:
        pass

    if input_array_dim_order == "tc2" and output_array_dim_order == "c2t" and motion == None:
        # order tc1 -> c1t
        array = np.transpose(array, axes=[1, 0])
    
    if input_array_dim_order == "c2" and output_array_dim_order == "c2ave" and motion == None:
        array = np.nanmean(array)
        array = np.expand_dims(array, axis=0)

    if input_array_dim_order == "tc2" and output_array_dim_order == "c2avet" and motion == None:
        # order tc1 -> c1t
        array = np.transpose(array, axes=[1, 0])
        array = np.nanmean(array, axis=0)
        array = np.expand_dims(array, axis=0)
    
    if input_array_dim_order == "c1" and output_array_dim_order == "c1ave" and motion == None:
        array = np.nanmean(array)
        array = np.expand_dims(array, axis=0)

    if input_array_dim_order == "tc1" and output_array_dim_order == "c1avet" and motion == None:
        # order tc1 -> c1t
        array = np.transpose(array, axes=[1, 0])
        array = np.nanmean(array, axis=0)
        array = np.expand_dims(array, axis=0)

    return array


# motion: all Nochange separate choose1 sum ave
# legend tc1c2: N s s
# curve tc1c2: sf N N
def reorganize_data_stress_strain(
        array
    ):
    # order tc1c2 -> (c1c2)t
    array = value_array_to_plot(array, "tc1c2", "(c1c2)t")

    # get coord
    char = pn.c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
    Coord = {}
    for key in ["Coord1", "Coord2"]:
        Coord[key] = pn.load_coord(lmp_path, char, key)
        Coord[key] = value_array_to_plot(Coord[key], "c1c2", "(c1c2)")
    # create legend label string
    stringlabel = {}
    for key in ["Coord1", "Coord2"]:
        stringlabel[key] = plot_func_input[key]["axis_label"]
    legend_label_string_list = [
        ",".join(list(c)) for c in zip(Coord["Coord1"], Coord["Coord2"])
    ]
    legend_label_string = stringlabel["Coord1"] + "," + stringlabel["Coord2"] + "=" + ",".join(legend_label_string_list)
    
    return (x_value_plot, y_value_plot, x_value_std_plot, y_value_std_plot, labelstring)

def select_coord(array, array_dim_order, coord1_index_array, coord2_index_array):
    if coord1_index_array is not None and coord2_index_array is not None: 
        if array_dim_order == 'tc1c2':
            if coord1_index_array == 'all':
                pass
            else:
                array = array[:, coord1_index_array, :]
            if coord2_index_array == 'all':
                pass
            else:
                array = array[:, :, coord2_index_array]
        elif array_dim_order == 'c1c2':
            if coord1_index_array == 'all':
                pass
            else:
                array = array[coord1_index_array, :]
            if coord2_index_array == 'all':
                pass
            else:
                array = array[:, coord2_index_array]
        elif array_dim_order == 't':
            pass
    elif coord1_index_array is not None and coord2_index_array is None:
        if array_dim_order == 'tc1':
            if coord1_index_array == 'all':
                pass
            else:
                array = array[:, coord1_index_array]
        elif array_dim_order == 'c1':
            if coord1_index_array == 'all':
                pass
            else:
                array = array[coord1_index_array]
        elif array_dim_order == 't':
            pass
    elif coord1_index_array is None and coord2_index_array is not None:
        if array_dim_order == 'tc2':
            if coord2_index_array == 'all':
                pass
            else:
                array = array[:, coord2_index_array]
        elif array_dim_order == 'c2':
            if coord2_index_array == 'all':
                pass
            else:
                array = array[coord2_index_array]
        elif array_dim_order == 't':
            pass

    return array



def plot_from_organized_data(
        fig, ax,
        if_on_paper,
        x_dict,
        y_dict,
        legend_label_list,
        legend_title,
        useerrorbar,
    ):

    if x_dict['legend_output_array_dim_order'] == 't' and y_dict['legend_output_array_dim_order'] in ['(c1c2)t', 'c1t', 'c2t', 'c1avet', 'c2avet']:
        # ax plot errorbar
        if useerrorbar:
            for i, legend_label in enumerate(legend_label_list):
                ax.errorbar(
                    x_dict['value'], y_dict['value'][i], xerr=x_dict['value_std'][i], yerr=y_dict['value_std'][i],
                    label=legend_label,
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                )
        else:
            for i, legend_label in enumerate(legend_label_list):
                ax.plot(
                    x_dict['value'], y_dict['value'][i],
                    label=legend_label,
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                )
    elif x_dict['legend_output_array_dim_order'] == 'c1c2ave' and y_dict['legend_output_array_dim_order'] in ['tc1c2ave']:
        # ax plot errorbar
        if useerrorbar:
            for i, legend_label in enumerate(legend_label_list):
                ax.errorbar(
                    x_dict['value'], y_dict['value'][i], xerr=x_dict['value_std'][i], yerr=y_dict['value_std'][i],
                    label=legend_label,
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                )
        else:
            for i, legend_label in enumerate(legend_label_list):
                ax.plot(
                    x_dict['value'], y_dict['value'][i],
                    label=legend_label,
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                )
    elif x_dict['legend_output_array_dim_order'] in ['(c1c2)t'] and y_dict['legend_output_array_dim_order'] in ['(c1c2)t']:
        # ax plot errorbar
        if useerrorbar:
            for i, legend_label in enumerate(legend_label_list):
                ax.errorbar(
                    x_dict['value'][i], y_dict['value'][i], xerr=x_dict['value_std'][i], yerr=y_dict['value_std'][i],
                    label=legend_label,
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                )
        else:
            for i, legend_label in enumerate(legend_label_list):
                ax.plot(
                    x_dict['value'][i], y_dict['value'][i],
                    label=legend_label,
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                )
    # legend x y
    if not if_on_paper:
        ax.legend(title=legend_title, loc='upper left', bbox_to_anchor=(1, 1))

    return (fig, ax)

def plotdata_figure_process_nochunk(
    lmp_path,
    if_on_paper,
    n_ave,
    inputstepsarray,
    x_name, x_scale_factor, fig_x_label,
    y_name, y_scale_factor, fig_y_label,
    useerrorbar,
    if_mv_over_m=False,
    if_include_0_y_axis=False,
    ):
    # set x y dict
    x_dict = {
        'name': x_name, 'scale_factor': x_scale_factor, 'fig_label': fig_x_label,
        'value': None, 'value_std': None,
    }
    y_dict = {
        'name': y_name, 'scale_factor': y_scale_factor, 'fig_label': fig_y_label,
        'value': None, 'value_std': None,
    }
    # get value for x y and rescale
    for d in [x_dict, y_dict]:
        # get value
        d['value'] = get_value_nochunk_include_ts_t_st(lmp_path, d['name'], n_ave, inputstepsarray, y_name)/d['scale_factor']
        
    # check if x-array y-array has the same shape
    xarray = x_dict['value']
    yarray = y_dict['value']
    if xarray.shape[-1] != yarray.shape[-1]:
        sys.exit("shape for x y are different")
    
    
    
    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    
    # xy scale label
    ax.set_xlabel(x_dict['fig_label'])
    ax.set_ylabel(y_dict['fig_label'])
    # rotate label
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    ax.plot(
        x_dict['value'], y_dict['value'],
        marker = ".",
        linestyle = 'None',
        markersize=12,
    )
    # legend x y
    if not if_on_paper:
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    # if_include_0_y_axis
    if if_include_0_y_axis:
        if ax.get_ylim()[0]*ax.get_ylim()[1] <= 0:
            pass
        else:
            if ax.get_ylim()[0] > 0:
                ax.set_ylim(bottom=0)
            else:
                ax.set_ylim(top=0)

    return (fig, ax)


# x-dict: name scale_factor fig_label
def plotdata_figure_process(
    lmp_path,
    if_on_paper,
    n_ave,
    inputstepsarray,
    coord1_index_array,
    coord2_index_array,
    x_name, x_scale_factor, fig_x_label, x_array_dim_order, x_array_legend_input_array_dim_order, x_array_legend_output_array_dim_order,
    y_name, y_scale_factor, fig_y_label, y_array_dim_order, y_array_legend_input_array_dim_order, y_array_legend_output_array_dim_order,
    useerrorbar,
    if_mv_over_m=False,
    ):
    # set x y dict
    x_dict = {
        'name': x_name, 'scale_factor': x_scale_factor, 'fig_label': fig_x_label,
        'array_dim_order': x_array_dim_order,
        'legend_input_array_dim_order': x_array_legend_input_array_dim_order, 'legend_output_array_dim_order': x_array_legend_output_array_dim_order,
        'value': None, 'value_std': None,
    }
    y_dict = {
        'name': y_name, 'scale_factor': y_scale_factor, 'fig_label': fig_y_label,
        'array_dim_order': y_array_dim_order,
        'legend_input_array_dim_order': y_array_legend_input_array_dim_order, 'legend_output_array_dim_order': y_array_legend_output_array_dim_order,
        'value': None, 'value_std': None,
    }
    if if_mv_over_m:
        if y_dict['name'] == "velocity_1":
            mv_name = 'mv_1'
        elif y_dict['name'] == "velocity_2":
            mv_name = 'mv_2'
        elif y_dict['name'] == "velocity_3":
            mv_name = 'mv_3'

        mv_dict = {
            'name': mv_name, 'scale_factor': y_scale_factor, 'fig_label': fig_y_label,
            'array_dim_order': y_array_dim_order,
            'legend_input_array_dim_order': y_array_legend_input_array_dim_order, 'legend_output_array_dim_order': y_array_legend_output_array_dim_order,
            'value': None, 'value_std': None,
        }
        m_dict = {
            'name': 'mass', 'scale_factor': y_scale_factor, 'fig_label': fig_y_label,
            'array_dim_order': y_array_dim_order,
            'legend_input_array_dim_order': y_array_legend_input_array_dim_order, 'legend_output_array_dim_order': y_array_legend_output_array_dim_order,
            'value': None, 'value_std': None,
        }

        for d in [mv_dict, m_dict]:
            # get value
            d['value'] = get_value_include_ts_t_st(lmp_path, d['name'], n_ave, inputstepsarray, y_name)
            # select corrd
            d['value'] = select_coord(d['value'], d['array_dim_order'], coord1_index_array, coord2_index_array)
            # reorganize data set for plot convenience
            d['value'] = value_array_to_plot(d['value'], d['legend_input_array_dim_order'], d["legend_output_array_dim_order"])
            
            # if need errorbar then get std value
            if useerrorbar:
                d['value_std'] = get_std_value_include_ts_t_st(lmp_path, d['name'] + "_std", n_ave, inputstepsarray, y_name)
                # select corrd
                d['value_std'] = select_coord(d['value_std'], d['array_dim_order'], coord1_index_array, coord2_index_array)
                # reorganize data set for plot convenience
                d['value_std'] = value_array_to_plot(d['value_std'], d['legend_input_array_dim_order'], d["legend_output_array_dim_order"])
        
        y_dict['value'] = mv_dict['value']/m_dict['value']/d['scale_factor']
        #y_dict['value_std'] = mv_dict['value_std']/m_dict['value_std']
        
        # get value for x y and rescale
        for d in [x_dict]:
            # get value
            d['value'] = get_value_include_ts_t_st(lmp_path, d['name'], n_ave, inputstepsarray, y_name)/d['scale_factor']
            # select corrd
            d['value'] = select_coord(d['value'], d['array_dim_order'], coord1_index_array, coord2_index_array)
            # reorganize data set for plot convenience
            d['value'] = value_array_to_plot(d['value'], d['legend_input_array_dim_order'], d["legend_output_array_dim_order"])
            
            # if need errorbar then get std value
            if useerrorbar:
                d['value_std'] = get_std_value_include_ts_t_st(lmp_path, d['name'] + "_std", n_ave, inputstepsarray, y_name)/d['scale_factor']
                # select corrd
                d['value_std'] = select_coord(d['value_std'], d['array_dim_order'], coord1_index_array, coord2_index_array)
                # reorganize data set for plot convenience
                d['value_std'] = value_array_to_plot(d['value_std'], d['legend_input_array_dim_order'], d["legend_output_array_dim_order"])
    else:
        # get value for x y and rescale
        for d in [x_dict, y_dict]:
            # get value
            d['value'] = get_value_include_ts_t_st(lmp_path, d['name'], n_ave, inputstepsarray, y_name)/d['scale_factor']
            # select corrd
            d['value'] = select_coord(d['value'], d['array_dim_order'], coord1_index_array, coord2_index_array)
            # reorganize data set for plot convenience
            d['value'] = value_array_to_plot(d['value'], d['legend_input_array_dim_order'], d["legend_output_array_dim_order"])
            
            # if need errorbar then get std value
            if useerrorbar:
                d['value_std'] = get_std_value_include_ts_t_st(lmp_path, d['name'] + "_std", n_ave, inputstepsarray, y_name)/d['scale_factor']
                # select corrd
                d['value_std'] = select_coord(d['value_std'], d['array_dim_order'], coord1_index_array, coord2_index_array)
                # reorganize data set for plot convenience
                d['value_std'] = value_array_to_plot(d['value_std'], d['legend_input_array_dim_order'], d["legend_output_array_dim_order"])
    
    # check if x-array y-array has the same shape
    xarray = x_dict['value']
    yarray = y_dict['value']
    if xarray.shape[-1] != yarray.shape[-1]:
        sys.exit("shape for x y are different")
    
    legend_dic = legend_dic_produce(y_dict, inputstepsarray, coord1_index_array, coord2_index_array)
    legend_label_list = legend_dic['label_list']
    legend_title = legend_dic['title']
    
    # plot ave_z velocity across y
    fig, ax = plt.subplots()
    
    # xy scale label
    ax.set_xlabel(x_dict['fig_label'])
    ax.set_ylabel(y_dict['fig_label'])
    # rotate label
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    (fig, ax) = plot_from_organized_data(
        fig, ax,
        if_on_paper,
        x_dict,
        y_dict,
        legend_label_list,
        legend_title,
        useerrorbar,
    )
    
    return (fig, ax)

# save figure file
def plotdata_figure_save(
        fig, ax,
        n_ave,
        x_name,
        y_name,
        inputstepsarray,
        coord1_index_array,
        coord2_index_array,
        y_array_legend_input_array_dim_order,
        y_array_legend_output_array_dim_order,
        if_on_paper,
    ):
    file_name = plotdata_file_name(
        if_on_paper, n_ave, x_name, y_name, inputstepsarray, coord1_index_array, coord2_index_array,
        y_array_legend_input_array_dim_order,
        y_array_legend_output_array_dim_order,
    )
    fig.savefig(
    "".join([
        plotdata_folder_path(if_on_paper, n_ave, x_name, y_name),
        file_name,
        ".png",
    ]),
    format="png",
    )
    # close figure after save
    plt.close('all')

# collection of inputs of variable in plot function
# v_name -> ifonpaper
# need output for figure name
plot_func_input = {
    "pressure": {"axis_label": "P",},
    "n_contact": {"axis_label": "Z",},
    "velocity_1": {"axis_label": "V",},
    "velocity_2": {"axis_label": "V",},
    "velocity_3": {"axis_label": "V",},
    "strain": {"axis_label": r'$\gamma$',},
    "inwall_force_1": {"axis_label": r'$F_x$',},
    "inwall_force_2": {"axis_label": r'$F_y$',},
    "inwall_force_3": {"axis_label": r'$F_z$',},
    "outwall_force_1": {"axis_label": r'$F_x$',},
    "outwall_force_2": {"axis_label": r'$F_y$',},
    "outwall_force_3": {"axis_label": r'$F_z$',},
    "zbottom_force_1": {"axis_label": r'$F_x$',},
    "zbottom_force_2": {"axis_label": r'$F_y$',},
    "zbottom_force_3": {"axis_label": r'$F_z$',},
    "strain_rate_21": {"axis_label": r'$\dot{\gamma_{21}}$',},
    "strain_rate_22": {"axis_label": r'$\dot{\gamma_{22}}$',},
    "strain_rate_23": {"axis_label": r'$\dot{\gamma_{23}}$',},
    "strain_rate_31": {"axis_label": r'$\dot{\gamma_{31}}$',},
    "strain_rate_32": {"axis_label": r'$\dot{\gamma_{32}}$',},
    "strain_rate_33": {"axis_label": r'$\dot{\gamma_{33}}$',},
    "strain_rate_21_middle": {"axis_label": r'$\dot{\gamma_{21}}$',},
    "strain_rate_22_middle": {"axis_label": r'$\dot{\gamma_{22}}$',},
    "strain_rate_23_middle": {"axis_label": r'$\dot{\gamma_{23}}$',},
    "strain_rate_31_middle": {"axis_label": r'$\dot{\gamma_{31}}$',},
    "strain_rate_32_middle": {"axis_label": r'$\dot{\gamma_{32}}$',},
    "strain_rate_33_middle": {"axis_label": r'$\dot{\gamma_{33}}$',},
    "stress_11": {"axis_label": r'$\sigma_{11}$',},
    "stress_22": {"axis_label": r'$\sigma_{22}$',},
    "stress_33": {"axis_label": r'$\sigma_{33}$',},
    "stress_12": {"axis_label": r'$\sigma_{12}$',},
    "stress_13": {"axis_label": r'$\sigma_{13}$',},
    "stress_23": {"axis_label": r'$\sigma_{23}$',},
    "mu_12": {"axis_label": r'$\mu_{12}$',},
    "mu_13": {"axis_label": r'$\mu_{13}$',},
    "mu_23": {"axis_label": r'$\mu_{23}$',},
    "mu_12_middle": {"axis_label": r'$\mu_{12}$',},
    "mu_13_middle": {"axis_label": r'$\mu_{13}$',},
    "mu_23_middle": {"axis_label": r'$\mu_{23}$',},
    "mu_tensor_12": {"axis_label": r'$\mu_{12}$',},
    "mu_tensor_13": {"axis_label": r'$\mu_{13}$',},
    "mu_tensor_23": {"axis_label": r'$\mu_{23}$',},
    "I_12": {"axis_label": r'$I_{12}$',},
    "I_13": {"axis_label": r'$I_{13}$',},
    "I_23": {"axis_label": r'$I_{23}$',},
    "I_tensor": {"axis_label": r'$I$',},
    "fraction": {"axis_label": r'$\phi$',},
    "inwall_stress_1": {"axis_label": r'$\sigma_{21}$',},
    "inwall_stress_2": {"axis_label": r'$\sigma_{22}$',},
    "inwall_stress_3": {"axis_label": r'$\sigma_{23}$',},
    "outwall_stress_1": {"axis_label": r'$\sigma_{21}$',},
    "outwall_stress_2": {"axis_label": r'$\sigma_{22}$',},
    "outwall_stress_3": {"axis_label": r'$\sigma_{23}$',},
    "zbottom_stress_1": {"axis_label": r'$\sigma_{31}$',},
    "zbottom_stress_2": {"axis_label": r'$\sigma_{32}$',},
    "zbottom_stress_3": {"axis_label": r'$\sigma_{33}$',},
    "Coord1": {"axis_label": r'y',},
    "Coord2": {"axis_label": r'z',},
}

def main():
    # set plot path 
    lmp_path = rr.folder_path_list_initial_to_last[-1]

    # plot nochunk
    plot_input_dict_list_nochunk = [
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_outwall_1', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{21}$' + "(static wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_outwall_2', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{22}$' + "(static wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_outwall_3', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{23}$' + "(static wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_inwall_1', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{21}$' + "(shearing wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_inwall_2', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{22}$' + "(shearing wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_inwall_3', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{23}$' + "(shearing wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_zbottom_1', "y_scale_factor": bottom_area*stress_scale_height, 'fig_y_label': r'$\sigma_{31}$' + "(bottom)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_zbottom_2', "y_scale_factor": bottom_area*stress_scale_height, 'fig_y_label': r'$\sigma_{32}$' + "(bottom)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5600000, 40000000, 500000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_zbottom_3', "y_scale_factor": bottom_area*stress_scale_height, 'fig_y_label': r'$\sigma_{33}$' + "(bottom)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
    ]
    plot_input_dict_list_nochunk_static = [
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_outwall_1', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{21}$' + "(static wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_outwall_2', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{22}$' + "(static wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_outwall_3', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{23}$' + "(static wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_inwall_1', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{21}$' + "(shearing wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_inwall_2', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{22}$' + "(shearing wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_inwall_3', "y_scale_factor": inwall_area*stress_scale_height, 'fig_y_label': r'$\sigma_{23}$' + "(shearing wall)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_zbottom_1', "y_scale_factor": bottom_area*stress_scale_height, 'fig_y_label': r'$\sigma_{31}$' + "(bottom)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_zbottom_2', "y_scale_factor": bottom_area*stress_scale_height, 'fig_y_label': r'$\sigma_{32}$' + "(bottom)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 1, 
            'inputstepsarray': np.arange(4800000, 5600000, 100000),
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$',
            'y_name': 'v_force_zbottom_3', "y_scale_factor": bottom_area*stress_scale_height, 'fig_y_label': r'$\sigma_{33}$' + "(bottom)",
            'useerrorbar': False,
            'if_include_0_y_axis': True,
        },
    ]
    for plot_input_dict in plot_input_dict_list_nochunk:
        plotdata_nochunk(**plot_input_dict)
    for plot_input_dict in plot_input_dict_list_nochunk_static:
        plotdata_nochunk(**plot_input_dict)
    # plot velocity-Coord1
    plot_input_dict_list_velocity = [
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.append(np.arange(4500000, 7000000, 500000), np.array([40000000])),
            'coord1_index_array': 'all',
            'coord2_index_array': 'all',
            'x_name': "Coord1", 'x_scale_factor': scale['Coord'], 'fig_x_label': 'y', 'x_array_dim_order': 'c1c2',
            'x_array_legend_input_array_dim_order': 'c1c2', 'x_array_legend_output_array_dim_order': 'c1c2ave',
            'y_name': 'velocity_1', "y_scale_factor": scale['velocity'], 'fig_y_label': 'V', 'y_array_dim_order': 'tc1c2',
            'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': 'tc1c2ave',
            'useerrorbar': False,
            'if_mv_over_m': True,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.append(np.arange(4500000, 7000000, 500000), np.array([40000000])),
            'coord1_index_array': 'all',
            'coord2_index_array': 'all',
            'x_name': "Coord1", 'x_scale_factor': scale['Coord'], 'fig_x_label': 'y', 'x_array_dim_order': 'c1c2',
            'x_array_legend_input_array_dim_order': 'c1c2', 'x_array_legend_output_array_dim_order': 'c1c2ave',
            'y_name': 'velocity_1', "y_scale_factor": scale['velocity'], 'fig_y_label': 'V', 'y_array_dim_order': 'tc1c2',
            'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': 'tc1c2ave',
            'useerrorbar': False,
            'if_mv_over_m': True,
        },
    ]
    
    for plot_input_dict in plot_input_dict_list_velocity:
        plotdata(**plot_input_dict)

    """
    # plot
    plot_input_dict_list = []
    plot_input_dict_list_wall = []
    
    plot_input_dict_list_inwall = [
        # wall stress average
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'inwall_stress_1', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{21}$', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'inwall_stress_1', "y_scale_factor": scale['stress'], 'fig_y_label': 'inwall_stress_1', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'inwall_stress_2', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{22}$', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'inwall_stress_2', "y_scale_factor": scale['stress'], 'fig_y_label': 'inwall_stress_2', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'inwall_stress_3', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{23}$', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'inwall_stress_3', "y_scale_factor": scale['stress'], 'fig_y_label': 'inwall_stress_3', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
    ]
    plot_input_dict_list_wall = plot_input_dict_list_wall + plot_input_dict_list_inwall
    
    plot_input_dict_list_outwall = [
        # wall stress average
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'outwall_stress_1', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{21}$', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'outwall_stress_1', "y_scale_factor": scale['stress'], 'fig_y_label': 'outwall_stress_1', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'outwall_stress_2', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{22}$', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'outwall_stress_2', "y_scale_factor": scale['stress'], 'fig_y_label': 'outwall_stress_2', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'outwall_stress_3', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{23}$', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': None,
            'coord2_index_array': 'all',
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'outwall_stress_3', "y_scale_factor": scale['stress'], 'fig_y_label': 'outwall_stress_3', 'y_array_dim_order': 'tc2',
            'y_array_legend_input_array_dim_order': 'tc2', 'y_array_legend_output_array_dim_order': 'c2avet',
            'useerrorbar': False,
        },
    ]
    plot_input_dict_list_wall = plot_input_dict_list_wall + plot_input_dict_list_outwall
    
    plot_input_dict_list_zbottom = [
        # wall stress average
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': 'all',
            'coord2_index_array': None,
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'zbottom_stress_1', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{31}$', 'y_array_dim_order': 'tc1',
            'y_array_legend_input_array_dim_order': 'tc1', 'y_array_legend_output_array_dim_order': 'c1avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': 'all',
            'coord2_index_array': None,
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'zbottom_stress_1', "y_scale_factor": scale['stress'], 'fig_y_label': 'zbottom_stress_1', 'y_array_dim_order': 'tc1',
            'y_array_legend_input_array_dim_order': 'tc1', 'y_array_legend_output_array_dim_order': 'c1avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': 'all',
            'coord2_index_array': None,
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'zbottom_stress_2', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{32}$', 'y_array_dim_order': 'tc1',
            'y_array_legend_input_array_dim_order': 'tc1', 'y_array_legend_output_array_dim_order': 'c1avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': 'all',
            'coord2_index_array': None,
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'zbottom_stress_2', "y_scale_factor": scale['stress'], 'fig_y_label': 'zbottom_stress_2', 'y_array_dim_order': 'tc1',
            'y_array_legend_input_array_dim_order': 'tc1', 'y_array_legend_output_array_dim_order': 'c1avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': True, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': 'all',
            'coord2_index_array': None,
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'zbottom_stress_3', "y_scale_factor": scale['stress'], 'fig_y_label': r'$\sigma_{33}$', 'y_array_dim_order': 'tc1',
            'y_array_legend_input_array_dim_order': 'tc1', 'y_array_legend_output_array_dim_order': 'c1avet',
            'useerrorbar': False,
        },
        {
            'lmp_path': lmp_path,
            'if_on_paper': False, 'n_ave': 51, 
            'inputstepsarray': np.arange(5100000, 40000000, 500000),
            'coord1_index_array': 'all',
            'coord2_index_array': None,
            'x_name': "strain", 'x_scale_factor': 1, 'fig_x_label': r'$\gamma$', 'x_array_dim_order': 't',
            'x_array_legend_input_array_dim_order': 't', 'x_array_legend_output_array_dim_order': 't',
            'y_name': 'zbottom_stress_3', "y_scale_factor": scale['stress'], 'fig_y_label': 'zbottom_stress_3', 'y_array_dim_order': 'tc1',
            'y_array_legend_input_array_dim_order': 'tc1', 'y_array_legend_output_array_dim_order': 'c1avet',
            'useerrorbar': False,
        },
    ]
    plot_input_dict_list_wall = plot_input_dict_list_wall + plot_input_dict_list_zbottom
    
    plot_input_dict_list = plot_input_dict_list + plot_input_dict_list_wall
    
    for plot_input_dict in plot_input_dict_list:
        plotdata(**plot_input_dict)
    """
    """
    

    # mu-I
    ## mu-I steady state
    for (coord1_index_array, coord2_index_array) in [('all', np.arange(15)), ('all',[0]),([0],np.arange(15)),([-1],np.arange(15))]:
        plot_input_dict_list_mu_I = [
            {
                'lmp_path': lmp_path,
                'if_on_paper': True, 'n_ave': 201, 
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_12", 'x_scale_factor': scale['I'], 'fig_x_label': r'$I_{12}$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_12_middle', "y_scale_factor": scale['mu'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
            {
                'lmp_path': lmp_path,
                'if_on_paper': False, 'n_ave': 201,
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_12", 'x_scale_factor': scale['I'], 'fig_x_label': r'$I_{12}$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_12_middle', "y_scale_factor": scale['mu'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
            {
                'lmp_path': lmp_path,
                'if_on_paper': True, 'n_ave': 201, 
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_tensor", 'x_scale_factor': scale['I_tensor'], 'fig_x_label': r'$I$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_tensor_12', "y_scale_factor": scale['mu_tensor'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
            {
                'lmp_path': lmp_path,
                'if_on_paper': False, 'n_ave': 201, 
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_tensor", 'x_scale_factor': scale['I_tensor'], 'fig_x_label': r'$I$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_tensor_12', "y_scale_factor": scale['mu_tensor'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
        ]
        for plot_input_dict in plot_input_dict_list_mu_I:
            plotdata(**plot_input_dict)
    
    for (coord1_index_array, coord2_index_array) in [([0,-1], np.arange(15))]:
        plot_input_dict_list_mu_I = [
            {
                'lmp_path': lmp_path,
                'if_on_paper': True, 'n_ave': 201, 
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_12", 'x_scale_factor': scale['I'], 'fig_x_label': r'$I_{12}$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_12_middle', "y_scale_factor": scale['mu'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
            {
                'lmp_path': lmp_path,
                'if_on_paper': False, 'n_ave': 201,
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_12", 'x_scale_factor': scale['I'], 'fig_x_label': r'$I_{12}$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_12_middle', "y_scale_factor": scale['mu'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
        ]
        for plot_input_dict in plot_input_dict_list_mu_I:
            plotdata1(**plot_input_dict)

    for (coord1_index_array, coord2_index_array) in [([6], np.arange(9,15))]:
        plot_input_dict_list_mu_I = [
            {
                'lmp_path': lmp_path,
                'if_on_paper': True, 'n_ave': 201, 
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_12", 'x_scale_factor': scale['I'], 'fig_x_label': r'$I_{12}$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_12_middle', "y_scale_factor": scale['mu'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
            {
                'lmp_path': lmp_path,
                'if_on_paper': False, 'n_ave': 201,
                'inputstepsarray': np.array([40000000]),
                'coord1_index_array': coord1_index_array,
                'coord2_index_array': coord2_index_array,
                'x_name': "I_12", 'x_scale_factor': scale['I'], 'fig_x_label': r'$I_{12}$', 'x_array_dim_order': 'tc1c2',
                'x_array_legend_input_array_dim_order': 'tc1c2', 'x_array_legend_output_array_dim_order': '(c1c2)t',
                'y_name': 'mu_12_middle', "y_scale_factor": scale['mu'], 'fig_y_label': r'$\mu_{12}$', 'y_array_dim_order': 'tc1c2',
                'y_array_legend_input_array_dim_order': 'tc1c2', 'y_array_legend_output_array_dim_order': '(c1c2)t',
                'useerrorbar': False,
                'if_mv_over_m': False,
            },
        ]
        for plot_input_dict in plot_input_dict_list_mu_I:
            plotdata2(**plot_input_dict)

    
    plot_quiver_from_chunk2D(
        lmp_path,
        201, "mu_12", "mu_12", np.arange(10000000,75000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=0.6,
    )

    plot_quiver_from_chunk2D(
        lmp_path,
        201, "I_12", "I_12", np.arange(10000000,75000000,5000000),
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=2*10**-3,
    )

    plot_quiver_from_chunk2D(
        lmp_path,
        201, "velocity_2", "velocity_3", np.arange(10000000,75000000,5000000),
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=0.0001, label_scale=0.0001,
        ifloglength=True,
        ifplotseparateupdown=True,
        quiver_scale_up=0.0001, label_scale_up=0.0001,
        quiver_scale_down=0.0005, label_scale_down=0.0005,
    )

    plot_quiver_from_chunk2D(
        lmp_path,
        201, "fraction", "fraction", np.arange(10000000,75000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0.4, vmax=0.7,
    )
    
    plot_quiver_from_chunk2D(
        lmp_path,
        51, "mu_12", "mu_12", np.arange(5260000,10000000,500000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=0.6,
    )

    plot_quiver_from_chunk2D(
        lmp_path,
        51, "I_12", "I_12", np.arange(5260000,10000000,500000),
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=2*10**-3,
    )

    plot_quiver_from_chunk2D(
        lmp_path,
        51, "velocity_2", "velocity_3", np.arange(5260000,10000000,500000),
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=0.0001, label_scale=0.0001,
        ifloglength=True,
        ifplotseparateupdown=True,
        quiver_scale_up=0.0001, label_scale_up=0.0001,
        quiver_scale_down=0.0005, label_scale_down=0.0005,
        contour_norm = "log",
        contour_v_min_max = "constant",
        vmin=10**-5, vmax=10**-1,
    )
    
    plot_quiver_from_chunk2D_fraction(
        lmp_path,
        51, "fraction", "fraction", np.arange(5260000,10000000,500000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0.4, vmax=0.7,
        ave_y=3, ave_z=3,
    )

    plot_quiver_from_chunk2D_fraction(
        lmp_path,
        201, "fraction", "fraction", np.arange(10000000,75000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0.4, vmax=0.7,
        ave_y=3, ave_z=3,
    )
    
    """
    # all -->   teach c1each c2sum     legend t each

    # 1Dplot from chunk2D
    # velocity

    # wall stress
    # contact number 
    #for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
    #    (np.arange(1000000,40000000,1000000),"strain", "fraction", 1, None, 1, None, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "n_contact", 1, None, 1, None, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"I_12", "mu_12_middle", I_scale, None, mu_scale, mu_scale_str, 'log'),
    #    (np.arange(1000000,40000000,1000000),"I_tensor", "mu_tensor_12", I_tensor_scale, I_tensor_scale_str, mu_tensor_scale, mu_tensor_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "I_12", 1, None, I_scale, None, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "mu_12_middle", 1, None, mu_scale, mu_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "I_tensor", 1, None, I_scale, None, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "mu_tensor_12", 1, None, mu_scale, mu_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "velocity_1", 1, None, velocity_scale, velocity_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "velocity_2", 1, None, velocity_scale, velocity_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "velocity_3", 1, None, velocity_scale, velocity_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "strain_rate_21_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "strain_rate_31_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "strain_rate_22_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "strain_rate_32_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "strain_rate_23_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "strain_rate_33_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "stress_11", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "stress_22", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "stress_33", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "stress_12", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "stress_13", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"strain", "stress_23", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"fraction", "mu_12_middle", 1, None, mu_scale, mu_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"fraction", "I_12", 1, None, I_scale, None, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"fraction", "mu_tensor_12", 1, None, mu_tensor_scale, mu_tensor_scale_str, 'linear'),
    #    (np.arange(1000000,40000000,1000000),"fraction", "I_tensor", 1, None, I_tensor_scale, I_tensor_scale_str, 'linear'),
    #]:
    
    """
    # contour for pressure and stress12 mu_12
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "strain_rate_21_middle", "strain_rate_21_middle", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*strain_rate_scale, x_scale_str=None, y_scale_factor=2**0.5*strain_rate_scale, y_scale_str=None,
        quiver_scale=strain_rate_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=5,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "mu_12", "mu_12", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=0.8,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "pressure", "pressure", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale_width, x_scale_str=None, y_scale_factor=2**0.5*stress_scale_width, y_scale_str=None,
        quiver_scale=stress_scale_width, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=10,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "stress_12", "stress_12", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale_width, x_scale_str=None, y_scale_factor=2**0.5*stress_scale_width, y_scale_str=None,
        quiver_scale=stress_scale_width, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=3,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        11, "fraction", "fraction", np.arange(5100000,6500000,100000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0.3, vmax=0.7,
    )
    

    # contour for pressure and stress12 mu_12
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "strain_rate_21_middle", "strain_rate_21_middle", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*strain_rate_scale, x_scale_str=None, y_scale_factor=2**0.5*strain_rate_scale, y_scale_str=None,
        quiver_scale=strain_rate_scale, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=5,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "mu_12", "mu_12", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*1, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=0.8,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "pressure", "pressure", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale_width, x_scale_str=None, y_scale_factor=2**0.5*stress_scale_width, y_scale_str=None,
        quiver_scale=stress_scale_width, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=10,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "stress_12", "stress_12", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5*stress_scale_width, x_scale_str=None, y_scale_factor=2**0.5*stress_scale_width, y_scale_str=None,
        quiver_scale=stress_scale_width, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0, vmax=3,
    )
    plot_quiver_from_chunk2D(
        lmp_path,
        301, "fraction", "fraction", np.arange(7000000,40000000,5000000),
        spaceave=None,
        x_scale_factor=2**0.5, x_scale_str=None, y_scale_factor=2**0.5*1, y_scale_str=None,
        quiver_scale=1, label_scale=1,
        ifloglength=False,
        ifstreamplot=False,
        ifcontour=True,
        contour_norm="linear",
        contour_v_min_max="constant", # or "min_to_max",
        vmin=0.3, vmax=0.7,
    )
    
    # 1Dplot from chunk2D
    # velocity
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.append(np.arange(4500000, 7000000, 500000), np.array([40000000])), "Coord1", "mv_1", coord_scale, coord_scale_str, velocity_scale, velocity_scale_str, 'log'),
        (np.append(np.arange(4500000, 7000000, 500000), np.array([40000000])), "Coord1", "mv_1", coord_scale, coord_scale_str, velocity_scale, velocity_scale_str, 'linear'),
    ]:
        coord1_index_array = np.arange(16)
        plot_1D_from_chunk2D(
            lmp_path,
            51, x_name, y_name, input_stepsarray, coord1_index_array, 'sumdividemass', legend_class = 't',
            x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
            x_scale=x_scale,
            useerrorbar=False,
        )
    
    

    # wall stress
    
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.arange(5000000,40000000,1000000),"strain", "inwall_stress_3", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "outwall_stress_1", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "outwall_stress_2", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "outwall_stress_3", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "zbottom_stress_1", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "zbottom_stress_2", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "zbottom_stress_3", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    ]:
        plot_1D_for_chunk1D_near_wall(
            lmp_path,
            51, x_name, y_name, input_stepsarray, coord_index_array=None, legend_class = 'c',
            x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
            x_scale=x_scale,
            useerrorbar=False,
            ifsum=True,
            showlegend=False,
        )
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.arange(5000000,40000000,1000000),"strain", "inwall_stress_1", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(5000000,40000000,1000000),"strain", "inwall_stress_2", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
    ]:
        plot_1D_for_chunk1D_near_wall(
            lmp_path,
            51, x_name, y_name, input_stepsarray, coord_index_array=None, legend_class = 'c',
            x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
            x_scale=x_scale,
            useerrorbar=False,
            ifsum=True,
            setymin = True,
            ylim_bottom = 0,
            showlegend=False,
        )
    # contact number
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale, ifdivideNcount) in [
        (np.arange(1000000,40000000,1000000),"strain", "n_contact", 1, None, 1, None, 'linear', True),
        (np.arange(1000000,40000000,1000000),"strain", "fraction", 1, None, 1, None, 'linear', False),
        ]:
        for (coord1_index_array, coord2_index_array, legend_class, ave_over_coord1, ave_over_coord2) in [
            ([0, 1, 2], 'ave', 'c1', True, False),
            ([-1, -2, -3], 'ave', 'c1', True, False),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                )
                plot_1D_from_chunk2D_mask(
                    lmp_path,
                    51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                    ifmaskstatic=False,
                    ifmasknonstatic=True,
                )
    # mu I
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, y_scale) in [
        (np.arange(6000000,40000000,1000000), "mu_12_middle", "I_12", mu_scale, mu_scale_str, I_scale, None, 'log'),
        (np.array([30000000]), "mu_12_middle", "I_12", mu_scale, mu_scale_str, I_scale, None, 'log'),
        # strain
        (np.array([30000000]), "strain", "pressure", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.array([30000000]), "strain", "stress_12", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.array([30000000]), "strain", "strain_rate_21_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
    ]:
        for (coord1_index_array, coord2_index_array) in [
                ([0], np.arange(0,16,2)),
                ([-1], np.arange(0,16,2)),
                (np.arange(0,15,2), [0]),
                ]:
            plot_1D_from_chunk2D_mask(
                lmp_path,
                51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                x_scale=x_scale,
                y_scale=y_scale,
                useerrorbar=useerrorbar,
                ifmaskstatic=False,
                ifmasknonstatic=True,
            )
            plot_1D_from_chunk2D(
                lmp_path,
                51, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                x_scale=x_scale,
                y_scale=y_scale,
                useerrorbar=useerrorbar,
            )
    
    
    # 1Dplot from chunk2D

    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.arange(1000000,40000000,1000000),"strain", "fraction", 1, None, 1, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "n_contact", 1, None, 1, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"I_12", "mu_12_middle", I_scale, None, mu_scale, mu_scale_str, 'log'),
        (np.arange(1000000,40000000,1000000),"I_tensor", "mu_tensor_12", I_tensor_scale, I_tensor_scale_str, mu_tensor_scale, mu_tensor_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "I_12", 1, None, I_scale, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "mu_12_middle", 1, None, mu_scale, mu_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "I_tensor", 1, None, I_scale, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "mu_tensor_12", 1, None, mu_scale, mu_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "velocity_1", 1, None, velocity_scale, velocity_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "velocity_2", 1, None, velocity_scale, velocity_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "velocity_3", 1, None, velocity_scale, velocity_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_21_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_31_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_22_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_32_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_23_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "strain_rate_33_middle", 1, None, strain_rate_scale, strain_rate_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_11", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_22", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_33", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_12", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_13", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "stress_23", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "mu_12_middle", 1, None, mu_scale, mu_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "I_12", 1, None, I_scale, None, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "mu_tensor_12", 1, None, mu_tensor_scale, mu_tensor_scale_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"fraction", "I_tensor", 1, None, I_tensor_scale, I_tensor_scale_str, 'linear'),
    ]:
        
        for (coord1_index_array, coord2_index_array) in [
            ([0], [0,1,5,7,8]),
            ([-1], [0,1,5,7,8]),
            ([0,1,2,3,4,5,6,-2,-1], [0]),
            ([0,1,2,3,4,5,6,-2,-1], [1]),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                )
                plot_1D_from_chunk2D_mask(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = 'c1c2',
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ifmaskstatic=False,
                    ifmasknonstatic=True,
                )
        for (coord1_index_array, coord2_index_array, legend_class) in [
            ([0,-1], 'ave', 'c1'),
            ('ave', [0], 'c2'),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                )
                #plot_1D_from_chunk2D_mask(
                #    lmp_path,
                #    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                #    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                #    x_scale=x_scale,
                #    useerrorbar=useerrorbar,
                #    ifmaskstatic=False,
                #    ifmasknonstatic=True,
                #)
    
    
    
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale, ifdivideNcount) in [
        (np.arange(1000000,40000000,1000000),"strain", "n_contact", 1, None, 1, None, 'linear', True),
        (np.arange(1000000,40000000,1000000),"strain", "fraction", 1, None, 1, None, 'linear', False),
        ]:
        for (coord1_index_array, coord2_index_array, legend_class, ave_over_coord1, ave_over_coord2) in [
            ([0, 1], 'ave', 'c1', True, False),
            ([-1, -2], 'ave', 'c1', True, False),
            ([0, 1, 2], 'ave', 'c1', True, False),
            ([1, 2, 3], 'ave', 'c1', True, False),
            ([2, 3, 4], 'ave', 'c1', True, False),
            ([3, 4, 5], 'ave', 'c1', True, False),
            ([4, 5, 6], 'ave', 'c1', True, False),
            ([5, 6, 7], 'ave', 'c1', True, False),
            ([-1, -2, -3], 'ave', 'c1', True, False),
            ('ave', [0], 'c2', False, False),
            ([0, 1], [0,1,2,3,4,5,6,7], None, True, True),
            ([-1, -2], [0,1,2,3,4,5,6,7], None, True, True),
            ([0, 1, 2], [0,1,2,3,4,5,6,7], None, True, True),
            ([1, 2, 3], [0,1,2,3,4,5,6,7], None, True, True),
            ([2, 3, 4], [0,1,2,3,4,5,6,7], None, True, True),
            ([3, 4, 5], [0,1,2,3,4,5,6,7], None, True, True),
            ([4, 5, 6], [0,1,2,3,4,5,6,7], None, True, True),
            ([5, 6, 7], [0,1,2,3,4,5,6,7], None, True, True),
            ([6, 7, 8], [0,1,2,3,4,5,6,7], None, True, True),
            ([7, 8, 9], [0,1,2,3,4,5,6,7], None, True, True),
            ([8, 9, 10], [0,1,2,3,4,5,6,7], None, True, True),
            ([9, 10, 11], [0,1,2,3,4,5,6,7], None, True, True),
            ([10, 11, 12], [0,1,2,3,4,5,6,7], None, True, True),
            ([11, 12, 13], [0,1,2,3,4,5,6,7], None, True, True),
            ([12, 13, 14], [0,1,2,3,4,5,6,7], None, True, True),
            ([-1, -2, -3], [0,1,2,3,4,5,6,7], None, True, True),
            ]:
            for useerrorbar in [True, False]:
                plot_1D_from_chunk2D(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                )
                plot_1D_from_chunk2D_mask(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord1_index_array, coord2_index_array, legend_class = legend_class,
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                    ave_over_coord1=ave_over_coord1, ave_over_coord2=ave_over_coord2,
                    ifdivideNcount=ifdivideNcount,
                    ifmaskstatic=False,
                    ifmasknonstatic=True,
                )
    
    
    # 1Dplot from chunk1D
    for (input_stepsarray, x_name, y_name, x_scale_factor, x_scale_str, y_scale_factor, y_scale_str, x_scale) in [
        (np.arange(1000000,40000000,1000000),"strain", "inwall_stress_1", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "inwall_stress_2", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "inwall_stress_3", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "outwall_stress_1", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "outwall_stress_2", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "outwall_stress_3", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "zbottom_stress_1", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "zbottom_stress_2", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        (np.arange(1000000,40000000,1000000),"strain", "zbottom_stress_3", 1, None, stress_scale_width, stress_scale_width_str, 'linear'),
        ]:
        for (coord_index_array) in [
            [0,1,2,3,4,5,-2,-1],
            ]:
            for useerrorbar in [True, False]:
                plot_1D_for_chunk1D_near_wall(
                    lmp_path,
                    n_ave, x_name, y_name, input_stepsarray, coord_index_array, legend_class = 'c',
                    x_scale_factor=x_scale_factor, x_scale_str=x_scale_str, y_scale_factor=y_scale_factor, y_scale_str=None,
                    x_scale=x_scale,
                    useerrorbar=useerrorbar,
                )
    
    

    # 2D plot (vector field, streamlines, etc)
    #plot_quiver
    for ifloglength in [True, False]:
        # mv
        
        plot_quiver_from_chunk2D(
            lmp_path,
            n_ave, "mv_2", "mv_3", np.arange(1000000,40000000,1000000),
            spaceave=None,
            x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
            quiver_scale=5*10**-9, label_scale=5*10**-9,
            ifloglength=ifloglength,
            ifplotseparateupdown=True,
            quiver_scale_up=5*10**-9, label_scale_up=5*10**-9,
            quiver_scale_down=5*10**-9, label_scale_down=5*10**-9,
        )
        plot_quiver_from_chunk2D(
            lmp_path,
            n_ave, "velocity_2", "velocity_3", np.arange(1000000,40000000,1000000),
            spaceave=None,
            x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
            quiver_scale=0.0001, label_scale=0.0001,
            ifloglength=ifloglength,
            ifplotseparateupdown=True,
            quiver_scale_up=0.0001, label_scale_up=0.0001,
            quiver_scale_down=0.0005, label_scale_down=0.0005,
        )

        plot_quiver_from_chunk2D(
            lmp_path,
            n_ave, "velocity_1", "velocity_1", np.arange(1000000,40000000,1000000),
            spaceave=None,
            x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
            quiver_scale=0.002, label_scale=0.002,
            ifloglength=ifloglength,
        )
    # fraction
    plot_quiver_from_chunk2D(
        lmp_path,
        n_ave, "fraction", "fraction", np.arange(1000000,40000000,1000000),
        spaceave=None,
        x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
        quiver_scale=0.2, label_scale=1,
        ifloglength=False,
        valueminus=0.5,
    )
    

    # plot shear rate for compare with quiver
    """
    
    
    print('\a')

# main exclusive
if __name__ == "__main__":
    
    main()