




# pics info
    # from simulation, geometry, size,
    # variable_component

# single pics, single plot,  1D-2D,   stress i j at x23 single timestep
    # n_ave, time
# multi pics, single plot,   1D-2D, stress i j at x23 multi timestep
    # n_ave, time
# single pics, single plot,  1D-1D,    stress i j at xk=k_index
    # pics info add: xk, k_index, time
# single pics, multi plot,      1D-1D,  stress i j at xk=k_index, time=timearray
    # pics info add: xk, k_index, time



# save data

# read data

# data for plot output in the same path as figure

# plot figure



# import
import os
from io import StringIO
import re
import time
from itertools import chain
from itertools import repeat
from itertools import islice
import pandas as pd
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import pickle
import functools
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D
# import module
import datapath as dp
import osmanage as om
import read_setting.read_setting as rr
# import calculate setting
import read_setting.calculate_setting as rc
# plot style

plt.style.use('classic')
plt.rcParams.update({'font.size': 16})

if "if_inwall_wall_gran" in rr.logfile.keys():
    if rr.logfile["if_inwall_wall_gran"] == "yes":
        if "wall_gran_type" in rr.logfile.keys():
            if rr.logfile["wall_gran_type"] == "1":
                ybottomwalltype = "rough (d=0.9)"
            elif rr.logfile["wall_gran_type"] == "2":
                ybottomwalltype = "rough (d=1)"
            elif rr.logfile["wall_gran_type"] == "3":
                ybottomwalltype = "rough (d=1.1)"
            else:
                sys.exit("can not get wall gran type")
        else:
            ybottomwalltype = "rough (d=1)"
    else:
        ybottomwalltype = "smooth"
else:
    ybottomwalltype = "smooth"
height = rr.logfile["z_length_create_dp_unit"]
width = rr.logfile["width_wall_dp_unit"]
periodlength = rr.logfile["x_period_dp_unit"]
labelstring_size_walltype = "L= " + periodlength + ", W= " + width + ", H= " + height + ", " + ybottomwalltype

# define function for extract data from fix txt to dataframe
if rr.logfile["shearwall"] == "zcylinder":
    chunk_method = 'rz'
if rr.logfile["shearwall"] == "yplane":
    chunk_method = 'yz'
# map dim index to coordinate
if chunk_method == "rz":
    map_dim_index_to_coordinate = ["t", "r", "z"]
elif chunk_method == "yz":
    map_dim_index_to_coordinate = ["x", "y", "z"]

map_dim_index_to_coordinate_reverse = {}

for i, co in enumerate(map_dim_index_to_coordinate):
    map_dim_index_to_coordinate_reverse[co]=i

diameter = float(rr.logfile['dp'])
width_wall_dp_unit = int(rr.logfile['width_wall_dp_unit'])

if chunk_method == "rz":
    ri = diameter*int(rr.logfile['ri_wall_dp_unit']) 
elif chunk_method == "yz":
    x_period = diameter*int(rr.logfile['x_period_dp_unit'])
else:
    sys.exit("chunk_method wrong")
g = float(rr.logfile['g'])
d_step = int(rr.logfile['freq_ave_chunk_momentum_mass_field'])
den = float(rr.logfile['den'])
height_dpunit = float(rr.logfile['zhi_chunk_dp_unit'])
stress_scale = den*g*height_dpunit*diameter

shear_rate_scale = float(rr.logfile['in_velocity'])/(float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp']))

if chunk_method == "rz":
    position_index_to_array_dim_index = {
                                    1: 1,
                                    2: 0,
                                    }
    n_r = int(rr.logfile['N_bin_r'])
    n_z = int(rr.logfile['N_bin_z'])
    n_1 = n_z
    n_2 = n_r
    n_12 = n_1*n_2
    
    dx = 1/n_r*width_wall_dp_unit
    dy = 1/n_z*height_dpunit
    x_array, y_array = np.meshgrid(
                                int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                )
    x_array = x_array.reshape((-1))
    y_array = y_array.reshape((-1))
    vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*diameter**3
elif chunk_method == "yz":
    n_y = int(rr.logfile['N_bin_y'])
    n_z = int(rr.logfile['N_bin_z'])
    if rr.logfile["chunk/atom"][1] == "y":
        position_index_to_array_dim_index = {
                                        1: 0,
                                        2: 1,
                                        }
        n_1 = n_y
        n_2 = n_z
    elif rr.logfile["chunk/atom"][1] == "z":
        position_index_to_array_dim_index = {
                                        2: 0,
                                        1: 1,
                                        }
        n_1 = n_z
        n_2 = n_y
    else:
        sys.exit("chunk_method wrong")
    n_12 = n_1*n_2
    dx = 1/n_y*width_wall_dp_unit
    dy = 1/n_z*height_dpunit
    vol_in_chunks = x_period*dx*dy*diameter**2
else:
    sys.exit("chunk_method wrong")


if chunk_method == "rz":
    # chink first dim unchange in the begining
    chunk_first_dim_coord = "z"
    chunk_second_dim_coord = "r"
elif chunk_method == "yz":
    if rr.logfile["chunk/atom"][1] == "y":
        chunk_first_dim_coord = "y"
        chunk_second_dim_coord = "z"
        xyztoCoor = {}
        xyztoCoor["y"] = "Coord1"
        xyztoCoor["z"] = "Coord2"
    elif rr.logfile["chunk/atom"][1] == "z":
        chunk_first_dim_coord = "z"
        chunk_second_dim_coord = "y"
        xyztoCoor["z"] = "Coord1"
        xyztoCoor["y"] = "Coord2"
    else:
        sys.exit("chunk_method wrong")
else:
    sys.exit("chunk_method wrong")

chunk_first_dim_number = [i for i, coord in enumerate(map_dim_index_to_coordinate) if coord == chunk_first_dim_coord][0]
chunk_second_dim_number = [i for i, coord in enumerate(map_dim_index_to_coordinate) if coord == chunk_second_dim_coord][0]

def lines_from_one_simu(lammps_directory):
    with open(lammps_directory + "output/stress/fix.stress.all") as f:
        lines = f.read().strip().split('\n')
    return lines

lines = lines_from_one_simu(dp.lammps_directory)
header = lines[2].split()[1:]
n_line_in_a_step = int(lines[3].split()[1])


def header_from_lines(lines):
    return lines[2].split()[1:]


def dic_index_of_variable_in_header(lines):
    header = header_from_lines(lines)
    dic = {}
    for i, name in enumerate(header):
        dic[name] = i
    return dic


def step_first_in_file(lines):
    return int(lines[3].split()[0])


def step_last_in_file(lines):
    return int(lines[-1 - n_line_in_a_step].split()[0])


def step_first_in_file_change_by_n_ave(n_ave,lines):
    return step_first_in_file(lines)+int((n_ave-1)/2*d_step)


def step_last_in_file_change_by_n_ave(n_ave,lines):
    return step_last_in_file(lines)-int((n_ave-1)/2*d_step)


def step2_fix_initial_to_last_func(index):
    if index == rr.n_loglist-1:
        lines = lines_from_one_simu(dp.lammps_directory)
        step2_fix = step_last_in_file(lines)
    else:
        lmp_dir_next = rc.folder_path_list_initial_to_last[index+1]
        lines_next = lines_from_one_simu(lmp_dir_next)
        step1_next = step_first_in_file(lines_next)

        lmp_dir = rc.folder_path_list_initial_to_last[index]
        lines = lines_from_one_simu(lmp_dir)
        logfile = rr.logfilelist_from_initial_to_lastest[index]
        freq = int(logfile["freq_ave_chunk_momentum_mass_field"])
        step1 = step_first_in_file(lines)
        step2 = step_last_in_file(lines)
        
        if step1_next <= step2:
            step2_fix = int(freq*int((step1_next - step1)/freq) + step1 - freq)
        else:
            step2_fix = int(step2)

    return step2_fix


def step_last_fix_change_by_n_ave(n_ave, index):
    return step2_fix_initial_to_last_func(index)-int((n_ave-1)/2*int(rr.logfilelist_from_initial_to_lastest[index]["freq_ave_chunk_momentum_mass_field"]))


def data_in_one_step(step, lines):

    n_line_0 = int(int(step - step_first_in_file(lines))/d_step)*(n_line_in_a_step+1) + 4
    n_line_1 = int(n_line_0 + n_line_in_a_step)
    ## select data
    data = [lines[t].split() for t in range(n_line_0, n_line_1)]
    data = np.array(data, dtype=np.float64)
    return data


def divide_zero(a,b):
    c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    return c


def value_in_a_step(step, variable_name,lines):
    return data_in_one_step(step,lines)[:,dic_index_of_variable_in_header(lines)[variable_name]]


def value_in_a_step_ave(step, variable_name, n_ave, lines):

    value = 0
    for i in range(n_ave):
        value += value_in_a_step(int(step - (n_ave-1)*d_step/2 + i*d_step), variable_name, lines)
    value /= n_ave
    return value


def time_in_a_step(step):
    return step*float(rr.logfile["ts"])


def time_in_a_step_from_start_rotate(step):
    return step*float(rr.logfile["ts"])-rc.rotate_start_time


def path_nve_subfolder_in_folder(n_ave, folder):
    subfolder = folder + "nve_" + str(n_ave) + "/"
    return subfolder


def add_nve_subfolder_in_folder(n_ave, folder):
    if not os.path.isdir(folder): 
        os.mkdir(folder)
    if not os.path.isdir(path_nve_subfolder_in_folder(n_ave, folder)): 
        os.mkdir(path_nve_subfolder_in_folder(n_ave, folder))


## empty grid
def if_grid_surround_empty(mass):
    if len(mass.shape) != 2:
        sys.exit("mass shape len is not 2")
    mass_n_1 = mass.shape[0]
    mass_n_2 = mass.shape[1]
    ifempty = (mass==0)
    expand_ifempty = np.empty([mass_n_1+2, mass_n_2+2], dtype=bool)
    expand_ifempty[:,0]  = False
    expand_ifempty[:,-1] = False
    expand_ifempty[0,:]  = False
    expand_ifempty[-1,:] = False
    expand_ifempty[1:-1, 1:-1] = ifempty
    
    center  = expand_ifempty[1:-1, 1:-1]
    left    = expand_ifempty[1:-1,  :-2]
    right   = expand_ifempty[1:-1, 2:  ]
    down    = expand_ifempty[:-2 , 1:-1]
    up      = expand_ifempty[2:  , 1:-1]

    return functools.reduce(np.logical_or, (center, left, right, down, up))


def if_grid_surround_not_empty(mass):
    return np.logical_not(if_grid_surround_empty(mass))


def initial_plot():
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    return (fig, ax)


def plot_quiver_position_label(fig, ax):
    
    if chunk_method == "rz":
        plt.xlabel('r')
        plt.ylabel('z')
    elif chunk_method == "yz":
        plt.xlabel('y')
        plt.ylabel('z')
    
    return (fig, ax)


def save_one_plot(fig, ax, foldersave, f_name, figformat="png", ifpickle=False):
    
    fig.savefig(foldersave + f_name + "." + figformat, format=figformat)
    if ifpickle:
        # Save figure handle to disk
        with open(foldersave + f_name + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
            pickle.dump(fig, f)
    plt.close('all')


class chunk(object):

    def __init__(self, n_ave, lmp_path):

        self.n_ave = n_ave
        self.lmp_path = lmp_path
        self.logfile = rr.read_log(self.lmp_path)
        self.lines = lines_from_one_simu(self.lmp_path)
        self.header = self.lines[2].split()[1:]
        self.n_line_in_a_step = int(self.lines[3].split()[1])
        self.step_first_in_file = int(self.lines[3].split()[0])
        self.step_second_in_file = int(self.lines[3 + self.n_line_in_a_step + 1].split()[0])
        self.step_last_in_file = int(self.lines[-1 - self.n_line_in_a_step].split()[0])
        self.d_step = self.step_second_in_file - self.step_first_in_file
        self.step_first_in_file_change_by_n_ave = self.step_first_in_file + int((self.n_ave-1)/2*self.d_step)
        self.step_last_in_file_change_by_n_ave = self.step_last_in_file - int((self.n_ave-1)/2*self.d_step)
        self.middle_step = int(
            int((self.step_last_in_file_change_by_n_ave - self.step_first_in_file_change_by_n_ave)/2/self.d_step)*self.d_step + self.step_first_in_file_change_by_n_ave
            )
        self.allsteps = np.arange(self.step_first_in_file_change_by_n_ave, self.step_last_in_file_change_by_n_ave, self.d_step)
        self.first_middle_last_steps = np.array([self.step_first_in_file_change_by_n_ave, self.middle_step, self.step_last_in_file_change_by_n_ave])
        self.extrasteps = np.array(
            [
                self.step_first_in_file_change_by_n_ave + 5*self.d_step*self.n_ave,
                ]
            )
        maskextra = np.logical_and(self.extrasteps > self.step_first_in_file_change_by_n_ave, self.extrasteps < self.step_last_in_file_change_by_n_ave)
        self.extrasteps = self.extrasteps[maskextra]
        self.first_extra_middle_last_steps = np.append(self.first_middle_last_steps, self.extrasteps)
        self.first_extra_middle_last_steps.sort()
        if "if_inwall_wall_gran" in rr.logfile.keys():
            if rr.logfile["if_inwall_wall_gran"] == "yes":
                if "wall_gran_type" in rr.logfile.keys():
                    if rr.logfile["wall_gran_type"] == "1":
                        self.ybottomwalltype = "rough (d=0.9)"
                    elif rr.logfile["wall_gran_type"] == "2":
                        self.ybottomwalltype = "rough (d=1)"
                    elif rr.logfile["wall_gran_type"] == "3":
                        self.ybottomwalltype = "rough (d=1.1)"
                    else:
                        sys.exit("can not get wall gran type")
                else:
                    self.ybottomwalltype = "rough (d=1)"
            else:
                self.ybottomwalltype = "smooth"
        else:
            self.ybottomwalltype = "smooth"

        self.height = rr.logfile["z_length_create_dp_unit"]
        self.width = rr.logfile["width_wall_dp_unit"]
        self.periodlength = rr.logfile["x_period_dp_unit"]
        self.labelstring_size_walltype = self.ybottomwalltype + "\n" + "L " + self.periodlength + "\n" + "W " + self.width + "\n" + "H " + self.height
        self.labelstring_size_walltype_one_line = self.ybottomwalltype + ", " + "L " + self.periodlength + ", " + "W " + self.width + ", " + "H " + self.height
        

    def data_in_one_step(self, step):

        n_line_0 = int(int(step - self.step_first_in_file)/self.d_step*(self.n_line_in_a_step+1) + 4)
        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
        ## select data
        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
        data = np.array(data, dtype=np.float64)
        
        return data

    def value_in_a_step(self, step, variable_name):
        
        return self.data_in_one_step(step)[:, self.header.index(variable_name)]

    def value_in_a_step_ave(self, step, variable_name):

        value = 0
        for i in range(self.n_ave):
            value += self.value_in_a_step(int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step), variable_name)
        value /= self.n_ave
        
        return value

    def stress_variable_name(self, i, j):

        if i == 1 and j == 1:
            variable_name = "c_stress[1]"
        elif i == 2 and j == 2:
            variable_name = "c_stress[2]"
        elif i == 3 and j == 3:
            variable_name = "c_stress[3]"
        elif i == 1 and j == 2:
            variable_name = "c_stress[4]"
        elif i == 1 and j == 3:
            variable_name = "c_stress[5]"
        elif i == 2 and j == 3:
            variable_name = "c_stress[6]"
        else:
            sys.exit("i j index not correct, only ij 11 22 33 12 13 23")

        return variable_name


    ########## plot 1D-2D stress-position ##########
    def datachunk_ave_one_step_stressij_x23(self, step, i, j):

        stress_array = value_in_a_step_ave(step, self.stress_variable_name(i,j), self.n_ave, self.lines)

        if chunk_method == "rz":
            r = value_in_a_step_ave(step, "v_r", self.n_ave, self.lines)
            z = value_in_a_step_ave(step, "v_z", self.n_ave, self.lines)
            x_array = r/diameter
            y_array = z/diameter
        elif chunk_method == "yz":
            y = value_in_a_step_ave(step, "Coord1", self.n_ave, self.lines)
            z = value_in_a_step_ave(step, "Coord2", self.n_ave, self.lines)
            x_array = y/diameter
            y_array = z/diameter
        else:
            sys.exit("chunk_method wrong")
        
        stress_array = stress_array/stress_scale/vol_in_chunks
        time = time_in_a_step_from_start_rotate(step)

        return [time, x_array, y_array, stress_array]



    ########## plot 1D-2D stress-position ##########
    def plotchunk_ave_one_step_stressij_x23(self, step, i, j, quiver_scale=0.1, label_scale=0.2):

        [time, x_array, y_array, stress_array] = self.datachunk_ave_one_step_stressij_x23(step, i, j)

        fig1, ax1 = plt.subplots()
        #fig1.figsize = [12.8, 9.6]
        plot_quiver_position_label(fig1, ax1)
        #ax1.set_title('velocity field r-z direction (average over theta)')
        Q = ax1.quiver(x_array, y_array, np.zeros_like(stress_array), stress_array,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                    )

        ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                    label = "equal stress = {:.2e}".format(label_scale) + ". At {:.2e} s".format(time),
                    labelpos='E',
                    coordinates='figure', angle=90)
        return (fig1, ax1)

    ########## save plot  ##########
    def save_stress_x23(self, stepsarray, i, j, figformat="png", ifpickle=False):
        savepath = dp.stress_folder_path(i,j)
        add_nve_subfolder_in_folder(self.n_ave, savepath)
        foldersave = path_nve_subfolder_in_folder(self.n_ave, savepath)
        for step in stepsarray:
            fig, ax = self.plotchunk_ave_one_step_stressij_x23(step, i, j)
            save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)

    def datachunk_stress_ij_fix_k_ave(self, step, i, j, k, k_index):

        fix_along_array_dim = position_index_to_array_dim_index[k]

        stress_array = value_in_a_step_ave(step, self.stress_variable_name(i,j), self.n_ave, self.lines)
        
        stress_array = stress_array/stress_scale/vol_in_chunks
        stress = np.take(stress_array, k_index, axis=fix_along_array_dim)
        
        if k == 2:
            anotherdim = 3
        elif k == 3:
            anotherdim = 2
        else:
            sys.exit("not 2 not 3")

        vector_origin = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[anotherdim-1]], self.n_ave, self.lines)
        vector = np.take(vector_origin, k_index, axis=fix_along_array_dim)

        time = time_in_a_step_from_start_rotate(step)
        
        return [time, vector, stress]

    ########## plot 1D-1D stress-position index change##########
    def plotchunk_stress_ij_fix_k_ave(self, step, i, j, k, k_index):

        [time, vector, stress] = self.datachunk_stress_ij_fix_k_ave(step, i, j, k, k_index)
        fig, ax = plt.subplots()
        ax.plot(vector, stress, label="{:.2e}".format(time),
                marker = ".",
                linestyle = 'None',
                markersize=16,
                )
        
        return (fig, ax)
    
    def save_stress_fix_k(self, stepsarray, i, j, k, k_index, figformat="png", ifpickle=False):
        savepath = dp.stress_fix_k_folder_path(i,j,k)
        om.create_directory(savepath)
        add_nve_subfolder_in_folder(self.n_ave, savepath)
        foldersave = path_nve_subfolder_in_folder(self.n_ave, savepath)
        for step in stepsarray:
            fig, ax = self.plotchunk_stress_ij_fix_k_ave(step, i, j, k, k_index)
            save_one_plot(fig, ax, foldersave, "step_" + str(int(step) + "k" + str(k_index)), figformat="png", ifpickle=False)


    def plotchunk_1(self, stepsarray, figformat="png", ifpickle=False):
        self.save_stress_x23(stepsarray, 1, 1, figformat="png", ifpickle=False)
        self.save_stress_x23(stepsarray, 1, 2, figformat="png", ifpickle=False)
        self.save_stress_x23(stepsarray, 1, 3, figformat="png", ifpickle=False)
        self.save_stress_x23(stepsarray, 2, 2, figformat="png", ifpickle=False)
        self.save_stress_x23(stepsarray, 2, 3, figformat="png", ifpickle=False)
        self.save_stress_x23(stepsarray, 3, 3, figformat="png", ifpickle=False)

    def plotchunk_1_step1to2(self, func, step1, step2, figformat="png", ifpickle=False):
        #if smallstep largestep
        # need combine
        stepsarray = np.arange(step1, step2, d_step)
        return func(stepsarray, figformat="png", ifpickle=False)

    def plotchunk_1_every(self, func, figformat="png", ifpickle=False):
        # need combine
        step1_default = step_first_in_file_change_by_n_ave(self.n_ave, self.lines)
        step2_default = step_last_in_file_change_by_n_ave(self.n_ave, self.lines)
        return self.plotchunk_1_step1to2(func, step1_default, step2_default, figformat="png", ifpickle=False)

    def plotchunk_1_list(self, func, stepslist, figformat="png", ifpickle=False):
        #if smallstep largestep
        # need combine
        stepsarray = np.array(stepslist)
        return func(stepsarray, figformat="png", ifpickle=False)

