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

if "if_ybottom_wall_gran" in rr.logfile.keys():
    if rr.logfile["if_ybottom_wall_gran"] == "yes":
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
velocity_scale = float(rr.logfile['in_velocity'])
if velocity_scale < 0:
    velocity_scale = -velocity_scale

if velocity_scale == 0:
    Sa_fake = 0.000002
    if chunk_method == "rz":
        velocity_scale = (Sa_fake*g*width_wall_dp_unit**3*diameter)**0.5
    elif chunk_method == "yz":
        velocity_scale = (Sa_fake*g*width_wall_dp_unit**3*diameter)**0.5

height_dpunit = float(rr.logfile['zhi_chunk_dp_unit'])
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
    with open(lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
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

    n_line_0 = int(int(step - step_first_in_file(lines))/d_step*(n_line_in_a_step+1) + 4)
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
        if "if_ybottom_wall_gran" in rr.logfile.keys():
            if rr.logfile["if_ybottom_wall_gran"] == "yes":
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


    ########## plot 2D-2D velocity-position ##########
    def plotchunk_ave_one_step_v23x23(self, step, figformat="png", ifpickle=False):
        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
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
        if chunk_method == "rz":
            mvr = value_in_a_step_ave(step, "v_mvr", self.n_ave, self.lines)
            mvz = value_in_a_step_ave(step, "v_mvz", self.n_ave, self.lines)
            vx_array = divide_zero(mvr,mass)/velocity_scale
            vy_array = divide_zero(mvz,mass)/velocity_scale
        elif chunk_method == "yz":
            mvy = value_in_a_step_ave(step, "v_mvy", self.n_ave, self.lines)
            mvz = value_in_a_step_ave(step, "v_mvz", self.n_ave, self.lines)
            vx_array = divide_zero(mvy,mass)/velocity_scale
            vy_array = divide_zero(mvz,mass)/velocity_scale
        else:
            sys.exit("chunk_method wrong")
        v_length_array = (vx_array**2+vy_array**2)**0.5
        max_v_length = np.amax(v_length_array)
        if dp.if_plot_velocity_field_scale_same == "yes":
            quiver_scale = dp.quiver_scale_velocity_xaxis_shearplanenormal_yaxis_z
        else:
            quiver_scale = max_v_length/2
        label_scale = quiver_scale

        fig1, ax1 = plt.subplots()
        plot_quiver_position_label(fig1, ax1)
        time = time_in_a_step_from_start_rotate(step)
            
        Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                    )

        ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                    " : {:.2e} wall velocity in 45 degree".format(label_scale)+ ". At {:.2e} s".format(time),
                    labelpos='E', coordinates='figure', angle=45)
        
        return (fig1, ax1)

    ########## plot 2D-2D velocity-position ##########
    def plotchunk_ave_one_step_v13x23(self, step ,figformat="png", ifpickle=False):
        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
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

        if chunk_method == "rz":
            mvt = value_in_a_step_ave(step, "v_mvt", self.n_ave, self.lines)
            mvz = value_in_a_step_ave(step, "v_mvz", self.n_ave, self.lines)
            vx_array = divide_zero(mvt,mass)/velocity_scale
            vy_array = divide_zero(mvz,mass)/velocity_scale
        elif chunk_method == "yz":
            mvx = value_in_a_step_ave(step, "v_mvx", self.n_ave, self.lines)
            mvz = value_in_a_step_ave(step, "v_mvz", self.n_ave, self.lines)
            vx_array = divide_zero(mvx,mass)/velocity_scale
            vy_array = divide_zero(mvz,mass)/velocity_scale
        else:
            sys.exit("chunk_method wrong")
        
        v_length_array = (vx_array**2+vy_array**2)**0.5
        max_v_length = np.amax(v_length_array)

        if dp.if_plot_velocity_field_scale_same == "yes":
            quiver_scale = dp.quiver_scale_velocity_xaxis_shearplanenormal_yaxis_z
        else:
            quiver_scale = max_v_length/2
        label_scale = quiver_scale

        fig1, ax1 = plt.subplots()
        plot_quiver_position_label(fig1, ax1)
        time = time_in_a_step_from_start_rotate(step)    
        Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                    )

        ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                    label = " : {:.2e} wall velocity in 45 degree".format(label_scale),
                    labelpos='E', coordinates='figure', angle=45)
        return (fig1, ax1)

    ########## plot 1D-2D volumn fraction-position ##########
    def plotchunk_ave_one_step_volumnfraction_x23(self, step, figformat="png", ifpickle=False):
        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
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
        quiver_scale = 0.2
        label_scale = 0.6
        vy_array = mass/float(rr.logfile['den'])/vol_in_chunks
        fig1, ax1 = plt.subplots()
        
        #fig1.figsize = [12.8, 9.6]
        plot_quiver_position_label(fig1, ax1)
        time = time_in_a_step_from_start_rotate(step)
        #ax1.set_title('velocity field r-z direction (average over theta)')
        Q = ax1.quiver(x_array, y_array, np.zeros_like(vy_array), vy_array,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                    )

        ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                    label = "equal solid fraction = {:.2e}".format(label_scale) + ". At {:.2e} s".format(time),
                    labelpos='E',
                    coordinates='figure', angle=90)
        return (fig1, ax1)

    ########## save plot 2D-2D velocity-position ##########
    def save_v23_x23(self, stepsarray, figformat="png", ifpickle=False):
        add_nve_subfolder_in_folder(self.n_ave, dp.f_momentum_mass_field_v23x23_path)
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_momentum_mass_field_v23x23_path)
        for step in stepsarray:
            fig, ax = self.plotchunk_ave_one_step_v23x23(step, figformat="png", ifpickle=False)
            save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)
            
    def save_v13_x23(self, stepsarray, figformat="png", ifpickle=False):
        add_nve_subfolder_in_folder(self.n_ave, dp.f_momentum_mass_field_v13x23_path)
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_momentum_mass_field_v13x23_path)
        for step in stepsarray:
            fig, ax = self.plotchunk_ave_one_step_v13x23(step, figformat="png", ifpickle=False)
            save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)
    
    ########## save plot 1D-2D volumn fraction-position ##########
    def save_fraction_x23(self, stepsarray, figformat="png", ifpickle=False):
        add_nve_subfolder_in_folder(self.n_ave, dp.f_momentum_mass_field_volumnfraction_x23_path)
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_momentum_mass_field_volumnfraction_x23_path)
        for step in stepsarray:
            fig, ax = self.plotchunk_ave_one_step_volumnfraction_x23(step, figformat="png", ifpickle=False)
            save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)

    ########## plot 1D-2D strain_rate-position ##########
    def plotchunk_strain_rate_i_j_x23(self, step, i, j, figformat="png", ifpickle=False):
        
        mass = np.resize(
            value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines),
            (n_1, n_2),
            )         
        vector_mv = np.resize(
            value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines),
            (n_1, n_2),
            ) 
        vector_x2 = np.resize(
            value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[1]], self.n_ave, self.lines),
            (n_1, n_2),
            ) 
        vector_x3 = np.resize(
            value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[2]], self.n_ave, self.lines),
            (n_1, n_2),
            ) 

        vector_vi = divide_zero(vector_mv, mass)
        
        if j == 1:
            vector_xi = vector_x2
        elif j == 2:
            vector_xi = vector_x3

        diff_along_array_dim = position_index_to_array_dim_index[j]
        if diff_along_array_dim == 0:
            middle_point_vector_x2 = (vector_x2[:-1,:] + vector_x2[1:,:])/2
            middle_point_vector_x3 = (vector_x3[:-1,:] + vector_x3[1:,:])/2
        elif diff_along_array_dim == 1:
            middle_point_vector_x2 = (vector_x2[:,:-1] + vector_x2[:,1:])/2
            middle_point_vector_x3 = (vector_x3[:,:-1] + vector_x3[:,1:])/2
        else:
            sys.exit("error")
        
        vector_vi_diff = np.diff(vector_vi,axis=diff_along_array_dim)
        vector_xi_diff = np.diff(vector_xi,axis=diff_along_array_dim)

        len_diff_axis = mass.shape[diff_along_array_dim]
        A = np.take(mass, np.arange(len_diff_axis-1), axis=diff_along_array_dim)
        B = np.take(mass, np.arange(1, len_diff_axis), axis=diff_along_array_dim)
        maskbothnotzero = np.logical_and(
            (A!=0), 
            (B!=0),
            )

        vector_vi_diff = vector_vi_diff*maskbothnotzero
    
        strain_rate = divide_zero(vector_vi_diff,vector_xi_diff)
        strain_rate /= shear_rate_scale
        
        middle_point_vector_x2 /= diameter
        middle_point_vector_x3 /= diameter

        quiver_scale = 0.5
        label_scale = quiver_scale
        fig1, ax1 = plt.subplots()
        
        plot_quiver_position_label(fig1, ax1)
        Q = ax1.quiver(middle_point_vector_x2, middle_point_vector_x3, np.zeros_like(vector_vi_diff), strain_rate,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                    )
        time = time_in_a_step_from_start_rotate(step)
        labelstring = (
            " : strain_rate_"
            + map_dim_index_to_coordinate[i]
            + map_dim_index_to_coordinate[j]
            + " = {:.2e}".format(label_scale)
            + ". At {:.2e} s".format(time)
        )

        ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                    label = labelstring,
                    labelpos='E',
                    coordinates='figure', angle=90)
        return (fig1, ax1)
        
    ########## save plot 1D-2D strain_rate-position ##########
    def save_strain_rate_ij_x23(self, stepsarray, i, j, figformat="png", ifpickle=False):
        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_strain_rate_i_j_x23(i,j))
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_strain_rate_i_j_x23(i,j))
        for step in stepsarray:
            fig, ax = self.plotchunk_strain_rate_i_j_x23(step, i, j, figformat="png", ifpickle=False)
            save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)

    ########## plot 1D-1D strain_rate-position index change##########
    def plotchunk_strain_rate_ij_ave_k_ave(self, step, fig, ax, i, j, k, figformat="png", ifpickle=False):

        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        mass = np.resize(mass, (n_1, n_2))
        
        vector_mv = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)
        vector_mv = np.resize(vector_mv, (n_1, n_2))
        vector_vi = divide_zero(vector_mv, mass)
        vector_j = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[j]], self.n_ave, self.lines)

        vector_j = np.resize(vector_j, (n_1, n_2))

        sum_along_array_dim = position_index_to_array_dim_index[k]

        diff_along_array_dim = position_index_to_array_dim_index[j]
        if diff_along_array_dim == 0:
            middle_point_vector_j = (vector_j[:-1,0] + vector_j[1:,0])/2
        elif diff_along_array_dim == 1:
            middle_point_vector_j = (vector_j[0,:-1] + vector_j[0,1:])/2
        else:
            sys.exit("error")

        len_in_diff_dim = mass.shape[diff_along_array_dim]
        diff_not_empty = np.logical_and(
            np.take(if_grid_surround_not_empty(mass), np.arange(len_in_diff_dim-1), axis=diff_along_array_dim), 
            np.take(if_grid_surround_not_empty(mass), np.arange(len_in_diff_dim-1)+1, axis=diff_along_array_dim),
            )
        count_not_empty = np.sum(diff_not_empty, axis=sum_along_array_dim)
        vector_vi_diff = divide_zero(np.sum(np.diff(vector_vi,axis=diff_along_array_dim)*diff_not_empty, axis=sum_along_array_dim), count_not_empty)
        vector_j_diff = np.sum(np.diff(vector_j,axis=diff_along_array_dim),axis=sum_along_array_dim)/vector_j.shape[sum_along_array_dim]
        strain_rate = divide_zero(vector_vi_diff,vector_j_diff)
        strain_rate /= shear_rate_scale
        middle_point_vector_j /= diameter

        time = time_in_a_step_from_start_rotate(step)
        ax.plot(middle_point_vector_j, strain_rate, label="{:.2e}".format(time),
                marker = ".",
                linestyle = 'None',
                markersize=16,
                )
        
        return (fig, ax)

    ########## plot 1D-1D strain_rate-position index change##########
    def plotchunk_strain_rate_ij_fix_k_ave(self, step, fig, ax, i, j, k, k_index, figformat="png", ifpickle=False):

        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        
        vector_mv = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)
        vector_vi = divide_zero(vector_mv, mass)
        vector_j = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[j]], self.n_ave, self.lines)

        vector_vi = np.resize(vector_vi, (n_1, n_2))
        vector_j = np.resize(vector_j, (n_1, n_2))

        fix_along_array_dim = position_index_to_array_dim_index[k]

        diff_along_array_dim = position_index_to_array_dim_index[j]
        if diff_along_array_dim == 0:
            middle_point_vector_j = (vector_j[:-1,0] + vector_j[1:,0])/2
        elif diff_along_array_dim == 1:
            middle_point_vector_j = (vector_j[0,:-1] + vector_j[0,1:])/2
        else:
            sys.exit("error")

        vector_vi_diff = np.take(np.diff(vector_vi,axis=diff_along_array_dim), k_index, axis=fix_along_array_dim)
        vector_j_diff = np.take(np.diff(vector_j,axis=diff_along_array_dim), k_index, axis=fix_along_array_dim)
        strain_rate = divide_zero(vector_vi_diff,vector_j_diff)
        strain_rate /= shear_rate_scale
        middle_point_vector_j /= diameter

        time = time_in_a_step_from_start_rotate(step)
        ax.plot(middle_point_vector_j, strain_rate, label="{:.2e}".format(time),
                marker = ".",
                linestyle = 'None',
                markersize=16,
                )
        
        return (fig, ax)
        
    ########## plot 1D-1D strain_rate-position index change##########
    def plotchunk_strain_rate_ij_ave_k_ave_manytime(self, stepsarray, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        if ifmanysimu:
            stepsarray.sort()
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(stepsarray, step2, side='right')
                stepsarray_index = stepsarray[:sort_final_index]

                for step in stepsarray_index:
                    fig, ax = self.plotchunk_strain_rate_ij_ave_k_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                stepsarray = stepsarray[sort_final_index:]
        else:
            for step in stepsarray:
                fig, ax = self.plotchunk_strain_rate_ij_ave_k_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
            
        ax.legend(
                title = "Time (s)",
                bbox_to_anchor=(1.04,1), 
                loc="upper left"
                )
        ax.set_xlabel(map_dim_index_to_coordinate[j])
        ax.set_ylabel("strain_rate" + map_dim_index_to_coordinate[i] + map_dim_index_to_coordinate[j] + "_average_" + map_dim_index_to_coordinate[k])
        plt.tight_layout()

        return (fig, ax)
    ########## save plot 1D-1D strain_rate-position index change##########
    def save_plotchunk_strain_rate_ij_ave_k_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=False):
        
        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_strain_rate_i_j_ave_k(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_strain_rate_i_j_ave_k(i,j,k))
        
        if inonefigure:
            
            fig, ax = self.plotchunk_strain_rate_ij_ave_k_ave_manytime(stepsarray, i, j, k, figformat="png", ifpickle=False)
            
            save_one_plot(
                fig, ax, foldersave,
                "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
                figformat="png", ifpickle=False,
            )

        else:
            
            for step in stepsarray:

                fig, ax = self.plotchunk_strain_rate_ij_ave_k_ave_manytime(np.array([step]), i, j, k, figformat="png", ifpickle=False)
                
                save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)

        plt.close('all')

    ########## plot 1D-1D velocity-position index change##########
    def plotchunk_velocity_i_ave_j_xk_ave(self, step, fig, ax, i, j, k, figformat="png", ifpickle=False):

        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        mvi = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)
        vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], self.n_ave, self.lines)
        mass = np.resize(mass, (n_1, n_2))
        mvi = np.resize(mvi, (n_1, n_2))
        vector = np.resize(vector, (n_1, n_2))
        
        sum_along_array_dim = position_index_to_array_dim_index[j]
        
        mvi = np.sum(mvi,axis=sum_along_array_dim)
        mass = np.sum(mass,axis=sum_along_array_dim)
        vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
        velocity = divide_zero(mvi, mass)
        velocity /= velocity_scale

        vector = vector/diameter
        time = time_in_a_step_from_start_rotate(step)
        ax.plot(vector, velocity, label="{:.2e}".format(time),
                marker = ".",
                linestyle = 'None',
                markersize=16,
                )

        return (fig, ax)
    
    ########## plot 1D-1D velocity-position index change##########
    def plotchunk_velocity_i_ave_j_xk_ave_manytime(self, stepsarray, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if ifmanysimu:
            stepsarray.sort()
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(stepsarray, step2, side='right')
                stepsarray_index = stepsarray[:sort_final_index]
                for step in stepsarray_index:
                    fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                
                stepsarray = stepsarray[sort_final_index:]
        else:
            for step in stepsarray:
                fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
        
        ax.legend(
            title = "Time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (Vwall)")
        ax.set_title(self.labelstring_size_walltype_one_line)
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D velocity-position index change##########
    def save_plotchunk_velocity_i_ave_j_xk_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_velocity_i_j_ave_k(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_velocity_i_j_ave_k(i,j,k))
        
        if inonefigure:
            
            fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave_manytime(stepsarray, i, j, k, figformat="png", ifpickle=False)

            save_one_plot(
                fig, ax, foldersave,
                "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
                figformat="png", ifpickle=False,
            )
        else:
            
            for step in stepsarray:

                fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave_manytime(np.array([step]), i, j, k, figformat="png", ifpickle=False)
                    
                save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)
    
    ########## plot 1D-1D ek-position index change##########
    def plotchunk_ekovermass_i_ave_j_xk_ave(self, step, fig, ax, i, j, k, figformat="png", ifpickle=False):

        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        eki = value_in_a_step_ave(step, "v_Ek" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)
        vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], self.n_ave, self.lines)
        mass = np.resize(mass, (n_1, n_2))
        eki = np.resize(eki, (n_1, n_2))
        vector = np.resize(vector, (n_1, n_2))
        
        sum_along_array_dim = position_index_to_array_dim_index[j]

        eki = np.sum(eki,axis=sum_along_array_dim)
        mass = np.sum(mass,axis=sum_along_array_dim)
        vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
        ekovermass = 2*divide_zero(eki, mass)
        
        ekovermass /= velocity_scale**2
        vector = vector/diameter
        time = time_in_a_step_from_start_rotate(step)
        ax.plot(vector, ekovermass, label="{:.2e}".format(time),
                marker = ".",
                linestyle = 'None',
                markersize=16,
                )

        return (fig, ax)
    
    ########## plot 1D-1D ek-position index change##########
    def plotchunk_ekovermass_i_ave_j_xk_ave_manytime(self, stepsarray, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if ifmanysimu:
            stepsarray.sort()
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(stepsarray, step2, side='right')
                stepsarray_index = stepsarray[:sort_final_index]
                for step in stepsarray_index:
                    fig, ax = self.plotchunk_ekovermass_i_ave_j_xk_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                
                stepsarray = stepsarray[sort_final_index:]
        else:
            for step in stepsarray:
                fig, ax = self.plotchunk_ekovermass_i_ave_j_xk_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
        
        ax.legend(
            title = "Time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_title(self.labelstring_size_walltype_one_line)
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("<E" + map_dim_index_to_coordinate[i] + ">" + "/" + "<m>")
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D ek-position index change##########
    def save_plotchunk_ekovermass_i_ave_j_xk_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekovermass_i_j_ave_k(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekovermass_i_j_ave_k(i,j,k))
        
        if inonefigure:
            
            fig, ax = self.plotchunk_ekovermass_i_ave_j_xk_ave_manytime(stepsarray, i, j, k, figformat="png", ifpickle=False)

            save_one_plot(
                fig, ax, foldersave,
                "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
                figformat="png", ifpickle=False,
            )
        else:
            
            for step in stepsarray:

                fig, ax = self.plotchunk_ekovermass_i_ave_j_xk_ave_manytime(np.array([step]), i, j, k, figformat="png", ifpickle=False)
                    
                save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)
    
    ##########  plot 1D-1D ekdiff-position index change##########
    def plotchunk_ekminusekaveovermass_i_ave_j_ave(self, step, fig, ax, i, j, k, figformat="png", ifpickle=False):

        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        ekminusekavi = divide_zero(
            2*value_in_a_step_ave(step, "v_Ek" + map_dim_index_to_coordinate[i],self.n_ave, self.lines) - divide_zero(value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)**2, mass),
            mass
            )
        vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], self.n_ave, self.lines)
        mass = np.resize(mass, (n_1, n_2))
        ekminusekavi = np.resize(ekminusekavi, (n_1, n_2))
        vector = np.resize(vector, (n_1, n_2))
        
        sum_along_array_dim = position_index_to_array_dim_index[j]

        ekminusekavi = np.sum(ekminusekavi,axis=sum_along_array_dim)
        mass = np.sum(mass,axis=sum_along_array_dim)
        vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
        ekminusekaveovermass = divide_zero(ekminusekavi, mass)
        
        vector = vector/diameter
        time = time_in_a_step_from_start_rotate(step)
        ax.plot(vector, ekminusekaveovermass, label="{:.2e}".format(time),
                marker = ".",
                linestyle = 'None',
                markersize=16,
                )

        return (fig, ax)
    
    ##########  plot 1D-1D ekdiff-position index change##########    
    def plotchunk_ekminusekaveovermass_i_ave_j_ave_manytime(self, stepsarray, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if ifmanysimu:
            stepsarray.sort()
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(stepsarray, step2, side='right')
                stepsarray_index = stepsarray[:sort_final_index]
                for step in stepsarray_index:
                    fig, ax = self.plotchunk_ekminusekaveovermass_i_ave_j_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                stepsarray = stepsarray[sort_final_index:]
        else:
            for step in stepsarray:
                fig, ax = self.plotchunk_ekminusekaveovermass_i_ave_j_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
        

        ax.legend(
                title = "Time (s)",
                bbox_to_anchor=(1.04,1),
                loc="upper left",
                )
        ax.set_title(self.labelstring_size_walltype_one_line)
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("<V" + map_dim_index_to_coordinate[i] + "^2>" + " - <V" + map_dim_index_to_coordinate[i] + ">^2")
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D ekdiff-position index change##########
    def save_plotchunk_ekminusekaveovermass_i_ave_j_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekminusekaveovermass_i_j_ave_k(i,j,k))
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekminusekaveovermass_i_j_ave_k(i,j,k))
        
        if inonefigure:
            
            fig, ax = self.plotchunk_ekminusekaveovermass_i_ave_j_ave_manytime(stepsarray, i, j, k, figformat="png", ifpickle=False)
            
            save_one_plot(
                fig, ax, foldersave,
                "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
                figformat="png", ifpickle=False,
            )

        else:
            
            for step in stepsarray:

                fig, ax = self.plotchunk_ekminusekaveovermass_i_ave_j_ave_manytime(np.array([step]), i, j, k, figformat="png", ifpickle=False)
                
                save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)

    def plotchunk_1(self, stepsarray, figformat="png", ifpickle=False):
        self.save_v23_x23(stepsarray, figformat="png", ifpickle=False)
        
        self.save_v13_x23(stepsarray, figformat="png", ifpickle=False)
        
        self.save_fraction_x23(stepsarray, figformat="png", ifpickle=False)
        
        for i in [0,1,2]:
            j=1
            k=2
            
            self.save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)
            
            j=2
            k=1
            
            self.save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)
            self.save_plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)
            self.save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
            self.save_plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
            self.save_plotchunk_velocity_i_ave_j_xk_ave(stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
        
        i=0
        j=1
        k=2
        self.save_plotchunk_strain_rate_ij_ave_k_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)
        self.save_plotchunk_strain_rate_ij_ave_k_ave(stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
        self.save_strain_rate_ij_x23(stepsarray, i, j, figformat="png", ifpickle=False)

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

    def plotchunk_velocity_fraction_step1to2(self, step1, step2, figformat="png", ifpickle=False):
        #if smallstep largestep
        # need combine
        self.plotchunk_1_step1to2(save_v23_x23, step1, step2, figformat="png", ifpickle=False)
        self.plotchunk_1_step1to2(save_v13_x23, step1, step2, figformat="png", ifpickle=False)
        self.plotchunk_1_step1to2(save_fraction_x23, step1, step2, figformat="png", ifpickle=False)

    def plotchunk_velocity_fraction_every(self, figformat="png", ifpickle=False):
        #if smallstep largestep
        # need combine
        self.plotchunk_1_every(save_v23_x23, figformat="png", ifpickle=False)
        self.plotchunk_1_every(save_v13_x23, figformat="png", ifpickle=False)
        self.plotchunk_1_every(save_fraction_x23, figformat="png", ifpickle=False)

    ########## get 1D time-1D fraction - 2D x23 ##########
    def manysteparray(self, steparray, v_name):

        n_step = len(steparray)

        initial_array = np.empty([n_line_in_a_step, n_step])
        
        for i, step in enumerate(steparray):
            
            initial_array[:,i] = value_in_a_step_ave(step, v_name, self.n_ave, self.lines)

        return initial_array

    ########## plot 1D volumn fraction-time for every grid ##########
    def check_volumnfraction_x23_increase_decrease_plot(self, steparray, figformat="png", ifpickle=False, ifmanysimu=True):
        
        om.create_directory(dp.f_fraction_check_everygrid)
        tolerence = 0.5

        if ifmanysimu:
            steparray.sort()
            mass = np.empty([n_line_in_a_step, 0])
            time = np.empty(0)
            steparray_new = np.empty(0)
            
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(steparray, step2, side='right')
                steparray_index = steparray[:sort_final_index]
                mass = np.append(mass, manysteparray(steparray_index, "c_m1", self.n_ave, self.lines), axis=-1)
                time = np.append(time, time_in_a_step(steparray_index), axis=-1)
                steparray_new = np.append(steparray_new, steparray_index, axis=-1)
                steparray = steparray[sort_final_index:]

            steparray = steparray_new
        else:
            mass = manysteparray(steparray, "c_m1", self.n_ave, self.lines)
            time = time_in_a_step(steparray)

        volumn_fraction = mass/float(rr.logfile['den'])/vol_in_chunks
        diff_volumn_fraction = np.diff(volumn_fraction, axis=1)
        diff_time = time_in_a_step(np.diff(steparray, axis=0))
        diff_t_volumn_fraction = diff_volumn_fraction/diff_time

        if chunk_method == "rz":
            r = value_in_a_step_ave(steparray[0], "v_r", self.n_ave, self.lines)
            z = value_in_a_step_ave(steparray[0], "v_z", self.n_ave, self.lines)
            x_array = r/diameter
            y_array = z/diameter
        elif chunk_method == "yz":
            y = value_in_a_step_ave(steparray[0], "Coord1", self.n_ave, self.lines)
            z = value_in_a_step_ave(steparray[0], "Coord2", self.n_ave, self.lines)
            x_array = y/diameter
            y_array = z/diameter
        else:
            sys.exit("chunk_method wrong")

        folder_increase = dp.f_fraction_check_everygrid + "increase/"
        om.create_directory(folder_increase)
        folder_decrease = dp.f_fraction_check_everygrid + "decrease/"
        om.create_directory(folder_decrease)
        folder_increase_decrease = dp.f_fraction_check_everygrid + "increase_decrease/"
        om.create_directory(folder_increase_decrease)

        for i in range(n_line_in_a_step):
            if not np.all(diff_t_volumn_fraction[i,:] < 0.1):
                if np.all(diff_t_volumn_fraction[i,:] > (0 - tolerence)) or np.all(diff_t_volumn_fraction[i,:] < (0 + tolerence)):
                        
                    if np.all(diff_t_volumn_fraction[i,:] > (0 - tolerence)) and np.all(diff_t_volumn_fraction[i,:] < (0 + tolerence)):
                        folder = folder_increase_decrease
                    elif np.all(diff_t_volumn_fraction[i,:] > (0 - tolerence)):
                        folder = folder_increase
                    elif np.all(diff_t_volumn_fraction[i,:] < (0 + tolerence)):
                        folder = folder_decrease
                    om.create_directory(folder)
                    if chunk_method == "rz":
                        labelstring = "r_ " + "{:.2e}".format(x_array[i]) + ", z_ " + "{:.2e}".format(y_array[i])
                    elif chunk_method == "yz":
                        labelstring = "y_ " + "{:.2e}".format(x_array[i]) + ", z_ " + "{:.2e}".format(y_array[i])
                    else:
                        sys.exit("chunk_method wrong")

                    fig = plt.figure()
                    ax = fig.add_subplot(111)

                    ax.plot(time, volumn_fraction[i,:], label=labelstring)
                    ax.legend(
                        title="position",
                        bbox_to_anchor=(1.04,1),
                        loc="upper left",
                    )
                    ax.set_xlabel("time (s)")
                    ax.set_ylabel("volumn fraction")
                    plt.tight_layout()

                    fig.savefig(folder + labelstring + "." + figformat, format=figformat)
                    if ifpickle:
                        # Save figure handle to disk
                        with open(folder + labelstring + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                            pickle.dump(fig, f)

                    plt.close('all')

    ########## plot 1D-1D velocity i ave j fix k - time index change##########
    def plotchunk_velocity_i_time_ave_j_fix_k_ave(self, steparray, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True):
        
        steparray.astype(int)

        if ifmanysimu:
            steparray.sort()
            mass = np.empty([n_line_in_a_step, 0])
            vector_mv = np.empty([n_line_in_a_step, 0])
            steparray_new = np.empty(0)
            
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(steparray, step2, side='right')
                steparray_index = steparray[:sort_final_index]
                mass = np.append(mass, manysteparray(steparray_index, "c_m1", self.n_ave, self.lines), axis=-1)
                vector_mv = np.append(vector_mv, manysteparray(steparray_index, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines), axis=-1)
                steparray_new = np.append(steparray_new, steparray_index, axis=-1)
                steparray = steparray[sort_final_index:]

            steparray = steparray_new
        else:
            mass = manysteparray(steparray, "c_m1", self.n_ave, self.lines)
            vector_mv = manysteparray(steparray_index, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)

        vector_vi_sum_array = np.empty(0)
        time_array = np.empty(0)

        for nn, step in enumerate(steparray):
            
            vector_mv_step = np.resize(vector_mv[:,nn], (n_1, n_2))
            mass_step = np.resize(mass[:,nn], (n_1, n_2))
            vector_vi = divide_zero(vector_mv_step, mass_step)
            fix_along_array_dim = position_index_to_array_dim_index[k]
            ave_along_array_dim = position_index_to_array_dim_index[j]
                    
            count_not_empty = np.sum(if_grid_surround_not_empty(mass_step), axis=ave_along_array_dim)
            vector_vi_sum = divide_zero(np.take(np.sum(vector_vi,axis=ave_along_array_dim), k_index, axis=fix_along_array_dim), np.take(count_not_empty, k_index, axis=fix_along_array_dim))

            time = time_in_a_step_from_start_rotate(step)
            vector_vi_sum_array= np.append(vector_vi_sum_array, vector_vi_sum)
            time_array= np.append(time_array, time)
        
        vector_vi_sum_array = vector_vi_sum_array/velocity_scale

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(time_array, vector_vi_sum_array, 
                label=self.labelstring_size_walltype,
                )
        plt.xticks(rotation=45)

        ax.legend(
            title="",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )

        return (fig, ax)########## plot 1D-1D omega i ave j fix k - time index change##########
    def plotchunk_omega_i_time_ave_j_fix_k_ave(self, steparray, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True):
        steparray.astype(int)

        if ifmanysimu:
            steparray.sort()
            mass = np.empty([n_line_in_a_step, 0])
            vector_mv = np.empty([n_line_in_a_step, 0])
            steparray_new = np.empty(0)
            
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(steparray, step2, side='right')
                steparray_index = steparray[:sort_final_index]

                mass = np.append(mass, manysteparray(steparray_index, "c_m1", self.n_ave, self.lines), axis=-1)
                vector_mv = np.append(vector_mv, manysteparray(steparray_index, "c_omega[1]", self.n_ave, self.lines), axis=-1)
                steparray_new = np.append(steparray_new, steparray_index, axis=-1)
                steparray = steparray[sort_final_index:]

            steparray = steparray_new
        else:
            mass = manysteparray(steparray, "c_m1", self.n_ave, self.lines)
            vector_mv = manysteparray(steparray_index, "c_omega[1]" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)

        vector_vi_sum_array = np.empty(0)
        time_array = np.empty(0)

        for nn, step in enumerate(steparray):
            
            vector_mv_step = np.resize(vector_mv[:,nn], (n_1, n_2))
            mass_step = np.resize(mass[:,nn], (n_1, n_2))
            vector_vi = divide_zero(vector_mv_step, mass_step)
            fix_along_array_dim = position_index_to_array_dim_index[k]
            ave_along_array_dim = position_index_to_array_dim_index[j]
                    
            count_not_empty = np.sum(if_grid_surround_not_empty(mass_step), axis=ave_along_array_dim)
            vector_vi_sum = np.take(np.sum(vector_vi,axis=ave_along_array_dim), k_index, axis=fix_along_array_dim)

            time = time_in_a_step_from_start_rotate(step)
            vector_vi_sum_array= np.append(vector_vi_sum_array, vector_vi_sum)
            time_array= np.append(time_array, time)
        
        vector_vi_sum_array = vector_vi_sum_array/velocity_scale

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(time_array, vector_vi_sum_array, 
                label=self.labelstring_size_walltype,
                )
        plt.xticks(rotation=45)

        ax.legend(
            title="",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )

        return (fig, ax)
    


    def save_plotchunk_velocity_i_time_ave_j_fix_k_ave(self, steparray, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True):

        om.create_directory(dp.f_velocity_i_time_ave_j_fix_k_ave_path)

        folder = dp.f_velocity_i_time_ave_j_fix_k_ave_path + "V_" + str(i) + "_ave_" + str(j) + "_fix_" + str(k) + "/"
        om.create_directory(folder)

        add_nve_subfolder_in_folder(self.n_ave, folder)
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, folder)
            
        fig, ax = self.plotchunk_velocity_i_time_ave_j_fix_k_ave(steparray, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True)
        
        ax.set_xlabel("time (s)")
        ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (Vwall)" + "_at " + map_dim_index_to_coordinate[k] + "=" + str(k_index))
        plt.tight_layout()

        figurepath_no_extension = foldersave + "combine" + str(int(steparray[0])) + "_to_" + str(int(steparray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')
    
    def save_plotchunk_omega_i_time_ave_j_fix_k_ave(self, steparray, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True):
        f_create = dp.diagram_path + "omega1/"
        om.create_directory(f_create)

        folder = f_create + "omega_" + 0 + "_ave_" + str(j) + "_fix_" + str(k) + "/"
        om.create_directory(folder)

        add_nve_subfolder_in_folder(self.n_ave, folder)
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, folder)
            
        fig, ax = self.plotchunk_omega_i_time_ave_j_fix_k_ave(steparray, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True)
        
        ax.set_xlabel("time (s)")
        ax.set_ylabel("Omega_x" + " (Vwall)" + "_at " + map_dim_index_to_coordinate[k] + "=" + str(k_index))
        plt.tight_layout()

        figurepath_no_extension = foldersave + "combine" + str(int(steparray[0])) + "_to_" + str(int(steparray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

    def save_plotchunk_velocity_i_time_near_wall_ave(self, steparray, i, figformat="png", ifpickle=False, ifmanysimu=True):

        om.create_directory(dp.f_velocity_i_time_ave_j_fix_k_ave_path)

        folder = dp.f_velocity_i_time_ave_j_fix_k_ave_path + "V_" + str(i) + "_near_wall" + "/"
        om.create_directory(folder)

        add_nve_subfolder_in_folder(self.n_ave, folder)
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, folder)
            
        fig, ax = self.plotchunk_velocity_i_time_ave_j_fix_k_ave(steparray, i, 2, 1, 0, figformat="png", ifpickle=False, ifmanysimu=True)

        ax.set_xlabel("time (s)")
        ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (V_wall)")
        plt.tight_layout()

        figurepath_no_extension = foldersave + "combine" + str(int(steparray[0])) + "_to_" + str(int(steparray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

    def save_plotchunk_omega_i_time_near_wall_ave(self, steparray, i, figformat="png", ifpickle=False, ifmanysimu=True):
        om.create_directory(dp.diagram_path + "omega_0_near_wall/")
        folder = dp.diagram_path + "omega_0_near_wall/" + "V_" + str(i) + "_near_wall" + "/"
        om.create_directory(folder)

        add_nve_subfolder_in_folder(self.n_ave, folder)
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, folder)
            
        fig, ax = self.plotchunk_omega_i_time_ave_j_fix_k_ave(steparray, i, 2, 1, 0, figformat="png", ifpickle=False, ifmanysimu=True)

        ax.set_xlabel("time (s)")
        ax.set_ylabel("Omega_x" + " (V_wall)")
        plt.tight_layout()

        figurepath_no_extension = foldersave + "combine" + str(int(steparray[0])) + "_to_" + str(int(steparray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')
        
    ########## plot 1D-1D fraction-ave over j position k index change##########
    def plotchunk_fraction_ave_j_ave(self, step, fig, ax, j, k, figformat="png", ifpickle=False, gridwidth=2):
        
        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], self.n_ave, self.lines)
        mass = np.resize(mass, (n_1, n_2))
        vector = np.resize(vector, (n_1, n_2))
        new_mass = 0
        new_vector = 0

        not_surroundemptyorselfempty = np.full((mass.shape[0]-gridwidth+1, mass.shape[1]), True)
        
        for nnn in range(gridwidth):
            new_mass += mass[nnn:mass.shape[0]-gridwidth+1+nnn, :]
            new_vector += vector[nnn:mass.shape[0]-gridwidth+1+nnn, :]
            not_surroundemptyorselfempty = np.logical_and(
                not_surroundemptyorselfempty,
                if_grid_surround_not_empty(mass[nnn:mass.shape[0]-gridwidth+1+nnn, :]),
                )
        mass = new_mass/gridwidth
        vector = new_vector/gridwidth
        fraction = mass/float(rr.logfile['den'])/vol_in_chunks
        
        sum_along_array_dim = position_index_to_array_dim_index[j]

        count_not_surroundemptyorselfempty = np.sum(not_surroundemptyorselfempty, axis=sum_along_array_dim)
        fraction_only_not_surround_empty = fraction*not_surroundemptyorselfempty

        fraction_ave = divide_zero(np.sum(fraction_only_not_surround_empty, axis=sum_along_array_dim), count_not_surroundemptyorselfempty)

        vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
        
        vector = vector/diameter

        time = time_in_a_step_from_start_rotate(step)

        ax.plot(vector, fraction_ave, label="{:.2e}".format(time))

        return (fig, ax)

    ########## plot 1D-1D fraction-ave over j position k index change##########
    def plotchunk_fraction_ave_j_ave_manytime(self, stepsarray, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if ifmanysimu:
            stepsarray.sort()
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(stepsarray, step2, side='right')
                stepsarray_index = stepsarray[:sort_final_index]
                for step in stepsarray_index:
                    fig, ax = self.plotchunk_fraction_ave_j_ave(step, fig, ax, j, k, figformat="png", ifpickle=False)
                
                stepsarray = stepsarray[sort_final_index:]
        else:
            for step in stepsarray:
                fig, ax = self.plotchunk_fraction_ave_j_ave(step, fig, ax, j, k, figformat="png", ifpickle=False)
        
        ax.legend(
            title = "time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("Volumn Fraction")
        ax.set_title(self.labelstring_size_walltype_one_line)
        plt.tight_layout()

        return (fig, ax)

    ########## save plot fraction-ave over j position k index change##########
    def save_plotchunk_fraction_ave_j_ave(self, stepsarray, j, k, figformat="png", ifpickle=False, inonefigure=False):
        
        f_path_fraction = dp.diagram_path + "fraction/"
        om.create_directory(f_path_fraction)

        f_path_ave_jk = f_path_fraction + "ave_" + str(j) + "x" + str(k) + "/"

        add_nve_subfolder_in_folder(self.n_ave, f_path_ave_jk)
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, f_path_ave_jk)
        
        if inonefigure:
            
            fig, ax = self.plotchunk_fraction_ave_j_ave_manytime(stepsarray, j, k, figformat="png", ifpickle=False)

            save_one_plot(
                fig, ax, foldersave,
                "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
                figformat="png", ifpickle=False,
            )
        else:
            
            for step in stepsarray:

                fig, ax = self.plotchunk_fraction_ave_j_ave_manytime(np.array([step]), j, k, figformat="png", ifpickle=False)
                    
                save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)
    
        plt.close('all')

    ########## plot 1D-1D velocity(no top)-position index change##########
    def plotchunk_velocity_i_ave_j_xk_ave_no_top(self, step, fig, ax, i, j, k, figformat="png", ifpickle=False):

        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        mvi = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)
        vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], self.n_ave, self.lines)
        mass = np.resize(mass, (n_1, n_2))
        mvi = np.resize(mvi, (n_1, n_2))
        vector = np.resize(vector, (n_1, n_2))

        mvi    = mvi*if_grid_surround_not_empty(mass)
        mass   = mass*if_grid_surround_not_empty(mass)

        sum_along_array_dim = position_index_to_array_dim_index[j]
        mvi = np.sum(mvi,axis=sum_along_array_dim)
        mass = np.sum(mass,axis=sum_along_array_dim)
        vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
        velocity = divide_zero(mvi, mass)
        velocity /= velocity_scale
        vector = vector/diameter
        time = time_in_a_step_from_start_rotate(step)
        ax.plot(vector, velocity, label="{:.2e}".format(time),
                marker = ".",
                linestyle = 'None',
                markersize=16,
                )

        return (fig, ax)
    
    ########## plot 1D-1D velocity(no top)-position index change##########
    def plotchunk_velocity_i_ave_j_xk_ave_no_top_manytime(self, stepsarray, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if ifmanysimu:
            stepsarray.sort()
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(stepsarray, step2, side='right')
                stepsarray_index = stepsarray[:sort_final_index]
                for step in stepsarray_index:
                    fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave_no_top(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                
                stepsarray = stepsarray[sort_final_index:]
        else:
            for step in stepsarray:
                fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave_no_top(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
        
        ax.legend(
            title = "Time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_xlabel(map_dim_index_to_coordinate[k] + " (d)")
        ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (V_wall)")
        ax.set_title(self.labelstring_size_walltype_one_line)
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D velocity(no top)-position index change##########
    def save_plotchunk_velocity_i_ave_j_xk_ave_no_top(self, stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_velocity_i_j_ave_k_no_top(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_velocity_i_j_ave_k_no_top(i,j,k))
        
        if inonefigure:
            
            fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave_no_top_manytime(stepsarray, i, j, k, figformat="png", ifpickle=False)

            save_one_plot(
                fig, ax, foldersave,
                "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
                figformat="png", ifpickle=False,
            )

        else:
            
            for step in stepsarray:

                fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave_no_top_manytime(np.array([step]), i, j, k, figformat="png", ifpickle=False)
                    
                save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)

    ########## plot 1D-1D ek(no top)-position index change##########
    def plotchunk_ek_i_ave_j_xk_ave_no_top(self, step, fig, ax, i, j, k, figformat="png", ifpickle=False, ifbottom=True):

        mass = value_in_a_step_ave(step, "c_m1", self.n_ave, self.lines)
        eki = value_in_a_step_ave(step, "v_Ek" + map_dim_index_to_coordinate[i], self.n_ave, self.lines)
        vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], self.n_ave, self.lines)
        mass = np.resize(mass, (n_1, n_2))
        eki = np.resize(eki, (n_1, n_2))
        vector = np.resize(vector, (n_1, n_2))
        if not ifbottom:
            eki    = eki*if_grid_surround_not_empty(mass)
            mass   = mass*if_grid_surround_not_empty(mass)
        if ifbottom:
            mass = mass[:, 0:5]
            eki = eki[:, 0:5]
            vector = vector[:, 0:5]


        sum_along_array_dim = position_index_to_array_dim_index[j]
        eki = np.sum(eki,axis=sum_along_array_dim)
        mass = np.sum(mass,axis=sum_along_array_dim)
        vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
        ekiovermass = divide_zero(eki, mass)

        ekiovermass /= velocity_scale**2
        vector = vector/diameter
        time = time_in_a_step_from_start_rotate(step)
        ax.plot(vector, ekiovermass, label="{:.2e}".format(time))

        return (fig, ax)
    
    ########## plot 1D-1D ek(no top)-position index change##########
    def plotchunk_ek_i_ave_j_xk_ave_no_top_manytime(self, stepsarray, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if ifmanysimu:
            stepsarray.sort()
            for index in range(rr.n_loglist):
                self.lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
                step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
                sort_final_index = np.searchsorted(stepsarray, step2, side='right')
                stepsarray_index = stepsarray[:sort_final_index]
                for step in stepsarray_index:
                    fig, ax = self.plotchunk_ek_i_ave_j_xk_ave_no_top(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                
                stepsarray = stepsarray[sort_final_index:]
        else:
            for step in stepsarray:
                fig, ax = self.plotchunk_ek_i_ave_j_xk_ave_no_top(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
        
        ax.legend(
            title = "Time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_xlabel(map_dim_index_to_coordinate[k] + " (d)")
        ax.set_ylabel("ek_" + map_dim_index_to_coordinate[i] + "/mass_" + map_dim_index_to_coordinate[j] + " (m^2/s^2)")
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D ek(no top)-position index change##########
    def save_plotchunk_ek_i_ave_j_xk_ave_no_top(self, stepsarray, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_ek_i_j_ave_k_no_top(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_ek_i_j_ave_k_no_top(i,j,k))
        
        if inonefigure:
            
            fig, ax = self.plotchunk_ek_i_ave_j_xk_ave_no_top_manytime(stepsarray, i, j, k, figformat="png", ifpickle=False)

            save_one_plot(
                fig, ax, foldersave,
                "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
                figformat="png", ifpickle=False,
            )

        else:
            
            for step in stepsarray:

                fig, ax = self.plotchunk_ek_i_ave_j_xk_ave_no_top_manytime(np.array([step]), i, j, k, figformat="png", ifpickle=False)
                    
                save_one_plot(fig, ax, foldersave, str(int(step)), figformat="png", ifpickle=False)


class chunk_include_pre(object):

    def __init__(self, n_ave, lmp_path):

        self.n_ave = n_ave
        self.lmp_path = lmp_path
        self.log_file_list_initial_to_last = rr.log_current_plus_previousfrom_initial_to_lastest(self.lmp_path)
        self.n_log_list = len(self.log_file_list_initial_to_last)
        
        first_rotate_index = 0
        for logfile in self.log_file_list_initial_to_last:
            if logfile["ifrotate"] == "yes" and float(logfile["Sa"]) != 0:
                break
            first_rotate_index += 1
        
        self.first_rotate_index = first_rotate_index
        
        self.folder_path_list_initial_to_last = rc.folder_path_list_initial_to_last
        self.first_rotate_step = chunk(self.n_ave, self.folder_path_list_initial_to_last[self.first_rotate_index]).step_first_in_file_change_by_n_ave
        self.end_step = chunk(self.n_ave, self.folder_path_list_initial_to_last[-1]).step_last_in_file_change_by_n_ave
        allsteps_since_rotate = np.empty(0)
        first_middle_last_steps_everysimu = np.empty(0)
        first_extra_middle_last_steps_everysimu = np.empty(0)
        for index in range(self.first_rotate_index, self.n_log_list):
            chunk_simu = chunk(self.n_ave, self.folder_path_list_initial_to_last[index])
            allsteps_since_rotate = np.append(allsteps_since_rotate, chunk_simu.allsteps)
            first_middle_last_steps_everysimu = np.append(
                first_middle_last_steps_everysimu, chunk_simu.first_middle_last_steps
                )
            first_extra_middle_last_steps_everysimu = np.append(
                first_extra_middle_last_steps_everysimu, chunk_simu.first_extra_middle_last_steps
                )
        self.allsteps_since_rotate = allsteps_since_rotate
        self.first_middle_last_steps_everysimu = first_middle_last_steps_everysimu
        self.first_extra_middle_last_steps_everysimu = first_extra_middle_last_steps_everysimu
    def current_and_remain_stepsarray(self, stepsarray, index):

        step2 = step_last_fix_change_by_n_ave(self.n_ave, index)
        sort_final_index = np.searchsorted(stepsarray, step2, side='right')
        stepsarray_current = stepsarray[:sort_final_index]
        stepsarray_remain = stepsarray[sort_final_index:]

        return [stepsarray_current, stepsarray_remain]

     
    ########## plot 1D-1D strain_rate-position index change##########
    def plotchunk_strain_rate_ij_ave_k_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        stepsarray_remain = stepsarray.copy()
        stepsarray_remain.sort()

        for index, lmp_path in enumerate(self.folder_path_list_initial_to_last):

            [stepsarray_current, stepsarray_remain] = self.current_and_remain_stepsarray(stepsarray_remain, index)

            chunkobject = chunk(self.n_ave, lmp_path)
            for step in stepsarray_current:
                chunkobject.plotchunk_strain_rate_ij_ave_k_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
            
        ax.legend(
                title = "Time (s)",
                bbox_to_anchor=(1.04,1), 
                loc="upper left"
                )
        ax.set_xlabel(map_dim_index_to_coordinate[j])
        ax.set_ylabel("strain_rate" + map_dim_index_to_coordinate[i] + map_dim_index_to_coordinate[j] + "_average_" + map_dim_index_to_coordinate[k])
        plt.tight_layout()

        return (fig, ax)
    ########## save plot 1D-1D strain_rate-position index change##########
    def save_plotchunk_strain_rate_ij_ave_k_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):
        
        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_strain_rate_i_j_ave_k(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_strain_rate_i_j_ave_k(i,j,k))
        
        fig, ax = self.plotchunk_strain_rate_ij_ave_k_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)
        
        save_one_plot(
            fig, ax, foldersave,
            "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
            figformat="png", ifpickle=False,
        )

        plt.close('all')

    ########## plot 1D-1D velocity-position index change##########
    def plotchunk_velocity_i_ave_j_xk_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        stepsarray_remain = stepsarray.copy()
        stepsarray_remain.sort()

        for index, lmp_path in enumerate(self.folder_path_list_initial_to_last):

            [stepsarray_current, stepsarray_remain] = self.current_and_remain_stepsarray(stepsarray_remain, index)

            chunkobject = chunk(self.n_ave, lmp_path)
            for step in stepsarray_current:
                chunkobject.plotchunk_velocity_i_ave_j_xk_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
  
        ax.legend(
            title = "Time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (Vwall)")
        ax.set_title(chunkobject.labelstring_size_walltype_one_line)
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D velocity-position index change##########
    def save_plotchunk_velocity_i_ave_j_xk_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_velocity_i_j_ave_k(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_velocity_i_j_ave_k(i,j,k))

        fig, ax = self.plotchunk_velocity_i_ave_j_xk_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)

        save_one_plot(
            fig, ax, foldersave,
            "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
            figformat="png", ifpickle=False,
        )     
    
    ########## plot 1D-1D ek-position index change##########
    def plotchunk_ekovermass_i_ave_j_xk_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        stepsarray_remain = stepsarray.copy()
        stepsarray_remain.sort()

        for index, lmp_path in enumerate(self.folder_path_list_initial_to_last):

            [stepsarray_current, stepsarray_remain] = self.current_and_remain_stepsarray(stepsarray_remain, index)

            chunkobject = chunk(self.n_ave, lmp_path)
            for step in stepsarray_current:
                chunkobject.plotchunk_ekovermass_i_ave_j_xk_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                
        ax.legend(
            title = "Time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("<E" + map_dim_index_to_coordinate[i] + ">" + "/" + "<m>")
        ax.set_title(chunkobject.labelstring_size_walltype_one_line)
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D ek-position index change##########
    def save_plotchunk_ekovermass_i_ave_j_xk_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekovermass_i_j_ave_k(i,j,k))
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekovermass_i_j_ave_k(i,j,k))
           
        fig, ax = self.plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)

        save_one_plot(
            fig, ax, foldersave,
            "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
            figformat="png", ifpickle=False,
        )
        
    
    ##########  plot 1D-1D ekdiff-position index change##########    
    def plotchunk_ekminusekaveovermass_i_ave_j_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        stepsarray_remain = stepsarray.copy()
        stepsarray_remain.sort()

        for index, lmp_path in enumerate(self.folder_path_list_initial_to_last):

            [stepsarray_current, stepsarray_remain] = self.current_and_remain_stepsarray(stepsarray_remain, index)

            chunkobject = chunk(self.n_ave, lmp_path)
            for step in stepsarray_current:
                chunkobject.plotchunk_ekminusekaveovermass_i_ave_j_ave(step, fig, ax, i, j, k, figformat="png", ifpickle=False)
                
        ax.legend(
                title = "Time (s)",
                bbox_to_anchor=(1.04,1),
                loc="upper left",
                )
        ax.set_title(chunkobject.labelstring_size_walltype_one_line)
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("<V" + map_dim_index_to_coordinate[i] + "^2>" + " - <V" + map_dim_index_to_coordinate[i] + ">^2")
        plt.tight_layout()

        return (fig, ax)

    ########## save plot 1D-1D ekdiff-position index change##########
    def save_plotchunk_ekminusekaveovermass_i_ave_j_ave(self, stepsarray, i, j, k, figformat="png", ifpickle=False):

        add_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekminusekaveovermass_i_j_ave_k(i,j,k))
        foldersave = path_nve_subfolder_in_folder(self.n_ave, dp.f_path_ekminusekaveovermass_i_j_ave_k(i,j,k))

        fig, ax = self.plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray, i, j, k, figformat="png", ifpickle=False)
        
        save_one_plot(
            fig, ax, foldersave,
            "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
            figformat="png", ifpickle=False,
        )
 
    ########## plot 1D-1D fraction-ave over j position k index change##########
    def plotchunk_fraction_ave_j_ave(self, stepsarray, j, k, figformat="png", ifpickle=False):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        stepsarray_remain = stepsarray.copy()
        stepsarray_remain.sort()
        for index, lmp_path in enumerate(self.folder_path_list_initial_to_last):

            [stepsarray_current, stepsarray_remain] = self.current_and_remain_stepsarray(stepsarray_remain, index)

            chunkobject = chunk(self.n_ave, lmp_path)
            for step in stepsarray_current:
                chunkobject.plotchunk_fraction_ave_j_ave(step, fig, ax, j, k, figformat="png", ifpickle=False)
                
        ax.legend(
            title = "time (s)",
            bbox_to_anchor=(1.04,1),
            loc="upper left",
            )
        ax.set_xlabel(map_dim_index_to_coordinate[k])
        ax.set_ylabel("Volumn Fraction")
        ax.set_title(chunkobject.labelstring_size_walltype_one_line)
        plt.tight_layout()

        return (fig, ax)

    ########## save plot fraction-ave over j position k index change##########
    def save_plotchunk_fraction_ave_j_ave(self, stepsarray, j, k, figformat="png", ifpickle=False):
        
        f_path_fraction = dp.diagram_path + "fraction/"
        om.create_directory(f_path_fraction)

        f_path_ave_jk = f_path_fraction + "ave_" + str(j) + "x" + str(k) + "/"

        add_nve_subfolder_in_folder(self.n_ave, f_path_ave_jk)
        
        foldersave = path_nve_subfolder_in_folder(self.n_ave, f_path_ave_jk)
        fig, ax = self.plotchunk_fraction_ave_j_ave(stepsarray, j, k, figformat="png", ifpickle=False)

        save_one_plot(
            fig, ax, foldersave,
            "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1])),
            figformat="png", ifpickle=False,
        )

    def plotchunk_all(self):
        stepsarray = self.first_extra_middle_last_steps_everysimu
        self.save_plotchunk_strain_rate_ij_ave_k_ave(stepsarray, 0, 1, 2, figformat="png", ifpickle=False)
        for i in range(3):
            self.save_plotchunk_velocity_i_ave_j_xk_ave(stepsarray, i, 2, 1, figformat="png", ifpickle=False)
            self.save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, i, 2, 1, figformat="png", ifpickle=False)
            self.save_plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray, i, 2, 1, figformat="png", ifpickle=False)
        self.save_plotchunk_fraction_ave_j_ave(stepsarray, 2, 1, figformat="png", ifpickle=False)
        



## class before ##


########## plot 2D-2D velocity-position ##########
def plotchunk_ave_one_step_v23x23(step, n_ave, lines, figformat="png", ifpickle=False):
    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    if chunk_method == "rz":
        r = value_in_a_step_ave(step, "v_r", n_ave, lines)
        z = value_in_a_step_ave(step, "v_z", n_ave, lines)
        x_array = r/diameter
        y_array = z/diameter
    elif chunk_method == "yz":
        y = value_in_a_step_ave(step, "Coord1", n_ave, lines)
        z = value_in_a_step_ave(step, "Coord2", n_ave, lines)
        x_array = y/diameter
        y_array = z/diameter
    else:
        sys.exit("chunk_method wrong")
    if chunk_method == "rz":
        mvr = value_in_a_step_ave(step, "v_mvr", n_ave, lines)
        mvz = value_in_a_step_ave(step, "v_mvz", n_ave, lines)
        vx_array = divide_zero(mvr,mass)/velocity_scale
        vy_array = divide_zero(mvz,mass)/velocity_scale
    elif chunk_method == "yz":
        mvy = value_in_a_step_ave(step, "v_mvy", n_ave, lines)
        mvz = value_in_a_step_ave(step, "v_mvz", n_ave, lines)
        vx_array = divide_zero(mvy,mass)/velocity_scale
        vy_array = divide_zero(mvz,mass)/velocity_scale
    else:
        sys.exit("chunk_method wrong")
    v_length_array = (vx_array**2+vy_array**2)**0.5
    max_v_length = np.amax(v_length_array)
    if dp.if_plot_velocity_field_scale_same == "yes":
        quiver_scale = dp.quiver_scale_velocity_xaxis_shearplanenormal_yaxis_z
    else:
        quiver_scale = max_v_length/2
    label_scale = quiver_scale

    fig1, ax1 = plt.subplots()
    plot_quiver_position_label(fig1, ax1)
        
    Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                )
    time = time_in_a_step_from_start_rotate(step)
    ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                " : {:.2e} wall velocity in 45 degree".format(label_scale) + ". At {:.2e} s".format(time),
                labelpos='E', coordinates='figure', angle=45)
    return (fig1, ax1)    

########## plot 2D-2D velocity-position ##########
def plotchunk_ave_one_step_v13x23(step, n_ave, lines,figformat="png", ifpickle=False):
    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    if chunk_method == "rz":
        r = value_in_a_step_ave(step, "v_r", n_ave, lines)
        z = value_in_a_step_ave(step, "v_z", n_ave, lines)
        x_array = r/diameter
        y_array = z/diameter
    elif chunk_method == "yz":
        y = value_in_a_step_ave(step, "Coord1", n_ave, lines)
        z = value_in_a_step_ave(step, "Coord2", n_ave, lines)
        x_array = y/diameter
        y_array = z/diameter
    else:
        sys.exit("chunk_method wrong")

    if chunk_method == "rz":
        mvt = value_in_a_step_ave(step, "v_mvt", n_ave, lines)
        mvz = value_in_a_step_ave(step, "v_mvz", n_ave, lines)
        vx_array = divide_zero(mvt,mass)/velocity_scale
        vy_array = divide_zero(mvz,mass)/velocity_scale
    elif chunk_method == "yz":
        mvx = value_in_a_step_ave(step, "v_mvx", n_ave, lines)
        mvz = value_in_a_step_ave(step, "v_mvz", n_ave, lines)
        vx_array = divide_zero(mvx,mass)/velocity_scale
        vy_array = divide_zero(mvz,mass)/velocity_scale
    else:
        sys.exit("chunk_method wrong")
    
    v_length_array = (vx_array**2+vy_array**2)**0.5
    max_v_length = np.amax(v_length_array)

    if dp.if_plot_velocity_field_scale_same == "yes":
        quiver_scale = dp.quiver_scale_velocity_xaxis_shearplanenormal_yaxis_z
    else:
        quiver_scale = max_v_length/2
    label_scale = quiver_scale

    fig1, ax1 = plt.subplots()
    plot_quiver_position_label(fig1, ax1)
        
    Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                )
    time = time_in_a_step_from_start_rotate(step)
    ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                label = " : {:.2e} wall velocity in 45 degree".format(label_scale) + ". At {:.2e} s".format(time),
                labelpos='E', coordinates='figure', angle=45)
    return (fig1, ax1)

########## plot 1D-2D volumn fraction-position ##########
def plotchunk_ave_one_step_volumnfraction_x23(step, n_ave, lines, figformat="png", ifpickle=False):
    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    if chunk_method == "rz":
        r = value_in_a_step_ave(step, "v_r", n_ave, lines)
        z = value_in_a_step_ave(step, "v_z", n_ave, lines)
        x_array = r/diameter
        y_array = z/diameter
    elif chunk_method == "yz":
        y = value_in_a_step_ave(step, "Coord1", n_ave, lines)
        z = value_in_a_step_ave(step, "Coord2", n_ave, lines)
        x_array = y/diameter
        y_array = z/diameter
    else:
        sys.exit("chunk_method wrong")
    quiver_scale = 0.2
    label_scale = 0.6
    vy_array = mass/float(rr.logfile['den'])/vol_in_chunks
    fig1, ax1 = plt.subplots()
    
    #fig1.figsize = [12.8, 9.6]
    plot_quiver_position_label(fig1, ax1)
    #ax1.set_title('velocity field r-z direction (average over theta)')
    Q = ax1.quiver(x_array, y_array, np.zeros_like(vy_array), vy_array,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                )
    time = time_in_a_step_from_start_rotate(step)
    ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                label = "equal solid fraction = {:.2e}".format(label_scale) + ". At {:.2e} s".format(time),
                labelpos='E',
                coordinates='figure', angle=90)
    return (fig1, ax1)

########## save plot 2D-2D velocity-position ##########
def save_v23_x23(stepsarray, n_ave, lines, figformat="png", ifpickle=False):
    add_nve_subfolder_in_folder(n_ave, dp.f_momentum_mass_field_v23x23_path)
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_momentum_mass_field_v23x23_path)
    for step in stepsarray:
        fig, ax = plotchunk_ave_one_step_v23x23(step, n_ave, lines, figformat="png", ifpickle=False)
        fig.savefig(foldersave + str(int(step)) + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(foldersave + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

def save_v13_x23(stepsarray, n_ave, lines, figformat="png", ifpickle=False):
    add_nve_subfolder_in_folder(n_ave, dp.f_momentum_mass_field_v13x23_path)
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_momentum_mass_field_v13x23_path)
    for step in stepsarray:
        fig, ax = plotchunk_ave_one_step_v13x23(step, n_ave, lines, figformat="png", ifpickle=False)
        fig.savefig(foldersave + str(int(step)) + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(foldersave + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')
########## save plot 1D-2D volumn fraction-position ##########
def save_fraction_x23(stepsarray, n_ave, lines, figformat="png", ifpickle=False):
    add_nve_subfolder_in_folder(n_ave, dp.f_momentum_mass_field_volumnfraction_x23_path)
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_momentum_mass_field_volumnfraction_x23_path)
    for step in stepsarray:
        fig, ax = plotchunk_ave_one_step_volumnfraction_x23(step, n_ave, lines, figformat="png", ifpickle=False)
        fig.savefig(foldersave + str(int(step)) + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(foldersave + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

########## plot 1D-2D strain_rate-position ##########
def plotchunk_strain_rate_i_j_x23(step, n_ave, lines, i, j, figformat="png", ifpickle=False):
    
    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    
    vector_mv = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines)
    vector_vi = divide_zero(vector_mv, mass)
    vector_x2 = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[1]], n_ave, lines)
    vector_x3 = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[2]], n_ave, lines)
    
    if j == 1:
        vector_xi = vector_x2
    elif j == 2:
        vector_xi = vector_x3

    vector_vi = np.resize(vector_vi, (n_1, n_2))
    vector_xi = np.resize(vector_xi, (n_1, n_2))
    vector_x2 = np.resize(vector_x2, (n_1, n_2))
    vector_x3 = np.resize(vector_x3, (n_1, n_2))

    diff_along_array_dim = position_index_to_array_dim_index[j]
    if diff_along_array_dim == 0:
        middle_point_vector_x2 = (vector_x2[:-1,:] + vector_x2[1:,:])/2
        middle_point_vector_x3 = (vector_x3[:-1,:] + vector_x3[1:,:])/2
    elif diff_along_array_dim == 1:
        middle_point_vector_x2 = (vector_x2[:,:-1] + vector_x2[:,1:])/2
        middle_point_vector_x3 = (vector_x3[:,:-1] + vector_x3[:,1:])/2
    else:
        sys.exit("error")
    
    vector_vi_diff = np.diff(vector_vi,axis=diff_along_array_dim)
    vector_xi_diff = np.diff(vector_xi,axis=diff_along_array_dim)
   
    strain_rate = divide_zero(vector_vi_diff,vector_xi_diff)
    strain_rate /= shear_rate_scale
    
    middle_point_vector_x2 /= diameter
    middle_point_vector_x3 /= diameter

    quiver_scale = 0.5
    label_scale = quiver_scale
    fig1, ax1 = plt.subplots()
    
    plot_quiver_position_label(fig1, ax1)
    Q = ax1.quiver(middle_point_vector_x2, middle_point_vector_x3, np.zeros_like(vector_vi_diff), strain_rate,
                units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                )
    
    time = time_in_a_step_from_start_rotate(step)

    labelstring = (
        " : strain_rate_"
        + map_dim_index_to_coordinate[i]
        + map_dim_index_to_coordinate[j]
        + " = {:.2e}".format(label_scale)
        + ". At {:.2e} s".format(time)
    )
    
    ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                label = labelstring,
                labelpos='E',
                coordinates='figure', angle=90)
    return (fig1, ax1)
    
########## save plot 1D-2D strain_rate-position ##########
def save_strain_rate_ij_x23(stepsarray, n_ave, lines, i, j, figformat="png", ifpickle=False):
    add_nve_subfolder_in_folder(n_ave, dp.f_path_strain_rate_i_j_x23(i,j))
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_path_strain_rate_i_j_x23(i,j))
    for step in stepsarray:
        fig, ax = plotchunk_strain_rate_i_j_x23(step, n_ave, lines, i, j, figformat="png", ifpickle=False)
        fig.savefig(foldersave + str(int(step)) + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(foldersave + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')



########## plot 1D-1D strain_rate-position index change##########
def plotchunk_strain_rate_ij_ave_k_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False):

    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    mass = np.resize(mass, (n_1, n_2))
    
    vector_mv = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines)
    vector_mv = np.resize(vector_mv, (n_1, n_2))
    vector_vi = divide_zero(vector_mv, mass)
    vector_j = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[j]], n_ave, lines)

    vector_j = np.resize(vector_j, (n_1, n_2))

    sum_along_array_dim = position_index_to_array_dim_index[k]

    diff_along_array_dim = position_index_to_array_dim_index[j]
    if diff_along_array_dim == 0:
        middle_point_vector_j = (vector_j[:-1,0] + vector_j[1:,0])/2
    elif diff_along_array_dim == 1:
        middle_point_vector_j = (vector_j[0,:-1] + vector_j[0,1:])/2
    else:
        sys.exit("error")

    len_in_diff_dim = mass.shape[diff_along_array_dim]
    diff_not_empty = np.logical_and(
        np.take(if_grid_surround_not_empty(mass), np.arange(len_in_diff_dim-1), axis=diff_along_array_dim), 
        np.take(if_grid_surround_not_empty(mass), np.arange(len_in_diff_dim-1)+1, axis=diff_along_array_dim),
        )
    count_not_empty = np.sum(diff_not_empty, axis=sum_along_array_dim)
    
    vector_vi_diff = np.sum(np.diff(vector_vi,axis=diff_along_array_dim)*diff_not_empty, axis=sum_along_array_dim)/count_not_empty
    vector_j_diff = np.sum(np.diff(vector_j,axis=diff_along_array_dim),axis=sum_along_array_dim)/vector_j.shape[sum_along_array_dim]
    strain_rate = divide_zero(vector_vi_diff,vector_j_diff)
    strain_rate /= shear_rate_scale
    middle_point_vector_j /= diameter

    time = time_in_a_step_from_start_rotate(step)
    ax.plot(middle_point_vector_j, strain_rate, label="{:.2e}".format(time))
    
    return (fig, ax)

########## plot 1D-1D strain_rate-position index change##########
def plotchunk_strain_rate_ij_fix_k_ave(step, fig, ax, n_ave, lines, i, j, k, k_index, figformat="png", ifpickle=False):

    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    
    vector_mv = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines)
    vector_vi = divide_zero(vector_mv, mass)
    vector_j = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[j]], n_ave, lines)

    vector_vi = np.resize(vector_vi, (n_1, n_2))
    vector_j = np.resize(vector_j, (n_1, n_2))

    fix_along_array_dim = position_index_to_array_dim_index[k]

    diff_along_array_dim = position_index_to_array_dim_index[j]
    if diff_along_array_dim == 0:
        middle_point_vector_j = (vector_j[:-1,0] + vector_j[1:,0])/2
    elif diff_along_array_dim == 1:
        middle_point_vector_j = (vector_j[0,:-1] + vector_j[0,1:])/2
    else:
        sys.exit("error")

    vector_vi_diff = np.take(np.diff(vector_vi,axis=diff_along_array_dim), k_index, axis=fix_along_array_dim)
    vector_j_diff = np.take(np.diff(vector_j,axis=diff_along_array_dim), k_index, axis=fix_along_array_dim)
    strain_rate = divide_zero(vector_vi_diff,vector_j_diff)
    strain_rate /= shear_rate_scale
    middle_point_vector_j /= diameter

    time = time_in_a_step_from_start_rotate(step)
    ax.plot(middle_point_vector_j, strain_rate, label="{:.2e}".format(time))
    
    return (fig, ax)
    
########## plot 1D-1D strain_rate-position index change##########
def plotchunk_strain_rate_ij_ave_k_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if ifmanysimu:
        stepsarray.sort()
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(stepsarray, step2, side='right')
            stepsarray_index = stepsarray[:sort_final_index]

            for step in stepsarray_index:
                fig, ax = plotchunk_strain_rate_ij_ave_k_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            stepsarray = stepsarray[sort_final_index:]
    else:
        for step in stepsarray:
            fig, ax = plotchunk_strain_rate_ij_ave_k_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
        
    ax.legend(
              title = "Time (s)",
              bbox_to_anchor=(1.04,1), 
              loc="upper left"
              )
    ax.set_xlabel(map_dim_index_to_coordinate[j])
    ax.set_ylabel("strain_rate" + map_dim_index_to_coordinate[i] + map_dim_index_to_coordinate[j] + "_average_" + map_dim_index_to_coordinate[k])
    plt.tight_layout()

    return (fig, ax)
########## save plot 1D-1D strain_rate-position index change##########
def save_plotchunk_strain_rate_ij_ave_k_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=False):
    
    add_nve_subfolder_in_folder(n_ave, dp.f_path_strain_rate_i_j_ave_k(i,j,k))
    
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_path_strain_rate_i_j_ave_k(i,j,k))
    
    if inonefigure:
        
        fig, ax = plotchunk_strain_rate_ij_ave_k_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
        
        figurepath_no_extension = foldersave + "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
    else:
        
        for step in stepsarray:

            fig, ax = plotchunk_strain_rate_ij_ave_k_ave_manytime(np.array([step]), n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            
            figurepath_no_extension = foldersave + str(int(step))
            
            fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig, f)

    plt.close('all')


########## plot 1D-1D velocity-position index change##########
def plotchunk_velocity_i_ave_j_xk_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False):

    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    mvi = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines)
    vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], n_ave, lines)
    mass = np.resize(mass, (n_1, n_2))
    mvi = np.resize(mvi, (n_1, n_2))
    vector = np.resize(vector, (n_1, n_2))
    
    sum_along_array_dim = position_index_to_array_dim_index[j]
    
    mvi = np.sum(mvi,axis=sum_along_array_dim)
    mass = np.sum(mass,axis=sum_along_array_dim)
    vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
    velocity = divide_zero(mvi, mass)
    velocity /= velocity_scale

    vector = vector/diameter
    time = time_in_a_step_from_start_rotate(step)
    ax.plot(vector, velocity, label="{:.2e}".format(time))

    return (fig, ax)
########## plot 1D-1D velocity-position index change##########
def plotchunk_velocity_i_ave_j_xk_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if ifmanysimu:
        stepsarray.sort()
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(stepsarray, step2, side='right')
            stepsarray_index = stepsarray[:sort_final_index]
            for step in stepsarray_index:
                fig, ax = plotchunk_velocity_i_ave_j_xk_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            
            stepsarray = stepsarray[sort_final_index:]
    else:
        for step in stepsarray:
            fig, ax = plotchunk_velocity_i_ave_j_xk_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
    
    ax.legend(
        title = "Time (s)",
        bbox_to_anchor=(1.04,1),
        loc="upper left",
        )
    ax.set_xlabel(map_dim_index_to_coordinate[k])
    ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (Vwall)")
    plt.tight_layout()

    return (fig, ax)

########## save plot 1D-1D velocity-position index change##########
def save_plotchunk_velocity_i_ave_j_xk_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

    add_nve_subfolder_in_folder(n_ave, dp.f_path_velocity_i_j_ave_k(i,j,k))
    
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_path_velocity_i_j_ave_k(i,j,k))
    
    if inonefigure:
        
        fig, ax = plotchunk_velocity_i_ave_j_xk_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)

        figurepath_no_extension = foldersave + "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')
    else:
        
        for step in stepsarray:

            fig, ax = plotchunk_velocity_i_ave_j_xk_ave_manytime(np.array([step]), n_ave, lines, i, j, k, figformat="png", ifpickle=False)
                  
            figurepath_no_extension = foldersave + str(int(step))
            
            fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig, f)
            plt.close('all')
 


########## plot 1D-1D ek-position index change##########
def plotchunk_ekovermass_i_ave_j_xk_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False):

    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    eki = value_in_a_step_ave(step, "v_Ek" + map_dim_index_to_coordinate[i], n_ave, lines)
    vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], n_ave, lines)
    mass = np.resize(mass, (n_1, n_2))
    eki = np.resize(eki, (n_1, n_2))
    vector = np.resize(vector, (n_1, n_2))
    
    sum_along_array_dim = position_index_to_array_dim_index[j]

    eki = np.sum(eki,axis=sum_along_array_dim)
    mass = np.sum(mass,axis=sum_along_array_dim)
    vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
    ekovermass = 2*divide_zero(eki, mass)
    
    ekovermass /= velocity_scale**2
    vector = vector/diameter
    time = time_in_a_step_from_start_rotate(step)
    ax.plot(vector, ekovermass, label="{:.2e}".format(time))

    return (fig, ax)
########## plot 1D-1D ek-position index change##########
def plotchunk_ekovermass_i_ave_j_xk_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if ifmanysimu:
        stepsarray.sort()
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(stepsarray, step2, side='right')
            stepsarray_index = stepsarray[:sort_final_index]
            for step in stepsarray_index:
                fig, ax = plotchunk_ekovermass_i_ave_j_xk_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            
            stepsarray = stepsarray[sort_final_index:]
    else:
        for step in stepsarray:
            fig, ax = plotchunk_ekovermass_i_ave_j_xk_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
    
    ax.legend(
        title = "Time (s)",
        bbox_to_anchor=(1.04,1),
        loc="upper left",
        )
    ax.set_xlabel(map_dim_index_to_coordinate[k])
    ax.set_ylabel("<E" + map_dim_index_to_coordinate[i] + ">" + "/" + "<m>")
    plt.tight_layout()

    return (fig, ax)

########## save plot 1D-1D ek-position index change##########
def save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

    add_nve_subfolder_in_folder(n_ave, dp.f_path_ekovermass_i_j_ave_k(i,j,k))
    
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_path_ekovermass_i_j_ave_k(i,j,k))
    
    if inonefigure:
        
        fig, ax = plotchunk_ekovermass_i_ave_j_xk_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)

        figurepath_no_extension = foldersave + "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')
    else:
        
        for step in stepsarray:

            fig, ax = plotchunk_ekovermass_i_ave_j_xk_ave_manytime(np.array([step]), n_ave, lines, i, j, k, figformat="png", ifpickle=False)
                  
            figurepath_no_extension = foldersave + str(int(step))
            
            fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig, f)
            plt.close('all')
 
    


##########  plot 1D-1D ekdiff-position index change##########
def plotchunk_ekminusekaveovermass_i_ave_j_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False):

    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    ekminusekavi = divide_zero(
            2*value_in_a_step_ave(step, "v_Ek" + map_dim_index_to_coordinate[i], n_ave, lines) - divide_zero(value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines)**2, mass),
            mass
            )
            
    vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], n_ave, lines)
    mass = np.resize(mass, (n_1, n_2))
    ekminusekavi = np.resize(ekminusekavi, (n_1, n_2))
    vector = np.resize(vector, (n_1, n_2))
    
    sum_along_array_dim = position_index_to_array_dim_index[j]

    ekminusekavi = np.sum(ekminusekavi,axis=sum_along_array_dim)
    mass = np.sum(mass,axis=sum_along_array_dim)
    vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
    ekminusekaveovermass = divide_zero(ekminusekavi, mass)
    
    ekminusekaveovermass /= velocity_scale**2
    vector = vector/diameter
    time = time_in_a_step_from_start_rotate(step)
    ax.plot(vector, ekminusekaveovermass, label="{:.2e}".format(time))

    return (fig, ax)
##########  plot 1D-1D ekdiff-position index change##########    
def plotchunk_ekminusekaveovermass_i_ave_j_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if ifmanysimu:
        stepsarray.sort()
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(stepsarray, step2, side='right')
            stepsarray_index = stepsarray[:sort_final_index]
            for step in stepsarray_index:
                fig, ax = plotchunk_ekminusekaveovermass_i_ave_j_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            stepsarray = stepsarray[sort_final_index:]
    else:
        for step in stepsarray:
            fig, ax = plotchunk_ekminusekaveovermass_i_ave_j_ave(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
    

    ax.legend(
              title = "Time (s)",
              bbox_to_anchor=(1.04,1),
              loc="upper left",
              )
    ax.set_xlabel(map_dim_index_to_coordinate[k])
    ax.set_ylabel("<V" + map_dim_index_to_coordinate[i] + "^2>" + " - <V" + map_dim_index_to_coordinate[i] + ">^2")
    plt.tight_layout()

    return (fig, ax)

########## save plot 1D-1D ekdiff-position index change##########
def save_plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

    add_nve_subfolder_in_folder(n_ave, dp.f_path_ekminusekaveovermass_i_j_ave_k(i,j,k))
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_path_ekminusekaveovermass_i_j_ave_k(i,j,k))
    
    if inonefigure:
        
        fig, ax = plotchunk_ekminusekaveovermass_i_ave_j_ave_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
        
        figurepath_no_extension = foldersave + "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

    else:
        
        for step in stepsarray:

            fig, ax = plotchunk_ekminusekaveovermass_i_ave_j_ave_manytime(np.array([step]), n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            
            figurepath_no_extension = foldersave + str(int(step))
            
            fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig, f)
            plt.close('all')
    plt.close('all')

def plotchunk_1(stepsarray, n_ave, lines, figformat="png", ifpickle=False):
    save_v23_x23(stepsarray, n_ave, lines, figformat="png", ifpickle=False)
    
    save_v13_x23(stepsarray, n_ave, lines, figformat="png", ifpickle=False)
    
    save_fraction_x23(stepsarray, n_ave, lines, figformat="png", ifpickle=False)
    
    for i in [0,1,2]:
        j=1
        k=2
        
        #save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
        
        j=2
        k=1
        
        #save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
        #save_plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
        #save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
        #save_plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
        save_plotchunk_velocity_i_ave_j_xk_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
    
    i=0
    j=1
    k=2
    save_plotchunk_strain_rate_ij_ave_k_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
    save_plotchunk_strain_rate_ij_ave_k_ave(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=True)
    save_strain_rate_ij_x23(stepsarray, n_ave, lines, i, j, figformat="png", ifpickle=False)


def plotchunk_1_step1to2(func, step1, step2, n_ave, lines, figformat="png", ifpickle=False):
    #if smallstep largestep
    # need combine
    stepsarray = np.arange(step1, step2, d_step)
    return func(stepsarray, n_ave, lines, figformat="png", ifpickle=False)

def plotchunk_1_every(func, n_ave, lines, figformat="png", ifpickle=False):
    # need combine
    step1_default = step_first_in_file_change_by_n_ave(n_ave, lines)
    step2_default = step_last_in_file_change_by_n_ave(n_ave, lines)
    return plotchunk_1_step1to2(func, step1_default, step2_default, n_ave, lines, figformat="png", ifpickle=False)

def plotchunk_1_list(func, stepslist, n_ave, lines, figformat="png", ifpickle=False):
    #if smallstep largestep
    # need combine
    stepsarray = np.array(stepslist)
    return func(stepsarray, n_ave, lines, figformat="png", ifpickle=False)

def plotchunk_velocity_fraction_step1to2(step1, step2, n_ave, lines, figformat="png", ifpickle=False):
    #if smallstep largestep
    # need combine
    plotchunk_1_step1to2(save_v23_x23, step1, step2, n_ave, lines, figformat="png", ifpickle=False)
    plotchunk_1_step1to2(save_v13_x23, step1, step2, n_ave, lines, figformat="png", ifpickle=False)
    plotchunk_1_step1to2(save_fraction_x23, step1, step2, n_ave, lines, figformat="png", ifpickle=False)

def plotchunk_velocity_fraction_every(n_ave, lines, figformat="png", ifpickle=False):
    #if smallstep largestep
    # need combine
    plotchunk_1_every(save_v23_x23, n_ave, lines, figformat="png", ifpickle=False)
    plotchunk_1_every(save_v13_x23, n_ave, lines, figformat="png", ifpickle=False)
    plotchunk_1_every(save_fraction_x23, n_ave, lines, figformat="png", ifpickle=False)


########## get 1D time-1D fraction - 2D x23 ##########

def manysteparray(steparray, v_name, n_ave, lines):

    n_step = len(steparray)

    initial_array = np.empty([n_line_in_a_step, n_step])
    
    for i, step in enumerate(steparray):
        
        initial_array[:,i] = value_in_a_step_ave(step, v_name, n_ave, lines)

    return initial_array


########## plot 1D volumn fraction-time for every grid ##########
def check_volumnfraction_x23_increase_decrease_plot(steparray, n_ave, lines, figformat="png", ifpickle=False, ifmanysimu=True):
    
    om.create_directory(dp.f_fraction_check_everygrid)
    tolerence = 0.5

    if ifmanysimu:
        steparray.sort()
        mass = np.empty([n_line_in_a_step, 0])
        time = np.empty(0)
        steparray_new = np.empty(0)
        
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(steparray, step2, side='right')
            steparray_index = steparray[:sort_final_index]
            mass = np.append(mass, manysteparray(steparray_index, "c_m1", n_ave, lines), axis=-1)
            time = np.append(time, time_in_a_step(steparray_index), axis=-1)
            steparray_new = np.append(steparray_new, steparray_index, axis=-1)
            steparray = steparray[sort_final_index:]

        steparray = steparray_new
    else:
        mass = manysteparray(steparray, "c_m1", n_ave, lines)
        time = time_in_a_step(steparray)

    volumn_fraction = mass/float(rr.logfile['den'])/vol_in_chunks
    diff_volumn_fraction = np.diff(volumn_fraction, axis=1)
    diff_time = time_in_a_step(np.diff(steparray, axis=0))
    diff_t_volumn_fraction = diff_volumn_fraction/diff_time

    if chunk_method == "rz":
        r = value_in_a_step_ave(steparray[0], "v_r", n_ave, lines)
        z = value_in_a_step_ave(steparray[0], "v_z", n_ave, lines)
        x_array = r/diameter
        y_array = z/diameter
    elif chunk_method == "yz":
        y = value_in_a_step_ave(steparray[0], "Coord1", n_ave, lines)
        z = value_in_a_step_ave(steparray[0], "Coord2", n_ave, lines)
        x_array = y/diameter
        y_array = z/diameter
    else:
        sys.exit("chunk_method wrong")

    folder_increase = dp.f_fraction_check_everygrid + "increase/"
    om.create_directory(folder_increase)
    folder_decrease = dp.f_fraction_check_everygrid + "decrease/"
    om.create_directory(folder_decrease)
    folder_increase_decrease = dp.f_fraction_check_everygrid + "increase_decrease/"
    om.create_directory(folder_increase_decrease)

    for i in range(n_line_in_a_step):
        if not np.all(diff_t_volumn_fraction[i,:] < 0.1):
            if np.all(diff_t_volumn_fraction[i,:] > (0 - tolerence)) or np.all(diff_t_volumn_fraction[i,:] < (0 + tolerence)):
                    
                if np.all(diff_t_volumn_fraction[i,:] > (0 - tolerence)) and np.all(diff_t_volumn_fraction[i,:] < (0 + tolerence)):
                    folder = folder_increase_decrease
                elif np.all(diff_t_volumn_fraction[i,:] > (0 - tolerence)):
                    folder = folder_increase
                elif np.all(diff_t_volumn_fraction[i,:] < (0 + tolerence)):
                    folder = folder_decrease
                om.create_directory(folder)
                if chunk_method == "rz":
                    labelstring = "r_ " + "{:.2e}".format(x_array[i]) + ", z_ " + "{:.2e}".format(y_array[i])
                elif chunk_method == "yz":
                    labelstring = "y_ " + "{:.2e}".format(x_array[i]) + ", z_ " + "{:.2e}".format(y_array[i])
                else:
                    sys.exit("chunk_method wrong")

                fig = plt.figure()
                ax = fig.add_subplot(111)

                ax.plot(time, volumn_fraction[i,:], label=labelstring)
                ax.legend(
                    title="position",
                    bbox_to_anchor=(1.04,1),
                    loc="upper left",
                )
                ax.set_xlabel("time (s)")
                ax.set_ylabel("volumn fraction")
                plt.tight_layout()

                fig.savefig(folder + labelstring + "." + figformat, format=figformat)
                if ifpickle:
                    # Save figure handle to disk
                    with open(folder + labelstring + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                        pickle.dump(fig, f)

                plt.close('all')


########## plot 1D-1D velocity i ave j fix k - time index change##########
def plotchunk_velocity_i_time_ave_j_fix_k_ave(steparray, n_ave, lines, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True):
    
    steparray.astype(int)

    if ifmanysimu:
        steparray.sort()
        mass = np.empty([n_line_in_a_step, 0])
        vector_mv = np.empty([n_line_in_a_step, 0])
        steparray_new = np.empty(0)
        
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(steparray, step2, side='right')
            steparray_index = steparray[:sort_final_index]

            mass = np.append(mass, manysteparray(steparray_index, "c_m1", n_ave, lines), axis=-1)
            vector_mv = np.append(vector_mv, manysteparray(steparray_index, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines), axis=-1)
            steparray_new = np.append(steparray_new, steparray_index, axis=-1)
            steparray = steparray[sort_final_index:]

        steparray = steparray_new
    else:
        mass = manysteparray(steparray, "c_m1", n_ave, lines)
        vector_mv = manysteparray(steparray_index, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines)

    vector_vi_sum_array = np.empty(0)
    time_array = np.empty(0)

    for nn, step in enumerate(steparray):
        
        vector_mv_step = np.resize(vector_mv[:,nn], (n_1, n_2))
        mass_step = np.resize(mass[:,nn], (n_1, n_2))
        vector_vi = divide_zero(vector_mv_step, mass_step)
        fix_along_array_dim = position_index_to_array_dim_index[k]
        ave_along_array_dim = position_index_to_array_dim_index[j]
                
        count_not_empty = np.sum(if_grid_surround_not_empty(mass_step), axis=ave_along_array_dim)
    
        vector_vi_sum = np.take(np.sum(vector_vi,axis=ave_along_array_dim), k_index, axis=fix_along_array_dim)/np.take(count_not_empty, k_index, axis=fix_along_array_dim)

        time = time_in_a_step_from_start_rotate(step)
        vector_vi_sum_array= np.append(vector_vi_sum_array, vector_vi_sum)
        time_array= np.append(time_array, time)
    
    vector_vi_sum_array = vector_vi_sum_array/velocity_scale

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(time_array, vector_vi_sum_array, 
            label=labelstring_size_walltype,
            )

    ax.legend(
        title="time (s)",
        bbox_to_anchor=(1.04,1),
        loc="upper left",
        )
    ax.set_xlabel("time (s)")
    ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (Vwall)" + " at " + map_dim_index_to_coordinate[k] + "=" + str(k_index))
    plt.tight_layout()

    return (fig, ax)


def save_plotchunk_velocity_i_time_ave_j_fix_k_ave(steparray, n_ave, lines, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True):

    om.create_directory(dp.f_velocity_i_time_ave_j_fix_k_ave_path)

    folder = dp.f_velocity_i_time_ave_j_fix_k_ave_path + "V_" + str(i) + "_ave_" + str(j) + "_fix_" + str(k) + "/"
    om.create_directory(folder)

    add_nve_subfolder_in_folder(n_ave, folder)
    
    foldersave = path_nve_subfolder_in_folder(n_ave, folder)
        
    fig, ax = plotchunk_velocity_i_time_ave_j_fix_k_ave(steparray, n_ave, lines, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True)

    figurepath_no_extension = foldersave + "combine" + str(int(steparray[0])) + "_to_" + str(int(steparray[-1]))

    fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
    if ifpickle:
        # Save figure handle to disk
        with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
            pickle.dump(fig, f)
    plt.close('all')
    


def save_plotchunk_velocity_i_time_near_wall_ave(steparray, n_ave, lines, i, figformat="png", ifpickle=False, ifmanysimu=True):
    j=2
    k=1
    k_index=0

    om.create_directory(dp.f_velocity_i_time_ave_j_fix_k_ave_path)

    folder = dp.f_velocity_i_time_ave_j_fix_k_ave_path + "V_" + str(i) + "_near_wall" + "/"
    om.create_directory(folder)

    add_nve_subfolder_in_folder(n_ave, folder)
    
    foldersave = path_nve_subfolder_in_folder(n_ave, folder)
        
    fig, ax = plotchunk_velocity_i_time_ave_j_fix_k_ave(steparray, n_ave, lines, i, j, k, k_index, figformat="png", ifpickle=False, ifmanysimu=True)

    figurepath_no_extension = foldersave + "combine" + str(int(steparray[0])) + "_to_" + str(int(steparray[-1]))

    fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
    if ifpickle:
        # Save figure handle to disk
        with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
            pickle.dump(fig, f)
    plt.close('all')
    

########## plot 1D-1D fraction-ave over j position k index change##########
def plotchunk_fraction_ave_j_ave(step, fig, ax, n_ave, lines, j, k, figformat="png", ifpickle=False, gridwidth=2):
    
    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], n_ave, lines)
    mass = np.resize(mass, (n_1, n_2))
    vector = np.resize(vector, (n_1, n_2))

    new_mass = 0
    new_vector = 0

    for nnn in range(gridwidth):
        new_mass += mass[nnn::gridwidth, :]
        new_vector += vector[nnn::gridwidth, :]

    mass = new_mass
    vector = new_vector
    fraction = mass/float(rr.logfile['den'])/(vol_in_chunks*gridwidth)

    not_surroundemptyorselfempty = if_grid_surround_not_empty(mass)

    sum_along_array_dim = position_index_to_array_dim_index[j]

    count_not_surroundemptyorselfempty = np.sum(not_surroundemptyorselfempty, axis=sum_along_array_dim)
    fraction_only_not_surround_empty = fraction*not_surroundemptyorselfempty

    fraction_ave = np.sum(fraction_only_not_surround_empty, axis=sum_along_array_dim)/count_not_surroundemptyorselfempty

    vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
    
    vector = vector/diameter

    time = time_in_a_step_from_start_rotate(step)

    ax.plot(vector, fraction_ave, label="{:.2e}".format(time))

    return (fig, ax)

########## plot 1D-1D fraction-ave over j position k index change##########
def plotchunk_fraction_ave_j_ave_manytime(stepsarray, n_ave, lines, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if ifmanysimu:
        stepsarray.sort()
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(stepsarray, step2, side='right')
            stepsarray_index = stepsarray[:sort_final_index]
            for step in stepsarray_index:
                fig, ax = plotchunk_fraction_ave_j_ave(step, fig, ax, n_ave, lines, j, k, figformat="png", ifpickle=False)
            
            stepsarray = stepsarray[sort_final_index:]
    else:
        for step in stepsarray:
            fig, ax = plotchunk_fraction_ave_j_ave(step, fig, ax, n_ave, lines, j, k, figformat="png", ifpickle=False)
    
    ax.legend(
        title = "time (s)",
        bbox_to_anchor=(1.04,1),
        loc="upper left",
        )
    ax.set_xlabel(map_dim_index_to_coordinate[k])
    ax.set_ylabel("Volumn Fraction")
    plt.tight_layout()

    return (fig, ax)

########## save plot fraction-ave over j position k index change##########
def save_plotchunk_fraction_ave_j_ave(stepsarray, n_ave, lines, j, k, figformat="png", ifpickle=False, inonefigure=False):
    
    f_path_fraction = dp.diagram_path + "fraction/"
    om.create_directory(f_path_fraction)

    f_path_ave_jk = f_path_fraction + "ave_" + str(j) + "_x" + str(k) + "/"

    add_nve_subfolder_in_folder(n_ave, f_path_ave_jk)
    
    foldersave = path_nve_subfolder_in_folder(n_ave, f_path_ave_jk)
    
    if inonefigure:
        
        fig, ax = plotchunk_fraction_ave_j_ave_manytime(stepsarray, n_ave, lines, j, k, figformat="png", ifpickle=False)

        figurepath_no_extension = foldersave + "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')
    else:
        
        for step in stepsarray:

            fig, ax = plotchunk_fraction_ave_j_ave_manytime(np.array([step]), n_ave, lines, j, k, figformat="png", ifpickle=False)
                  
            figurepath_no_extension = foldersave + str(int(step))
            
            fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig, f)
            plt.close('all')
 
    plt.close('all')


def plotVymax_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False):
    add_nve_subfolder_in_folder(n_ave, dp.f_max_velocity_near_wall)

    with open(dp.lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:  
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        n_line_in_a_step = int(lines[3].split()[1])

    def plotVymax_1(step1_1, step2_1, figformat="png", ifpickle=False):
        data0 = np.array([lines[t].split() for t in range(4, 4+int(n_line_in_a_step))], dtype=np.float64)
        
        if chunk_method == "rz":
            x_array, y_array = np.meshgrid(
                                        int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                        (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                        )
            x_array = x_array.reshape((-1))
            x_array_to_movingwall_dp_unit = x_array - rr.logfile['ri_wall_dp_unit']

        elif chunk_method == "yz":
            if rr.logfile["chunk/atom"][1] == "y":
                x_array = data0[:,dic_index_of_variable_in_header(lines)["Coord1"]]
                x_array = x_array/diameter
                
            elif rr.logfile["chunk/atom"][1] == "z":
                x_array = data0[:,dic_index_of_variable_in_header(lines)["Coord2"]]
                x_array = x_array/diameter
            else:
                sys.exit("wrong")
            x_array_to_movingwall_dp_unit = x_array
        else:
            sys.exit("chunk_method wrong")

        index_select_array = (x_array_to_movingwall_dp_unit<=2.5)

        n_line_in_a_step_new = np.sum(1*index_select_array)
        step_array = np.arange(step1_1, step2_1, d_step)
        data_array = np.empty([step_array.shape[0], n_line_in_a_step_new, len(header)])
        for index in range(step_array.shape[0]):
            data_array[index,:,:] = data_in_one_step(step_array[index],lines)

        step_sum_all = 0
        for index in range(n_ave):
            step_sum_all += step_array[index:step_array.shape[0]+1-n_ave+index]
        step_mean_all = step_sum_all/n_ave
        data_sum_all = 0
        for index in range(n_ave):
            data_sum_all += data_array[index:step_array.shape[0]+1-n_ave+index]
        data_mean_all = data_sum_all/n_ave

        max_vy = np.empty(step_mean_all.shape[0])
        ave_vy = np.empty(step_mean_all.shape[0])
        for index in range(step_mean_all.shape[0]):
 
            data = data_mean_all[index]

            if chunk_method == "rz":
                x_array, y_array = np.meshgrid(
                                            int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                            (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                            )
                x_array = x_array.reshape((-1))
                y_array = y_array.reshape((-1))
                vy_array = divide_zero(data[:,dic_index_of_variable_in_header(lines)["v_mvz"]],data[:,dic_index_of_variable_in_header(lines)["c_m1"]])/velocity_scale
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = data[:,dic_index_of_variable_in_header(lines)["Coord1"]]
                    y_array = data[:,dic_index_of_variable_in_header(lines)["Coord2"]]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = data[:,dic_index_of_variable_in_header(lines)["Coord2"]]
                    y_array = data[:,dic_index_of_variable_in_header(lines)["Coord1"]]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
                vy_array = divide_zero(data[:,dic_index_of_variable_in_header(lines)["v_mvz"]],data[:,dic_index_of_variable_in_header(lines)["c_m1"]])/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            kk = np.argmax(np.absolute(vy_array))
            max_vy[index] = vy_array[kk]
            ave_vy[index] = np.mean(vy_array)
        
        time_mean_all = time_in_a_step_from_start_rotate(step_mean_all)
        
        fig_handle = plt.figure()

        plt.xlabel('time')
        plt.ylabel('maxVz')
        plt.plot(time_mean_all, max_vy)
        plt.tight_layout()
        
        fig_handle.savefig(dp.f_max_velocity_near_wall + "maxVy" + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(dp.f_max_velocity_near_wall + "maxVy" + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig_handle, f)
        plt.close('all')
    
    if if_plot_to_last:
        step1_default = step_first_in_file_change_by_n_ave(n_ave,lines)
        step2_default = step_last_in_file_change_by_n_ave(n_ave,lines)
        plotVymax_1(step1_default, step2_default, figformat="png", ifpickle=ifpickle)
    else:
        plotVymax_1(step1, step2, figformat="png", ifpickle=ifpickle)



########## plot 1D-1D velocity(no top)-position index change##########
def plotchunk_velocity_i_ave_j_xk_ave_no_top(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False):

    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    mvi = value_in_a_step_ave(step, "v_mv" + map_dim_index_to_coordinate[i], n_ave, lines)
    vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], n_ave, lines)
    mass = np.resize(mass, (n_1, n_2))
    mvi = np.resize(mvi, (n_1, n_2))
    vector = np.resize(vector, (n_1, n_2))

    mvi    = mvi*if_grid_surround_not_empty(mass)
    mass   = mass*if_grid_surround_not_empty(mass)

    sum_along_array_dim = position_index_to_array_dim_index[j]
    mvi = np.sum(mvi,axis=sum_along_array_dim)
    mass = np.sum(mass,axis=sum_along_array_dim)
    vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
    velocity = divide_zero(mvi, mass)
    velocity /= velocity_scale
    vector = vector/diameter
    time = time_in_a_step_from_start_rotate(step)
    ax.plot(vector, velocity, label="{:.2e}".format(time))

    return (fig, ax)
########## plot 1D-1D velocity(no top)-position index change##########
def plotchunk_velocity_i_ave_j_xk_ave_no_top_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if ifmanysimu:
        stepsarray.sort()
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(stepsarray, step2, side='right')
            stepsarray_index = stepsarray[:sort_final_index]
            for step in stepsarray_index:
                fig, ax = plotchunk_velocity_i_ave_j_xk_ave_no_top(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            
            stepsarray = stepsarray[sort_final_index:]
    else:
        for step in stepsarray:
            fig, ax = plotchunk_velocity_i_ave_j_xk_ave_no_top(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
    
    ax.legend(
        title = "Time (s)",
        bbox_to_anchor=(1.04,1),
        loc="upper left",
        )
    ax.set_xlabel(map_dim_index_to_coordinate[k] + " (d)")
    ax.set_ylabel("V" + map_dim_index_to_coordinate[i] + " (V_wall)")
    plt.tight_layout()

    return (fig, ax)

########## save plot 1D-1D velocity(no top)-position index change##########
def save_plotchunk_velocity_i_ave_j_xk_ave_no_top(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

    add_nve_subfolder_in_folder(n_ave, dp.f_path_velocity_i_j_ave_k_no_top(i,j,k))
    
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_path_velocity_i_j_ave_k_no_top(i,j,k))
    
    if inonefigure:
        
        fig, ax = plotchunk_velocity_i_ave_j_xk_ave_no_top_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)

        figurepath_no_extension = foldersave + "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

    else:
        
        for step in stepsarray:

            fig, ax = plotchunk_velocity_i_ave_j_xk_ave_no_top_manytime(np.array([step]), n_ave, lines, i, j, k, figformat="png", ifpickle=False)
                  
            figurepath_no_extension = foldersave + str(int(step))
            
            fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig, f)
            plt.close('all')
 




########## plot 1D-1D ek(no top)-position index change##########
def plotchunk_ek_i_ave_j_xk_ave_no_top(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False, ifbottom=True):

    mass = value_in_a_step_ave(step, "c_m1", n_ave, lines)
    eki = value_in_a_step_ave(step, "v_Ek" + map_dim_index_to_coordinate[i], n_ave, lines)
    vector = value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], n_ave, lines)
    mass = np.resize(mass, (n_1, n_2))
    eki = np.resize(eki, (n_1, n_2))
    vector = np.resize(vector, (n_1, n_2))
    if not ifbottom:
        eki    = eki*if_grid_surround_not_empty(mass)
        mass   = mass*if_grid_surround_not_empty(mass)
    if ifbottom:
        mass = mass[:, 0:5]
        eki = eki[:, 0:5]
        vector = vector[:, 0:5]


    sum_along_array_dim = position_index_to_array_dim_index[j]
    eki = np.sum(eki,axis=sum_along_array_dim)
    mass = np.sum(mass,axis=sum_along_array_dim)
    vector = np.sum(vector,axis=sum_along_array_dim)/vector.shape[sum_along_array_dim]
    ekiovermass = divide_zero(eki, mass)
    
    ekiovermass /= velocity_scale**2
    vector = vector/diameter
    time = time_in_a_step_from_start_rotate(step)
    ax.plot(vector, ekiovermass, label="{:.2e}".format(time))

    return (fig, ax)
########## plot 1D-1D ek(no top)-position index change##########
def plotchunk_ek_i_ave_j_xk_ave_no_top_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, ifmanysimu=True):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if ifmanysimu:
        stepsarray.sort()
        for index in range(rr.n_loglist):
            lines = lines_from_one_simu(rc.folder_path_list_initial_to_last[index])
            step2 = step_last_fix_change_by_n_ave(n_ave, index)
            sort_final_index = np.searchsorted(stepsarray, step2, side='right')
            stepsarray_index = stepsarray[:sort_final_index]
            for step in stepsarray_index:
                fig, ax = plotchunk_ek_i_ave_j_xk_ave_no_top(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
            
            stepsarray = stepsarray[sort_final_index:]
    else:
        for step in stepsarray:
            fig, ax = plotchunk_ek_i_ave_j_xk_ave_no_top(step, fig, ax, n_ave, lines, i, j, k, figformat="png", ifpickle=False)
    
    ax.legend(
        title = "Time (s)",
        bbox_to_anchor=(1.04,1),
        loc="upper left",
        )
    ax.set_xlabel(map_dim_index_to_coordinate[k] + " (d)")
    ax.set_ylabel("ek_" + map_dim_index_to_coordinate[i] + "/mass_" + map_dim_index_to_coordinate[j] + " (V_wall^2)")
    plt.tight_layout()

    return (fig, ax)

########## save plot 1D-1D ek(no top)-position index change##########
def save_plotchunk_ek_i_ave_j_xk_ave_no_top(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False, inonefigure=False):

    add_nve_subfolder_in_folder(n_ave, dp.f_path_ek_i_j_ave_k_no_top(i,j,k))
    
    foldersave = path_nve_subfolder_in_folder(n_ave, dp.f_path_ek_i_j_ave_k_no_top(i,j,k))
    
    if inonefigure:
        
        fig, ax = plotchunk_ek_i_ave_j_xk_ave_no_top_manytime(stepsarray, n_ave, lines, i, j, k, figformat="png", ifpickle=False)

        figurepath_no_extension = foldersave + "combine" + str(int(stepsarray[0])) + "_to_" + str(int(stepsarray[-1]))

        fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

    else:
        
        for step in stepsarray:

            fig, ax = plotchunk_ek_i_ave_j_xk_ave_no_top_manytime(np.array([step]), n_ave, lines, i, j, k, figformat="png", ifpickle=False)
                  
            figurepath_no_extension = foldersave + str(int(step))
            
            fig.savefig(figurepath_no_extension + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(figurepath_no_extension + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig, f)
            plt.close('all')
