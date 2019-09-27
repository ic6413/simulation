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
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D
# import module
import datapath as dp
import read_setting.read_setting as rr
# define function for extract data from fix txt to dataframe
if rr.logfile["shearwall"] == "zcylinder":
    chunk_method = 'rz'
    
if rr.logfile["shearwall"] == "yplane":
    chunk_method = 'yz'
    wallfile1 = "force_zbottom_to_particle.allstep"
    wallfile2 = "force_y_top_to_particle.allstep"
    wallfile3 = "force_y_bottom_to_particle.allstep"

diameter = float(rr.logfile['dp'])
width_wall_dp_unit = int(rr.logfile['width_wall_dp_unit'])
if chunk_method == "rz":
    ri = diameter*int(rr.logfile['ri_wall_dp_unit']) 
elif chunk_method == "yz":
    x_period = diameter*int(rr.logfile['x_period_dp_unit'])
else:
    sys.exit("chunk_method wrong")
g = float(rr.logfile['g'])
d_step = int(rr.logfile['freq_ave_wall'])
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

if chunk_method == "rz":
    n_1 = int(rr.logfile['N_bin_r'])
    n_2 = int(rr.logfile['N_bin_z'])
elif chunk_method == "yz":
    n_1 = int(rr.logfile['N_bin_y'])
    n_2 = int(rr.logfile['N_bin_z'])
else:
    sys.exit("chunk_method wrong")

n_12 = n_1*n_2
if chunk_method == "rz":
    x_array, y_array = np.meshgrid(
                                int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                )
elif chunk_method == "yz":
    if rr.logfile["chunk/atom"][1] == "y":
        y_array, x_array = np.meshgrid(
                                       (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                       (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                    )
    elif rr.logfile["chunk/atom"][1] == "z":
        x_array, y_array = np.meshgrid(
                                       (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                       (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                    )
    else:
        sys.exit("wrong")
    
else:
    sys.exit("chunk_method wrong")

dx = 1/n_1*width_wall_dp_unit
dy = 1/n_2*height_dpunit
x_array = x_array.reshape((-1))
y_array = y_array.reshape((-1))
if chunk_method == "rz":
    vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*diameter**3
elif chunk_method == "yz":
    vol_in_chunks = x_period*dx*dy*diameter**2
else:
    sys.exit("chunk_method wrong")


def plot_wall_force(if_plot_to_last, step1, step2):
    for wallfile in [wallfile1, wallfile2, wallfile3]:    
        with open(dp.lammps_directory + "output/wall/" + wallfile) as f:
            
            lines = f.read().strip().split('\n')
            
            header = lines[1].split()[1:]
            step1_default = int(lines[2].split()[0])
            step2_default = int(lines[-1].split()[0])
            
        def plot_wall_force_1(step1_1, step2_1):
        
            n_line_0 = (step1_1 - step1_default)/d_step + 2
            n_line_1 = (step2_1 - step1_default)/d_step + 2
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            
            ## attach data
            df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            ## repeat timestep

            variable1 = 'v_t'
            
            for variable2 in header:

                x_array = df[variable1].values
                y_array = df[variable2].values

                fig1, ax1 = plt.subplots()

                plt.xlabel(variable1)
                plt.ylabel(variable2)
                
                plt.plot(x_array, y_array)
                
                fig1.savefig(dp.f_wall_force_plot_path + variable2 + "_" + str(step1_1) + "_" + str(step2_1))
                plt.close('all')
    
        if if_plot_to_last:
            plot_wall_force_1(step1_default, step2_default)
        else:
            plot_wall_force_1(step1, step2)


def plot_wall_force_ave(if_plot_to_last, step1, step2, n_ave):
    f_wall_force_plot_path_nve = dp.f_wall_force_plot_path + "nve_" + str(n_ave) + "/"

    if not os.path.isdir(f_wall_force_plot_path_nve): 
            os.mkdir(f_wall_force_plot_path_nve)
    
    for wallfile in [wallfile1, wallfile2, wallfile3]: 
        
        with open(dp.lammps_directory + "output/wall/" + wallfile) as f:
            
            lines = f.read().strip().split('\n')
            
            header = lines[1].split()[1:]
            step1_default = int(lines[2].split()[0])
            step2_default = int(lines[-1].split()[0])
            
        def plot_wall_force_1(step1_1, step2_1):
        
            n_line_0 = (step1_1 - step1_default)/d_step + 2
            n_line_1 = (step2_1 - step1_default)/d_step + 2
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            
            ## attach data
            df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            ## repeat timestep

            variable1 = 'v_t'
            
            for variable2 in header:

                x_array = df[variable1].values
                y_array = df[variable2].values

                def ave_over(array, n):
                    length = len(array)
                    for i in range(n)
                    array[0:length-i]
                x_array = x_array[0:-1]+x_array[0:-1]

                fig1, ax1 = plt.subplots()

                plt.xlabel(variable1)
                plt.ylabel(variable2)
                
                plt.plot(x_array, y_array)
                
                fig1.savefig(dp.f_wall_force_plot_path + variable2 + "_" + str(step1_1) + "_" + str(step2_1))
                plt.close('all')
    
        
        if if_plot_to_last:
            plot_wall_force_1(step1_default, step2_default)
        else:
            plot_wall_force_1(step1, step2)

