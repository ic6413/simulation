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
import pickle
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D
# import module
import datapath as dp
import osmanage as om
import read_setting.read_setting as rr
# import calculate setting
import read_setting.calculate_setting as rc
# import combine
import combine_other_run.combine_previous_run as cc
# plot style

plt.style.use('classic')

# define function for extract data from fix txt to dataframe
if rr.logfile["shearwall"] == "zcylinder":
    wallfiles = [
                "force_zbottom_to_particle.allstep",
                "force_outwall_to_particle.allstep",
                "force_inwall_to_particle.allstep",
    ]
if rr.logfile["shearwall"] == "yplane":
    wallfiles = [
                "force_zbottom_to_particle.allstep",
                "force_y_top_to_particle.allstep",
                "force_y_bottom_to_particle.allstep",
    ]
    
    

diameter = float(rr.logfile['dp'])
width_wall_dp_unit = int(rr.logfile['width_wall_dp_unit'])

g = float(rr.logfile['g'])
d_step = int(rr.logfile['freq_ave_wall'])


def plot_wall_force(if_plot_to_last, step1, step2, figformat="png", ifpickle=False, ifplotfrominitial=False, ifplotfromrotate=True):
    f_wall_force_plot_path = dp.f_wall_force_plot_path
    om.create_directory(f_wall_force_plot_path)
    for wallfile in wallfiles:
        if ifplotfrominitial or ifplotfromrotate:
            with open(dp.lammps_directory + "output/wall/" + wallfile) as f:
                lines = f.read().strip().split('\n')
            header = lines[1].split()[1:]
            df_full = cc.combine_previous_wall_data(wallfile, ifplotfrominitial, ifplotfromrotate)
        else:
            with open(dp.lammps_directory + "output/wall/" + wallfile) as f:
                lines = f.read().strip().split('\n')
            header = lines[1].split()[1:]
            data = [lines[t].split() for t in range(2, len(lines))]
            df_full = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        
        step1_default = int(df_full['TimeStep'].min())
        step2_default = int(df_full['TimeStep'].max())
        
        def plot_wall_force_1(step1_1, step2_1):
            df = df_full[(df_full['TimeStep']<=step2_1) & (df_full['TimeStep']>=step1_1)]
            
            ## repeat timestep
            
            for variable2 in header:
                fig_handle = plt.figure()
                x_array = df["v_t"].values - rc.rotate_start_time
                y_array = df[variable2].values
                
                plt.xlabel("time(s)")
                if "force" in variable2:
                    label_y = variable2 + " (N)"
                else:
                    label_y = variable2
                plt.ylabel(label_y)
                
                plt.plot(x_array, y_array)
                plt.tight_layout()
                fig_handle.savefig(f_wall_force_plot_path + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + "." + figformat, format=figformat)
                if ifpickle:
                    # Save figure handle to disk
                    with open(f_wall_force_plot_path + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                        pickle.dump(fig_handle, f)

                plt.close('all')
    
        if if_plot_to_last:
            plot_wall_force_1(step1_default, step2_default)
        else:
            plot_wall_force_1(step1, step2)


def plot_wall_force_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False, ifplotfrominitial=False, ifplotfromrotate=True):
    f_wall_force_plot_path = dp.f_wall_force_plot_path
    om.create_directory(f_wall_force_plot_path)
    f_wall_force_plot_path_nve = f_wall_force_plot_path + "nve_" + str(n_ave) + "/"

    om.create_directory(f_wall_force_plot_path_nve)
    
    for wallfile in wallfiles:
        if ifplotfrominitial or ifplotfromrotate:
            with open(dp.lammps_directory + "output/wall/" + wallfile) as f:
                lines = f.read().strip().split('\n')
            header = lines[1].split()[1:]
            df_full = cc.combine_previous_wall_data(wallfile, ifplotfrominitial, ifplotfromrotate)
        else:
            with open(dp.lammps_directory + "output/wall/" + wallfile) as f:
                lines = f.read().strip().split('\n')
            header = lines[1].split()[1:]
            data = [lines[t].split() for t in range(2, len(lines))]
            df_full = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        
        def plot_wall_force_1(step1_1, step2_1):
        
            df = df_full[(df_full['TimeStep']<=step2_1) & (df_full['TimeStep']>=step1_1)]
            ## repeat timestep
            
            for variable2 in header:
                fig_handle = plt.figure()
                x_array = df["v_t"].values - rc.rotate_start_time
                y_array = df[variable2].values

                def ave_over(array, n):
                    length = len(array)
                    anwser_array = 0
                    for i in range(n):
                        anwser_array = anwser_array + array[i: i+length-n_ave+1]
                    anwser_array = anwser_array/n
                    return anwser_array

                x_array = ave_over(x_array, n_ave)
                y_array = ave_over(y_array, n_ave)

                plt.xlabel("time(s)")
                if "force" in variable2:
                    label_y = variable2 + " (N)"
                else:
                    label_y = variable2
                plt.ylabel(label_y)

                plt.plot(x_array, y_array)
                plt.tight_layout()


                fig_handle.savefig(f_wall_force_plot_path_nve + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + "." + figformat, format=figformat)
                if ifpickle:
                    # Save figure handle to disk
                    with open(f_wall_force_plot_path_nve + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                        pickle.dump(fig_handle, f)

                plt.close('all')


        step1_default = int(df_full['TimeStep'].min())
        step2_default = int(df_full['TimeStep'].max())
        
        if if_plot_to_last:
            plot_wall_force_1(step1_default, step2_default)
        else:
            plot_wall_force_1(step1, step2)

