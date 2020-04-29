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
import os
import read_setting as rr
# import calculate setting

# import combine
import combine_other_run.combine_previous_run as cc
# plot style

plt.style.use('classic')
plt.rcParams.update({'font.size': 16})


if "if_inwall_wall_gran" in rr.logfile.keys():
    if rr.logfile["if_inwall_wall_gran"] == "yes":
        if "wall_gran_type" in rr.logfile.keys():
            if rr.logfile["wall_gran_type"] == "1":
                ybottomwalltype = "rough (d=1.1)"
            elif rr.logfile["wall_gran_type"] == "2":
                ybottomwalltype = "rough (d=1)"
            elif rr.logfile["wall_gran_type"] == "3":
                ybottomwalltype = "rough (d=0.9)"
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
labelstring_size_walltype = "L: " + periodlength + "\n" + "W: " + width + "\n" + "H: " + height + "\n" + ybottomwalltype

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

force_scale = 0.6*float(rr.logfile['den'])*g*float(height)**2*float(periodlength)/2*diameter**3


def plot_wall_force(if_plot_to_last, step1, step2, figformat="png", ifpickle=False, ifplotfrominitial=False, ifplotfromrotate=True):
    f_wall_force_plot_path = dp.f_wall_force_plot_path
    os.makedirs(f_wall_force_plot_path, exist_ok=True)
    for wallfile in wallfiles:
        if ifplotfrominitial or ifplotfromrotate:
            with open(rr.lammps_directory + "output/wall/" + wallfile) as f:
                lines = f.read().strip().split('\n')
            header = lines[1].split()[1:]
            df_full = cc.combine_previous_wall_data(wallfile, ifplotfrominitial, ifplotfromrotate)
        else:
            with open(rr.lammps_directory + "output/wall/" + wallfile) as f:
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
                fig = plt.figure()
                ax = fig.add_subplot(111)
                x_array = df["v_t"].values - rr.rotate_start_time
                y_array = df[variable2].values
                if "force" in variable2:
                    y_array /= force_scale
                
                ax.set_xlabel("time (s)")
                
                if "force" in variable2:
                    if "y_top" in variable2:
                        wallstring = " on static wall"
                    
                    elif "y_bottom" in variable2:
                        wallstring = " on moving wall"
                    
                    elif "zbottom" in variable2:
                        wallstring = " on ground"

                    force_string = "F" + variable2[-1]

                    label_y = force_string + wallstring  + " (normalized)"
                else:
                    label_y = variable2
                ax.set_ylabel(label_y)
                
                ax.plot(x_array, y_array, label=labelstring_size_walltype)
                ax.legend(
                    title="",
                    bbox_to_anchor=(1.04,1),
                    loc="upper left",
                    )
                plt.xticks(rotation=45)
                plt.tight_layout()
                fig.savefig(f_wall_force_plot_path + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + "." + figformat, format=figformat)
                if ifpickle:
                    # Save figure handle to disk
                    with open(f_wall_force_plot_path + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                        pickle.dump(fig, f)

                plt.close('all')
    
        if if_plot_to_last:
            plot_wall_force_1(step1_default, step2_default)
        else:
            plot_wall_force_1(step1, step2)


def plot_wall_force_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False, ifplotfrominitial=False, ifplotfromrotate=True):
    f_wall_force_plot_path = dp.f_wall_force_plot_path
    os.makedirs(f_wall_force_plot_path, exist_ok=True)
    f_wall_force_plot_path_nve = f_wall_force_plot_path + "nve_" + str(n_ave) + "/"

    os.makedirs(f_wall_force_plot_path_nve, exist_ok=True)
    
    for wallfile in wallfiles:
        if ifplotfrominitial or ifplotfromrotate:
            with open(rr.lammps_directory + "output/wall/" + wallfile) as f:
                lines = f.read().strip().split('\n')
            header = lines[1].split()[1:]
            df_full = cc.combine_previous_wall_data(wallfile, ifplotfrominitial, ifplotfromrotate)
        else:
            with open(rr.lammps_directory + "output/wall/" + wallfile) as f:
                lines = f.read().strip().split('\n')
            header = lines[1].split()[1:]
            data = [lines[t].split() for t in range(2, len(lines))]
            df_full = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        
        def plot_wall_force_1(step1_1, step2_1):
        
            df = df_full[(df_full['TimeStep']<=step2_1) & (df_full['TimeStep']>=step1_1)]
            ## repeat timestep
            
            for variable2 in header:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                x_array = df["v_t"].values - rr.rotate_start_time
                y_array = df[variable2].values
                if "force" in variable2:
                    y_array /= force_scale

                def ave_over(array, n):
                    length = len(array)
                    anwser_array = 0
                    for i in range(n):
                        anwser_array = anwser_array + array[i: i+length-n_ave+1]
                    anwser_array = anwser_array/n
                    return anwser_array

                x_array = ave_over(x_array, n_ave)
                y_array = ave_over(y_array, n_ave)

                ax.set_xlabel("time (s)")
                if "force" in variable2:
                    if "y_top" in variable2:
                        wallstring = " on static wall"
                    
                    elif "y_bottom" in variable2:
                        wallstring = " on moving wall"
                    
                    elif "zbottom" in variable2:
                        wallstring = " on ground"

                    force_string = "F" + variable2[-1]

                    label_y = force_string + wallstring  + " (normalized)"
                else:
                    label_y = variable2

                ax.set_ylabel(label_y)

                ax.plot(x_array, y_array, label=labelstring_size_walltype)
                ax.legend(
                    title="",
                    bbox_to_anchor=(1.04,1),
                    loc="upper left",
                    )
                plt.xticks(rotation=45)
                plt.tight_layout()


                fig.savefig(f_wall_force_plot_path_nve + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + "." + figformat, format=figformat)
                if ifpickle:
                    # Save figure handle to disk
                    with open(f_wall_force_plot_path_nve + variable2 + "_" + str(step1_1) + "_" + str(step2_1) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                        pickle.dump(fig, f)

                plt.close('all')


        step1_default = int(df_full['TimeStep'].min())
        step2_default = int(df_full['TimeStep'].max())
        
        if if_plot_to_last:
            plot_wall_force_1(step1_default, step2_default)
        else:
            plot_wall_force_1(step1, step2)

