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
import read_setting.read_setting as rr
# plot style

plt.style.use('classic')
# define function for extract data from fix txt to dataframe
if rr.logfile["shearwall"] == "zcylinder":
    chunk_method = 'rz'
if rr.logfile["shearwall"] == "yplane":
    chunk_method = 'yz'

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




if chunk_method == "rz":
    n_1 = int(rr.logfile['N_bin_r'])
    n_2 = int(rr.logfile['N_bin_z'])
    n_12 = n_1*n_2
    
    dx = 1/n_1*width_wall_dp_unit
    dy = 1/n_2*height_dpunit
    x_array, y_array = np.meshgrid(
                                int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                )
    x_array = x_array.reshape((-1))
    y_array = y_array.reshape((-1))
    vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*diameter**3
elif chunk_method == "yz":
    n_1 = int(rr.logfile['N_bin_y'])
    n_2 = int(rr.logfile['N_bin_z'])
    n_12 = n_1*n_2
    dx = 1/n_1*width_wall_dp_unit
    dy = 1/n_2*height_dpunit
    vol_in_chunks = x_period*dx*dy*diameter**2
else:
    sys.exit("chunk_method wrong")


def plotchunk(if_plot_to_last, step1, step2, figformat="png", ifpickle=False):
    
    with open(dp.lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
        
        lines = f.read().strip().split('\n')
        
        header = lines[2].split()[1:]
        n_line_in_a_step = int(lines[3].split()[1])
        step1_default = int(lines[3].split()[0])
        step2_default = int(lines[-1 - n_line_in_a_step].split()[0])

    def plotchunk_1(step1_1, step2_1, figformat="png", ifpickle=False):
        
        for step in range(step1_1, step2_1, d_step):
            n_line_0 = (step - step1_1)/d_step*(n_line_in_a_step+1) + 4
            n_line_1 = n_line_0 + n_line_in_a_step
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            
            ## attach data
            df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            ## repeat timestep

            def divide_zero(a,b):
                c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
                return c

            if chunk_method == "rz":
                x_array, y_array = np.meshgrid(
                                            int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                            (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                            )
                x_array = x_array.reshape((-1))
                y_array = y_array.reshape((-1))
                vx_array = divide_zero(df['v_mvr'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = df['Coord1'].values
                    y_array = df['Coord2'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = df['Coord2'].values
                    y_array = df['Coord1'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
                vx_array = divide_zero(df['v_mvy'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
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
            if chunk_method == "rz":
                plt.xlabel('r')
                plt.ylabel('z')
            elif chunk_method == "yz":
                plt.xlabel('y')
                plt.ylabel('z')
            #ax1.set_title('velocity field r-z direction (average over theta)')
            Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                        units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                        )

            ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                        "This arrow present {:.2E} wall velocity in 45 degree".format(label_scale),
                        labelpos='E', coordinates='figure', angle=45)

            fig1.savefig(dp.f_momentum_mass_field_v23x23_path + str(step), format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(dp.f_momentum_mass_field_v23x23_path + str(step) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')

            if chunk_method == "rz":
                x_array, y_array = np.meshgrid(
                                            int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                            (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                            )
                x_array = x_array.reshape((-1))
                y_array = y_array.reshape((-1))
                vx_array = divide_zero(df['v_mvt'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = df['Coord1'].values
                    y_array = df['Coord2'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = df['Coord2'].values
                    y_array = df['Coord1'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
                vx_array = divide_zero(df['v_mvx'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            v_length_array = (vx_array**2+vy_array**2)**0.5
            max_v_length = np.amax(v_length_array)
            if dp.if_plot_velocity_field_scale_same == "yes":
                quiver_scale = dp.quiver_scale_velocity_xaxis_shearplaneshear_yaxis_z
            else:
                quiver_scale = max_v_length/2
            label_scale = quiver_scale
            fig1, ax1 = plt.subplots()
            if chunk_method == "rz":
                plt.xlabel('r(position), theta(velocity)')
                plt.ylabel('z')
            elif chunk_method == "yz":
                plt.xlabel('y(position), x(velocity)')
                plt.ylabel('z')
            else:
                sys.exit("chunk_method wrong")
            
            #ax1.set_title('velocity field r-z direction (average over theta)')
            Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                        units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                        )
            ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                        "This arrow present {:.2E} wall velocity in 45 degree".format(label_scale),
                        labelpos='E', coordinates='figure', angle=45)

            fig1.savefig(dp.f_momentum_mass_field_v13x23_path + str(step), format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(dp.f_momentum_mass_field_v13x23_path + str(step) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')


            quiver_scale = 0.3
            label_scale = quiver_scale
            vx_array = 0*df['c_m1'].values
            vy_array = df['c_m1'].values/float(rr.logfile['den'])/vol_in_chunks
            fig1, ax1 = plt.subplots()
            
            #fig1.figsize = [12.8, 9.6]
            if chunk_method == "rz":
                plt.xlabel('r')
                plt.ylabel('z')
            elif chunk_method == "yz":
                plt.xlabel('y')
                plt.ylabel('z')
            #ax1.set_title('velocity field r-z direction (average over theta)')
            Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                        units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                        )

            ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                        "This arrow present {:.2E}".format(label_scale) + 'solid fraction',
                        labelpos='E',
                        coordinates='figure', angle=90)

            fig1.savefig(dp.f_momentum_mass_field_density_x23_path + str(step), format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(dp.f_momentum_mass_field_density_x23_path + str(step) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')
    
    if if_plot_to_last:
        plotchunk_1(step1_default, step2_default, figformat="png", ifpickle=False)
    else:
        plotchunk_1(step1, step2, figformat="png", ifpickle=False)


def plotchunk_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False):
    f_momentum_mass_field_v13x23_path_nve = dp.f_momentum_mass_field_v13x23_path + "nve_" + str(n_ave) + "/"
    f_momentum_mass_field_v23x23_path_nve = dp.f_momentum_mass_field_v23x23_path + "nve_" + str(n_ave) + "/"
    f_momentum_mass_field_density_x23_path_nve = dp.f_momentum_mass_field_density_x23_path + "nve_" + str(n_ave) + "/"
    post_process_folder_paths = [
        f_momentum_mass_field_v13x23_path_nve,
        f_momentum_mass_field_v23x23_path_nve,
        f_momentum_mass_field_density_x23_path_nve,
    ]
    for post_process_folder_path in post_process_folder_paths:
        if not os.path.isdir(post_process_folder_path): 
                os.mkdir(post_process_folder_path)

    with open(dp.lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
        
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        n_line_in_a_step = int(lines[3].split()[1])
        step1_default = int(lines[3].split()[0])
        step2_default = int(lines[-1 - n_line_in_a_step].split()[0])-n_ave*d_step
        
    def plotchunk_1(step1_1, step2_1, figformat="png", ifpickle=False):
        def data_inloop(step_smallloop):
            n_line_0 = (step_smallloop - step1_1)/d_step*(n_line_in_a_step+1) + 4
            n_line_1 = n_line_0 + n_line_in_a_step
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            data = np.array(data, dtype=np.float64)
            return data
        for step_in in range(step1_1, step2_1, d_step):
            step = 0
            data = 0
            for step_smallloop in range(step_in, step_in+n_ave*d_step, d_step):
                step += step_smallloop
                data += data_inloop(step_smallloop)
                
            step = step/n_ave
            data = data/n_ave
            ## attach data
            df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            ## repeat timestep

            def divide_zero(a,b):
                c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
                return c

            if chunk_method == "rz":
                x_array, y_array = np.meshgrid(
                                            int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                            (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                            )
                x_array = x_array.reshape((-1))
                y_array = y_array.reshape((-1))
                vx_array = divide_zero(df['v_mvr'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = df['Coord1'].values
                    y_array = df['Coord2'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = df['Coord2'].values
                    y_array = df['Coord1'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
                vx_array = divide_zero(df['v_mvy'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
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
            if chunk_method == "rz":
                plt.xlabel('r')
                plt.ylabel('z')
            elif chunk_method == "yz":
                plt.xlabel('y')
                plt.ylabel('z')
            #ax1.set_title('velocity field r-z direction (average over theta)')
            Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                        units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                        )

            ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                        "This arrow present {:.2E} wall velocity in 45 degree".format(label_scale),
                        labelpos='E', coordinates='figure', angle=45)
            fig1.savefig(f_momentum_mass_field_v23x23_path_nve + str(int(step)), format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(f_momentum_mass_field_v23x23_path_nve + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')

            if chunk_method == "rz":
                x_array, y_array = np.meshgrid(
                                            int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit,
                                            (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                            )
                x_array = x_array.reshape((-1))
                y_array = y_array.reshape((-1))
                vx_array = divide_zero(df['v_mvt'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = df['Coord1'].values
                    y_array = df['Coord2'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = df['Coord2'].values
                    y_array = df['Coord1'].values
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
                vx_array = divide_zero(df['v_mvx'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            v_length_array = (vx_array**2+vy_array**2)**0.5
            max_v_length = np.amax(v_length_array)
            if dp.if_plot_velocity_field_scale_same == "yes":
                quiver_scale = dp.quiver_scale_velocity_xaxis_shearplaneshear_yaxis_z
            else:
                quiver_scale = max_v_length/2
            label_scale = quiver_scale
            fig1, ax1 = plt.subplots()
            if chunk_method == "rz":
                plt.xlabel('r(position), theta(velocity)')
                plt.ylabel('z')
            elif chunk_method == "yz":
                plt.xlabel('y(position), x(velocity)')
                plt.ylabel('z')
            else:
                sys.exit("chunk_method wrong")
            
            #ax1.set_title('velocity field r-z direction (average over theta)')
            Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                        units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                        )
            ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                        "This arrow present {:.2E} wall velocity in 45 degree".format(label_scale),
                        labelpos='E', coordinates='figure', angle=45)
            fig1.savefig(f_momentum_mass_field_v13x23_path_nve + str(int(step)), format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(f_momentum_mass_field_v13x23_path_nve + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')


            quiver_scale = 0.2
            label_scale = 0.6
            vx_array = 0*df['c_m1'].values
            vy_array = df['c_m1'].values/float(rr.logfile['den'])/vol_in_chunks
            fig1, ax1 = plt.subplots()
            
            #fig1.figsize = [12.8, 9.6]
            if chunk_method == "rz":
                plt.xlabel('r')
                plt.ylabel('z')
            elif chunk_method == "yz":
                plt.xlabel('y')
                plt.ylabel('z')
            #ax1.set_title('velocity field r-z direction (average over theta)')
            Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                        units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                        )

            ax1.quiverkey(Q, 0.2, 0.95,
                        label_scale, label = "{:.2E}".format(label_scale) + "quiver_scale" + str(quiver_scale) + 'solid fraction',
                        labelpos='E',
                        coordinates='figure', angle=90)

            fig1.savefig(f_momentum_mass_field_density_x23_path_nve + str(int(step)), format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(f_momentum_mass_field_density_x23_path_nve + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')
    
    if if_plot_to_last:
        plotchunk_1(step1_default, step2_default, figformat="png", ifpickle=False)
    else:
        plotchunk_1(step1, step2, figformat="png", ifpickle=False)

