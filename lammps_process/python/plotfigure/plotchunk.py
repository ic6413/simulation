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

if chunk_method == "rz":
    # chink first dim unchange in the begining
    chunk_first_dim = "z"
    chunk_second_dim = "r"
elif chunk_method == "yz":
    if rr.logfile["chunk/atom"][1] == "y":
        chunk_first_dim = "y"
        chunk_second_dim = "z"
    elif rr.logfile["chunk/atom"][1] == "z":
        chunk_first_dim = "z"
        chunk_second_dim = "y"
    else:
        sys.exit("chunk_method wrong")
else:
    sys.exit("chunk_method wrong")


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
        def index_of_variable_in_header(variable_name_in_header):
            return [n for n in range(0,len(header)) if header[n]==variable_name_in_header][0]
    def plotchunk_1(step1_1, step2_1, figformat="png", ifpickle=False):
        def data_inloop(step_smallloop):
            n_line_0 = (step_smallloop - step1_1)/d_step*(n_line_in_a_step+1) + 4
            n_line_1 = n_line_0 + n_line_in_a_step
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            data = np.array(data, dtype=np.float64)
            return data
        
        
        step_array = np.arange(step1_1, step2_1, d_step)
        data_array = np.empty([step_array.shape[0], n_line_in_a_step, len(header)])
        for index in range(step_array.shape[0]):
            data_array[index,:,:] = data_inloop(step_array[index])

        step_sum_all = 0
        for index in range(n_ave):
            step_sum_all += step_array[index:step_array.shape[0]+1-n_ave+index]
        step_mean_all = step_sum_all/n_ave
        data_sum_all = 0
        for index in range(n_ave):
            data_sum_all += data_array[index:step_array.shape[0]+1-n_ave+index]
        data_mean_all = data_sum_all/n_ave

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
        elif chunk_method == "yz":
            if rr.logfile["chunk/atom"][1] == "y":
                x_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                y_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                x_array = x_array/diameter
                y_array = y_array/diameter
                
            elif rr.logfile["chunk/atom"][1] == "z":
                x_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                y_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                x_array = x_array/diameter
                y_array = y_array/diameter
            else:
                sys.exit("wrong")
        else:
            sys.exit("chunk_method wrong")

        if chunk_method == "rz":
            vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvr")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
            vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
        elif chunk_method == "yz":
            vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvy")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
            vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
        else:
            sys.exit("chunk_method wrong")


        for index in range(step_mean_all.shape[0]):    
            step = step_mean_all[index]
            data = data_mean_all[index]
            vx_array = vx_array_all[index]
            vy_array = vy_array_all[index]
            
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
            fig1.savefig(f_momentum_mass_field_v23x23_path_nve + str(int(step)) + "." + figformat, format=figformat)
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
                vx_array = divide_zero(data[:,index_of_variable_in_header("v_mvt")],data[:,index_of_variable_in_header("c_m1")])/velocity_scale
                vy_array = divide_zero(data[:,index_of_variable_in_header("v_mvz")],data[:,index_of_variable_in_header("c_m1")])/velocity_scale
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = data[:,index_of_variable_in_header("Coord1")]
                    y_array = data[:,index_of_variable_in_header("Coord2")]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = data[:,index_of_variable_in_header("Coord2")]
                    y_array = data[:,index_of_variable_in_header("Coord1")]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
                vx_array = divide_zero(data[:,index_of_variable_in_header("v_mvx")],data[:,index_of_variable_in_header("c_m1")])/velocity_scale
                vy_array = divide_zero(data[:,index_of_variable_in_header("v_mvz")],data[:,index_of_variable_in_header("c_m1")])/velocity_scale
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
            fig1.savefig(f_momentum_mass_field_v13x23_path_nve + str(int(step)) + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(f_momentum_mass_field_v13x23_path_nve + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')


            quiver_scale = 0.2
            label_scale = 0.6
            vx_array = 0*data[:,index_of_variable_in_header("c_m1")]
            vy_array = data[:,index_of_variable_in_header("c_m1")]/float(rr.logfile['den'])/vol_in_chunks
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

            fig1.savefig(f_momentum_mass_field_density_x23_path_nve + str(int(step)) + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(f_momentum_mass_field_density_x23_path_nve + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig1, f)
            plt.close('all')
    
    if if_plot_to_last:
        plotchunk_1(step1_default, step2_default, figformat="png", ifpickle=ifpickle)
    else:
        plotchunk_1(step1, step2, figformat="png", ifpickle=ifpickle)


def plotVymax_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False):
    f_max_velocity_near_wall = dp.f_max_velocity_near_wall + "nve_" + str(n_ave) + "/"
    if not os.path.isdir(f_max_velocity_near_wall): 
            os.mkdir(f_max_velocity_near_wall)
    with open(dp.lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
        
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        n_line_in_a_step = int(lines[3].split()[1])
        step1_default = int(lines[3].split()[0])
        step2_default = int(lines[-1 - n_line_in_a_step].split()[0])-n_ave*d_step
        def index_of_variable_in_header(variable_name_in_header):
            return [n for n in range(0,len(header)) if header[n]==variable_name_in_header][0]
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
                x_array = data0[:,index_of_variable_in_header("Coord1")]
                x_array = x_array/diameter
                
            elif rr.logfile["chunk/atom"][1] == "z":
                x_array = data0[:,index_of_variable_in_header("Coord2")]
                x_array = x_array/diameter
            else:
                sys.exit("wrong")
            x_array_to_movingwall_dp_unit = x_array
        else:
            sys.exit("chunk_method wrong")

        index_select_array = (x_array_to_movingwall_dp_unit<=2.5)
        def data_inloop_select(step_smallloop):
            n_line_0 = (step_smallloop - step1_1)/d_step*(n_line_in_a_step+1) + 4
            n_line_1 = n_line_0 + n_line_in_a_step
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            data = np.array(data, dtype=np.float64)
            data = data[index_select_array,:]
            return data
        n_line_in_a_step_new = np.sum(1*index_select_array)
        step_array = np.arange(step1_1, step2_1, d_step)
        data_array = np.empty([step_array.shape[0], n_line_in_a_step_new, len(header)])
        for index in range(step_array.shape[0]):
            data_array[index,:,:] = data_inloop_select(step_array[index])

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
                vy_array = divide_zero(data[:,index_of_variable_in_header("v_mvz")],data[:,index_of_variable_in_header("c_m1")])/velocity_scale
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = data[:,index_of_variable_in_header("Coord1")]
                    y_array = data[:,index_of_variable_in_header("Coord2")]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = data[:,index_of_variable_in_header("Coord2")]
                    y_array = data[:,index_of_variable_in_header("Coord1")]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
                vy_array = divide_zero(data[:,index_of_variable_in_header("v_mvz")],data[:,index_of_variable_in_header("c_m1")])/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            kk = np.argmax(np.absolute(vy_array))
            max_vy[index] = vy_array[kk]
            ave_vy[index] = np.mean(vy_array)
        
        time_mean_all = step_mean_all*float(rr.logfile["ts"])
        
        fig_handle = plt.figure()

        plt.xlabel('time')
        plt.ylabel('maxVz')
        plt.plot(time_mean_all, max_vy)
        plt.tight_layout()
        
        fig_handle.savefig(f_max_velocity_near_wall + "maxVy" + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(f_max_velocity_near_wall + "maxVy" + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig_handle, f)
        plt.close('all')
    
    if if_plot_to_last:
        plotVymax_1(step1_default, step2_default, figformat="png", ifpickle=ifpickle)
    else:
        plotVymax_1(step1, step2, figformat="png", ifpickle=ifpickle)



def plotchunk_shearratexyaveoverz_y_ave(if_plot_to_last, step1, step2, n_ave, picksteplist, figformat="png", ifpickle=False, filepath=None):
    shear_rate_scale = float(rr.logfile['in_velocity'])/(float(rr.logfile['width_wall_dp_unit'])*float(rr.logfile['dp']))
    f_shearratexyaveoverz_y_ave_path_nve = dp.f_shearratexyaveoverz_y_ave_path + "nve_" + str(n_ave) + "/"
    post_process_folder_paths = [
        f_shearratexyaveoverz_y_ave_path_nve
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
        def index_of_variable_in_header(variable_name_in_header):
            return [n for n in range(0,len(header)) if header[n]==variable_name_in_header][0]
    def plotchunk_1_shearratexyaveoverz_y_ave(step1_1, step2_1, figformat="png", ifpickle=False):
        
        def data_inloop(step_smallloop):
            n_line_0 = (step_smallloop - step1_1)/d_step*(n_line_in_a_step+1) + 4
            n_line_1 = n_line_0 + n_line_in_a_step
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            data = np.array(data, dtype=np.float64)
            return data
        
        
        step_array = np.arange(step1_1, step2_1, d_step)
        data_array = np.empty([step_array.shape[0], n_line_in_a_step, len(header)])
        for index in range(step_array.shape[0]):
            data_array[index,:,:] = data_inloop(step_array[index])

        step_sum_all = 0
        for index in range(n_ave):
            step_sum_all += step_array[index:step_array.shape[0]+1-n_ave+index]
        step_mean_all = step_sum_all/n_ave
        data_sum_all = 0
        for index in range(n_ave):
            data_sum_all += data_array[index:step_array.shape[0]+1-n_ave+index]
        data_mean_all = data_sum_all/n_ave

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
        elif chunk_method == "yz":
            if rr.logfile["chunk/atom"][1] == "y":
                x_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                y_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                x_array = x_array/diameter
                y_array = y_array/diameter
                
            elif rr.logfile["chunk/atom"][1] == "z":
                x_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                y_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                x_array = x_array/diameter
                y_array = y_array/diameter
            else:
                sys.exit("wrong")
        else:
            sys.exit("chunk_method wrong")

        if chunk_method == "rz":
            vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvt")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/shear_rate_scale
            vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/shear_rate_scale
        elif chunk_method == "yz":
            vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvx")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/shear_rate_scale
            vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/shear_rate_scale
        else:
            sys.exit("chunk_method wrong")
        n_y = np.sum(x_array == x_array[0])
        
        x_array = x_array.reshape(-1, n_y)
        x_array = x_array[:,0]
        x_array_diff = np.diff(x_array)
        x_array = (x_array[0:-1] + x_array[1:])/2
        if picksteplist=="all":    
            for index in range(step_mean_all.shape[0]):    
                step = step_mean_all[index]
                vx_array = vx_array_all[index]
                vy_array = vy_array_all[index]
                vx_array_diff = np.diff(vx_array.reshape(-1, n_y).sum(axis=1)/n_y)
                vy_array_diff = np.diff(vy_array.reshape(-1, n_y).sum(axis=1)/n_y)

                vx_array = vx_array_diff/(x_array_diff*diameter)
                vy_array = vy_array_diff/(x_array_diff*diameter)      
                
                fig_handle = plt.figure()

                plt.xlabel('y')
                plt.ylabel('shear_rate_xy_average_over_z')
                plt.plot(x_array, vx_array)
                plt.tight_layout()

                
                
                fig_handle.savefig(f_shearratexyaveoverz_y_ave_path_nve + str(int(step)) + "." + figformat, format=figformat)
                if ifpickle:
                    # Save figure handle to disk
                    with open(f_shearratexyaveoverz_y_ave_path_nve + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                        pickle.dump(fig_handle, f)
                plt.close('all')
        else:
            fig_handle = plt.figure()
            for step in picksteplist:
                index = int((step-step_mean_all[0])/(step_mean_all[1]-step_mean_all[0]))
                step = step_mean_all[index]
                vx_array = vx_array_all[index]
                vy_array = vy_array_all[index]
                vx_array_diff = np.diff(vx_array.reshape(-1, n_y).sum(axis=1)/n_y)
                vy_array_diff = np.diff(vy_array.reshape(-1, n_y).sum(axis=1)/n_y)

                vx_array = vx_array_diff/(x_array_diff*diameter)
                vy_array = vy_array_diff/(x_array_diff*diameter)
                time = step*float(rr.logfile["ts"])
                plt.plot(x_array, vx_array, label=" time="+ "{:.2E}".format(time))
            plt.title("height="+str(rr.logfile["z_length_create_dp_unit"]))
            plt.legend()
            plt.xlabel('y')
            plt.ylabel('shear_rate_xy_average_over_z')
            plt.plot(x_array, vx_array)
            plt.tight_layout()
             
            fig_handle.savefig(f_shearratexyaveoverz_y_ave_path_nve + "from" + str(int(picksteplist[0])) + "to" + str(int(picksteplist[-1])) + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(f_shearratexyaveoverz_y_ave_path_nve + "from" + str(int(picksteplist[0])) + "to" + str(int(picksteplist[-1])) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig_handle, f)
            plt.close('all')
    
    if if_plot_to_last:
        plotchunk_1_shearratexyaveoverz_y_ave(step1_default, step2_default, figformat="png", ifpickle=ifpickle)
    else:
        plotchunk_1_shearratexyaveoverz_y_ave(step1, step2, figformat="png", ifpickle=ifpickle)


def plotchunk_vxaveoverz_y_ave(if_plot_to_last, step1, step2, n_ave, picksteplist, figformat="png", ifpickle=False):
    f_vxaveoverz_y_ave_path_nve = dp.f_vxaveoverz_y_ave_path + "nve_" + str(n_ave) + "/"
    if not os.path.isdir(f_vxaveoverz_y_ave_path_nve): 
        os.mkdir(f_vxaveoverz_y_ave_path_nve)
    index=0
    with open(rc.folder_path_list_last_to_initial[index] + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
        lines = f.read().strip().split('\n')
    n_line_in_a_step = int(lines[3].split()[1])
    header = lines[2].split()[1:]
    step2_default = int(lines[-1 - n_line_in_a_step].split()[0])-n_ave*d_step
    step1_default = 0
    def index_of_variable_in_header(variable_name_in_header):
        return [n for n in range(0,len(header)) if header[n]==variable_name_in_header][0]
    def plotchunk_1_vxaveoverz_y_ave(step1_1, step2_1, figformat="png", ifpickle=False):
        for index in range(rc.n_log_list):


            with open(rc.folder_path_list_last_to_initial[index] + "output/momentum_mass_field/fix.momentum_mass_field.all")as f:
                lines = f.read().strip().split('\n') 
            n_line_in_a_step = int(lines[3].split()[1])
            step_begin = int(lines[3].split()[0])
            step_last = int(lines[-1 - n_line_in_a_step].split()[0])-n_ave*d_step
            if index>=1:
                step_last = next_step1
            if step1_1 < step_begin:
                step1_select = step_begin
            else:
                step1_select = step1_1

            if step2_1 > step_last:
                step2_select = step_last
            else:
                step2_select = step2_1
            step_array = np.arange(step1_select, step2_select, d_step)

            def data_inloop(step_smallloop):
                n_line_0 = (step_smallloop - step_begin)/d_step*(n_line_in_a_step+1) + 4
                n_line_1 = n_line_0 + n_line_in_a_step
                ## select data
                data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
                
                data = np.array(data, dtype=np.float64)
                return data

            data_array = np.empty([step_array.shape[0], n_line_in_a_step, len(header)])
            for k in range(step_array.shape[0]):
                data_array[k,:,:] = data_inloop(step_array[k])
            next_step1 = step1_select
            
            if index == 0:
                step_array_out = step_array
                data_array_out = data_array
            else:
                step_array_out = np.concatenate([step_array,step_array_out])
                data_array_out = np.concatenate([data_array,data_array_out])

            if step1_1>=step_begin:
                break 

        step_sum_all = 0
        for k in range(n_ave):
            step_sum_all += step_array_out[k:step_array_out.shape[0]+1-n_ave+k]
        step_mean_all = step_sum_all/n_ave
        data_sum_all = 0
        for k in range(n_ave):
            data_sum_all += data_array_out[k:step_array_out.shape[0]+1-n_ave+k]
        data_mean_all = data_sum_all/n_ave

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
        elif chunk_method == "yz":
            if rr.logfile["chunk/atom"][1] == "y":
                x_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                y_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                x_array = x_array/diameter
                y_array = y_array/diameter
                
            elif rr.logfile["chunk/atom"][1] == "z":
                x_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                y_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                x_array = x_array/diameter
                y_array = y_array/diameter
            else:
                sys.exit("wrong")
        else:
            sys.exit("chunk_method wrong")

        if chunk_method == "rz":
            vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvt")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
            vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
        elif chunk_method == "yz":
            vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvx")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
            vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
        else:
            sys.exit("chunk_method wrong")
        n_y = np.sum(x_array == x_array[0])
        
        x_array = x_array.reshape(-1, n_y)
        x_array = x_array[:,0]
        
        if picksteplist=="all":
            for index in range(step_mean_all.shape[0]):    
                step = step_mean_all[index]
                vx_array = vx_array_all[index]
                vx_array = vx_array.reshape(-1, n_y).sum(axis=1)/n_y
                fig_handle = plt.figure()

                plt.xlabel('y')
                plt.ylabel('Vx_average_over_z')
                plt.plot(x_array, vx_array)
                plt.tight_layout()
                
                fig_handle.savefig(f_vxaveoverz_y_ave_path_nve + str(int(step)) + "." + figformat, format=figformat)
                if ifpickle:
                    # Save figure handle to disk
                    with open(f_vxaveoverz_y_ave_path_nve + str(int(step)) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                        pickle.dump(fig_handle, f)
                plt.close('all')

        else:
            fig_handle = plt.figure()
            for step in picksteplist:
                index = int((step-step_mean_all[0])/(step_mean_all[1]-step_mean_all[0]))
                step = step_mean_all[index]
                vx_array = vx_array_all[index]
                vx_array = vx_array.reshape(-1, n_y).sum(axis=1)/n_y
                time = (step)*float(rr.logfile["ts"])
                plt.plot(x_array, vx_array, label=" time=" + "{:.2E}".format(time))
            plt.title("height="+str(rr.logfile["z_length_create_dp_unit"]))
            plt.legend()
            plt.xlabel('y')
            plt.ylabel('Vx_average_over_z')
            plt.tight_layout()  
            fig_handle.savefig(f_vxaveoverz_y_ave_path_nve + "from" + str(int(picksteplist[0])) + "to" + str(int(picksteplist[-1])) + "." + figformat, format=figformat)
            if ifpickle:
                # Save figure handle to disk
                with open(f_vxaveoverz_y_ave_path_nve + "from" + str(int(picksteplist[0])) + "to" + str(int(picksteplist[-1])) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                    pickle.dump(fig_handle, f)
            plt.close('all')

    if if_plot_to_last:
        plotchunk_1_vxaveoverz_y_ave(step1_default, step2_default, figformat="png", ifpickle=ifpickle)
    else:
        plotchunk_1_vxaveoverz_y_ave(step1, step2, figformat="png", ifpickle=ifpickle)


def plotchunk_ek_d1_aveover_d2_d3_ave(if_plot_to_last, step1, step2, n_ave, picksteplist, figformat="png", ifpickle=False, inputfilepath=None, d1="y",d2="z",d3="y"):

    f_ek_d1_aveover_d2_d3_ave_path = dp.diagram_path + "ek" + d1 + "aveover" + d2 + "_" + d3 + "_ave/"
    om.create_directory(f_ek_d1_aveover_d2_d3_ave_path)
    f_ek_d1_aveover_d2_d3_ave_path_nve = f_ek_d1_aveover_d2_d3_ave_path + "nve_" + str(n_ave) + "/"
    if not os.path.isdir(f_ek_d1_aveover_d2_d3_ave_path_nve): 
        os.mkdir(f_ek_d1_aveover_d2_d3_ave_path_nve)
    index=0
    if inputfilepath == None:
        with open(rc.folder_path_list_last_to_initial[index] + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
            lines = f.read().strip().split('\n')
    else:
        with open(inputfilepath) as f:
            lines = f.read().strip().split('\n')
    
    n_line_in_a_step = int(lines[3].split()[1])
    header = lines[2].split()[1:]
    step2_default = int(lines[-1 - n_line_in_a_step].split()[0])-n_ave*d_step
    step1_default = 0
    def index_of_variable_in_header(variable_name_in_header):
        return [n for n in range(0,len(header)) if header[n]==variable_name_in_header][0]
    def plotchunk_1_ek_d1_aveover_d2_d3_ave(step1_1, step2_1, figformat="png", ifpickle=False):
        
        if inputfilepath == None:
            with open(rc.folder_path_list_last_to_initial[0] + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
                lines = f.read().strip().split('\n')
        else:
            with open(inputfilepath) as f:
                lines = f.read().strip().split('\n')
        n_lines = len(lines)

        n_line_in_a_step = int(lines[3].split()[1])

        n_timestep = int((n_lines-3)/(n_line_in_a_step+1))

        step_array = np.array([int(lines[3+n*(n_line_in_a_step+1)].split()[0]) for n in range(n_timestep)])

        def data_inloop(k):
            n_line_0 = k*(n_line_in_a_step+1) + 4
            n_line_1 = n_line_0 + n_line_in_a_step
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            
            data = np.array(data, dtype=np.float64)
            return data

        data_array = np.empty([step_array.shape[0], n_line_in_a_step, len(header)])
        
        for k in range(step_array.shape[0]):
            data_array[k,:,:] = data_inloop(k)

        step_array_out = step_array
        data_array_out = data_array
        

        step_sum_all = 0
        for k in range(n_ave):
            step_sum_all += step_array_out[k:step_array_out.shape[0]+1-n_ave+k]
        step_mean_all = (step_sum_all/n_ave)
        data_sum_all = 0
        for k in range(n_ave):
            data_sum_all += data_array_out[k:step_array_out.shape[0]+1-n_ave+k]
        data_mean_all = data_sum_all/n_ave

        if chunk_method == "rz":
            coord1_array = int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit
            n_grid_coord1 = len(coord1_array)
            coord2_array = np.arange(n_2)+0.5/n_2*height_dpunit
            n_grid_coord2 = len(coord2_array)
        elif chunk_method == "yz":
            coord1_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
            n_grid_coord1 = np.sum(coord1_array == coord1_array[0])
            coord1_array = coord1_array.reshape(-1, n_grid_coord1)
            coord1_array = coord1_array[:,0]
            n_grid_coord2 = len(coord1_array)
        else:
            sys.exit("chunk_method wrong")

        if chunk_method == "yz":
            dic_map_coordinate_to_header={}
            if rr.logfile["chunk/atom"][1] == "y":
                dic_map_coordinate_to_header["y"] = "Coord1"
                dic_map_coordinate_to_header["z"] = "Coord2"
            elif rr.logfile["chunk/atom"][1] == "z":
                dic_map_coordinate_to_header["y"] = "Coord2"
                dic_map_coordinate_to_header["z"] = "Coord1"
            else:
                sys.exit("chunk_method wrong")

        axis_number_for_d2 = {}
        if chunk_method == "rz":
            axis_number_for_d2["r"] = 0
            axis_number_for_d2["r"] = 1
        elif chunk_method == "yz":
            if rr.logfile["chunk/atom"][1] == "y":
                axis_number_for_d2["y"] = 0
                axis_number_for_d2["z"] = 1
            elif rr.logfile["chunk/atom"][1] == "z":
                axis_number_for_d2["y"] = 1
                axis_number_for_d2["z"] = 0
            else:
                sys.exit("chunk_method wrong")
        else:
            sys.exit("chunk_method wrong")



        if chunk_method == "rz":
            grid_r = int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_wall_dp_unit
            grid_z = np.arange(n_2)+0.5/n_2*height_dpunit
            if d3 == "r":
                x_array = grid_r
            elif d3 == "z":
                x_array = grid_z
            else:
                sys.exit("chunk_method wrong")
        elif chunk_method == "yz":
            x_array = x_array/diameter
        else:
            sys.exit("chunk_method wrong")

        n_grid_normal_to_shear_plane = n_grid_coord1
        n_grid_z = n_grid_coord2
        
        y_array_all = data_mean_all[:,:,index_of_variable_in_header("v_Ek"+d1)]
        
        
        
        fig_handle = plt.figure()

        for step in picksteplist:
            
            for index_2 in range(len(step_mean_all)):
                if step == step_mean_all[index_2]:
                    break
            if not step == step_mean_all[index_2]:
                breakpoint()
                sys.exit("no match step")
            index = index_2
            step = step_mean_all[index]
            y_array = y_array_all[index]
            number_sum = axis_number_for_d2[d2]
            y_array = y_array.reshape(-1, n_grid_coord1).sum(axis=axis_number_for_d2[d2])/number_sum
            time = (step)*float(rr.logfile["ts"])
            plt.plot(x_array[0:9], y_array[0:9], label=" time=" + "{:.2E}".format(time), linestyle="None", marker=11, markersize=15)
        plt.title("height="+str(rr.logfile["z_length_create_dp_unit"]))
        plt.legend()
        plt.xlabel(d3)
        plt.ylabel('ek'+d1)
        plt.tight_layout()  
        fig_handle.savefig(f_ek_d1_aveover_d2_d3_ave_path_nve + "from" + str(int(picksteplist[0])) + "to" + str(int(picksteplist[-1])) + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(f_ek_d1_aveover_d2_d3_ave_path_nve + "from" + str(int(picksteplist[0])) + "to" + str(int(picksteplist[-1])) + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig_handle, f)
        plt.close('all')

    if if_plot_to_last:
        plotchunk_1_ek_d1_aveover_d2_d3_ave(step1_default, step2_default, figformat="png", ifpickle=ifpickle)
    else:
        plotchunk_1_ek_d1_aveover_d2_d3_ave(step1, step2, figformat="png", ifpickle=ifpickle)


def plotchunk_ekyaveoverz_y_ave(if_plot_to_last, step1, step2, n_ave, picksteplist, figformat="png", ifpickle=False, inputfilepath=None):
    plotchunk_ek_d1_aveover_d2_d3_ave(if_plot_to_last, step1, step2, n_ave, picksteplist, figformat="png", ifpickle=False, inputfilepath=None)



def plotchunk_vxaveoverz_y_ave_combine_diff_simu(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False):
    lmp_path_list = [
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_240_Sa_2e-6_his_yes_xmu_5e-1_noshear_2e7/f_2e7_run_5e7/",
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_120_Sa_2e-6_his_yes_xmu_5e-1_noshear_2e7/f_2e7_run_5e7/",
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_60_Sa_2e-6_his_yes_xmu_5e-1_noshear_until_echeck_1e-17/f_5e6_run_1e7/",
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_30_Sa_2e-6_his_yes_xmu_5e-1_noshearfor_5e6_and_run_1e7/f_5e6_run_1e7/",
]
    start_step = [200000000,20000000,5000000,5000000]
    step1_list = [step+9000000 for step in start_step]
    step2_list = [step+9010000 for step in start_step]
    height_list = [240,120,60,30]

    fig_handle = plt.figure()
    for k, lmp_path in enumerate(lmp_path_list):
        picksteplist = [step1_list[k]]
        f_vxaveoverz_y_ave_path_nve = dp.lammps_directory
        post_process_folder_paths = [
            f_vxaveoverz_y_ave_path_nve
        ]
        for post_process_folder_path in post_process_folder_paths:
            if not os.path.isdir(post_process_folder_path): 
                    os.mkdir(post_process_folder_path)

        with open(lmp_path + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
            
            lines = f.read().strip().split('\n')
            header = lines[2].split()[1:]
            n_line_in_a_step = int(lines[3].split()[1])
            step1_default = int(lines[3].split()[0])
            step2_default = int(lines[-1 - n_line_in_a_step].split()[0])-n_ave*d_step
            def index_of_variable_in_header(variable_name_in_header):
                return [n for n in range(0,len(header)) if header[n]==variable_name_in_header][0]
        def plotchunk_1_vxaveoverz_y_ave(step1_1, step2_1, figformat="png", ifpickle=False):
            
            def data_inloop(step_smallloop):
                n_line_0 = (step_smallloop - step1_1)/d_step*(n_line_in_a_step+1) + 4
                n_line_1 = n_line_0 + n_line_in_a_step
                ## select data
                data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
                data = np.array(data, dtype=np.float64)
                return data
            
            
            step_array = np.arange(step1_1, step2_1, d_step)
            data_array = np.empty([step_array.shape[0], n_line_in_a_step, len(header)])
            for index in range(step_array.shape[0]):
                data_array[index,:,:] = data_inloop(step_array[index])

            step_sum_all = 0
            for index in range(n_ave):
                step_sum_all += step_array[index:step_array.shape[0]+1-n_ave+index]
            step_mean_all = step_sum_all/n_ave
            data_sum_all = 0
            for index in range(n_ave):
                data_sum_all += data_array[index:step_array.shape[0]+1-n_ave+index]
            data_mean_all = data_sum_all/n_ave

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
            elif chunk_method == "yz":
                if rr.logfile["chunk/atom"][1] == "y":
                    x_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                    y_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                    
                elif rr.logfile["chunk/atom"][1] == "z":
                    x_array = data_mean_all[0][:,index_of_variable_in_header("Coord2")]
                    y_array = data_mean_all[0][:,index_of_variable_in_header("Coord1")]
                    x_array = x_array/diameter
                    y_array = y_array/diameter
                else:
                    sys.exit("wrong")
            else:
                sys.exit("chunk_method wrong")

            if chunk_method == "rz":
                vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvt")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
                vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
            elif chunk_method == "yz":
                vx_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvx")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
                vy_array_all = divide_zero(data_mean_all[:,:,index_of_variable_in_header("v_mvz")],data_mean_all[:,:,index_of_variable_in_header("c_m1")])/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            n_y = np.sum(x_array == x_array[0])
            
            x_array = x_array.reshape(-1, n_y)
            x_array = x_array[:,0]
            
            
                
            step = step1_list[k]
            index = int((step-step_mean_all[0])/(step_mean_all[1]-step_mean_all[0]))
            step = step_mean_all[index]
            vx_array = vx_array_all[index]
            vx_array = vx_array.reshape(-1, n_y).sum(axis=1)/n_y
            time = (step-start_step[k])*float(rr.logfile["ts"])
            plt.plot(x_array, vx_array, label="height="+str(height_list[k]))
            plt.title("time="+ "{:.2E}".format(time))
               
        plotchunk_1_vxaveoverz_y_ave(step1_list[k], step2_list[k], figformat="png", ifpickle=ifpickle)
    
    plt.legend()
    plt.tight_layout()     
    plt.xlabel('y')
    plt.ylabel('Vx_average_over_z')

    if picksteplist=="all":
        fig_handle.savefig(f_vxaveoverz_y_ave_path_nve + "combine_diff_height" + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(f_vxaveoverz_y_ave_path_nve + "combine_diff_height" + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig_handle, f)
        plt.close('all')

    else:
        
        fig_handle.savefig(f_vxaveoverz_y_ave_path_nve + "combine_diff_height" + "." + figformat, format=figformat)
        if ifpickle:
            # Save figure handle to disk
            with open(f_vxaveoverz_y_ave_path_nve + "combine_diff_height" + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig_handle, f)
    plt.close('all')