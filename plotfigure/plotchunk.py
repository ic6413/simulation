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
import rename_variable as rv
import datapath as dp
import read_setting.read_setting as rr
# define function for extract data from fix txt to dataframe

if rr.logfile["chunk/atom"][0] == "chunk_r_z":
    chunk_method = 'rz'
if rr.logfile["chunk/atom"][0] == "chunk_y_z":
    chunk_method = 'yz'
    
diameter = float(rr.logfile['dp'])
width_dp_unit = int(rr.logfile['width_dp_unit'])
if chunk_method == "rz":
    ri = diameter*int(rr.logfile['ri_dp_unit']) 
elif chunk_method == "yz":
    x_period = diameter*int(rr.logfile['x_period_dp_unit'])
else:
    sys.exit("chunk_method wrong")
g = float(rr.logfile['g'])
d_step = int(rr.logfile['freq_ave_chunk_momentum_mass_field'])
velocity_scale = float(rr.logfile['in_velocity'])
if velocity_scale < 0:
    sys.exit("velocity scale smaller than zero")

if velocity_scale == 0:
    Sa_fake = 0.000002
    if chunk_method == "rz":
        velocity_scale = ri*(Sa_fake*g*width_dp_unit**3*diameter/ri**2)**0.5
    elif chunk_method == "yz":
        velocity_scale = (Sa_fake*g*width_dp_unit**3*diameter)**0.5

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
                                int(rr.logfile['ri_dp_unit']) + (np.arange(n_1)+0.5)/n_1*width_dp_unit,
                                (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                )
elif chunk_method == "yz":
    if rr.logfile["chunk/atom"][1] == "y":
        y_array, x_array = np.meshgrid(
                                       (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                       (np.arange(n_1)+0.5)/n_1*width_dp_unit,
                                    )
    elif rr.logfile["chunk/atom"][1] == "z":
        x_array, y_array = np.meshgrid(
                                       (np.arange(n_1)+0.5)/n_1*width_dp_unit,
                                       (np.arange(n_2)+0.5)/n_2*height_dpunit,
                                    )
    else:
        sys.exit("wrong")
    
else:
    sys.exit("chunk_method wrong")

dx = 1/n_1*width_dp_unit
dy = 1/n_2*height_dpunit
x_array = x_array.reshape((-1))
y_array = y_array.reshape((-1))
if chunk_method == "rz":
    vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*diameter**3
elif chunk_method == "yz":
    vol_in_chunks = x_period*dx*dy*diameter**2
else:
    sys.exit("chunk_method wrong")


def plotchunk(if_plot_to_last, step1, step2):

    with open(dp.lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
        
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        step1_default = int(lines[3].split()[0])
        step2_default = int(lines[-1 - n_12].split()[0])
        
    def plotchunk_1(step1_1, step2_1):
        for step in range(step1_1, step2_1, d_step):
            n_line_0 = (step - step1_1)/d_step*(n_12+1) + 4
            n_line_1 = n_line_0 + n_12
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            ## attach data
            df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            ## repeat timestep

            def divide_zero(a,b):
                c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
                return c

            if chunk_method == "rz":
                vx_array = divide_zero(df['v_mvr'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                vx_array = divide_zero(df['v_mvy'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            
            v_length_array = (vx_array**2+vy_array**2)**0.5
            max_v_length = np.amax(v_length_array)
            quiver_scale = max_v_length/2
            label_scale = max_v_length/2
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

            ax1.quiverkey(Q, 0.2, 0.9, label_scale,
                        label = "labescale" + str(label_scale) + "quiver_scale" + str(quiver_scale) + '*arrow length of wall velocity in 45 degree, velocity field r-z direction',
                        labelpos='E', coordinates='figure', angle=45)

            fig1.savefig(dp.f_momentum_mass_field_rz_path + str(step))
            plt.close('all')

            if chunk_method == "rz":
                vx_array = divide_zero(df['v_mvt'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                vx_array = divide_zero(df['v_mvx'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            v_length_array = (vx_array**2+vy_array**2)**0.5
            max_v_length = np.amax(v_length_array)
            quiver_scale = max_v_length/2
            label_scale = max_v_length/2
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
            ax1.quiverkey(Q, 0.2, 0.9, label_scale,
                        label = "labescale" + str(label_scale) + "quiver_scale" + str(quiver_scale) + '*arrow length of wall velocity in 45 degree, velocity field r-z direction',
                        labelpos='E', coordinates='figure', angle=45)

            fig1.savefig(dp.f_momentum_mass_field_rtheta_path + str(step))
            plt.close('all')


            quiver_scale = 0.2
            label_scale = 0.6
            vx_array = 0*df['c_m1'].values
            vy_array = df['c_m1'].values/float(rr.logfile['den'])/vol_in_chunks
            fig1, ax1 = plt.subplots()
            fig1.set_size_inches(12.8, 9.6)
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

            ax1.quiverkey(Q, 0.2, 0.9,
                        label_scale, label = "labescale" + str(label_scale) + "quiver_scale" + str(quiver_scale) + 'solid fraction',
                        labelpos='E',
                        coordinates='figure', angle=90)

            fig1.savefig(dp.f_momentum_mass_field_density_path + str(step))
            plt.close('all')
    
    if if_plot_to_last:
        plotchunk_1(step1_default, step2_default)
    else:
        plotchunk_1(step1, step2)


def plotchunk_ave(if_plot_to_last, step1, step2, n_ave):
    f_momentum_mass_field_rtheta_path_nve = dp.f_momentum_mass_field_rtheta_path + "nve_" + str(n_ave) + "/"
    f_momentum_mass_field_rz_path_nve = dp.f_momentum_mass_field_rz_path + "nve_" + str(n_ave) + "/"
    f_momentum_mass_field_density_path_nve = dp.f_momentum_mass_field_density_path + "nve_" + str(n_ave) + "/"
    post_process_folder_paths = [
        f_momentum_mass_field_rtheta_path_nve,
        f_momentum_mass_field_rz_path_nve,
        f_momentum_mass_field_density_path_nve,
    ]
    for post_process_folder_path in post_process_folder_paths:
        if not os.path.isdir(post_process_folder_path): 
                os.mkdir(post_process_folder_path)

    with open(dp.lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all") as f:
        
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        step1_default = int(lines[3].split()[0])
        step2_default = int(lines[-1 - n_12].split()[0])-n_ave*d_step
        
    def plotchunk_1(step1_1, step2_1):
        def data_inloop(step_smallloop):
            n_line_0 = (step_smallloop - step1_1)/d_step*(n_12+1) + 4
            n_line_1 = n_line_0 + n_12
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
                vx_array = divide_zero(df['v_mvr'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                vx_array = divide_zero(df['v_mvy'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            
            v_length_array = (vx_array**2+vy_array**2)**0.5
            max_v_length = np.amax(v_length_array)
            quiver_scale = max_v_length/2
            label_scale = max_v_length/2
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

            ax1.quiverkey(Q, 0.2, 0.9, label_scale,
                        label = "labescale" + str(label_scale) + "quiver_scale" + str(quiver_scale) + '*arrow length of wall velocity in 45 degree, velocity field r-z direction',
                        labelpos='E', coordinates='figure', angle=45)
            fig1.savefig(f_momentum_mass_field_rz_path_nve + str(int(step)))
            plt.close('all')

            if chunk_method == "rz":
                vx_array = divide_zero(df['v_mvt'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            elif chunk_method == "yz":
                vx_array = divide_zero(df['v_mvx'].values,df['c_m1'].values)/velocity_scale
                vy_array = divide_zero(df['v_mvz'].values,df['c_m1'].values)/velocity_scale
            else:
                sys.exit("chunk_method wrong")
            v_length_array = (vx_array**2+vy_array**2)**0.5
            max_v_length = np.amax(v_length_array)
            quiver_scale = max_v_length/2
            label_scale = max_v_length/2
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
            ax1.quiverkey(Q, 0.2, 0.9, label_scale,
                        label = "labescale" + str(label_scale) + "quiver_scale" + str(quiver_scale) + '*arrow length of wall velocity in 45 degree, velocity field r-z direction',
                        labelpos='E', coordinates='figure', angle=45)

            fig1.savefig(f_momentum_mass_field_rtheta_path_nve + str(int(step)))
            plt.close('all')


            quiver_scale = 0.2
            label_scale = 0.6
            vx_array = 0*df['c_m1'].values
            vy_array = df['c_m1'].values/float(rr.logfile['den'])/vol_in_chunks
            fig1, ax1 = plt.subplots()
            fig1.set_size_inches(12.8, 9.6)
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

            ax1.quiverkey(Q, 0.2, 0.9,
                        label_scale, label = "labescale" + str(label_scale) + "quiver_scale" + str(quiver_scale) + 'solid fraction',
                        labelpos='E',
                        coordinates='figure', angle=90)

            fig1.savefig(f_momentum_mass_field_density_path_nve + str(int(step)))
            plt.close('all')
    
    if if_plot_to_last:
        plotchunk_1(step1_default, step2_default)
    else:
        plotchunk_1(step1, step2)


def chunkfile_to_dataframe(file):

    with open(file) as f:
        lines = f.read().strip().split('\n')
        id_line_timestep = [n for n, line in enumerate(lines) if line.startswith('# Timestep')]
        n_chunks = int(lines[id_line_timestep[0]+2].split()[2])
        header = lines[id_line_timestep[0]+1].split()[1:]
        id_line_timestep.append(len(lines))
        iter = chain.from_iterable(range(id + 3, id_line_timestep[i + 1] - 1) for i, id in enumerate(id_line_timestep[0: -1]))
        ## select data
        data = [lines[t].split() for t in iter]
        ## attach data
        df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        ## repeat timestep
        steps = list(
            chain.from_iterable(
                repeat(
                    lines[id + 2].split()[0], id_line_timestep[i + 1] - -id - 4) for i, id in enumerate(id_line_timestep[0: -1]
                    )
                )
            )
        steps = np.asarray(steps,dtype=np.float64)
        ## insert timesteps
        df.insert(1, 'step', steps)
        
    return df

