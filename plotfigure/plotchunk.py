# import
import os.path
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

velocity_scale = dp.omega_in*dp.r_in

width_dpunit = int(rr.logfile['ro_dp_unit'])-int(rr.logfile['ri_dp_unit'])
height_dpunit = float(rr.logfile['zhi_chunk_dp_unit'])
diameter = float(rr.logfile['dp'])

step_0 = 10000
d_step = 10000

n_r = dp.N_bin_r
n_z = dp.N_bin_z
x_array, y_array = np.meshgrid(
                               int(rr.logfile['ri_dp_unit']) + (np.arange(n_r)+0.5)/n_r*width_dpunit,
                               (np.arange(n_z)+0.5)/n_z*height_dpunit,
                               )
dx = 1/n_r*width_dpunit
dy = 1/n_z*height_dpunit
x_array = x_array.reshape((-1))
y_array = y_array.reshape((-1))

vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*diameter**3


def plotchunk(step, file):
    
    n_line_0 = (step - step_0)/d_step*(n_r*n_z+1) + 4
    n_line_1 = n_line_0 + n_r*n_z
    

    with open(file) as f:
        
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        ## select data
        data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
        ## attach data
        df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        ## repeat timestep

    quiver_scale = 10
    label_scale = 10
    vx_array = np.divide(df['v_mvr'].values, df['c_m1'].values)
    vy_array = np.divide(df['v_mvz'].values, df['c_m1'].values)
    fig1, ax1 = plt.subplots()
    plt.xlabel('r')
    plt.ylabel('z')
    #ax1.set_title('velocity field r-z direction (average over theta)')
    Q = ax1.quiver(x_array, y_array, vx_array/velocity_scale, vy_array/velocity_scale,
                   units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                   )

    ax1.quiverkey(Q, 0.2, 0.9, label_scale, label = str(label_scale)+'*arrow length of wall velocity in 45 degree, velocity field r-z direction', labelpos='E',
                   coordinates='figure', angle=45)

    fig1.savefig(dp.f_momentum_mass_field_rz_path + str(step))
    plt.close('all')

    quiver_scale = 100
    label_scale = 100

    vx_array = np.divide(df['v_mvt'].values, df['c_m1'].values)
    vy_array = np.divide(df['v_mvz'].values, df['c_m1'].values)
    fig1, ax1 = plt.subplots()
    plt.xlabel('r(position), theta(velocity)')
    plt.ylabel('z')
    #ax1.set_title('velocity field r-z direction (average over theta)')
    Q = ax1.quiver(x_array, y_array, vx_array/velocity_scale, vy_array/velocity_scale,
                   units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                   )
    ax1.quiverkey(Q, 0.2, 0.9, label_scale, label = str(label_scale)+'*arrow length of wall velocity in 45 degree, velocity field r-theta direction', labelpos='E',
                   coordinates='figure', angle=45)

    fig1.savefig(dp.f_momentum_mass_field_rtheta_path + str(step))
    plt.close('all')


    quiver_scale = 0.2
    label_scale = 1
    vx_array = 0*df['c_m1'].values
    vy_array = df['c_m1'].values/float(rr.logfile['den'])/vol_in_chunks
    fig1, ax1 = plt.subplots()
    plt.xlabel('r')
    plt.ylabel('z')
    #ax1.set_title('velocity field r-z direction (average over theta)')
    Q = ax1.quiver(x_array, y_array, vx_array, vy_array,
                   units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                   )

    ax1.quiverkey(Q, 0.2, 0.9, label_scale, label = 'solid fraction', labelpos='E',
                   coordinates='figure', angle=90)

    fig1.savefig(dp.f_momentum_mass_field_density_path + str(step))
    plt.close('all')


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

