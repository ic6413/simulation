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
# define function for extract data from fix txt to dataframe

velocity_scale = dp.omega_in*dp.r_in

width_dpunit = 17 #(dp.r_out - dp.r_in)/dp.dp0
height_dpunit = 37

step_0 = 10000
d_step = 10000

n_r = dp.N_bin_r
n_z = dp.N_bin_z



def plotchunk(step, file):
    quiver_scale = 10
    label_scale = 10
    n_line_0 = (step - step_0)/d_step*(n_r*n_z+1) + 4
    n_line_1 = n_line_0 + n_r*n_z
    x_array, y_array = np.meshgrid(np.arange(n_r)/n_r*width_dpunit, np.arange(n_z)/n_z*height_dpunit)

    with open(file) as f:
        
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        ## select data
        data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
        ## attach data
        df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        ## repeat timestep
    label_x = 'v_vr'
    label_y = 'vz'
    vx_array = df[label_x].values
    vy_array = df[label_y].values
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
    label_x = 'v_vt'
    label_y = 'vz'
    vx_array = df[label_x].values
    vy_array = df[label_y].values
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

