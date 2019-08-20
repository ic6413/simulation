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

step_0 = 0
d_step = 100
n_r = 17
n_z = 10
label_x = 'v_vr'
label_y = 'vz'

def plotchunk(step, file):
    n_line_0 = (step - step_0)/d_step*(n_r*n_z+1) + 4
    n_line_1 = n_line_0 + n_r*n_z
    x_array, y_array = np.meshgrid(np.arange(n_r), np.arange(n_z))

    with open(file) as f:
        
        lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        ## select data
        data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
        ## attach data
        df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        ## repeat timestep

    vx_array = df[label_x].values
    vy_array = df[label_y].values
    fig1, ax1 = plt.subplots()
    plt.xlabel('r')
    plt.ylabel('z')
    ax1.set_title('velocity field in r-z direction (average over theta)')
    Q = ax1.quiver(x_array, y_array, vx_array, vy_array, units='width')
    fig1.savefig(str(step), bbox_inches='tight')
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

