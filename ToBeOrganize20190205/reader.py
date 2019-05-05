import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

filename = sys.argv[1]
f = open(filename)
lines = f.read().strip().split('\n')
header = lines[2][2:].split(' ')
lines = lines[3:]
idx = 0
while idx < len(lines):
    data = []
    time = lines[idx].split(' ')[0]
    rownum = int(lines[idx].split(' ')[1])
    idx += 1
    for j in range(rownum):
        data.append(lines[idx][2:].split(' '))
        idx += 1

    datapd = pd.DataFrame(data=data, columns=header)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make the grid
    x = datapd['Coord1'].values.astype('float64')
    y = datapd['Coord2'].values.astype('float64')
    z = datapd['Coord3'].values.astype('float64')

    # Make the direction data for the arrows
    u = datapd['vx'].values.astype('float64')
    v = datapd['vy'].values.astype('float64')
    w = datapd['vz'].values.astype('float64')

    ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)

    fig.savefig(sys.argv[2] + time + '.png')   # save the figure to file
    plt.close(fig)    # close the figure


