# import site package
import os
import sys
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datapath as dp
# import my package
import read_setting as rr
diameter = float(rr.log_variable["dp"])
lmp_path = rr.folder_path_list_initial_to_last[-1]
# python read line
def extractdata(timestep, variablename):
    filepath_rela_lmp_path = "output/single_all/dump.all.single." + str(timestep)
    
    with open(lmp_path + filepath_rela_lmp_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    n_line_in_file = len(lines)
    first_9_lines = lines[0: 9]
    # header
    header = first_9_lines[8].split()[2:]

    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[9:]

    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    index_id = header.index('id')
    index_variable = header.index(variablename)

    id = lines[:, index_id]
    ind = np.argsort(id)
    value = lines[:, index_variable]
    value = np.take_along_axis(value, ind, axis=0)
    return value

# max distance
def distance(s1, s2):
    dx = extractdata(s1, 'x') - extractdata(s2, 'x')
    dy = extractdata(s1, 'y') - extractdata(s2, 'y')
    dz = extractdata(s1, 'z') - extractdata(s2, 'z')
    dyz = (dy**2+dz**2)**0.5
    return [extractdata(s1, 'y'), extractdata(s1, 'z'), dyz]

def plotstaticzone(s1, s2, rate_lessthanave=0.01):
    [y, z, dyz] = distance(s1, s2)
    ave = np.average(dyz)
    ifstatic = (dyz < rate_lessthanave*ave)
    nonstatic = np.logical_not(ifstatic)
    static_y = ifstatic*y
    static_z = ifstatic*z
    nonstatic_y = nonstatic*y
    nonstatic_z = nonstatic*z
    # normalize
    static_y = static_y/diameter
    static_z = static_z/diameter
    nonstatic_y = nonstatic_y/diameter
    nonstatic_z = nonstatic_z/diameter
    
    fig, ax = plt.subplots()
    ax.set_title(
        "staticzone \n(" + "moving distance smaller than " + str(rate_lessthanave) + "*average)" + "\n step " + str(s1) + " to " + str(s2)
    )
    ax.plot(
        static_y, static_z,
        marker = ".",
        linestyle = 'None',
        markersize=1,
    )
    
    ax.plot(
        nonstatic_y, nonstatic_z,
        marker = ".",
        linestyle = 'None',
        markersize=1,
    )
    
    subfoldername = "staticzone/"
    os.makedirs(dp.diagram_path + subfoldername, exist_ok=True)
    fig.savefig(
    "_".join([
        dp.diagram_path + subfoldername + 'step',
        str(s1),
        str(s2),
    ]),
    format="png",
    )

# define static region
mask = (y-9) - (9-0.5)/(15.5-10.5)*(x-15.5) < 0

# define nonstatic region
mask = (y-25) - (25-0.5)/(15.5-5)*(x-5) > 0

# calculate variable for static region and nonstatic region

