import sys
import pprint
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import simulation module
import datapath as dp

# import plot module
import plotfigure.inputvariable as pi
# inputvariable
step1 = pi.step1
step2 = pi.step2

import plot

#=====
def plot_v(array):
	x = np.arange(step1, step2+1)
	plt.plot(x, array)
	x_label = 'Step'
	plt.xlabel(x_label)
	plt.xscale('linear')
	plt.yscale('linear', linthreshy=1)
	plt.xticks(rotation=20, fontsize=14)
	plt.yticks(fontsize=14)
#=====

outputfolder = dp.debug_print_path + 'fig/'
f_read = dp.f_thermo + '.h5'
variable_name_list = 'all' # 'all'

df = pd.read_hdf(f_read)

if variable_name_list == 'all':
    variable_name_list = list(df)


df_step = df.loc[df['Step'].isin(list(range(step1, step2+1)))]

for variable_name in variable_name_list:

    y = df_step[variable_name].values
    fig = plt.figure()
    plot_v(y)
    plt.ylabel(variable_name)
    outputfig = outputfolder + variable_name + '_step_' + str(step1) + '_' + str(step2) + '.png'
    fig.savefig(outputfig, bbox_inches='tight')   # save the figure to file
