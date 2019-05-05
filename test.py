#%%

import sys
import pprint
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# import project module
import check_simulation as cs
import output_control as oc
import datapath as dp

outputfolder = dp.debug_print_path + 'fig/'
f_read = dp.f_custom + '.h5'
id_i = oc.id_i

method = 1

#[number_contact_total, contact_id_collection_no_dup] = cs.number_contact_total_id_collection(f_read, id_i, oc.step1, oc.step2, oc.error_tolerence)

def plot_v(array):

    x = np.arange(oc.step1, oc.step2+1)
    #ax = plt.semilogy(x, array)
    fig, ax = plt.subplots()
    ax.plot(x, array)
    x_label = 'Step'
    plt.xlabel(x_label)
    plt.xscale('linear')
    plt.yscale('symlog')
    plt.xticks(rotation=20, fontsize=14)
    plt.yticks(fontsize=14)

def add_id_title(idj_or_idw):
	plt.title('idj = {}'.format(str(idj_or_idw)))
	
def add_y_label(ylabel_name):
	plt.ylabel(ylabel_name)

def outputfig_path(ylabel_name, idj_or_idw):
	outputfig = outputfolder + ylabel_name + '_i_' +str(id_i) + '_j_' + str(idj_or_idw) + '_step_' + str(oc.step1) + '_' + str(oc.step2) + '.png'
	return outputfig



def plot_nth_column_v(array2D, n):
	array = array2D[:,n]
	plot_v(array)


def plot_length_v(array2D):
	array = cs.length(array2D)[:,0]
	plot_v(array)


def plot_1_id(idj_or_idw):
	
	# velocity
	[vijt, vijn] = cs.vij_many_steps(id_i, idj_or_idw, f_read, oc.step1, oc.step2, method)

	fig = plt.figure()
	plot_length_v(vijt)
	add_id_title(idj_or_idw)
	add_y_label('vijt')

	fig.savefig(outputfig_path('vijt', idj_or_idw), bbox_inches='tight')   # save the figure to file




idj_or_idw=15583
# velocity
[vijt, vijn] = cs.vij_many_steps(id_i, idj_or_idw, f_read, oc.step1, oc.step2, method)

array = cs.length(vijt)[:,0]

x = np.arange(oc.step1, oc.step2+1)
#ax = plt.semilogy(x, array)
fig, ax = plt.subplots()
ax.plot(x, array)
x_label = 'Step'

plt.xlabel(x_label)
plt.xscale('linear')
plt.yscale('symlog', linthreshy=1e-20)
plt.xticks(rotation=20, fontsize=14)
plt.yticks(fontsize=14)

ylim = ax.get_ylim()

#plt.yticks([ylim[1]/2,ylim[1]])
# from matplotlib.ticker import MultipleLocator


#ml = MultipleLocator(1)

#ax.yaxis.set_major_locator(MultipleLocator(0.01))
#ax.yaxis.set_minor_locator(ml)


#ax.axis([ 50000, 300000,0, 0.01])

add_id_title(idj_or_idw)
add_y_label('vijt')

fig.savefig(outputfig_path('vijt', idj_or_idw), bbox_inches='tight')   # save the figure to file


#%%
from plot_vij_overlapij_ftn_all_cntact import plot_1_id
plot_1_id(15583)
#%%
