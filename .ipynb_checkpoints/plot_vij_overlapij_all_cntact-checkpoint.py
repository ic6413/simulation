import sys
import pprint
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import project module
import check_simulation as cs
import output_control as oc
import datapath as dp

outputfolder = dp.debug_print_path + 'fig/'
f_read = dp.f_custom + '.h5'
id_i = oc.id_i

method = 1

[number_contact_total, contact_id_collection_no_dup] = cs.number_contact_total_id_collection(f_read, id_i, oc.step1, oc.step2, oc.error_tolerence)

def plot_v(array):
	x = np.arange(oc.step1, oc.step2+1)
	plt.plot(x, array)
	x_label = 'Step'
	plt.xlabel(x_label)
	plt.xscale('linear')
	plt.yscale('linear')
	plt.xticks(rotation=20, fontsize=14)
	plt.yticks(fontsize=14)
	plt.title(str(idj_or_idw))


def plot_nth_column_v(array2D, n):
	array = array2D[:,n]
	plot_v(array)


def plot_length_v(array2D):
	array = cs.length(array2D)[:,0]
	plot_v(array)



for idj_or_idw in contact_id_collection_no_dup:
	[vijt, vijn] = cs.vij_many_steps(id_i, idj_or_idw, f_read, oc.step1, oc.step2, method)

	fig = plt.figure()
	plot_length_v(vijt)
	y_label = 'vijt'
	plt.ylabel(y_label)
	outputfig = outputfolder + 'vijt_' + 'i_' +str(id_i) + '_j_' + str(idj_or_idw) + '_step_' + str(oc.step1) + '_' + str(oc.step2) + '.png'
	fig.savefig(outputfig, bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_length_v(vijn)
	y_label = 'vijn'
	plt.ylabel(y_label)
	outputfig = outputfolder + 'vijn_' + 'i_' +str(id_i) + '_j_' + str(idj_or_idw) + '_step_' + str(oc.step1) + '_' + str(oc.step2) + '.png'
	fig.savefig(outputfig, bbox_inches='tight')   # save the figure to file

	
	overlap_length = cs.overlap_length_many_steps(id_i, idj_or_idw, f_read, oc.step1, oc.step2, method)

	fig = plt.figure()
	plot_length_v(overlap_length)
	y_label = 'overlap_length'
	plt.ylabel(y_label)
	outputfig = outputfolder + 'ol_' + 'i_' +str(id_i) + '_j_' + str(idj_or_idw) + '_step_' + str(oc.step1) + '_' + str(oc.step2) + '.png'
	fig.savefig(outputfig, bbox_inches='tight')   # save the figure to file


