import sys
import pprint
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import simulation module
import datapath as dp

# import plot module
import plotfigure.input as pi

# import calculate module
import calculate.checksimulation as cs


# input
id_i = pi.id_i
step1 = pi.step1
step2 = pi.step2
error_tolerence = pi.error_tolerence
# begin


outputfolder = dp.debug_print_path + 'fig/'
f_read = dp.f_custom + '.h5'
id_i = id_i

method = 1

#[number_contact_total, contact_id_collection_no_dup] = cs.number_contact_total_id_collection(f_read, id_i, step1, step2, error_tolerence)

def plot_v(array):
	x = np.arange(step1, step2+1)
	plt.plot(x, array)
	x_label = 'Step'
	plt.xlabel(x_label)
	plt.xscale('linear')
	plt.yscale('symlog', linthreshy=1e-20)
	plt.xticks(rotation=20, fontsize=14)
	plt.yticks(fontsize=14)

def add_id_title(idj_or_idw):
	plt.title('idj = {}'.format(str(idj_or_idw)))
	
def add_y_label(ylabel_name):
	plt.ylabel(ylabel_name)

def outputfig_path(ylabel_name, idj_or_idw):
	outputfig = outputfolder + ylabel_name + '_i_' +str(id_i) + '_j_' + str(idj_or_idw) + '_step_' + str(step1) + '_' + str(step2) + '.png'
	return outputfig



def plot_nth_column_v(array2D, n):
	array = array2D[:,n]
	plot_v(array)


def plot_length_v(array2D):
	array = cs.length(array2D)[:,0]
	plot_v(array)


def plot_1_id(idj_or_idw):
	
	# velocity
	[vijt, vijn] = cs.vij_many_steps(id_i, idj_or_idw, f_read, step1, step2, method)

	fig = plt.figure()
	plot_length_v(vijt)
	add_id_title(idj_or_idw)
	add_y_label('vijt')
	fig.savefig(outputfig_path('vijt', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_length_v(vijn)
	add_id_title(idj_or_idw)
	add_y_label('vijn')
	fig.savefig(outputfig_path('vijn', idj_or_idw), bbox_inches='tight')   # save the figure to file

	# overlap
	overlap_length = cs.overlap_length_many_steps(id_i, idj_or_idw, f_read, step1, step2, method)

	fig = plt.figure()
	plot_length_v(overlap_length)
	add_id_title(idj_or_idw)
	add_y_label('overlap_length')
	fig.savefig(outputfig_path('overlap_length', idj_or_idw), bbox_inches='tight')   # save the figure to file

	# force
	[fnk,fngamma,ftk_include_his,ftgamma] = cs.f_many_steps(f_read, id_i, idj_or_idw, step1, step2, error_tolerence, method)
	fn = fnk + fngamma
	ft = ftk_include_his + ftgamma

	
	fig = plt.figure()
	plot_length_v(fnk)
	add_id_title(idj_or_idw)
	add_y_label('fnk')
	fig.savefig(outputfig_path('fnk', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_length_v(fngamma)
	add_id_title(idj_or_idw)
	add_y_label('fngamma')
	fig.savefig(outputfig_path('fngamma', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_length_v(ftk_include_his)
	add_id_title(idj_or_idw)
	add_y_label('ftk_include_his')
	fig.savefig(outputfig_path('ftk_include_his', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_length_v(ftgamma)
	add_id_title(idj_or_idw)
	add_y_label('ftgamma')
	fig.savefig(outputfig_path('ftgamma', idj_or_idw), bbox_inches='tight')   # save the figure to file

	# force
	fig = plt.figure()
	plot_length_v(fn)
	add_id_title(idj_or_idw)
	add_y_label('fn')
	fig.savefig(outputfig_path('fn', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_length_v(ft)
	add_id_title(idj_or_idw)
	add_y_label('ft')
	fig.savefig(outputfig_path('ft', idj_or_idw), bbox_inches='tight')   # save the figure to file

	# work
	[work_ft, work_fn] = cs.work_ftfn_many_steps(f_read, id_i, idj_or_idw, step1, step2, error_tolerence, method)
	fig = plt.figure()
	plot_v(work_ft)
	add_id_title(idj_or_idw)
	add_y_label('work_ft')
	fig.savefig(outputfig_path('work_ft', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_v(work_fn)
	add_id_title(idj_or_idw)
	add_y_label('work_fn')
	fig.savefig(outputfig_path('work_fn', idj_or_idw), bbox_inches='tight')   # save the figure to file

	[work_ftgamma, work_fngamma] = cs.work_ftgammafngamma_many_steps(f_read, id_i, idj_or_idw, step1, step2, error_tolerence, method)
	fig = plt.figure()
	plot_v(work_ftgamma)
	add_id_title(idj_or_idw)
	add_y_label('work_ftgamma')
	fig.savefig(outputfig_path('work_ftgamma', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_v(work_fngamma)
	add_id_title(idj_or_idw)
	add_y_label('work_fngamma')
	fig.savefig(outputfig_path('work_fngamma', idj_or_idw), bbox_inches='tight')   # save the figure to file


	[sum_ftwork, sum_fnwork] = cs.cumsum_work(f_read, id_i, idj_or_idw, step1, step2, error_tolerence, method)
	fig = plt.figure()
	plot_v(sum_ftwork)
	add_id_title(idj_or_idw)
	add_y_label('sum_ftwork')
	fig.savefig(outputfig_path('sum_ftwork', idj_or_idw), bbox_inches='tight')   # save the figure to file

	fig = plt.figure()
	plot_v(sum_fnwork)
	add_id_title(idj_or_idw)
	add_y_label('sum_fnwork')
	fig.savefig(outputfig_path('sum_fnwork', idj_or_idw), bbox_inches='tight')   # save the figure to file

	plt.close('all')



#for idj_or_idw in contact_id_collection_no_dup:
#	fig = plt.figure()	
#	plot_1_id(idj_or_idw)
#	add_id_title(idj_or_idw)
