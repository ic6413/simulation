import sys
import pprint
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import simulation module
import osmanage as om
import datapath as dp

# import plot module
import plotfigure.inputvariable as pi

# import calculate module
import calculate.checksimulation as cs

om.create_directory(dp.debug_print_path)
om.create_directory(dp.debug_fig_path)
om.create_directory(dp.debug_fig_thermo_path)
om.create_directory(dp.debug_fig_oneatom_path)
om.create_directory(dp.debug_fig_atomij_path)


def create_figclass_i_jorw(step1, step2, array_x, xlabel, array_y, ylabel, id_i, id_jorw):
	if id_jorw > 0:
		result = lammp_figure_atom_ij(step1, step2, array_x, xlabel, array_y, ylabel, id_i, id_jorw)
	elif id_jorw < 0:
		result = lammp_figure_atomi_wall(step1, step2, array_x, xlabel, array_y, ylabel, id_i, id_jorw)
	else:
		sys.exit('idjorw = 0')
	return result

# plot for 1D array_y or length of 2D array_y
class lammp_figure(object):

    def __init__(self, step1, step2, array_x, xlabel, array_y, ylabel):
        self.array_x = array_x
        self.xlabel = xlabel
        self.array_y = array_y
        self.ylabel = ylabel
        self.step1 = step1
        self.step2 = step2
        

    def create_lammp_figure(self):
        
        fig = plt.figure()
        
        def array_to1D(array):
            array_dim = array.ndim
            if array_dim == 1:
                pass
            elif array_dim == 2:
                if array.shape[1] == 1:
                    array = array[:,0]
                else:
                    array = cs.length(array)[:,0]
            else:
                sys.exit('array dim not 1 not 2')
            return array
        
        array_x = array_to1D(self.array_x)
        array_y = array_to1D(self.array_y)

        plt.plot(array_x, array_y)

        plt.xscale('linear')
        plt.yscale('linear')
        #plt.yscale('symlog', linthreshy=1e-20)
        plt.xticks(rotation=20, fontsize=14)
        plt.yticks(fontsize=14)

        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        self.add_title()
        
        return fig
        
    def add_title(self):
        pass

    def filename(self):
        pass

    def create_and_save(self, outputfolder):
        fig = self.create_lammp_figure()
        fig.savefig(outputfolder + self.filename(), bbox_inches='tight')   # save the figure to file

class lammp_figure_thermo(lammp_figure):
    
    def __init__(self, step1, step2, array_x, xlabel, array_y, ylabel):
        super().__init__(step1, step2, array_x, xlabel, array_y, ylabel)
    
    def add_title(self):
        plt.title('thermo')

    def filename(self):
        filename = self.ylabel + '_' + self.xlabel + '_step_' + str(self.step1) + '_' + str(self.step2) + '.png'
        return filename

class lammp_figure_atom_single(lammp_figure):
    
    def __init__(self, step1, step2, array_x, xlabel, array_y, ylabel, atomid):
        super().__init__(step1, step2, array_x, xlabel, array_y, ylabel)
        self.atomid = atomid
    
    def add_title(self):
        plt.title('atom id = {}'.format(str(self.atomid)))

    def filename(self):
        filename = self.ylabel + '_' + self.xlabel + '_id_' + str(self.atomid) + '_step_' + str(self.step1) + '_' + str(self.step2) + '.png'
        return filename

class lammp_figure_atom_ij(lammp_figure):
    
    def __init__(self, step1, step2, array_x, xlabel, array_y, ylabel, id_i, id_j):
        super().__init__(step1, step2, array_x, xlabel, array_y, ylabel)
        self.id_i = id_i
        self.id_j = id_j

    def add_title(self):
        plt.title('id_i = {id_i} id_j = {id_j}'.format(id_i=str(self.id_i), id_j=str(self.id_j)))

    def filename(self):
        filename = self.ylabel + '_' + self.xlabel + '_idi_' + str(self.id_i) + '_idj_' + str(self.id_j) + '_step_' + str(self.step1) + '_' + str(self.step2) + '.png'
        return filename

class lammp_figure_atomi_wall(lammp_figure):
    
    def __init__(self, step1, step2, array_x, xlabel, array_y, ylabel, id_i, id_w):
        super().__init__(step1, step2, array_x, xlabel, array_y, ylabel)
        self.id_i = id_i
        self.id_j = id_w
    
    def add_title(self):
        plt.title('id_i = {id_i} wallid = {id_j}'.format(id_i=str(self.id_i), id_j=str(self.id_j)))

    def filename(self):
        filename = self.ylabel + '_' + self.xlabel + '_idi_' + str(self.id_i) + '_idwall_' + str(self.id_j) + '_step_' + str(self.step1) + '_' + str(self.step2) + '.png'
        return filename
