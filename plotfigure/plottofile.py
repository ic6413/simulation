import sys
import pprint
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
# import simulation module
import osmanage as om
import datapath as dp
# import calculate module
import calculate.checksimulation as cs

# import datatofile module
import createdata.datatofile as cd

# create folder
om.create_directory(dp.debug_fig_atomij_path)
om.create_directory(dp.debug_fig_thermo_path)
om.create_directory(dp.debug_fig_oneatom_path)
# define method
method = 1

# =========format===========
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
        self.fig = plt.figure()
        plt.style.use('classic')
        #plt.autoscale(enable=True, axis='y', tight=True)
        plt.xscale('linear')
        plt.yscale('linear')
        #plt.yscale('linear')
        #plt.yscale('symlog', linthreshy=1e-20)
        

    def create_lammp_figure(self):
        

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

        
        plt.xticks(rotation=20, fontsize=14)
        plt.yticks(fontsize=14)

        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        self.add_title()

        
    def add_title(self):
        pass

    def filename(self):
        pass

    def create_and_save(self, outputfolder):
        self.create_lammp_figure()
        self.fig.savefig(outputfolder + self.filename(), bbox_inches='tight')   # save the figure to file

class lammp_figure_thermo(lammp_figure):
    
    def __init__(self, step1, step2, array_x, xlabel, array_y, ylabel):
        super().__init__(step1, step2, array_x, xlabel, array_y, ylabel)
    
    def add_title(self):
        plt.title('thermo')

    def filename(self):
        filename = self.ylabel + '_' + self.xlabel + '_step_' + str(self.step1) + '_' + str(self.step2) + '.png'
        return filename

class lammp_figure_max(lammp_figure):
    
    def __init__(self, step1, step2, array_x, xlabel, array_y, ylabel, maxlabel):
        super().__init__(step1, step2, array_x, xlabel, array_y, ylabel)
        self.maxlabel = maxlabel
    
    def add_title(self):
        plt.title('max' + self.maxlabel)

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

class lammp_3Dtrajfigure(lammp_figure):

    def __init__(self, step1, step2, array_x, array_y, array_z, xlabel,ylabel,zlabel):
        super().__init__(step1, step2, array_x, xlabel, array_y, ylabel)
        self.array_z = array_z
        self.zlabel = zlabel

    def skipsteps(self, step):
        self.array_x = self.array_x[::step]
        self.array_y = self.array_y[::step]
        self.array_z = self.array_z[::step]
        

    def create_3Dlammp_figure(self):
        
        self.skipsteps(1000)
        
        ax = self.fig.add_subplot(111, projection='3d')
        vx = np.diff(self.array_x)
        vy = np.diff(self.array_y)
        vz = np.diff(self.array_z)

        ax.quiver(self.array_x[:-1], self.array_y[:-1], self.array_z[:-1], vx, vy, vz, normalize=True, length = 0.0002)

        plt.xticks(rotation=20, fontsize=14)
        plt.yticks(fontsize=14)

        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.ylabel(self.zlabel)

        self.add_title()

        
    def add_title(self):
        plt.title('3D')

    def filename(self):
        filename = self.xlabel + '_' + self.ylabel + '_' + self.zlabel + '_step_' + str(self.step1) + '_' + str(self.step2) + '.png'
        return filename

    def create_and_save(self, outputfolder):
        self.create_3Dlammp_figure()
        self.fig.savefig(outputfolder + self.filename(), bbox_inches='tight')   # save the figure to file

# =======format end ======

# ========begin==========

class plotclass(object):

    def __init__(self):
        om.create_directory(dp.debug_print_path)
        om.create_directory(dp.debug_fig_path)
        om.create_directory(dp.debug_fig_thermo_path)
        om.create_directory(dp.debug_fig_oneatom_path)
        om.create_directory(dp.debug_fig_atomij_path)
        

class plotfromcustom(plotclass):

    def __init__(self, step1, step2, fromtraceorall):
        super().__init__()
        self.step1 = step1
        self.step2 = step2
        self.fromtraceorall = fromtraceorall
        print("plot step from {step1} to {step2}".format(step1=step1, step2=step2))
        cd.dumptofile(fromtraceorall).dump_custom()

    def plotij_singlej(self, id_i, id_jorw):

        step1 = self.step1
        step2 = self.step2

        if id_jorw > 0:
            manystepsclass = cs.manysteps_idj(dp.f_custom, id_i, id_jorw, step1, step2, method)
        elif id_jorw < 0:
            manystepsclass = cs.manysteps_wall(dp.f_custom, id_i, id_jorw, step1, step2, method)
        else:
            sys.exit("idj = 0")

        steps = np.arange(step1, step2)

        ifoverlap = manystepsclass.ifoverlap()

        # xij
        xij = manystepsclass.xij()
        x_array = steps
        x_label = 'step'
        y_array = xij
        y_label = 'xij'

        def plot_1D_xarray_1D_yarray(x_array1D, x_label1D, y_array1D, y_label1D):

            figclass = create_figclass_i_jorw(step1, step2, x_array1D, x_label1D, y_array1D, y_label1D, id_i, id_jorw)
            figclass.create_and_save(dp.debug_fig_atomij_path)
            plt.close('all')
            

        if y_array.ndim == 2 and y_array.shape[-1] == 3:

            def splitxyzarray(array):
                #transform (n,3) array to x y z component
                array1 = array[:,0]
                array2 = array[:,1]
                array3 = array[:,2]
                return [array1, array2, array3]

            def addxyzlabel(label):
                #add xyz label
                xyzlabel = ['_x', '_y', '_z']
                label1 = label + xyzlabel[0]
                label2 = label + xyzlabel[1]
                label3 = label + xyzlabel[2]
                return [label1, label2, label3]

                
            [y_label1, y_label2, y_label3] = addxyzlabel(y_label)
            [y_array1, y_array2, y_array3] = splitxyzarray(y_array)

        ## velocity
        #vijt = manystepsclass.vijt()
        #vijn = manystepsclass.vijn()
        #vijt_contactpoint = manystepsclass.vijt_contactpoint()
        #vijt_half_pre = manystepsclass.vijt_half_pre()
        ## overlap
        #overlap_length = manystepsclass.overlap_length()
        # force
        [fnk,fngamma,ftk_include_his,ftgamma] = manystepsclass.f_many_steps()
        fn = fnk + fngamma
        ft = ftk_include_his + ftgamma
        # work
        #[work_ft, work_fn] = manystepsclass.work_ftfn_many_steps()
        #[work_ftgamma, work_fngamma] = manystepsclass.work_ftgammafngamma_many_steps()
        #[work_ftk, work_fnk] = manystepsclass.work_ftkfnk_many_steps()
        #[sum_ftwork, sum_fnwork] = manystepsclass.cumsum_work()
        #sum_ftkwork = np.cumsum(work_ftk, axis=0)
        #sum_fnkwork = np.cumsum(work_fnk, axis=0)
        # gravity
        gravity = manystepsclass.gravity()
        x_array_label_list = [
            [steps,'step']
        ]
        y_array_label_list = [
            [cs.projection_scalar(fn,xij),'fprojectxij'],
            #[xij, 'xij'],
            #[vijt,'vijt'],
            #[vijn,'vijn'],
            #[vijt_half_pre,'vijt_half_pre'],
            #[vijt_contactpoint,'vijt_contactpoint'],
            #[overlap_length,'overlap_length'],
            #[fnk,'fnk'],
            #[fngamma,'fngamma'],
            #[ftk_include_his,'ftk_include_his'],
            #[ftgamma,'ftgamma'],
            #[fn,'fn'],
            #[ft,'ft'],
            #[fn+ft,'fn+ft'],
            #[work_ft,'work_ft'],
            #[work_fn,'work_fn'],
            #[work_ftgamma,'work_ftgamma'],
            #[work_fngamma,'work_fngamma'],
            #[work_ftk,'work_ftk'],
            #[work_fnk,'work_fnk'],
            #[sum_ftwork,'sum_ftwork'],
            #[sum_fnwork,'sum_fnwork'],
            #[sum_ftkwork,'sum_ftkwork'],
            #[sum_fnkwork,'sum_fnkwork'],
            #[gravity,'gravity'],
        ]
        # plot x y z component
        for variable_string in y_array_label_list:
            variable = variable_string[0]
            if variable.ndim == 2 and variable.shape[-1] == 3:
                addxyz = ['_x', '_y', '_z']
                for i in range(3):
                    y_array_label_list.append([variable[:,i], variable_string[1]+addxyz[i]])


        for xy_array_label in itertools.product(x_array_label_list, y_array_label_list):

            x_array = xy_array_label[0][0]
            x_label = xy_array_label[0][1]
            y_array = xy_array_label[1][0]
            y_label = xy_array_label[1][1]

            figclass = create_figclass_i_jorw(step1, step2, x_array, x_label, y_array, y_label, id_i, id_jorw)
            figclass.create_and_save(dp.debug_fig_atomij_path)
            plt.close('all')

        print("finish plot ij or iw")


    def plotij(self, id_i):

        step1 = self.step1
        step2 = self.step2
        [step1_contactchange_ncjis1_or_ncwall1, step2_contactchange_ncjis1_or_ncwall1] = cs.steps_n_c_wall_less_1_n_c_j_less_1(
            dp.f_custom, id_i, step1, step2
        )

        [n_c_atom, idlist_atom] = cs.number_contact_atom_id_collection(
            dp.f_custom, id_i, step1, step2
        )
        [n_c_wall, idlist_wall] = cs.number_contact_wall_id_collection(
            dp.f_custom, id_i, step1, step2
        )
        contactid_list = np.append(idlist_atom, idlist_wall)

        for id_jorw in contactid_list:
            self.plotij_singlej(id_i, id_jorw)

        print("finish plot ij")


    def plotijtest(self, id_i):

        step1 = self.step1
        step2 = self.step2

        for id_jorw in [15583]:

            if id_jorw > 0:
                manystepsclass = cs.manysteps_idj(dp.f_custom, id_i, id_jorw, step1, step2, method)
            elif id_jorw < 0:
                manystepsclass = cs.manysteps_wall(dp.f_custom, id_i, id_jorw, step1, step2, method)
            else:
                sys.exit("idj = 0")

            steps = np.arange(step1, step2)

            [work_ftk, work_fnk] = manystepsclass.work_ftkfnk_many_steps()
            sum_ftkwork = np.cumsum(work_ftk, axis=0)
            sum_fnkwork = np.cumsum(work_fnk, axis=0)
            [work_ftgamma, work_fngamma] = manystepsclass.work_ftgammafngamma_many_steps()
            sum_ftgammawork = np.cumsum(work_ftgamma, axis=0)
            sum_fngammawork = np.cumsum(work_fngamma, axis=0)
            [fnk,fngamma,ftk_include_his,ftgamma] = manystepsclass.f_many_steps()
            fn = fnk + fngamma
            ft = ftk_include_his + ftgamma
            vijt = manystepsclass.vijt()
            vijn = manystepsclass.vijn()
            x_array_label_list = [
                    [steps,'step']
                ]
            y_array_label_list = [
                [ft,'ft'],
                [fn,'fn'],
                [vijt,'vijt'],
                [vijn,'vijn'],
                [cs.length(ft)/cs.length(fn),'ftoverfn'],
                [ftgamma,'ftgamma'],
                [ftk_include_his,'ftk_include_his'],
                [work_ftk,'work_ftk'],
                [work_fnk,'work_fnk'],
                [sum_ftkwork,'sum_ftkwork'],
                [sum_fnkwork,'sum_fnkwork'],
                [sum_ftgammawork,'sum_ftgammawork'],
                [sum_fngammawork,'sum_fngammawork'],  
                ]
            # plot x y z component
            for variable_string in y_array_label_list:
                variable = variable_string[0]
                if variable.ndim == 2 and variable.shape[-1] == 3:
                    addxyz = ['_x', '_y', '_z']
                    for i in range(3):
                        y_array_label_list.append([variable[:,i], variable_string[1]+addxyz[i]])


            for xy_array_label in itertools.product(x_array_label_list, y_array_label_list):

                x_array = xy_array_label[0][0]
                x_label = xy_array_label[0][1]
                y_array = xy_array_label[1][0]
                y_label = xy_array_label[1][1]

                figclass = create_figclass_i_jorw(step1, step2, x_array, x_label, y_array, y_label, id_i, id_jorw)
                figclass.create_and_save(dp.debug_fig_atomij_path)
                plt.close('all')

        print("finish plot ij test")

        
    def plotmaxKE_everystep(self, maxlabel, variable_name_list):

        cd.dumptofile(self.fromtraceorall).dump_custom_max(maxlabel)
        debug_fig_max_path = dp.debug_fig_path + 'max_' + maxlabel + '/'
        om.create_directory(debug_fig_max_path)

        df = pd.read_hdf(dp.put_maxlabel_on_file(maxlabel, dp.f_custom))

        if variable_name_list == 'all':
            new_variable_name_list = list(df)
        else:
            new_variable_name_list = variable_name_list

        df_step = cs.extract_dataframe(df, self.step1, self.step2)

        steps = np.arange(self.step1, self.step2)

        x_array_label_list = [
            [steps,'step']
        ]

        y_array_label_list = [
            [df_step[variable_name].values, variable_name] for variable_name in new_variable_name_list
            ]

        for xy_array_label in itertools.product(x_array_label_list, y_array_label_list):

            x_array = xy_array_label[0][0]
            x_label = xy_array_label[0][1]
            y_array = xy_array_label[1][0]
            y_label = xy_array_label[1][1]

            figclass = lammp_figure_max(self.step1, self.step2, x_array, x_label, y_array, y_label, maxlabel)
            figclass.create_and_save(debug_fig_max_path)
            plt.close('all')

        print("finish plot max")


class plotfromcustomselect(plotclass):

    def __init__(self, step1, step2, id_i, fromtraceorall):
        super().__init__()
        self.step1 = step1
        self.step2 = step2
        self.id_i = id_i
        self.id_i_list = [id_i]
        self.fromtraceorall = fromtraceorall
        print("plot step from {step1} to {step2}".format(step1=step1, step2=step2))
        cd.dumptofile(fromtraceorall).dump_custom_select(self.id_i_list)
        self.f_read_custom_single = dp.put_id_on_file(self.id_i_list, dp.f_custom)

    def plotsingle(self, variable_name_list):
        df = pd.read_hdf(self.f_read_custom_single)

        if variable_name_list == 'all':
            new_variable_name_list = list(df)
        else:
            new_variable_name_list = variable_name_list

        df_step = cs.extract_dataframe(df, self.step1, self.step2)

        steps = np.arange(self.step1, self.step2)
        
        x_array_label_list = [
            [steps,'step']
        ]

        y_array_label_list = [
            [df_step[variable_name].values, variable_name] for variable_name in new_variable_name_list
            ]

        for y_array_label in y_array_label_list:
            if y_array_label[1]=='x':
                x_component = y_array_label[0]
            elif y_array_label[1]=='y':
                y_component = y_array_label[0]
            elif y_array_label[1]=='z':
                z_component = y_array_label[0]
            else:
                pass
        
        def d_zaxis_theta(x_component,y_component):
            d_zaxis = (x_component**2+y_component**2)**0.5
            theta = np.arctan(y_component/x_component)
            return [d_zaxis,theta]

        [d_zaxis,theta] = d_zaxis_theta(x_component,y_component)
        y_array_label_list.append([d_zaxis,'d_zaxis'])
        y_array_label_list.append([theta,'theta'])

        for xy_array_label in itertools.product(x_array_label_list, y_array_label_list):

            x_array = xy_array_label[0][0]
            x_label = xy_array_label[0][1]
            y_array = xy_array_label[1][0]
            y_label = xy_array_label[1][1]

            figclass = lammp_figure_atom_single(self.step1, self.step2, x_array, x_label, y_array, y_label, self.id_i)
            figclass.create_and_save(dp.debug_fig_oneatom_path)
            plt.close('all')

            
        x_array = theta
        x_label = 'theta'
        y_array = z_component
        y_label = 'z'

        figclass = lammp_figure_atom_single(self.step1, self.step2, x_array, x_label, y_array, y_label, self.id_i)
        figclass.create_and_save(dp.debug_fig_oneatom_path)
        plt.close('all')

        print("finish plot single i")


    def plot3Dtraj(self):

        df = pd.read_hdf(self.f_read_custom_single)

        df_step = cs.extract_dataframe(df, self.step1, self.step2)

        array_x = df_step['x'].values
        array_y = df_step['y'].values
        array_z = df_step['z'].values
        label_x = 'x'
        label_y = 'y'
        label_z = 'z'

        figclass = lammp_3Dtrajfigure(self.step1, self.step2, array_x, array_y, array_z, label_x,label_y,label_z)
        figclass.create_and_save(dp.debug_fig_oneatom_path)
        plt.close('all')

        print("finish plot 3Dtraj")


class plotfromthermo(plotclass):

    def __init__(self, step1, step2):
        super().__init__()
        self.step1 = step1
        self.step2 = step2
        print("plot step from {step1} to {step2}".format(step1=step1, step2=step2))
        cd.thermo_hdf5_csv()

    def plotthermo(self, variable_name_list):

        df = pd.read_hdf(dp.f_thermo)

        if variable_name_list == 'all':
            new_variable_name_list = list(df)
        else:
            new_variable_name_list = variable_name_list

        df_step = cs.extract_dataframe(df, self.step1, self.step2)

        steps = np.arange(self.step1, self.step2)

        x_array_label_list = [
            [steps,'step']
        ]

        y_array_label_list = [
            [df_step[variable_name].values, variable_name] for variable_name in new_variable_name_list
            ]

        for xy_array_label in itertools.product(x_array_label_list, y_array_label_list):

            x_array = xy_array_label[0][0]
            x_label = xy_array_label[0][1]
            y_array = xy_array_label[1][0]
            y_label = xy_array_label[1][1]

            figclass = lammp_figure_thermo(self.step1, self.step2, x_array, x_label, y_array, y_label)
            figclass.create_and_save(dp.debug_fig_thermo_path)

        plt.close('all')

        print("finish plot thermo")


class plotfromtraceprint_idi(plotclass):

    def __init__(self, id_i):
        super().__init__()
        self.id_i = id_i
        

    def traceprinttexttoarray(self, file):
        with open(file) as f:
            array = np.loadtxt(f, skiprows=1)
        return array

    def labeltoarray(self, label):
        file = dp.add_label_put_id_tracepath(label, self.id_i)
        return self.traceprinttexttoarray(file)

    def labelorstepto1Darray(self, label, step1, step2):
        if label == "step":
            array = np.arange(step1,step2)
        else:
            array = self.labeltoarray(label)
        
        if step1 < 1+dp.startstep:
            sys.exit("step1 should >= 1+startstep, step1={step1}, startstep={startstep}".format(step1=step1, startstep = dp.startstep))
        if step2 > dp.startstep+1+len(array):
            sys.exit("step2 should <= startstep+1+len(array)")

        row1 = step1 - 1 - dp.startstep
        row2 = step2 - 1 - dp.startstep

        return array[row1:row2]
        
        
    def plotsingle(self, x_label, y_label, step1, step2):

        if step1 < dp.startstep+1:
            sys.exit("step1 in trace print should be larger than startstep, step1={step1}, startstep={startstep}".format(step1=step1, startstep = dp.startstep))
        x_array = self.labelorstepto1Darray(x_label, step1, step2)
        y_array = self.labelorstepto1Darray(y_label, step1, step2)

        figclass = lammp_figure_atom_single(step1, step2, x_array, x_label, y_array, y_label, self.id_i)
        figclass.create_and_save(dp.debug_fig_oneatom_path)
        plt.close('all')

        print("finish plot single i for" + y_label + " v.s " + x_label)

    def plotsingle_multifigure(self, x_label_list, y_label_list, step1, step2):
        
        for x_label in x_label_list:
            for y_label in y_label_list:
                self.plotsingle(x_label, y_label, step1, step2)


class plotfromtraceprint_max(plotclass):

    def __init__(self, maxlabel):
        super().__init__()
        self.maxlabel = maxlabel
        self.debug_fig_max_path = dp.debug_fig_path + 'max_' + self.maxlabel + '/'
        om.create_directory(self.debug_fig_max_path)
        
        
    def traceprinttexttoarray(self, file):
        with open(file) as f:
            array = np.loadtxt(f, skiprows=1)
        return array

    def labeltoarray(self, label):
        file = dp.add_label_max_tracepath(label)
        return self.traceprinttexttoarray(file)

    def labelorstepto1Darray(self, label, step1, step2):
        if label == "step":
            array = np.arange(step1,step2)
        else:
            array = self.labeltoarray(label)
        if step1 < 1+dp.startstep:
            sys.exit("step1 should >= 1+startstep, step1={step1}, startstep={startstep}".format(step1=step1, startstep = dp.startstep))
        if step2 > dp.startstep+1+len(array):
            sys.exit("step2 should <= startstep+1+len(array)")

        row1 = step1 - 1 - dp.startstep
        row2 = step2 - 1 - dp.startstep

        return array[row1:row2]
        
        
    def plotsingle(self, x_label, y_label, step1, step2):

        if step1 < dp.startstep+1:
            sys.exit("step1 in trace print should be larger than startstep, step1={step1}, startstep={startstep}".format(step1=step1, startstep = dp.startstep))

        x_array = self.labelorstepto1Darray(x_label, step1, step2)
        y_array = self.labelorstepto1Darray(y_label, step1, step2)
        figclass = lammp_figure_max(step1, step2, x_array, x_label, y_array, y_label, self.maxlabel)
        figclass.create_and_save(self.debug_fig_max_path)
        plt.close('all')

        print("finish plot single i for" + y_label + " v.s " + x_label)

    def plotsingle_multifigure(self, x_label_list, y_label_list, step1, step2):
        
        for x_label in x_label_list:
            for y_label in y_label_list:
                self.plotsingle(x_label, y_label, step1, step2)