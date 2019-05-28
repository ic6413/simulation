import sys
import pprint
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
# import simulation module
import osmanage as om
import datapath as dp

# import plot module
import plotfigure.plotlammpmodule as pp

# import calculate module
import calculate.checksimulation as cs

# inputvariable
f_read_custom = dp.f_custom
f_read_thermo = dp.f_thermo
outputfolder_atomij = dp.debug_fig_atomij_path
outputfolder_thermo = dp.debug_fig_thermo_path
outputfolder_oneatom = dp.debug_fig_oneatom_path
om.create_directory(outputfolder_atomij)
om.create_directory(outputfolder_thermo)
om.create_directory(outputfolder_oneatom)
method = 1

# ========begin==========

class plotclass(object):

    def __init__(self, step1, step2):
        self.step1 = step1
        self.step2 = step2
        print("plot step from {step1} to {step2}".format(step1=step1, step2=step2))

    def plotsingle(self, id_i, variable_name_list):
        f_read_custom_single = dp.put_id_on_file([id_i], dp.f_custom)
        df = pd.read_hdf(f_read_custom_single)

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
        
        def radius_theta(x_component,y_component):
            radius = (x_component**2+y_component**2)**0.5
            theta = np.arctan(y_component/x_component)
            return [radius,theta]

        [radius,theta] = radius_theta(x_component,y_component)
        y_array_label_list.append([radius,'radius'])
        y_array_label_list.append([theta,'theta'])

        for xy_array_label in itertools.product(x_array_label_list, y_array_label_list):

            x_array = xy_array_label[0][0]
            x_label = xy_array_label[0][1]
            y_array = xy_array_label[1][0]
            y_label = xy_array_label[1][1]

            figclass = pp.lammp_figure_atom_single(self.step1, self.step2, x_array, x_label, y_array, y_label, id_i)
            figclass.create_and_save(outputfolder_oneatom)
            plt.close('all')

            
        x_array = theta
        x_label = 'theta'
        y_array = z_component
        y_label = 'z'

        figclass = pp.lammp_figure_atom_single(self.step1, self.step2, x_array, x_label, y_array, y_label, id_i)
        figclass.create_and_save(outputfolder_oneatom)
        plt.close('all')

        print("finish plot single i")


    def plotij_singlej(self, id_i, id_jorw):

        step1 = self.step1
        step2 = self.step2

        if id_jorw > 0:
            manystepsclass = cs.manysteps_idj(f_read_custom, id_i, id_jorw, step1, step2, method)
        elif id_jorw < 0:
            manystepsclass = cs.manysteps_wall(f_read_custom, id_i, id_jorw, step1, step2, method)
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

            figclass = pp.create_figclass_i_jorw(step1, step2, x_array1D, x_label1D, y_array1D, y_label1D, id_i, id_jorw)
            figclass.create_and_save(outputfolder_atomij)
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

            figclass = pp.create_figclass_i_jorw(step1, step2, x_array, x_label, y_array, y_label, id_i, id_jorw)
            figclass.create_and_save(outputfolder_atomij)
            plt.close('all')

        print("finish plot ij or iw")


    def plotij(self, id_i):

        step1 = self.step1
        step2 = self.step2
        [step1_contactchange_ncjis1_or_ncwall1, step2_contactchange_ncjis1_or_ncwall1] = cs.steps_n_c_wall_less_1_n_c_j_less_1(
            f_read_custom, id_i, step1, step2
        )

        [n_c_atom, idlist_atom] = cs.number_contact_atom_id_collection(
            f_read_custom, id_i, step1, step2
        )
        [n_c_wall, idlist_wall] = cs.number_contact_wall_id_collection(
            f_read_custom, id_i, step1, step2
        )
        contactid_list = np.append(idlist_atom, idlist_wall)

        for id_jorw in contactid_list:
            self.plotij_singlej(id_i, id_jorw)

        print("finish plot ij")


    def plot3Dtraj(self, id_i):
        try:
            f_read_custom_single = dp.put_id_on_file([id_i], dp.f_custom)
        except FileNotFoundError:
            import createdata.datatofile as cd
            cd.dump_custom_select([id_i])
        df = pd.read_hdf(f_read_custom_single)

        df_step = cs.extract_dataframe(df, self.step1, self.step2)

        array_x = df_step['x'].values
        array_y = df_step['y'].values
        array_z = df_step['z'].values
        label_x = 'x'
        label_y = 'y'
        label_z = 'z'

        figclass = pp.lammp_3Dtrajfigure(self.step1, self.step2, array_x, array_y, array_z, label_x,label_y,label_z)
        figclass.create_and_save(outputfolder_oneatom)
        plt.close('all')

        print("finish plot 3Dtraj")


    def plotijtest(self, id_i):

        step1 = self.step1
        step2 = self.step2

        for id_jorw in [15583]:

            if id_jorw > 0:
                manystepsclass = cs.manysteps_idj(f_read_custom, id_i, id_jorw, step1, step2, method)
            elif id_jorw < 0:
                manystepsclass = cs.manysteps_wall(f_read_custom, id_i, id_jorw, step1, step2, method)
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

                figclass = pp.create_figclass_i_jorw(step1, step2, x_array, x_label, y_array, y_label, id_i, id_jorw)
                figclass.create_and_save(outputfolder_atomij)
                plt.close('all')

        print("finish plot ij test")

    def plotthermo(self, variable_name_list):

        df = pd.read_hdf(f_read_thermo)

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

            figclass = pp.lammp_figure_thermo(self.step1, self.step2, x_array, x_label, y_array, y_label)
            figclass.create_and_save(outputfolder_thermo)

        plt.close('all')

        print("finish plot thermo")
