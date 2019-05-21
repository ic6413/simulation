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
import plotfigure.inputvariable as pi

# import calculate module
import calculate.checksimulation as cs
import createdata.inputvariable as cdi

# inputvariable
f_read_custom = pi.f_custom
f_read_thermo = pi.f_thermo
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


    def plotij(self, id_i):

        [step1_contactchange_ncjis1_or_ncwall1, step2_contactchange_ncjis1_or_ncwall1] = cs.steps_n_c_wall_less_1_n_c_j_less_1(
            f_read_custom, id_i, self.step1, self.step2
        )

        for index, step1 in enumerate(step1_contactchange_ncjis1_or_ncwall1):
            step2 = step2_contactchange_ncjis1_or_ncwall1[index]

            [n_c_atom, idlist_atom] = cs.number_contact_atom_id_collection(
                f_read_custom, id_i, step1, step2
            )

            [n_c_wall, idlist_wall] = cs.number_contact_wall_id_collection(
                f_read_custom, id_i, step1, step2
            )
            for id_jorw in (np.append(idlist_atom, idlist_wall)):

                if id_jorw > 0:
                    manystepsclass = cs.manysteps_idj(f_read_custom, id_i, id_jorw, step1, step2, method)
                elif id_jorw < 0:
                    manystepsclass = cs.manysteps_wall(f_read_custom, id_i, id_jorw, step1, step2, method)
                else:
                    sys.exit("idj = 0")

                steps = np.arange(step1, step2+1)

                # xij
                xij = manystepsclass.xij()

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
#
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

            steps = np.arange(step1, step2+1)

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

    def plotthermo(self, variable_name_list):

        df = pd.read_hdf(f_read_thermo)

        if variable_name_list == 'all':
            new_variable_name_list = list(df)
        else:
            new_variable_name_list = variable_name_list

        df_step = cs.extract_dataframe(df, self.step1, self.step2)

        steps = np.arange(self.step1, self.step2+1)

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


    def plotsingle(self, id_i, variable_name_list):
        f_read_custom_single = dp.hdf5_csv_path + cdi.put_id_on_file([id_i], 'custom_id') + '.h5'
        df = pd.read_hdf(f_read_custom_single)

        if variable_name_list == 'all':
            new_variable_name_list = list(df)
        else:
            new_variable_name_list = variable_name_list

        df_step = cs.extract_dataframe(df, self.step1, self.step2)

        steps = np.arange(self.step1, self.step2+1)
        
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

            figclass = pp.lammp_figure_atom_single(self.step1, self.step2, x_array, x_label, y_array, y_label, id_i)
            figclass.create_and_save(outputfolder_oneatom)
            plt.close('all')
        