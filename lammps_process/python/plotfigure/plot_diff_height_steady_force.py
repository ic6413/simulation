# import
import os
from io import StringIO
import re
import time
from itertools import chain
from itertools import repeat
from itertools import islice
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import pickle
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D
# import module
import datapath as dp
import read_setting.read_setting as rr
# import calculate setting
import read_setting.calculate_setting as rc
n_lastpoints = 250

# plot style

plt.style.use('classic')


# define function for extract data from fix txt to dataframe
if rr.logfile["shearwall"] == "zcylinder":
    wallfiles = [
                "force_zbottom_to_particle.allstep",
    ]
if rr.logfile["shearwall"] == "yplane":
    wallfiles = [
                "force_zbottom_to_particle.allstep",
    ]
    
    

diameter = float(rr.logfile['dp'])
width_wall_dp_unit = int(rr.logfile['width_wall_dp_unit'])

g = float(rr.logfile['g'])
d_step = int(rr.logfile['freq_ave_wall'])

lmp_path_list = [
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_240_Sa_2e-6_his_yes_xmu_5e-1_noshear_2e7/f_2e7_run_5e7/f_55e6_run_5e7/",
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_120_Sa_2e-6_his_yes_xmu_5e-1_noshear_2e7/f_2e7_run_5e7/",
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_60_Sa_2e-6_his_yes_xmu_5e-1_noshear_until_echeck_1e-17/f_5e6_run_1e7/f_15e6_run_15e7_not_complete/",
    "/home/ic6413/lmp_run/block_xp_50_w_16_h_30_Sa_2e-6_his_yes_xmu_5e-1_noshearfor_5e6_and_run_1e7/f_5e6_run_1e7/f_15e6_run_15e7/",
]

height_list = [
    240,
    120,
    60,
    30,
]
def plot_wall_diff_height_steady_time(if_plot_to_last, step1, step2, figformat="png", ifpickle=False, ifplotfrominitial=True):
    fig_handle = plt.figure()
    total_time_array = np.empty(len(height_list))
    for k, lammps_directory in enumerate(lmp_path_list):
        
        for wallfile in wallfiles:
            if ifplotfrominitial:
                logfilelist_from_initial_to_lastest = rr.log_current_plus_previousfrom_initial_to_lastest(lammps_directory)
                for index in range(len(logfilelist_from_initial_to_lastest)):
                    n_log_list = len(logfilelist_from_initial_to_lastest)
                    folder_path_list_initial_to_last = [lammps_directory+"../"*(n_log_list-1-n) for n in range(n_log_list)]
                    folder_path_list_last_to_initial = [folder_path_list_initial_to_last[n_log_list-1-n] for n in range(n_log_list)]

                    with open(folder_path_list_last_to_initial[index] + "output/wall/" + wallfile) as f:
                        lines = f.read().strip().split('\n')
                    header = lines[1].split()[1:]
                    
                    step1_default = int(lines[2].split()[0])
                    ## select data
                    data = [lines[t].split() for t in range(2, len(lines))]
                    ## attach data
                    df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
                    if index >= 1:
                        df = df.loc[df['TimeStep']<next_step1]
                    array = df['v_t'].values
                    calculate_setting_diclist_from_initial_to_last=[]
                    log_path_list_last_to_initial = [lammps_directory+"../"*n+'log.lammps' for n in range(n_log_list)]

                    log_path_list_initial_to_last = [log_path_list_last_to_initial[n_log_list-1-n] for n in range(n_log_list)]

                    def count_restart_time(index):
                        restart_time = 0
                        for n in range(index+1):
                            if n == 0:
                                restart_time_from_last_to_current = 0
                            else:
                                logfile = logfilelist_from_initial_to_lastest[n]
                                rst_from_current = str(logfile['rst_from'])
                                last_log_path = log_path_list_initial_to_last[n-1]
                            
                                if os.path.isfile(last_log_path):
                                    with open(last_log_path, mode='r') as f:
                                        lines = f.read().strip().split('\n')
                                else:
                                    sys.exit("file not exist")

                                for n_line1, line in enumerate(lines):
                                    if len(line) >= 2:
                                        if line.split()[0]=="Step" and line.split()[1]=="Time":
                                            break
                                
                                for n_line2, line in enumerate(lines):
                                    if len(line) >= 2:
                                        if line.split()[0]=="Loop" and line.split()[1]=="time":           
                                            break

                                for n_line, line in enumerate(lines):
                                    if n_line>n_line1 and n_line<n_line2:
                                        if line.split()[0]==rst_from_current:
                                            restart_time_from_last_to_current = float(line.split()[1])                        
                                            break
                            restart_time += restart_time_from_last_to_current
                        return restart_time
                    def calculate_setting_dic(index):
                        dic = {}
                        dic["previous_time"] = count_restart_time(index)
                        return dic                    
                    for n_log in range(n_log_list):
                        calculate_setting_diclist_from_initial_to_last.append(calculate_setting_dic(n_log))
                    array += calculate_setting_diclist_from_initial_to_last[n_log_list-1-index]["previous_time"]
                    df['v_t'] = array
                    next_step1 = step1_default
                    
                    if index == 0:
                        df_out = df
                    else:
                        df_out = pd.concat([df,df_out])

                df_full = df_out
            else:
                with open(dp.lammps_directory + "output/wall/" + wallfile) as f:
                    lines = f.read().strip().split('\n')
                header = lines[1].split()[1:]
                
                data = [lines[t].split() for t in range(2, len(lines))]
                df_full = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            
            step1_default = int(df_full['TimeStep'].min())
            step2_default = int(df_full['TimeStep'].max())
        
        def total_time(step1_1, step2_1):
            df = df_full[(df_full['TimeStep']<=step2_1) & (df_full['TimeStep']>=step1_1)]
            
            ## repeat timestep

            variable1 = 'v_t'
            
            for variable2 in ["v_force_zbottom_x"]:
                
                x_array = df[variable1].values
                y_array = df[variable2].values
                ave_yarray = np.sum(y_array[-n_lastpoints-1:-1])/n_lastpoints
                if ave_yarray > 0:
                    if_larger90 = (y_array>(ave_yarray*0.9))
                    first_index_larger90 = [n for n in range(len(if_larger90)) if if_larger90[n] == True][0]
                    if_larger10 = (y_array>(ave_yarray*0.1))
                    first_index_larger10 = [n for n in range(len(if_larger10)) if if_larger10[n] == True][0]
                elif ave_yarray < 0:
                    if_larger90 = (y_array<(ave_yarray*0.9))
                    first_index_larger90 = [n for n in range(len(if_larger90)) if if_larger90[n] == True][0]
                    if_larger10 = (y_array<(ave_yarray*0.1))
                    first_index_larger10 = [n for n in range(len(if_larger10)) if if_larger10[n] == True][0]
                else:
                    sys.exit("average equal 0")
                totaltime = x_array[first_index_larger90]-x_array[first_index_larger10]
            return totaltime

        if if_plot_to_last:
            total_time_array[k] = total_time(step1_default, step2_default)
        else:
            total_time_array[k] = total_time(step1, step2)
    plt.plot(np.array(height_list), total_time_array) 
    plt.xlabel("height")
    plt.ylabel("time_from10to90")
    plt.tight_layout()
    fig_handle.savefig(dp.lammps_directory + "walls" + "." + figformat, format=figformat)
    if ifpickle:
        # Save figure handle to disk
        with open(dp.lammps_directory + "walls" + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
            pickle.dump(fig_handle, f)

    plt.close('all')