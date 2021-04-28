# import
import os
from io import StringIO
import re
from itertools import chain
from itertools import repeat
from itertools import islice
import pandas as pd
import numpy as np
# import module
import datapath as dp
import read_setting as rr

def combine_previous_wall_data(wallfile_name, ifplotfrominitial=False, ifplotfromrotate=True):
    for index in range(rr.n_simu_total):

        if ifplotfromrotate:
            dic = rr.log_variable_dic_list[::-1][index]
            if dic["ifrotate"] == "no" or float(dic["Sa"]) == 0:
                break

        with open(rr.folder_path_list_last_to_initial[index] + "output/wall/" + wallfile_name) as f:
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
        array += rr.log_variable_dic_list[rr.n_simu_total-1-index]["previous_time"]
        df['v_t'] = array
        next_step1 = step1_default
        
        if index == 0:
            df_out = df
        else:
            df_out = pd.concat([df,df_out])

    df_full = df_out

    return df_full


def combine_previous_wall_data_and_save(wallfile_name, ifplotfrominitial=False, ifplotfromrotate=True):
    combine_previous_wall_data(wallfile_name, ifplotfrominitial, ifplotfromrotate).to_hdf(dp.combine_previous_from_output + "wall/" + wallfile_name, key='df', mode='w', format='fixed')

def combine_previous_allwall_data_and_save():
    if rr.log_variable["shearwall"] == "zcylinder":
        wallfiles_name = [
                    "force_zbottom_to_particle.allstep",
                    "force_outwall_to_particle.allstep",
                    "force_inwall_to_particle.allstep",
        ]
    if rr.log_variable["shearwall"] == "yplane":
        wallfiles_name = [
                    "force_zbottom_to_particle.allstep",
                    "force_y_top_to_particle.allstep",
                    "force_y_bottom_to_particle.allstep",
        ]
    for wallfile_name in wallfiles_name:
        combine_previous_wall_data_and_save(wallfile_name, ifplotfrominitial=False, ifplotfromrotate=True)


def combine_previous_chunk_data():
    for index in range(rr.n_simu_total):
        log_variable_inthis_index = rr.log_variable_dic_list[index]
        d_step = int(log_variable_inthis_index['freq_ave_chunk_momentum_mass_field'])
        with open(rr.folder_path_list_last_to_initial[index] + "output/momentum_mass_field/" + "fix.momentum_mass_field.all") as f:
            lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        n_line_in_a_step = int(lines[3].split()[1])
        step1_default = int(lines[3].split()[0])
        step2_default = int(lines[-1 - n_line_in_a_step].split()[0])-d_step
        def index_of_variable_in_header(variable_name_in_header):
            return [n for n in range(0,len(header)) if header[n]==variable_name_in_header][0]
        def data_inloop(step_smallloop):
            n_line_0 = (step_smallloop - step1_default)/d_step*(n_line_in_a_step+1) + 4
            n_line_1 = n_line_0 + n_line_in_a_step
            ## select data
            data = [lines[t].split() for t in range(int(n_line_0), int(n_line_1))]
            data = np.array(data, dtype=np.float64)
            return data
        step_array = np.arange(step1_default, step2_default + d_step, d_step)
        time_array = step_array*log_variable_inthis_index["ts"] + rr.log_variable_dic_list[rr.n_simu_total-1-index]["previous_time"]
        data_array = np.empty([step_array.shape[0], n_line_in_a_step, len(header)])
        for n1 in range(step_array.shape[0]):
            data_array[n1,:,:] = data_inloop(step_array[n1])
        
        ## select data
        data = [lines[t].split() for t in range(2, len(lines))]
        ## attach data
        if index == 0:
            df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        if index >= 1:
            df = df.loc[df['TimeStep']<next_step1]
        array = df['v_t'].values
        array += rr.log_variable_dic_list[rr.n_simu_total-1-index]["previous_time"]
        df['v_t'] = array
        next_step1 = step1_default
        
        if index == 0:
            df_out = df
        else:
            df_out = pd.concat([df,df_out])

    df_full = df_out

    return df_full