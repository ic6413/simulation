# import
import os
import re
import time
from itertools import chain
from itertools import repeat
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D
import rename_variable as rv
# =======================================

# Change sign when switch i j
header_change_sign = ['fjinx', 'fjiny', 'fjinz', 'fjitx', 'fjity', 'fjitz', 'vijnx', 'vijny', 'vijnz', 'vijtx', 'vijty', 'vijtz']
# Switch i j for id and type
header_change_from = ['id_i', 'id_j', 'type_i', 'type_j']
header_change_to = ['id_j', 'id_i', 'type_j', 'type_i']

# ====================================== thermo =========================

def thermofile_to_dataframe(file):
    with open(file) as f:
        lines = f.read().strip().split('\n')
    # index of first row of every block
    id_headers = [n for n, line in enumerate(lines) if line.startswith('Step Atoms')]
    if len(id_headers) > 1:
        string = 'number of header = {n}'.format(n=len(id_headers))
        sys.exit(string)
    id_header = id_headers[0]
    # index of data end
    id_ends = [n for n, line in enumerate(lines) if line.startswith('Loop time of')]
    id_end = id_ends[0]
    # header
    header = lines[id_header].split()[0:]
    # select data
    id_data_begin = id_header + 1
    data = [lines[t].split() for t in range(id_data_begin, id_end)]
    # attach data
    df = pd.DataFrame(data = data, columns = header, dtype = 'float64')

    return df


def thermo_file_to_h5_csv(f_thermo, f_output_name, override ='no'):

    if os.path.isfile(f_output_name + '.h5') and (override != 'yes'):
        sys.exit('h5 file exist')
    else:
        df = thermofile_to_dataframe(f_thermo)
        df.to_hdf(f_output_name + '.h5', key='df', mode='w')
        df.to_csv(f_output_name + '.csv', encoding='utf-8')


# ====================================== pari custom =====================

# define function for extract data from dump txt to dataframe
def dumpfile_to_dataframe(file):
	# read data as lines
	with open(file) as f:
		lines = f.read().strip().split('\n')
	# index of first row of every block	
	id_ITEMTIMESTEPs = [n for n, line in enumerate(lines) if line.startswith('ITEM: TIMESTEP')]
	# header
	header = lines[id_ITEMTIMESTEPs[0] + 8].split()[2:]
	# index of row timestep, start, end
	id_ITEMTIMESTEPs.append(len(lines))
	# combine ranges of blocks
	iter = chain.from_iterable(range(id+9, id_ITEMTIMESTEPs[i + 1]) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1]))
	# select data
	data = [lines[t].split() for t in iter]
	# attach data
	df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
	# repeat timestep
	steps = list(chain.from_iterable(repeat(lines[id + 1], id_ITEMTIMESTEPs[i + 1] - (id + 9)) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1])))
	steps = np.asarray(steps,dtype=np.float64)
	# insert timesteps
	df.insert(1, 'step', steps)
	return df

# extract custom for specific id
def select_custom(df_custom, select_id_list):
    if select_id_list == 'all':
        return df_custom
    else:
        # select traceid
        df_select_custom = df_custom.loc[df_custom['id'].isin(select_id_list)]
        # check id appear in every steps
        steps_array = df_select_custom['step'].values
        check_first_step_is_zero = (steps_array[0: len(select_id_list)] == 0).all()
        check_other_steps = (steps_array[0: -len(select_id_list)] + 1 == steps_array[len(select_id_list): ]).all()
        check_all = (check_other_steps and check_first_step_is_zero)
        if not check_all:
            firstwrong = np.where(np.logical_not((steps_array[0: -len(select_id_list)] + 1 == steps_array[len(select_id_list): ])))[0][0]
            print('one of select atom not in step' + repr(steps_array[firstwrong]))
        return df_select_custom

# create ji pair from ij pair and combine together
def pair_double(df_pair, header_change_from, header_change_to, header_change_sign):
    # copy
    df_switch = df_pair.copy()
    # change sign
    df_switch[header_change_sign] = -df_switch[header_change_sign].values
    # change header names
    mapper = {v: header_change_to[i] for i, v in enumerate(header_change_from)}
    df_switch = df_switch.rename(columns=mapper)
    # reorder columns
    df_switch = df_switch[list(df_pair)]
    # combine df_pair and df_switch
    dfp_data_double = np.empty((2*df_pair.shape[0], df_pair.shape[1]))
    dfp_data_double[::2] = df_pair.values
    dfp_data_double[1::2] = df_switch.values
    df_double = pd.DataFrame(data = dfp_data_double, columns = list(df_pair))
    return df_double

# merge pairdouble and select_custom, output dfcip
def merge_ci_p(dfc, dfp2):
    # add suffix on custom
    dfc = dfc.add_suffix('_i')
    dfc = dfc.rename(columns={'step_i': 'step'})
    dfcip = pd.merge(dfc, dfp2, how='left', on=['step', 'id_i'])
    # check if type_i_x = type_i_y except nan (means no pair for the timestep). If yes then delete type_i_y
    type_notconsistent = np.where(dfcip['type_i_x'].values != dfcip['type_i_y'].values)[0]
    where_type_i_y_nan = np.where(np.isnan(dfcip['type_i_y'].values))[0]
    check_only_nan_inconsistent = (type_notconsistent == where_type_i_y_nan).all()
    if ~check_only_nan_inconsistent:
        sys.exit('type_i not consistent when merge custom and pair')
    else:
        # delete rename type_i_x, delete type_i_y
        dfcip = dfcip.drop('type_i_y', axis=1)
        dfcip = dfcip.rename(columns={'type_i_x': 'type_i'})
    return dfcip

# merge dfcip and custom_j, output dfcipcj
def merge_cip_cj(dfcip, dfc):
    # add suffix on custom
    dfc = dfc.add_suffix('_j')
    dfc = dfc.rename(columns={'step_j': 'step'})
    dfcipcj = pd.merge(dfcip, dfc, how='left', on=['step', 'id_j'])
    # check if type_j_x = type_j_y except nan (means no pair for the timestep). If yes then delete type_j_y
    type_notconsistent = np.where(dfcipcj['type_j_x'].values != dfcipcj['type_j_y'].values)[0]
    where_type_j_y_nan = np.where(np.isnan(dfcipcj['type_j_y'].values))[0]
    check_only_nan_inconsistent = (type_notconsistent == where_type_j_y_nan).all()
    if ~check_only_nan_inconsistent:
        sys.exit('type_j not consistent when merge custom and pair')
    else:
        # delete rename type_j_x, delete type_j_y
        dfcipcj = dfcipcj.drop('type_j_y', axis=1)
        dfcipcj = dfcipcj.rename(columns={'type_j_x': 'type_j'})
    return dfcipcj

# merge pairdouble, custom_i and custom_j. output dfcipcj 
def merge_ci_p_cj(dfci_select, dfp2, dfc):
    dfcip = merge_ci_p(dfci_select, dfp2)
    dfcipcj = merge_cip_cj(dfcip, dfc)
    return dfcipcj

def file_to_h5_csv(f_dumppair, f_dumpcustom, id_i_list, f_output_name, override ='no'):
    
    if f_dumppair != None:
        # if file exist and override not yes, then not to write files
        if os.path.isfile(f_output_name + '.h5') and (override != 'yes'):
            sys.exit('h5 file exist')
        else:
            dfp = dumpfile_to_dataframe(f_dumppair)
            dfp = dfp.rename(columns=rv.rename_mapper_pair)
            dfc = dumpfile_to_dataframe(f_dumpcustom)
            dfc_select = select_custom(dfc, id_i_list)
            dfp2 = pair_double(dfp, header_change_from, header_change_to, header_change_sign)
            df_select_cipcj = merge_ci_p_cj(dfc_select, dfp2, dfc)
            df_select_cipcj.to_hdf(f_output_name + '.h5', key='df', mode='w')
            df_select_cipcj.to_csv(f_output_name + '.csv', encoding='utf-8')

    else:
        if os.path.isfile(f_output_name + '.h5') and (override != 'yes'):
            sys.exit('h5 file exist')
        else:
            dfc = dumpfile_to_dataframe(f_dumpcustom)
            dfc_select = select_custom(dfc, id_i_list)
            dfc_select.to_hdf(f_output_name + '.h5', key='df', mode='w')
            dfc_select.to_csv(f_output_name + '.csv', encoding='utf-8')


def create_directory(post_process_folder_path):
    # define the name of the directory to be created
    if not os.path.isdir(post_process_folder_path):
        try:  
            os.mkdir(post_process_folder_path)
        except OSError:  
            print ("Creation of the directory %s failed" % post_process_folder_path)
        else:  
            print ("Successfully created the directory %s " % post_process_folder_path)