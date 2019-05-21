# import
import os.path
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
# import module in simulation folder
import osmanage as om
import datapath as dp
# import module in createdata folder
import createdata.dataarrange as da
import createdata.inputvariable as cdi
# ======================================

# inputvariable
post_process_path = dp.post_process_path
hdf5_csv_path = dp.hdf5_csv_path
thermo_path = dp.thermo_path
custom_path = dp.custom_path
f_thermo = cdi.f_thermo
f_custom = cdi.f_custom
f_cipcj = cdi.f_cipcj
# current module inputvariable
custom_id_i_list = cdi.custom_id_i_list
combine_id_i_list = cdi.combine_id_i_list
pair_path = dp.pair_path
# end inputvariable

om.create_directory(post_process_path)
om.create_directory(hdf5_csv_path)

def thermo_hdf5_csv():
    da.thermo_file_to_h5_csv(thermo_path, f_thermo, override ='no')

def dump_custom():
    da.file_to_h5_csv(None, custom_path, custom_id_i_list, f_custom, override ='no')

def dump_custom_pair_combine():
    da.file_to_h5_csv(pair_path, custom_path, combine_id_i_list, f_cipcj, override ='no')

def dump_custom_select(id_i_list):
    df = pd.read_hdf(f_custom + '.h5', 'df')
    df = da.select_custom(df, id_i_list)
    f_output_name = dp.hdf5_csv_path + cdi.put_id_on_file(id_i_list, 'custom_id')
    df.to_hdf(f_output_name + '.h5', key='df', mode='w')
    df.to_csv(f_output_name + '.csv', encoding='utf-8')
    
    