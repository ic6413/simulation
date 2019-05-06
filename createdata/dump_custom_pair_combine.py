# import
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
import createdata.data_arrange as da
import createdata.inputvariable as cdi
# ======================================

# inputvariable begin
f_cipcj = cdi.f_cipcj
combine_id_i_list = cdi.combine_id_i_list
custom_path = dp.custom_path
pair_path = dp.pair_path
hdf5_csv_path = dp.hdf5_csv_path
post_process_path = dp.post_process_path
# inputvariable end

# create folder
om.create_directory(post_process_path)
om.create_directory(hdf5_csv_path)
da.file_to_h5_csv(pair_path, custom_path, combine_id_i_list, f_cipcj, override ='yes')
