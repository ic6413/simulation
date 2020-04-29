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
# import project module
import dataarrange as da
import datapath as dp
import output_control as oc
import os
# ======================================

custom_id_i_list = 'all'

# create folder
os.makedirs(dp.post_process_path, exist_ok=True)
os.makedirs(dp.hdf5_csv_path, exist_ok=True)
da.file_to_h5_csv(None, dp.custom_near_trace_path, custom_id_i_list, dp.f_custom, override ='yes')
