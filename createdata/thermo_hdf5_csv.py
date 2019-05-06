# thermo hdf5
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
# ======================================

# create folder
om.create_directory(dp.post_process_path)
om.create_directory(dp.hdf5_csv_path)
da.thermo_file_to_h5_csv(dp.thermo_path, dp.f_thermo, override ='yes')