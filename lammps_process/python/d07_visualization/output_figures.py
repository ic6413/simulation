#!/usr/bin/env python3
import os
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data as dd
import d00_utils.plot as dp
import numpy as np

lmp_folder_path = di.lmp_folder_path

# pickle variable
(log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)


# plot velocity 1
(fig, ax) = dp.plot_velocity_1_ave_y(
    log_variable_dic_list,
    51,
    np.arange(34000000,40000000,1000000),
    float(log_variable['dp']), "y",
    float(log_variable['in_velocity']), "V",
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
)
dp.save_close(fig, ax, di.figures_file_path('velocity_1'))
# plot streamline

# plot wall force

# plot wall stress

# plot mu-I

# plot fraction