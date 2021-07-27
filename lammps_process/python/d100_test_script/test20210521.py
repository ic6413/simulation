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

# create folder
os.makedirs(di.plots_for_view_folder(log_variable_dic_list), exist_ok=True)
os.makedirs(di.plots_for_paper_folder(log_variable_dic_list), exist_ok=True)

for if_on_paper in [False]:
    
    # plot fraction for average region
    dp.plot_ave_value_for_select_region_yz_strain(
        log_variable_dic_list,
        201,
        np.array([35000000]),
        1, r'$\gamma$',
        1, 'fraction',
        'fraction',
        0, 16,
        0, 4,
        if_on_paper=if_on_paper,
        ifrotate_tick=True,
        ifshrink=False,
        fixtimeave_id_fortime='avspatial_ave',
        fixtimeave_id_forcoord='avspatial_ave',
        fixtimeave_id_name=None,
        is_std=False,
        is_calculated_v=True,
    )
