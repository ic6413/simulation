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

# plot velocity 1
(fig, ax) = dp.plot_velocity_1_ave_y(
    log_variable_dic_list,
    51,
    np.append(np.arange(4500000, 7000000, 500000), np.array([40000000])),
    float(log_variable['dp']), "y",
    abs(float(log_variable['in_velocity'])), "V",
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
)
dp.save_close(fig, ax, di.plots_for_view_file_path(log_variable_dic_list, 'velocity_1.png'))

# plot wall average wall stress
stress_scale_width = float(log_variable['den'])*float(log_variable['g'])*float(log_variable['width_wall_dp_unit'])*float(log_variable['dp'])
stress_scale_width_str = r'$\rho_s g w $'
stress_scale_height = float(log_variable['den'])*float(log_variable['g'])*float(log_variable['z_length_create_dp_unit'])*float(log_variable['dp'])
stress_scale_height_str = r'$\rho_s g h $'
inwall_area = (
    float(log_variable['x_period_dp_unit'])
    *float(log_variable['z_length_create_dp_unit'])
    *float(log_variable['dp'])**2
)
bottom_area = (
    float(log_variable['width_wall_dp_unit'])
    *float(log_variable['x_period_dp_unit'])
    *float(log_variable['dp'])**2
)
sidewall_force_scale = inwall_area*stress_scale_height
bottomewall_force_scale = bottom_area*stress_scale_height
wallinputstepsarray = np.arange(5000000,45000000,1000000)
"""
this is for note not on paper
for (force_v_name, fixtimeave_id_name, force_scale_factor, fig_y_label, path) in [
    ('inwall_force_1', 'timeav_inwall_force', sidewall_force_scale, r'$\sigma_{12}$', 'inwall_force_1.png'),
    ('inwall_force_2', 'timeav_inwall_force', sidewall_force_scale, r'$\sigma_{22}$', 'inwall_force_2.png'),
    ('inwall_force_3', 'timeav_inwall_force', sidewall_force_scale, r'$\sigma_{32}$', 'inwall_force_3.png'),
    ('outwall_force_1', 'timeav_outwall_force', sidewall_force_scale, r'$\sigma_{12}$', 'outwall_force_1.png'),
    ('outwall_force_2', 'timeav_outwall_force', sidewall_force_scale, r'$\sigma_{22}$', 'outwall_force_2.png'),
    ('outwall_force_3', 'timeav_outwall_force', sidewall_force_scale, r'$\sigma_{32}$', 'outwall_force_3.png'),
    ('zbottom_force_1', 'timeav_zbottom_force', bottomewall_force_scale, r'$\sigma_{13}$', 'zbottom_force_1.png'),
    ('zbottom_force_2', 'timeav_zbottom_force', bottomewall_force_scale, r'$\sigma_{23}$', 'zbottom_force_2.png'),
    ('zbottom_force_3', 'timeav_zbottom_force', bottomewall_force_scale, r'$\sigma_{33}$', 'zbottom_force_3.png'),
]:
"""
for (force_v_name, fixtimeave_id_name, force_scale_factor, fig_y_label, path) in [
    ('inwall_force_1', 'timeav_inwall_force', sidewall_force_scale, 'shearing wall ' + r'$\sigma_{12}$', 'inwall_force_1.png'),
    ('inwall_force_2', 'timeav_inwall_force', sidewall_force_scale, 'shearing wall ' + r'$\sigma_{22}$', 'inwall_force_2.png'),
    ('inwall_force_3', 'timeav_inwall_force', sidewall_force_scale, 'shearing wall ' + r'$\sigma_{32}$', 'inwall_force_3.png'),
    ('outwall_force_1', 'timeav_outwall_force', sidewall_force_scale, 'static wall ' + r'$\sigma_{12}$', 'outwall_force_1.png'),
    ('outwall_force_2', 'timeav_outwall_force', sidewall_force_scale, 'static wall ' + r'$\sigma_{22}$', 'outwall_force_2.png'),
    ('outwall_force_3', 'timeav_outwall_force', sidewall_force_scale, 'static wall ' + r'$\sigma_{32}$', 'outwall_force_3.png'),
    ('zbottom_force_1', 'timeav_zbottom_force', bottomewall_force_scale, 'bottom ' + r'$\sigma_{13}$', 'zbottom_force_1.png'),
    ('zbottom_force_2', 'timeav_zbottom_force', bottomewall_force_scale, 'bottom ' + r'$\sigma_{23}$', 'zbottom_force_2.png'),
    ('zbottom_force_3', 'timeav_zbottom_force', bottomewall_force_scale, 'bottom ' + r'$\sigma_{33}$', 'zbottom_force_3.png'),
]:


    (fig, ax) = dp.plot_total_wall_force_strain(
        log_variable_dic_list,
        51,
        wallinputstepsarray,
        1, r'$\gamma$',
        force_scale_factor, fig_y_label,
        force_v_name,
        if_on_paper=False,
        if_include_0_y_axis=True,
        fixtimeave_id_name=fixtimeave_id_name,
        ifrotate_tick=True,
        ifshrink=False,
    )
    dp.save_close(fig, ax, di.plots_for_view_file_path(log_variable_dic_list, path))

# plot streamline
dp.save_plot_velocity_streamplot_contour(
    log_variable_dic_list,
    51,
    np.array([6000000, 15000000]),
    'Coord1', float(log_variable['dp']), 'y',
    'Coord2', float(log_variable['dp']), 'z',
    'velocity_2', abs(float(log_variable['in_velocity'])),
    'velocity_3', abs(float(log_variable['in_velocity'])),
    0.001,
    0.001,
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name=None,
    is_calculated_v=True,
    contour_norm = 'log',
    contour_v_min_max = "min_to_max", # or "min_to_max",
    vmin = 10**-5,
    vmax = 10**0,
    ifrotate_tick=True,
    ifshrink=False,
)
# plot fraction

# plot wall stress

# plot mu-I