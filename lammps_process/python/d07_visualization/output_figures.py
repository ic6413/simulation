#!/usr/bin/env python3
import os
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data as dd
import d00_utils.plot as dp
import d00_utils.data_for_plot as ddfp
import numpy as np

lmp_folder_path = di.lmp_folder_path

# pickle variable
(log_variable_dic_list, n_simu_total, log_variable, folder_path_list_initial_to_last) = dr.dump_variable(lmp_folder_path)

# create folder
os.makedirs(di.plots_for_view_folder(log_variable_dic_list), exist_ok=True)
os.makedirs(di.plots_for_paper_folder(log_variable_dic_list), exist_ok=True)
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
for if_on_paper in [False, True]:
    """
    # plot fraction for strain = 0
    dp.save_plot_fraction_contour(
        log_variable_dic_list,
        1,
        np.array([5000000]),
        'Coord1', float(log_variable['dp']), 'y',
        'Coord2', float(log_variable['dp']), 'z',
        'fraction', 1,
        if_on_paper=if_on_paper,
        fixtimeave_id_fortime='avspatial_ave',
        fixtimeave_id_forcoord='avspatial_ave',
        fixtimeave_id_name=None,
        is_calculated_v=True,
        contour_norm = 'linear',
        contour_v_min_max = "constant", # or "min_to_max",
        vmin = 10**-5,
        vmax = 10**0,
        ifrotate_tick=True,
        ifshrink=False,
        add_pre_filename_string='f0001',
    )

    # plot fraction
    dp.save_plot_fraction_contour(
        log_variable_dic_list,
        201,
        np.array([6000000, 15000000, 35000000]),
        'Coord1', float(log_variable['dp']), 'y',
        'Coord2', float(log_variable['dp']), 'z',
        'fraction', 1,
        if_on_paper=if_on_paper,
        fixtimeave_id_fortime='avspatial_ave',
        fixtimeave_id_forcoord='avspatial_ave',
        fixtimeave_id_name=None,
        is_calculated_v=True,
        contour_norm = 'linear',
        contour_v_min_max = "constant", # or "min_to_max",
        vmin = 10**-5,
        vmax = 10**0,
        ifrotate_tick=True,
        ifshrink=False,
        add_pre_filename_string='f0002',
    )

    # plot fraction and for average region
    add_pre_filename_string = 'f0003'
    for (n_y_0, d_n_y, n_z_0, d_n_z, filename) in [
        (0, 16, 0, 4, add_pre_filename_string + 'fraction_region.png'),
        (0, 8, 0, 4, add_pre_filename_string + 'fraction_region_near_shearing_1.png'),
        (0, 2, 0, 4, add_pre_filename_string + 'fraction_region_near_shearing_2.png'),
        (8, 8, 0, 4, add_pre_filename_string + 'fraction_region_static_1.png'),
        (14, 2, 0, 4, add_pre_filename_string + 'fraction_region_static_2.png'),
    ]:
        dp.save_plot_ave_value_for_select_region_yz_strain(
            log_variable_dic_list,
            201,
            np.arange(4000000,45000000,1000000),
            1, r'$\gamma$',
            1, 'fraction',
            'fraction',
            n_y_0, d_n_y,
            n_z_0, d_n_z,
            if_on_paper=if_on_paper,
            ifrotate_tick=True,
            ifshrink=False,
            fixtimeave_id_fortime='avspatial_ave',
            fixtimeave_id_forcoord='avspatial_ave',
            fixtimeave_id_name=None,
            is_std=False,
            is_calculated_v=True,
            filename=filename,
        )

    # plot contact number for average region
    add_pre_filename_string = 'f0004'
    for (n_y_0, d_n_y, n_z_0, d_n_z, filename) in [
        (0, 16, 0, 4, add_pre_filename_string + 'contact_region.png'),
        (0, 8, 0, 4, add_pre_filename_string + 'contact_region_near_shearing_1.png'),
        (0, 2, 0, 4, add_pre_filename_string + 'contact_region_near_shearing_2.png'),
        (8, 8, 0, 4, add_pre_filename_string + 'contact_region_static_1.png'),
        (14, 2, 0, 4, add_pre_filename_string + 'contact_region_static_2.png'),
    ]:
        dp.save_plot_ave_ratio_value_for_select_region_yz_strain(
            log_variable_dic_list,
            201,
            np.arange(4000000,45000000,1000000),
            1, r'$\gamma$',
            1, 'contact number',
            'n_contact', 'Ncount',
            n_y_0, d_n_y,
            n_z_0, d_n_z,
            if_on_paper=if_on_paper,
            ifrotate_tick=True,
            ifshrink=False,
            fixtimeave_id_fortime='avspatial_ave',
            fixtimeave_id_forcoord='avspatial_ave',
            fixtimeave_id_name1='contact_ave',
            fixtimeave_id_name2='avspatial_ave',
            is_calculated_v1=False,
            is_calculated_v2=False,
            filename=filename,
        )

    # plot streamline
    dp.save_plot_velocity_streamplot_contour(
        log_variable_dic_list,
        201,
        np.array([6000000, 15000000, 35000000]),
        'Coord1', float(log_variable['dp']), 'y',
        'Coord2', float(log_variable['dp']), 'z',
        'velocity_2', abs(float(log_variable['in_velocity'])),
        'velocity_3', abs(float(log_variable['in_velocity'])),
        0.001,
        0.001,
        if_on_paper=if_on_paper,
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
        add_pre_filename_string='f0005',
    )

    # plot velocity 1, number of average is 1
    dp.save_plot_velocity_1_ave_y(
        log_variable_dic_list,
        1,
        np.append(np.arange(4990000, 5050000, 10000), np.array([15000000, 35000000])),
        float(log_variable['dp']), "y",
        abs(float(log_variable['in_velocity'])), "Vx",
        if_on_paper=if_on_paper,
        if_include_0_y_axis=True,
        n_sum_over_axis = 2,
        fixtimeave_id_name='avspatial_ave',
        ifrotate_tick=True,
        ifshrink=False,
        filename='f0006' + 'velocity_1.png',
    )

    # plot velocity 1
    dp.save_plot_velocity_1_ave_y(
        log_variable_dic_list,
        51,
        np.append(np.arange(4500000, 7000000, 500000), np.array([15000000, 35000000])),
        float(log_variable['dp']), "y",
        abs(float(log_variable['in_velocity'])), "Vx",
        if_on_paper=if_on_paper,
        if_include_0_y_axis=True,
        n_sum_over_axis = 2,
        fixtimeave_id_name='avspatial_ave',
        ifrotate_tick=True,
        ifshrink=False,
        filename='f0007' + 'velocity_1.png',
    )

    velocityinputstepsarray = np.arange(5000000,45000000,1000000)
    # plot velocity near shearing wall
    dp.save_plot_velocity_ave_z_strain(
        log_variable_dic_list,
        51,
        velocityinputstepsarray,
        0,
        1, r'$\gamma$',
        abs(float(log_variable['in_velocity'])), "Vx near shear wall",
        mv_v_name="mv_1",
        if_on_paper=False,
        if_include_0_y_axis=True,
        n_sum_over_axis = 2,
        fixtimeave_id_name='avspatial_ave',
        ifrotate_tick=True,
        ifshrink=False,
        filename='f0007' + 'velocity_near_shear_wall.png',
    )
    # plot velocity near static wall
    dp.save_plot_velocity_ave_z_strain(
        log_variable_dic_list,
        51,
        velocityinputstepsarray,
        -1,
        1, r'$\gamma$',
        abs(float(log_variable['in_velocity'])), "Vx near static wall",
        mv_v_name="mv_1",
        if_on_paper=False,
        if_include_0_y_axis=True,
        n_sum_over_axis = 2,
        fixtimeave_id_name='avspatial_ave',
        ifrotate_tick=True,
        ifshrink=False,
        filename='f0007' + 'velocity_near_static_wall.png',
    )
    dp.save_plot_velocity_select_z_vs_y(
        log_variable_dic_list,
        201,
        np.array([velocityinputstepsarray[-1]]),
        [0,7,14,15,16],
        float(log_variable['dp']), 'y',
        abs(float(log_variable['in_velocity'])), 'Vx',
        mv_v_name="mv_1",
        if_on_paper=False,
        if_include_0_y_axis=True,
        fixtimeave_id_name='avspatial_ave',
        ifrotate_tick=True,
        ifshrink=False,
        filename='f0007' + 'velocity_select_z_vs_y1.png',
    )
    dp.save_plot_velocity_select_z_vs_y(
        log_variable_dic_list,
        201,
        np.array([velocityinputstepsarray[-1]]),
        [0,7,14,15,16],
        float(log_variable['dp']), 'y',
        abs(float(log_variable['in_velocity'])), 'Vy',
        mv_v_name="mv_2",
        if_on_paper=False,
        if_include_0_y_axis=True,
        fixtimeave_id_name='avspatial_ave',
        ifrotate_tick=True,
        ifshrink=False,
        filename='f0007' + 'velocity_select_z_vs_y2.png',
    )
    dp.save_plot_velocity_select_z_vs_y(
        log_variable_dic_list,
        201,
        np.array([velocityinputstepsarray[-1]]),
        [0,7,14,15,16],
        float(log_variable['dp']), 'y',
        abs(float(log_variable['in_velocity'])), 'Vz',
        mv_v_name="mv_3",
        if_on_paper=False,
        if_include_0_y_axis=True,
        fixtimeave_id_name='avspatial_ave',
        ifrotate_tick=True,
        ifshrink=False,
        filename='f0007' + 'velocity_select_z_vs_y3.png',
    )
    # plot wall average wall stress
    
    add_pre_filename_string = 'f0008'
    for (force_v_name, fixtimeave_id_name, force_scale_factor, fig_y_label, filename) in [
        ('inwall_force_1', 'timeav_inwall_force', sidewall_force_scale, 'shearing wall ' + r'$\sigma_{21}$', add_pre_filename_string + 'inwall_force_1.png'),
        ('inwall_force_2', 'timeav_inwall_force', sidewall_force_scale, 'shearing wall ' + r'$\sigma_{22}$', add_pre_filename_string + 'inwall_force_2.png'),
        ('inwall_force_3', 'timeav_inwall_force', sidewall_force_scale, 'shearing wall ' + r'$\sigma_{23}$', add_pre_filename_string + 'inwall_force_3.png'),
        ('outwall_force_1', 'timeav_outwall_force', sidewall_force_scale, 'static wall ' + r'$\sigma_{21}$', add_pre_filename_string + 'outwall_force_1.png'),
        ('outwall_force_2', 'timeav_outwall_force', sidewall_force_scale, 'static wall ' + r'$\sigma_{22}$', add_pre_filename_string + 'outwall_force_2.png'),
        ('outwall_force_3', 'timeav_outwall_force', sidewall_force_scale, 'static wall ' + r'$\sigma_{23}$', add_pre_filename_string + 'outwall_force_3.png'),
        ('zbottom_force_1', 'timeav_zbottom_force', bottomewall_force_scale, 'bottom ' + r'$\sigma_{31}$', add_pre_filename_string + 'zbottom_force_1.png'),
        ('zbottom_force_2', 'timeav_zbottom_force', bottomewall_force_scale, 'bottom ' + r'$\sigma_{32}$', add_pre_filename_string + 'zbottom_force_2.png'),
        ('zbottom_force_3', 'timeav_zbottom_force', bottomewall_force_scale, 'bottom ' + r'$\sigma_{33}$', add_pre_filename_string + 'zbottom_force_3.png'),
    ]:
        dp.save_plot_total_wall_force_strain(
            log_variable_dic_list,
            51,
            wallinputstepsarray,
            1, r'$\gamma$',
            force_scale_factor, fig_y_label,
            force_v_name,
            if_on_paper=if_on_paper,
            if_include_0_y_axis=True,
            fixtimeave_id_name=fixtimeave_id_name,
            ifrotate_tick=True,
            ifshrink=False,
            filename=filename
        )
    
    # plot force ratio wall
    add_pre_filename_string = 'f0009'
    for (force_v1_name, force_v2_name, fixtimeave_id_name, force_scale_factor, fig_y_label, filename) in [
        ('inwall_force_1', 'inwall_force_2', 'timeav_inwall_force', 1, 'shearing wall ' + r'$\mu_{12}$', add_pre_filename_string + 'inwall_force_ratio_12.png'),
        ('outwall_force_1', 'outwall_force_2', 'timeav_outwall_force', 1, 'static wall ' + r'$\mu_{12}$', add_pre_filename_string + 'outwall_force_ratio_12.png'),
        ('zbottom_force_1', 'zbottom_force_3', 'timeav_zbottom_force', 1, 'bottom ' + r'$\mu_{13}$', add_pre_filename_string + 'zbottom_force_ratio_13.png'),
    ]:
        dp.save_plot_total_wall_force_ratio_strain(
            log_variable_dic_list,
            51,
            wallinputstepsarray,
            1, r'$\gamma$',
            force_scale_factor, fig_y_label,
            force_v1_name, force_v2_name,
            if_on_paper=if_on_paper,
            if_include_0_y_axis=True,
            fixtimeave_id_name=fixtimeave_id_name,
            ifrotate_tick=True,
            ifshrink=False,
            filename=filename
        )
    
    # plot wall stress vs coord
    add_pre_filename_string = 'f0010'
    for (v_name, fixtimeave_id_name, fig_y_label, filename) in [
        ('chunk_inwall_force_1', 'ave_std_inwall', 'inwall_stress ' + r'$\sigma_{21}$', add_pre_filename_string + 'inwall_stress_21_coord.png'),
        ('chunk_inwall_force_2', 'ave_std_inwall', 'inwall_stress ' + r'$\sigma_{22}$', add_pre_filename_string + 'inwall_stress_22_coord.png'),
        ('chunk_inwall_force_3', 'ave_std_inwall', 'inwall_stress ' + r'$\sigma_{23}$', add_pre_filename_string + 'inwall_stress_23_coord.png'),
        ('chunk_outwall_force_1', 'ave_std_outwall', 'outwall_stress ' + r'$\sigma_{21}$', add_pre_filename_string + 'outwall_stress_21_coord.png'),
        ('chunk_outwall_force_2', 'ave_std_outwall', 'outwall_stress ' + r'$\sigma_{22}$', add_pre_filename_string + 'outwall_stress_22_coord.png'),
        ('chunk_outwall_force_3', 'ave_std_outwall', 'outwall_stress ' + r'$\sigma_{23}$', add_pre_filename_string + 'outwall_stress_23_coord.png'),
        ]:
        dp.save_plot_wall_stress_coord(
            log_variable_dic_list,
            201,
            np.arange(3000000,45000000,10000000),
            v_name,
            'Coord2',
            fixtimeave_id_name,
            'z',
            fig_y_label,
            float(log_variable['dp']),
            stress_scale_height,
            if_on_paper=if_on_paper,
            ifrotate_tick=True,
            ifshrink=False,
            filename = filename,
        )
    for (v_name, fixtimeave_id_name, fig_y_label, filename) in [
        ('chunk_zbottom_force_1', 'ave_std_zbottom', 'zbottom_stress ' + r'$\sigma_{31}$', add_pre_filename_string + 'zbottom_stress_31_coord.png'),
        ('chunk_zbottom_force_2', 'ave_std_zbottom', 'zbottom_stress ' + r'$\sigma_{32}$', add_pre_filename_string + 'zbottom_stress_32_coord.png'),
        ('chunk_zbottom_force_3', 'ave_std_zbottom', 'zbottom_stress ' + r'$\sigma_{33}$', add_pre_filename_string + 'zbottom_stress_33_coord.png'),
        ]:
        dp.save_plot_wall_stress_coord(
            log_variable_dic_list,
            201,
            np.arange(3000000,45000000,10000000),
            v_name,
            'Coord1',
            fixtimeave_id_name,
            'y',
            fig_y_label,
            float(log_variable['dp']),
            stress_scale_height,
            if_on_paper=if_on_paper,
            ifrotate_tick=True,
            ifshrink=False,
            filename = filename,
        )
    # plot wall stress ratio vs coord
    add_pre_filename_string = 'f0011'
    for (v_name_1, v_name_2, fixtimeave_id_name, fig_y_label, filename) in [
        ('chunk_inwall_force_1', 'chunk_inwall_force_2', 'ave_std_inwall', 'shearing wall ' + r'$\mu_{12}$', add_pre_filename_string + 'inwall_stress_ratio_12_coord.png'),
        ('chunk_inwall_force_3', 'chunk_inwall_force_2', 'ave_std_inwall', 'shearing wall ' + r'$\mu_{32}$', add_pre_filename_string + 'inwall_stress_ratio_32_coord.png'),
        ('chunk_outwall_force_1', 'chunk_outwall_force_2', 'ave_std_outwall', 'static wall ' + r'$\mu_{12}$', add_pre_filename_string + 'outwall_stress_ratio_12_coord.png'),
        ('chunk_outwall_force_3', 'chunk_outwall_force_2', 'ave_std_outwall', 'static wall ' + r'$\mu_{32}$', add_pre_filename_string + 'outwall_stress_ratio_32_coord.png'),
        ]:
        dp.save_plot_wall_stress_ratio_coord(
            log_variable_dic_list,
            201,
            np.arange(3000000,45000000,10000000),
            v_name_1,
            v_name_2,
            'Coord2',
            fixtimeave_id_name,
            'z',
            fig_y_label,
            float(log_variable['dp']),
            1,
            if_on_paper=if_on_paper,
            ifrotate_tick=True,
            ifshrink=False,
            filename = filename,
        )
    for (v_name_1, v_name_2, fixtimeave_id_name, fig_y_label, filename) in [
        ('chunk_zbottom_force_1', 'chunk_zbottom_force_3', 'ave_std_zbottom', 'zbottom' + r'$\mu_{13}$', add_pre_filename_string + 'zbottom_stress_ratio_31_coord.png'),
        ('chunk_zbottom_force_2', 'chunk_zbottom_force_3', 'ave_std_zbottom', 'zbottom' + r'$\mu_{23}$', add_pre_filename_string + 'zbottom_stress_ratio_32_coord.png'),
        ]:
        dp.save_plot_wall_stress_ratio_coord(
            log_variable_dic_list,
            201,
            np.arange(3000000,45000000,10000000),
            v_name_1,
            v_name_2,
            'Coord1',
            fixtimeave_id_name,
            'y',
            fig_y_label,
            float(log_variable['dp']),
            1,
            if_on_paper=if_on_paper,
            ifrotate_tick=True,
            ifshrink=False,
            filename = filename,
        )
    
    # plot wall stress

    # plot mu-I
    coord1 = ddfp.get_coord_by_variable(len(log_variable_dic_list)-1, 'Coord1', log_variable_dic_list, 'avspatialstress_ave')
    d_coord1 = coord1[1,0] - coord1[0,0]
    coord2 = ddfp.get_coord_by_variable(len(log_variable_dic_list)-1, 'Coord2', log_variable_dic_list, 'avspatialstress_ave')
    d_coord2 = coord2[0,1] - coord2[0,0]
    binvolume = (
        d_coord1
        *d_coord2
        *float(log_variable_dic_list[-1]['x_period_dp_unit'])
        *float(log_variable_dic_list[-1]['dp'])
    )
    stress_scale_width = float(log_variable['den'])*float(log_variable['g'])*float(log_variable['width_wall_dp_unit'])*float(log_variable['dp'])
    
    for (inputstepsarray, add_pre_filename_string, array_index_y) in [
        (np.array([35000000]), 'f0012', np.array([0])),
        (np.array([3000000]), 'f0013', np.array([0])),
        (np.array([35000000]), 'f0014', np.array([1])),
        (np.array([3000000]), 'f0015', np.array([1])),
    ]:
        for (v_name, filename, fig_y_label) in [
            ('stress_multiply_binvolume_11', add_pre_filename_string + 'stress_multiply_binvolume_11_near_shearing', r'$\sigma_{11}$' + 'near_shearing'),
            ('stress_multiply_binvolume_22', add_pre_filename_string + 'stress_multiply_binvolume_22_near_shearing', r'$\sigma_{22}$' + 'near_shearing'),
            ('stress_multiply_binvolume_33', add_pre_filename_string + 'stress_multiply_binvolume_33_near_shearing', r'$\sigma_{33}$' + 'near_shearing'),
            ('stress_multiply_binvolume_12', add_pre_filename_string + 'stress_multiply_binvolume_12_near_shearing', r'$\sigma_{12}$' + 'near_shearing'),
            ('stress_multiply_binvolume_13', add_pre_filename_string + 'stress_multiply_binvolume_13_near_shearing', r'$\sigma_{13}$' + 'near_shearing'),
            ('stress_multiply_binvolume_23', add_pre_filename_string + 'stress_multiply_binvolume_23_near_shearing', r'$\sigma_{23}$' + 'near_shearing'),
        ]:
            v_name_x_axis = 'Coord2'
            x_scale_factor = float(log_variable['dp'])
            fig_x_label = 'z'
            array_index_z = np.arange(10)

            dp.save_plot_variable_vs_y_or_z_or_strain(
                log_variable_dic_list,
                201,
                inputstepsarray,
                x_scale_factor, fig_x_label,
                stress_scale_height*binvolume, fig_y_label,
                v_name, v_name_x_axis,
                array_index_y,
                array_index_z,
                if_on_paper=False,
                ifrotate_tick=True,
                ifshrink=False,
                fixtimeave_id_fortime='avspatialstress_ave',
                fixtimeave_id_forcoord='avspatialstress_ave',
                fixtimeave_id_name='avspatialstress_ave',
                is_std=False,
                is_calculated_v=False,
                filename=filename,
            )
    for (inputstepsarray, add_pre_filename_string, array_index_y) in [
        (np.array([35000000]), 'f0016', np.array([-1])),
        (np.array([3000000]), 'f0017', np.array([-1])),
        (np.array([35000000]), 'f0018', np.array([-2])),
        (np.array([3000000]), 'f0019', np.array([-2])),
    ]:
        for (v_name, filename, fig_y_label) in [
            ('stress_multiply_binvolume_11', add_pre_filename_string + 'stress_multiply_binvolume_11_near_static', r'$\sigma_{11}$' + 'near_static'),
            ('stress_multiply_binvolume_22', add_pre_filename_string + 'stress_multiply_binvolume_22_near_static', r'$\sigma_{22}$' + 'near_static'),
            ('stress_multiply_binvolume_33', add_pre_filename_string + 'stress_multiply_binvolume_33_near_static', r'$\sigma_{33}$' + 'near_static'),
            ('stress_multiply_binvolume_12', add_pre_filename_string + 'stress_multiply_binvolume_12_near_static', r'$\sigma_{12}$' + 'near_static'),
            ('stress_multiply_binvolume_13', add_pre_filename_string + 'stress_multiply_binvolume_13_near_static', r'$\sigma_{13}$' + 'near_static'),
            ('stress_multiply_binvolume_23', add_pre_filename_string + 'stress_multiply_binvolume_23_near_static', r'$\sigma_{23}$' + 'near_static'),
        ]:
            v_name_x_axis = 'Coord2'
            x_scale_factor = float(log_variable['dp'])
            fig_x_label = 'z'
            array_index_z = np.arange(10)

            dp.save_plot_variable_vs_y_or_z_or_strain(
                log_variable_dic_list,
                201,
                inputstepsarray,
                x_scale_factor, fig_x_label,
                stress_scale_height*binvolume, fig_y_label,
                v_name, v_name_x_axis,
                array_index_y,
                array_index_z,
                if_on_paper=False,
                ifrotate_tick=True,
                ifshrink=False,
                fixtimeave_id_fortime='avspatialstress_ave',
                fixtimeave_id_forcoord='avspatialstress_ave',
                fixtimeave_id_name='avspatialstress_ave',
                is_std=False,
                is_calculated_v=False,
                filename=filename,
            )
    
    for (inputstepsarray, add_pre_filename_string, array_index_y) in [
        (np.array([50200000]), 'f0020', np.array([0])),
        (np.array([5050000]), 'f0021', np.array([0])),
        (np.array([5100000]), 'f0022', np.array([0])),
        (np.array([4000000]), 'f0023', np.array([0])),
        (np.array([50200000]), 'f0024', np.array([1])),
        (np.array([5050000]), 'f0025', np.array([1])),
        (np.array([5100000]), 'f0026', np.array([1])),
        (np.array([4000000]), 'f0027', np.array([1])),
    ]:
        for (v_name, filename, fig_y_label) in [
            ('omega_1', add_pre_filename_string + 'omega_1_near_shearing', r'$\omega_{1}$' + 'near_shearing'),
            ('omega_2', add_pre_filename_string + 'omega_2_near_shearing', r'$\omega_{2}$' + 'near_shearing'),
            ('omega_3', add_pre_filename_string + 'omega_3_near_shearing', r'$\omega_{3}$' + 'near_shearing'),
        ]:
            v_name_x_axis = 'Coord2'
            x_scale_factor = float(log_variable['dp'])
            fig_x_label = 'z'
            array_index_z = np.arange(10)

            dp.save_plot_variable_vs_y_or_z_or_strain(
                log_variable_dic_list,
                1,
                inputstepsarray,
                x_scale_factor, fig_x_label,
                1, fig_y_label,
                v_name, v_name_x_axis,
                array_index_y,
                array_index_z,
                if_on_paper=False,
                ifrotate_tick=True,
                ifshrink=False,
                fixtimeave_id_fortime='avspatial_omega_ave',
                fixtimeave_id_forcoord='avspatial_omega_ave',
                fixtimeave_id_name='avspatial_omega_ave',
                is_std=False,
                is_calculated_v=False,
                filename=filename,
            )
    for (inputstepsarray, add_pre_filename_string, array_index_y) in [
        (np.array([50200000]), 'f0028', np.array([-1])),
        (np.array([5050000]), 'f0029', np.array([-1])),
        (np.array([5100000]), 'f0030', np.array([-1])),
        (np.array([4000000]), 'f0031', np.array([-1])),
        (np.array([50200000]), 'f0032', np.array([-2])),
        (np.array([5050000]), 'f0033', np.array([-2])),
        (np.array([5100000]), 'f0034', np.array([-2])),
        (np.array([4000000]), 'f0035', np.array([-2])),
    ]:
        for (v_name, filename, fig_y_label) in [
            ('omega_1', add_pre_filename_string + 'omega_1_near_static', r'$\omega_{1}$' + 'near_static'),
            ('omega_2', add_pre_filename_string + 'omega_2_near_static', r'$\omega_{2}$' + 'near_static'),
            ('omega_3', add_pre_filename_string + 'omega_3_near_static', r'$\omega_{3}$' + 'near_static'),
        ]:
            v_name_x_axis = 'Coord2'
            x_scale_factor = float(log_variable['dp'])
            fig_x_label = 'z'
            array_index_z = np.arange(10)

            dp.save_plot_variable_vs_y_or_z_or_strain(
                log_variable_dic_list,
                1,
                inputstepsarray,
                x_scale_factor, fig_x_label,
                1, fig_y_label,
                v_name, v_name_x_axis,
                array_index_y,
                array_index_z,
                if_on_paper=False,
                ifrotate_tick=True,
                ifshrink=False,
                fixtimeave_id_fortime='avspatial_omega_ave',
                fixtimeave_id_forcoord='avspatial_omega_ave',
                fixtimeave_id_name='avspatial_omega_ave',
                is_std=False,
                is_calculated_v=False,
                filename=filename,
            )
    """
    # new plot mu and mu for middle region
    # whole time
    dp.save_plot_mu(
        51, 8, 8, np.arange(5300000,74000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu$', 'mu-gamma-later.png',
        if_on_paper=if_on_paper,
    )
    dp.save_plot_mu(
        51, 8, 8, np.arange(5300000,74000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu$', 'mu-ignoredia-gamma-later.png',
        if_on_paper=if_on_paper,
        ignore_diagonal=True,
    )
    dp.save_plot_I(
        51, 8, 8, np.arange(5300000,74000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$I$', 'I-gamma-later.png',
        if_on_paper=if_on_paper,
    )
    dp.save_plot_mu_ij(
        3, 3, 51, 8, 8, np.arange(5300000,74000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu_{33}$', 'mu33-gamma.png',
        if_on_paper=if_on_paper,
    )
    dp.save_plot_mu_ij(
        2, 2, 51, 8, 8, np.arange(5300000,74000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu_{22}$', 'mu22-gamma.png',
        if_on_paper=if_on_paper,
    )
    dp.save_plot_mu_ij(
        1, 1, 51, 8, 8, np.arange(5300000,74000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu_{11}$', 'mu11-gamma.png',
        if_on_paper=if_on_paper,
    )
    # begin
    dp.save_plot_mu(
        1, 8, 8, np.arange(5000000,6000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu$', 'mu-gamma-begin.png',
        if_on_paper=if_on_paper,
    )
    dp.save_plot_mu(
        1, 8, 8, np.arange(5000000,6000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu$', 'mu-ignoredia-gamma-begin.png',
        if_on_paper=if_on_paper,
        ignore_diagonal=True,
    )
    dp.save_plot_I(
        1, 8, 8, np.arange(5000000,6000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$I$', 'I-gamma-begin.png',
        if_on_paper=if_on_paper,
    )
    # plot mu diagonal at begin
    dp.save_plot_mu_ij(
        3, 3, 1, 8, 8, np.arange(5000000,6000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu_{33}$', 'mu33-gamma-begin.png',
        if_on_paper=if_on_paper,
    )
    dp.save_plot_mu_ij(
        2, 2, 1, 8, 8, np.arange(5000000,6000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu_{22}$', 'mu22-gamma-begin.png',
        if_on_paper=if_on_paper,
    )
    dp.save_plot_mu_ij(
        1, 1, 1, 8, 8, np.arange(5000000,6000000,10000), log_variable_dic_list, np.s_[:,3,3],
        r'$\gamma$', r'$\mu_{11}$', 'mu11-gamma-begin.png',
        if_on_paper=if_on_paper,
    )