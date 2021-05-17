import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data_for_plot as ddfp
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def include_0_y_axis(fig, ax):
    if ax.get_ylim()[0]*ax.get_ylim()[1] <= 0:
        pass
    else:
        if ax.get_ylim()[0] > 0:
            ax.set_ylim(bottom=0)
        else:
            ax.set_ylim(top=0)
    return (fig, ax)

def show_legend(fig, ax):
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    return (fig, ax)

def action_for_not_on_paper(fig, ax):
    (fig, ax) = show_legend(fig, ax)
    return (fig, ax)

def rotate_ticklabel_shrink(fig, ax):
    # rotate ticklabel
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    return (fig, ax)

def save_close(fig, ax, filepath):
    fig.savefig(
        filepath,
        format="png",
    )
    # close figure after save
    plt.close('all')

def api_ax_dim_xy_11(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        label = "",
        marker = ".",
        linestyle = 'None',
        markersize=12,
    ):
    ax.plot(
        x_value/x_scale_factor, y_value/y_scale_factor,
        label = label,
        marker = marker,
        linestyle = linestyle,
        markersize=markersize,
    )
    return (fig, ax)

def api_ax_dim_xy_12(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        labels_list,
        marker = ".",
        linestyle = 'None',
        markersize=12,
    ):
    if y_value.shape[0] != len(labels_list):
        sys.exit('y_value.shape[0] != len(labels_list)')

    for n, label in enumerate(labels_list):
        (fig, ax) = api_ax_dim_xy_11(
            fig, ax,
            x_value, x_scale_factor,
            y_value[n], y_scale_factor,
            label = label,
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    return (fig, ax)

def api_ax_quiver_xyqv_2222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        q_value, q_scale_factor,
        v_value, v_scale_factor,
        quiver_scale,
        label_scale,
        label="",
    ):
    Q = ax.quiver(
        x_value/x_scale_factor, y_value/y_scale_factor, q_value/q_scale_factor, v_value/v_scale_factor,
        units='width',angles='xy', scale_units='xy', scale=quiver_scale,
        cmap='hsv',
    )
    ax.quiverkey(
        Q, 0.1, 0.95, label_scale,
        label = label,
        labelpos='E',
        coordinates='figure', angle=45,
    )
    return (fig, ax)

def api_ax_streamplot_xyqv_2222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        q_value, q_scale_factor,
        v_value, v_scale_factor,
    ):
    strm = ax.streamplot(
        x_value[:,0]/x_scale_factor, y_value[0,:]/y_scale_factor,
        np.transpose(q_value/q_scale_factor), np.transpose(v_value/v_scale_factor),
        linewidth=1, color='k',
        density=[0.8, 0.8],
    )
    return (fig, ax)

def api_ax_contour_xyv_222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        vector_value, vector_scale_factor,
        contour_norm = 'linear',
        contour_v_min_max = "constant", # or "min_to_max",
        vmin = 10**-5,
        vmax = 10**0,
    ):

    d_x_value = x_value[:,0][1]-x_value[:,0][0]
    d_y_value = y_value[0,:][1]-y_value[0,:][0]
    x_value_expand = x_value[:,0] - d_x_value/2
    x_value_expand = np.append(x_value_expand, (x_value_expand[-1]+d_x_value))
    y_value_expand = y_value[0,:] - d_y_value/2
    y_value_expand = np.append(y_value_expand, (y_value_expand[-1]+d_y_value))
    # plot

    if contour_norm == 'linear' and contour_v_min_max == "min_to_max":
        colornorm = colors.Normalize(vmin=np.transpose(vector_value/vector_scale_factor).min(), vmax=np.transpose(vector_value/vector_scale_factor).max())
    elif contour_norm == 'log' and contour_v_min_max == "min_to_max":
        colornorm = colors.LogNorm(vmin=np.transpose(vector_value/vector_scale_factor).min(), vmax=np.transpose(vector_value/vector_scale_factor).max())
    elif contour_norm == 'linear' and contour_v_min_max == "constant":
        colornorm = colors.Normalize(vmin=vmin, vmax=vmax)
    elif contour_norm == 'log' and contour_v_min_max == "constant":
        colornorm = colors.LogNorm(vmin=vmin, vmax=vmax)

    contour = ax.pcolor(
        x_value_expand/x_scale_factor, y_value_expand/y_scale_factor, np.transpose(vector_value/vector_scale_factor),
        norm=colornorm,
        cmap='coolwarm',
    )
    fig.colorbar(contour, ax=ax, extend='max')
    return (fig, ax)

def api_plot_1D_1D_dim_xy_11(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    ):
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    (fig, ax) = api_ax_dim_xy_11(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        label = "",
        marker = ".",
        linestyle = 'None',
        markersize=12,
    )
    return (fig, ax)

def api_plot_1D_1D_dim_xy_12(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    labels_list,
    ):
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    (fig, ax) = api_ax_dim_xy_12(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        labels_list,
        marker = ".",
        linestyle = 'None',
        markersize=12,
    )
    return (fig, ax)

def api_quiver_2D_2D_dim_xyqv_2222(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    q_value, q_scale_factor,
    v_value, v_scale_factor,
    quiver_scale,
    label_scale,
    label="",
    ):
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    (fig, ax) = api_ax_quiver_xyqv_2222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        q_value, q_scale_factor,
        v_value, v_scale_factor,
        quiver_scale,
        label_scale,
        label="",
    )
    (fig, ax) = rotate_ticklabel_shrink(fig, ax)
    return (fig, ax)

def api_streamplot_2D_2D_dim_xyqv_2222(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    q_value, q_scale_factor,
    v_value, v_scale_factor,
    ):
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    # plot
    (fig, ax) = api_ax_streamplot_xyqv_2222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        q_value, q_scale_factor,
        v_value, v_scale_factor,
    )
    return (fig, ax)

def api_contour_2D_2D_dim_xyv_222(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    vector_value, vector_scale_factor,
    contour_norm = 'linear',
    contour_v_min_max = "constant", # or "min_to_max",
    vmin = 10**-5,
    vmax = 10**0,
    ):
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    
    (fig, ax) = api_ax_contour_xyv_222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        vector_value, vector_scale_factor,
        contour_norm = contour_norm,
        contour_v_min_max = contour_v_min_max, # or "min_to_max",
        vmin = vmin,
        vmax = vmax,
    )
    return (fig, ax)

def api_velocity_streamplot_contour_xyqv_2222(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    q_value, q_scale_factor,
    v_value, v_scale_factor,
    quiver_scale,
    label_scale,
    if_on_paper=False,
    if_include_0_y_axis=True,
    contour_norm = 'linear',
    contour_v_min_max = "constant", # or "min_to_max",
    vmin = 10**-5,
    vmax = 10**0,
    ):
    
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)

    (fig, ax) = api_ax_streamplot_xyqv_2222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        q_value, q_scale_factor,
        v_value, v_scale_factor,
    )
    vector_value = np.sqrt(q_value**2 + v_value**2)
    vector_value = np.ma.masked_where(np.logical_not(vector_value > 0), vector_value)
    (fig, ax) = api_ax_contour_xyv_222(
        fig, ax,
        x_value, x_scale_factor,
        y_value, y_scale_factor,
        vector_value, q_scale_factor,
        contour_norm = contour_norm,
        contour_v_min_max = contour_v_min_max, # or "min_to_max",
        vmin = vmin,
        vmax = vmax,
    )

    (fig, ax) = rotate_ticklabel_shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)

    return (fig, ax)

def plot_velocity_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_1",
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ):
    (coord, velocity, labels_list) = ddfp.velocity_ave_y(
        mv_v_name,
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        n_sum_over_axis = n_sum_over_axis,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    (fig, ax) = api_plot_1D_1D_dim_xy_12(
        coord, coord_scale_factor, fig_x_label,
        velocity, velocity_scale_factor, fig_y_label,
        labels_list,
    )
    (fig, ax) = rotate_ticklabel_shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
    return (fig, ax)

def plot_total_wall_force_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    force_scale_factor, fig_y_label,
    force_v_name,
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_name='timeav_inwall_force',
    ):
    (strain, total_wall_force) = ddfp.total_wall_force_strain(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        force_v_name,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    (fig, ax) = api_plot_1D_1D_dim_xy_11(
        strain, strain_scale_factor, fig_x_label,
        total_wall_force, force_scale_factor, fig_y_label,
    )
    (fig, ax) = rotate_ticklabel_shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
    return (fig, ax)

def plot_velocity_1_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ):
    (fig, ax) =  plot_velocity_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_1",
    if_on_paper=if_on_paper,
    if_include_0_y_axis=if_include_0_y_axis,
    n_sum_over_axis = n_sum_over_axis,
    fixtimeave_id_name=fixtimeave_id_name,
    )
    return (fig, ax)

def plot_velocity_2_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ):
    (fig, ax) =  plot_velocity_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_2",
    if_on_paper=if_on_paper,
    if_include_0_y_axis=if_include_0_y_axis,
    n_sum_over_axis = n_sum_over_axis,
    fixtimeave_id_name=fixtimeave_id_name,
    )
    return (fig, ax)

def plot_velocity_3_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ):
    (fig, ax) =  plot_velocity_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_3",
    if_on_paper=if_on_paper,
    if_include_0_y_axis=if_include_0_y_axis,
    n_sum_over_axis = n_sum_over_axis,
    fixtimeave_id_name=fixtimeave_id_name,
    )
    return (fig, ax)


def save_plot_velocity_streamplot_contour(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_x, x_scale_factor, fig_x_label,
    v_name_y, y_scale_factor, fig_y_label,
    v_name_q, q_scale_factor,
    v_name_v, v_scale_factor,
    quiver_scale,
    label_scale,
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name=None,
    is_calculated_v=True,
    contour_norm = 'linear',
    contour_v_min_max = "constant", # or "min_to_max",
    vmin = 10**-5,
    vmax = 10**0,
    ):
    
    (coord_1, coord_2, V_1, V_2) = ddfp.quiver_data(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_q,
    v_name_v,
    v_name_x,
    v_name_y,
    fixtimeave_id_fortime=fixtimeave_id_fortime,
    fixtimeave_id_forcoord=fixtimeave_id_forcoord,
    fixtimeave_id_name=fixtimeave_id_name,
    is_calculated_v=is_calculated_v,
    )
    for n, step in enumerate(inputstepsarray):
        (x_value, y_value, q_value, v_value) = (coord_1, coord_2, V_1[n], V_2[n])
        # plot ave_z velocity across y
        (fig, ax) = api_velocity_streamplot_contour_xyqv_2222(
            x_value, x_scale_factor, fig_x_label,
            y_value, y_scale_factor, fig_y_label,
            q_value, q_scale_factor,
            v_value, v_scale_factor,
            quiver_scale,
            label_scale,
            if_on_paper=if_on_paper,
            if_include_0_y_axis=if_include_0_y_axis,
            contour_norm = contour_norm,
            contour_v_min_max = contour_v_min_max, # or "min_to_max",
            vmin = vmin,
            vmax = vmax,
        )
        os.makedirs(di.plots_for_view_folder(log_variable_dic_list), exist_ok=True)
        os.makedirs(di.plots_for_paper_folder(log_variable_dic_list), exist_ok=True)
        filepath = di.plots_for_view_file_path(log_variable_dic_list, 'stream_' + str(step) + '.png')
        save_close(fig, ax, filepath)

