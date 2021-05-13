import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data_for_plot as ddfp
import matplotlib.pyplot as plt

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

def api_plot_1D_1D_dim_xy_11(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    if_on_paper=False,
    if_include_0_y_axis=False,
    ):
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    ax.plot(
        x_value/x_scale_factor, y_value/y_scale_factor,
        marker = ".",
        linestyle = 'None',
        markersize=12,
    )
    (fig, ax) = rotate_ticklabel_shrink(fig, ax)
    if not if_on_paper:
        action_for_not_on_paper(fig, ax)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
    return (fig, ax)

def api_plot_1D_1D_dim_xy_12(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    labels_list,
    if_on_paper=False,
    if_include_0_y_axis=False,
    ):
    # plot ave_z velocity across y
    fig, ax = plt.subplots()

    # xy label
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    for n, label in enumerate(labels_list):
        ax.plot(
            x_value/x_scale_factor, y_value[n]/y_scale_factor,
            label = label,
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    (fig, ax) = rotate_ticklabel_shrink(fig, ax)
    if not if_on_paper:
        action_for_not_on_paper(fig, ax)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
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
        if_on_paper=if_on_paper,
        if_include_0_y_axis=if_include_0_y_axis,
    )
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
    fixtimeave_id_name='avspatial_ave',
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
        if_on_paper=if_on_paper,
        if_include_0_y_axis=if_include_0_y_axis,
    )
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
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
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
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
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
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    )
    return (fig, ax)