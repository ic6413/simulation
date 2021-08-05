import os
import sys
import numpy as np
import d00_utils.input_text as di
import d00_utils.read_log as dr
import d00_utils.data_for_plot as ddfp
import d00_utils.calculate_new_variable as dcn
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# setting matplotlib
# Here are some matplotlib style settings that I (JSW) find useful for publication-quality plots
# http://tobackgroup.physics.tamu.edu/useful/matplotlib-formatting.txt
# First, we set the width of a figure to be exactly the width of a single column in PRL or PRD (in
# inches), and set the aspect ratio to be 4:3.

plt.rc('figure', figsize=(3.375, 3.375*(.75)), titlesize=10)

# This line sets the resolution (number of dots per inch) for saving figures to raster formats like
# PNG.  It's mainly useful to set the size of the figures displayed in the webbrowser in the Jupyter
# or IPython notebook
plt.rc('savefig', dpi=180)
#plt.rc('axes.formatter', limits=(0,0), use_mathtext=True)
plt.rc('savefig', bbox='tight')
plt.ticklabel_format(style='sci', axis='both')

# The default for scatter plots is to put three points in each legend entry.  Usually one is
# sufficient.
#plt.rc('legend', scatterpoints=1, fontsize=10, frameon=True)

# We would like the text to be created using TeX, so we get full support for LaTeX math syntax.
#plt.rc('text', usetex=True)

# We want all the text to be Computer Modern, and 10pt.  This way, combined with the correct setting
# of the size of the figure, the text in the figure will exactly match the text in the rest of the
# manuscript.

plt.rc('font', size=10, weight='normal', family='sans-serif', serif=['Helvetica', 'Arial'])

plt.rc('axes', labelsize=12, titlesize=8, labelweight='normal')

# http://aeturrell.com/2018/01/31/publication-quality-plots-in-python/
#plt.rc('xtick', labelsize=10)
#plt.rc('ytick', labelsize=10)
plt.rc('figure', autolayout=False)
plt.rc('lines', linewidth=2)
plt.rc('lines', markersize=4)
plt.rc('mathtext', fontset="stix")

# more setting
# http://physicalmodelingwithpython.blogspot.com/2015/06/making-plots-for-publication.html

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

def rotate_ticklabel(fig, ax):
    # rotate ticklabel
    plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
    plt.setp(ax.get_xticklabels(), rotation=30)
    return (fig, ax)

def shrink(fig, ax):
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    return (fig, ax)

def save_close(fig, ax, filepath, ifaddpng=True):
    if ifaddpng:
        if filepath[-3:] != 'png':
            filepath = filepath + '.png'
    fig.savefig(
        filepath,
        format="png",
    )
    # close figure after save
    plt.close('all')
    
def create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=False,
    ):

    if if_on_paper:
        os.makedirs(di.plots_for_paper_folder(log_variable_dic_list), exist_ok=True)
        filepath = di.plots_for_paper_file_path(log_variable_dic_list, filename)
        os.makedirs(di.plots_for_sync_paper_folder(log_variable_dic_list), exist_ok=True)
        filepath_dropbox = di.dropbox_plots_for_paper_file_path(log_variable_dic_list, filename)
    else:
        os.makedirs(di.plots_for_view_folder(log_variable_dic_list), exist_ok=True)
        filepath = di.plots_for_view_file_path(log_variable_dic_list, filename)
        os.makedirs(di.plots_for_sync_view_folder(log_variable_dic_list), exist_ok=True)
        filepath_dropbox = di.dropbox_plots_for_view_file_path(log_variable_dic_list, filename)
    save_close(fig, ax, filepath)
    save_close(fig, ax, filepath_dropbox)

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
        errorlimit=0.01,
    ):
    # check if y_value alomost equally spaced
    y_space = y_value[0,1:]-y_value[0,:-1]
    first_y_space = y_space[0]
    error = np.abs((y_space - first_y_space)/y_space)
    iferror_array = (error > errorlimit)
    iferror = np.any(iferror_array)
    if iferror:
        sys.exit('y_value is not alomost equally spaced')
    else:
        fixed_y = y_value[0, 0] + first_y_space*np.arange(y_value.shape[1])
    strm = ax.streamplot(
        x_value[:,0]/x_scale_factor, fixed_y/y_scale_factor,
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
    ifrotate_tick=True,
    ifshrink=False,
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
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    return (fig, ax)

def api_streamplot_2D_2D_dim_xyqv_2222(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    q_value, q_scale_factor,
    v_value, v_scale_factor,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    titlelabel="",
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
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        plt.title(titlelabel)
    return (fig, ax)

def api_contour_2D_2D_dim_xyv_222(
    x_value, x_scale_factor, fig_x_label,
    y_value, y_scale_factor, fig_y_label,
    vector_value, vector_scale_factor,
    contour_norm = 'linear',
    contour_v_min_max = "constant", # or "min_to_max",
    vmin = 10**-5,
    vmax = 10**0,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    titlelabel="",
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
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        plt.title(titlelabel)
    
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
    ifrotate_tick=True,
    ifshrink=False,
    titlelabel="",
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

    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        plt.title(titlelabel)
    return (fig, ax)

def plot_variable_vs_y_or_z_or_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    x_scale_factor, fig_x_label,
    y_scale_factor, fig_y_label,
    v_name, v_name_x_axis,
    array_index_y,
    array_index_z,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name=None,
    is_std=False,
    is_calculated_v=False,
    ):
    # get 3D value
    value = ddfp.get_ave_value(n_ave, fixtimeave_id_fortime, v_name, inputstepsarray, log_variable_dic_list, fixtimeave_id_name=fixtimeave_id_name, is_std=is_std, is_calculated_v=is_calculated_v,)
    # get 1D value
    value = value[np.arange(len(inputstepsarray)) , array_index_y, array_index_z]
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    coord1 = ddfp.get_coord_by_variable(len(log_variable_dic_list)-1, "Coord1", log_variable_dic_list, fixtimeave_id_name)[:, 0][array_index_y]
    coord2 = ddfp.get_coord_by_variable(len(log_variable_dic_list)-1, "Coord2", log_variable_dic_list, fixtimeave_id_name)[0, :][array_index_z]
    if len(inputstepsarray) != 1 and len(array_index_y) == 1 and len(array_index_z) == 1 and v_name_x_axis == 'strain':
        x_value = strain
    if len(inputstepsarray) == 1 and len(array_index_y) != 1 and len(array_index_z) == 1 and v_name_x_axis == 'Coord1':
        x_value = coord1
    if len(inputstepsarray) == 1 and len(array_index_y) == 1 and len(array_index_z) != 1 and v_name_x_axis == 'Coord2':
        x_value = coord2
    else:
        sys.exit('len of array wrong')
    
    (fig, ax) = api_plot_1D_1D_dim_xy_11(
        x_value, x_scale_factor, fig_x_label,
        value, y_scale_factor, fig_y_label,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        diameter = float(log_variable['dp'])
        if len(inputstepsarray) != 1 and len(array_index_y) == 1 and len(array_index_z) == 1 and v_name_x_axis == 'strain':
            title_1 = "y = {:.2f}, z = {:.2f}".format(coord1[0]/diameter, coord2[0]/diameter)
        if len(inputstepsarray) == 1 and len(array_index_y) != 1 and len(array_index_z) == 1 and v_name_x_axis == 'Coord1':
            title_1 = r'$\gamma = $' + '{:.3e}'.format(strain[0]) + "z = {:.2f}".format(coord2[0]/diameter)
        if len(inputstepsarray) == 1 and len(array_index_y) == 1 and len(array_index_z) != 1 and v_name_x_axis == 'Coord2':
            title_1 = r'$\gamma = $' + '{:.3e}'.format(strain[0]) + "y = {:.2f}".format(coord1[0]/diameter)
        else:
            sys.exit('len of array wrong')
        titlelabel = (
            title_1
            + "\n each point average over strain = " + '{:.4e}'.format(
                n_ave*ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_fortime, log_variable_dic_list)*float(log_variable["ts"])*shear_rate
            )
        )
        plt.title(titlelabel)

    return (fig, ax)

def plot_velocity_ave_z(
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
    ifrotate_tick=True,
    ifshrink=False,
    ):
    (coord, velocity, labels_list) = ddfp.velocity_ave_z(
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
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        titlelabel = (
            "each point average over strain = " + '{:.4e}'.format(
                n_ave
                *ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_name, log_variable_dic_list)
                *float(log_variable["ts"])
                *shear_rate
            )
        )
        plt.title(titlelabel)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
    return (fig, ax)

def plot_velocity_ave_z_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    index_coord1,
    coord_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_1",
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
    ):
    (coord, velocity, labels_list) = ddfp.velocity_ave_z(
        mv_v_name,
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        n_sum_over_axis = n_sum_over_axis,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    strain_array = ddfp.strain_from_rotate_start(inputstepsarray, log_variable_dic_list[-1])
    (fig, ax) = api_plot_1D_1D_dim_xy_11(
        strain_array, 1, fig_x_label,
        velocity[:, index_coord1], velocity_scale_factor, fig_y_label,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        titlelabel = (
            "each point average over strain = " + '{:.4e}'.format(
                n_ave
                *ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_name, log_variable_dic_list)
                *float(log_variable["ts"])
                *shear_rate
            )
        )
        plt.title(titlelabel)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
    return (fig, ax)

def plot_velocity_select_z_vs_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    array_index_coord2,
    coord_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_1",
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
    ):
    velocity = ddfp.velocity(
        mv_v_name,
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    velocity = velocity[:, :, array_index_coord2]
    velocity = velocity.reshape((-1, len(array_index_coord2)))
    velocity = np.transpose(velocity)
    strain_array = ddfp.strain_from_rotate_start(inputstepsarray, log_variable_dic_list[-1])
    str_gamma_list = [r'$\gamma = $' + '{:.3e}'.format(strain) for strain in strain_array]
    coord2_array = ddfp.get_coord_by_variable(len(log_variable_dic_list)-1, "Coord2", log_variable_dic_list, fixtimeave_id_name)[0, array_index_coord2]
    str_coord2_list = ['height' + '{:.3e}'.format(coord2/float(log_variable_dic_list[-1]['dp'])) for coord2 in coord2_array]
    labels_list=[]
    coord1 = ddfp.get_coord_by_variable(len(log_variable_dic_list)-1, "Coord1", log_variable_dic_list, fixtimeave_id_name)[:, 0]
    for str1 in str_gamma_list:
        for str2 in str_coord2_list:
            labels_list.append(str1 + ",_" + str2)
    (fig, ax) = api_plot_1D_1D_dim_xy_12(
        coord1, coord_scale_factor, fig_x_label,
        velocity, velocity_scale_factor, fig_y_label,
        labels_list,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        titlelabel = (
            "each point average over strain = " + '{:.4e}'.format(
                n_ave
                *ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_name, log_variable_dic_list)
                *float(log_variable["ts"])
                *shear_rate
            )
        )
        plt.title(titlelabel)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
    return (fig, ax)

def plot_wall_stress_coord(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name,
    coord_name,
    fixtimeave_id_name,
    fig_x_label,
    fig_y_label,
    x_scale_factor,
    y_scale_factor,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    ):
    (stress, coord) = ddfp.chunk_wall_stress_coord(
        log_variable_dic_list,
        fixtimeave_id_name,
        v_name,
        n_ave,
        inputstepsarray,
        coord_name,
    )
    strain_array = ddfp.strain_from_rotate_start(inputstepsarray, log_variable_dic_list[-1])
    labels_list = [r'$\gamma = $' + '{:.3e}'.format(strain) for strain in strain_array]
    (fig, ax) = api_plot_1D_1D_dim_xy_12(
        coord, x_scale_factor, fig_x_label,
        stress, y_scale_factor, fig_y_label,
        labels_list,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))

        titlelabel = (
            "each curve result from average over strain = " + '{:.4e}'.format(
                n_ave
                *ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_name, log_variable_dic_list)
                *float(log_variable["ts"])
                *shear_rate
            )
        )
        plt.title(titlelabel)
    return (fig, ax)

def plot_wall_stress_ratio_coord(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_1,
    v_name_2,
    coord_name,
    fixtimeave_id_name,
    fig_x_label,
    fig_y_label,
    x_scale_factor,
    y_scale_factor,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    ifabs=True,
    ):
    (stress_ratio, coord) = ddfp.chunk_wall_stress_ratio_coord(
        log_variable_dic_list,
        fixtimeave_id_name,
        v_name_1,
        v_name_2,
        n_ave,
        inputstepsarray,
        coord_name,
    )
    strain_array = ddfp.strain_from_rotate_start(inputstepsarray, log_variable_dic_list[-1])
    labels_list = [r'$\gamma = $' + '{:.2e}'.format(strain) for strain in strain_array]
    if ifabs:
        stress_ratio = np.abs(stress_ratio)
    (fig, ax) = api_plot_1D_1D_dim_xy_12(
        coord, x_scale_factor, fig_x_label,
        stress_ratio, y_scale_factor, fig_y_label,
        labels_list,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))

        titlelabel = (
            "each curve result from average over strain = " + '{:.4e}'.format(
                n_ave
                *ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_name, log_variable_dic_list)
                *float(log_variable["ts"])
                *shear_rate
            )
        )
        plt.title(titlelabel)
    return (fig, ax)

def plot_ave_value_for_select_region_yz_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    y_scale_factor, fig_y_label,
    v_name,
    n_y_0, d_n_y,
    n_z_0, d_n_z,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name=None,
    is_std=False,
    is_calculated_v=False,
    ):
    (value, first_coord_1, last_coord_1, first_coord_2, last_coord_2) = ddfp.api_ave_value_for_select_region_yz(
        v_name,
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        n_y_0, d_n_y,
        n_z_0, d_n_z,
        fixtimeave_id_fortime=fixtimeave_id_fortime,
        fixtimeave_id_forcoord=fixtimeave_id_forcoord,
        fixtimeave_id_name=fixtimeave_id_name,
        is_std=is_std,
        is_calculated_v=is_calculated_v,
    )
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    (fig, ax) = api_plot_1D_1D_dim_xy_11(
        strain, strain_scale_factor, fig_x_label,
        value, y_scale_factor, fig_y_label,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        diameter = float(log_variable['dp'])
        titlelabel = (
            "average over y = {:.2f} to {:.2f}, z = {:.2f} to {:.2f}".format(first_coord_1/diameter, last_coord_1/diameter, first_coord_2/diameter, last_coord_2/diameter)
            + "\n each point average over strain = " + '{:.4e}'.format(
                n_ave*ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_fortime, log_variable_dic_list)*float(log_variable["ts"])*shear_rate
            )
        )
        plt.title(titlelabel)

    return (fig, ax)

def plot_ave_ratio_value_for_select_region_yz_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    y_scale_factor, fig_y_label,
    v1_name, v2_name,
    n_y_0, d_n_y,
    n_z_0, d_n_z,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name1=None, fixtimeave_id_name2=None,
    is_calculated_v1=False, is_calculated_v2=False,
    ):
    (value, first_coord_1, last_coord_1, first_coord_2, last_coord_2) = ddfp.api_ave_ratio_value_for_select_region_yz(
        v1_name, v2_name,
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        n_y_0, d_n_y,
        n_z_0, d_n_z,
        fixtimeave_id_fortime=fixtimeave_id_fortime,
        fixtimeave_id_forcoord=fixtimeave_id_forcoord,
        fixtimeave_id_name1=fixtimeave_id_name1,
        fixtimeave_id_name2=fixtimeave_id_name2,
        is_calculated_v1=is_calculated_v1,
        is_calculated_v2=is_calculated_v2,
    )
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    (fig, ax) = api_plot_1D_1D_dim_xy_11(
        strain, strain_scale_factor, fig_x_label,
        value, y_scale_factor, fig_y_label,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        diameter = float(log_variable['dp'])
        titlelabel = (
            "average over y = {:.2f} to {:.2f}, z = {:.2f} to {:.2f}".format(first_coord_1/diameter, last_coord_1/diameter, first_coord_2/diameter, last_coord_2/diameter)
            + "\n each point average over strain = " + '{:.4e}'.format(
                n_ave*ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_fortime, log_variable_dic_list)*float(log_variable["ts"])*shear_rate
            )
        )
        plt.title(titlelabel)

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
    ifrotate_tick=True,
    ifshrink=False,
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
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        titlelabel = (
            "each point average over strain = " + '{:.4e}'.format(
                n_ave*ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_name, log_variable_dic_list)*float(log_variable["ts"])*shear_rate
            )
        )
        plt.title(titlelabel)
    if if_include_0_y_axis:
        (fig, ax) = include_0_y_axis(fig, ax)
    return (fig, ax)

def plot_total_wall_force_ratio_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    force_ratio_scale_factor, fig_y_label,
    force_v1_name, force_v2_name,
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_name='timeav_inwall_force',
    ifrotate_tick=True,
    ifshrink=False,
    ifabs=True,
    ):
    (strain, total_wall_force_ratio) = ddfp.total_wall_force_ratio_strain(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        force_v1_name, force_v2_name,
        fixtimeave_id_name=fixtimeave_id_name,
    )
    if ifabs:
        total_wall_force_ratio = np.abs(total_wall_force_ratio)
    (fig, ax) = api_plot_1D_1D_dim_xy_11(
        strain, strain_scale_factor, fig_x_label,
        total_wall_force_ratio, force_ratio_scale_factor, fig_y_label,
    )
    if ifrotate_tick:
        (fig, ax) = rotate_ticklabel(fig, ax)
    if ifshrink:
        (fig, ax) = shrink(fig, ax)
    if not if_on_paper:
        (fig, ax) = action_for_not_on_paper(fig, ax)
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        titlelabel = (
            "each point average over strain = " + '{:.4e}'.format(
                n_ave*ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_name, log_variable_dic_list)*float(log_variable["ts"])*shear_rate
            )
        )
        plt.title(titlelabel)
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
    ifrotate_tick=True,
    ifshrink=False,
    ):
    (fig, ax) =  plot_velocity_ave_z(
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
    ifrotate_tick=ifrotate_tick,
    ifshrink=ifshrink,
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
    ifrotate_tick=True,
    ifshrink=False,
    ):
    (fig, ax) =  plot_velocity_ave_z(
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
    ifrotate_tick=ifrotate_tick,
    ifshrink=ifshrink,
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
    ifrotate_tick=True,
    ifshrink=False,
    ):
    (fig, ax) =  plot_velocity_ave_z(
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
    ifrotate_tick=ifrotate_tick,
    ifshrink=ifshrink,
    )
    return (fig, ax)

def save_plot_ave_value_for_select_region_yz_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    y_scale_factor, fig_y_label,
    v_name,
    n_y_0, d_n_y,
    n_z_0, d_n_z,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name=None,
    is_std=False,
    is_calculated_v=False,
    filename='fraction_region',
    ):
    (fig, ax) = plot_ave_value_for_select_region_yz_strain(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        strain_scale_factor, fig_x_label,
        y_scale_factor, fig_y_label,
        v_name,
        n_y_0, d_n_y,
        n_z_0, d_n_z,
        if_on_paper=if_on_paper,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
        fixtimeave_id_fortime=fixtimeave_id_fortime,
        fixtimeave_id_forcoord=fixtimeave_id_forcoord,
        fixtimeave_id_name=fixtimeave_id_name,
        is_std=is_std,
        is_calculated_v=is_calculated_v,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
def save_plot_ave_ratio_value_for_select_region_yz_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    y_scale_factor, fig_y_label,
    v1_name, v2_name,
    n_y_0, d_n_y,
    n_z_0, d_n_z,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name1=None, fixtimeave_id_name2=None,
    is_calculated_v1=False, is_calculated_v2=False,
    filename='fraction_region',
    ):
    (fig, ax) = plot_ave_ratio_value_for_select_region_yz_strain(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        strain_scale_factor, fig_x_label,
        y_scale_factor, fig_y_label,
        v1_name, v2_name,
        n_y_0, d_n_y,
        n_z_0, d_n_z,
        if_on_paper=if_on_paper,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
        fixtimeave_id_fortime=fixtimeave_id_fortime,
        fixtimeave_id_forcoord=fixtimeave_id_forcoord,
        fixtimeave_id_name1=fixtimeave_id_name1,
        fixtimeave_id_name2=fixtimeave_id_name2,
        is_calculated_v1=is_calculated_v1,
        is_calculated_v2=is_calculated_v2,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_total_wall_force_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    force_scale_factor, fig_y_label,
    force_v_name,
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_name='timeav_inwall_force',
    ifrotate_tick=True,
    ifshrink=False,
    filename='wall',
    ):
    (fig, ax) = plot_total_wall_force_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    force_scale_factor, fig_y_label,
    force_v_name,
    if_on_paper=if_on_paper,
    if_include_0_y_axis=if_include_0_y_axis,
    fixtimeave_id_name=fixtimeave_id_name,
    ifrotate_tick=ifrotate_tick,
    ifshrink=ifshrink,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_total_wall_force_ratio_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    force_scale_factor, fig_y_label,
    force_v1_name, force_v2_name,
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_name='timeav_inwall_force',
    ifrotate_tick=True,
    ifshrink=False,
    filename='wall',
    ifabs=True,
    ):
    (fig, ax) = plot_total_wall_force_ratio_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    strain_scale_factor, fig_x_label,
    force_scale_factor, fig_y_label,
    force_v1_name, force_v2_name,
    if_on_paper=if_on_paper,
    if_include_0_y_axis=if_include_0_y_axis,
    fixtimeave_id_name=fixtimeave_id_name,
    ifrotate_tick=ifrotate_tick,
    ifshrink=ifshrink,
    ifabs=ifabs,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_velocity_1_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
    filename='velocity_1.png',
    ):
    (fig, ax) = plot_velocity_1_ave_y(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        coord1_scale_factor, fig_x_label,
        velocity_scale_factor, fig_y_label,
        if_on_paper=if_on_paper,
        if_include_0_y_axis=if_include_0_y_axis,
        n_sum_over_axis = n_sum_over_axis,
        fixtimeave_id_name=fixtimeave_id_name,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_velocity_2_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
    filename='velocity_2.png',
    ):
    (fig, ax) = plot_velocity_2_ave_y(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        coord1_scale_factor, fig_x_label,
        velocity_scale_factor, fig_y_label,
        if_on_paper=if_on_paper,
        if_include_0_y_axis=if_include_0_y_axis,
        n_sum_over_axis = n_sum_over_axis,
        fixtimeave_id_name=fixtimeave_id_name,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
def save_plot_velocity_3_ave_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
    filename='velocity_3.png',
    ):
    (fig, ax) = plot_velocity_3_ave_y(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        coord1_scale_factor, fig_x_label,
        velocity_scale_factor, fig_y_label,
        if_on_paper=if_on_paper,
        if_include_0_y_axis=if_include_0_y_axis,
        n_sum_over_axis = n_sum_over_axis,
        fixtimeave_id_name=fixtimeave_id_name,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
def save_plot_variable_vs_y_or_z_or_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    x_scale_factor, fig_x_label,
    y_scale_factor, fig_y_label,
    v_name, v_name_x_axis,
    array_index_y,
    array_index_z,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    fixtimeave_id_fortime='avspatial_ave',
    fixtimeave_id_forcoord='avspatial_ave',
    fixtimeave_id_name=None,
    is_std=False,
    is_calculated_v=False,
    filename='variable.png',
    ):
    (fig, ax) = plot_variable_vs_y_or_z_or_strain(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        x_scale_factor, fig_x_label,
        y_scale_factor, fig_y_label,
        v_name, v_name_x_axis,
        array_index_y,
        array_index_z,
        if_on_paper=if_on_paper,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
        fixtimeave_id_fortime=fixtimeave_id_fortime,
        fixtimeave_id_forcoord=fixtimeave_id_forcoord,
        fixtimeave_id_name=fixtimeave_id_name,
        is_std=is_std,
        is_calculated_v=is_calculated_v,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
def save_plot_velocity_ave_z_strain(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    index_coord1,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_1",
    if_on_paper=False,
    if_include_0_y_axis=True,
    n_sum_over_axis = 2,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
    filename='velocity_1.png',
    ):
    (fig, ax) = plot_velocity_ave_z_strain(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        index_coord1,
        coord1_scale_factor, fig_x_label,
        velocity_scale_factor, fig_y_label,
        mv_v_name=mv_v_name,
        if_on_paper=if_on_paper,
        if_include_0_y_axis=if_include_0_y_axis,
        n_sum_over_axis = n_sum_over_axis,
        fixtimeave_id_name=fixtimeave_id_name,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
def save_plot_velocity_select_z_vs_y(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    array_index_coord2,
    coord1_scale_factor, fig_x_label,
    velocity_scale_factor, fig_y_label,
    mv_v_name="mv_1",
    if_on_paper=False,
    if_include_0_y_axis=True,
    fixtimeave_id_name='avspatial_ave',
    ifrotate_tick=True,
    ifshrink=False,
    filename='velocity_select_z_vs_y.png',
    ):
    (fig, ax) = plot_velocity_select_z_vs_y(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        array_index_coord2,
        coord1_scale_factor, fig_x_label,
        velocity_scale_factor, fig_y_label,
        mv_v_name=mv_v_name,
        if_on_paper=if_on_paper,
        if_include_0_y_axis=if_include_0_y_axis,
        fixtimeave_id_name=fixtimeave_id_name,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
def save_plot_wall_stress_coord(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name,
    coord_name,
    fixtimeave_id_name,
    fig_x_label,
    fig_y_label,
    x_scale_factor,
    y_scale_factor,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    filename='stress_coord.png',
    ):

    (fig,ax) = plot_wall_stress_coord(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        v_name,
        coord_name,
        fixtimeave_id_name,
        fig_x_label,
        fig_y_label,
        x_scale_factor,
        y_scale_factor,
        if_on_paper=if_on_paper,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
        )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
def save_plot_wall_stress_ratio_coord(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_1,
    v_name_2,
    coord_name,
    fixtimeave_id_name,
    fig_x_label,
    fig_y_label,
    x_scale_factor,
    y_scale_factor,
    if_on_paper=False,
    ifrotate_tick=True,
    ifshrink=False,
    filename='stress_ratio_coord.png',
    ifabs=True,
    ):

    (fig,ax) = plot_wall_stress_ratio_coord(
        log_variable_dic_list,
        n_ave,
        inputstepsarray,
        v_name_1,
        v_name_2,
        coord_name,
        fixtimeave_id_name,
        fig_x_label,
        fig_y_label,
        x_scale_factor,
        y_scale_factor,
        if_on_paper=if_on_paper,
        ifrotate_tick=ifrotate_tick,
        ifshrink=ifshrink,
        ifabs=ifabs,
        )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )
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
    ifrotate_tick=True,
    ifshrink=False,
    filenamelist='bystep',
    add_pre_filename_string=None,
    ):

    if filenamelist=='bystep':
        filenamelist = ['stream_' + str(step) + '.png' for n, step in enumerate(inputstepsarray)]
    else:
        filenamelist=filenamelist
    if add_pre_filename_string is not None:
        filenamelist = [add_pre_filename_string + '_' + filename for filename in filenamelist]
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
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        titlelabel = (
            r'$\gamma = $' + '{:.2e}'.format(ddfp.strain_from_rotate_start(step, log_variable_dic_list[-1]))
            + "\n"
            + "each figure average over strain = " + '{:.4e}'.format(
                n_ave*ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_fortime, log_variable_dic_list)*float(log_variable["ts"])*shear_rate
            )
        )
        
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
            ifrotate_tick=ifrotate_tick,
            ifshrink=ifshrink,
            titlelabel=titlelabel,
        )
        create_folder_save(
            log_variable_dic_list,
            filenamelist[n],
            fig, ax,
            if_on_paper=if_on_paper,
        )

def save_plot_fraction_contour(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_x, x_scale_factor, fig_x_label,
    v_name_y, y_scale_factor, fig_y_label,
    v_name_q, q_scale_factor,
    if_on_paper=False,
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
    filenamelist='bystep',
    add_pre_filename_string=None,
    ):

    if filenamelist=='bystep':
        filenamelist = ['fraction_' + str(step) + '.png' for n, step in enumerate(inputstepsarray)]
    else:
        filenamelist=filenamelist
    if add_pre_filename_string is not None:
        filenamelist = [add_pre_filename_string + '_' + filename for filename in filenamelist]
    (coord_1, coord_2, Q) = ddfp.contour_data(
    log_variable_dic_list,
    n_ave,
    inputstepsarray,
    v_name_q,
    v_name_x,
    v_name_y,
    fixtimeave_id_fortime=fixtimeave_id_fortime,
    fixtimeave_id_forcoord=fixtimeave_id_forcoord,
    fixtimeave_id_name=fixtimeave_id_name,
    is_calculated_v=is_calculated_v,
    )
    Q[Q == 0] = np.nan
    for n, step in enumerate(inputstepsarray):
        (x_value, y_value, q_value) = (coord_1, coord_2, Q[n])
        
        # plot ave_z velocity across y
        log_variable = log_variable_dic_list[-1]
        shear_rate = float(log_variable['in_velocity'])/(float(log_variable['width_wall_dp_unit'])*float(log_variable['dp']))
        titlelabel = (
            r'$\gamma = $' + '{:.2e}'.format(ddfp.strain_from_rotate_start(step, log_variable_dic_list[-1]))
            + "\n"
            + "each figure average over strain = " + '{:.4e}'.format(
                n_ave*ddfp.get_d_step(len(log_variable_dic_list)-1, fixtimeave_id_fortime, log_variable_dic_list)*float(log_variable["ts"])*shear_rate
            )
        )
        (fig, ax) = api_contour_2D_2D_dim_xyv_222(
            x_value, x_scale_factor, fig_x_label,
            y_value, y_scale_factor, fig_y_label,
            q_value, q_scale_factor,
            if_on_paper=if_on_paper,
            contour_norm = contour_norm,
            contour_v_min_max = contour_v_min_max, # or "min_to_max",
            vmin = vmin,
            vmax = vmax,
            ifrotate_tick=ifrotate_tick,
            ifshrink=ifshrink,
            titlelabel=titlelabel,
        )
        create_folder_save(
            log_variable_dic_list,
            filenamelist[n],
            fig, ax,
            if_on_paper=if_on_paper,
        )

# define plot for mu I
def plot_mu_ij_I_ij(
    i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    ):
    mu_ij = dcn.mu_ij(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    I_ij = dcn.I_ij(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    fig, ax = plt.subplots()
    
    ax.plot(
        np.log10(np.abs(I_ij)),
        np.abs(mu_ij),
        marker = ".",
        linestyle = 'None',
        markersize=12,
        )
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)

    return (fig, ax)

def plot_mu_I(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    ignore_diagonal=False,
    ):
    mu = dcn.mu(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=ignore_diagonal, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    I = dcn.I(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=ignore_diagonal, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    fig, ax = plt.subplots()
    ax.plot(
        np.log10(I),
        mu,
        marker = ".",
        linestyle = 'None',
        markersize=12,
        )
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    return (fig, ax)

def plot_mu(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    ignore_diagonal=False,
    ):
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    mu = dcn.mu(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=ignore_diagonal, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    fig, ax = plt.subplots()
    ax.plot(
        strain,
        mu,
        marker = ".",
        linestyle = 'None',
        markersize=12,
        )
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    return (fig, ax)

def plot_I(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    ignore_diagonal=False,
    ):
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    I = dcn.I(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ignore_diagonal=ignore_diagonal, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    fig, ax = plt.subplots()
    ax.plot(
        strain,
        np.log10(I),
        marker = ".",
        linestyle = 'None',
        markersize=12,
        )
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    return (fig, ax)

def plot_mu_ij(
    i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    ):
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    mu_ij = dcn.mu_ij(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    fig, ax = plt.subplots()
    ax.plot(
        strain,
        mu_ij,
        marker = ".",
        linestyle = 'None',
        markersize=12,
        )
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    return (fig, ax)

def plot_I_ij(
    i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    iflog=True,
    ):
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    I_ij = dcn.I_ij(i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    fig, ax = plt.subplots()
    if iflog:
        ax.plot(
            strain,
            np.log10(np.abs(I_ij)),
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    else:
        ax.plot(
            strain,
            I_ij,
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    return (fig, ax)

def plot_trace_I(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    iflog=True,
    ):
    strain = ddfp.strain_from_rotate_start(np.array(inputstepsarray), log_variable_dic_list[-1])
    I_ij = dcn.trace_I(n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, ifcorrect_by_ave_velocity=True, ifcorrect_by_wall=True)[index]
    fig, ax = plt.subplots()
    if iflog:
        ax.plot(
            strain,
            np.log10(np.abs(I_ij)),
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    else:
        ax.plot(
            strain,
            I_ij,
            marker = ".",
            linestyle = 'None',
            markersize=12,
        )
    ax.set_xlabel(fig_x_label)
    ax.set_ylabel(fig_y_label)
    return (fig, ax)

def save_plot_mu_ij_I_ij(
    i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label, filename,
    if_on_paper=False,
    ):
    (fig, ax) = plot_mu_ij_I_ij(
    i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_mu_I(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label, filename,
    if_on_paper=False,
    ignore_diagonal=False,
    ):
    fig, ax = plot_mu_I(
        n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
        fig_x_label, fig_y_label,
        ignore_diagonal=ignore_diagonal,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_mu(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label, filename,
    if_on_paper=False,
    ignore_diagonal=False,
    ):
    fig, ax = plot_mu(
        n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
        fig_x_label, fig_y_label,
        ignore_diagonal=ignore_diagonal,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_I(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label, filename,
    if_on_paper=False,
    ignore_diagonal=False,
    ):
    fig, ax = plot_I(
        n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
        fig_x_label, fig_y_label,
        ignore_diagonal=ignore_diagonal,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_mu_ij(
    i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label, filename,
    if_on_paper=False,
    ):
    fig, ax = plot_mu_ij(
        i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
        fig_x_label, fig_y_label,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_I_ij(
    i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label, filename,
    if_on_paper=False,
    iflog=True,
    ):
    fig, ax = plot_I_ij(
        i, j, n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
        fig_x_label, fig_y_label,
        iflog=iflog,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

def save_plot_trace_I(
    n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
    fig_x_label, fig_y_label, filename,
    if_on_paper=False,
    iflog=True,
    ):
    fig, ax = plot_trace_I(
        n_ave, n_ave_coord1, n_ave_coord2, inputstepsarray, log_variable_dic_list, index,
        fig_x_label, fig_y_label,
        iflog=iflog,
    )
    create_folder_save(
        log_variable_dic_list,
        filename,
        fig, ax,
        if_on_paper=if_on_paper,
    )

# combine fig from different simu
def combine_figure(list_of_fig_ax):
    n = len(list_of_fig_ax)
    fig, axs = plt.subplots(1, n)
    for i in range(n):
        fig_sub, ax = list_of_fig_ax[i]
        axs[i] = ax
        axs[i].set_title(fig_sub.title)