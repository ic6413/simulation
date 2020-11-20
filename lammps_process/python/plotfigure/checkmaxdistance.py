# import site package
import os
import sys
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datapath as dp
# import my package
import read_setting as rr
diameter = float(rr.logfile["dp"])
lmp_path = rr.folder_path_list_initial_to_last[-1]
# python read line
def extractdata(timestep, variablename):
    filepath_rela_lmp_path = "output/single_all/dump.all.single." + str(timestep)
    
    with open(lmp_path + filepath_rela_lmp_path) as f:
        # line list
        lines = f.read().strip().split('\n')
    n_line_in_file = len(lines)
    first_9_lines = lines[0: 9]
    # header
    header = first_9_lines[8].split()[2:]

    # calculate index of line to delete
    # start from timestep row
    lines[:] = lines[9:]

    # creating a list
    lines = [line.split() for line in lines]
    lines = np.asarray(lines, dtype=np.float64, order='F')
    index_id = header.index('id')
    index_variable = header.index(variablename)

    id = lines[:, index_id]
    ind = np.argsort(id)
    value = lines[:, index_variable]
    value = np.take_along_axis(value, ind, axis=0)
    return value

# max distance
def distance(s1, s2):
    dx = extractdata(s1, 'x') - extractdata(s2, 'x')
    dy = extractdata(s1, 'y') - extractdata(s2, 'y')
    dz = extractdata(s1, 'z') - extractdata(s2, 'z')
    dyz = (dy**2+dz**2)**0.5
    return [extractdata(s1, 'y'), extractdata(s1, 'z'), dyz]

def plotstaticzone(s1, s2, rate_lessthanave=0.01):
    [y, z, dyz] = distance(s1, s2)
    ave = np.average(dyz)
    ifstatic = (dyz < rate_lessthanave*ave)
    nonstatic = np.logical_not(ifstatic)
    static_y = ifstatic*y
    static_z = ifstatic*z
    nonstatic_y = nonstatic*y
    nonstatic_z = nonstatic*z
    # normalize
    static_y = static_y/diameter
    static_z = static_z/diameter
    nonstatic_y = nonstatic_y/diameter
    nonstatic_z = nonstatic_z/diameter
    
    fig, ax = plt.subplots()
    ax.set_title(
        "staticzone \n(" + "moving distance smaller than " + str(rate_lessthanave) + "*average)" + "\n step " + str(s1) + " to " + str(s2)
    )
    ax.plot(
        static_y, static_z,
        marker = ".",
        linestyle = 'None',
        markersize=1,
    )
    
    ax.plot(
        nonstatic_y, nonstatic_z,
        marker = ".",
        linestyle = 'None',
        markersize=1,
    )
    
    subfoldername = "staticzone/"
    os.makedirs(dp.diagram_path + subfoldername, exist_ok=True)
    fig.savefig(
    "_".join([
        dp.diagram_path + subfoldername + 'step',
        str(s1),
        str(s2),
    ]),
    format="png",
    )

# define static region
mask = (y-9) - (9-0.5)/(15.5-10.5)*(x-15.5) < 0

# define nonstatic region
mask = (y-25) - (25-0.5)/(15.5-5)*(x-5) > 0

# calculate variable for static region and nonstatic region

    def plot_1D_from_chunk2D(
            n_ave, x_name, y_name, inputstepsarray, coord1_index_array, coord2_index_array,
            ave_over_coord1=False, ave_over_coord2=False,
            figure_class=None, legend_class=None, spaceave=None,
            x_scale_factor=1, x_scale_str=None, y_scale_factor=1, y_scale_str=None,
            x_scale = 'linear', y_scale = 'linear',
            useerrorbar = True,
            ifdivideNcount = False,
            ifmaskstatic=False,
            ifmasknonstatic=False,
        ):
        char = c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
        Coord1 = load_coord(char, "Coord1")
        Coord2 = load_coord(char, "Coord2")
        x = Coord1/diameter
        y = Coord2/diameter
        masknonstatic = (y-9) - (9-0.5)/(15.5-10.5)*(x-15.5) < 0
        maskstatic = (y-25) - (25-0.5)/(15.5-5)*(x-5) > 0
        if ifmaskstatic:
            maskstaticornot = maskstatic
        if ifmasknonstatic:
            maskstaticornot = masknonstatic
        # str for variable used in figure label legend
        x_name_label_in_figure = transfer_v_name_to_label_str_in_figure(x_name).replace('_middle', '')
        y_name_label_in_figure = transfer_v_name_to_label_str_in_figure(y_name).replace('_middle', '')
        # legend_class: time, coord1, coord2
        # spaceave: coord1, coord2
        
        time_array = time_from_start_rotate(inputstepsarray)
        if x_name == "timestep":
            x_value = inputstepsarray
            x_value_std = 0
        elif x_name == "time":
            x_value = time_array
            x_value_std = 0
        else:
            x_value = get_value(lmp_path, x_name, n_ave, inputstepsarray)
            x_value_std = get_value(lmp_path, x_name + "_std", n_ave, inputstepsarray)

        if y_name == "timestep":
            y_value = inputstepsarray
            y_value_std = 0
        elif y_name == "time":
            y_value = time_array
            y_value_std = 0
        else:
            y_value = get_value(lmp_path, y_name, n_ave, inputstepsarray)
            y_value_std = get_value(lmp_path, y_name + "_std", n_ave, inputstepsarray)
        if ifdivideNcount:
            Ncount_value = get_value(lmp_path, 'Ncount', n_ave, inputstepsarray)
        # scale factor
        x_value = x_value/x_scale_factor
        x_value_std = x_value_std/x_scale_factor
        y_value = y_value/y_scale_factor
        y_value_std = y_value_std/y_scale_factor

        # subfolder name
        subfoldername = y_name + "_" + x_name + "/"
        os.makedirs(dp.diagram_path + subfoldername, exist_ok=True)
        if useerrorbar:
            os.makedirs(dp.diagram_path + subfoldername + "errorbar/", exist_ok=True)

        # plot ave_z velocity across y
        fig, ax = plt.subplots()
        # title

        if x_scale_str is None:
            x_label_str = x_name_label_in_figure
        else:
            x_label_str = x_name_label_in_figure + " (" + x_scale_str + ")"
        ax.set_xlabel(x_label_str)
        
        if y_scale_str is None:
            y_label_str = y_name_label_in_figure
        else:
            y_label_str = y_name_label_in_figure + " (" + y_scale_str + ")"
        ax.set_ylabel(y_label_str)
        if ifdivideNcount:
            ax.set_ylabel(y_label_str + '/N')
        
        ax.set_xscale(x_scale)
        ax.set_yscale(y_scale)
        char = c_r_npyfilepath_coord_char[y_name]["coordinate_characteristic"]
        Coord1 = load_coord(char, "Coord1")
        Coord2 = load_coord(char, "Coord2")
        # legend_class: time, coord1, coord2
        if legend_class == 'tc1c2':
            for indexstep, step in enumerate(inputstepsarray):
                for coord1_index in coord1_index_array:
                    for coord2_index in coord2_index_array:
                        time = time_from_start_rotate(step)
                        label = " ".join([
                            "t={:.2e} s".format(time),
                            "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                            "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                        ])

                        if x_name == "timestep" or x_name == "time":
                            x_value_plot = x_value[indexstep]
                            x_value_std_plot = 0
                        else:
                            x_value_plot = x_value[indexstep, coord1_index, coord2_index]
                            x_value_std_plot = x_value_std[indexstep, coord1_index, coord2_index]
                        
                        y_value_plot = y_value[indexstep, coord1_index, coord2_index]
                        y_value_std_plot = y_value_std[indexstep, coord1_index, coord2_index]
                        if ifdivideNcount:
                            Ncount_value_plot = Ncount_value[indexstep, coord1_index, coord2_index]
                            y_value_plot = y_value_plot/Ncount_value_plot
                            y_value_std_plot = y_value_std_plot/Ncount_value_plot
                        if useerrorbar:
                            ax.errorbar(
                                x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                                label=label,
                                marker = ".",
                                linestyle = 'None',
                                markersize=12,
                            )
                        else:
                            ax.plot(
                                x_value_plot, y_value_plot,
                                label=label,
                                marker = ".",
                                linestyle = 'None',
                                markersize=12,
                            )
        elif legend_class == 'c1c2':
            for coord1_index in coord1_index_array:
                for coord2_index in coord2_index_array:
                    
                    label=" ".join([
                        "t={:.2e} to {:.2e}s".format(time_array[0], time_array[-1]),
                        "y={:.1f} ({})".format(Coord1[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                        "z={:.1f} ({})".format(Coord2[coord1_index, coord2_index]/coord_scale, coord_scale_str),
                    ])

                    if x_name == "timestep" or x_name == "time":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, coord1_index, coord2_index]
                        x_value_std_plot = x_value_std[:, coord1_index, coord2_index]
                    
                    y_value_plot = y_value[:, coord1_index, coord2_index]
                    y_value_std_plot = y_value_std[:, coord1_index, coord2_index]
                    
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, coord1_index, coord2_index]
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
        elif legend_class == 'c1':
            if coord2_index_array=='ave':
                if ave_over_coord1:
                    maskstaticornot_plot = maskstaticornot[coord1_index_array, :]
                    if x_name == "timestep" or x_name == "time":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, coord1_index_array, :]
                        x_value_std_plot = x_value_std[:, coord1_index_array, :]
                        x_value_plot = np.mean(x_value_plot*maskstaticornot_plot, axis=(1,2))
                        x_value_std_plot = np.mean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_plot = y_value[:, coord1_index_array, :]
                    y_value_std_plot = y_value_std[:, coord1_index_array, :]
                    y_value_plot = np.mean(y_value_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_std_plot = np.mean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))
                    
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, coord1_index_array, :]
                        Ncount_value_plot = np.mean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    ystring = "_".join(
                        ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                    )
                    label=" ".join([
                        "t={:.2e} to {:.2e}s".format(time_array[0], time_array[-1]),
                        "average y=",
                        ystring,
                        "({})".format(coord_scale_str),
                        "z=average",
                    ])
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                else:
                    for coord1_index in coord1_index_array:
                        label=" ".join([
                            "t={:.2e} to {:.2e}s".format(time_array[0], time_array[-1]),
                            "y={:.1f} ({})".format(Coord1[coord1_index, 0]/coord_scale, coord_scale_str),
                            "z=average",
                        ])
                        maskstaticornot_plot = maskstaticornot[coord1_index, :]
                        if x_name == "timestep" or x_name == "time":
                            x_value_plot = x_value[:]
                            x_value_std_plot = 0
                        else:
                            x_value_plot = x_value[:, coord1_index, :]
                            x_value_std_plot = x_value_std[:, coord1_index, :]
                            #breakpoint()
                            x_value_plot = np.mean(x_value_plot*maskstaticornot_plot, axis=1)
                            x_value_std_plot = np.mean(x_value_std_plot*maskstaticornot_plot, axis=1)
                        y_value_plot = y_value[:, coord1_index, :]
                        y_value_std_plot = y_value_std[:, coord1_index, :]
                        y_value_plot = np.mean(y_value_plot*maskstaticornot_plot, axis=1)
                        y_value_std_plot = np.mean(y_value_std_plot*maskstaticornot_plot, axis=1)
                        if ifdivideNcount:
                            Ncount_value_plot = Ncount_value[:, coord1_index, :]
                            Ncount_value_plot = np.mean(Ncount_value_plot*maskstaticornot_plot, axis=1)
                            y_value_plot = y_value_plot/Ncount_value_plot
                            y_value_std_plot = y_value_std_plot/Ncount_value_plot
                        if useerrorbar:
                            ax.errorbar(
                                x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                                label=label,
                                marker = ".",
                                linestyle = 'None',
                                markersize=12,
                            )
                        else:
                            ax.plot(
                                x_value_plot, y_value_plot,
                                label=label,
                                marker = ".",
                                linestyle = 'None',
                                markersize=12,
                            )
        elif legend_class == 'c2':
            if coord1_index_array=='ave':
                if ave_over_coord2:
                    maskstaticornot_plot = maskstaticornot[:, coord2_index_array]
                    if x_name == "timestep" or x_name == "time":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, :, coord2_index_array]
                        x_value_std_plot = x_value_std[:, :, coord2_index_array]
                        x_value_plot = np.mean(x_value_plot*maskstaticornot_plot, axis=(1,2))
                        x_value_std_plot = np.mean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_plot = y_value[:, :, coord2_index_array]
                    y_value_std_plot = y_value_std[:, :, coord2_index_array]
                    y_value_plot = np.mean(y_value_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_std_plot = np.mean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, :, coord2_index_array]
                        Ncount_value_plot = np.mean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    zstring = "_".join(
                        ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                    )
                    label=" ".join([
                        "t={:.2e} to {:.2e}s".format(time_array[0], time_array[-1]),
                        "y=average",
                        "average z=",
                        zstring,
                        "({})".format(coord_scale_str),
                    ])
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                else:
                    for coord2_index in coord2_index_array:
                        
                        label=" ".join([
                            "t={:.2e} to {:.2e}s".format(time_array[0], time_array[-1]),
                            "y=average",
                            "z={:.1f} ({})".format(Coord2[0, coord2_index]/coord_scale, coord_scale_str),
                        ])
                        maskstaticornot_plot = maskstaticornot[:, coord2_index]
                        if x_name == "timestep" or x_name == "time":
                            x_value_plot = x_value[:]
                            x_value_std_plot = 0
                        else:
                            x_value_plot = x_value[:, :, coord2_index]
                            x_value_std_plot = x_value_std[:, :, coord2_index]
                            x_value_plot = np.mean(x_value_plot*maskstaticornot_plot, axis=1)
                            x_value_std_plot = np.mean(x_value_std_plot*maskstaticornot_plot, axis=1)
                        y_value_plot = y_value[:, :, coord2_index]
                        y_value_std_plot = y_value_std[:, :, coord2_index]
                        y_value_plot = np.mean(y_value_plot*maskstaticornot_plot, axis=1)
                        y_value_std_plot = np.mean(y_value_std_plot*maskstaticornot_plot, axis=1)
                        if ifdivideNcount:
                            Ncount_value_plot = Ncount_value[:, :, coord2_index]
                            Ncount_value_plot = np.mean(Ncount_value_plot*maskstaticornot_plot, axis=1)
                            y_value_plot = y_value_plot/Ncount_value_plot
                            y_value_std_plot = y_value_std_plot/Ncount_value_plot
                        if useerrorbar:
                            ax.errorbar(
                                x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                                label=label,
                                marker = ".",
                                linestyle = 'None',
                                markersize=12,
                            )
                        else:
                            ax.plot(
                                x_value_plot, y_value_plot,
                                label=label,
                                marker = ".",
                                linestyle = 'None',
                                markersize=12,
                            )
        elif legend_class == None:
            if ave_over_coord1:
                if ave_over_coord2:
                    maskstaticornot_plot = maskstaticornot[coord1_index_array, :][:, coord2_index_array]
                    if x_name == "timestep" or x_name == "time":
                        x_value_plot = x_value[:]
                        x_value_std_plot = 0
                    else:
                        x_value_plot = x_value[:, coord1_index_array, :][:, :, coord2_index_array]
                        x_value_std_plot = x_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                        x_value_plot = np.mean(x_value_plot*maskstaticornot_plot, axis=(1,2))
                        x_value_std_plot = np.mean(x_value_std_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_plot = y_value[:, coord1_index_array, :][:, :, coord2_index_array]
                    y_value_std_plot = y_value_std[:, coord1_index_array, :][:, :, coord2_index_array]
                    y_value_plot = np.mean(y_value_plot*maskstaticornot_plot, axis=(1,2))
                    y_value_std_plot = np.mean(y_value_std_plot*maskstaticornot_plot, axis=(1,2))
                    if ifdivideNcount:
                        Ncount_value_plot = Ncount_value[:, coord1_index_array, :][:, :, coord2_index_array]
                        Ncount_value_plot = np.mean(Ncount_value_plot*maskstaticornot_plot, axis=(1,2))
                        y_value_plot = y_value_plot/Ncount_value_plot
                        y_value_std_plot = y_value_std_plot/Ncount_value_plot
                    zstring = "_".join(
                        ["{:.1f}".format(Coord2[0, coord2_index]/coord_scale) for coord2_index in coord2_index_array]
                    )
                    ystring = "_".join(
                        ["{:.1f}".format(Coord1[coord1_index, 0]/coord_scale) for coord1_index in coord1_index_array]
                    )
                    label=" ".join([
                        "t={:.2e} to {:.2e}s".format(time_array[0], time_array[-1]),
                        "average y=",
                        ystring,
                        "average z=",
                        zstring,
                        "({})".format(coord_scale_str),
                    ])
                    
                    if useerrorbar:
                        ax.errorbar(
                            x_value_plot, y_value_plot, xerr=x_value_std_plot, yerr=y_value_std_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                    else:
                        ax.plot(
                            x_value_plot, y_value_plot,
                            label=label,
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                        )
                
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.setp(ax.xaxis.get_minorticklabels(), rotation=30)
        plt.setp(ax.get_xticklabels(), rotation=30)
        if useerrorbar:
            fig.savefig(
            "_".join([
                dp.diagram_path + subfoldername + "errorbar/" + 'step',
                transfer_time_to_str(inputstepsarray),
                'coord1',
                transfer_coor_to_str(coord1_index_array, ave_over_coord1),
                'coord2',
                transfer_coor_to_str(coord2_index_array, ave_over_coord2),
                ".png",
            ]),
            format="png",
            )
        else:
            fig.savefig(
                "_".join([
                    dp.diagram_path + subfoldername + 'step',
                    transfer_time_to_str(inputstepsarray),
                    'coord1',
                    transfer_coor_to_str(coord1_index_array, ave_over_coord1),
                    'coord2',
                    transfer_coor_to_str(coord2_index_array, ave_over_coord2),
                    ".png",
                ]),
                format="png",
            )
        # close figure after save
        plt.close('all')

