#!/usr/bin/env python
import numpy as np
import subprocess
from pylatex import Document, Section, Subsection, Figure, SubFigure, Package, NewPage
from pylatex.utils import italic
from pylatex.utils import NoEscape
import os
import os.path
import re
import sys
import d00_utils.input_text as di
import d00_utils.read_log as dr

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def plots_path_list(log_variable_dic_list, walkpath):
    pic_path_list = []
    for dirname, dirnames, filenames in os.walk(walkpath):
        # print path to all subdirectories first.
        for subdirname in dirnames:
            print(os.path.join(dirname, subdirname))

        # print path to all filenames.
        for filename in filenames:
            pic_path_list.append(os.path.join(dirname, filename))

    pic_path_list = natural_sort(pic_path_list)
    return pic_path_list

def plots_rela_path_list(log_variable_dic_list, walkpath):
    plots_rela_path_list = []
    for dirname, dirnames, filenames in os.walk(walkpath):
        # print path to all subdirectories first.
        for subdirname in dirnames:
            print(os.path.join(dirname, subdirname))

        # print path to all filenames.
        for filename in filenames:
            path = os.path.join(dirname, filename)
            relative_path = os.path.relpath(path, walkpath)
            plots_rela_path_list.append(relative_path)

    plots_rela_path_list = natural_sort(plots_rela_path_list)
    return plots_rela_path_list

def produce_latex(log_variable_dic_list, pic_path_list, default_filepath):
    
    os.makedirs(di.latex_folder_path(log_variable_dic_list), exist_ok=True)
    doc = Document(default_filepath=default_filepath)
    doc.packages.append(Package('placeins', "section"))

    log_variable = log_variable_dic_list[-1]
    with doc.create(Section("Simulation Plot")):
        doc.append("Width= " + str(log_variable["width_wall_dp_unit"]) + " diameter")
        doc.append("\nHeight= " + str(log_variable["z_length_create_dp_unit"]) + " diameter")
        doc.append("\nPeriodic length in shear direction = " + str(log_variable["x_period_dp_unit"]) + " diameter")
        doc.append("\nSavage_number = " + str(log_variable["Sa"]))
        doc.append("\nV_wall = " + str(float(log_variable["in_velocity"])/float(log_variable["dp"])) + " diameter/second")
        for n, pic_path in enumerate(pic_path_list[:49]):
            with doc.create(Figure(position='!ht')) as pic:
                pic.add_image(pic_path)
                capstring = r'\label{fig:' + str(n) + '} This is a figure.'
                pic.add_caption(NoEscape(capstring))
    with doc.create(Section("Simulation Plot2")):
        for n, pic_path in enumerate(pic_path_list[49:]):
            with doc.create(Figure(position='!ht')) as pic:
                pic.add_image(pic_path)
                capstring = r'\label{fig:' + str(n+49) + '} This is a figure.'
                pic.add_caption(NoEscape(capstring))
        doc.append(NoEscape("\\FloatBarrier"))
        
    doc.generate_pdf(clean_tex=False)

def produce_latex_view(log_variable_dic_list):
    default_filepath=di.latex_for_view_file_path(log_variable_dic_list)
    walkpath = di.plots_for_view_folder(log_variable_dic_list)
    pic_path_list = plots_path_list(log_variable_dic_list, walkpath)
    produce_latex(log_variable_dic_list, pic_path_list, default_filepath)

def produce_latex_paper(log_variable_dic_list):
    default_filepath=di.latex_for_paper_file_path(log_variable_dic_list)
    walkpath = di.plots_for_paper_folder(log_variable_dic_list)
    pic_path_list = plots_path_list(log_variable_dic_list, walkpath)
    produce_latex(log_variable_dic_list, pic_path_list, default_filepath)


def produce_latex_side_by_side(list_of_log_variable_dic_list, pic_rela_path_list, default_filepath, ifpaper):
    walkpath_list = []
    for log_variable_dic_list in list_of_log_variable_dic_list:
        if ifpaper:
            walkpath_list.append(di.plots_for_paper_folder(log_variable_dic_list))
        else:
            walkpath_list.append(di.plots_for_view_folder(log_variable_dic_list))
    os.makedirs(di.latex_folder_path(log_variable_dic_list), exist_ok=True)
    doc = Document(default_filepath=default_filepath)
    doc.packages.append(Package('placeins', "section"))

    log_variable = log_variable_dic_list[-1]
    with doc.create(Section("Simulation Plot")):
        doc.append("Width= " + str(log_variable["width_wall_dp_unit"]) + " diameter")
        doc.append("\nHeight= " + str(log_variable["z_length_create_dp_unit"]) + " diameter")
        doc.append("\nPeriodic length in shear direction = " + str(log_variable["x_period_dp_unit"]) + " diameter")
        doc.append("\nSavage_number = " + str(log_variable["Sa"]))
        doc.append("\nV_wall = " + str(float(log_variable["in_velocity"])/float(log_variable["dp"])) + " diameter/second")
        for n, pic_rela_path in enumerate(pic_rela_path_list):

            with doc.create(Figure(position='h!')) as pic:

                with doc.create(SubFigure(
                        position='b',
                        width=NoEscape(r'0.45\linewidth'))) as subfig1:
                    pic_path = os.path.join(walkpath_list[0], pic_rela_path)
                    subfig1.add_image(pic_path, width=NoEscape(r'\linewidth'))
                    capstring = r'\label{fig:' + str(n) + '} This is a figure.'
                    subfig1.add_caption(NoEscape(capstring))
                with doc.create(SubFigure(
                        position='b',
                        width=NoEscape(r'0.45\linewidth'))) as subfig2:
                    pic_path = os.path.join(walkpath_list[1], pic_rela_path)
                    subfig2.add_image(pic_path, width=NoEscape(r'\linewidth'))
                    capstring = r'\label{fig:' + str(n) + '} This is a figure.'
                    subfig2.add_caption(NoEscape(capstring))
                with doc.create(SubFigure(
                        position='b',
                        width=NoEscape(r'0.45\linewidth'))) as subfig3:
                    pic_path = os.path.join(walkpath_list[2], pic_rela_path)
                    subfig3.add_image(pic_path, width=NoEscape(r'\linewidth'))
                    capstring = r'\label{fig:' + str(n) + '} This is a figure.'
                    subfig3.add_caption(NoEscape(capstring))
        
    doc.generate_pdf(clean_tex=False)