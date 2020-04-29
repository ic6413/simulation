#!/usr/bin/env python
import sys
import numpy as np
import subprocess
from pylatex import Document, Section, Subsection, Figure, Package
from pylatex.utils import italic
from pylatex.utils import NoEscape
import os
import os.path
import read_setting as rr
import datapath as dp
import re
import os

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

n_nve = int(sys.argv[1])

# create folders
for root, dirnames, filenames in os.walk(dp.diagram_path):
    for dirname in dirnames:
        folderpath = os.path.join(root, dirname)
        createfolderpath = folderpath.replace(dp.diagram_path, dp.latex_pics_path)
        os.makedirs(createfolderpath, exist_ok=True)
# copy diagram wall    
subprocess.run(["cp", "-r", dp.f_wall_force_plot_path, dp.latex_pics_path])
subprocess.run(["cp", "-r", dp.f_wall_force_plot_from_initial_path, dp.latex_pics_path])

# copy max velocity near wall   
subprocess.run(["cp", "-r", dp.f_max_velocity_near_wall, dp.latex_pics_path])

# copy diagram velocity only for begin and last
for root, dirnames, filenames in os.walk(dp.f_momentum_mass_field_path):
    # if not empty
    if filenames:
        sorted_filenames = natural_sort(filenames)
        first_filename = sorted_filenames[0]
        last_filename = sorted_filenames[-1]
        copy_filenames = [first_filename, last_filename]
        for filename in copy_filenames:
            filepath = os.path.join(root, filename)
            f_copy_path = filepath.replace(dp.diagram_path, dp.latex_pics_path)
            subprocess.run(["cp", filepath, f_copy_path])


pics_list = []
for root, dirnames, filenames in os.walk(dp.latex_path):
    if n_nve != 1:
        re_key = re.compile("/nve_" + str(n_nve) + "/")
        if re_key.search(root+"/"):
            for filename in filenames:
                relapath_root = root.replace(dp.latex_path, '')
                pics_list.append(os.path.join(relapath_root, filename))
    elif n_nve == 1:
        re_key = re.compile("/nve_[0-9]*/")
        if not re_key.search(root+"/"):
            for filename in filenames:
                relapath_root = root.replace(dp.latex_path, '')
                pics_list.append(os.path.join(relapath_root, filename))
        re_key = re.compile("/nve_" + "1" + "/")
        if re_key.search(root+"/"):
            for filename in filenames:
                relapath_root = root.replace(dp.latex_path, '')
                pics_list.append(os.path.join(relapath_root, filename))
os.chdir(dp.latex_path)
doc = Document(default_filepath=dp.latex_path+"plots_" + "nve_" + str(n_nve))
doc.packages.append(Package('placeins',"section"))


with doc.create(Section("Simulation Plot")):
    doc.append("Width= " + str(rr.logfile["width_wall_dp_unit"]) + " diameter")
    doc.append("\nHeight= " + str(rr.logfile["z_length_create_dp_unit"]) + " diameter")
    doc.append("\nPeriodic length in shear direction = " + str(rr.logfile["x_period_dp_unit"]) + " diameter")
    doc.append("\nSavage_number = " + str(rr.logfile["Sa"]))
    doc.append("\nV_wall = " + str(float(rr.logfile["in_velocity"])/float(rr.logfile["dp"])) + " diameter/second")
    pic_path_wall_from_initial = natural_sort([pic for pic in pics_list if ("wall_force_from_initial" in pic)])
    pic_path_wall = natural_sort([pic for pic in pics_list if ("wall_force" in pic) and ("wall_force_from_initial" not in pic)])
    pic_path_vfield = natural_sort([pic for pic in pics_list if ("momentum_mass_field" in pic)])
    pic_path_maxvy = natural_sort([pic for pic in pics_list if ("max_velocity_near_wall" in pic)])
    with doc.create(Subsection('Wall')):
        for pic_path in pic_path_wall:
            with doc.create(Figure(position='h!')) as pic:
                pic.add_image(pic_path)
                pic.add_caption(NoEscape("\label{fig:epsart}"))
    doc.append(NoEscape("\\FloatBarrier"))
    with doc.create(Subsection('Velocity Field and Density')):
        for pic_path in pic_path_vfield:
            with doc.create(Figure(position='h!')) as pic:
                pic.add_image(pic_path)
                pic.add_caption(NoEscape("\label{fig:epsart}"))
    doc.append(NoEscape("\\FloatBarrier"))
    with doc.create(Subsection('Wall from Initial')):
        for pic_path in pic_path_wall_from_initial:
            with doc.create(Figure(position='h!')) as pic:
                pic.add_image(pic_path)
                pic.add_caption(NoEscape("\label{fig:epsart}"))
    doc.append(NoEscape("\\FloatBarrier"))
    with doc.create(Subsection('Maximum of Vy near wall')):
        for pic_path in pic_path_maxvy:
            with doc.create(Figure(position='h!')) as pic:
                pic.add_image(pic_path)
                pic.add_caption(NoEscape("\label{fig:epsart}"))

doc.generate_pdf(clean_tex=False)

os.chdir(rr.lammps_directory)