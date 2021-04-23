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

def produce_latex(subfolder_under_diagram, latexpdfname):
    pic_path_list = []
    for dirname, dirnames, filenames in os.walk(dp.diagram_path + subfolder_under_diagram):
        # print path to all subdirectories first.
        for subdirname in dirnames:
            print(os.path.join(dirname, subdirname))

        # print path to all filenames.
        for filename in filenames:
            pic_path_list.append(os.path.join(dirname, filename))

    pic_path_list = natural_sort(pic_path_list)
    os.chdir(dp.latex_path)
    doc = Document(default_filepath=dp.latex_path+latexpdfname)
    doc.packages.append(Package('placeins', "section"))

    #breakpoint()

    with doc.create(Section("Simulation Plot")):
        doc.append("Width= " + str(rr.logfile["width_wall_dp_unit"]) + " diameter")
        doc.append("\nHeight= " + str(rr.logfile["z_length_create_dp_unit"]) + " diameter")
        doc.append("\nPeriodic length in shear direction = " + str(rr.logfile["x_period_dp_unit"]) + " diameter")
        doc.append("\nSavage_number = " + str(rr.logfile["Sa"]))
        doc.append("\nV_wall = " + str(float(rr.logfile["in_velocity"])/float(rr.logfile["dp"])) + " diameter/second")
        
           
        for pic_path in pic_path_list:
            with doc.create(Figure(position='h!')) as pic:
                #breakpoint()
                pic.add_image(pic_path)
                pic.add_caption(NoEscape("\label{fig:epsart}"))
        
        doc.append(NoEscape("\\FloatBarrier"))
        


    doc.generate_pdf(clean_tex=False)

    os.chdir(rr.lammps_directory)
# main exclusive
def main():
    produce_latex('onpaper/', "plots_onpaper")
    
    produce_latex('notonpaper/', "plots_notonpaper")
# main exclusive
if __name__ == "__main__":
    main()