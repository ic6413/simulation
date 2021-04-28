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

# main exclusive
def main():
    
    # copy diagram
    subprocess.run(["cp", "-r", dp.diagram_path, dp.latex_pics_path])

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

    os.chdir(dp.latex_path)
    doc = Document(default_filepath=dp.latex_path+"plots")
    doc.packages.append(Package('placeins', "section"))

    #breakpoint()

    with doc.create(Section("Simulation Plot")):
        doc.append("Width= " + str(rr.log_variable["width_wall_dp_unit"]) + " diameter")
        doc.append("\nHeight= " + str(rr.log_variable["z_length_create_dp_unit"]) + " diameter")
        doc.append("\nPeriodic length in shear direction = " + str(rr.log_variable["x_period_dp_unit"]) + " diameter")
        doc.append("\nSavage_number = " + str(rr.log_variable["Sa"]))
        doc.append("\nV_wall = " + str(float(rr.log_variable["in_velocity"])/float(rr.log_variable["dp"])) + " diameter/second")
        list_plot_no_errorbar = [
            ("Coord1", "mv_1"),
        ]
        list_plot_errorbar = [
            ("I_12", "mu_12_middle"),
            ("I_tensor", "mu_tensor_12"),
            ("strain", "I_12"),
            ("strain", "mu_12_middle"),
            ("strain", "I_tensor"),
            ("strain", "mu_tensor_12"),
            ("strain", "fraction"),
            ("strain", "velocity_1"),
            ("strain", "velocity_2"),
            ("strain", "velocity_3"),
            ("strain", "strain_rate_21_middle"),
            ("strain", "strain_rate_31_middle"),
            ("strain", "strain_rate_22_middle"),
            ("strain", "strain_rate_32_middle"),
            ("strain", "strain_rate_23_middle"),
            ("strain", "strain_rate_33_middle"),
            ("strain", "stress_11"),
            ("strain", "stress_22"),
            ("strain", "stress_33"),
            ("strain", "stress_12"),
            ("strain", "stress_13"),
            ("strain", "stress_23"),
            ("fraction", "mu_12_middle"),
            ("fraction", "I_12"),
            ("fraction", "mu_tensor_12"),
            ("fraction", "I_tensor"),
            ("strain", "fraction"),
            ("strain", "inwall_stress_1"),
            ("strain", "inwall_stress_2"),
            ("strain", "inwall_stress_3"),
            ("strain", "outwall_stress_1"),
            ("strain", "outwall_stress_2"),
            ("strain", "outwall_stress_3"),
            ("strain", "zbottom_stress_1"),
            ("strain", "zbottom_stress_2"),
            ("strain", "zbottom_stress_3"),
        ]
        for (x_name, y_name) in list_plot_no_errorbar:
            folder_name =  y_name + "_" + x_name + "/"
            folder_path = dp.latex_pics_path + "diagram/n_ave_51/" + folder_name
            onlyfilepaths = [
                os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))
            ]
            pic_path_list = natural_sort(onlyfilepaths)

            subsection_title = y_name + " vs. " + x_name
            #breakpoint()
            with doc.create(
                Subsection(
                    subsection_title.replace('_middle', '')
                )
            ):
                for pic_path in pic_path_list:
                    with doc.create(Figure(position='h!')) as pic:
                        #breakpoint()
                        pic.add_image(pic_path)
                        pic.add_caption(NoEscape("\label{fig:epsart}"))
                
            doc.append(NoEscape("\\FloatBarrier"))
        

        for (x_name, y_name) in list_plot_errorbar:
            folder_name =  y_name + "_" + x_name + "/"
            folder_path = dp.latex_pics_path + "diagram/n_ave_51/" + folder_name
            onlyfilepaths = [
                os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))
            ]
            pic_path_list = natural_sort(onlyfilepaths)

            folder_path = dp.latex_pics_path + "diagram/n_ave_51/" + folder_name + "errorbar/"
            onlyfilepaths = [
                os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))
            ]
            pic_path_list_errorbar = natural_sort(onlyfilepaths)

            subsection_title = y_name + " vs. " + x_name
            #breakpoint()
            with doc.create(
                Subsection(
                    subsection_title.replace('_middle', '')
                )
            ):
                for pic_path in pic_path_list:
                    with doc.create(Figure(position='h!')) as pic:
                        #breakpoint()
                        pic.add_image(pic_path)
                        pic.add_caption(NoEscape("\label{fig:epsart}"))
                
            doc.append(NoEscape("\\FloatBarrier"))
        
            with doc.create(
                Subsection(
                    subsection_title.replace('_middle', '') + " (with errorbar)"
                )
            ):
                for pic_path in pic_path_list_errorbar:
                    with doc.create(Figure(position='h!')) as pic:
                        pic.add_image(pic_path)
                        pic.add_caption(NoEscape("\label{fig:epsart}"))
                
            doc.append(NoEscape("\\FloatBarrier"))


    doc.generate_pdf(clean_tex=False)

    os.chdir(rr.lammps_directory)
    
# main exclusive
if __name__ == "__main__":
    main()