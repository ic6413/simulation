import numpy as np
import subprocess

from pylatex import Document, Section, Subsection, Figure
from pylatex.utils import italic
from pylatex.utils import NoEscape
import os
import os.path

import read_setting.read_setting as rr
import datapath as dp

"""
subprocess.run(["cp", "-r", dp.f_wall_force_plot_path, dp.latex_pics_path])
"""

for field_folder in [
    dp.f_momentum_mass_field_v23x23_path,
    dp.f_momentum_mass_field_v13x23_path,
    dp.f_momentum_mass_field_density_x23_path,]:
    if os.path.isdir(dp.f_momentum_mass_field_path):
        if dp.f_momentum_mass_field_path in field_folder:
            string1=field_folder.replace(dp.f_momentum_mass_field_path, '')
    if os.path.isdir(dp.f_momentum_mass_field_samescale_path):
        if dp.f_momentum_mass_field_samescale_path in field_folder:
            string1=field_folder.replace(dp.f_momentum_mass_field_samescale_path, '')
    subprocess.run(["mkdir", dp.latex_pics_path + string1])
    first_step = int(rr.logfile['rst_from']) + int(rr.logfile['freq_ave_chunk_momentum_mass_field'])
    end_step = int(rr.logfile['rst_from']) + int(rr.logfile['runstep']) - int(rr.logfile['freq_ave_chunk_momentum_mass_field'])
    pics_first = field_folder + str(first_step) + ".png"
    pics_end = field_folder + str(end_step) + ".png"
    if os.path.exists(pics_first):
        subprocess.run(["cp", pics_first, dp.latex_pics_path + string1])
    if os.path.exists(pics_end):
        subprocess.run(["cp", pics_end, dp.latex_pics_path + string1])

pics_list = []
for root, dirnames, filenames in os.walk(dp.latex_path):
    for filename in [f for f in filenames if f.endswith(".png")]:
        relapath_root = root.replace(dp.latex_path, '')
        pics_list.append(os.path.join(relapath_root, filename))

if __name__ == '__main__':
    os.chdir(dp.latex_path)
    doc = Document(default_filepath=dp.latex_path+"full")

    with doc.create(Section('The fancy stuff')):

        with doc.create(Subsection('pictures')):
            for pic_path in pics_list:
                with doc.create(Figure(position='h!')) as pic:
                    pic.add_image(pic_path)
                    pic.add_caption(NoEscape("\label{fig:epsart} A figure caption. The figure captions are automatically numbered."))

    doc.generate_pdf('full', clean_tex=False)

os.chdir(dp.lammps_directory)