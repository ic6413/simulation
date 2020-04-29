#!/usr/bin/env python
import sys
n_ave = int(sys.argv[1])

import numpy as np

import plotfigure.plotchunk as pp
import plotfigure.plotchunkstress as pps
import read_setting as rr

###chunk current
chunkobject = pp.chunk(n_ave, rr.lammps_directory)
current_steps = chunkobject.first_extra_middle_last_steps

"""
chunkobject.save_v23_x23(current_steps)
chunkobject.save_v13_x23(current_steps)
chunkobject.save_fraction_x23(current_steps)
chunkobject.save_strain_rate_ij_x23(current_steps, 0, 1)
chunkobject.save_strain_rate_ij_x23(current_steps, 2, 1)
"""
###chunk include pre
chunk_all = pp.chunk_include_pre(n_ave, rr.lammps_directory)
allsteps_since_rotate = chunk_all.allsteps_since_rotate
stepsarray_all = chunk_all.first_extra_middle_last_steps_everysimu
#stepsarray_all_tmp = stepsarray_all[np.array([0,1,4,7])]
## use allsteps since rotate to plot near wall velocity
"""
chunkobject.save_plotchunk_velocity_i_time_near_wall_ave(allsteps_since_rotate, 0, figformat="png", ifpickle=False, ifmanysimu=True)
chunkobject.save_plotchunk_velocity_i_time_near_wall_ave(allsteps_since_rotate, 1, figformat="png", ifpickle=False, ifmanysimu=True)
chunkobject.save_plotchunk_velocity_i_time_near_wall_ave(allsteps_since_rotate, 2, figformat="png", ifpickle=False, ifmanysimu=True)

###chunk include pre

chunk_all.save_plotchunk_strain_rate_ij_ave_k_ave(stepsarray_all, 0, 1, 2, figformat="png", ifpickle=False)
chunk_all.save_plotchunk_fraction_ave_j_ave(stepsarray_all, 2, 1, figformat="png", ifpickle=False)
"""
for i in range(3):
    chunk_all.save_plotchunk_velocity_i_ave_j_xk_ave(stepsarray_all, i, 2, 1, figformat="png", ifpickle=False)
"""
for i in range(3):
    chunk_all.save_plotchunk_ekovermass_i_ave_j_xk_ave(stepsarray_all_tmp, i, 2, 1, figformat="png", ifpickle=False)
    chunk_all.save_plotchunk_ekminusekaveovermass_i_ave_j_ave(stepsarray_all_tmp, i, 2, 1, figformat="png", ifpickle=False)
"""


###chunk stress current
chunkobject = pps.chunk(n_ave, rr.lammps_directory)
current_steps = chunkobject.first_extra_middle_last_steps
chunkobject.plotchunk_1(current_steps)
