#!/usr/bin/env python
import sys
import plotfigure.plotchunk as pp
import datapath as dp
if_plot_to_last = bool(int(sys.argv[1]))
step1           = int(sys.argv[2])
step2           = int(sys.argv[3])
n_ave           = int(sys.argv[4])

def plotchunk_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False):
    lines = pp.lines_from_one_simu(dp.lammps_directory)
    if if_plot_to_last:
        pp.plotchunk_velocity_fraction_every(n_ave, lines, figformat="png", ifpickle=False)
    else:
        pp.plotchunk_velocity_fraction_step1to2(step1, step2, n_ave, lines, figformat="png", ifpickle=False)

plotchunk_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False)