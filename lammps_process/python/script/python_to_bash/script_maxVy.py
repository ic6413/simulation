#!/usr/bin/env python
import sys
import plotfigure.plotchunk as pc
if_plot_to_last = bool(int(sys.argv[1]))
step1           = int(sys.argv[2])
step2           = int(sys.argv[3])
n_ave           = int(sys.argv[4])
pc.plotVymax_ave(if_plot_to_last, step1, step2, n_ave, figformat="png", ifpickle=False)