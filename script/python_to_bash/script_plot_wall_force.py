#!/usr/bin/env python
import sys
import plotfigure.plot_wall_force as ppwf
if_plot_to_last = bool(int(sys.argv[1]))
step1           = int(sys.argv[2])
step2           = int(sys.argv[3])
n_ave           = int(sys.argv[4])
if n_ave == 1:
    ppwf.plot_wall_force(if_plot_to_last, step1, step2)
else:
    ppwf.plot_wall_force_ave(if_plot_to_last, step1, step2, n_ave)