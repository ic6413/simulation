#!/usr/bin/env python
import sys
import datapath as dp
import createdata.datatofile as cd
import calculate.checkoutputfile as cco
import plotfigure.plotchunk as pc

if_plot_to_last = bool(int(sys.argv[1]))
step2           = int(sys.argv[3])
step1           = int(sys.argv[2])
pc.plotchunk(if_plot_to_last, step1, step2)