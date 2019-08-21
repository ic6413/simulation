#!/usr/bin/env python
import sys
import datapath as dp
import createdata.datatofile as cd
import calculate.checkoutputfile as cco
import plotfigure.plotchunk as pc
for i in range(20000,100,20200):
    pc.plotchunk(i, dp.lammps_directory + "output/velocity_field/fix.velocity_field.all")