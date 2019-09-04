#!/usr/bin/env python
import sys
import datapath as dp
import createdata.datatofile as cd
import calculate.checkoutputfile as cco
import plotfigure.plotchunk as pc
for i in range(10000,3000000,10000):
    pc.plotchunk(i, dp.lammps_directory + "output/momentum_mass_field/fix.momentum_mass_field.all")