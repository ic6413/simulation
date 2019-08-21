#!/usr/bin/env python
import sys
import datapath as dp
import createdata.datatofile as cd
import calculate.checkoutputfile as cco
import plotfigure.plottofile as ppf


# thermo
cd.thermo_hdf5_csv()

cd.dumptofile(fromtraceorall='trace').dump_custom()
cd.dumptofile(fromtraceorall='all').dump_custom()
step1 = 20000 #23730
step2 = 20645
id_i = 131923

error_tolerence = 1e-15
method_list = [2, 0, 1] # method 2 has smaller error than method 0, 1

cco.checkoverlap(id_i, 20000, 20645).checkprint()
cco.checkoverlap(109753, 20000, 20645).checkprint()

"""
cco.checkforce(id_i, 20000, 21000, error_tolerence, method_list).checkprint()

cd.dumptofile('trace').dump_custom_max("c_KEt_atom")

cco.checkmaxid("c_KEt_atom", step1, step2).checkprint()

cco.checkforce(109753, 20980, 21000, error_tolerence, method_list).checkprint()

# check j force, must have store pairforce when running lammps
#cco.check_ft_1j_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()
# check wall force, must have store wallforce when running lammps
#cco.check_ft_1w_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()

"""