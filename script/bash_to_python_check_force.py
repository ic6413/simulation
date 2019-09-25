#!/usr/bin/env python
import sys
import datapath as dp
import createdata.datatofile as cd
import calculate.checkoutputfile as cco
import plotfigure.plottofile as ppf
import read_setting.read_setting as rr

id_i = int(sys.argv[1])
step1 = int(sys.argv[2])
step2 = int(sys.argv[3])
fromtraceorall_input = sys.argv[4]
error_tolerence = 1e-7
method_list = list(range(0, 3))

cd.dumptofile(fromtraceorall=fromtraceorall_input).dump_custom()
# cehck overlap
cco.checkoverlap(id_i, step1, step2).checkprint()
# check j force, must have store pairforce when running lammps
cco.check_ft_1j_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()
# check wall force, must have store wallforce when running lammps
cco.check_ft_1w_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()
# check force for all contact
cco.checkforce(id_i, step1, step2, error_tolerence, method_list).checkprint()
