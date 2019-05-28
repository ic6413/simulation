#!/usr/bin/env python
import getattribute.create_attributes as gc
gc.define_attribute_dict()

import createdata.datatofile as cd
cd.thermo_hdf5_csv()
cd.dump_custom()
id_i = 15556
step1 = 57860
step2 = 57900

import calculate.checkoutputfile as cco
cco.checkoverlap(id_i, step1, step2).checkprint()
error_tolerence = 1e-7
method_list = list(range(0, 3))
cco.checkforce(id_i, step1, step2, error_tolerence, method_list).checkprint()
cco.check_ft_1j_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()
cco.check_ft_1w_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()

import plotfigure.plottofile as ppf
ppf.plotclass(step1, step2).plotij(id_i)

variable_name_list = 'all'
ppf.plotclass(step1, step2).plotthermo(variable_name_list)


id_i = 15556

ppf.plotclass(0, 100000).plot3Dtraj(id_i)
ppf.plotclass(0, 10000).plot3Dtraj(id_i)

for i in range(0,800000,200000):
    ppf.plotclass(i, i+200000).plot3Dtraj(id_i)
cd.dump_custom_select([id_i])
y_variables = [
        'fx',
        'fy',
        'fz',
        'x',
        'y',
        'z',
        'vx',
        'vy',
        'vz',
        'omegax',
        'omegay',
        'omegaz',
        'v_KEtr_atom',
        'c_KEr_atom',
        'c_KEt_atom'
    ]
variable_name_list = 'all'

try:
    ppf.plotclass(step1, step2).plotsingle(id_i, variable_name_list)
except FileNotFoundError:
    cd.dump_custom_select([id_i])

