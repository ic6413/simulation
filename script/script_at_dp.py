#!/usr/bin/env python
import sys
import datapath as dp
import createdata.datatofile as cd
import calculate.checkoutputfile as cco
import plotfigure.plottofile as ppf


# attribute
#gc.define_attribute_dict()

# thermo
cd.thermo_hdf5_csv()

step1 = 0
step2 = 100

# thermo plot
variable_name_list = 'all'

ppf.plotfromthermo(step1, step2).plotthermo(variable_name_list)

cd.dump_custom_max("c_KEt_atom")

cco.checkmaxid("c_KEt_atom", step1, step2).checkprint()

variable_name_list2 = ['c_KEt_atom','c_KEr_atom','v_KEtr_atom','id',]

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
        'c_KEr_atom',
        'c_KEt_atom'
    ]
ppf.plotfromtraceprint_max("c_KEt_atom").plotsingle_multifigure(["step"], y_variables, 1, 101)
"""
ppf.plotfromcustom(step1, step2).plotmaxKE_everystep(
    'c_KEt_atom', 'all'
    )
"""

"""
# custom

cd.dump_custom()
id_i = 6109
step1 = 4120248
step2 = 6000000


cco.checkoverlap(id_i, step1, step2).checkprint()

error_tolerence = 1e-7
method_list = list(range(0, 3))
#cco.checkforce(id_i, step1, step2, error_tolerence, method_list).checkprint()
cco.check_ft_1j_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()
#cco.check_ft_1w_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()

step1 = 3060000
step2 = 4000000
#ppf.plotfromcustom(step1, step2).plotij(id_i)

"""

"""
id_i = 15556

ppf.plotfromcustomselect(0, 100000, id_i).plot3Dtraj()
ppf.plotfromcustomselect(0, 10000, id_i).plot3Dtraj()

for i in range(0,800000,200000):
    ppf.plotfromcustomselect(i, i+200000, id_i).plot3Dtraj()
"""
"""
id_i = 6109
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
    ppf.plotfromcustomselect(step1, step2, id_i).plotsingle(variable_name_list)
except FileNotFoundError:
    cd.dump_custom_select([id_i])
"""
"""
id_i = 15583
try:
    ppf.plotfromcustomselect(200000, 300000, id_i).plotsingle(variable_name_list)
except FileNotFoundError:
    cd.dump_custom_select([id_i])
"""

