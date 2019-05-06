import pprint
import datetime
import numpy as np
# import module in simulation folder
import osmanage as om
import datapath as dp
# import module in calculate folder
import calculate.checksimulation as cs
import calculate.inputvariable as ci

# input
id_i = ci.id_i
f_custom_path = ci.f_custom_path

check_f_function = cs.fjwi_plus_check_multistep_multicontact_fromcustom  # or cs.fjwi_plus_check_multistep_multicontact_fromcustom_v1

method_list = list(range(0, 3))

file = dp.debug_print_path + ci.f_check_f_name

np.set_printoptions(threshold = 100000)
np.set_printoptions(precision=6)

om.create_directory(dp.debug_print_path)
pprint.pprint(datetime.datetime.now(), open(file, 'a'))
pprint.pprint("print error larger than {error_tolerence}".format(error_tolerence=ci.error_tolerence) , open(file, 'a'))
pprint.pprint("step {first} to {second}".format(first=ci.step1, second=ci.step2) , open(file, 'a'))

for method_i in method_list:

	[f_step_error_array, fi_cal_in_error_step, fi_plus_in_error_step] = check_f_function(f_custom_path + '.h5', id_i, ci.step1, ci.step2, ci.error_tolerence, method=method_i)
	pprint.pprint("method {method_i}".format(method_i=method_i) , open(file, 'a'))
	pprint.pprint("f_step_error_array", open(file, 'a'))
	pprint.pprint(f_step_error_array, open(file, 'a'))
	pprint.pprint("fi_plus_in_error_step", open(file, 'a'))
	pprint.pprint(fi_plus_in_error_step, open(file, 'a'))
	pprint.pprint("fi_cal_in_error_step", open(file, 'a'))
	pprint.pprint(fi_cal_in_error_step, open(file, 'a'))
	pprint.pprint('===', open(file, 'a'))
