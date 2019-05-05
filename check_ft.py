import pprint
import datetime
import numpy as np
# import project module
import check_simulation as cs
import output_control as oc
import datapath as dp
import data_arrange as da

check_f_function = cs.fjwi_plus_check_multistep_multicontact_fromcustom  # or cs.fjwi_plus_check_multistep_multicontact_fromcustom_v1

method_list = list(range(0, 3))

file = dp.debug_print_path + oc.f_name_check_f

np.set_printoptions(threshold = 100000)
np.set_printoptions(precision=6)

da.create_directory(dp.debug_print_path)
pprint.pprint(datetime.datetime.now(), open(file, 'a'))
pprint.pprint("print error larger than {error_tolerence}".format(error_tolerence=oc.error_tolerence) , open(file, 'a'))
pprint.pprint("step {first} to {second}".format(first=oc.step1, second=oc.step2) , open(file, 'a'))

for method_i in method_list:

	[f_step_error_array, fi_cal_in_error_step, fi_plus_in_error_step] = check_f_function(dp.f_custom + '.h5', oc.id_i, oc.step1, oc.step2, oc.error_tolerence, method=method_i)
	pprint.pprint("method {method_i}".format(method_i=method_i) , open(file, 'a'))
	pprint.pprint("f_step_error_array", open(file, 'a'))
	pprint.pprint(f_step_error_array, open(file, 'a'))
	pprint.pprint("fi_plus_in_error_step", open(file, 'a'))
	pprint.pprint(fi_plus_in_error_step, open(file, 'a'))
	pprint.pprint("fi_cal_in_error_step", open(file, 'a'))
	pprint.pprint(fi_cal_in_error_step, open(file, 'a'))
	pprint.pprint('===', open(file, 'a'))
