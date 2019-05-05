import sys
import pprint
import datetime
import numpy as np
# import project module
import check_simulation as cs
import output_control as oc
import datapath as dp
import data_arrange as da


method_list = list(range(0, 3))  # method = 0, 1, 2 are correct. method = 3~7(8) wrong

worj=sys.argv[1]

f_name_check_f_1contact = 'check_contact_force_' + worj + '_id_' + str(oc.id_i) +'step_' + str(oc.step1) + '_' + str(oc.step2)

file = dp.debug_print_path + f_name_check_f_1contact


with open (file, 'a') as f:

	np.set_printoptions(threshold = 100000)
	np.set_printoptions(precision=6)

	da.create_directory(dp.debug_print_path)
	pprint.pprint(datetime.datetime.now(), f)
	pprint.pprint("print error larger than {error_tolerence}".format(error_tolerence=oc.error_tolerence), f)
	pprint.pprint("step {first} to {second}".format(first=oc.step1, second=oc.step2), f)
	pprint.pprint("(wall_id, wallname) = {first}".format(first=cs.wall_list_id_name), f)
	

	for method_i in method_list:
		[f_step_error_array, fji_cal_in_error_step, fji_plus_in_error_step, id_j] = cs.fjwi_plus_check_multistep_1contact_fromcustom(dp.f_custom + '.h5', oc.id_i, oc.step1, oc.step2, oc.error_tolerence, method_i, worj)
		pprint.pprint("id_j {id_j}".format(id_j=id_j), f)
		pprint.pprint ("method {method_i}".format(method_i=method_i), f)
		pprint.pprint("f_step_error_array", f)
		pprint.pprint(f_step_error_array, f)
		pprint.pprint("fji_plus_in_error_step", f)
		pprint.pprint(fji_plus_in_error_step, f)
		pprint.pprint("fji_cal_in_error_step", f)
		pprint.pprint(fji_cal_in_error_step, f)
		pprint.pprint('===', f)
