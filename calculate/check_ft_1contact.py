import sys
import pprint
import datetime
import numpy as np
# import module in simulation folder
import osmanage as om
import datapath as dp
# import module in calculate folder
import calculate.checksimulation as cs
import calculate.input as ci

# ====== input begin ========
method_list = ci.method_list  # method = 0, 1, 2 are correct. method = 3~7(8) wrong

worj = ci.worj

id_i = ci.id_i

step1 = ci.step1

step2 = ci.step2

error_tolerence = ci.error_tolerence

# ======= input end ==============

f_name_check_f_1contact = 'check_contact_force_' + worj + '_id_' + str(ci.id_i) +'step_' + str(ci.step1) + '_' + str(ci.step2)

file = dp.debug_print_path + f_name_check_f_1contact


with open (file, 'a') as f:

	np.set_printoptions(threshold = 100000)
	np.set_printoptions(precision=6)

	om.create_directory(dp.debug_print_path)
	pprint.pprint(datetime.datetime.now(), f)
	pprint.pprint("print error larger than {error_tolerence}".format(error_tolerence=ci.error_tolerence), f)
	pprint.pprint("step {first} to {second}".format(first=ci.step1, second=ci.step2), f)
	pprint.pprint("(wall_id, wallname) = {first}".format(first=cs.wall_list_id_name), f)
	

	for method_i in method_list:
		[f_step_error_array, fji_cal_in_error_step, fji_plus_in_error_step, id_j] = cs.fjwi_plus_check_multistep_1contact_fromcustom(dp.f_custom + '.h5', ci.id_i, ci.step1, ci.step2, ci.error_tolerence, method_i, worj)
		pprint.pprint("id_j {id_j}".format(id_j=id_j), f)
		pprint.pprint ("method {method_i}".format(method_i=method_i), f)
		pprint.pprint("f_step_error_array", f)
		pprint.pprint(f_step_error_array, f)
		pprint.pprint("fji_plus_in_error_step", f)
		pprint.pprint(fji_plus_in_error_step, f)
		pprint.pprint("fji_cal_in_error_step", f)
		pprint.pprint(fji_cal_in_error_step, f)
		pprint.pprint('===', f)
