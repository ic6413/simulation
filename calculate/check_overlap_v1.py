# v1 after  write class
import pprint
import datetime
import numpy as np
# import module in simulation folder
import osmanage as om
import datapath as dp
# import module in calculate folder
import calculate.checksimulation as cs
import calculate.inputvariable as ci

f_custom_path = ci.f_custom_path

file = dp.debug_print_path + ci.f_check_overlap_name

[step_id_ifover_diffnext, ifoverlap_ij_iw_array] = cs.contact_check_multistep_v1(f_custom_path + '.h5', ci.id_i, ci.step1, ci.step2, ci.error_tolerence)

om.create_directory(dp.debug_print_path)
pprint.pprint(datetime.datetime.now(), open(file, 'a'))
pprint.pprint("step {first} to {second}".format(first=ci.step1, second=ci.step2) , open(file, 'a'))
pprint.pprint("print step, id, and ifcontact. that ifcontact changed from last step to current step", open(file, 'a'))
np.set_printoptions(threshold = 1000000)
pprint.pprint(step_id_ifover_diffnext, open(file, 'a'))
pprint.pprint('===', open(file, 'a'))