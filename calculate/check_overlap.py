import pprint
import datetime
import numpy as np
# import module in simulation folder
import osmanage as om
import datapath as dp
# import module in calculate folder
import calculate.checksimulation as cs
import calculate.inputvariable as ci

file = dp.debug_print_path + ci.f_inputvariable_name

[step_id_ifover_diffnext, ifoverlap_ij_iw_array, initial_overlap_id] = cs.contact_check_multistep(dp.f_custom + '.h5', ci.id_i, ci.step1, ci.step2, ci.error_tolerence)

np.set_printoptions(threshold = 100000)

om.create_directory(dp.debug_print_path)
pprint.pprint(datetime.datetime.now(), open(file, 'a'))
pprint.pprint("step {first} to {second}".format(first=ci.step1, second=ci.step2) , open(file, 'a'))
pprint.pprint("(wall_id, wallname) = {first}".format(first=cs.wall_list_id_name), open(file, 'a'))
pprint.pprint("initial overlap id are {first}".format(first=initial_overlap_id), open(file, 'a'))
pprint.pprint("print step, id, and ifcontact. that ifcontact changed from last step to current step", open(file, 'a'))
pprint.pprint(step_id_ifover_diffnext, open(file, 'a'))
pprint.pprint('===', open(file, 'a'))
