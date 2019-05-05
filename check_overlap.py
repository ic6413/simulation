import pprint
import datetime
import numpy as np
# import project module
import check_simulation as cs
import output_control as oc
import datapath as dp
import data_arrange as da

file = dp.debug_print_path + oc.f_name_check_c

[step_id_ifover_diffnext, ifoverlap_ij_iw_array, initial_overlap_id] = cs.contact_check_multistep(dp.f_custom + '.h5', oc.id_i, oc.step1, oc.step2, oc.error_tolerence)

np.set_printoptions(threshold = 100000)

da.create_directory(dp.debug_print_path)
pprint.pprint(datetime.datetime.now(), open(file, 'a'))
pprint.pprint("step {first} to {second}".format(first=oc.step1, second=oc.step2) , open(file, 'a'))
pprint.pprint("(wall_id, wallname) = {first}".format(first=cs.wall_list_id_name), open(file, 'a'))
pprint.pprint("initial overlap id are {first}".format(first=initial_overlap_id), open(file, 'a'))
pprint.pprint("print step, id, and ifcontact. that ifcontact changed from last step to current step", open(file, 'a'))
pprint.pprint(step_id_ifover_diffnext, open(file, 'a'))
pprint.pprint('===', open(file, 'a'))
