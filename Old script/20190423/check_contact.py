import pprint
# import project module
import check_simulation as cs
import output_control as oc
import datapath as dp
import data_arrange as da

result = cs.contact_check_multistep_singlecontact(dp.f_cipcj + '.h5', oc.id_i, oc.id_j, oc.step1, oc.step2)
da.create_directory(dp.debug_print_path)
pprint.pprint(result, open(dp.debug_print_path + f_name_check_c, 'a'))