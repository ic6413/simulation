import createdata.inputvariable
method_list = list(range(0, 3))
worj = 'j'
id_i = 15556
step1 = 0
step2 = 930000
error_tolerence = 1e-8
f_custom_path = createdata.inputvariable.f_custom + '.h5'

f_check_f_name = 'check_force_' + 'id_' + str(id_i) +'step_' + str(step1) + '_' + str(step2)
f_check_overlap_name = 'check_overlap_' + 'id_' + str(id_i) +'step_' + str(step1) + '_' + str(step2)
f_check_f_1contact_name = 'check_1contact_force_' + worj + '_id_' + str(id_i) +'step_' + str(step1) + '_' + str(step2)