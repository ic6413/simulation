# import
import os
# function

def put_id_on_file(id_list, f_name_without_id):
    
    if custom_id_i_list != 'all':
        f_name_add_id = (f_name_without_id + '_' + '_'.join(str(i) for i in custom_id_i_list))
    else:
        f_name_add_id = (f_name_without_id + '_' + 'all')

    return f_name_add_id 

#=========== lammmps folder path ====================

folder_path = os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')
#=========== data arrange ========================

# select combine id_i
combine_id_i_list = 'all'
custom_id_i_list = 'all' # id_list or 'all'
# fcipcj name
f_cipcj = put_id_on_file(combine_id_i_list, 'cipcj_id')
# fcustom name
f_custom = put_id_on_file(custom_id_i_list, 'custom_id')

#============ check ft overlap ========================

# select id_i to check ft and overlap
id_i = 15556

# step
[step1, step2] = [173000, 300000]   # 57858 ~ 57861 no contact, 57862 contact wall. wrong at 57861~59780   [57830,57860] [27940, 27962]
f_name_check_f = 'check_force_' + 'id_' + str(id_i) +'step_' + str(step1) + '_' + str(step2)
f_name_check_c = 'check_overlap' + '_' + str(step1) + '_' + str(step2)
# error tolerence
error_tolerence = 1e-6

# overlap tolerence
overlap_tolerence = 0.0

#==================== end ===========================
#======= 1contact ====================

