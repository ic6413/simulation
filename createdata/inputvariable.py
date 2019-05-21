import datapath as dp

# select combine id_i
combine_id_i_list = 'all'
custom_id_i_list = 'all' # id_list or 'all'

def put_id_on_file(id_list, f_name_without_id):
    
    if id_list != 'all':
        f_name_add_id = (f_name_without_id + '_' + '_'.join(str(i) for i in id_list))
    else:
        f_name_add_id = (f_name_without_id + '_' + 'all')

    return f_name_add_id

f_cipcj = dp.hdf5_csv_path + put_id_on_file(combine_id_i_list, 'cipcj_id')
f_custom = dp.hdf5_csv_path + put_id_on_file(custom_id_i_list, 'custom_id')
f_thermo = dp.hdf5_csv_path + 'thermo'