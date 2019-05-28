# import module in simulation folder
import osmanage as om
import datapath as dp
# import module in createdata folder
import createdata.dataarrange as da
# ======================================

om.create_directory(dp.post_process_path)
om.create_directory(dp.hdf5_csv_path)

override='no'
ifcsv='no'

def thermo_hdf5_csv():
    da.handlelog(dp.thermo_path, override, ifcsv).tohdf5()

def dump_custom():
    da.handledumpcustom(dp.custom_path, override, ifcsv).tohdf5()

def dump_custom_select(id_i_list):
    da.handlecustomselect(dp.custom_path, id_i_list, override, ifcsv)

def dump_custom_pair_combine(id_i_list):
    da.handle_merge_custom_pair(dp.custom_path, dp.pair_path, id_i_list, override, ifcsv)
    
    