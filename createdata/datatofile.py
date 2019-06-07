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
    da.handlelog(override, ifcsv).tohdf5()

def dump_custom():
    da.handledumpcustom(override, ifcsv).tohdf5()

def dump_custom_select(id_i_list):
    da.handlecustomselect(id_i_list, override, ifcsv).tohdf5()

def dump_custom_pair_combine(id_i_list):
    da.handle_merge_custom_pair(id_i_list, override, ifcsv).tohdf5()

def dump_custom_max(maxlabel):
    da.handlecustom_max_everysteps(maxlabel, override, ifcsv).tohdf5()

