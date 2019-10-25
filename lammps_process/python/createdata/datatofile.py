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

def chunk_hdf5_csv():
    da.handlechunk(override, ifcsv).tohdf5()

class dumptofile:

    def __init__(self, fromtraceorall):
        self.fromtraceorall = fromtraceorall

    def dump_custom(self):
        da.handledumpcustom(override, ifcsv, self.fromtraceorall).tohdf5()

    def dump_custom_select(self, id_i_list):
        da.handlecustomselect(id_i_list, override, ifcsv, self.fromtraceorall).tohdf5()

    def dump_custom_pair_combine(self, id_i_list):
        da.handle_merge_custom_pair(id_i_list, override, ifcsv, self.fromtraceorall).tohdf5()

    def dump_custom_max(self, maxlabel):
        da.handlecustom_max_everysteps(maxlabel, override, ifcsv, self.fromtraceorall).tohdf5()

