# import
import os.path
from io import StringIO
import re
import time
from itertools import chain
from itertools import repeat
from itertools import islice
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D
# import module
import rename_variable as rv
import datapath as dp
# =======================================

# Change sign when switch i j
header_change_sign = ['fjinx', 'fjiny', 'fjinz', 'fjitx', 'fjity', 'fjitz', 'vijnx', 'vijny', 'vijnz', 'vijtx', 'vijty', 'vijtz']
# Switch i j for id and type
header_change_from = ['id_i', 'id_j', 'type_i', 'type_j']
header_change_to = ['id_j', 'id_i', 'type_j', 'type_i']
# Rename variable in pair
rename_mapper_pair = {
	'c_p_idtype_trace[1]':'id_i',
	'c_p_idtype_trace[2]':'id_j',
	'c_p_idtype_trace[3]':'type_i',
	'c_p_idtype_trace[4]':'type_j',
	'c_p_fv_trace[1]':'distance',
	'c_p_fv_trace[2]':'fjinx',
	'c_p_fv_trace[3]':'fjiny',
	'c_p_fv_trace[4]':'fjinz',
	'c_p_fv_trace[5]':'fjin',
	'c_p_fv_trace[6]':'fjitx',
	'c_p_fv_trace[7]':'fjity',
	'c_p_fv_trace[8]':'fjitz',
	'c_p_fv_trace[9]':'fjit',
	'c_p_fv_trace[10]':'vijnx',
	'c_p_fv_trace[11]':'vijny',
	'c_p_fv_trace[12]':'vijnz',
	'c_p_fv_trace[13]':'vijtx',
	'c_p_fv_trace[14]':'vijty',
	'c_p_fv_trace[15]':'vijtz',
	'c_p_idtype_all[1]':'id_i',
	'c_p_idtype_all[2]':'id_j',
	'c_p_idtype_all[3]':'type_i',
	'c_p_idtype_all[4]':'type_j',
	'c_p_fv_all[1]':'distance',
	'c_p_fv_all[2]':'fjinx',
	'c_p_fv_all[3]':'fjiny',
	'c_p_fv_all[4]':'fjinz',
	'c_p_fv_all[5]':'fjin',
	'c_p_fv_all[6]':'fjitx',
	'c_p_fv_all[7]':'fjity',
	'c_p_fv_all[8]':'fjitz',
	'c_p_fv_all[9]':'fjit',
	'c_p_fv_all[10]':'vijnx',
	'c_p_fv_all[11]':'vijny',
	'c_p_fv_all[12]':'vijnz',
	'c_p_fv_all[13]':'vijtx',
	'c_p_fv_all[14]':'vijty',
	'c_p_fv_all[15]':'vijtz',
	}
# ====================================== thermo =========================

def save_dataframe(df, h5filename, ifcsv):

    df.to_hdf(h5filename, key='df', mode='w', format='fixed')
    print ("h5 saved")
    if ifcsv == "yes":
        df.to_csv(h5filename[:-3]+".csv", encoding='utf-8')
        print ("csv saved")
    else:
        pass

def thermofile_to_dataframe(file):
    with open(file) as f:
        lines = f.read().strip().split('\n')
    # index of first row of every block
    id_headers = [n for n, line in enumerate(lines) if line.startswith('Step Atoms')]
    if len(id_headers) > 1:
        string = 'log file has {n} sections, only extract the last section'.format(n=len(id_headers))
        print(string)
    # choose last section of thermo output
    id_header = id_headers[-1]
    # index of data end
    id_ends = [n for n, line in enumerate(lines) if line.startswith('Loop time of')]
    # find data end
    if len(id_ends) < len(id_headers):
        id_end = len(lines)-3
        laststep = lines[id_end].split()[0]
        print ("simulation not complete, thermo only run to step {laststep}, use the 3rd-last line in log file as end of data".format(laststep=laststep))
    else:
        id_end = id_ends[0]
    # header
    header = lines[id_header].split()[0:]
    # select data
    id_data_begin = id_header + 1
    data = [lines[t].split() for t in range(id_data_begin, id_end)]
    # attach data
    df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
    return df

class handlelammpfile(object):
    def __init__(self, f_out, override='no', ifcsv='no', fromtraceorall='trace'):
        self.f_out = f_out
        self.override = override
        self.ifcsv = ifcsv
        self.ifoutputexist = os.path.isfile(self.f_out)

    def todataframe(self):
        pass
    
    def tohdf5(self):
        if self.ifoutputexist and (self.override == 'no'):
            print ('h5 file exist, so not create again')
        else:
            print ("creating h5")
            print ("extract dataframe")
            save_dataframe(self.todataframe(), self.f_out, self.ifcsv)

class handlelog(handlelammpfile):

    def __init__(self, override, ifcsv):
        f_out = dp.f_thermo
        super().__init__(f_out, override, ifcsv)
        self.f_log_in = dp.thermo_path
        print("handle log")

    def todataframe(self):

        if os.path.isfile(self.f_out) and (self.override == 'no'):
            print ('h5 file exist, so create dataframe from previous saved file')
            df = pd.read_hdf(self.f_out)
        else:
            print ("reading f_log_in and creating h5")

            df = thermofile_to_dataframe(self.f_log_in)
        return df

    def tohdf5(self):
        super().tohdf5()

class handlechunk(handlelammpfile):

    def __init__(self, override, ifcsv):
        f_out = dp.f_chunk
        super().__init__(f_out, override, ifcsv)
        self.f_chunk_in = dp.chunk_path
        print("handle chunk")

    def todataframe(self):

        if os.path.isfile(self.f_out) and (self.override == 'no'):
            print ('h5 file exist, so create dataframe from previous saved file')
            df = pd.read_hdf(self.f_out)
        else:
            print ("reading f_chunk_in and creating h5")

            df = chunkfile_to_dataframe(self.f_chunk_in)
        return df

    def tohdf5(self):
        super().tohdf5()

class handledump(handlelammpfile):

    def __init__(self, f_in, f_out, override, ifcsv, fromtraceorall):

        super().__init__(f_out, override, ifcsv, fromtraceorall)
        self.f_in = f_in
        print("handle dump")

    def todataframe(self):

        if os.path.isfile(self.f_out) and (self.override == 'no'):
            print ('h5 file exist, so create dataframe from previous saved file')
            df = pd.read_hdf(self.f_out)
        else:
            print ("reading f_in and creating h5")

            with open(self.f_in) as f:
                # if first string in the first line is digit, than see the whole file as a csv without any string
                firststring = f.readline().split()[0]
            
                if firststring.isdigit():
                    with open(dp.custom_all_firststep_path) as f_first_for_header:
                        # header
                        for i in range(9):
                            header = f_first_for_header.readline().split()[2:]
                    header.insert(0, "step")
                    df = pd.read_csv(self.f_in, names=header, delim_whitespace=True)

                else:
                    f.seek(0)
                    # header
                    for i in range(9):
                        header = f.readline().split()[2:]
                    n_col = len(header)

                    # timestep
                    pattern_TIMESTEP = re.compile("TIMESTEP\n(.*?)IT", re.DOTALL)
                    f.seek(0)
                    steps = pattern_TIMESTEP.findall(f.read())
                    steps = "".join(steps)
                    steps = np.fromstring(steps, dtype=np.float64 ,sep=' ')

                    # data
                    # (?:IT|\Z) means end with IT or fileend but not capture the group in this parentheses
                    pattern_data = re.compile(header[-1] + " \n(.*?)(?:IT|\Z)", re.DOTALL)
                    f.seek(0)
                    data = pattern_data.findall(f.read())

                    # check len data and steps
                    if (len(data)!=steps.shape[0]):
                        sys.exit('len(data) not equal len of steps')
                    
                    # count number of row in each timestep
                    count_rows = np.asarray([onestepdata.count('\n') for onestepdata in data])
                    steps = np.repeat(steps, count_rows)

                    # data string to numpy array
                    data_combine = "".join(data)
                    data_array = np.fromstring(data_combine, dtype=np.float64 ,sep=' ')
                    data_array = data_array.reshape(-1, n_col)

                    # create dataframe
                    df = pd.DataFrame(data = data_array, columns = header, dtype = np.float64)
                    
                    # insert timesteps
                    df.insert(0, 'step', steps)

                #qqq = re.search(r'^ITEM\:\s([^\n]+)', f_read, re.MULTILINE).group(1)
            #    #ccci = [b.start(0) for b in aaa]
            #    #ccc = [b.group(0) for b in aaa]
            #    
            #
            #    timeend = time.time()
            #    print(timeend - timestart)
            #    breakpoint()
            #
            #    timestart = time.time()
            #    id_ITEMTIMESTEPs = []
            #    for n, line in enumerate(f):
            #        if line.startswith('ITEM: TIMESTEP'):
            #            id_ITEMTIMESTEPs.append(n)
            #    timeend = time.time()
            #    print(timeend - timestart)
            #    breakpoint()
            #    f.seek(0)
            #    num_lines = sum(1 for i, line in enumerate(f))
            #    breakpoint()
            #    f.seek(0)
            #    # index of row timestep, start, end
            #    id_ITEMTIMESTEPs.append(num_lines)
            #    id_ITEMTIMESTEPs = np.asarray(id_ITEMTIMESTEPs, dtype=int)
            #    id_steps = id_ITEMTIMESTEPs+1
            #    steps = [line for t, line in enumerate(f) if t in id_steps]  #for i,elm in islice(enumerate(some_list),7,40):
            #    f.seek(0)
            #    steps = np.asarray(steps,dtype=int)
            #    # repeat timestep
            #    num_repeats = id_ITEMTIMESTEPs[1:] - (id_ITEMTIMESTEPs[:-1] + 9)
            #    steps = np.repeat(steps, num_repeats)
            #    
            #    # header
            #    for i in range(id_ITEMTIMESTEPs[0] + 8):
            #        header = f.readline().split()[2:]
            #    f.seek(0)
            #    # combine ranges of blocks
            #    iter = chain.from_iterable(range(id+9, id_ITEMTIMESTEPs[i + 1]) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1]))
            #    # select data
            #    data = [line.split() for t, line in enumerate(f) if t in iter]
            #
            ## attach data
            #df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            ## insert timesteps
            #df.insert(1, 'step', steps)
            
            
            
            #with open(file) as f:
            #    lines = f.read().strip().split('\n')
            ## index of first row of every block	
            #id_ITEMTIMESTEPs = [n for n, line in enumerate(lines) if line.startswith('ITEM: TIMESTEP')]
            ## header
            #header = lines[id_ITEMTIMESTEPs[0] + 8].split()[2:]
            ## index of row timestep, start, end
            #id_ITEMTIMESTEPs.append(len(lines))
            ## combine ranges of blocks
            #iter = chain.from_iterable(range(id+9, id_ITEMTIMESTEPs[i + 1]) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1]))
            ## select data
            #data = [lines[t].split() for t in iter]
            ## attach data
            #df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
            ## repeat timestep
            #steps = list(chain.from_iterable(repeat(lines[id + 1], id_ITEMTIMESTEPs[i + 1] - (id + 9)) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1])))
            #steps = np.asarray(steps,dtype=np.float64)
            ## insert timesteps
            #df.insert(1, 'step', steps)
        return df

    def tohdf5(self):
        super().tohdf5()

class handledumpcustom(handledump):

    def __init__(self, override, ifcsv, fromtraceorall):
        f_out = dp.f_custom
        if fromtraceorall == 'trace':
            f_in = dp.custom_near_trace_path
        elif fromtraceorall == 'all':
            f_in = dp.custom_all_path
        else:
            sys.exit("fromtraceorall not trace not all, please select one of them")
        super().__init__(f_in, f_out, override, ifcsv, fromtraceorall)
        print("handle dumpcustom")

class handledumppair(handledump):

    def __init__(self, override, ifcsv, fromtraceorall):
        f_out = dp.f_pair
        f_in = dp.pair_all_path
        super().__init__(f_in, f_out, override, ifcsv, fromtraceorall)
        print("handle dumppair")


class handlecustomselect(handlelammpfile):

    def __init__(self, id_i_list, override, ifcsv, fromtraceorall):
        f_out = dp.put_id_on_file(id_i_list, dp.f_custom)
        super().__init__(f_out, override, ifcsv, fromtraceorall)
        self.id_i_list = id_i_list
        self.f_custom_all_out = dp.f_custom
        self.dumpcustomall = handledumpcustom(override, ifcsv, fromtraceorall)
        print("handle customselect")
   

    def todataframe(self):

        if os.path.isfile(self.f_out) and (self.override == 'no'):
            print ('h5 file exist, so create dataframe from previous saved file')
            dfc_select = pd.read_hdf(self.f_out)
        else:
            print ("reading f_in_custom and creating h5")
        
            self.dumpcustomall.tohdf5()
            dfc = pd.read_hdf(self.f_custom_all_out)
            dfc_select = select_custom(dfc, self.id_i_list)

        return dfc_select

    def tohdf5(self):
        super().tohdf5()

class handlecustom_max_everysteps(handlelammpfile):

    def __init__(self, maxlabel, override, ifcsv, fromtraceorall):
        f_out = dp.put_maxlabel_on_file(maxlabel, dp.f_custom)
        super().__init__(f_out, override, ifcsv, fromtraceorall)
        self.maxlabel = maxlabel
        self.f_custom_all_out = dp.f_custom
        self.dumpcustomall = handledumpcustom(override, ifcsv, fromtraceorall)
        print("handle custom_max_everysteps")

    def todataframe(self):

        if os.path.isfile(self.f_out) and (self.override == 'no'):
            print ('h5 file exist, so create dataframe from previous saved file')
            df_max = pd.read_hdf(self.f_out)
        else:
            print ("reading f_in_custom and creating h5")
            self.dumpcustomall.tohdf5()
            dfc = pd.read_hdf(self.f_custom_all_out)
            # find maxKE in everysteps
            df_max = dfc.loc[dfc.groupby(["step"])[self.maxlabel].idxmax()]

        return df_max

    def tohdf5(self):
        super().tohdf5()


class handle_merge_custom_pair(handlelammpfile):

    def __init__(self, id_i_list, override, ifcsv, fromtraceorall):
        f_out = dp.put_id_on_file(id_i_list, dp.f_cipcj)
        super().__init__(f_out, override, ifcsv, fromtraceorall)
        self.id_i_list = id_i_list
        self.f_custom_all_out = dp.f_custom
        self.f_pair_all_out = dp.f_pair
        self.dumpcustomall = handledumpcustom(override, ifcsv, fromtraceorall)
        self.dumpcustomselect = handlecustomselect(id_i_list, override, ifcsv, fromtraceorall)
        self.dumppairall = handledumppair(override, ifcsv, fromtraceorall)
        print("handle merge_custom_pair")


    def todataframe(self):

        if os.path.isfile(self.f_out) and (self.override == 'no'):
            print ('h5 file exist, so create dataframe from previous saved file')
            df_select_cipcj = pd.read_hdf(self.f_out)
        else:
            print ("reading f_in_custom f_in_pair and creating h5")

            self.dumpcustomall.tohdf5()
            dfc = pd.read_hdf(self.f_custom_all_out)

            self.dumpcustomselect.tohdf5()
            dfc_select = pd.read_hdf(self.dumpcustomselect.f_out)

            self.dumppairall.tohdf5()
            dfp = pd.read_hdf(self.f_pair_all_out)
            dfp = dfp.rename(columns=rename_mapper_pair)

            
            dfp2 = pair_double(dfp, header_change_from, header_change_to, header_change_sign)
            df_select_cipcj = merge_ci_p_cj(dfc_select, dfp2, dfc)
        
        return df_select_cipcj

    def tohdf5(self):
        super().tohdf5()


def thermo_file_to_h5_csv(f_thermo, f_output_name, override ='no'):

    if os.path.isfile(f_output_name + '.h5') and (override != 'yes'):
        print ('h5 file exist, so not create again')
    else:
        df = thermofile_to_dataframe(f_thermo)
        df.to_hdf(f_output_name + '.h5', key='df', mode='w')
        df.to_csv(f_output_name + '.csv', encoding='utf-8')


# ====================================== pari custom =====================

# define function for extract data from dump txt to dataframe
def dumpfile_to_dataframe(file, selectid='all'):
	
    with open(file) as f:

        # header
        for i in range(9):
            header = f.readline().split()[2:]
        n_col = len(header)

        # timestep
        pattern_TIMESTEP = re.compile("TIMESTEP\n(.*?)IT", re.DOTALL)
        f.seek(0)
        steps = pattern_TIMESTEP.findall(f.read())
        steps = "".join(steps)
        steps = np.fromstring(steps, dtype=np.float64 ,sep=' ')

        # data
        # (?:IT|\Z) means end with IT or fileend but not capture the group in this parentheses
        pattern_data = re.compile(header[-1] + " \n(.*?)(?:IT|\Z)", re.DOTALL)
        f.seek(0)
        data = pattern_data.findall(f.read())

    # check len data and steps
    if (len(data)!=steps.shape[0]):
        sys.exit('len(data) not equal len of steps')
    
    # count number of row in each timestep
    count_rows = np.asarray([onestepdata.count('\n') for onestepdata in data])
    steps = np.repeat(steps, count_rows)

    # data string to numpy array
    data_combine = "".join(data)
    data_array = np.fromstring(data_combine, dtype=np.float64 ,sep=' ')
    data_array = data_array.reshape(-1, n_col)

    # create dataframe
    df = pd.DataFrame(data = data_array, columns = header, dtype = np.float64)
    
    # insert timesteps
    df.insert(0, 'step', steps)

        #qqq = re.search(r'^ITEM\:\s([^\n]+)', f_read, re.MULTILINE).group(1)
    #    #ccci = [b.start(0) for b in aaa]
    #    #ccc = [b.group(0) for b in aaa]
    #    
    #
    #    timeend = time.time()
    #    print(timeend - timestart)
    #    breakpoint()
    #
    #    timestart = time.time()
    #    id_ITEMTIMESTEPs = []
    #    for n, line in enumerate(f):
    #        if line.startswith('ITEM: TIMESTEP'):
    #            id_ITEMTIMESTEPs.append(n)
    #    timeend = time.time()
    #    print(timeend - timestart)
    #    breakpoint()
    #    f.seek(0)
    #    num_lines = sum(1 for i, line in enumerate(f))
    #    breakpoint()
    #    f.seek(0)
    #    # index of row timestep, start, end
    #    id_ITEMTIMESTEPs.append(num_lines)
    #    id_ITEMTIMESTEPs = np.asarray(id_ITEMTIMESTEPs, dtype=int)
    #    id_steps = id_ITEMTIMESTEPs+1
    #    steps = [line for t, line in enumerate(f) if t in id_steps]  #for i,elm in islice(enumerate(some_list),7,40):
    #    f.seek(0)
    #    steps = np.asarray(steps,dtype=int)
    #    # repeat timestep
    #    num_repeats = id_ITEMTIMESTEPs[1:] - (id_ITEMTIMESTEPs[:-1] + 9)
    #    steps = np.repeat(steps, num_repeats)
    #    
    #    # header
    #    for i in range(id_ITEMTIMESTEPs[0] + 8):
    #        header = f.readline().split()[2:]
    #    f.seek(0)
    #    # combine ranges of blocks
    #    iter = chain.from_iterable(range(id+9, id_ITEMTIMESTEPs[i + 1]) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1]))
    #    # select data
    #    data = [line.split() for t, line in enumerate(f) if t in iter]
    #
    ## attach data
    #df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
    ## insert timesteps
    #df.insert(1, 'step', steps)
    
    
    
    #with open(file) as f:
    #    lines = f.read().strip().split('\n')
    ## index of first row of every block	
    #id_ITEMTIMESTEPs = [n for n, line in enumerate(lines) if line.startswith('ITEM: TIMESTEP')]
    ## header
    #header = lines[id_ITEMTIMESTEPs[0] + 8].split()[2:]
    ## index of row timestep, start, end
    #id_ITEMTIMESTEPs.append(len(lines))
    ## combine ranges of blocks
    #iter = chain.from_iterable(range(id+9, id_ITEMTIMESTEPs[i + 1]) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1]))
    ## select data
    #data = [lines[t].split() for t in iter]
    ## attach data
    #df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
    ## repeat timestep
    #steps = list(chain.from_iterable(repeat(lines[id + 1], id_ITEMTIMESTEPs[i + 1] - (id + 9)) for i, id in enumerate(id_ITEMTIMESTEPs[0: -1])))
    #steps = np.asarray(steps,dtype=np.float64)
    ## insert timesteps
    #df.insert(1, 'step', steps)

    return df

# extract custom for specific id
def select_custom(df_custom, select_id_list):

    if select_id_list == 'all':
        return df_custom
    else:
        # number of step of df_custom
        steps_of_df_custom = df_custom['step'].values
        firststep = steps_of_df_custom[0]
        secondstep = steps_of_df_custom[steps_of_df_custom!=firststep][0]
        d_step = secondstep - firststep
        laststep = steps_of_df_custom[-1]
        n_step_df_custom = (laststep - firststep)/d_step + 1
        # number of select id
        n_select = len(select_id_list)
        # select traceid for dataframe
        df_select_custom = df_custom.loc[df_custom['id'].isin(select_id_list)]
        # number of row of df_select_custom
        n_row_df_select_custom = len(df_select_custom.index)
        # check len of step array is n_select*n_step_df_custom
        if n_row_df_select_custom < n_select*n_step_df_custom:
            sys.exit('len of steps selected less than number_of_atom*number of step, maybe one of the select atom not appear in all the steps of df_custom')
        elif n_row_df_select_custom > n_select*n_step_df_custom:
            sys.exit('len of steps selected larger than expected')
        else:
            pass

        steps_array_expected = np.repeat(np.arange(firststep, laststep+1, d_step), n_select)
        steps_df_select_custom = df_select_custom['step'].values
        check_steps = (steps_array_expected == steps_df_select_custom)
        if not np.all(check_steps):
            id_step_firstwrong = np.nonzero(~check_steps)[0][0]
            step_firstwrong = steps_array_expected[id_step_firstwrong]
            sys.exit('step ' + repr(step_firstwrong) + 'has less or more atom than expected')
        else:
            pass
        
        return df_select_custom

# create ji pair from ij pair and combine together
def pair_double(df_pair, header_change_from, header_change_to, header_change_sign):
    # copy
    df_switch = df_pair.copy()
    # change sign
    df_switch[header_change_sign] = -df_switch[header_change_sign].values
    # change header names
    mapper = {v: header_change_to[i] for i, v in enumerate(header_change_from)}
    df_switch = df_switch.rename(columns=mapper)
    # reorder columns
    df_switch = df_switch[df_pair.columns.values.tolist()]
    # combine df_pair and df_switch
    dfp_data_double = np.empty((2*df_pair.shape[0], df_pair.shape[1]))
    dfp_data_double[::2] = df_pair.values
    dfp_data_double[1::2] = df_switch.values
    df_double = pd.DataFrame(data = dfp_data_double, columns = df_pair.columns.values.tolist())
    return df_double

# merge pairdouble and select_custom, output dfcip
def merge_ci_p(dfc, dfp2):
    # add suffix on custom
    dfc = dfc.add_suffix('_i')
    dfc = dfc.rename(columns={'step_i': 'step'})
    dfcip = pd.merge(dfc, dfp2, how='left', on=['step', 'id_i'])
    # check if type_i_x = type_i_y except nan (means no pair for the timestep). If yes then delete type_i_y
    type_notconsistent = np.where(dfcip['type_i_x'].values != dfcip['type_i_y'].values)[0]
    where_type_i_y_nan = np.where(np.isnan(dfcip['type_i_y'].values))[0]
    check_only_nan_inconsistent = np.all(type_notconsistent == where_type_i_y_nan)
    if ~check_only_nan_inconsistent:
        sys.exit('type_i not consistent when merge custom and pair')
    else:
        # delete rename type_i_x, delete type_i_y
        dfcip = dfcip.drop('type_i_y', axis=1)
        dfcip = dfcip.rename(columns={'type_i_x': 'type_i'})
    return dfcip

# merge dfcip and custom_j, output dfcipcj
def merge_cip_cj(dfcip, dfc):
    # add suffix on custom
    dfc = dfc.add_suffix('_j')
    dfc = dfc.rename(columns={'step_j': 'step'})
    dfcipcj = pd.merge(dfcip, dfc, how='left', on=['step', 'id_j'])
    # check if type_j_x = type_j_y except nan (means no pair for the timestep). If yes then delete type_j_y
    type_notconsistent = np.where(dfcipcj['type_j_x'].values != dfcipcj['type_j_y'].values)[0]
    where_type_j_y_nan = np.where(np.isnan(dfcipcj['type_j_y'].values))[0]
    check_only_nan_inconsistent = np.all(type_notconsistent == where_type_j_y_nan)
    if ~check_only_nan_inconsistent:
        sys.exit('type_j not consistent when merge custom and pair')
    else:
        # delete rename type_j_x, delete type_j_y
        dfcipcj = dfcipcj.drop('type_j_y', axis=1)
        dfcipcj = dfcipcj.rename(columns={'type_j_x': 'type_j'})
    return dfcipcj

# merge pairdouble, custom_i and custom_j. output dfcipcj 
def merge_ci_p_cj(dfci_select, dfp2, dfc):
    dfcip = merge_ci_p(dfci_select, dfp2)
    dfcipcj = merge_cip_cj(dfcip, dfc)
    return dfcipcj

def file_to_h5_csv(f_dumppair, f_dumpcustom, id_i_list, f_output_name, override ='no'):
    
    if f_dumppair != None:
        # if file exist and override not yes, then not to write files
        if os.path.isfile(f_output_name + '.h5') and (override != 'yes'):
            print ('h5 file exist, so not create again')
        else:
            dfp = dumpfile_to_dataframe(f_dumppair)
            dfp = dfp.rename(columns=rename_mapper_pair)
            dfc = dumpfile_to_dataframe(f_dumpcustom)
            dfc_select = select_custom(dfc, id_i_list)
            dfp2 = pair_double(dfp, header_change_from, header_change_to, header_change_sign)
            df_select_cipcj = merge_ci_p_cj(dfc_select, dfp2, dfc)
            df_select_cipcj.to_hdf(f_output_name + '.h5', key='df', mode='w')
            df_select_cipcj.to_csv(f_output_name + '.csv', encoding='utf-8')

    else:
        if os.path.isfile(f_output_name + '.h5') and (override != 'yes'):
            print ('h5 file exist, so not create again')
        else:
            dfc = dumpfile_to_dataframe(f_dumpcustom)
            dfc_select = select_custom(dfc, id_i_list)
            dfc_select.to_hdf(f_output_name + '.h5', key='df', mode='w')
            dfc_select.to_csv(f_output_name + '.csv', encoding='utf-8')

# ====================================== fix chunk =====================

# define function for extract data from fix txt to dataframe
def chunkfile_to_dataframe(file):

    with open(file) as f:
        lines = f.read().strip().split('\n')
        id_line_timestep = [n for n, line in enumerate(lines) if line.startswith('# Timestep')]
        n_chunks = int(lines[id_line_timestep[0]+2].split()[2])
        header = lines[id_line_timestep[0]+1].split()[1:]
        id_line_timestep.append(len(lines))
        iter = chain.from_iterable(range(id + 3, id_line_timestep[i + 1] - 1) for i, id in enumerate(id_line_timestep[0: -1]))
        ## select data
        data = [lines[t].split() for t in iter]
        ## attach data
        df = pd.DataFrame(data = data, columns = header, dtype = 'float64')
        ## repeat timestep
        steps = list(
            chain.from_iterable(
                repeat(
                    lines[id + 2].split()[0], id_line_timestep[i + 1] - -id - 4) for i, id in enumerate(id_line_timestep[0: -1]
                    )
                )
            )
        steps = np.asarray(steps,dtype=np.float64)
        ## insert timesteps
        df.insert(1, 'step', steps)
        
    return df

