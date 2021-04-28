# import modules
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import pickle
# import mymodules
import read_setting as rr
import datapath as dp
# chunk method
if rr.log_variable["shearwall"] == "zcylinder":
    chunk_method = 'rz'
if rr.log_variable["shearwall"] == "yplane":
    chunk_method = 'yz'

# map dim index to coordinate
if rr.log_variable["shearwall"] == "zcylinder":
    map_dim_index_to_coordinate = ["t", "r", "z"]
elif rr.log_variable["shearwall"] == "yplane":
    map_dim_index_to_coordinate = ["x", "y", "z"]

# xyztoCoor
if rr.log_variable["shearwall"] == "zcylinder":
    # chink first dim unchange in the begining
    chunk_first_dim_coord = "z"
    chunk_second_dim_coord = "r"
elif rr.log_variable["shearwall"] == "yplane":
    if "chunk/atom 23" in rr.log_variable.keys():
        if rr.log_variable["chunk/atom 23"][1] == "y":
            chunk_first_dim_coord = "y"
            chunk_second_dim_coord = "z"
            xyztoCoor = {}
            xyztoCoor["y"] = "Coord1"
            xyztoCoor["z"] = "Coord2"
        elif rr.log_variable["chunk/atom 23"][1] == "z":
            chunk_first_dim_coord = "z"
            chunk_second_dim_coord = "y"
            xyztoCoor["z"] = "Coord1"
            xyztoCoor["y"] = "Coord2"
        else:
            sys.exit("chunk_method wrong")
    else:
        chunk_first_dim_coord = "y"
        chunk_second_dim_coord = "z"
        xyztoCoor = {}
        xyztoCoor["y"] = "Coord1"
        xyztoCoor["z"] = "Coord2"
else:
    sys.exit("chunk_method wrong")

# n_1 n_2 position_index_to_array_dim_index
if rr.log_variable["shearwall"] == "zcylinder":
    position_index_to_array_dim_index = {
                                    1: 1,
                                    2: 0,
                                    }
    n_r = int(rr.log_variable['N_bin_r'])
    n_z = int(rr.log_variable['N_bin_z'])
    n_1 = n_z
    n_2 = n_r
    n_12 = n_1*n_2
    
    dx = 1/n_r*int(rr.log_variable['width_wall_dp_unit'])
    dy = 1/n_z*float(rr.log_variable['zhi_chunk_dp_unit'])
    x_array, y_array = np.meshgrid(
                                int(rr.log_variable['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*int(rr.log_variable['width_wall_dp_unit']),
                                (np.arange(n_2)+0.5)/n_2*float(rr.log_variable['zhi_chunk_dp_unit']),
                                )
    x_array = x_array.reshape((-1))
    y_array = y_array.reshape((-1))
    vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*float(rr.log_variable['dp'])**3
elif rr.log_variable["shearwall"] == "yplane":
    n_y = int(rr.log_variable['N_bin_y'])
    n_z = int(rr.log_variable['N_bin_z'])
    if "chunk/atom 23" in rr.log_variable.keys():
        if rr.log_variable["chunk/atom 23"][1] == "y":
            position_index_to_array_dim_index = {
                                            1: 0,
                                            2: 1,
                                            }
            n_1 = n_y
            n_2 = n_z
        elif rr.log_variable["chunk/atom 23"][1] == "z":
            position_index_to_array_dim_index = {
                                            2: 0,
                                            1: 1,
                                            }
            n_1 = n_z
            n_2 = n_y
        else:
            sys.exit("chunk_method wrong")
    else:
        position_index_to_array_dim_index = {
                                        1: 0,
                                        2: 1,
                                        }
        n_1 = n_y
        n_2 = n_z

    n_12 = n_1*n_2
    dx = 1/n_y*int(rr.log_variable['width_wall_dp_unit'])
    dy = 1/n_z*float(rr.log_variable['zhi_chunk_dp_unit'])
    vol_in_chunks = float(rr.log_variable['dp'])*int(rr.log_variable['x_period_dp_unit'])*dx*dy*float(rr.log_variable['dp'])**2
else:
    sys.exit("chunk_method wrong")
# time count from step 0
def time_from_step_0(step):
    return step*float(rr.log_variable["ts"])

# time count from rotate started
def time_from_start_rotate(step):
    return time_from_step_0(step)-rr.log_variable["rotate_start_time"]

def extend_chunk_object_to_multisimu(class_name, *arg):
    ob1 = class_name(*arg)
    
    # check if header same
    # if header not the same, output simu_index, and do not append line
    # if header the same, combine lines
    def info_dic_for_all_simu_has_sameheader(self):
        index_list_same_header = []
        for i in rr.n_simu_total - 1:
            header = chunk.infofrom_filepath_n_ave_dic(
                rr.folder_path_list_initial_to_last[i] + self.f_path_rela_lmpfolder,
                self.n_ave,
            )
            if header == self.header:
                index_list_same_header.append(i)
        return index_list_same_header
    
    @staticmethod
    def separate_steps_to_different_simulation(stepslist):
        stepslist.sort()
        stepslist_in_each_simu_list = [] # [first_simu_steps_list, second_simu_steps_list.....]
        stepslist_in_each_simu = []
        k = 0
        if rr.n_simu_total >= 2:
            for i in range(rr.n_simu_total):
                for step in stepslist[k:]:
                    if i == (rr.n_simu_total - 1):
                        stepslist_in_each_simu_list.append(stepslist[k:])
                    else:
                        if step > int(rr.log_variable_dic_list[i]["rst_from"]) and step <= int(rr.log_variable_dic_list[i+1]["rst_from"]):
                            stepslist_in_each_simu.append(step)
                            k = k+1
                        elif step > int(rr.log_variable_dic_list[i+1]["rst_from"]):
                            stepslist_in_each_simu_list.append(stepslist_in_each_simu)
                            stepslist_in_each_simu = []
                            break
                        else:
                            breakpoint()
                            sys.exit("step < rst_from")
        return stepslist_in_each_simu_list



# class chunk common
class chunk(object):

    if "if_inwall_wall_gran" in rr.log_variable.keys():
        if rr.log_variable["if_inwall_wall_gran"] == "yes":
            if "wall_gran_type" in rr.log_variable.keys():
                if rr.log_variable["wall_gran_type"] == "1":
                    ybottomwalltype = "rough (d=0.9)"
                elif rr.log_variable["wall_gran_type"] == "2":
                    ybottomwalltype = "rough (d=1)"
                elif rr.log_variable["wall_gran_type"] == "3":
                    ybottomwalltype = "rough (d=1.1)"
                else:
                    sys.exit("can not get wall gran type")
            else:
                ybottomwalltype = "rough (d=1)"
        else:
            ybottomwalltype = "smooth"
    else:
        ybottomwalltype = "smooth"

    height = rr.log_variable["z_length_create_dp_unit"]
    width = rr.log_variable["width_wall_dp_unit"]
    periodlength = rr.log_variable["x_period_dp_unit"]
    labelstring_size_walltype = ybottomwalltype + "\n" + "L " + periodlength + "\n" + "W " + width + "\n" + "H " + height
    labelstring_size_walltype_one_line = ybottomwalltype + ", " + "L " + periodlength + ", " + "W " + width + ", " + "H " + height

    def __init__(self, n_ave, lmp_path, f_path_rela_lmpfolder):
        self.n_ave = n_ave
        self.lmp_path = lmp_path
        self.f_path_rela_lmpfolder = f_path_rela_lmpfolder
        self.infodic_fromfile = chunk.infofrom_filepath_n_ave_dic(lmp_path + f_path_rela_lmpfolder, n_ave)
        # infofromlines dictionary
        self.lines = self.infodic_fromfile["lines"]
        self.header = self.infodic_fromfile["header"]
        self.n_line_in_a_step = self.infodic_fromfile["n_line_in_a_step"]
        self.step_first_in_file = self.infodic_fromfile["step_first_in_file"]
        self.step_second_in_file = self.infodic_fromfile["step_second_in_file"]
        self.step_last_in_file = self.infodic_fromfile["step_last_in_file"]
        self.d_step = self.infodic_fromfile["d_step"]
        self.step_first_in_file_change_by_n_ave = self.infodic_fromfile["step_first_in_file_change_by_n_ave"]
        self.step_last_in_file_change_by_n_ave = self.infodic_fromfile["step_last_in_file_change_by_n_ave"]
        self.middle_step = self.infodic_fromfile["middle_step"]
        self.allsteps = self.infodic_fromfile["allsteps"]
        self.first_middle_last_steps = self.infodic_fromfile["first_middle_last_steps"]
        self.extrasteps = self.infodic_fromfile["extrasteps"]
        self.first_extra_middle_last_steps = self.infodic_fromfile["first_extra_middle_last_steps"]
        self.index_wasread_from = rr.n_simu_total-1
    # get info dictionary from filepath and n_ave
    @staticmethod
    def infofrom_filepath_n_ave_dic(filepath, n_ave):
        # get lines from files
        with open(filepath) as f:
            lines = f.read().strip().split('\n')
        header = lines[2].split()[1:]
        n_line_in_a_step = int(lines[3].split()[1])
        step_first_in_file = int(lines[3].split()[0])
        step_second_in_file = int(lines[3 + n_line_in_a_step + 1].split()[0])
        step_last_in_file = int(lines[-1 - n_line_in_a_step].split()[0])
        d_step = step_second_in_file - step_first_in_file
        step_first_in_file_change_by_n_ave = step_first_in_file + (n_ave-1)/2*d_step
        step_last_in_file_change_by_n_ave = step_last_in_file - (n_ave-1)/2*d_step
        if step_first_in_file_change_by_n_ave > step_last_in_file_change_by_n_ave:
            sys.exit("error: step_first_in_file_change_by_n_ave > step_last_in_file_change_by_n_ave, n_ave too large") 
        middle_step = int(
            int((step_last_in_file_change_by_n_ave - step_first_in_file_change_by_n_ave)/2/d_step)*d_step + step_first_in_file_change_by_n_ave
            )
        allsteps = np.arange(step_first_in_file_change_by_n_ave, step_last_in_file_change_by_n_ave, d_step)
        first_middle_last_steps = np.array([step_first_in_file_change_by_n_ave, middle_step, step_last_in_file_change_by_n_ave])
        extrasteps = np.array(
            [
                step_first_in_file_change_by_n_ave + 5*d_step*n_ave,
                ]
            )
        maskextra = np.logical_and(extrasteps > step_first_in_file_change_by_n_ave, extrasteps < step_last_in_file_change_by_n_ave)
        extrasteps = extrasteps[maskextra]
        first_extra_middle_last_steps = np.append(first_middle_last_steps, extrasteps)
        first_extra_middle_last_steps.sort()
        infofrom_filepath_n_ave_dic = {
            "lines": lines,
            "header": header,
            "n_line_in_a_step": n_line_in_a_step,
            "step_first_in_file": step_first_in_file,
            "step_second_in_file": step_second_in_file,
            "step_last_in_file": step_last_in_file,
            "d_step": d_step,
            "step_first_in_file_change_by_n_ave": step_first_in_file_change_by_n_ave,
            "step_last_in_file_change_by_n_ave": step_last_in_file_change_by_n_ave,
            "middle_step": middle_step,
            "allsteps": allsteps,
            "first_middle_last_steps": first_middle_last_steps,
            "extrasteps": extrasteps,
            "first_extra_middle_last_steps": first_extra_middle_last_steps,
        }
        return infofrom_filepath_n_ave_dic
    
    def infofrom_filepath_n_ave_dic_by_stepsarray(self, stepsarray):
        
        # get smallest steps in stepsarray
        min_step = np.min(stepsarray)
        # if min_step <= rst_from then read previous

        # check index_read_from
        for index_read_from in range(rr.n_simu_total):
            if min_step > int(rr.log_variable_dic_list[index_read_from]["rst_from"]):
                break
            else:
                pass
        if index_read_from < self.index_wasread_from:
            for j in range(self.index_wasread_from-1, index_read_from, -1):
                # infofromlines dictionary
                infodic_fromfile = chunk.infofrom_filepath_n_ave_dic(rr.folder_path_list_initial_to_last[j] + self.f_path_rela_lmpfolder, self.n_ave)
                # revise lines and corelatted attribute from infofromlines dictionary
                # check if header same
                if self.header != infodic_fromfile["header"]:
                    breakpoint()
                    sys.exit("header not same ")
                elif self.n_line_in_a_step != infodic_fromfile["n_line_in_a_step"]:
                    sys.exit("n_line_in_a_step not same")
                elif self.d_step != infodic_fromfile["d_step"]:
                    sys.exit("d_step not same")
                else:
                    self.lines = infodic_fromfile["lines"] + self.lines[3:]
                    self.step_first_in_file = infodic_fromfile["step_first_in_file"]
                    self.step_second_in_file = infodic_fromfile["step_second_in_file"]
                    self.step_last_in_file = infodic_fromfile["step_last_in_file"]
                    self.step_first_in_file_change_by_n_ave = infodic_fromfile["step_first_in_file_change_by_n_ave"]
                    self.step_last_in_file_change_by_n_ave = infodic_fromfile["step_last_in_file_change_by_n_ave"]
                    self.middle_step = infodic_fromfile["middle_step"]
                    self.allsteps = np.arange(self.step_first_in_file_change_by_n_ave, self.step_last_in_file_change_by_n_ave, self.d_step)
                    self.first_middle_last_steps = np.append(infodic_fromfile["first_middle_last_steps"], self.first_middle_last_steps)
                    self.extrasteps = infodic_fromfile["extrasteps"]
                    self.first_extra_middle_last_steps = np.append(infodic_fromfile["first_extra_middle_last_steps"], self.first_extra_middle_last_steps)
                    # reset self.index_wasread_from
            self.index_wasread_from = index_read_from
        
    
    def get_coord(self, variable_name):
        pass
    def firstdata(self, variable_name):
        pass
    def value_in_a_step_ave(self, step, variable_name, index_n_1=None, index_n_2=None, ddof=1):
        pass
    def collectdataforplot(self, stepsarray, variable_name, index_n_1=None, index_n_2=None, ddof=1):
        pass
    # n_ave path
    def path_nve_subfolder_in_folder(self, folder):
        subfolder = folder + "nve_" + str(self.n_ave) + "/"
        return subfolder
    # add n_ave subfolder
    def add_nve_subfolder_in_folder(self, folder):
        if not os.path.isdir(folder): 
            os.mkdir(folder)
        if not os.path.isdir(self.path_nve_subfolder_in_folder(folder)): 
            os.mkdir(self.path_nve_subfolder_in_folder(folder))

    @staticmethod
    def save_one_plot(fig, ax, savepath, f_name, figformat="png", ifpickle=False, bbox_inches = None):
        # using bbox_inches = 'tight' to ensure legend in plot and outside axe
        fig.savefig(savepath + f_name + "." + figformat, format=figformat, bbox_inches = bbox_inches)
        if ifpickle:
            # Save figure handle to disk
            with open(savepath + f_name + ".pickle", 'wb') as f: # should be 'wb' rather than 'w'
                pickle.dump(fig, f)
        plt.close('all')

    @staticmethod
    def plot_quiver_position_label(fig, ax):
    
        if rr.log_variable["shearwall"] == "zcylinder":
            plt.xlabel('r')
            plt.ylabel('z')
        elif rr.log_variable["shearwall"] == "yplane":
            plt.xlabel('y')
            plt.ylabel('z')
        
        return (fig, ax)

class chunk1D(chunk):
    def __init__(self, n_ave, lmp_path, f_path_rela_lmpfolder, coodfilepath_rela_lmpfolder):
        super().__init__(n_ave, lmp_path, f_path_rela_lmpfolder)
        self.coodfilepath_rela_lmpfolder = coodfilepath_rela_lmpfolder
        # coordinate in chunk
        try:
            self.chunk_coor_array = self.get_coord("Coord2")
        except:
            # old version chunk
            if rr.log_variable["shearwall"] == "zcylinder":
                self.chunk_coor_array = self.firstdata("v_z")
            else:
                sys.exit("chunk_method wrong")
        self.chunk_coor_array_label = "z"
        
    def get_coord(self, variable_name):
        if variable_name in self.header:
            coordanswer = self.firstdata(variable_name)
        else:
            with open(self.lmp_path + self.coodfilepath_rela_lmpfolder) as f:
                lines = f.read().strip().split('\n')
            n_line_0 = 4
            n_line_1 = n_line_0 + self.n_line_in_a_step
            header = lines[2].split()[1:]
            ## select data
            data = [lines[t].split() for t in range(n_line_0, n_line_1)]
            coordanswer = np.asarray(data, dtype=np.float64, order='F')[:, header.index(variable_name)]
        return coordanswer

    def firstdata(self, variable_name):
        n_line_0 = 4
        n_line_1 = n_line_0 + self.n_line_in_a_step
        ## select data
        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
        return (
            np.asarray(data, dtype=np.float64, order='F')[:, self.header.index(variable_name)]
        )

    def value_in_a_step_ave(self, step, variable_name, index_n=None, ddof=1):
        # mean value
        value = 0
        # if lmp output std
        sq_variable_name = variable_name + "_sq"
        if sq_variable_name in self.header:
            ave_square_2 = 0
            if index_n == None:
                for i in range(self.n_ave):
                    oldvalue = value
                    step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                    n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                    n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                    ## select data
                    data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                    data = np.asarray(data, dtype=np.float64, order='F')
                    realchunk = (data[:, self.header.index('Ncount')] > 10**-9)
                    x = data[:, self.header.index(variable_name)]
                    value = value + x*realchunk
                    # get std from lmp output
                    ave_square_2_plus = data[:, self.header.index(sq_variable_name)]
                    ave_square_2 = ave_square_2 + ave_square_2_plus*realchunk
            else:
                for i in range(self.n_ave):
                    oldvalue = value
                    step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                    n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                    n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                    ## select data
                    data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                    data = np.asarray(data, dtype=np.float64, order='F')
                    realchunk = (data[:, self.header.index('Ncount')] > 10**-9)[index_n]
                    x = data[:, self.header.index(variable_name)][index_n]
                    value = value + x*realchunk
                    # get std from lmp output
                    ave_square_2_plus = data[:, self.header.index(sq_variable_name)][index_n]
                    ave_square_2 = ave_square_2 + ave_square_2_plus*realchunk
            value = value/self.n_ave
            ave_square_2 = ave_square_2/self.n_ave
            # calculate std
            totaln = self.n_ave*int(rr.log_variable["repeat_ave_chunk_wallforce"])
            if (ddof == 1 and totaln == 1):
                std = 0
            else:
                std2 = (ave_square_2 - value**2)*totaln/(totaln-ddof)
                std = std2**0.5
        else:
            # sum_diff_square_2 = sum of (xi-x_ave)**2
            sum_diff_square_2 = 0
            if index_n == None:
                for i in range(self.n_ave):
                    oldvalue = value
                    step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                    n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                    n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                    ## select data
                    data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                    data = np.asarray(data, dtype=np.float64, order='F')
                    x = data[:, self.header.index(variable_name)]
                    value = value + (
                        x - value
                        )/(i + 1)
                    sum_diff_square_2 += (x-value)*(x-oldvalue)

            else:
                for i in range(self.n_ave):
                    oldvalue = value
                    step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                    n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                    n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                    ## select data
                    data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                    data = np.asarray(data, dtype=np.float64, order='F')
                    x = data[:, self.header.index(variable_name)][index_n]
                    value = value + (
                        x - value
                        )/(i + 1)
                    sum_diff_square_2 += (x-value)*(x-oldvalue)
            # calculate std
            if (ddof == 1 and self.n_ave == 1):
                std = 0
            else:
                std = (
                    sum_diff_square_2/(self.n_ave-ddof)
                )**0.5
        return (value, std)

    def collectdataforplot(self, stepsarray, variable_name, index_n=None, ddof=1):
        # check array size by using first step
        onesteparrayshape = self.value_in_a_step_ave(stepsarray[0], variable_name, index_n, ddof)[0].shape
        # total shape
        totalshape = (len(stepsarray),) + onesteparrayshape
        # initialize array
        value_array = np.empty(totalshape)
        std_array = np.empty(totalshape)
        for i, step in enumerate(stepsarray):
            (value, std) = self.value_in_a_step_ave(stepsarray[i], variable_name, index_n, ddof)
            value_array[i] = value
            std_array[i] = std
        return (value_array, std_array)

    # data 1D-1D
    def datachunk_ave_one_step_XY(
        self,
        step, Y_name,
        Y_scale,
        X_scale=float(rr.log_variable['dp']),
        ):

        # reduce to 1D X_vector Y_vector
        X_vector = self.chunk_coor_array
        X_label = self.chunk_coor_array_label
        Y_vector, Y_std_vector = self.value_in_a_step_ave(step, Y_name)

        # rescale
        X_vector = X_vector/X_scale
        (Y_vector, Y_std_vector) = (Y_vector/Y_scale, Y_std_vector/Y_scale)

        # time from rotate
        time = time_from_start_rotate(step)

        return (time, X_vector, X_label, Y_vector, Y_std_vector)

    # plot 1D-1D
    def plot_XY(
        self,
        stepsarray,
        Y_name,
        Y_scale,
        Ylabel,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        fig, ax = plt.subplots()
        for step in stepsarray:
            (time, X_vector, Xlabel, Y_vector, Y_std_vector) = (
                self.datachunk_ave_one_step_XY(
                    step, Y_name,
                    Y_scale,
                    X_scale=float(rr.log_variable['dp']),
                )
            )

            # plot
            ax.errorbar(X_vector, Y_vector, yerr=Y_std_vector,
                    label="t={:.2e} s".format(time),
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                    )
        
        ax.set_xlabel(Xlabel + " ({:3g})".format(X_scale))
        ax.set_ylabel(Ylabel + " ({:3g})".format(Y_scale))
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)
    
    # save X-Y
    def save_XY_fix_k(
        self,
        stepsarray,
        Y_name,
        Y_scale,
        Ylabel,
        savepath,
        X_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        if onefigure:
            fig, ax = self.plot_XY(
                                    stepsarray,
                                    Y_name,
                                    Y_scale,                                        
                                    Ylabel,
                                    X_scale,
                                    )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + "-".join([str(number) for number in stepsarray]),
                figformat, ifpickle, bbox_inches,
                )
        else:
            for step in stepsarray:
                fig, ax = self.plot_XY(
                                        [step],
                                        Y_name,
                                        Y_scale,
                                        Ylabel,
                                        X_scale,
                                        )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + str(int(step)),
                figformat, ifpickle, bbox_inches,
                )

    # data at 1 points-time
    def datachunk_ave_one_step_time_variable(
        self,
        stepsarray, Y_name, Y_scale,
        index_n,
        ):
        
        Y_vector = np.empty(len(stepsarray), dtype=np.float64)
        Y_std_vector = np.empty(len(stepsarray), dtype=np.float64)
        # Y_vector
        for i, step in enumerate(stepsarray):
            (Y_vector[i], Y_std_vector[i]) = self.value_in_a_step_ave(step, Y_name, index_n=index_n)
            
        # scale
        Y_vector = Y_vector/Y_scale
        Y_std_vector = Y_std_vector/Y_scale
        # time from rotate
        time = time_from_start_rotate(stepsarray)
        return (time, Y_vector, Y_std_vector)

    # plot 1D-1D
    def plot_time_variable(
        self,
        stepsarray, Y_name, Y_scale,
        Ylabel,
        index_n,
        ):
        
        fig, ax = plt.subplots()
        (time, Y_vector, Y_std_vector) = self.datachunk_ave_one_step_time_variable(
            stepsarray, Y_name, Y_scale,
            index_n,
        )
        value_coor_dpunit = self.chunk_coor_array[index_n]/float(rr.log_variable['dp'])
        # plot
        ax.errorbar(time, Y_vector, yerr=Y_std_vector,
                label=self.chunk_coor_array_label + "={:.1f}, ".format(value_coor_dpunit),
                marker = ".",
                linestyle = 'None',
                markersize=12,
                )
        
        ax.set_xlabel("time")
        ax.set_ylabel(Ylabel + " ({:3g})".format(Y_scale))
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)

    # save time vs variable
    def save_time_variable(
        self,
        stepsarray, Y_name, Y_scale,
        Ylabel,
        index_n,
        savepath,
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        fig, ax = self.plot_time_variable(
            stepsarray, Y_name, Y_scale,
            Ylabel,
            index_n,
        )
        chunk.save_one_plot(
            fig, ax, savepath,
            "step_" + str(stepsarray[0]) + "_to_" + str(stepsarray[-1]) + "_n_" + str(index_n),
            figformat, ifpickle, bbox_inches,
            )

class chunkinwallforce(chunk1D):
    coodfilepath_rela_lmpfolder = "output/wall/chunk/coord_inwall"
    call_header_by_bettername = {
            "F_inwall_0": "v_inwall_per_atom_1", 
            "F_inwall_1": "v_inwall_per_atom_2", 
            "F_inwall_2": "v_inwall_per_atom_3",
            "F_inwall_0_sq": "v_inwall_per_atom_1_sq", 
            "F_inwall_1_sq": "v_inwall_per_atom_2_sq", 
            "F_inwall_2_sq": "v_inwall_per_atom_3_sq",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/wall/chunk/inwallforcefile", chunkinwallforce.coodfilepath_rela_lmpfolder)

# plot 1D-1D
    def plot_XFricCoeffi(
        self,
        stepsarray,
        Y_name,
        Y_scale,
        Ylabel,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        fig, ax = plt.subplots()
        for step in stepsarray:
            (time, X_vector, Xlabel, F0_vector, F0_std_vector) = (
                self.datachunk_ave_one_step_XY(
                    step, "v_inwall_per_atom_1",
                    Y_scale,
                    X_scale=float(rr.log_variable['dp']),
                )
            )

            (time, X_vector, Xlabel, F1_vector, F1_std_vector) = (
                self.datachunk_ave_one_step_XY(
                    step, "v_inwall_per_atom_2",
                    Y_scale,
                    X_scale=float(rr.log_variable['dp']),
                )
            )

            Y_vector = F0_vector/F1_vector
            Y_std_vector = 0
            # plot
            ax.plot(X_vector, Y_vector,
                    label="t={:.2e} s".format(time),
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                    )
        
        ax.set_xlabel(Xlabel + " ({:3g})".format(X_scale))
        ax.set_ylabel("Stress_ratio" + " ({:3g})".format(Y_scale))
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)
    
    # save X-Y
    def save_XFricCoeffi_fix_k(
        self,
        stepsarray,
        Y_name,
        Y_scale,
        Ylabel,
        savepath,
        X_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        if onefigure:
            fig, ax = self.plot_XFricCoeffi(
                                    stepsarray,
                                    Y_name,
                                    Y_scale,                                        
                                    Ylabel,
                                    X_scale,
                                    )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + "-".join([str(number) for number in stepsarray]),
                figformat, ifpickle, bbox_inches,
                )
        else:
            for step in stepsarray:
                fig, ax = self.plot_XFricCoeffi(
                                        [step],
                                        Y_name,
                                        Y_scale,
                                        Ylabel,
                                        X_scale,
                                        )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + str(int(step)),
                figformat, ifpickle, bbox_inches,
                )
class chunkoutwallforce(chunk1D):
    coodfilepath_rela_lmpfolder = "output/wall/chunk/coord_outwall"
    call_header_by_bettername = {
            "F_outwall_0": "v_outwall_per_atom_1", 
            "F_outwall_1": "v_outwall_per_atom_2", 
            "F_outwall_2": "v_outwall_per_atom_3",
            "F_outwall_0_sq": "v_outwall_per_atom_1_sq", 
            "F_outwall_1_sq": "v_outwall_per_atom_2_sq", 
            "F_outwall_2_sq": "v_outwall_per_atom_3_sq",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/wall/chunk/outwallforcefile", chunkoutwallforce.coodfilepath_rela_lmpfolder)
 
class chunkzbottomwallforce(chunk1D):
    coodfilepath_rela_lmpfolder = "output/wall/chunk/coord_zbottom"
    call_header_by_bettername = {
            "F_zbottom_0": "v_zbottom_per_atom_1", 
            "F_zbottom_1": "v_zbottom_per_atom_2", 
            "F_zbottom_2": "v_zbottom_per_atom_3",
            "F_zbottom_0_sq": "v_zbottom_per_atom_1_sq",
            "F_zbottom_1_sq": "v_zbottom_per_atom_2_sq",
            "F_zbottom_2_sq": "v_zbottom_per_atom_3_sq",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/wall/chunk/zbottomforcefile", chunkzbottomwallforce.coodfilepath_rela_lmpfolder)
        # coordinate in chunk
        try:
            self.chunk_coor_array = self.get_coord("Coord1")
        except:
            # old version chunk
            if rr.log_variable["shearwall"] == "zcylinder":
                self.chunk_coor_array = self.firstdata("v_r")
            else:
                sys.exit("chunk_method wrong")
        
        # chunk_coor_array_label
        if rr.log_variable["shearwall"] == "zcylinder":
            self.chunk_coor_array_label = "r"
        elif rr.log_variable["shearwall"] == "yplane":
            self.chunk_coor_array_label = "y"
        else:
            sys.exit("chunk_method wrong")


class chunk2D(chunk):
    coodfilepath_rela_lmpfolder = "output/wall/chunk/coord_chunk_2_3"
    def __init__(self, n_ave, lmp_path, f_path_rela_lmpfolder):
        super().__init__(n_ave, lmp_path, f_path_rela_lmpfolder)
        # coordinate in chunk
        try:
            self.chunk_x_array = self.get_coord("Coord1")
            self.chunk_y_array = self.get_coord("Coord2")
        except:
            # old version chunk
            if rr.log_variable["shearwall"] == "zcylinder":
                self.chunk_x_array = self.firstdata("v_r")
                self.chunk_y_array = self.firstdata("v_z")
            else:
                sys.exit("chunk_method wrong")
        # x-y coordinate in chunk
        if rr.log_variable["shearwall"] == "zcylinder":
            self.chunk_x_array_label = "r"
            self.chunk_y_array_label = "z"
        elif rr.log_variable["shearwall"] == "yplane":
            self.chunk_x_array_label = "y"
            self.chunk_y_array_label = "z"
        else:
            sys.exit("chunk_method wrong")

    def get_coord(self, variable_name):
        if variable_name in self.header:
            coordanswer = self.firstdata(variable_name)
        else:
            with open(self.lmp_path + chunk2D.coodfilepath_rela_lmpfolder) as f:
                
                lines = f.read().strip().split('\n')
                n_line_0 = 4
                n_line_1 = n_line_0 + self.n_line_in_a_step
                header = lines[2].split()[1:]
                ## select data
                data = [lines[t].split() for t in range(n_line_0, n_line_1)]
                coordanswer = np.asarray(data, dtype=np.float64, order='F')[:, header.index(variable_name)].reshape(n_1, n_2)
        return coordanswer

    def firstdata(self, variable_name):
        n_line_0 = 4
        n_line_1 = n_line_0 + self.n_line_in_a_step
        ## select data
        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
        return (
            np.asarray(data, dtype=np.float64, order='F')[:, self.header.index(variable_name)].reshape(n_1, n_2)
        )

    def value_in_a_step_ave(self, step, variable_name, index_n_1=None, index_n_2=None, ddof=1):
        # mean value
        value = 0
        # if lmp output std
        sq_variable_name = variable_name + "_sq"
        if sq_variable_name in self.header:
            ave_square_2 = 0

            if index_n_1 == None:
                if index_n_2 == None:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        realchunk = (data[:, self.header.index('Ncount')] > 10**-9).reshape(n_1, n_2)
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)
                        value = value + x*realchunk
                        ave_square_2_plus = data[:, self.header.index(sq_variable_name)].reshape(n_1, n_2)
                        ave_square_2 = ave_square_2 + ave_square_2_plus*realchunk
                            
                else:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        realchunk = (data[:, self.header.index('Ncount')] > 10**-9).reshape(n_1, n_2)[:, index_n_2]
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)[:, index_n_2]
                        value = value + x*realchunk
                        ave_square_2_plus = data[:, self.header.index(sq_variable_name)].reshape(n_1, n_2)[:, index_n_2]
                        ave_square_2 = ave_square_2 + ave_square_2_plus*realchunk
            else:
                if index_n_2 == None:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        realchunk = (data[:, self.header.index('Ncount')] > 10**-9).reshape(n_1, n_2)[index_n_1, :]
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)[index_n_1, :]
                        value = value + x*realchunk
                        ave_square_2_plus = data[:, self.header.index(sq_variable_name)].reshape(n_1, n_2)[index_n_1, :]
                        ave_square_2 = ave_square_2 + ave_square_2_plus*realchunk
                else:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        realchunk = (data[:, self.header.index('Ncount')] > 10**-9).reshape(n_1, n_2)[index_n_1, index_n_2]
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)[index_n_1, index_n_2]
                        value = value + x*realchunk
                        ave_square_2_plus = data[:, self.header.index(sq_variable_name)].reshape(n_1, n_2)[index_n_1, index_n_2]
                        ave_square_2 = ave_square_2 + ave_square_2_plus*realchunk
            # calculate std
            value = value/self.n_ave
            ave_square_2 = ave_square_2/self.n_ave
            totaln = self.n_ave*int(rr.log_variable["repeat_ave_chunk_momentum_mass_field"])
            if (ddof == 1 and totaln == 1):
                std = 0
            else:
                std2 = (ave_square_2 - value**2)*totaln/(totaln-ddof)
                std = std2**0.5
        
        else:    
            # sum_diff_square_2 = sum of (xi-x_ave)**2
            sum_diff_square_2 = 0

            if index_n_1 == None:
                if index_n_2 == None:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)
                        value = value + (
                            x - value
                            )/(i + 1)
                        sum_diff_square_2 += (x-value)*(x-oldvalue)
                else:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)[:, index_n_2]
                        value = value + (
                            x - value
                            )/(i + 1)
                        sum_diff_square_2 += (x-value)*(x-oldvalue)
                        
            else:
                if index_n_2 == None:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)[index_n_1, :]
                        value = value + (
                            x - value
                            )/(i + 1)
                        sum_diff_square_2 += (x-value)*(x-oldvalue)
                else:
                    for i in range(self.n_ave):
                        oldvalue = value
                        step_inloop = int(step - (self.n_ave-1)*self.d_step/2 + i*self.d_step)
                        n_line_0 = int(int(step_inloop - self.step_first_in_file)/self.d_step)*(self.n_line_in_a_step+1) + 4
                        n_line_1 = int(n_line_0 + self.n_line_in_a_step)
                        ## select data
                        data = [self.lines[t].split() for t in range(n_line_0, n_line_1)]
                        data = np.asarray(data, dtype=np.float64, order='F')
                        x = data[:, self.header.index(variable_name)].reshape(n_1, n_2)[index_n_1, index_n_2]
                        value = value + (
                            x - value
                            )/(i + 1)
                        sum_diff_square_2 += (x-value)*(x-oldvalue)

            if (ddof == 1 and self.n_ave == 1):
                std = 0
            else:
                std = (
                    sum_diff_square_2/(self.n_ave-ddof)
                )**0.5
        return (value, std)

    def collectdataforplot(self, stepsarray, variable_name, index_n_1=None, index_n_2=None, ddof=1):
        # check array size by using first step
        onesteparrayshape = self.value_in_a_step_ave(stepsarray[0], variable_name, index_n_1, index_n_2, ddof)[0].shape
        # total shape
        totalshape = (len(stepsarray),) + onesteparrayshape
        # initialize array
        value_array = np.empty(totalshape)
        std_array = np.empty(totalshape)
        for i, step in enumerate(stepsarray):
            (value, std) = self.value_in_a_step_ave(stepsarray[i], variable_name, index_n_1, index_n_2, ddof)
            value_array[i] = value
            std_array[i] = std
        return (value_array, std_array)

    # data quiver 2D 2D
    def datachunk_ave_one_step_quiver_x23(
        self,
        step, Q_name, V_name, Q_scale, V_scale,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        ):

        # V array
        V_array = self.value_in_a_step_ave(step, V_name, index_n_1=None, index_n_2=None)[0]

        # Q array, 0 if Q_name None
        if Q_name == None:
            Q_array = np.zeros_like(V_array)
        else:
            Q_array = self.value_in_a_step_ave(step, Q_name, index_n_1=None, index_n_2=None)[0]

        # rescale
        x_array = self.chunk_x_array/x_scale
        y_array = self.chunk_y_array/y_scale
        Q_array = Q_array/Q_scale
        V_array = V_array/V_scale

        # time from rotate
        time = time_from_start_rotate(step)

        return [time, x_array, y_array, Q_array, V_array]

    # plot quiver
    def plotquiver(
        self,
        time, x_array, y_array, Q_array, V_array,
        step, Q_name, V_name, Q_scale, V_scale,
        quiver_scale=0.1, label_scale=0.2,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        ):
        fig1, ax1 = plt.subplots()
        #fig1.figsize = [12.8, 9.6]
        chunk.plot_quiver_position_label(fig1, ax1)
        #ax1.set_title('velocity field r-z direction (average over theta)')
        Q = ax1.quiver(x_array, y_array, Q_array, V_array,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                    )

        
        ax1.set_xlabel(map_dim_index_to_coordinate[1] + " ({:3g})".format(x_scale))
        ax1.set_ylabel(map_dim_index_to_coordinate[2] + " ({:3g})".format(y_scale))
        self.fig1 = fig1
        self.ax1 = ax1
        self.quiverobject = Q

    def quiveraddkey(
        self, keylabel, label_scale,
        ):
        self.ax1.quiverkey(
            self.quiverobject, 0.2, 0.95, label_scale,
            label = keylabel,
            labelpos='E',
            coordinates='figure', angle=90
            )

    # plot quiver
    def plotquiver_x23(
        self,
        step, Q_name, V_name, Q_scale, V_scale,
        quiver_scale=0.1, label_scale=0.2,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        ):

        [time, x_array, y_array, Q_array, V_array] = (
            self.datachunk_ave_one_step_quiver_x23(
            step, Q_name, V_name, Q_scale, V_scale,
            x_scale, y_scale,
            )
        )
        self.plotquiver(
        time, x_array, y_array, Q_array, V_array,
        step, Q_name, V_name, Q_scale, V_scale,
        quiver_scale=0.1, label_scale=0.2,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        )
        self.quiveraddkey("equal (" + Q_name + ", " + V_name + ") = {:.2e}".format(label_scale) + ". At {:.2e} s".format(time), label_scale)
    
    # save quiver
    def savequiver_x23(
        self,
        stepsarray,
        Q_name, V_name, Q_scale, V_scale,
        savepath,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False,
        quiver_scale=0.1, label_scale=0.2,
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        for step in stepsarray:
            self.plotquiver_x23(
                step, Q_name, V_name, Q_scale, V_scale,
                quiver_scale, label_scale,
                x_scale, y_scale,
            )
            chunk.save_one_plot(self.fig1, self.ax1, savepath, str(int(step)), figformat, ifpickle)

    # data diff quiver 2D 2D
    def datachunk_ave_one_step_quiver_diff_x23(
        self,
        step, Q_name, V_name, Q_scale, V_scale,
        diff_position_index,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        ):
        
        [time, x_array, y_array, Q_array, V_array] = self.datachunk_ave_one_step_quiver_x23(
                                                                                            step, Q_name, V_name, Q_scale, V_scale,
                                                                                            x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
                                                                                            )
        diff_along_array_dim = position_index_to_array_dim_index[diff_position_index]
        Ncount_array = self.value_in_a_step_ave(step, 'Ncount', index_n_1=None, index_n_2=None)[0]
        mask_Ncount_small_larger_1_array = (Ncount_array>=1)
        if diff_along_array_dim == 0:
            middle_point_x_array = (x_array[:-1,:] + x_array[1:,:])/2
            middle_point_y_array = (y_array[:-1,:] + y_array[1:,:])/2
            vector_diffQ = np.diff(Q_array, axis=diff_along_array_dim)/np.diff(x_array, axis=diff_along_array_dim)
            vector_diffV = np.diff(V_array, axis=diff_along_array_dim)/np.diff(x_array, axis=diff_along_array_dim)
            mask_array = np.logical_and(mask_Ncount_small_larger_1_array[:-1,:],mask_Ncount_small_larger_1_array[1:,:])
        elif diff_along_array_dim == 1:
            middle_point_x_array = (x_array[:,:-1] + x_array[:,1:])/2
            middle_point_y_array = (y_array[:,:-1] + y_array[:,1:])/2
            vector_diffQ = np.diff(Q_array, axis=diff_along_array_dim)/np.diff(y_array, axis=diff_along_array_dim)
            vector_diffV = np.diff(V_array, axis=diff_along_array_dim)/np.diff(y_array, axis=diff_along_array_dim)
            mask_array = np.logical_and(mask_Ncount_small_larger_1_array[:,:-1],mask_Ncount_small_larger_1_array[:,1:])

        middle_point_x_array = mask_array*middle_point_x_array
        middle_point_y_array = mask_array*middle_point_y_array
        vector_diffQ = mask_array*vector_diffQ
        vector_diffV = mask_array*vector_diffV

        return [time, middle_point_x_array, middle_point_y_array, vector_diffQ, vector_diffV]

    # plot quiver diff quiver 2D 2D
    def plotquiver_diff_x23(
        self,
        step, Q_name, V_name, Q_scale, V_scale,
        diff_position_index,
        quiver_scale=0.1, label_scale=0.2,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        ):

        [time, x_array, y_array, Q_array, V_array] = (
            self.datachunk_ave_one_step_quiver_diff_x23(
            step, Q_name, V_name, Q_scale, V_scale,
            diff_position_index,
            x_scale, y_scale,
            )
        )
        self.plotquiver(
        time, x_array, y_array, Q_array, V_array,
        step, Q_name, V_name, Q_scale, V_scale,
        quiver_scale=0.1, label_scale=0.2,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        )
        self.quiveraddkey("equal (diff_" + Q_name + ", diff_" + V_name + ") = {:.2e}".format(label_scale) + ". At {:.2e} s".format(time), label_scale)
    
    # save quiver diff quiver 2D 2D
    def savequiver_diff_x23(
        self,
        stepsarray,
        Q_name, V_name, Q_scale, V_scale,
        diff_position_index,
        savepath,
        x_scale=float(rr.log_variable['dp']), y_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False,
        quiver_scale=0.1, label_scale=0.2,
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        for step in stepsarray:
            self.plotquiver_diff_x23(
                step, Q_name, V_name, Q_scale, V_scale,
                diff_position_index,
                quiver_scale, label_scale,
                x_scale, y_scale,
            )
            chunk.save_one_plot(self.fig1, self.ax1, savepath, str(int(step)), figformat, ifpickle)

    # data 1D-1D diff
    def datachunk_ave_one_step_XY_diff_fix_k(
        self,
        step, coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index,
        diff_coord_index,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        # coord_index_horizental_in_plot should not equal k
        if coord_index_horizental_in_plot == k:
            sys.exit("coord_index_horizental_in_plot should not equal k")

        Y_array, Y_std_array = self.value_in_a_step_ave(step, Y_name, index_n_1=None, index_n_2=None)

        if position_index_to_array_dim_index[coord_index_horizental_in_plot] == 0:
            X_array = self.chunk_x_array
            X_label = self.chunk_x_array_label
        elif position_index_to_array_dim_index[coord_index_horizental_in_plot] == 1:
            X_array = self.chunk_y_array
            X_label = self.chunk_y_array_label
        
        if position_index_to_array_dim_index[diff_coord_index] == 0:
            diff_label = self.chunk_x_array_label
            diff_coord_array = self.chunk_x_array
            middle_X_array = ((X_array[:-1,:-1] + X_array[1:,:-1])/2 + (X_array[:-1,1:] + X_array[1:,1:])/2)/2
            d_diff_coord_array = (diff_coord_array[:-1,:-1] - diff_coord_array[1:,:-1] + diff_coord_array[:-1,1:] - diff_coord_array[1:,1:])/2
            d_Y_array = (Y_array[:-1,:-1] - Y_array[1:,:-1] + Y_array[:-1,1:] - Y_array[1:,1:])/2
            d_Y_std_array = (
                (Y_std_array[:-1,:-1])**2 + (Y_std_array[1:,:-1])**2 + (Y_std_array[:-1,1:])**2 + (Y_std_array[1:,1:])**2
                )**0.5
            
        elif position_index_to_array_dim_index[diff_coord_index] == 1:
            diff_label = self.chunk_y_array_label
            diff_coord_array = self.chunk_y_array
            middle_X_array = ((X_array[:-1, :-1] + X_array[:-1, 1:])/2 + (X_array[1:, :-1] + X_array[1:, 1:])/2)/2
            d_diff_coord_array = (diff_coord_array[:-1, :-1] - diff_coord_array[:-1, 1:] + diff_coord_array[1:, :-1] - diff_coord_array[1:, 1:])/2
            d_Y_array = (Y_array[:-1, :-1] - Y_array[:-1, 1:] + Y_array[1:, :-1] - Y_array[1:, 1:])/2
            d_Y_std_array = (
                (Y_std_array[:-1,:-1])**2 + (Y_std_array[:-1,1:])**2 + (Y_std_array[1:,:-1])**2 + (Y_std_array[1:,1:])**2
                )**0.5
            
        if position_index_to_array_dim_index[k] == 0:
            middle_X_vector = middle_X_array[k_index, :]
            d_Y_vector = d_Y_array[k_index, :]
            d_diff_coord_vector = d_diff_coord_array[k_index, :]
            d_Y_std_vector = d_Y_std_array[k_index, :]
        elif position_index_to_array_dim_index[k] == 1:
            middle_X_vector = middle_X_array[:, k_index]
            d_Y_vector = d_Y_array[:, k_index]
            d_diff_coord_vector = d_diff_coord_array[:, k_index]
            d_Y_std_vector = d_Y_std_array[:, k_index]
        if np.any(d_diff_coord_vector==0):
            breakpoint()
        d_Y_d_diff_coord_vector = d_Y_vector/d_diff_coord_vector
        d_Y_d_diff_coord_std_vector = d_Y_std_vector/d_diff_coord_vector

        # rescale
        middle_X_vector = middle_X_vector/X_scale
        (d_Y_d_diff_coord_vector, d_Y_d_diff_coord_std_vector) = (d_Y_d_diff_coord_vector/Y_scale*X_scale, d_Y_d_diff_coord_std_vector/Y_scale*X_scale)

        # time from rotate
        time = time_from_start_rotate(step)

        return (time, middle_X_vector, X_label, d_Y_d_diff_coord_vector, d_Y_d_diff_coord_std_vector, diff_label)
    # plot 1D-1D diff
    def plot_XY_diff_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index_array,
        Ylabel,
        diff_coord_index,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        fig, ax = plt.subplots()
        for k_index in k_index_array:
            for step in stepsarray:
                (time, X_vector, Xlabel, Y_vector, Y_std_vector, diff_label) = (
                    self.datachunk_ave_one_step_XY_diff_fix_k(
                        step, coord_index_horizental_in_plot, Y_name,
                        Y_scale,
                        k, k_index,
                        diff_coord_index,
                        X_scale=float(rr.log_variable['dp']),
                    )
                )
                
                # k_value
                if position_index_to_array_dim_index[k] == 0:
                    k_value = self.chunk_x_array[k_index, 0]
                    
                elif position_index_to_array_dim_index[k] == 1:
                    k_value = self.chunk_y_array[0, k_index]

                k_value = k_value/float(rr.log_variable['dp'])
                # plot
                ax.errorbar(X_vector, Y_vector, yerr=Y_std_vector,
                        label="t={:.2e} s".format(time) + ", " + map_dim_index_to_coordinate[k] + "=" + "{:.2f}".format(k_value),
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                        )

        ax.set_xlabel(Xlabel + " ({:3g})".format(X_scale))
        ax.set_ylabel("d(" + Ylabel + ")/d(" + diff_label + ")" + " ({:3g})".format(Y_scale/X_scale))
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)
    # save X-Y diff
    def save_XY_diff_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index_array,
        Ylabel,
        diff_coord_index,
        savepath,
        X_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        if onefigure:
            fig, ax = self.plot_XY_diff_fix_k(
                                        stepsarray,
                                        coord_index_horizental_in_plot, Y_name,
                                        Y_scale,
                                        k, k_index_array,
                                        Ylabel,
                                        diff_coord_index,
                                        X_scale,
                                        )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + "-".join([str(number) for number in stepsarray]) + map_dim_index_to_coordinate[k] + "_" + "-".join([str(number) for number in k_index_array]),
                figformat, ifpickle, bbox_inches,
                )
        else:
            for k_index in k_index_array:
                for step in stepsarray:
                    fig, ax = self.plot_XY_diff_fix_k(
                                                [step],
                                                coord_index_horizental_in_plot, Y_name,
                                                Y_scale,
                                                k, [k_index],
                                                Ylabel,
                                                diff_coord_index,
                                                X_scale,
                                                )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + str(int(step)) + map_dim_index_to_coordinate[k] + "_" + str(k_index),
                figformat, ifpickle, bbox_inches,
                )

    # data 1D-1D
    def datachunk_ave_one_step_XY_fix_k(
        self,
        step, coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        # coord_index_horizental_in_plot should not equal k
        if coord_index_horizental_in_plot == k:
            sys.exit("coord_index_horizental_in_plot should not equal k")

        # reduce to 1D X_vector Y_vector
        if position_index_to_array_dim_index[k] == 0:
            X_vector = self.chunk_y_array[k_index, :]
            X_label = self.chunk_y_array_label
            Y_vector, Y_std_vector = self.value_in_a_step_ave(step, Y_name, index_n_1=k_index, index_n_2=None)
        elif position_index_to_array_dim_index[k] == 1:
            X_vector = self.chunk_x_array[:, k_index]
            X_label = self.chunk_x_array_label
            Y_vector, Y_std_vector = self.value_in_a_step_ave(step, Y_name, index_n_1=None, index_n_2=k_index)

        # rescale
        X_vector = X_vector/X_scale
        (Y_vector, Y_std_vector) = (Y_vector/Y_scale, Y_std_vector/Y_scale)

        # time from rotate
        time = time_from_start_rotate(step)

        return (time, X_vector, X_label, Y_vector, Y_std_vector)

    # plot 1D-1D
    def plot_XY_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index_array,
        Ylabel,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        fig, ax = plt.subplots()
        for k_index in k_index_array:
            for step in stepsarray:

                (time, X_vector, Xlabel, Y_vector, Y_std_vector) = (
                    self.datachunk_ave_one_step_XY_fix_k(
                        step, coord_index_horizental_in_plot, Y_name,
                        Y_scale,
                        k, k_index,
                        X_scale=float(rr.log_variable['dp']),
                    )
                )
                
                # k_value
                if position_index_to_array_dim_index[k] == 0:
                    k_value = self.chunk_x_array[k_index, 0]
                    
                elif position_index_to_array_dim_index[k] == 1:
                    k_value = self.chunk_y_array[0, k_index]

                k_value = k_value/float(rr.log_variable['dp'])
                # plot
                ax.errorbar(X_vector, Y_vector, yerr=Y_std_vector,
                        label="t={:.2e} s".format(time) + ", " + map_dim_index_to_coordinate[k] + "=" + "{:.2f}".format(k_value),
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                        )
        ax.set_xlabel(Xlabel + " ({:3g})".format(X_scale))
        ax.set_ylabel(Ylabel + " ({:3g})".format(Y_scale))
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)
    
    # save X-Y
    def save_XY_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index_array,
        Ylabel,
        savepath,
        X_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        if onefigure:
            fig, ax = self.plot_XY_fix_k(
                                        stepsarray,
                                        coord_index_horizental_in_plot, Y_name,
                                        Y_scale,
                                        k, k_index_array,
                                        Ylabel,
                                        X_scale,
                                        )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + "-".join([str(number) for number in stepsarray]) + map_dim_index_to_coordinate[k] + "_" + "-".join([str(number) for number in k_index_array]),
                figformat, ifpickle, bbox_inches,
                )
        else:
            for k_index in k_index_array:
                for step in stepsarray:
                    fig, ax = self.plot_XY_fix_k(
                                                [step],
                                                coord_index_horizental_in_plot, Y_name,
                                                Y_scale,
                                                k, [k_index],
                                                Ylabel,
                                                X_scale,
                                                )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + str(int(step)) + map_dim_index_to_coordinate[k] + "_" + str(k_index),
                figformat, ifpickle, bbox_inches,
                )

    # data at 1 points-time
    def datachunk_ave_one_step_time_variable(
        self,
        stepsarray, Y_name, Y_scale,
        index_n_1, index_n_2,
        ):
        
        Y_vector = np.empty(len(stepsarray), dtype=np.float64)
        Y_std_vector = np.empty(len(stepsarray), dtype=np.float64)
        # Y_vector
        for i, step in enumerate(stepsarray):
            (Y_vector[i], Y_std_vector[i]) = self.value_in_a_step_ave(step, Y_name, index_n_1=index_n_1, index_n_2=index_n_2)
            
        # scale
        Y_vector = Y_vector/Y_scale
        Y_std_vector = Y_std_vector/Y_scale
        # time from rotate
        time = time_from_start_rotate(stepsarray)
        return (time, Y_vector, Y_std_vector)

    # plot 1D-1D
    def plot_time_variable(
        self,
        stepsarray, Y_name, Y_scale,
        Ylabel,
        index_n_1, index_n_2,
        ):
        
        fig, ax = plt.subplots()
        (time, Y_vector, Y_std_vector) = self.datachunk_ave_one_step_time_variable(
            stepsarray, Y_name, Y_scale,
            index_n_1, index_n_2,
        )
        value_x_dpunit = self.chunk_x_array[index_n_1, index_n_2]/float(rr.log_variable['dp'])
        value_y_dpunit = self.chunk_y_array[index_n_1, index_n_2]/float(rr.log_variable['dp'])
        # plot
        ax.errorbar(time, Y_vector, yerr=Y_std_vector,
                label=self.chunk_x_array_label + "={:.1f}, ".format(value_x_dpunit) + self.chunk_y_array_label + "={:.1f}".format(value_y_dpunit),
                marker = ".",
                linestyle = 'None',
                markersize=12,
                )
        
        ax.set_xlabel("time")
        ax.set_ylabel(Ylabel + " ({:3g})".format(Y_scale))
        
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)

    # save time vs variable
    def save_time_variable(
        self,
        stepsarray, Y_name, Y_scale,
        Ylabel,
        index_n_1, index_n_2,
        savepath,
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        fig, ax = self.plot_time_variable(
            stepsarray, Y_name, Y_scale,
            Ylabel,
            index_n_1, index_n_2,
        )
        chunk.save_one_plot(
            fig, ax, savepath,
            "step_" + str(stepsarray[0]) + "_to_" + str(stepsarray[-1]) + "_n_1_" + str(index_n_1) + "_n_2_" + str(index_n_2),
            figformat, ifpickle, bbox_inches,
            )

# chunk stress
class chunkstress(chunk2D):

    call_header_by_bettername = {
            "stress_00": "c_stress[1]",
            "stress_11": "c_stress[2]",
            "stress_22": "c_stress[3]",
            "stress_01": "c_stress[4]",
            "stress_02": "c_stress[5]",
            "stress_12": "c_stress[6]",
            "stress_00_sq": "c_stress[1]_sq",
            "stress_11_sq": "c_stress[2]_sq",
            "stress_22_sq": "c_stress[3]_sq",
            "stress_01_sq": "c_stress[4]_sq",
            "stress_02_sq": "c_stress[5]_sq",
            "stress_12_sq": "c_stress[6]_sq",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/stress/fix.stress.all")
        # stress_variable_name dictionary

    # data 1D-1D
    def datachunk_ave_one_step_pressure_fix_k(
        self,
        step, coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index,
        X_scale=float(rr.log_variable['dp']),
        ):
        (time, X_vector, X_label, stress_00_vector, stress_00_std_vector) = self.datachunk_ave_one_step_XY_fix_k(
        step, coord_index_horizental_in_plot, self.call_header_by_bettername["stress_00"],
        Y_scale,
        k, k_index,
        X_scale=float(rr.log_variable['dp']),
        )
        (time, X_vector, X_label, stress_11_vector, stress_11_std_vector) = self.datachunk_ave_one_step_XY_fix_k(
        step, coord_index_horizental_in_plot, self.call_header_by_bettername["stress_11"],
        Y_scale,
        k, k_index,
        X_scale=float(rr.log_variable['dp']),
        )
        (time, X_vector, X_label, stress_22_vector, stress_22_std_vector) = self.datachunk_ave_one_step_XY_fix_k(
        step, coord_index_horizental_in_plot, self.call_header_by_bettername["stress_22"],
        Y_scale,
        k, k_index,
        X_scale=float(rr.log_variable['dp']),
        )

        pressure_vector = -1/3*(stress_00_vector + stress_11_vector + stress_22_vector)
        pressure_std_vector = 0

        return (time, X_vector, X_label, pressure_vector, pressure_std_vector)        

# chunk momentum
class chunkmomentum(chunk2D):

    call_header_by_bettername = {
            "mv_1": "v_mv2",
            "mv_0": "v_mv1",
            "mv_2": "v_mv3",
            "mass": "c_m1",
            "mv_1_sq": "v_mv2_sq",
            "mv_0_sq": "v_mv1_sq",
            "mv_2_sq": "v_mv3_sq",
            "mass_sq": "c_m1_sq",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/momentum_mass_field/fix.momentum_mass_field.all")

    # data 1D-1D diff
    def datachunk_ave_one_step_XY_shear_rate_fix_k(
        self,
        step, coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index,
        diff_coord_index,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        # coord_index_horizental_in_plot should not equal k
        if coord_index_horizental_in_plot == k:
            sys.exit("coord_index_horizental_in_plot should not equal k")

        mass_array, mass_std_array = self.value_in_a_step_ave(step, chunkmomentum.call_header_by_bettername["mass"], index_n_1=None, index_n_2=None)

        Y_array, Y_std_array = self.value_in_a_step_ave(step, Y_name, index_n_1=None, index_n_2=None)

        Y_array = Y_array/mass_array
        Y_std_array = Y_std_array/mass_array

        if position_index_to_array_dim_index[coord_index_horizental_in_plot] == 0:
            X_array = self.chunk_x_array
            X_label = self.chunk_x_array_label
        elif position_index_to_array_dim_index[coord_index_horizental_in_plot] == 1:
            X_array = self.chunk_y_array
            X_label = self.chunk_y_array_label
        if position_index_to_array_dim_index[diff_coord_index] == 0:
            diff_label = self.chunk_x_array_label
            diff_coord_array = self.chunk_x_array
            middle_X_array = ((X_array[:-1,:-1] + X_array[1:,:-1])/2 + (X_array[:-1,1:] + X_array[1:,1:])/2)/2
            d_diff_coord_array = (diff_coord_array[:-1,:-1] - diff_coord_array[1:,:-1] + diff_coord_array[:-1,1:] - diff_coord_array[1:,1:])/2
            d_Y_array = (Y_array[:-1,:-1] - Y_array[1:,:-1] + Y_array[:-1,1:] - Y_array[1:,1:])/2

            d_Y_std_array = (
                (Y_std_array[:-1,:-1])**2 + (Y_std_array[1:,:-1])**2 + (Y_std_array[:-1,1:])**2 + (Y_std_array[1:,1:])**2
                )**0.5
            
        elif position_index_to_array_dim_index[diff_coord_index] == 1:
            diff_label = self.chunk_y_array_label
            diff_coord_array = self.chunk_y_array
            middle_X_array = ((X_array[:-1, :-1] + X_array[:-1, 1:])/2 + (X_array[1:, :-1] + X_array[1:, 1:])/2)/2
            d_diff_coord_array = (diff_coord_array[:-1, :-1] - diff_coord_array[:-1, 1:] + diff_coord_array[1:, :-1] - diff_coord_array[1:, 1:])/2
            d_Y_array = (Y_array[:-1, :-1] - Y_array[:-1, 1:] + Y_array[1:, :-1] - Y_array[1:, 1:])/2
            d_Y_std_array = (
                (Y_std_array[:-1,:-1])**2 + (Y_std_array[:-1,1:])**2 + (Y_std_array[1:,:-1])**2 + (Y_std_array[1:,1:])**2
                )**0.5
            
        if position_index_to_array_dim_index[k] == 0:
            middle_X_vector = middle_X_array[k_index, :]
            d_Y_vector = d_Y_array[k_index, :]
            d_diff_coord_vector = d_diff_coord_array[k_index, :]
            d_Y_std_vector = d_Y_std_array[k_index, :]
        elif position_index_to_array_dim_index[k] == 1:
            middle_X_vector = middle_X_array[:, k_index]
            d_Y_vector = d_Y_array[:, k_index]
            d_diff_coord_vector = d_diff_coord_array[:, k_index]
            d_Y_std_vector = d_Y_std_array[:, k_index]
        if np.any(d_diff_coord_vector==0):
            breakpoint()
        d_Y_d_diff_coord_vector = d_Y_vector/d_diff_coord_vector
        d_Y_d_diff_coord_std_vector = d_Y_std_vector/d_diff_coord_vector

        # rescale
        middle_X_vector = middle_X_vector/X_scale
        (d_Y_d_diff_coord_vector, d_Y_d_diff_coord_std_vector) = (d_Y_d_diff_coord_vector/Y_scale*X_scale, d_Y_d_diff_coord_std_vector/Y_scale*X_scale)

        # time from rotate
        time = time_from_start_rotate(step)

        return (time, middle_X_vector, X_label, d_Y_d_diff_coord_vector, d_Y_d_diff_coord_std_vector, diff_label)

    # plot 1D-1D diff
    def plot_XY_shear_rate_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index_array,
        Ylabel,
        diff_coord_index,
        X_scale=float(rr.log_variable['dp']),
        ):
        
        fig, ax = plt.subplots()
        for k_index in k_index_array:
            for step in stepsarray:
                (time, X_vector, Xlabel, Y_vector, Y_std_vector, diff_label) = (
                    self.datachunk_ave_one_step_XY_shear_rate_fix_k(
                        step, coord_index_horizental_in_plot, Y_name,
                        Y_scale,
                        k, k_index,
                        diff_coord_index,
                        X_scale=float(rr.log_variable['dp']),
                    )
                )
                
                # k_value
                if position_index_to_array_dim_index[k] == 0:
                    k_value = self.chunk_x_array[k_index, 0]
                    
                elif position_index_to_array_dim_index[k] == 1:
                    k_value = self.chunk_y_array[0, k_index]

                k_value = k_value/float(rr.log_variable['dp'])
                # plot
                ax.errorbar(X_vector, Y_vector, yerr=Y_std_vector,
                        label="t={:.2e} s".format(time) + ", " + map_dim_index_to_coordinate[k] + "=" + "{:.2f}".format(k_value),
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                        )

        ax.set_xlabel(Xlabel + " ({:3g})".format(X_scale))
        ax.set_ylabel("strain_rate_" + Ylabel + "_" + diff_label + " ({:3g})".format(Y_scale/X_scale))
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)
    # save X-Y diff
    def save_XY_shear_rate_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index_array,
        Ylabel,
        diff_coord_index,
        savepath,
        X_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        if onefigure:
            fig, ax = self.plot_XY_shear_rate_fix_k(
                                        stepsarray,
                                        coord_index_horizental_in_plot, Y_name,
                                        Y_scale,
                                        k, k_index_array,
                                        Ylabel,
                                        diff_coord_index,
                                        X_scale,
                                        )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + "-".join([str(number) for number in stepsarray]) + map_dim_index_to_coordinate[k] + "_" + "-".join([str(number) for number in k_index_array]),
                figformat, ifpickle, bbox_inches,
                )
        else:
            for k_index in k_index_array:
                for step in stepsarray:
                    fig, ax = self.plot_XY_shear_rate_fix_k(
                                                [step],
                                                coord_index_horizental_in_plot, Y_name,
                                                Y_scale,
                                                k, [k_index],
                                                Ylabel,
                                                diff_coord_index,
                                                X_scale,
                                                )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + str(int(step)) + map_dim_index_to_coordinate[k] + "_" + str(k_index),
                figformat, ifpickle, bbox_inches,
                )

class chunkmuI(chunk2D):

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/momentum_mass_field/fix.momentum_mass_field.all")
    
    def data_to_plot_mu_I(self, step, coord_index_horizental_in_plot,
        i, j,
        k, k_index,
        X_scale=float(rr.log_variable['dp']),
        ):
        if position_index_to_array_dim_index[k] + 1 == 1:
            n_k = n_1
        elif position_index_to_array_dim_index[k] + 1 == 2:
            n_k = n_2
        else:
            sys.exit("k is wrong")
        k_index_middle = np.arange(n_k-1)[k_index]
        stress_scale = float(rr.log_variable['den'])*float(rr.log_variable['g'])*float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['dp'])

        stress_data = chunkstress(self.n_ave, self.lmp_path)
        pressure_vector= 1/4*(
            stress_data.datachunk_ave_one_step_pressure_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle,
                X_scale=float(rr.log_variable['dp']),
                )[3][0:-1]
            +stress_data.datachunk_ave_one_step_pressure_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle,
                X_scale=float(rr.log_variable['dp']),
                )[3][1:]
            +stress_data.datachunk_ave_one_step_pressure_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle+1,
                X_scale=float(rr.log_variable['dp']),
                )[3][0:-1]
            +stress_data.datachunk_ave_one_step_pressure_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle+1,
                X_scale=float(rr.log_variable['dp']),
                )[3][1:]
            )
        stress_ij = 1/4*(
            stress_data.datachunk_ave_one_step_XY_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle,
                X_scale=float(rr.log_variable['dp']),
            )[3][0:-1]
            +stress_data.datachunk_ave_one_step_XY_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle,
                X_scale=float(rr.log_variable['dp']),
            )[3][1:]
            +stress_data.datachunk_ave_one_step_XY_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle+1,
                X_scale=float(rr.log_variable['dp']),
            )[3][0:-1]
            +stress_data.datachunk_ave_one_step_XY_fix_k(
                step, coord_index_horizental_in_plot, chunkstress.call_header_by_bettername["stress_" + str(i) + str(j)],
                stress_scale,
                k, k_index_middle+1,
                X_scale=float(rr.log_variable['dp']),
            )[3][1:]
        )

        velocity_scale = float(rr.log_variable['in_velocity'])
        if velocity_scale < 0:
            velocity_scale = -velocity_scale
        shear_rate_scale = velocity_scale/(float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['dp']))
        momentum_data = chunkmomentum(self.n_ave, self.lmp_path)
        shear_rate_abs_2 = 0
        for (i_loop, j_loop) in [
            (0, 1),
            (0, 2),
            (1, 2),
            (2, 1),
            ]:

            (time, middle_X_vector, X_label, d_Y_d_diff_coord_vector, d_Y_d_diff_coord_std_vector, diff_label) = momentum_data.datachunk_ave_one_step_XY_shear_rate_fix_k(
                step, coord_index_horizental_in_plot, chunkmomentum.call_header_by_bettername["mv_" + str(i_loop)],
                shear_rate_scale,
                k, k_index_middle,
                j_loop,
                X_scale=float(rr.log_variable['dp']),
                )
            #if len(d_Y_d_diff_coord_vector) == 12:
            #    breakpoint()
            if i_loop == i and j_loop == j:
                shear_rate_ij = d_Y_d_diff_coord_vector
                shear_rate_ij_std = d_Y_d_diff_coord_std_vector

            shear_rate_abs_2 = shear_rate_abs_2 + d_Y_d_diff_coord_vector**2
        #breakpoint()
        shear_rate_abs = shear_rate_abs_2**0.5
        mu = np.abs(stress_ij)/pressure_vector #*shear_rate_abs/np.abs(shear_rate_ij)
        #mu = np.abs(stress_ij)/pressure_vector*shear_rate_abs/np.abs(shear_rate_ij)
        inertia = shear_rate_abs*float(rr.log_variable["dp"])/(pressure_vector/float(rr.log_variable["den"]))**0.5
        # mask too small pressure
        maskoutsmallpressure = (pressure_vector > 5*10**-9)
        #if not np.all(maskoutsmallpressure):
        #    breakpoint()
        mu = mu[maskoutsmallpressure]
        inertia = inertia[maskoutsmallpressure]
        #if np.any(mu>0.3):
        #    breakpoint()
        #if np.any(inertia>1):
        #    breakpoint()  
        return (time, mu, inertia)

    def plot_mu_I_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot,
        i, j,
        k, k_index_array,
        X_scale=float(rr.log_variable['dp']),
        ifgroupk_indexortimeorz="time"
        ):
        fig, ax = plt.subplots()
        if ifgroupk_indexortimeorz=="time":
            for k_index in k_index_array:
                # k_value
                if position_index_to_array_dim_index[k] == 0:
                    k_value = self.chunk_x_array[k_index, 0]
                    
                elif position_index_to_array_dim_index[k] == 1:
                    k_value = self.chunk_y_array[0, k_index]
                
                k_value = k_value/float(rr.log_variable['dp'])

                (time, mu ,inertia) = self.data_to_plot_mu_I(
                        stepsarray[0], coord_index_horizental_in_plot,
                        i, j,
                        k, k_index,
                        X_scale=float(rr.log_variable['dp']),
                    )

                def adddim(array_to_expand):
                    # check array size by using first step
                    onesteparrayshape = array_to_expand.shape
                    # total shape
                    totalshape = (len(stepsarray),) + onesteparrayshape
                    # initialize array
                    array_expanded_empty = np.empty(totalshape)
                    return array_expanded_empty

                time_array = adddim(time)
                mu_array = adddim(mu)
                inertia_array = adddim(inertia)

                for i, step in enumerate(stepsarray):
                    (time, mu ,inertia) = self.data_to_plot_mu_I(
                        step, coord_index_horizental_in_plot,
                        i, j,
                        k, k_index,
                        X_scale=float(rr.log_variable['dp']),
                    )
                    time_array[i] = time 
                    mu_array[i] = mu 
                    inertia_array[i] = inertia 

                n_mu_difftime = time_array.shape[1]   
                for j in range(n_mu_difftime):
                    # plot
                    ax.errorbar(inertia[:,j], mu[:,j], yerr=0,
                        label="t={:.2e}".format(time[0]) + "to {:.2e} s".format(time[-1]) + ", " + map_dim_index_to_coordinate[k] + "=" + "{:.2f}".format(k_value),
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                        )
        elif ifgroupk_indexortimeorz=="z":
            for step in stepsarray:
                for k_index in k_index_array:
                    (time, mu ,inertia) = self.data_to_plot_mu_I(
                        step, coord_index_horizental_in_plot,
                        i, j,
                        k, k_index,
                        X_scale=float(rr.log_variable['dp']),
                    )
                    
                    # k_value
                    if position_index_to_array_dim_index[k] == 0:
                        k_value = self.chunk_x_array[k_index, 0]
                        
                    elif position_index_to_array_dim_index[k] == 1:
                        k_value = self.chunk_y_array[0, k_index]

                    k_value = k_value/float(rr.log_variable['dp'])
                    # plot
                    ax.errorbar(inertia, mu, yerr=0,
                            label="t={:.2e} s".format(time) + ", " + map_dim_index_to_coordinate[k] + "=" + "{:.2f}".format(k_value),
                            marker = ".",
                            linestyle = 'None',
                            markersize=12,
                            )
        else:
            sys.exit("ifgroupk_indexortimeorz not defined well")
        ax.set_xlabel("I")
        ax.set_ylabel("mu")
        ax.set_xscale('log')
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # title
        ax.set_title("# of snapshot {:g}".format(self.n_ave))

        return (fig, ax)
    
    # save X-Y
    def save_mu_I_fix_k(
        self,
        stepsarray,
        coord_index_horizental_in_plot,
        i, j,
        k, k_index_array,
        savepath,
        X_scale=float(rr.log_variable['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):
        self.infofrom_filepath_n_ave_dic_by_stepsarray(stepsarray)
        if onefigure:
            fig, ax = self.plot_mu_I_fix_k(
                                        stepsarray,
                                        coord_index_horizental_in_plot,
                                        i,j,
                                        k, k_index_array,
                                        X_scale,
                                        )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + "-".join([str(number) for number in stepsarray]) + map_dim_index_to_coordinate[k] + "_" + "-".join([str(number) for number in k_index_array]),
                figformat, ifpickle, bbox_inches,
                )
        else:
            for k_index in k_index_array:
                for step in stepsarray:
                    fig, ax = self.plot_mu_I_fix_k(
                                                [step],
                                                coord_index_horizental_in_plot,
                                                i,j,
                                                k, [k_index],
                                                X_scale,
                                                )
            chunk.save_one_plot(
                fig, ax, savepath,
                "step_" + str(int(step)) + map_dim_index_to_coordinate[k] + "_" + str(k_index),
                figformat, ifpickle, bbox_inches,
                )

# chunk kinetic energy
class chunkKE(chunk2D):

    call_header_by_bettername = {
            "Ek_0": "v_Ek1",
            "Ek_1": "v_Ek2",
            "Ek_2": "v_Ek3",
            "mass": "c_m1",
            "Ek_0_sq": "v_Ek1_sq",
            "Ek_1_sq": "v_Ek2_sq",
            "Ek_2_sq": "v_Ek3_sq",
            "mass_sq": "c_m1_sq",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/momentum_mass_field/fix.momentum_mass_field.all")

# chunk omega
class chunkomega(chunk2D):

    call_header_by_bettername = {
            "omega_0": "c_omega[1]",
            "omega_1": "c_omega[2]",
            "omega_2": "c_omega[3]",
            "omega_0_sq": "c_omega[1]_sq",
            "omega_1_sq": "c_omega[2]_sq",
            "omega_2_sq": "c_omega[3]_sq",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/momentum_mass_field/fix.omega.all")
        # stress_variable_name dictionary    

# chunk from many simu

# check 

def checkNchunkint():
    if int(rr.log_variable["freq_ave_chunk_momentum_mass_field"]) != 0:
        ob1 = chunkstress(1, rr.lammps_directory)
        for step_inloop in ob1.allsteps:
            n_line_0 = int(int(step_inloop - ob1.step_first_in_file)/ob1.d_step)*(ob1.n_line_in_a_step+1) + 4
            n_line_1 = int(n_line_0 + ob1.n_line_in_a_step)
            ## select data
            Nchunklist = [ob1.lines[t].split()[ob1.header.index('Ncount')] for t in range(n_line_0, n_line_1)]
            for x in Nchunklist:
                try: 
                    int(x)
                except ValueError:
                    print("one of chunk in step {step_inloop} contain {number} Ncount".format(step_inloop=step_inloop, number=x))
                    if float(x) < 0.1:
                        Nchunkarray = np.asarray(Nchunklist, dtype=np.float64, order='F')
                        totalNchunk = np.sum(Nchunkarray)
                        print(totalNchunk)
                        sys.exit("Ncount (number of atom in chunk) not int")
        print("all chunk has integer number of atoms")


def run_main_by_stepsarray(number_average, stepsarray):
    # execute only if run as a script
    ob3 = chunkmomentum(number_average, rr.lammps_directory)
    toppath = dp.diagram_path + "Ncount/"
    os.makedirs(toppath, exist_ok=True)
    ob3.add_nve_subfolder_in_folder(toppath)
    savepath = ob3.path_nve_subfolder_in_folder(toppath)
    for index_n_2 in range(5):
        ob3.save_time_variable(
            stepsarray, "Ncount", 1,
            "Ncount",
            1, index_n_2,
            savepath
        )

    obmui = chunkmuI(number_average, rr.lammps_directory)
    for (i,j) in [
        (0, 1),
        (0, 2),
        (1, 2),
    ]:
        toppath = dp.f_mu_I_path + str(i) + str(j) + "/"
        os.makedirs(toppath, exist_ok=True)
        obmui.add_nve_subfolder_in_folder(toppath)
        savepath = obmui.path_nve_subfolder_in_folder(toppath)
        obmui.save_mu_I_fix_k(
            stepsarray,
            2,
            i, j,
            1, [0,1,2,5,10],
            savepath,
        )
    # stress
    stress_scale = float(rr.log_variable['den'])*float(rr.log_variable['g'])*float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['dp'])
    bin_volume = float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['bin_y_dp_unit_approximate'])*float(rr.log_variable['bin_z_dp_unit_approximate'])*float(rr.log_variable['dp'])**3
    stress_volume_scale = stress_scale*bin_volume
    if int(rr.log_variable["freq_ave_chunk_momentum_mass_field"]) != 0:
        ob1 = chunkstress(number_average, rr.lammps_directory)
        # plot near vertical wall and bottom
        for (coord_index_horizental_in_plot, k, index_of_k) in [
            (2, 1, 0),
            (2, 1, -1),
            (1, 2, 0),
            ]:
            for key in [key for key in chunkstress.call_header_by_bettername.keys() if "_sq" not in key]:
                toppath = dp.f_stress_field_samescale_path + key + "/"
                os.makedirs(toppath, exist_ok=True)
                ob1.add_nve_subfolder_in_folder(toppath)
                savepath = ob1.path_nve_subfolder_in_folder(toppath)
                ob1.save_XY_fix_k(
                    stepsarray,
                    coord_index_horizental_in_plot, chunkstress.call_header_by_bettername[key],
                    stress_volume_scale,
                    k, [index_of_k],
                    key,
                    savepath,
                )

    # execute only if run as a script
    # chunk wall force
    eachbin2D_on_vertical_area = float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['bin_z_dp_unit_approximate'])*float(rr.log_variable['dp'])**2
    eachbin2D_on_bottom_area = float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['bin_y_dp_unit_approximate'])*float(rr.log_variable['dp'])**2
    
    vertical_wallforce_scale = stress_scale*eachbin2D_on_vertical_area
    horizental_wallforce_scale = stress_scale*eachbin2D_on_bottom_area

    ob1 = chunkinwallforce(number_average, rr.lammps_directory)
    toppath = dp.diagram_path + "stress_ratio/" + key + "/"
    os.makedirs(toppath, exist_ok=True)
    ob1.add_nve_subfolder_in_folder(toppath)
    savepath = ob1.path_nve_subfolder_in_folder(toppath)
    ob1.save_XFricCoeffi_fix_k(
        stepsarray,
        None,
        vertical_wallforce_scale,
        key + "/(bin_area_on_wall)/(stress_scale)",
        savepath,
    )

    for ob1 in [
        chunkinwallforce(number_average, rr.lammps_directory),
        chunkoutwallforce(number_average, rr.lammps_directory),
    ]:
        for key in [key for key in ob1.call_header_by_bettername.keys() if "_sq" not in key]:
            toppath = dp.diagram_path + "wallforce/" + key + "/"
            os.makedirs(toppath, exist_ok=True)
            ob1.add_nve_subfolder_in_folder(toppath)
            savepath = ob1.path_nve_subfolder_in_folder(toppath)
            ob1.save_XY_fix_k(
                stepsarray,
                ob1.call_header_by_bettername[key],
                vertical_wallforce_scale,
                key + "/(bin_area_on_wall)/(stress_scale)",
                savepath,
            )
            ob1.save_time_variable(
                                stepsarray, 
                                ob1.call_header_by_bettername[key],
                                vertical_wallforce_scale,
                                key + "/(bin_area_on_wall)/(stress_scale)",
                                0,
                                savepath,
                                )
        del ob1
    
    ob2 = chunkzbottomwallforce(number_average, rr.lammps_directory)
    for key in [key for key in ob2.call_header_by_bettername.keys() if "_sq" not in key]:
        toppath = dp.diagram_path + "wallforce/" + key + "/"
        os.makedirs(toppath, exist_ok=True)
        ob2.add_nve_subfolder_in_folder(toppath)
        savepath = ob2.path_nve_subfolder_in_folder(toppath)
        ob2.save_XY_fix_k(
            stepsarray,
            ob2.call_header_by_bettername[key],
            horizental_wallforce_scale,
            key + "/(bin_area_on_wall)/(stress_scale)",
            savepath,
        )
        ob2.save_time_variable(
                            stepsarray,
                            ob2.call_header_by_bettername[key],
                            horizental_wallforce_scale,
                            key + "/(bin_area_on_wall)/(stress_scale)",
                            0,
                            savepath,
                            )
    
    # chunkmomentum
    ob3 = chunkmomentum(number_average, rr.lammps_directory)
    velocity_scale = float(rr.log_variable['in_velocity'])
    if velocity_scale < 0:
        velocity_scale = -velocity_scale
    for key in ["mv_0", "mv_1", "mv_2"]:
        toppath = dp.diagram_path + "velocity/" + key + "/"
        os.makedirs(toppath, exist_ok=True)
        ob3.add_nve_subfolder_in_folder(toppath)
        savepath = ob3.path_nve_subfolder_in_folder(toppath)
        ob3.save_XY_fix_k(
            stepsarray,
            2,
            ob3.call_header_by_bettername[key],
            velocity_scale,
            1, [0,-1],
            key,
            savepath,
        )
        
        for index_n_1 in [0, -1]:
            for index_n_2 in range(n_2-3):
                ob3.save_time_variable(
                                    stepsarray,
                                    ob3.call_header_by_bettername[key],
                                    velocity_scale,
                                    key,
                                    index_n_1, index_n_2,
                                    savepath,
                                    )

        for diff_index in [1,2]:
            toppath = dp.diagram_path + "velocity_diff/" + key + "_d" + str(diff_index) + "/"
            os.makedirs(toppath, exist_ok=True)
            ob3.add_nve_subfolder_in_folder(toppath)
            savepath = ob3.path_nve_subfolder_in_folder(toppath)
            for k_index in [0,2,5,10,-1]:
                ob3.save_XY_shear_rate_fix_k(
                    stepsarray,
                    2,
                    ob3.call_header_by_bettername[key],
                    velocity_scale,
                    1, [k_index],
                    key,
                    diff_index,
                    savepath,
                )

def run_main_first_middle_last_steps(number_average):
    # execute only if run as a script
    
    ob1 = chunkstress(number_average, rr.lammps_directory)
    stepsarray = ob1.first_middle_last_steps
    run_main_by_stepsarray(number_average, stepsarray)

def run_main_by_step(number_average, step_begin, step_end_not_include, delta_step):
    stepsarray = np.arange(step_begin, step_end_not_include, delta_step)
    run_main_by_stepsarray(number_average, stepsarray)

def run_main_all_steps(number_average):
    # execute only if run as a script
    
    # stress
    stress_scale = float(rr.log_variable['den'])*float(rr.log_variable['g'])*float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['dp'])
    bin_volume = float(rr.log_variable['width_wall_dp_unit'])*float(rr.log_variable['bin_y_dp_unit_approximate'])*float(rr.log_variable['bin_z_dp_unit_approximate'])*float(rr.log_variable['dp'])**3
    stress_volume_scale = stress_scale*bin_volume
    if int(rr.log_variable["freq_ave_chunk_momentum_mass_field"]) != 0:
        ob1 = chunkstress(number_average, rr.lammps_directory)
        # plot near vertical wall and bottom
        for index_n_1 in [0, -1]:
            for index_n_2 in range(n_2-3):
                ob2 = chunkstress(1, rr.lammps_directory)
                for key in [key for key in chunkstress.call_header_by_bettername.keys() if "_sq" not in key]:
                    toppath = dp.f_stress_field_samescale_path + key + "_time/"
                    os.makedirs(toppath, exist_ok=True)
                    ob2.add_nve_subfolder_in_folder(toppath)
                    savepath = ob2.path_nve_subfolder_in_folder(toppath)

                    ob2.save_time_variable(
                        ob2.allsteps, chunkstress.call_header_by_bettername[key],
                        stress_scale,
                        key,
                        index_n_1, index_n_2,
                        savepath,
                    )

# main exclusive
if __name__ == "__main__":
    # input as n_ave, step1, step2 (not include), d_step. 
    #       or n_ave, 'all'
    #       or n_ave, 'middle'
    number_average = int(sys.argv[1])
    if len(sys.argv) == 3:
        if sys.argv[2] == 'all':
            run_main_all_steps(number_average)
        elif sys.argv[2] == 'middle':
            run_main_first_middle_last_steps(number_average)
    elif len(sys.argv) == 5:
        run_main_by_step(number_average, int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    else:
        sys.exit("len(sys.argv) is not correct")
