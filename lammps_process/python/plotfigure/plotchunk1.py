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
if rr.logfile["shearwall"] == "zcylinder":
    chunk_method = 'rz'
if rr.logfile["shearwall"] == "yplane":
    chunk_method = 'yz'

# map dim index to coordinate
if rr.logfile["shearwall"] == "zcylinder":
    map_dim_index_to_coordinate = ["t", "r", "z"]
elif rr.logfile["shearwall"] == "yplane":
    map_dim_index_to_coordinate = ["x", "y", "z"]

# xyztoCoor
if rr.logfile["shearwall"] == "zcylinder":
    # chink first dim unchange in the begining
    chunk_first_dim_coord = "z"
    chunk_second_dim_coord = "r"
elif rr.logfile["shearwall"] == "yplane":
    if "chunk/atom 23" in rr.logfile.keys():
        if rr.logfile["chunk/atom 23"][1] == "y":
            chunk_first_dim_coord = "y"
            chunk_second_dim_coord = "z"
            xyztoCoor = {}
            xyztoCoor["y"] = "Coord1"
            xyztoCoor["z"] = "Coord2"
        elif rr.logfile["chunk/atom 23"][1] == "z":
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
if rr.logfile["shearwall"] == "zcylinder":
    position_index_to_array_dim_index = {
                                    1: 1,
                                    2: 0,
                                    }
    n_r = int(rr.logfile['N_bin_r'])
    n_z = int(rr.logfile['N_bin_z'])
    n_1 = n_z
    n_2 = n_r
    n_12 = n_1*n_2
    
    dx = 1/n_r*int(rr.logfile['width_wall_dp_unit'])
    dy = 1/n_z*float(rr.logfile['zhi_chunk_dp_unit'])
    x_array, y_array = np.meshgrid(
                                int(rr.logfile['ri_wall_dp_unit']) + (np.arange(n_1)+0.5)/n_1*int(rr.logfile['width_wall_dp_unit']),
                                (np.arange(n_2)+0.5)/n_2*float(rr.logfile['zhi_chunk_dp_unit']),
                                )
    x_array = x_array.reshape((-1))
    y_array = y_array.reshape((-1))
    vol_in_chunks = np.pi*((x_array+0.5*dx)**2-(x_array-0.5*dx)**2)*(y_array+0.5*dy-(y_array-0.5*dy))*float(rr.logfile['dp'])**3
elif rr.logfile["shearwall"] == "yplane":
    n_y = int(rr.logfile['N_bin_y'])
    n_z = int(rr.logfile['N_bin_z'])
    if "chunk/atom 23" in rr.logfile.keys():
        if rr.logfile["chunk/atom 23"][1] == "y":
            position_index_to_array_dim_index = {
                                            1: 0,
                                            2: 1,
                                            }
            n_1 = n_y
            n_2 = n_z
        elif rr.logfile["chunk/atom 23"][1] == "z":
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
    dx = 1/n_y*int(rr.logfile['width_wall_dp_unit'])
    dy = 1/n_z*float(rr.logfile['zhi_chunk_dp_unit'])
    vol_in_chunks = float(rr.logfile['dp'])*int(rr.logfile['x_period_dp_unit'])*dx*dy*float(rr.logfile['dp'])**2
else:
    sys.exit("chunk_method wrong")
# time count from step 0
def time_from_step_0(step):
    return step*float(rr.logfile["ts"])

# time count from rotate started
def time_from_start_rotate(step):
    return time_from_step_0(step)-rr.logfile["rotate_start_time"]

# class chunk common
class chunk(object):

    def __init__(self, n_ave, lmp_path, f_path_rela_lmpfolder):
        self.n_ave = n_ave
        self.lmp_path = lmp_path
        self.f_path_rela_lmpfolder = f_path_rela_lmpfolder
        # get lines from files
        with open(self.lmp_path + f_path_rela_lmpfolder) as f:
            self.lines = f.read().strip().split('\n')
        self.logfile = rr.logfile
        self.header = self.lines[2].split()[1:]
        self.n_line_in_a_step = int(self.lines[3].split()[1])
        self.step_first_in_file = int(self.lines[3].split()[0])
        self.step_second_in_file = int(self.lines[3 + self.n_line_in_a_step + 1].split()[0])
        self.step_last_in_file = int(self.lines[-1 - self.n_line_in_a_step].split()[0])
        self.d_step = self.step_second_in_file - self.step_first_in_file
        self.step_first_in_file_change_by_n_ave = self.step_first_in_file + int((self.n_ave)/2*self.d_step)
        self.step_last_in_file_change_by_n_ave = self.step_last_in_file - int((self.n_ave)/2*self.d_step)
        self.middle_step = int(
            int((self.step_last_in_file_change_by_n_ave - self.step_first_in_file_change_by_n_ave)/2/self.d_step)*self.d_step + self.step_first_in_file_change_by_n_ave
            )
        self.allsteps = np.arange(self.step_first_in_file_change_by_n_ave, self.step_last_in_file_change_by_n_ave, self.d_step)
        self.first_middle_last_steps = np.array([self.step_first_in_file_change_by_n_ave, self.middle_step, self.step_last_in_file_change_by_n_ave])
        self.extrasteps = np.array(
            [
                self.step_first_in_file_change_by_n_ave + 5*self.d_step*self.n_ave,
                ]
            )
        maskextra = np.logical_and(self.extrasteps > self.step_first_in_file_change_by_n_ave, self.extrasteps < self.step_last_in_file_change_by_n_ave)
        self.extrasteps = self.extrasteps[maskextra]
        self.first_extra_middle_last_steps = np.append(self.first_middle_last_steps, self.extrasteps)
        self.first_extra_middle_last_steps.sort()
        if "if_inwall_wall_gran" in rr.logfile.keys():
            if rr.logfile["if_inwall_wall_gran"] == "yes":
                if "wall_gran_type" in rr.logfile.keys():
                    if rr.logfile["wall_gran_type"] == "1":
                        self.ybottomwalltype = "rough (d=0.9)"
                    elif rr.logfile["wall_gran_type"] == "2":
                        self.ybottomwalltype = "rough (d=1)"
                    elif rr.logfile["wall_gran_type"] == "3":
                        self.ybottomwalltype = "rough (d=1.1)"
                    else:
                        sys.exit("can not get wall gran type")
                else:
                    self.ybottomwalltype = "rough (d=1)"
            else:
                self.ybottomwalltype = "smooth"
        else:
            self.ybottomwalltype = "smooth"

        self.height = rr.logfile["z_length_create_dp_unit"]
        self.width = rr.logfile["width_wall_dp_unit"]
        self.periodlength = rr.logfile["x_period_dp_unit"]
        self.labelstring_size_walltype = self.ybottomwalltype + "\n" + "L " + self.periodlength + "\n" + "W " + self.width + "\n" + "H " + self.height
        self.labelstring_size_walltype_one_line = self.ybottomwalltype + ", " + "L " + self.periodlength + ", " + "W " + self.width + ", " + "H " + self.height

        

    def firstdata(self, variable_name):
        pass
    def value_in_a_step_ave(self, step, variable_name, index_n_1=None, index_n_2=None, ddof=1):
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
    
        if rr.logfile["shearwall"] == "zcylinder":
            plt.xlabel('r')
            plt.ylabel('z')
        elif rr.logfile["shearwall"] == "yplane":
            plt.xlabel('y')
            plt.ylabel('z')
        
        return (fig, ax)

class chunk1D(chunk):
    def __init__(self, n_ave, lmp_path, f_path_rela_lmpfolder):
        super().__init__(n_ave, lmp_path, f_path_rela_lmpfolder)
        # coordinate in chunk
        if rr.logfile["shearwall"] == "zcylinder":
            self.chunk_coor_array = self.firstdata("v_z")
            self.chunk_coor_array_label = "z"
        elif rr.logfile["shearwall"] == "yplane":
            self.chunk_coor_array = self.firstdata("Coord2")
            self.chunk_coor_array_label = "z"
        else:
            sys.exit("chunk_method wrong")

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
        # sum_square_2 = sum of (xi-x_ave)**2
        sum_square_2 = 0

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
                sum_square_2 += (x-value)*(x-oldvalue)

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
                sum_square_2 += (x-value)*(x-oldvalue)
        
        if (ddof == 1 and self.n_ave == 1):
            std = 0
        else:
            std = (
                sum_square_2/(self.n_ave-ddof)
            )**0.5
        return (value, std)
    # data 1D-1D
    def datachunk_ave_one_step_XY(
        self,
        step, Y_name,
        Y_scale,
        X_scale=float(rr.logfile['dp']),
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
        X_scale=float(rr.logfile['dp']),
        ):
        
        fig, ax = plt.subplots()
        for step in stepsarray:
            (time, X_vector, Xlabel, Y_vector, Y_std_vector) = (
                self.datachunk_ave_one_step_XY(
                    step, Y_name,
                    Y_scale,
                    X_scale=float(rr.logfile['dp']),
                )
            )

            # plot
            ax.errorbar(X_vector, Y_vector, yerr=Y_std_vector,
                    label="t={:.2e} s".format(time),
                    marker = ".",
                    linestyle = 'None',
                    markersize=12,
                    )
        
        ax.set_xlabel(Xlabel)
        ax.set_ylabel(Ylabel)
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
        X_scale=float(rr.logfile['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):

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
        value_coor_dpunit = self.chunk_coor_array[index_n]/float(rr.logfile['dp'])
        # plot
        ax.errorbar(time, Y_vector, yerr=Y_std_vector,
                label=self.chunk_coor_array_label + "={:.1f}, ".format(value_coor_dpunit),
                marker = ".",
                linestyle = 'None',
                markersize=12,
                )
        
        ax.set_xlabel("time")
        ax.set_ylabel(Ylabel)
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
    call_header_by_bettername = {
            "F_inwall_0": "v_inwall_per_atom_1", 
            "F_inwall_1": "v_inwall_per_atom_2", 
            "F_inwall_2": "v_inwall_per_atom_3",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/wall/chunk/inwallforcefile")
        # coordinate in chunk
        if rr.logfile["shearwall"] == "zcylinder":
            self.chunk_coor_array = super().firstdata("v_z")
            self.chunk_coor_array_label = "z"
        elif rr.logfile["shearwall"] == "yplane":
            self.chunk_coor_array = super().firstdata("Coord2")
            self.chunk_coor_array_label = "z"
        else:
            sys.exit("chunk_method wrong")

class chunkoutwallforce(chunk1D):
    call_header_by_bettername = {
            "F_outwall_0": "v_outwall_per_atom_1", 
            "F_outwall_1": "v_outwall_per_atom_2", 
            "F_outwall_2": "v_outwall_per_atom_3",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/wall/chunk/outwallforcefile")
        # coordinate in chunk
        if rr.logfile["shearwall"] == "zcylinder":
            self.chunk_coor_array = super().firstdata("v_z")
            self.chunk_coor_array_label = "z"
        elif rr.logfile["shearwall"] == "yplane":
            self.chunk_coor_array = super().firstdata("Coord2")
            self.chunk_coor_array_label = "z"
        else:
            sys.exit("chunk_method wrong")
        
class chunkzbottomwallforce(chunk1D):
    call_header_by_bettername = {
            "F_zbottom_0": "v_zbottom_per_atom_1", 
            "F_zbottom_1": "v_zbottom_per_atom_2", 
            "F_zbottom_2": "v_zbottom_per_atom_3",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/wall/chunk/zbottomforcefile")
        # coordinate in chunk
        if rr.logfile["shearwall"] == "zcylinder":
            self.chunk_coor_array = super().firstdata("v_r")
            self.chunk_coor_array_label = "r"
        elif rr.logfile["shearwall"] == "yplane":
            self.chunk_coor_array = super().firstdata("Coord1")
            self.chunk_coor_array_label = "y"
        else:
            sys.exit("chunk_method wrong")


class chunk2D(chunk):

    def __init__(self, n_ave, lmp_path, f_path_rela_lmpfolder):
        super().__init__(n_ave, lmp_path, f_path_rela_lmpfolder)
        # x-y coordinate in chunk
        if rr.logfile["shearwall"] == "zcylinder":
            self.chunk_x_array = self.firstdata("v_r")
            self.chunk_y_array = self.firstdata("v_z")
            self.chunk_x_array_label = "r"
            self.chunk_y_array_label = "z"
        elif rr.logfile["shearwall"] == "yplane":
            self.chunk_x_array = self.firstdata("Coord1")
            self.chunk_y_array = self.firstdata("Coord2")
            self.chunk_x_array_label = "y"
            self.chunk_y_array_label = "z"
        else:
            sys.exit("chunk_method wrong")

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
        # sum_square_2 = sum of (xi-x_ave)**2
        sum_square_2 = 0

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
                    sum_square_2 += (x-value)*(x-oldvalue)
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
                    sum_square_2 += (x-value)*(x-oldvalue)
                    
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
                    sum_square_2 += (x-value)*(x-oldvalue)
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
                    sum_square_2 += (x-value)*(x-oldvalue)

        if (ddof == 1 and self.n_ave == 1):
            std = 0
        else:
            std = (
                sum_square_2/(self.n_ave-ddof)
            )**0.5
        return (value, std)

    # data quiver 2D 2D
    def datachunk_ave_one_step_quiver_x23(
        self,
        step, Q_name, V_name, Q_scale, V_scale,
        x_scale=float(rr.logfile['dp']), y_scale=float(rr.logfile['dp']),
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
    def plotquiver_x23(
        self,
        step, Q_name, V_name, Q_scale, V_scale,
        quiver_scale=0.1, label_scale=0.2,
        x_scale=float(rr.logfile['dp']), y_scale=float(rr.logfile['dp']),
        ):

        [time, x_array, y_array, Q_array, V_array] = (
            self.datachunk_ave_one_step_quiver_x23(
            step, Q_name, V_name, Q_scale, V_scale,
            x_scale, y_scale,
            )
        )

        fig1, ax1 = plt.subplots()
        #fig1.figsize = [12.8, 9.6]
        chunk.plot_quiver_position_label(fig1, ax1)
        #ax1.set_title('velocity field r-z direction (average over theta)')
        Q = ax1.quiver(x_array, y_array, Q_array, V_array,
                    units='width',angles='xy', scale_units='xy', scale=quiver_scale,
                    )

        ax1.quiverkey(Q, 0.2, 0.95, label_scale,
                    label = "equal stress = {:.2e}".format(label_scale) + ". At {:.2e} s".format(time),
                    labelpos='E',
                    coordinates='figure', angle=90)
        return (fig1, ax1)
    
    # save quiver
    def savequiver_x23(
        self,
        stepsarray,
        Q_name, V_name, Q_scale, V_scale,
        savepath,
        x_scale=float(rr.logfile['dp']), y_scale=float(rr.logfile['dp']),
        figformat="png", ifpickle=False,
        quiver_scale=0.1, label_scale=0.2,
        ):
        
        for step in stepsarray:
            fig1, ax1 = (
                self.plotquiver_x23(
                    step, Q_name, V_name, Q_scale, V_scale,
                    quiver_scale, label_scale,
                    x_scale, y_scale,
                )
            )
            chunk.save_one_plot(fig1, ax1, savepath, str(int(step)), figformat, ifpickle)

    # data 1D-1D
    def datachunk_ave_one_step_XY_fix_k(
        self,
        step, coord_index_horizental_in_plot, Y_name,
        Y_scale,
        k, k_index,
        X_scale=float(rr.logfile['dp']),
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
        X_scale=float(rr.logfile['dp']),
        ):
        
        fig, ax = plt.subplots()
        for k_index in k_index_array:
            for step in stepsarray:

                (time, X_vector, Xlabel, Y_vector, Y_std_vector) = (
                    self.datachunk_ave_one_step_XY_fix_k(
                        step, coord_index_horizental_in_plot, Y_name,
                        Y_scale,
                        k, k_index,
                        X_scale=float(rr.logfile['dp']),
                    )
                )
                
                # k_value
                if position_index_to_array_dim_index[k] == 0:
                    k_value = self.value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], index_n_1=k_index, index_n_2=0)[0]
                    
                elif position_index_to_array_dim_index[k] == 1:
                    k_value = self.value_in_a_step_ave(step, xyztoCoor[map_dim_index_to_coordinate[k]], index_n_1=0, index_n_2=k_index)[0]

                k_value = k_value/float(rr.logfile['dp'])
                # plot
                ax.errorbar(X_vector, Y_vector, yerr=Y_std_vector,
                        label="t={:.2e} s".format(time) + ", " + map_dim_index_to_coordinate[k] + "=" + "{:.2f}".format(k_value),
                        marker = ".",
                        linestyle = 'None',
                        markersize=12,
                        )
        
        ax.set_xlabel(Xlabel)
        ax.set_ylabel(Ylabel)
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
        X_scale=float(rr.logfile['dp']),
        figformat="png", ifpickle=False, bbox_inches="tight", onefigure=True
        ):

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
        value_x_dpunit = self.chunk_x_array[index_n_1, index_n_2]/float(rr.logfile['dp'])
        value_y_dpunit = self.chunk_y_array[index_n_1, index_n_2]/float(rr.logfile['dp'])
        # plot
        ax.errorbar(time, Y_vector, yerr=Y_std_vector,
                label=self.chunk_x_array_label + "={:.1f}, ".format(value_x_dpunit) + self.chunk_y_array_label + "={:.1f}".format(value_y_dpunit),
                marker = ".",
                linestyle = 'None',
                markersize=12,
                )
        
        ax.set_xlabel("time")
        ax.set_ylabel(Ylabel)
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
            "F_inwall_0": "v_inwall_per_atom_1", 
            "F_inwall_1": "v_inwall_per_atom_2", 
            "F_inwall_2": "v_inwall_per_atom_3", 
            "F_outwall_0": "v_outwall_per_atom_1", 
            "F_outwall_1": "v_outwall_per_atom_2", 
            "F_outwall_2": "v_outwall_per_atom_3", 
            "F_zbottom_0": "v_zbottom_per_atom_1", 
            "F_zbottom_1": "v_zbottom_per_atom_2", 
            "F_zbottom_2": "v_zbottom_per_atom_3",
            }

    def __init__(self, n_ave, lmp_path):
        super().__init__(n_ave, lmp_path, "output/stress/fix.stress.all")
        # stress_variable_name dictionary
        

# chunk momentum
class chunkmomentum(chunk2D):

    def __init__(self, n_ave, lmp_path, f_path_rela_lmpfolder):
        super().__init__(n_ave, lmp_path, "output/momentum_mass_field/fix.momentum_mass_field.all")

    
# chunk from many simu



def run_main(number_average):
    # execute only if run as a script
    
    stress_scale = float(rr.logfile['den'])*float(rr.logfile['g'])*float(rr.logfile['zhi_chunk_dp_unit'])*float(rr.logfile['dp'])
    

    
    ob1 = chunkstress(number_average, rr.lammps_directory)
    for key in chunkstress.call_header_by_bettername.keys():

        toppath = dp.f_stress_field_samescale_path + key + "/"
        os.makedirs(toppath, exist_ok=True)
        ob1.add_nve_subfolder_in_folder(toppath)
        savepath = ob1.path_nve_subfolder_in_folder(toppath)
        ob1.save_XY_fix_k(
            ob1.first_middle_last_steps,
            2, chunkstress.call_header_by_bettername[key],
            stress_scale,
            1, [0,-1],
            key,
            savepath,
        )

    

    for index_n_1 in [0, -1]:
        for index_n_2 in range(n_2-3):
            ob2 = chunkstress(1, rr.lammps_directory)
            for key in chunkstress.call_header_by_bettername.keys():
                toppath = dp.f_stress_field_samescale_path + key + "_time/"
                os.makedirs(toppath, exist_ok=True)
                ob2.add_nve_subfolder_in_folder(toppath)
                savepath = ob2.path_nve_subfolder_in_folder(toppath)

                ob2.save_time_variable(
                    ob2.allsteps, chunkstress.call_header_by_bettername[key],
                    stress_scale,
                    key,
                    index_n_1, index_n_2,
                    savepath
                )

    # execute only if run as a script
    stress_scale = float(rr.logfile['den'])*float(rr.logfile['g'])*float(rr.logfile['zhi_chunk_dp_unit'])*float(rr.logfile['dp'])
    real_scale = stress_scale*float(rr.logfile['x_period_dp_unit'])*float(rr.logfile['bin_z_dp_unit_approximate'])*float(rr.logfile['dp'])**2

    for ob1 in [
        chunkinwallforce(number_average, rr.lammps_directory),
        chunkoutwallforce(number_average, rr.lammps_directory),
    ]:
        for key in ob1.call_header_by_bettername.keys():
            toppath = dp.diagram_path + "wallforce/" + key + "/"
            os.makedirs(toppath, exist_ok=True)
            ob1.add_nve_subfolder_in_folder(toppath)
            savepath = ob1.path_nve_subfolder_in_folder(toppath)
            ob1.save_XY_fix_k(
                ob1.first_middle_last_steps,
                ob1.call_header_by_bettername[key],
                real_scale,
                key,
                savepath,
            )
            stepsarray = np.arange(ob1.first_middle_last_steps[0], ob1.first_middle_last_steps[0] + 5*number_average, number_average)
            ob1.save_time_variable(
                                stepsarray, 
                                ob1.call_header_by_bettername[key],
                                real_scale,
                                key,
                                0,
                                savepath,
                                )
        del ob1

    ob2 = chunkzbottomwallforce(number_average, rr.lammps_directory)
    for key in ob2.call_header_by_bettername.keys():
        toppath = dp.diagram_path + "wallforce/" + key + "/"
        os.makedirs(toppath, exist_ok=True)
        ob2.add_nve_subfolder_in_folder(toppath)
        savepath = ob2.path_nve_subfolder_in_folder(toppath)
        ob2.save_XY_fix_k(
            ob2.first_middle_last_steps,
            ob2.call_header_by_bettername[key],
            real_scale,
            key,
            savepath,
        )
        stepsarray = np.arange(ob2.first_middle_last_steps[0], ob2.first_middle_last_steps[0] + 5*number_average, number_average)
        ob2.save_time_variable(
                            stepsarray, 
                            ob2.call_header_by_bettername[key],
                            real_scale,
                            key,
                            0,
                            savepath,
                            )
# main exclusive
if __name__ == "__main__":
    
    run_main(5000)
