import sys

import read_setting as rr

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

def data_from_multi_simu(classobject, datamethodtouse, *args, **kwargs):
    classobject.datamethodtouse(step, *args, **kwargs)


class chunkmulti(object):
    def __init__(self, n_ave, lmp_path_list, f_path_rela_lmpfolder):
        self.n_ave = n_ave
        self.lmp_path_list = lmp_path_list
        self.f_path_rela_lmpfolder = f_path_rela_lmpfolder
        
        
        # get lines from files
        for i in rr.n_simu_total:
        with open(self.lmp_path + f_path_rela_lmpfolder) as f:
            self.lines = f.read().strip().split('\n')
            header = self.lines[2].split()[1:]

        # check if header same, 
        # if header not the same, output simu_index, and do not append line
        # if header the same, combine lines

        self.log_variable = rr.log_variable
        self.header = self.lines[2].split()[1:]
        self.n_line_in_a_step = int(self.lines[3].split()[1])
        self.step_first_in_file = int(self.lines[3].split()[0])
        self.step_second_in_file = int(self.lines[3 + self.n_line_in_a_step + 1].split()[0])
        self.step_last_in_file = int(self.lines[-1 - self.n_line_in_a_step].split()[0])
        self.d_step = self.step_second_in_file - self.step_first_in_file
        self.step_first_in_file_change_by_n_ave = self.step_first_in_file + (self.n_ave-1)/2*self.d_step
        self.step_last_in_file_change_by_n_ave = self.step_last_in_file - (self.n_ave-1)/2*self.d_step
        if self.step_first_in_file_change_by_n_ave > self.step_last_in_file_change_by_n_ave:
            sys.exit("error: step_first_in_file_change_by_n_ave > step_last_in_file_change_by_n_ave, n_ave too large") 
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
        if "if_inwall_wall_gran" in rr.log_variable.keys():
            if rr.log_variable["if_inwall_wall_gran"] == "yes":
                if "wall_gran_type" in rr.log_variable.keys():
                    if rr.log_variable["wall_gran_type"] == "1":
                        self.ybottomwalltype = "rough (d=0.9)"
                    elif rr.log_variable["wall_gran_type"] == "2":
                        self.ybottomwalltype = "rough (d=1)"
                    elif rr.log_variable["wall_gran_type"] == "3":
                        self.ybottomwalltype = "rough (d=1.1)"
                    else:
                        sys.exit("can not get wall gran type")
                else:
                    self.ybottomwalltype = "rough (d=1)"
            else:
                self.ybottomwalltype = "smooth"
        else:
            self.ybottomwalltype = "smooth"

        self.height = rr.log_variable["z_length_create_dp_unit"]
        self.width = rr.log_variable["width_wall_dp_unit"]
        self.periodlength = rr.log_variable["x_period_dp_unit"]
        self.labelstring_size_walltype = self.ybottomwalltype + "\n" + "L " + self.periodlength + "\n" + "W " + self.width + "\n" + "H " + self.height
        self.labelstring_size_walltype_one_line = self.ybottomwalltype + ", " + "L " + self.periodlength + ", " + "W " + self.width + ", " + "H " + self.height
    