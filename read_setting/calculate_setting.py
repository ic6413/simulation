import sys
import os
import read_setting.read_setting as rr



lammps_directory = os.getcwd() + '/'  #os.path.expanduser('~/lammps_simulation/run/0103_mu05_nofreeze_traceid15556/')

n_log_list = len(rr.logfilelist_from_initial_to_lastest)

log_path_list_last_to_initial = [lammps_directory+"../"*n+'log.lammps' for n in range(n_log_list)]
log_path_list_initial_to_last = [log_path_list_last_to_initial[n_log_list-1-n] for n in range(n_log_list)]
calculate_setting_diclist_from_initial_to_last = []

folder_path_list_initial_to_last = [lammps_directory+"../"*(n_log_list-1-n) for n in range(n_log_list)]
folder_path_list_last_to_initial = [folder_path_list_initial_to_last[n_log_list-1-n] for n in range(n_log_list)]

def count_restart_time(index):
    restart_time = 0
    for n in range(index+1):
        if n == 0:
            restart_time_from_last_to_current = 0
        else:
            logfile = rr.logfilelist_from_initial_to_lastest[n]
            rst_from_current = str(logfile['rst_from'])
            last_log_path = log_path_list_initial_to_last[n-1]
        
            if os.path.isfile(last_log_path):
                with open(last_log_path, mode='r') as f:
                    lines = f.read().strip().split('\n')
            else:
                sys.exit("file not exist")

            for n_line1, line in enumerate(lines):
                if len(line) >= 2:
                    if line.split()[0]=="Step" and line.split()[1]=="Time":
                        break
            
            for n_line2, line in enumerate(lines):
                if len(line) >= 2:
                    if line.split()[0]=="Loop" and line.split()[1]=="time":           
                        break

            for n_line, line in enumerate(lines):
                if n_line>n_line1 and n_line<n_line2:
                    if line.split()[0]==rst_from_current:
                        restart_time_from_last_to_current = float(line.split()[1])                        
                        break
        restart_time += restart_time_from_last_to_current
    return restart_time


def calculate_setting_dic(index):
    dic = {}
    dic["previous_time"] = count_restart_time(index)
    return dic

for index in range(n_log_list):
    calculate_setting_diclist_from_initial_to_last.append(calculate_setting_dic(index))

calculate_setting_dic = calculate_setting_diclist_from_initial_to_last[-1]