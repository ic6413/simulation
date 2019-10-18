#!/usr/bin/env python
import sys
import os

import read_setting.read_setting as rr
import datapath as dp

rst_from_current = int(rr.logfile['rst_from'])

# === current module inputvariable ===

def restart_time_and_last_rst_from(currentfolderpath):
    last_folder_path = currentfolderpath + '../'
    last_log_path = currentfolderpath + '../log.lammps'

    if os.path.isfile(last_log_path):
            
        with open(last_log_path, mode='r') as f:
            lines = f.read().strip().split('\n')

    else:
        sys.exit("file not exist")

    n_lines_Step_Time = [(n, line) for n, line in enumerate(lines) if line.split()[0]=="Step" and line.split()[1]=="Time"]
    (n_line1, line1) = n_lines_Step_Time[0]

    n_line_Loop_time = [(n, line) for n, line in enumerate(lines) if line.split()[0]=="Loop" and line.split()[1]=="time"]
    (n_line2, line2) = n_line_Loop_time[0]

    n_lines_pick = [(n, line) for n, line in enumerate(lines) if n>n_line1 and n<n_line2 and int(line.split()[0])==rst_from_current]
    (n_line3, line3) = n_lines_pick[0]

    restart_time = line3[1]

    last_rst_from = int([line.split()[3] for line in lines if line.startswith("variable") and line.split()[1] == "rst_from" and line.split()[2] == "index"][0])

    return [restart_time, last_rst_from, last_folder_path]

def get_restart_time():
    total_previous_time = 0

    if rst_from_current != 0:

        [restart_time, last_rst_from, last_folder_path] = restart_time_and_last_rst_from(dp.lammps_directory)
        total_previous_time += float(restart_time)

        while last_rst_from != 0:
            [restart_time, last_rst_from, last_folder_path] = restart_time_and_last_rst_from(last_folder_path)
            total_previous_time += float(restart_time)
 
    return total_previous_time