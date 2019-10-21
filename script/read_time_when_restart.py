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
    
    for n_line1, line in enumerate(lines):
        if len(line) >= 2:
            if line.split()[0]=="Step" and line.split()[1]=="Time":
                break
    
    for n_line2, line in enumerate(lines):
        if len(line) >= 2:
            if line.split()[0]=="Loop" and line.split()[1]=="time":           
                break

    for n, line in enumerate(lines):
        if n>n_line1 and n<n_line2:
            if int(line.split()[0])==rst_from_current:
                restart_time = line.split()[1]
                break

    for n, line in enumerate(lines):
        if len(line) >= 2:
            if line.startswith("variable") and line.split()[1] == "rst_from" and line.split()[2] == "index":
                last_rst_from = int(line.split()[3])
                break

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

print(get_restart_time())
