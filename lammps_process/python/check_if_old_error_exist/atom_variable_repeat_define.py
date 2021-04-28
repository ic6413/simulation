import sys
import os

# check if log_variable output has atom variable define two times
def if_atom_variable_define_notonly_1_times(log_variablepath):
    if os.path.isfile(log_variablepath):
        with open(log_variablepath, mode='r') as f:
            lines = f.read().strip().split('\n')
    else:
        sys.exit("log file not exist")
    # variable lines
    lines_start_variable = [line for line in lines if line.startswith("variable")]
    variable_names = [line.split()[1] for line in lines_start_variable]
    variable_names_set = set(variable_names)
    for v_name in variable_names_set:
        satisfy_lines = [line for line in lines_start_variable if line.split()[1] == v_name]
        if len(satisfy_lines) >= 2:
            # check if satisfy lines contain atom variable
            for satisfy_line in satisfy_lines:
                if satisfy_line.split()[2] == "atom":
                    define_2_times_ok_list = [
                        'check_nbid1',
                        'check_nbid2',
                        'check_nbid3',
                        'check_nbmaxKEt',
                        'check_nbmaxKEr',
                        'check_nbmaxKEtr',
                    ]
                    if v_name not in define_2_times_ok_list:
                        sys.exit("error: atom variable {v_name} define two times".format(v_name=v_name))
                    else:
                        print("atom variable {v_name} define two times".format(v_name=v_name))
# main exclusive
if __name__ == "__main__":
    lammps_directory = os.getcwd() + '/'
    if_atom_variable_define_notonly_1_times(lammps_directory + '/log.lammps')
