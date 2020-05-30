import os
import check_if_old_error_exist.atom_variable_repeat_define as ca

# main exclusive
if __name__ == "__main__":
    lammps_directory = os.getcwd() + '/'
    ca.if_atom_variable_define_notonly_1_times(lammps_directory + '/log.lammps')