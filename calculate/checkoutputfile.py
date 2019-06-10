import sys
import os.path
import pprint
import datetime
import numpy as np
# import module in simulation folder
import osmanage as om
import datapath as dp
# import module in calculate folder
import calculate.checksimulation as cs

# input
f_custom_path = dp.f_custom
debug_print_path = dp.debug_print_path

class checkfile(object):

    def __init__(self, step1, step2):
        self.step1 = step1
        self.step2 = step2
        # begin
        om.create_directory(dp.debug_print_path)
        np.set_printoptions(threshold = 100000)
        np.set_printoptions(precision=6)

    def f_name(self):
        pass

    def filepath(self):
        return debug_print_path + self.f_name()

    def stringstep(self):
        return "step_{first}_{second}".format(first=self.step1, second=self.step2)

    def string_wallname_id(self):
        return "(wall_id, wallname) = {first}".format(first=cs.wall_list_id_name)

    def print_string_list(self, string_list):
        pprint.pprint ('\n'.join(string_list), open(self.filepath(), 'a'))


class checkfile_look_idi(checkfile):

    def __init__(self, id_i, step1, step2):
        super().__init__(step1, step2)
        self.id_i = id_i

    def f_name(self):
        pass

    def string_idi_step(self):
        return 'idi_' + str(self.id_i) + '_' + self.stringstep()

class checkforce(checkfile_look_idi):

    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2)
        self.error_tolerence = error_tolerence
        self.method_list = method_list # method = 0, 1, 2 are correct. method = 3~7(8) wrong
        print("checking force")

    def f_name(self):
        return 'check_force_' + self.string_idi_step()

    def string_errortolerence(self):
        return "print error larger than error_tolerence = {error_tolerence}".format(error_tolerence=self.error_tolerence)

    def checkprint(self):

        file = self.filepath()

        if os.path.isfile(file):
            print ('file exist, calculate again and append to the end')
        else:
            print ("create file and print to file")


        check_f_function = cs.fjwi_plus_check_multistep_multicontact_fromcustom  # or cs.fjwi_plus_check_multistep_multicontact_fromcustom_v1

        self.print_string_list([
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.string_errortolerence(), self.stringstep(),
        ])

        for method_i in self.method_list:

            [f_step_error_array, fi_cal_in_error_step, fi_plus_in_error_step] = check_f_function(f_custom_path, self.id_i, self.step1, self.step2, self.error_tolerence, method=method_i)
            self.print_string_list([
                "method {method_i}".format(method_i=method_i),
                "f_step_error_array",
            ])
            pprint.pprint(f_step_error_array, open(file, 'a'))
            pprint.pprint("fi_plus_in_error_step", open(file, 'a'))
            pprint.pprint(fi_plus_in_error_step, open(file, 'a'))
            pprint.pprint("fi_cal_in_error_step", open(file, 'a'))
            pprint.pprint(fi_cal_in_error_step, open(file, 'a'))
            pprint.pprint('===========', open(file, 'a'))
        
        print("finish printing to file")

class check_ft_1_contact(checkforce):
    
    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2, error_tolerence, method_list)

    def f_name(self):
        pass

    def number_contact(self):
        pass

    def check_f(self, method_i):
        pass
    

class check_ft_1j_contact(check_ft_1_contact):
    
    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2, error_tolerence, method_list)
        print("check contact with single atom j")

    def f_name(self):
        return 'check_1contact_force_' + 'j' + '_' + self.string_idi_step()

    def number_contact(self):
        [n, id_collection] = cs.number_contact_atom_id_collection(f_custom_path, self.id_i, self.step1, self.step2)
        return n

    def check_f(self, method_i):
        result = cs.fjwi_plus_check_multistep_1contact_fromcustom(f_custom_path, self.id_i, self.step1, self.step2, self.error_tolerence, method_i, 'j')
        return result

    def checkprint(self):
        file = self.filepath()
        if os.path.isfile(file):
            print ('file exist, calculate again and append to the end')
        else:
            print ("create file and print to file")
        
        if self.number_contact()==1:

            self.print_string_list([
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.string_errortolerence(), self.stringstep(),
                self.string_wallname_id()
            ])
            for method_i in self.method_list:
                [f_step_error_array, fji_cal_in_error_step, fji_plus_in_error_step, id_j] = self.check_f(method_i)
                
                pprint.pprint("id_j {id_j}".format(id_j=id_j), open(file, 'a'))
                self.print_string_list([
                    "method {method_i}".format(method_i=method_i),
                    "f_step_error_array",
                ])
                pprint.pprint(f_step_error_array, open(file, 'a'))
                pprint.pprint("fji_plus_in_error_step", open(file, 'a'))
                pprint.pprint(fji_plus_in_error_step, open(file, 'a'))
                pprint.pprint("fji_cal_in_error_step", open(file, 'a'))
                pprint.pprint(fji_cal_in_error_step, open(file, 'a'))
                pprint.pprint('=========', open(file, 'a'))
        else:
            print("contact is {n_contact} not 1 so not print".format(n_contact=self.number_contact()))
            pprint.pprint("contact is {n_contact} not 1 so not print".format(n_contact=self.number_contact()), open(file, 'a'))
            pprint.pprint('=========', open(file, 'a'))
        print("finish printing to file")

class check_ft_1w_contact(check_ft_1_contact):

    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2, error_tolerence, method_list)
        print("check contact with single wall")

    def f_name(self):
        return 'check_1contact_force_' + 'w' + '_' + self.string_idi_step()

    def number_contact(self):
        [n, id_collection] = cs.number_contact_wall_id_collection(f_custom_path, self.id_i, self.step1, self.step2)
        return n

    def check_f(self, method_i):
        result = cs.fjwi_plus_check_multistep_1contact_fromcustom(f_custom_path, self.id_i, self.step1, self.step2, self.error_tolerence, method_i, 'w')
        return result

    def checkprint(self):
        file = self.filepath()
        if os.path.isfile(file):
            print ('file exist, calculate again and append to the end')
        else:
            print ("create file and print to file")
        
        if self.number_contact()==1:

            self.print_string_list([
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.string_errortolerence(), self.stringstep(),
                self.string_wallname_id()
            ])
            for method_i in self.method_list:
                [f_step_error_array, fji_cal_in_error_step, fji_plus_in_error_step, id_j] = self.check_f(method_i)
                
                pprint.pprint("id_j {id_j}".format(id_j=id_j), open(file, 'a'))
                self.print_string_list([
                    "method {method_i}".format(method_i=method_i),
                    "f_step_error_array",
                ])
                pprint.pprint(f_step_error_array, open(file, 'a'))
                pprint.pprint("fji_plus_in_error_step", open(file, 'a'))
                pprint.pprint(fji_plus_in_error_step, open(file, 'a'))
                pprint.pprint("fji_cal_in_error_step", open(file, 'a'))
                pprint.pprint(fji_cal_in_error_step, open(file, 'a'))
                pprint.pprint('=========', open(file, 'a'))
        else:
            print("contact is {n_contact} not 1 so not print".format(n_contact=self.number_contact()))
            pprint.pprint("contact is {n_contact} not 1 so not print".format(n_contact=self.number_contact()), open(file, 'a'))
            pprint.pprint('=========', open(file, 'a'))

        print("finish printing to file")

        
class checkoverlap(checkfile_look_idi):
    
    def __init__(self, id_i, step1, step2):
        super().__init__(id_i, step1, step2)
        print ('creating overlap')
    
    def f_name(self):
        return 'check_overlap_' + self.string_idi_step()
    
    def contact_check_multistep(self):
        result_contact_check_multistep = cs.contact_check_multistep(f_custom_path, self.id_i, self.step1, self.step2)
        return result_contact_check_multistep
    
    def checkprint(self):
        file = self.filepath()
        if os.path.isfile(file):
            print ('file exist, calculate again and append to the end')
        else:
            print ("create file and print to file")
        [step_id_ifover_diffnext, ifoverlap_ij_iw_array, initial_overlap_id] = self.contact_check_multistep()
        self.print_string_list([
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.stringstep(),
            self.string_wallname_id(),
            "initial overlap id are {first}".format(first=initial_overlap_id),
            "print step, id, and ifcontact. that ifcontact changed from last step to current step",
        ])
        pprint.pprint(step_id_ifover_diffnext, open(file, 'a'))
        pprint.pprint('==========', open(file, 'a'))
        print ("file printed")


# v1 means contact_check_multistep_v1
class checkoverlap_v1(checkoverlap):
    
    def __init__(self, id_i, step1, step2):
        super().__init__(id_i, step1, step2)

    def contact_check_multistep(self):
        result_contact_check_multistep = cs.contact_check_multistep_v1(f_custom_path, self.id_i, self.step1, self.step2)
        return result_contact_check_multistep



class checkmaxid(checkfile):
    
    def __init__(self, maxlabel, step1, step2):
        super().__init__(step1, step2)
        self.maxlabel = maxlabel
        self.maxh5file = dp.put_maxlabel_on_file(maxlabel, dp.f_custom)
        print ('creating maxid')
    
    def f_name(self):
        return 'maxlabel_' + self.maxlabel + "_step_" + str(self.step1) + "_" + str(self.step2)
    
    def checkprint(self):
        file = self.filepath()
        if os.path.isfile(file):
            print ('file exist, calculate again and append to the end')
        else:
            print ("create file and print to file")
        id_step = cs.max_id_step(self.maxh5file, self.step1, self.step2)
        self.print_string_list([
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.stringstep(),
            "step, id for maximum",
        ])
        pprint.pprint(id_step, open(file, 'a'))
        pprint.pprint('==========', open(file, 'a'))
        print ("file printed")