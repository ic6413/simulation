import createdata.datatofile as cd
import calculate.checkoutputfile as cco
import plotfigure.plottofile as ppf
import getattribute.create_attributes as gc



class scripttofile(object):
    def __init__(self):
        pass

    def startstring(self):
        pass

    def tofile(self):
        pass

    def executive(self):
        print(self.startstring())
        self.tofile()
        print('finished' + self.startstring())


class script_attribute(scripttofile):
    def __init__(self):
        super().__init__()
    
    def startstring(self):
        return 'creating attribute'

    def tofile(self):
        gc.define_attribute_dict()
        

class script_thermo(scripttofile):
    def __init__(self):
        super().__init__()
    
    def startstring(self):
        return 'creating thermo h5'

    def tofile(self):
        cd.thermo_hdf5_csv()


class script_custom(scripttofile):
    def __init__(self):
        super().__init__()
    
    def startstring(self):
        return 'creating custom h5'

    def tofile(self):
        cd.dump_custom()

class script_custom_single(scripttofile):
    def __init__(self, id_i):
        self.id_i = id_i
    
    def startstring(self):
        return 'creating custom_single h5 for idi = {idi}'.format(idi=self.id_i)

    def tofile(self):
        cd.dump_custom_select([self.id_i])

class script_steps(scripttofile):
    def __init__(self, step1, step2):
        self.step1 = step1
        self.step2 = step2


class script_idi_steps(script_steps):
    def __init__(self, id_i, step1, step2):
        super().__init__(step1, step2)
        self.id_i = id_i
        

class script_checkoverlap(script_idi_steps):
    def __init__(self, id_i, step1, step2):
        super().__init__(id_i, step1, step2)

    def startstring(self):
        return 'creating overlap'

    def tofile(self):
        cco.checkoverlap(self.id_i, self.step1, self.step2).checkprint()

class script_checkforce(script_idi_steps):
    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2)
        self.error_tolerence = error_tolerence
        self.method_list = method_list

class script_checkforce_all(script_checkforce):
    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2, error_tolerence, method_list)

    def startstring(self):
        return 'creating check force'

    def tofile(self):
        cco.checkforce(self.id_i, self.step1, self.step2, self.error_tolerence, self.method_list).checkprint()

class script_checkforce_1contactatom(script_checkforce):
    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2, error_tolerence, method_list)

    def startstring(self):
        return 'creating check 1contactatom force'

    def tofile(self):
        cco.check_ft_1j_contact(self.id_i, self.step1, self.step2, self.error_tolerence, self.method_list).checkprint()

class script_checkforce_1contactwall(script_checkforce):
    def __init__(self, id_i, step1, step2, error_tolerence, method_list):
        super().__init__(id_i, step1, step2, error_tolerence, method_list)

    def startstring(self):
        return 'creating check 1contactwall force'

    def tofile(self):
        cco.check_ft_1w_contact(self.id_i, self.step1, self.step2, self.error_tolerence, self.method_list).checkprint()


class script_plotij(script_idi_steps):
    def __init__(self, id_i, step1, step2):
        super().__init__(id_i, step1, step2)
    
    def startstring(self):
        return 'plotting ij, idi = {idi}'.format(idi=self.id_i)

    def tofile(self):
        ppf.plotclass(self.step1, self.step2).plotij(self.id_i)


class script_plotthermo(script_steps):
    def __init__(self, step1, step2, variable_name_list):
        super().__init__(step1, step2)
        self.variable_name_list = variable_name_list
    
    def startstring(self):
        return 'plotting thermo'

    def tofile(self):
        ppf.plotclass(self.step1, self.step2).plotthermo(self.variable_name_list)


class script_plotsingle(script_idi_steps):
    def __init__(self, id_i, step1, step2, variable_name_list):
        super().__init__(id_i, step1, step2)
        self.variable_name_list = variable_name_list
    
    def startstring(self):
        return 'plotting id {id}'.format(id=self.id_i)

    def tofile(self):
        try:
            ppf.plotclass(self.step1, self.step2).plotsingle(self.id_i, self.variable_name_list)
        except FileNotFoundError:
            script_custom_single(self.id_i).executive()

class script_plot3D(script_idi_steps):
    def __init__(self, id_i, step1, step2):
        super().__init__(id_i, step1, step2)
    
    def startstring(self):
        return 'plotting 3D id {id}'.format(id=self.id_i)

    def tofile(self):
        try:
            ppf.plotclass(self.step1, self.step2).plot3Dtraj(self.id_i)
        except FileNotFoundError:
            script_custom_single(self.id_i).executive()