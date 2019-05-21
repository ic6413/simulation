#!/usr/bin/env python


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
        import getattribute.create_attributes
script_attribute().executive()
import script.scriptclasses as ss
ss.script_thermo().executive()
ss.script_custom().executive()
id_i = 1
step1 = 900001
step2 = 1377826
ss.script_checkoverlap(id_i, step1, step2).executive()
error_tolerence = 1e-7
method_list = list(range(0, 3))
#ss.script_checkforce_all(id_i, step1, step2, error_tolerence, method_list).executive()
ss.script_checkforce_1contactatom(id_i, step1, step2, error_tolerence, method_list).executive()
ss.script_checkforce_1contactwall(id_i, step1, step2, error_tolerence, method_list).executive()
ss.script_plotij(id_i, step1, step2,).executive()
ss.script_plotthermo(step1, step2, 'all').executive()


ss.script_custom_single(id_i).executive()
y_variables = [
        'fx',
        'fy',
        'fz',
        'x',
        'y',
        'z',
        'vx',
        'vy',
        'vz',
        'omegax',
        'omegay',
        'omegaz',
        'v_KEtr_atom',
        'c_KEr_atom',
        'c_KEt_atom'
    ]
y_variables = 'all'

ss.script_plotsingle(id_i, step1, step2, y_variables).executive()







def script_skip(attrfile=1, hdf5file=1, overlapfile=1, forcefile=1, contact1frocefile=1, plotfile=1):
    
    if attrfile == 1:       
        import getattribute.create_attributes
    
    if hdf5file == 1:  
        import createdata.datatofile as cd
        print ('start create thermo')
        cd.thermo_hdf5_csv()
        print ('finished create thermo')
        print ('start create custom')
        cd.dump_custom()
        print ('finished create custom')
    
    print ('start print to file')
    import calculate.checkoutputfile as cco
    id_i = 15556
    step1 = 51914
    step2 = 870000
    error_tolerence = 1e-7
    method_list = list(range(1, 3))  # method = 1, 2 are correct. method 0 almost correct. method = 3~7(8) wrong
    if overlapfile == 1:    
        print ('start calculate overlap')
        cco.checkoverlap(id_i, step1, step2).checkprint()
        print ('end calculate overlap')
    if forcefile == 1:
        print ('start calculate ft')
        cco.checkforce(id_i, step1, step2, error_tolerence, method_list).checkprint()
    if contact1frocefile==1:    
        print ('start calculate ft for 1 j contact')
        cco.check_ft_1j_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()
        print ('start calculate ft for 1 wall contact')
        cco.check_ft_1w_contact(id_i, step1, step2, error_tolerence, method_list).checkprint()

    if plotfile == 1:
        print ('start plot')
        import plotfigure.plottofile as ppf
        id_i = 15556
        step1 = 0
        step2 = 50000
        ppf.plotclass(step1, step2).plotij(id_i)
        ppf.plotclass(step1, step2).plotthermo('all')

#script_skip(
#    attrfile=0,
#    hdf5file=0,
#    overlapfile=0,
#    forcefile=0,
#    contact1frocefile=1,
#    plotfile=1
#)