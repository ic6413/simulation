# data_post_process_LAMMPS
# import
import os.path
import json
import re
import pprint
import time
from itertools import chain
from itertools import repeat
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D
# import module in simulation folder
import os
import datapath as dp
import read_setting as rr


# ==== inputvariable varables ====
# current module variable
overlap_tolerence = 0
# see reading value as zero
abs_error_tolerence = dp.abs_error_tolerence
# ==== inputvariable varables end ====
# timestep
ts = float(rr.logfile['ts']) 
# atom radius
dp0 = float(rr.logfile['dp']) 
density = flaot(rr.logfile["den"])
# intersection pointmu
r_in = float(rr.logfile["r_in"]) 
r_out = float(rr.logfile["r_out"]) 
# gravity
g = float(rr.logfile["g"])  
# parameter
mu = float(rr.logfile["mu"]) 
kn = float(rr.logfile["kn"])  
kt = float(rr.logfile["kt"])  
gamma_n = float(rr.logfile["gamma_n"]) 
gamma_t = float(rr.logfile["gamma_t"]) 
type_radius_array = int(rr.logfile["type_radius_array"]) 
zlo_wall = rr.logfile["zlo_wall"]
walls = rr.logfile["walls"]
# ====================================== end import attribute


# transform list to numpy array
g = np.asarray(g)


def reset_nan_to_zero(A):
    where_are_NaNs = np.isnan(A)
    A[where_are_NaNs] = 0
    return A


# define extract dataframe by step1 step2
def extract_dataframe(df, step1, step2):

    def test1(df, step1, step2):

        if 'Step' in df.columns.values.tolist():
            steps = df['Step'].values
        elif 'step' in df.columns.values.tolist():
            steps = df['step'].values
            
        else:
            sys.exit('no step in header')
        
        firststep = steps[0]
        laststep = steps[-1]
        if step1 < firststep:
            sys.exit("step1 should be larger or equal to first step in array, step1={step1}, firststep={firststep}".format(step1=step1,firststep=firststep))
        if step2 > laststep+1:
            sys.exit("step2 should be smaller or equal to laststep+1 in array, step2={step2}, laststep={laststep}".format(step2=step2,laststep=laststep))

        row1 = steps.searchsorted(step1-0.5)
        row2 = steps.searchsorted(step2-0.5)
        df_step = df.iloc[row1:row2]
        return df_step

    #def test2(df, step1, step2):
    #    if 'Step' in df.columns.values.tolist():
    #        df_step = df.loc[df['Step'].isin(list(range(step1, step2)))]
    #    elif 'step' in df.columns.values.tolist():
    #        df_step = df.loc[df['step'].isin(list(range(step1, step2)))]
    #    else:
    #        sys.exit('no step in header')
    #    return df_step
    
    return test1(df, step1, step2)


# define function for get radius of each atom
def radius_by_type(type): 
    type_expand = type[..., np.newaxis]
    r = np.sum((type_expand == type_radius_array[0])*type_radius_array[1], axis=-1)
    return r

def mass_by_type(type, density):
    return 4/3*np.pi*radius_by_type(type)**3*density


# np.square np.inner cross
def length(vector):
	return np.sqrt(np.sum(vector**2, axis=-1, keepdims=True))


def unit(vector):
    unit = np.zeros_like(vector)
    index_vector_not_zero_not_nan = np.logical_and(np.any(vector != 0, axis=-1), ~np.any(np.isnan(vector), axis=-1))
    vector_select = vector[index_vector_not_zero_not_nan]
    unit_select = vector_select/length(vector_select)
    unit[index_vector_not_zero_not_nan] = unit_select
    return unit


def projection_scalar(v, n):
	return np.sum(v*unit(n), axis=-1, keepdims=True)


def projection_vector(v, n):
	return unit(n)*projection_scalar(v, n)


def verticle_vector(v, n):
	return v - projection_vector(v, n)


def costheta(v, n):
    return projection_scalar(v, n)/length(v)


# r_p_wall_vector, wall_position is radius or plane height
def riw_p_vector(xi, plane_point, plane_normal):
    return projection_vector(xi - plane_point, plane_normal)

def riw_cy_vector(xi, r, center_point, axis_vector):
    i_axis = verticle_vector(xi - center_point, axis_vector)
    return (length(i_axis) - r)*unit(i_axis)

class single_atom(object):

    def __init__(self, type, x, v, f, om, tq):
        self.type = type
        self.x = x
        self.f = f
        self.v = v
        self.om = om
        self.tq = tq

    def radius(self):
        return radius_by_type(self.type)

    def mass(self):
        return 4/3*np.pi*self.radius()**3*density

    def acceleration(self):
        return self.f/self.mass()

    def v_half(self):
        return self.v + self.acceleration()*ts/2

    def rotational_inertia(self):
        return 2/5*self.mass()*self.radius()**2

    def om_half(self):
        return self.om + self.tq/self.rotational_inertia()*ts/2


class overlap_to_i(object):
    
    def __init__(self, typei, xi):
        self.typei = typei
        self.xi = xi

    def xij(self):
        pass
    
    def overlap_length(self):
        pass

    def ifoverlap(self):
        return (self.overlap_length() >= overlap_tolerence)

    def overlapij_vector(self):
        return -self.overlap_length()*unit(self.xij())*self.ifoverlap()

class j_overlap_class(overlap_to_i):

    def __init__(
        self,
        typei, xi,
        typej, xj
    ):
        super().__init__(self, typei, xi)
        self.typej = typej
        self.xj = xj

    def xij(self):
        xij = self.xi - self.xj
        return xij
        
    def overlap_length(self):
        return radius_by_type(self.typei) + radius_by_type(self.typej) - length(self.xij())

class wall_overlap_class(overlap_to_i):

    def __init__(self,
                typei, xi,
                typew, xw,
                ):
        super().__init__(typei, xi)
        self.typew = typew
        self.xw = np.asarray(xw)

    def nearest_point_fun(self, xi):
	    pass

    def nearest_point(self):
	    return self.nearest_point_fun(self.xi)

    def xij(self):
        return self.xi - self.xw

    def overlap_length(self):
        pass

    def ifoverlap(self):
        return (self.overlap_length() >= overlap_tolerence)

    def overlapij_vector(self):
        pass

    def overlapiw(self):
        result1 = overlapij(radius_by_type(self.typei), 0, self.xi, self.nearest_point())
        ifo = result1[0]
        o_vector = result1[1]
        return [ifo, o_vector]

class wall_p_overlap_class(wall_overlap_class):

    def __init__(self,
                typei, xi,
                center_point, plane_normal,
                typew, xw,
                ):
        super().__init__(
            typei, xi,
            typew, xw,
            )
        self.center_point = np.asarray(center_point)
        self.plane_normal = np.asarray(plane_normal)
        self.r = 0

    def nearest_point_fun(self, xi):
	    return verticle_vector(xi - self.center_point, self.plane_normal) + self.center_point

    def overlap_length(self):
        d_to_plane = np.abs(projection_scalar(self.xi - self.center_point, self.plane_normal))
        return radius_by_type(self.typei) - d_to_plane

    def overlapij_vector(self):
        ps = projection_scalar(self.xi - self.center_point, self.plane_normal)
        return -self.overlap_length()*unit(self.plane_normal)*self.ifoverlap()*ps/np.abs(ps)

class wall_cy_overlap_class(wall_overlap_class):
    
    def __init__(self,
                typei, xi,
                center_point, axis_vector, r,
                typew, xw,
                ):
        super().__init__(
            typei, xi,
            typew, xw,
            )
        self.center_point = np.asarray(center_point)
        self.axis_vector = np.asarray(axis_vector)
        self.r = r

    def nearest_point_fun(self, xi):
        axispoint = projection_vector(xi - self.center_point, self.axis_vector) + self.center_point
        return self.r*unit(xi - axispoint) + axispoint

    def vector_to_axis_fun(self, xi):
        return verticle_vector(xi - self.center_point, self.axis_vector)
    
    def vector_to_axis(self):
        return self.vector_to_axis_fun(self.xi)

    def d_axis(self):
        return length(self.vector_to_axis())

    def d_to_surface(self):
        return self.d_axis() - self.r

    def overlap_length(self):
        return radius_by_type(self.typei) - np.abs(self.d_to_surface())

    def overlapij_vector(self):
        n = unit(self.vector_to_axis())*self.d_to_surface()/np.abs(self.d_to_surface())
        return -self.overlap_length()*n*self.ifoverlap()


class intersection_to_i(object):

    def __init__(self, typei, xi, vi, fi, omi, tqi, history_t_k, total_length, total_displacement, method):
        self.typei = typei
        self.xi = xi
        self.vi = vi
        self.fi = fi
        self.omi = omi
        self.tqi = tqi
        self.atomi = single_atom(self.typei, self.xi, self.vi, self.fi, self.omi, self.tqi)
        self.method = method
        self.history_t_k = history_t_k
        self.total_length = total_length
        self.total_displacement = total_displacement

    def vi_half(self):
        vi_half_xi_plus = update_v_x(self.xi, self.vi, self.fi/self.atomi.mass())
        vi_half = vi_half_xi_plus[0]
        return vi_half

    def xi_plus(self):
        vi_half_xi_plus = update_v_x(self.xi, self.vi, self.fi/self.atomi.mass())
        xi_plus = vi_half_xi_plus[1]
        return xi_plus

    def xij(self):
        pass

    def overlap_length(self):
        pass

    def ifoverlap(self):
        return (self.overlap_length() >= overlap_tolerence)

    def overlapij_vector(self):
        return -self.overlap_length()*unit(self.xij())*self.ifoverlap()
    
    def vij_half(self):
        pass

    def xij_plus(self):
        return self.xij()+self.vij_half()*ts

    def overlap_length_plus(self):
        pass

    def ifoverlap_plus(self):
        return (self.overlap_length_plus() >= overlap_tolerence)

    def overlapij_vector_plus(self):
        return -self.overlap_length_plus()*unit(self.xij_plus())*self.ifoverlap_plus()

    def vijn_half(self):
        return projection_vector(self.vij_half(), self.xij_plus())

    def vijt_half(self):
        pass
        
    def meff(self):
        pass

    def fji_n_k_plus(self):
        return -kn*self.overlapij_vector_plus()*self.ifoverlap_plus()

    def fji_n_gamma_plus(self):
        return -self.meff()*gamma_n*self.vijn_half()*self.ifoverlap_plus()
    
    def fji_n_plus(self):
        return self.fji_n_k_plus() + self.fji_n_gamma_plus()

    def fji_t_k_plus_one_step(self):
        return -kt*(self.vijt_half()*ts)*self.ifoverlap_plus()
    
    def fji_t_gamma_plus(self):
        return -self.meff()*gamma_t*self.vijt_half()*self.ifoverlap_plus()

    def displacement_ji_singlestep(self):
        return self.vijt_half()*ts*self.ifoverlap_plus()

    def fji_t_plus(self):
        
        method = self.method

        history_force_out = history_force(method, self.xij(), self.xij_plus(), self.vijt_half(), self.history_t_k, self.total_displacement, self.total_length)

        length_singlestep = length(self.displacement_ji_singlestep())
        fji_t_plus = self.fji_t_gamma_plus() + self.fji_t_k_plus_one_step() + history_force_out
        
        truncated = length(fji_t_plus) > mu*length(self.fji_n_plus())

        fji_t_plus = truncated*mu*length(self.fji_n_plus())*unit(fji_t_plus) + (1 - truncated)*fji_t_plus

        return fji_t_plus
    
    def history_t_k_plus(self):
        return self.fji_t_plus() - self.fji_t_gamma_plus()
    
    def total_displacement_plus(self):
        return -1/kt*self.history_t_k_plus()

    def total_length_plus(self):
        return length(-1/kt*self.history_t_k_plus())

    def tqji_plus(self):
        return np.cross(self.atomi.radius()*unit(-self.xij_plus()), self.fji_t_plus())


class j_class(intersection_to_i):

    def __init__(self,
                typei, xi, vi, fi, omi, tqi,
                typej, xj, vj, fj, omj, tqj, 
                history_t_k, total_length, total_displacement, method
    ):
        super().__init__(typei, xi, vi, fi, omi, tqi, history_t_k, total_length, total_displacement, method)
        self.typej = typej
        self.xj = xj
        self.vj = vj
        self.fj = fj
        self.omj = omj
        self.tqj = tqj
        self.atomj = single_atom(self.typej, self.xj, self.vj, self.fj, self.omj, self.tqj)

    def xij(self):
        xij = self.xi - self.xj
        return xij
        
    def overlap_length(self):
        return self.atomi.radius() + self.atomj.radius() - length(self.xij())

    def overlap_length_plus(self):
        return self.atomi.radius() + self.atomj.radius() - length(self.xij_plus())

    def vij_half(self):
        return self.atomi.v_half() - self.atomj.v_half()
    
    def vijt_half(self):
        vijt_half = (
            verticle_vector(self.vij_half(), self.xij_plus())
            +self.atomi.radius()*np.cross(self.atomi.om_half(), unit(-self.xij_plus()))
            -self.atomj.radius()*np.cross(self.atomj.om_half(), unit(self.xij_plus()))
        )
        return vijt_half
    
    def meff(self):
        mi = self.atomi.mass()
        mj = self.atomj.mass()
        return mi*mj/(mi + mj)


class wall_class(intersection_to_i):

    def __init__(self,
                typei, xi, vi, fi, omi, tqi,
                history_t_k, total_length, total_displacement, method,
                typew, xw, vw, aw, omw, alphaw,
                ):
        super().__init__(typei, xi, vi, fi, omi, tqi, history_t_k, total_length, total_displacement, method)
        self.typew = typew
        self.xw = np.asarray(xw)
        self.vw = np.asarray(vw)
        self.aw = np.asarray(aw)
        self.omw = np.asarray(omw)
        self.alphaw = np.asarray(alphaw)

    def nearest_point_fun(self, xi):
	    pass

    def nearest_point(self):
	    return self.nearest_point_fun(self.xi)
    
    def nearest_point_plus(self):
        return self.nearest_point_fun(self.xi_plus())

    def xij(self):
        return self.xi - self.xw

    def overlap_length(self):
        pass

    def overlap_length_plus(self):
        pass

    def vij_half(self):
        return self.atomi.v_half()

    def vijt_half(self):
        vijt_half = (
            verticle_vector(self.vij_half(), self.xij_plus())
            +self.atomi.radius()*np.cross(self.atomi.om_half(), unit(-self.xij_plus()))
        )
        return vijt_half

    def meff(self):
        mi = self.atomi.mass()
        return mi

    def ifoverlap(self):
        return (self.overlap_length() >= overlap_tolerence)

    def ifoverlap_plus(self):
        return (self.overlap_length_plus() >= overlap_tolerence)

    def overlapij_vector(self):
        pass

    def overlapij_vector_plus(self):
        pass

    def overlapiw(self):
        result2 = overlapij(self.atomi.radius(), 0, self.xi, self.nearest_point())
        ifo = result2[0]
        o_vector = result2[1]
        return [ifo, o_vector]

    def overlapiw_plus(self):
        result3 = overlapij(self.atomi.radius(), 0, self.xi_plus(), self.nearest_point_plus())
        ifo = result3[0]
        o_vector = result3[1]
        return [ifo, o_vector]


class wall_p_class(wall_class):

    def __init__(self,
                typei, xi, vi, fi, omi, tqi,
                center_point, plane_normal, history_t_k, total_length, total_displacement, method,
                typew, xw, vw, aw, omw, alphaw,
                ):
        super().__init__(
            typei, xi, vi, fi, omi, tqi, history_t_k, total_length, total_displacement, method,
            typew, xw, vw, aw, omw, alphaw,
            )
        self.center_point = np.asarray(center_point)
        self.plane_normal = np.asarray(plane_normal)
        self.r = 0

    def nearest_point_fun(self, xi):
	    return verticle_vector(xi - self.center_point, self.plane_normal) + self.center_point

    def overlap_length(self):
        d_to_plane = np.abs(projection_scalar(self.xi - self.center_point, self.plane_normal))
        return self.atomi.radius() - d_to_plane

    def overlap_length_plus(self):
        d_to_plane = np.abs(projection_scalar(self.xi_plus() - self.center_point, self.plane_normal))
        return self.atomi.radius() - d_to_plane 

    def overlapij_vector(self):
        ps = projection_scalar(self.xi - self.center_point, self.plane_normal)
        return -self.overlap_length()*unit(self.plane_normal)*self.ifoverlap()*ps/np.abs(ps)

    def overlapij_vector_plus(self):
        ps = projection_scalar(self.xi_plus() - self.center_point, self.plane_normal)
        return -self.overlap_length_plus()*unit(self.plane_normal)*self.ifoverlap_plus()*ps/np.abs(ps)


class wall_cy_class(wall_class):
    
    def __init__(self,
                typei, xi, vi, fi, omi, tqi,
                center_point, axis_vector, r, history_t_k, total_length, total_displacement, method,
                typew, xw, vw, aw, omw, alphaw,
                ):
        super().__init__(
            typei, xi, vi, fi, omi, tqi, history_t_k, total_length, total_displacement, method,
            typew, xw, vw, aw, omw, alphaw,
            )
        self.center_point = np.asarray(center_point)
        self.axis_vector = np.asarray(axis_vector)
        self.r = r

    def nearest_point_fun(self, xi):
        axispoint = projection_vector(xi - self.center_point, self.axis_vector) + self.center_point
        return self.r*unit(xi - axispoint) + axispoint


    def vector_to_axis_fun(self, xi):
        return verticle_vector(xi - self.center_point, self.axis_vector)
    
    def vector_to_axis(self):
        return self.vector_to_axis_fun(self.xi)

    def vector_to_axis_plus(self):
        return self.vector_to_axis_fun(self.xi_plus())

    def d_axis(self):
        return length(self.vector_to_axis())

    def d_axis_plus(self):
        return length(self.vector_to_axis_plus())

    def d_to_surface(self):
        return self.d_axis() - self.r

    def d_to_surface_plus(self):
        return self.d_axis_plus() - self.r

    def overlap_length(self):
        return self.atomi.radius() - np.abs(self.d_to_surface())

    def overlap_length_plus(self):
        return self.atomi.radius() - np.abs(self.d_to_surface_plus())

    def overlapij_vector(self):
        n = unit(self.vector_to_axis())*self.d_to_surface()/np.abs(self.d_to_surface())
        return -self.overlap_length()*n*self.ifoverlap()

    def overlapij_vector_plus(self):
        n = unit(self.vector_to_axis_plus())*self.d_to_surface_plus()/np.abs(self.d_to_surface_plus())
        return -self.overlap_length_plus()*n*self.ifoverlap_plus()

def setnocontacttonan(array, ifoverlapcheck):
    array[~ifoverlapcheck[:,0],:]=np.nan
    return array

def overlapij(ri, rj, xi, xj):
    xij = (xi - xj)
    overlap_length = ri + rj - length(xij)
    ifoverlap = (overlap_length >= overlap_tolerence)
    overlapij_vector = -overlap_length*unit(xij)*ifoverlap
    overlap_length = setnocontacttonan(overlap_length, ifoverlap)
    overlapij_vector = setnocontacttonan(overlapij_vector, ifoverlap)
    return [ifoverlap, overlapij_vector, overlap_length]


def contact_ids_inonestep(df_onestep, id_i):
    
    dfi_onestep = df_onestep.loc[df_onestep['id']==id_i]
    
    if len(dfi_onestep.index) != 1:
        sys.exit('row number of dfi_onestep is not 1')
    else:
        pass

    typei = dfi_onestep[['type']].values
    ri = radius_by_type(typei)
    xi = dfi_onestep[['x','y','z']].values

    dfjs_onestep = df_onestep.loc[df_onestep['id']!=id_i]
    typej = dfjs_onestep[['type']].values
    rj = radius_by_type(typej)
    xj = dfjs_onestep[['x','y','z']].values
    [ifoverlap_ij_onestep, overlapij_vector_onestep, overlapij_length_onestep] = overlapij(ri, rj, xi, xj)

    contact_idj = dfjs_onestep[['id']].values[ifoverlap_ij_onestep]

    idwalls = np.empty(len(walls), dtype=int)
    ifoverlap_iwall_onestep = np.empty(len(walls), dtype=bool)
    
    for n, walllist in enumerate(walls):

        idwalls[n] = -n-1

        if walllist[0] == 'p':
            wall = wall_p_overlap_class(typei, xi, walllist[1], walllist[2], walllist[3], walllist[4])
        elif walllist[0] == 'cy':
            wall = wall_cy_overlap_class(typei, xi, walllist[1], walllist[2], walllist[3], walllist[4], walllist[5])
        else:
            print('walltype not p not cy')

        ifoverlap_iwall_onestep[n] = wall.ifoverlap()
    contact_idwall = idwalls[ifoverlap_iwall_onestep]

    contact_idj_wall = np.append(contact_idj,contact_idwall)

    return contact_idj_wall

def airviscous(v):
    viscous = rr.logfile["gamma_air"]*rr.logfile["mp"]
    return -viscous*v

def gravity(mi):
    return mi*g

def force_sum_except_contact(mi, v):
    force_sum_except_contact = 0
    force_sum_except_contact += gravity(mi)
    if rr.logfile["ifairviscous"]=="yes":
        force_sum_except_contact += airviscous(v)
    return force_sum_except_contact

def history_force(method, xij_cal, xij_plus_cal, vijt_half, previous_ktforce, total_displacement, total_length):
    
    ps = projection_scalar(total_displacement, vijt_half)
    project_scalar_unit = np.divide(ps, abs(ps), out=np.zeros_like(ps), where=ps!=0) 

    if method == 0:
        # covariant  =============================================
        n1 = unit(xij_cal)
        n2 = unit(xij_plus_cal)
        dn = n2 - n1
        t1 = unit(verticle_vector(n2, n1))
        t2 = unit(verticle_vector(-n1, n2))
        history_force_out = previous_ktforce + projection_scalar(previous_ktforce, t1)*(t2-t1)

    elif method == 1:
        # covariant  =============================================
        n1 = unit(xij_cal)
        n2 = unit(xij_plus_cal)
        dn = n2 - n1
        dn_t1 = verticle_vector(dn, n1)
        dn_t2 = verticle_vector(dn, n2)
        history_ktforce_para = projection_scalar(previous_ktforce, dn_t1)*unit(dn_t2)
        history_ktforce_verti = verticle_vector(previous_ktforce, dn_t1)
        history_force_out = history_ktforce_para + history_ktforce_verti

    elif method == 2:
        # correct================================================================================
        history_force_out = -kt*verticle_vector(total_displacement, xij_plus_cal)

    elif method == 3:
        # maybe correct================================================================================
        history_force_out = -kt*length(total_displacement)*unit(vijt_half)*project_scalar_unit

    elif method == 4:
        # maybe correct================================================================================
        history_force_out = -kt*projection_scalar(total_displacement, vijt_half)*unit(vijt_half)

    # wrong when direction change#===============
    elif method == 5:
        history_force_out = -kt*total_length*unit(vijt_half)
        
    elif method == 6:
        history_force_out = -kt*total_length*unit(vijt_half)*project_scalar_unit

    elif method == 7:
        history_force_out = -kt*length(total_displacement)*unit(vijt_half)

    elif method == 8:
        history_force_out = -kt*total_length*unit(vijt_half)*project_scalar_unit

    else:
        sys.exit('method not defined')

    return history_force_out    


def error_ratio(error, denominator):
    return np.divide(error, denominator, out=np.zeros_like(error), where=(denominator!=0))


def add_step_to_array(steps_array, array):
    
    n_row = array.shape[0]
    n_column = array.shape[1]+1
    array_addstep = np.empty([n_row, n_column])
    array_addstep[:,0:1] = steps_array
    array_addstep[:,1:] = array
    
    return array_addstep


def step_error_array_with_index(step1, step2, fi_cal, fi, error_tolerence):

    #f_diff = (fi_cal - fi)
    #fji_error_steps = np.divide(f_diff, fi, out=np.zeros_like(f_diff), where=fi!=0)
    
    #errornorm = np.nansum(abs(fji_error_steps), axis=-1)
    #errorindex = (errornorm > error_tolerence)

    errorindex = np.logical_not(
        np.isclose(fi_cal, fi, rtol=error_tolerence, atol=abs_error_tolerence).all(axis=-1)
    )

    fji_error_steps = np.absolute((fi_cal[errorindex,:]-fi[errorindex,:])/fi[errorindex,:])
    steps_array = np.arange(step1+1, step2)[errorindex].reshape((-1,1))
    step_error_array = add_step_to_array(steps_array, fji_error_steps)
    return [step_error_array, errorindex]


def update_v_x(x,v,a):
    v_half = v + a*ts/2
    x_plus = x + (v + a*ts/2)*ts
    return [v_half, x_plus]


def update_v_x_pre(x,v,a):
    v_half_pre = v + a*(-ts)/2
    x_minus = x + (v + a*(-ts)/2)*(-ts)
    return [v_half_pre, x_minus]


def calculate_m_r_I_vh_xp_omh(type, density, f, tq, x, v, om):
    
    m = mass_by_type(type, density)
    r = radius_by_type(type)
    I = 2/5*m*r**2
    a = f/m
    alpha = tq/I
    [v_half, x_plus_cal]=update_v_x(x,v,a)
    result4=update_v_x(0,om,alpha)
    om_half=result4[0]
    return [m, r, I, v_half, x_plus_cal, om_half]


def calculate_m_r_I_vh_xp_omh_pre(type, density, f, tq, x, v, om):
    
    m = mass_by_type(type, density)
    r = radius_by_type(type)
    I = 2/5*m*r**2
    a = f/m
    alpha = tq/I
    [v_half_pre, x_minus_cal]=update_v_x_pre(x,v,a)
    result5=update_v_x_pre(0,om,alpha)
    om_half_pre = result5[0]
    return [m, r, I, v_half_pre, x_minus_cal, om_half_pre]


def calculate_wall_r_vh_xp_omh(wall):
    
    [v_half, x_plus_cal]=update_v_x(wall.xw,wall.vw,wall.aw)
    result6=update_v_x(np.zeros(3),wall.omw,wall.alphaw)
    om_half = result6[0]
    r = wall.r
    
    return [r, v_half, x_plus_cal, om_half]


def calculate_wall_r_vh_xp_omh_pre(wall):
    
    [v_half, x_plus_cal]=update_v_x_pre(wall.xw,wall.vw,wall.aw)
    [om_half, d_theta]=update_v_x_pre(0,wall.omw,wall.alphaw)
    r = wall.r 
    
    return [r, v_half, x_plus_cal, om_half]



def create_wall_class_from_walllist(
    walllist,type,x,v,f,om,tq, ftk_include_his, total_length, total_displacement, method,
):
    
    if walllist[0] == 'p':
        wall = wall_p_class(
                            type,x,v,f,om,tq,
                            walllist[1], walllist[2], ftk_include_his, total_length, total_displacement, method,
                            walllist[3], walllist[4], walllist[5], walllist[6], walllist[7], walllist[8],
                            )
    elif walllist[0] == 'cy':
        wall = wall_cy_class(
                            type,x,v,f,om,tq,
                            walllist[1], walllist[2], walllist[3], ftk_include_his, total_length, total_displacement,method,
                            walllist[4], walllist[5], walllist[6], walllist[7], walllist[8], walllist[9]
                            )
    else:
        print('walltype not p not cy')

    return wall


def reindex_by_step(df, step1, step2):
    df = df.set_index('step')
    df = df.reindex(np.arange(step1, step2))
    return df


def get_type_x_v_f_om_tq_from_df(dfi, step1, step2):
    
    type = dfi[['type']].values
    x = dfi[['x','y','z']].values
    v = dfi[['vx','vy','vz']].values
    f = dfi[['fx','fy','fz']].values
    om = dfi[['omegax','omegay','omegaz']].values
    tq = dfi[['tqx','tqy','tqz']].values

    return [type, x, v, f, om, tq]


def cal_f_no_history(overlapij_vector_plus_cal, ifoverlap_plus_cal, meff, vijn_half, vijt_half):

    fnk = -kn*overlapij_vector_plus_cal*ifoverlap_plus_cal
    fnk = reset_nan_to_zero(fnk)
    fngamma = -meff*gamma_n*vijn_half*ifoverlap_plus_cal
    fngamma = reset_nan_to_zero(fngamma)
    ftk = -kt*(vijt_half*ts)*ifoverlap_plus_cal
    ftk = reset_nan_to_zero(ftk)
    ftgamma = -meff*gamma_t*vijt_half*ifoverlap_plus_cal
    ftgamma = reset_nan_to_zero(ftgamma)

    return [fnk,fngamma,ftk,ftgamma]


def fji_plus_cal(typei, typej, xi, xj, vi, vj, fi, fj, previous_ktforce, total_length, total_displacement, omi, omj, tqi, tqj, method):
    
    [mi, ri, I_i, vi_half, xi_plus_cal, omi_half] = calculate_m_r_I_vh_xp_omh(typei, density, fi, tqi, xi, vi, omi)
    [mj, rj, I_j, vj_half, xj_plus_cal, omj_half] = calculate_m_r_I_vh_xp_omh(typej, density, fj, tqj, xj, vj, omj)

    meff = mi*mj/(mi + mj)
    vij_half = vi_half - vj_half
    xij_cal = xi - xj
    xij_plus_cal = xi_plus_cal - xj_plus_cal
    [ifoverlap_plus_cal, overlapij_vector_plus_cal, overlapij_length_plus_cal] = overlapij(ri, rj, xi_plus_cal, xj_plus_cal)
    vijn_half = projection_vector(vij_half, xij_plus_cal)
    vijt_half_no_om = verticle_vector(vij_half, xij_plus_cal)
    vijt_half = vijt_half_no_om + ri*np.cross(omi_half, unit(-xij_plus_cal)) - rj*np.cross(omj_half, unit(xij_plus_cal))
    
    [fnk,fngamma,ftk,ftgamma] = cal_f_no_history(overlapij_vector_plus_cal, ifoverlap_plus_cal, meff, vijn_half, vijt_half)
    fjin_plus_cal = fnk + fngamma
    
    displacement_singlestep = vijt_half*ts*ifoverlap_plus_cal
    displacement_singlestep = reset_nan_to_zero(displacement_singlestep)

    history_force_out = history_force(method, xij_cal, xij_plus_cal, vijt_half, previous_ktforce, total_displacement, total_length)
    
    history_force_out = history_force_out*ifoverlap_plus_cal
    history_force_out = reset_nan_to_zero(history_force_out)
    total_displacement = reset_nan_to_zero(total_displacement)
    total_length = reset_nan_to_zero(total_length)

    length_singlestep = length(displacement_singlestep)
    total_displacement = total_displacement + displacement_singlestep
    total_length = total_length + length_singlestep

    fjit_plus_cal = ftgamma + ftk + history_force_out
    truncated = length(fjit_plus_cal) > mu*length(fjin_plus_cal)
    fjit_plus_cal = truncated*mu*length(fjin_plus_cal)*unit(fjit_plus_cal) + (1 - truncated)*fjit_plus_cal
    ftk_include_his_out = fjit_plus_cal - ftgamma
    total_displacement = -1/kt*ftk_include_his_out
    total_length = length(-1/kt*ftk_include_his_out)
    
    tqji_plus_cal = np.cross(ri*unit(-xij_plus_cal),fjit_plus_cal)

    return [ftk_include_his_out, total_displacement, fjit_plus_cal, fjin_plus_cal, tqji_plus_cal, total_length]


def fwi_plus_cal(typei, wall, xi, vi, fi, previous_ktforce, total_length, total_displacement, omi, tqi, method):
    
    [mi, ri, I_i, vi_half, xi_plus_cal, omi_half] = calculate_m_r_I_vh_xp_omh(typei, density, fi, tqi, xi, vi, omi)
    [rj, vj_half, xj_plus_cal, omj_half] = calculate_wall_r_vh_xp_omh(wall)

    meff = wall.meff()
    xij_cal = xi - wall.nearest_point()
    xij_plus_cal = xi_plus_cal - wall.nearest_point_plus()
    
    ifoverlap_plus_cal = wall.ifoverlap_plus()
    overlapij_vector_plus_cal = wall.overlapij_vector_plus()

    vij_half = vi_half - vj_half
    vijn_half = projection_vector(vij_half, xij_plus_cal)
    vijt_half_no_om = verticle_vector(vij_half, xij_plus_cal)
    vijt_half = vijt_half_no_om + ri*np.cross(omi_half, unit(-xij_plus_cal)) + rj*np.cross(omj_half, unit(xij_plus_cal))

    [fnk,fngamma,ftk,ftgamma] = cal_f_no_history(overlapij_vector_plus_cal, ifoverlap_plus_cal, meff, vijn_half, vijt_half)
    fjin_plus_cal = fnk + fngamma

    displacement_singlestep = vijt_half*ts*ifoverlap_plus_cal
    displacement_singlestep = reset_nan_to_zero(displacement_singlestep)

    history_force_out = history_force(method, xij_cal, xij_plus_cal, vijt_half, previous_ktforce, total_displacement, total_length)

    history_force_out = history_force_out*ifoverlap_plus_cal
    history_force_out = reset_nan_to_zero(history_force_out)
    total_displacement = reset_nan_to_zero(total_displacement)
    total_length = reset_nan_to_zero(total_length)
    
    length_singlestep = length(displacement_singlestep)
    total_displacement = total_displacement + displacement_singlestep
    total_length = total_length + length_singlestep

    fjit_plus_cal = ftgamma + ftk + history_force_out
    
    truncated = length(fjit_plus_cal) > mu*length(fjin_plus_cal)
    fjit_plus_cal = truncated*mu*length(fjin_plus_cal)*unit(fjit_plus_cal) + (1 - truncated)*fjit_plus_cal
    ftk_include_his_out = fjit_plus_cal - ftgamma
    total_displacement = -1/kt*ftk_include_his_out
    total_length = length(-1/kt*ftk_include_his_out)
    
    tqji_plus_cal = np.cross(ri*unit(-xij_plus_cal),fjit_plus_cal)
        
    return [ftk_include_his_out, total_displacement, fjit_plus_cal, fjin_plus_cal, tqji_plus_cal, total_length]


def fjwi_plus_cal_multistep_multicontact_fromcustom(f_read, id_i, step1, step2, method, version=0):

    if version==0:
        df = pd.read_hdf(f_read, 'df')
        
        df_step = extract_dataframe(df, step1, step2)

        df_step_id_i = df_step.loc[df_step['id']==id_i]
        n_step1tostep2_id_i = len(df_step_id_i.index)
        n_step = step2 - step1
        if n_step1tostep2_id_i != n_step:
            steps_expected = np.arange(step1, step2)
            # only compare for same length part before array end
            min_len = min(n_step1tostep2_id_i, n_step)
            checksteps = df_step_id_i['step'].values[0:min_len] == steps_expected[0:min_len]
            if checksteps.all():
                id_first_wrong_step = min_len
            else:
                id_first_wrong_step = np.nonzero(~checksteps)[0][0]

            first_step_wrong = steps_expected[id_first_wrong_step]
            sys.exit('id_i is not in df read on step = {first_step_wrong}'.format(first_step_wrong=first_step_wrong))
        else:
            pass

        df_firststep = extract_dataframe(df, step1, step1+1)
            
        fisststep_contact_ids = contact_ids_inonestep(df_firststep, id_i)

        if len(fisststep_contact_ids) != 0:
            sys.exit('there exist contact in the begin step. can not get history. contact id list = {fisststep_contact_ids}'.format(fisststep_contact_ids=fisststep_contact_ids))
        else:
            pass

        groups_byid = df_step.groupby(['id'])

        dfi = reindex_by_step(groups_byid.get_group(id_i), step1, step2)
        
        typei = dfi[['type']].values[:-1]
        xi = dfi[['x','y','z']].values[:-1]
        vi = dfi[['vx','vy','vz']].values[:-1]
        fi = dfi[['fx','fy','fz']].values[:-1]
        omi = dfi[['omegax','omegay','omegaz']].values[:-1]
        tqi = dfi[['tqx','tqy','tqz']].values[:-1]
        xi_plus = dfi[['x','y','z']].values[1:]
        fi_plus = dfi[['fx','fy','fz']].values[1:]
        tqi_plus = dfi[['tqx','tqy','tqz']].values[1:]

        sum_fjit_plus_cal_steps = 0
        sum_fjin_plus_cal_steps = 0
        sum_tqji_plus_cal_steps = 0

        id_list = groups_byid.groups.keys()

        id_j_list = [i for i in id_list if i !=id_i]

        for id_j in id_j_list:

            dfj = reindex_by_step(groups_byid.get_group(id_j), step1, step2)

            typej = dfj[['type']].values[:-1]  
            xj = dfj[['x','y','z']].values[:-1]
            vj = dfj[['vx','vy','vz']].values[:-1]
            fj = dfj[['fx','fy','fz']].values[:-1]
            omj = dfj[['omegax','omegay','omegaz']].values[:-1]
            tqj = dfj[['tqx','tqy','tqz']].values[:-1]
            

            fjit_plus_cal_steps = np.empty([step2-step1-1, 3])
            fjin_plus_cal_steps = np.empty([step2-step1-1, 3])
            tqji_plus_cal_steps = np.empty([step2-step1-1, 3])

            fjit_plus_cal = np.zeros((1,3))
            fjin_plus_cal = np.zeros((1,3))
            tqji_plus_cal = np.zeros((1,3))
            ftk_include_his = np.zeros((1,3))
            total_displacement = np.zeros((1,3))
            total_length = np.zeros((1,1))
            

            for step in range(step1+1, step2):

                k = step - (step1+1)

                [ftk_include_his, total_displacement, fjit_plus_cal, fjin_plus_cal, tqji_plus_cal, total_length] = (
                    fji_plus_cal(typei[k:k+1], typej[k:k+1], xi[k:k+1], xj[k:k+1], vi[k:k+1], vj[k:k+1], fi[k:k+1], fj[k:k+1], ftk_include_his, total_length, total_displacement, omi[k:k+1], omj[k:k+1], tqi[k:k+1], tqj[k:k+1], method=method)
                )
                fjit_plus_cal_steps[k:k+1] = fjit_plus_cal
                fjin_plus_cal_steps[k:k+1] = fjin_plus_cal
                tqji_plus_cal_steps[k:k+1] = tqji_plus_cal
                

            sum_fjit_plus_cal_steps = sum_fjit_plus_cal_steps + reset_nan_to_zero(fjit_plus_cal_steps)
            sum_fjin_plus_cal_steps = sum_fjin_plus_cal_steps + reset_nan_to_zero(fjin_plus_cal_steps)
            sum_tqji_plus_cal_steps = sum_tqji_plus_cal_steps + reset_nan_to_zero(tqji_plus_cal_steps)

        for walllist in walls:

            fjit_plus_cal_steps = np.empty([step2-step1-1, 3])
            fjin_plus_cal_steps = np.empty([step2-step1-1, 3])
            tqji_plus_cal_steps = np.empty([step2-step1-1, 3])

            fjit_plus_cal = np.zeros((1,3))
            fjin_plus_cal = np.zeros((1,3))
            tqji_plus_cal = np.zeros((1,3))
            ftk_include_his = np.zeros((1,3))
            total_displacement = np.zeros((1,3))
            total_length = np.zeros((1,1))
            total_project_length = np.zeros((1,1))

            xj = 0*xi
            vj = 0*xi
            omj = 0*xi
            xj_plus = 0*xi_plus
            
            for step in range(step1+1, step2):

                k = step - (step1+1)

                wall = create_wall_class_from_walllist(
                            walllist,
                            typei[k:k+1],xi[k:k+1],vi[k:k+1],fi[k:k+1],omi[k:k+1],tqi[k:k+1],
                            ftk_include_his, total_length, total_displacement, method,
                            )


                [ftk_include_his, total_displacement, fjit_plus_cal, fjin_plus_cal, tqji_plus_cal, total_length] = (
                    fwi_plus_cal(typei[k:k+1], wall, xi[k:k+1], vi[k:k+1], fi[k:k+1], ftk_include_his, total_length, total_displacement, omi[k:k+1], tqi[k:k+1], method=method)
                )
                fjit_plus_cal_steps[k:k+1] = fjit_plus_cal
                fjin_plus_cal_steps[k:k+1] = fjin_plus_cal
                tqji_plus_cal_steps[k:k+1] = tqji_plus_cal

            sum_fjit_plus_cal_steps = sum_fjit_plus_cal_steps + reset_nan_to_zero(fjit_plus_cal_steps)
            sum_fjin_plus_cal_steps = sum_fjin_plus_cal_steps + reset_nan_to_zero(fjin_plus_cal_steps)
            sum_tqji_plus_cal_steps = sum_tqji_plus_cal_steps + reset_nan_to_zero(tqji_plus_cal_steps)
        
        

    elif version==1:
        df = pd.read_hdf(f_read, 'df')
        df_step = extract_dataframe(df, step1, step2)

        df_step_id_i = df_step.loc[df_step['id']==id_i]
        n_step1tostep2_id_i = len(df_step_id_i.index)
        n_step = step2 - step1
        if n_step1tostep2_id_i != n_step:
            steps_expected = np.arange(step1, step2)
            # only compare for same length part before array end
            min_len = min(n_step1tostep2_id_i, n_step)
            checksteps = df_step_id_i['step'].values[0:min_len] == steps_expected[0:min_len]
            if checksteps.all():
                id_first_wrong_step = min_len
            else:
                id_first_wrong_step = np.nonzero(~checksteps)[0][0]

            first_step_wrong = steps_expected[id_first_wrong_step]
            sys.exit('id_i is not in df read on step = {first_step_wrong}'.format(first_step_wrong=first_step_wrong))
        else:
            pass
            
        df_firststep = extract_dataframe(df, step1, step1+1)

        fisststep_contact_ids = contact_ids_inonestep(df_firststep, id_i)

        if len(fisststep_contact_ids) != 0:
            sys.exit('there exist contact in the begin step. can not get history. contact id list = {fisststep_contact_ids}'.format(fisststep_contact_ids=fisststep_contact_ids))
        else:
            pass

        groups_byid = df_step.groupby(['id'])
        
        dfi = (groups_byid.get_group(id_i))
        
        typei = dfi[['type']].values[:-1]
        xi = dfi[['x','y','z']].values[:-1]
        vi = dfi[['vx','vy','vz']].values[:-1]
        fi = dfi[['fx','fy','fz']].values[:-1]
        omi = dfi[['omegax','omegay','omegaz']].values[:-1]
        tqi = dfi[['tqx','tqy','tqz']].values[:-1]
        fi_plus = dfi[['fx','fy','fz']].values[1:]
        tqi_plus = dfi[['tqx','tqy','tqz']].values[1:]

        sum_fjit_plus_cal_steps = 0
        sum_fjin_plus_cal_steps = 0
        sum_tqji_plus_cal_steps = 0

        id_list = groups_byid.groups.keys()

        id_j_list = [i for i in id_list if i !=id_i]

        for id_j in id_j_list:

            dfj = (groups_byid.get_group(id_j))

            index_j_exist = (dfj[['step']].values).astype(int)[:,0].T - step1
            fjit_plus_cal_steps = np.empty([step2-step1+1, 3])
            fjin_plus_cal_steps = np.empty([step2-step1+1, 3])
            tqji_plus_cal_steps = np.empty([step2-step1+1, 3])

            history_t_k_many_steps = np.zeros((step2-step1+1,3))
            total_displacement_many_steps = np.zeros((step2-step1+1,3))
            total_length_many_steps = np.zeros((step2-step1+1,1))

            for k, n in enumerate(index_j_exist):

                typei = dfi[['type']].values[n:n+1]
                xi = dfi[['x','y','z']].values[n:n+1]
                vi = dfi[['vx','vy','vz']].values[n:n+1]
                fi = dfi[['fx','fy','fz']].values[n:n+1]
                omi = dfi[['omegax','omegay','omegaz']].values[n:n+1]
                tqi = dfi[['tqx','tqy','tqz']].values[n:n+1]
                #xi_plus = dfi[['x','y','z']].values[1:]
                #fi_plus = dfi[['fx','fy','fz']].values[1:]
                #tqi_plus = dfi[['tqx','tqy','tqz']].values[1:]


                typej = dfj[['type']].values[k:k+1]
                xj = dfj[['x','y','z']].values[k:k+1]
                vj = dfj[['vx','vy','vz']].values[k:k+1]
                fj = dfj[['fx','fy','fz']].values[k:k+1]
                omj = dfj[['omegax','omegay','omegaz']].values[k:k+1]
                tqj = dfj[['tqx','tqy','tqz']].values[k:k+1]
                #xj_plus = dfj[['x','y','z']].values[n+1:n+2]
                #fj_plus = dfj[['fx','fy','fz']].values[n+1:n+2]
                #tqj_plus = dfj[['tqx','tqy','tqz']].values[n+1:n+2]



                jtoi = j_class(
                    typei, xi, vi, fi, omi, tqi,
                    typej, xj, vj, fj, omj, tqj,
                    history_t_k_many_steps[n], total_length_many_steps[n], total_displacement_many_steps[n],method
                    )
                
                fjit_plus_cal_steps[n:n+1] = jtoi.fji_t_plus()
                fjin_plus_cal_steps[n:n+1] = jtoi.fji_n_plus()
                tqji_plus_cal_steps[n:n+1] = jtoi.tqji_plus()
                history_t_k_many_steps[n:n+1] = jtoi.history_t_k_plus()
                total_displacement_many_steps[n:n+1] = jtoi.total_displacement_plus()
                total_length_many_steps[n:n+1] = jtoi.total_length_plus()

            sum_fjit_plus_cal_steps = sum_fjit_plus_cal_steps + fjit_plus_cal_steps
            sum_fjin_plus_cal_steps = sum_fjin_plus_cal_steps + fjin_plus_cal_steps
            sum_tqji_plus_cal_steps = sum_tqji_plus_cal_steps + tqji_plus_cal_steps



        for walllist in walls:

            fjit_plus_cal_steps = np.empty([step2-step1+1, 3])
            fjin_plus_cal_steps = np.empty([step2-step1+1, 3])
            tqji_plus_cal_steps = np.empty([step2-step1+1, 3])

            history_t_k_many_steps = np.zeros((step2-step1+1,3))
            total_displacement_many_steps = np.zeros((step2-step1+1,3))
            total_length_many_steps = np.zeros((step2-step1+1,1))

            for n in range(0, step2-step1):

                typei = dfi[['type']].values[n:n+1]
                xi = dfi[['x','y','z']].values[n:n+1]
                vi = dfi[['vx','vy','vz']].values[n:n+1]
                fi = dfi[['fx','fy','fz']].values[n:n+1]
                omi = dfi[['omegax','omegay','omegaz']].values[n:n+1]
                tqi = dfi[['tqx','tqy','tqz']].values[n:n+1]
                #xi_plus = dfi[['x','y','z']].values[1:]
                #fi_plus = dfi[['fx','fy','fz']].values[1:]
                #tqi_plus = dfi[['tqx','tqy','tqz']].values[1:]

                wtoi = create_wall_class_from_walllist(
                    walllist,typei, xi, vi, fi, omi, tqi,
                    history_t_k_many_steps[n], total_length_many_steps[n], total_displacement_many_steps[n],method,)



                fjit_plus_cal_steps[n:n+1] = wtoi.fji_t_plus()
                fjin_plus_cal_steps[n:n+1] = wtoi.fji_n_plus()
                tqji_plus_cal_steps[n:n+1] = wtoi.tqji_plus()
                history_t_k_many_steps[n:n+1] = wtoi.history_t_k_plus()
                total_displacement_many_steps[n:n+1] = wtoi.total_displacement_plus()
                total_length_many_steps[n:n+1] = wtoi.total_length_plus()

            sum_fjit_plus_cal_steps = sum_fjit_plus_cal_steps + (fjit_plus_cal_steps)
            sum_fjin_plus_cal_steps = sum_fjin_plus_cal_steps + (fjin_plus_cal_steps)
            sum_tqji_plus_cal_steps = sum_tqji_plus_cal_steps + (tqji_plus_cal_steps)
        
    # other force, gravity, airviscous
    typei = dfi[['type']].values[:-1]
    mi = mass_by_type(typei, density)
    fi_only_contact_cal_plus = sum_fjit_plus_cal_steps + sum_fjin_plus_cal_steps
    if int(rr.logfile['freq_ave_wall']) != 0:
        fi_only_contact_plus = dfi[['f_sum_pairforce[1]','f_sum_pairforce[2]','f_sum_pairforce[3]']].values[1:] 
    else:
        fi_only_contact_plus = fi_plus - force_sum_except_contact(mi, vi+0.5*fi/mi*ts)
    fi_cal_plus = fi_only_contact_cal_plus + force_sum_except_contact(mi, vi+0.5*fi/mi*ts)
    tqi_cal_plus = sum_tqji_plus_cal_steps

    return [fi_cal_plus, fi_plus, tqi_cal_plus, tqi_plus, fi_only_contact_cal_plus,fi_only_contact_plus]


def fjwi_plus_check_multistep_multicontact_fromcustom_inputvariablefunc(f_read, id_i, step1, step2, error_tolerence, method, version=0):

    [fi_cal_plus, fi_plus, tqi_cal_plus, tqi_plus, fi_only_contact_cal_plus, fi_only_contact_plus] = fjwi_plus_cal_multistep_multicontact_fromcustom(f_read, id_i, step1, step2, method, version)
    
    [fi_step_error_array, fi_errorindex] = step_error_array_with_index(step1, step2, fi_cal_plus, fi_plus, error_tolerence)
    [tqi_step_error_array, tqi_errorindex] = step_error_array_with_index(step1, step2, tqi_cal_plus, tqi_plus, error_tolerence)
    [fi_only_contact_step_error_array, fi_only_contact_errorindex] = step_error_array_with_index(step1, step2, fi_only_contact_cal_plus, fi_only_contact_plus, error_tolerence)
    
    steps_array = np.arange(step1+1, step2)
    steps_array = steps_array.reshape((-1,1))

    fi_cal_plus_in_error_step = (add_step_to_array(steps_array, fi_cal_plus))[fi_errorindex]
    fi_plus_in_error_step = (add_step_to_array(steps_array, fi_plus))[fi_errorindex]

    tqi_cal_plus_in_error_step = (add_step_to_array(steps_array, tqi_cal_plus))[tqi_errorindex]
    tqi_plus_in_error_step = (add_step_to_array(steps_array, tqi_plus))[tqi_errorindex]

    fi_only_contact_cal_plus_in_error_step = (add_step_to_array(steps_array, fi_only_contact_cal_plus))[fi_only_contact_errorindex]
    fi_only_contact_plus_in_error_step = (add_step_to_array(steps_array, fi_only_contact_plus))[fi_only_contact_errorindex]

    anwser = [fi_step_error_array,              fi_cal_plus_in_error_step,              fi_plus_in_error_step,
              tqi_step_error_array,             tqi_cal_plus_in_error_step,             tqi_plus_in_error_step,
              fi_only_contact_step_error_array, fi_only_contact_cal_plus_in_error_step, fi_only_contact_plus_in_error_step,]
    
    return anwser


def contact_check_multistep(f_read, id_i, step1, step2):
    df = pd.read_hdf(f_read, 'df')
    
    df_step = extract_dataframe(df, step1, step2)[['step', 'id', 'type', 'x', 'y', 'z',]]

    df_step_id_i = df_step.loc[df_step['id']==id_i]
    n_step1tostep2_id_i = len(df_step_id_i.index)
    n_step = step2 - step1
    if n_step1tostep2_id_i != n_step:
        steps_expected = np.arange(step1, step2)
        # only compare for same length part before array end
        min_len = min(n_step1tostep2_id_i, n_step)

        checksteps = df_step_id_i['step'].values[0:min_len] == steps_expected[0:min_len]
        if checksteps.all():
            id_first_wrong_step = min_len
        else:
            id_first_wrong_step = np.nonzero(~checksteps)[0][0]

        first_step_wrong = steps_expected[id_first_wrong_step]
        sys.exit('id_i is not in df read on step = {first_step_wrong}'.format(first_step_wrong=first_step_wrong))
    else:
        pass

    groups_byid = df_step.groupby(['id'])
    
    dfi = reindex_by_step(groups_byid.get_group(id_i), step1, step2)
    
    typei = dfi[['type']].values
    ri = radius_by_type(typei)
    xi = dfi[['x','y','z']].values

    id_list = groups_byid.groups.keys()

    id_j_list = [i for i in id_list if i !=id_i]

    id_j_list = np.asarray(id_j_list)

    number_idj = len(id_j_list)
    number_wall = len(walls)
    # overlap array time in first dim. idj in second dim

    ifoverlap_ij_array = np.empty((step2 - step1, number_idj), dtype=bool)
    ifoverlap_iw_array = np.empty((step2 - step1, number_wall), dtype=bool)

    for n, id_j in enumerate(id_j_list):

        dfj = reindex_by_step(groups_byid.get_group(id_j), step1, step2)
        typej = dfj[['type']].values
        rj = radius_by_type(typej)
        xj = dfj[['x','y','z']].values
        
        [ifoverlap_ij, overlapij_vector, overlapij_length] = overlapij(ri, rj, xi, xj)
        ifoverlap_ij_array[:,n:n+1] = ifoverlap_ij
    
    id_j_walls = id_j_list

    for n, walllist in enumerate(walls):

        id_j_walls = np.append(id_j_walls, -n-1)

        if walllist[0] == 'p':
            wall = wall_p_overlap_class(typei,xi,walllist[1], walllist[2], walllist[3], walllist[4])
        elif walllist[0] == 'cy':
            wall = wall_cy_overlap_class(typei,xi,walllist[1], walllist[2], walllist[3], walllist[4], walllist[5])
        else:
            print('walltype not p not cy')

        ifoverlap_iw = wall.ifoverlap()
        ifoverlap_iw_array[:,n:n+1] = ifoverlap_iw
    
    ifoverlap_ij_iw_array = np.concatenate((ifoverlap_ij_array, ifoverlap_iw_array), axis=1)
    ifoverlap_nantozero = reset_nan_to_zero(ifoverlap_ij_iw_array)
    diff_next = (ifoverlap_nantozero[:-1] != ifoverlap_nantozero[1:])
    index_diff_next = np.nonzero(diff_next)

    step_id_ifover_diffnext = np.empty((3, len(index_diff_next[0])), dtype=int)
    step_id_ifover_diffnext[2:3, :] = ifoverlap_nantozero[index_diff_next[0], index_diff_next[1]]
    step_id_ifover_diffnext[0:1, :] = index_diff_next[0] + step1
    step_id_ifover_diffnext[1:2, :] = id_j_walls[index_diff_next[1]]
    step_id_ifover_diffnext = np.transpose(step_id_ifover_diffnext)

    step_id_ifover_difflast = step_id_ifover_diffnext
    # change diffnext step to difflast step
    step_id_ifover_difflast[:,0] = step_id_ifover_difflast[:,0] + 1
    # change ifcontact in last step to ifcontact in next step
    step_id_ifover_difflast[:,2] = 1-step_id_ifover_difflast[:,2]

    # add initial contact
    initial_overlap = ifoverlap_nantozero[0:1] == True
    index_initial_overlap = np.nonzero(initial_overlap)
    id_j_wall_initial_overlap = id_j_walls[index_initial_overlap[1]].astype(int)
    step_id_ifover_initial = np.empty([len(id_j_wall_initial_overlap), 3], dtype=int)
    step_id_ifover_initial[:, 0:1] = step1
    step_id_ifover_initial[:, 1] = id_j_wall_initial_overlap
    step_id_ifover_initial[:, 2] = 1
    ifoverlap_ij_iw_array = ifoverlap_ij_iw_array.astype(int)

    return [step_id_ifover_difflast, ifoverlap_ij_iw_array, step_id_ifover_initial]


def contact_check_multistep_v1(f_read, id_i, step1, step2):

    df = pd.read_hdf(f_read, 'df')
    
    df_step = extract_dataframe(df, step1, step2)[['step', 'id', 'type', 'x', 'y', 'z',]]

    df_step_id_i = df_step.loc[df_step['id']==id_i]
    n_step1tostep2_id_i = len(df_step_id_i.index)
    n_step = step2 - step1
    if n_step1tostep2_id_i != n_step:
        steps_expected = np.arange(step1, step2)
        # only compare for same length part before array end
        min_len = min(n_step1tostep2_id_i, n_step)
        checksteps = df_step_id_i['step'].values[0:min_len] == steps_expected[0:min_len]
        if checksteps.all():
            id_first_wrong_step = min_len
        else:
            id_first_wrong_step = np.nonzero(~checksteps)[0][0]

        first_step_wrong = steps_expected[id_first_wrong_step]
        sys.exit('id_i is not in df read on step = {first_step_wrong}'.format(first_step_wrong=first_step_wrong))
    else:
        pass

    groups_byid = df_step.groupby(['id'])
    
    dfi = (groups_byid.get_group(id_i))
    
    typei = dfi[['type']].values
    xi = dfi[['x','y','z']].values

    id_list = groups_byid.groups.keys()

    id_j_list = [i for i in id_list if i !=id_i]

    id_j_list = np.asarray(id_j_list)

    number_idj = len(id_j_list)
    number_wall = len(walls)
    # overlap array time in first dim. idj in second dim

    ifoverlap_ij_array = np.full((step2 - step1, number_idj), False)
    ifoverlap_iw_array = np.full((step2 - step1, number_wall), False)

    for n, id_j in enumerate(id_j_list):

        dfj = (groups_byid.get_group(id_j))
        
        def selectindexfori(dfj):
            indexjfori = (dfj[['step']].values-step1).astype(int)
            indexjfori = indexjfori[:,0].T
            return indexjfori
        
            
        indexjfori = selectindexfori(dfj)
        typei = dfi[['type']].values[indexjfori]
        xi = dfi[['x','y','z']].values[indexjfori]

        typej = dfj[['type']].values
        xj = dfj[['x','y','z']].values

        jtoi = j_overlap_class(typei, xi, typej, xj)

        ifoverlap_ij = jtoi.ifoverlap()
        

        
        ifoverlap_ij_array[selectindexfori(dfj),n:n+1] = ifoverlap_ij
    
    id_j_walls = id_j_list

    for n, walllist in enumerate(walls):

        id_j_walls = np.append(id_j_walls, -n-1)
        if walllist[0] == 'p':
            wtoi = wall_p_overlap_class(typei,xi,walllist[1], walllist[2], walllist[3], walllist[4])
        elif walllist[0] == 'cy':
            wtoi = wall_cy_overlap_class(typei,xi,walllist[1], walllist[2], walllist[3], walllist[4], walllist[5])
        else:
            print('walltype not p not cy')

        ifoverlap_iw = wtoi.ifoverlap()
        

        ifoverlap_iw_array[:,n:n+1] = ifoverlap_iw

    
    ifoverlap_ij_iw_array = np.concatenate((ifoverlap_ij_array, ifoverlap_iw_array), axis=1)
    ifoverlap_nantozero = reset_nan_to_zero(ifoverlap_ij_iw_array)
    diff_next = (ifoverlap_ij_iw_array[:-1] != ifoverlap_ij_iw_array[1:])

    index_diff_next = np.nonzero(diff_next)

    step_id_ifover_diffnext = np.empty((3, len(index_diff_next[0])), dtype=int)
    step_id_ifover_diffnext[2:3, :] = ifoverlap_ij_iw_array[index_diff_next[0], index_diff_next[1]]
    step_id_ifover_diffnext[0:1, :] = index_diff_next[0] + step1
    step_id_ifover_diffnext[1:2, :] = id_j_walls[index_diff_next[1]]
    step_id_ifover_diffnext = np.transpose(step_id_ifover_diffnext)

    step_id_ifover_difflast = step_id_ifover_diffnext
    step_id_ifover_difflast[:,0] = step_id_ifover_difflast[:,0] + 1
    step_id_ifover_difflast[:,2] = 1-step_id_ifover_difflast[:,2]

    # add initial contact
    initial_overlap = ifoverlap_nantozero[0:1] == True
    index_initial_overlap = np.nonzero(initial_overlap)
    id_j_wall_initial_overlap = id_j_walls[index_initial_overlap[1]]
    step_id_ifover_initial = np.empty([len(id_j_wall_initial_overlap), 3], dtype=int)
    step_id_ifover_initial[:, 0:1] = step1
    step_id_ifover_initial[:, 1] = id_j_wall_initial_overlap
    step_id_ifover_initial[:, 2] = 1
    ifoverlap_ij_iw_array = ifoverlap_ij_iw_array.astype(int)

    return [step_id_ifover_difflast, ifoverlap_ij_iw_array, step_id_ifover_initial]


def steps_n_c_wall_less_1_n_c_j_less_1(f_read, id_i, step1, step2):
    
    [step_id_ifover_difflast, ifoverlap_ij_iw_array, step_id_ifover_initial] = contact_check_multistep(f_read, id_i, step1, step2)
    step_id_ifover = np.concatenate((step_id_ifover_initial, step_id_ifover_difflast), axis=0)
    steps = step_id_ifover[:,0]
    ids = step_id_ifover[:,1]
    ifoverlap = step_id_ifover[:,2]
    # ifoverlap_j change 0 to -1, see wall as no contact
    ifoverlap_j = 2*ifoverlap-1
    ifoverlap_j[ids<0] = 0
    # ifoverlap_wall change 0 to -1, see j as no contact
    ifoverlap_wall = 2*ifoverlap-1
    ifoverlap_wall[ids>0] = 0
    
    ncj = np.cumsum(ifoverlap_j, dtype=int) 
    ncwall = np.cumsum(ifoverlap_wall, dtype=int)
    
    contactchange_ncjis1_or_ncwall1 = np.append((steps[:-1] != steps[1:]), np.array([True])) & ((ncj==1) | (ncwall==1))
    step1_contactchange_ncjis1_or_ncwall1 = steps[contactchange_ncjis1_or_ncwall1]
    index = np.nonzero(contactchange_ncjis1_or_ncwall1)[0]
    if contactchange_ncjis1_or_ncwall1[-1]==True:
        # ignore last index and append step2
        step2_contactchange_ncjis1_or_ncwall1 = steps[index[:-1]+1]
        step2_contactchange_ncjis1_or_ncwall1 = np.append(step2_contactchange_ncjis1_or_ncwall1,(np.array([step2])))
    else:
        step2_contactchange_ncjis1_or_ncwall1 = steps[index+1]
    
    return [step1_contactchange_ncjis1_or_ncwall1, step2_contactchange_ncjis1_or_ncwall1]


def collect_contact_id_no_dup(f_read, id_i, step1, step2):
    [step_id_ifover_difflast, ifoverlap_ij_iw_array, step_id_ifover_initial] = contact_check_multistep(f_read, id_i, step1, step2)
    step_id_ifover = np.concatenate((step_id_ifover_initial, step_id_ifover_difflast), axis=0)
    contact_id_collection = (step_id_ifover[:,1])[step_id_ifover[:,2]==1]
    unique_ids = np.unique(contact_id_collection)
    return unique_ids


def number_contact_atom_id_collection(f_read, id_i, step1, step2):
    contact_id_collection_no_dup = collect_contact_id_no_dup(f_read, id_i, step1, step2)
    contact_atom_id_collection_no_dup = contact_id_collection_no_dup[contact_id_collection_no_dup > 0]
    number_contact_atom = len(contact_atom_id_collection_no_dup)
    return [number_contact_atom, contact_atom_id_collection_no_dup]


def number_contact_wall_id_collection(f_read, id_i, step1, step2):
    contact_id_collection_no_dup = collect_contact_id_no_dup(f_read, id_i, step1, step2)
    contact_wall_id_collection_no_dup = contact_id_collection_no_dup[contact_id_collection_no_dup < 0]
    number_contact_wall = len(contact_wall_id_collection_no_dup)
    
    return [number_contact_wall, contact_wall_id_collection_no_dup]


def number_contact_total_id_collection(f_read, id_i, step1, step2):
    
    contact_id_collection_no_dup = collect_contact_id_no_dup(f_read, id_i, step1, step2)
    number_contact_total = len(contact_id_collection_no_dup)
    
    return [number_contact_total, contact_id_collection_no_dup]


def fjwi_plus_cal_multistep_1contact_fromcustom(f_read, id_i, idj_or_idw, step1, step2, method):

    df = pd.read_hdf(f_read, 'df')
    
    df_step = extract_dataframe(df, step1, step2)

    df_step_id_i = df_step.loc[df_step['id']==id_i]
    n_step1tostep2_id_i = len(df_step_id_i.index)
    n_step = step2 - step1
    if n_step1tostep2_id_i != n_step:
        steps_expected = np.arange(step1, step2)
        # only compare for same length part before array end
        min_len = min(n_step1tostep2_id_i, n_step)
        checksteps = df_step_id_i['step'].values[0:min_len] == steps_expected[0:min_len]
        if checksteps.all():
            id_first_wrong_step = min_len
        else:
            id_first_wrong_step = np.nonzero(~checksteps)[0][0]
        first_step_wrong = steps_expected[id_first_wrong_step]
        sys.exit('id_i is not in df read on step = {first_step_wrong}'.format(first_step_wrong=first_step_wrong))
    else:
        pass

    groups_byid = df_step.groupby(['id'])
    
    dfi = reindex_by_step(groups_byid.get_group(id_i), step1, step2)
    
    typei = dfi[['type']].values[:-1]
    xi = dfi[['x','y','z']].values[:-1]
    vi = dfi[['vx','vy','vz']].values[:-1]
    fi = dfi[['fx','fy','fz']].values[:-1]
    fpair = dfi[['f_sum_pairforce[1]','f_sum_pairforce[2]','f_sum_pairforce[3]']].values[:-1]
    omi = dfi[['omegax','omegay','omegaz']].values[:-1]
    tqi = dfi[['tqx','tqy','tqz']].values[:-1]
    
    fi_plus = dfi[['fx','fy','fz']].values[1:]
    fpair_plus = dfi[['f_sum_pairforce[1]','f_sum_pairforce[2]','f_sum_pairforce[3]']].values[1:]
    
    tqi_plus = dfi[['tqx','tqy','tqz']].values[1:]

    [mi, ri, I_i, vi_half_pre, xi_minus_cal, omi_half_pre] = calculate_m_r_I_vh_xp_omh_pre(typei, density, fi, tqi, xi, vi, omi)

    force_not_contact = force_sum_except_contact(mi, vi_half_pre)
    force_not_contact_plus = force_sum_except_contact(mi, vi+0.5*fi/mi*ts)

    fjtoi = fpair
    fwtoi = fi - (fpair + force_not_contact)

    fjtoi_plus = fpair_plus
    fwtoi_plus = fi_plus - (fpair_plus + force_not_contact_plus)

    if idj_or_idw > 0:

        id_j = idj_or_idw

        dfj = reindex_by_step(groups_byid.get_group(id_j), step1, step2)

        typej = dfj[['type']].values[:-1]  
        xj = dfj[['x','y','z']].values[:-1]
        vj = dfj[['vx','vy','vz']].values[:-1]
        fj = dfj[['fx','fy','fz']].values[:-1]
        omj = dfj[['omegax','omegay','omegaz']].values[:-1]
        tqj = dfj[['tqx','tqy','tqz']].values[:-1]
        xj_plus = dfj[['x','y','z']].values[1:]
        fj_plus = dfj[['fx','fy','fz']].values[1:]
        tqj_plus = dfj[['tqx','tqy','tqz']].values[1:]

        [mj, rj, I_j, vj_half_pre, xj_minus_cal, omj_half_pre] = calculate_m_r_I_vh_xp_omh_pre(typej, density, fj, tqj, xj, vj, omj)

        meff = mi*mj/(mi + mj)
        f_jorw_to_i_total = fjtoi
        f_jorw_to_i_total_plus = fjtoi_plus

        xij = xi - xj
        [ifoverlap, overlapij_vector, overlapij_length] = overlapij(ri, rj, xi, xj)

    if idj_or_idw < 0:

        idw = -idj_or_idw-1
        walllist = walls[idw]

        wall = create_wall_class_from_walllist(
                walllist,
                typei,xi,vi,fi,omi,tqi, 
                np.zeros(3), np.zeros(3), np.zeros(3), method,
                )

        [rj, vj_half_pre, xj_minus_cal, omj_half_pre] = calculate_wall_r_vh_xp_omh_pre(wall)
        meff = wall.meff()
        f_jorw_to_i_total = fwtoi
        f_jorw_to_i_total_plus = fwtoi_plus
        xij = xi - wall.nearest_point()
        ifoverlap = wall.ifoverlap()
        
        overlapij_vector = wall.overlapij_vector()

    vij_half_pre = vi_half_pre - vj_half_pre
    vijn_half_pre = projection_vector(vij_half_pre, xij)
    vijt_half_no_om_pre = verticle_vector(vij_half_pre, xij)
    vijt_half_pre = vijt_half_no_om_pre + ri*np.cross(omi_half_pre, unit(-xij)) - rj*np.cross(omj_half_pre, unit(xij))

    [fnk,fngamma,ftk,ftgamma] = cal_f_no_history(overlapij_vector, ifoverlap, meff, vijn_half_pre, vijt_half_pre)
    fjin = fnk + fngamma
    friction_ratio = np.divide(length(f_jorw_to_i_total - fjin), length(fjin), out=np.zeros_like(length(fjin)), where=length(fjin)!=0)
    ftk_include_his = f_jorw_to_i_total - fjin - ftgamma
    total_displacement = -ftk_include_his/kt
    total_length = length(total_displacement)

    if idj_or_idw > 0:
        [ftk_include_his_plus_cal, total_displacement_plus_cal, fjit_plus_cal, fjin_plus_cal, tqji_plus_cal, total_length_plus_cal] = (
            fji_plus_cal(typei, typej, xi, xj, vi, vj, fi, fj, ftk_include_his, total_length, total_displacement, omi, omj, tqi, tqj, method=method)
        )

    if idj_or_idw < 0:
        [ftk_include_his_plus_cal, total_displacement_plus_cal, fjit_plus_cal, fjin_plus_cal, tqji_plus_cal, total_length_plus_cal] = (
            fwi_plus_cal(typei, wall, xi, vi, fi, ftk_include_his, total_length, total_displacement, omi, tqi, method=method)
        )
    
    fji_cal = fjit_plus_cal + fjin_plus_cal

    return [fji_cal, f_jorw_to_i_total_plus, tqji_plus_cal, tqi_plus, fjit_plus_cal, fjin_plus_cal]

    
def fjwi_plus_check_multistep_1contact_fromcustom(f_read, id_i, step1, step2, error_tolerence, method, jorw):

    if jorw == 'j':
        [n, id_collection] = number_contact_atom_id_collection(f_read, id_i, step1, step2)
    elif jorw == 'w':
        [n, id_collection] = number_contact_wall_id_collection(f_read, id_i, step1, step2)
    else:
        sys.exit("jorw not j not w")

    if n > 1:
        sys.exit("number of contact is {n} for ".format(n=n) + jorw + ' should be 0 or 1')
    elif n == 1:
        [fji_cal, fji_plus, tqji_plus_cal, tqi_plus, fjit_plus_cal, fjin_plus_cal] = fjwi_plus_cal_multistep_1contact_fromcustom(f_read, id_i, id_collection[0], step1, step2, method)
        [f_step_error_array, f_errorindex] = step_error_array_with_index(step1, step2, fji_cal, fji_plus, error_tolerence)
        [tq_step_error_array, tq_errorindex] = step_error_array_with_index(step1, step2, tqji_plus_cal, tqi_plus, error_tolerence)
        
        steps_array = np.arange(step1+1, step2)
        steps_array = steps_array.reshape((-1,1))

        step_fji_cal = add_step_to_array(steps_array, fji_cal)
        step_fji_plus = add_step_to_array(steps_array, fji_plus)
        fji_cal_in_error_step = step_fji_cal[f_errorindex]
        fji_plus_in_error_step = step_fji_plus[f_errorindex]
    elif n==0:
        print("number of contact is {n} for ".format(n=n) + jorw + ' so not check')
    else:
        sys.exit("number of contact is {n} for ".format(n=n) + jorw + ' should >=0')

    return [f_step_error_array, fji_cal_in_error_step, fji_plus_in_error_step, id_collection[0]]


def get_vijt(typei, typej, xi, xj, vi, vj, omi, omj):
    
    ri = radius_by_type(typei)
    rj = radius_by_type(typej)
    xij = xi - xj
    vij = vi - vj
    vijt_no_om = verticle_vector(vij, xij)

    return vijt_no_om + ri*np.cross(omi, unit(-xij)) - rj*np.cross(omj, unit(xij))


def get_vijn(xi, xj, vi, vj):

    xij = xi - xj
    vij = vi - vj

    return projection_vector(vij, xij)


def get_viwt(typei, wall, xi, vi, omi):
    
    ri = radius_by_type(typei)
    rj = wall.r
    vj = wall.vw
    omj = wall.omw
    xij = xi - wall.nearest_point()
    vij = vi - vj
    
    vijt_no_om = verticle_vector(vij, xij)
    
    return vijt_no_om + ri*np.cross(omi, unit(-xij)) - rj*np.cross(omj, unit(xij))


def get_viwn(typei, wall, xi, vi):

    vj = wall.vw
    xij = xi - wall.nearest_point()
    vij = vi - vj

    return projection_vector(vij, xij)



class manysteps(object):
    
    def __init__(self, f_read, id_i, idj_or_idw, step1, step2, method):
        self.f_read = f_read
        self.id_i = id_i
        self.idj_or_idw = idj_or_idw
        self.step1 = step1
        self.step2 = step2
        self.method = method
        df = pd.read_hdf(self.f_read, 'df')
        df_step = extract_dataframe(df, self.step1, self.step2)
        self.groups_byid = df_step.groupby(['id'])
        self.dfi = reindex_by_step(self.groups_byid.get_group(self.id_i), self.step1, self.step2)
        [self.typei, self.xi, self.vi, self.fi, self.omi, self.tqi] = get_type_x_v_f_om_tq_from_df(self.dfi, self.step1, self.step2)
    
    def gravity(self):
        mi = mass_by_type(self.typei, density)
        return gravity(mi)

    def force_sum_except_contact(self):
        mi = mass_by_type(self.typei, density)
        return force_sum_except_contact(mi, self.vi+0.5*self.fi/mi*ts)


    def vh_omh_pre(self, type, density, f, tq, x, v, om):
    
        m = mass_by_type(type, density)
        r = radius_by_type(type)
        I = 2/5*m*r**2
        a = f/m
        alpha = tq/I
        [v_half_pre, x_minus_cal]=update_v_x_pre(x,v,a)
        result5=update_v_x_pre(0,om,alpha)
        om_half_pre = result5[0]

        return [v_half_pre, om_half_pre]


    def work_ftfn_many_steps(self, vijt, vijn, f_many_steps):

        [fnk,fngamma,ftk_include_his,ftgamma] = f_many_steps
        fn = fnk + fngamma
        ft = ftk_include_his + ftgamma
        
        def innerproduct(a,b):
            return np.sum(a*b, axis=-1)

        work_ft = innerproduct(ft, vijt)*ts
        work_fn = innerproduct(fn, vijn)*ts

        return [work_ft, work_fn]

    def work_ftgammafngamma_many_steps(self, vijt, vijn, f_many_steps):

        [fnk,fngamma,ftk_include_his,ftgamma] = f_many_steps

        def innerproduct(a,b):
            return np.sum(a*b, axis=-1)

        work_ftgamma = innerproduct(ftgamma,vijt)*ts
        work_fngamma = innerproduct(fngamma,vijn)*ts
        return [work_ftgamma, work_fngamma]

    def work_ftkfnk_many_steps(self, vijt, vijn, f_many_steps):

        [fnk,fngamma,ftk_include_his,ftgamma] = f_many_steps

        def innerproduct(a,b):
            return np.sum(a*b, axis=-1)

        work_ftk = innerproduct(ftk_include_his,vijt)*ts
        work_fnk = innerproduct(fnk,vijn)*ts
        return [work_ftk, work_fnk]


class manysteps_idj(manysteps):

    def __init__(self, f_read, id_i, idj_or_idw, step1, step2, method):
        super().__init__(f_read, id_i, idj_or_idw, step1, step2, method)
        self.fj = reindex_by_step(self.groups_byid.get_group(self.idj_or_idw), self.step1, self.step2)
        [self.typej, self.xj, self.vj, self.fj, self.omj, self.tqj] = get_type_x_v_f_om_tq_from_df(self.fj, self.step1, self.step2)
    def vijt(self):
        return get_vijt(self.typei, self.typej, self.xi, self.xj, self.vi, self.vj, self.omi, self.omj)

    def vijt_contactpoint(self):
        ri = radius_by_type(self.typei)
        rj = radius_by_type(self.typej)
        [ifoverlap, overlapij_vector, overlap_length] = overlapij(ri, rj, self.xi, self.xj)
        return ifoverlap*(rj/(ri+rj)*get_vijt(self.typei, self.typej, self.xi, self.xj, self.vi, self.vj, np.zeros_like(self.vi), np.zeros_like(self.vj)))

    def vijn(self):
        return get_vijn(self.xi, self.xj, self.vi, self.vj)
    
    def vijt_half_pre(self):
        [mi, ri, I_i, vi_half_pre, xi_minus_cal, omi_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typei, density, self.fi, self.tqi, self.xi, self.vi, self.omi)
        [mj, rj, I_j, vj_half_pre, xj_minus_cal, omj_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typej, density, self.fj, self.tqj, self.xj, self.vj, self.omj)
        return get_vijt(self.typei, self.typej, self.xi, self.xj, vi_half_pre, vj_half_pre, omi_half_pre, omj_half_pre)

    def vijn_half_pre(self):
        [mi, ri, I_i, vi_half_pre, xi_minus_cal, omi_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typei, density, self.fi, self.tqi, self.xi, self.vi, self.omi)
        [mj, rj, I_j, vj_half_pre, xj_minus_cal, omj_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typej, density, self.fj, self.tqj, self.xj, self.vj, self.omj)
        return get_vijn(self.xi, self.xj, vi_half_pre, vj_half_pre)
    
    def ifoverlap(self):
        ri = radius_by_type(self.typei)
        rj = radius_by_type(self.typej)
        [ifoverlap, overlapij_vector, overlap_length] = overlapij(ri, rj, self.xi, self.xj)
        return ifoverlap

    def overlap_length(self):
        ri = radius_by_type(self.typei)
        rj = radius_by_type(self.typej)
        [ifoverlap, overlapij_vector, overlap_length] = overlapij(ri, rj, self.xi, self.xj)
        return overlap_length

    def f_jorw_to_i_total(self):
        # only 1 contact can use this function
        fpair = self.dfi[['f_sum_pairforce[1]','f_sum_pairforce[2]','f_sum_pairforce[3]']].values
        fjtoi = fpair
        return fjtoi

    def xij(self):
        return self.xi - self.xj

    def f_many_steps(self):
        # only 1 contact can use this function

        [mi, ri, I_i, vi_half_pre, xi_minus_cal, omi_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typei, density, self.fi, self.tqi, self.xi, self.vi, self.omi)
        id_j = self.idj_or_idw
        [mj, rj, I_j, vj_half_pre, xj_minus_cal, omj_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typej, density, self.fj, self.tqj, self.xj, self.vj, self.omj)
        meff = mi*mj/(mi + mj)
        [ifoverlap, overlapij_vector, overlapij_length] = overlapij(ri, rj, self.xi, self.xj)

        [fnk,fngamma,ftk,ftgamma] = cal_f_no_history(overlapij_vector, ifoverlap, meff, self.vijn_half_pre(), self.vijt_half_pre())
        fjin = fnk + fngamma
        friction_ratio = np.divide(length(self.f_jorw_to_i_total() - fjin), length(fjin), out=np.zeros_like(length(fjin)), where=length(fjin)!=0)
        ftk_include_his = self.f_jorw_to_i_total() - fjin - ftgamma

        return [fnk,fngamma,ftk_include_his,ftgamma]
    
    def work_ftfn_many_steps(self):
        return super().work_ftfn_many_steps(self.vijt(), self.vijn(), self.f_many_steps())

    def work_ftgammafngamma_many_steps(self):
        return super().work_ftgammafngamma_many_steps(self.vijt(), self.vijn(), self.f_many_steps())

    def work_ftkfnk_many_steps(self):
        return super().work_ftkfnk_many_steps(self.vijt(), self.vijn(), self.f_many_steps())

    def cumsum_work(self):
        [work_ft, work_fn] = self.work_ftfn_many_steps()
        sum_ftwork = np.cumsum(work_ft, axis=0)
        sum_fnwork = np.cumsum(work_fn, axis=0)
        return [sum_ftwork, sum_fnwork]

class manysteps_wall(manysteps):

    def __init__(self, f_read, id_i, idj_or_idw, step1, step2, method):
        super().__init__(f_read, id_i, idj_or_idw, step1, step2, method)

        id_walls = -self.idj_or_idw-1
        walllist = walls[id_walls]

        self.wall = create_wall_class_from_walllist(
                walllist,
                self.typei,self.xi,self.vi,self.fi,self.omi,self.tqi, 
                np.zeros(3), np.zeros(3), np.zeros(3), self.method,
                )
    def vijt(self):
        return get_viwt(self.typei, self.wall, self.xi, self.vi, self.omi)

    def vijt_contactpoint(self):
        ifoverlap = self.wall.ifoverlap()
        return ifoverlap*(get_viwt(self.typei, self.wall, self.xi, self.vi, np.zeros_like(self.vi)))


    def vijn(self):
        return get_viwn(self.typei, self.wall, self.xi, self.vi)

    def vijt_half_pre(self):
        [mi, ri, I_i, vi_half_pre, xi_minus_cal, omi_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typei, density, self.fi, self.tqi, self.xi, self.vi, self.omi)
        return get_viwt(self.typei, self.wall, self.xi, vi_half_pre, omi_half_pre)
    def vijn_half_pre(self):    
        [mi, ri, I_i, vi_half_pre, xi_minus_cal, omi_half_pre] = calculate_m_r_I_vh_xp_omh_pre(self.typei, density, self.fi, self.tqi, self.xi, self.vi, self.omi)
        return get_viwn(self.typei, self.wall, self.xi, vi_half_pre)

    def ifoverlap(self):
        return self.wall.ifoverlap()

    def overlap_length(self):
        return self.wall.overlap_length()

    def f_jorw_to_i_total(self):
        # only 1 contact can use this function
        fpair = self.dfi[['f_sum_pairforce[1]','f_sum_pairforce[2]','f_sum_pairforce[3]']].values
        mi = mass_by_type(self.typei, density)
        force_not_contact = force_sum_except_contact(mi, self.vi+0.5*self.fi/mi*ts)
        fwtoi = self.fi - (fpair + force_not_contact)
        return fwtoi

    def xij(self):
        return self.xi - self.wall.nearest_point()


    def f_many_steps(self):
        # only 1 contact can use this function

        
        meff = self.wall.meff()
        ifoverlap = self.wall.ifoverlap()
        
        overlapij_vector = self.wall.overlapij_vector()

        [fnk,fngamma,ftk,ftgamma] = cal_f_no_history(overlapij_vector, ifoverlap, meff, self.vijn_half_pre(), self.vijt_half_pre())
        fjin = fnk + fngamma
        friction_ratio = np.divide(length(self.f_jorw_to_i_total() - fjin), length(fjin), out=np.zeros_like(length(fjin)), where=length(fjin)!=0)
        ftk_include_his = self.f_jorw_to_i_total() - fjin - ftgamma

        return [fnk,fngamma,ftk_include_his,ftgamma]

    def work_ftfn_many_steps(self):
        return super().work_ftfn_many_steps(self.vijt(), self.vijn(), self.f_many_steps())

    def work_ftgammafngamma_many_steps(self):
        return super().work_ftgammafngamma_many_steps(self.vijt(), self.vijn(), self.f_many_steps())

    def work_ftkfnk_many_steps(self):
        return super().work_ftkfnk_many_steps(self.vijt(), self.vijn(), self.f_many_steps())

    def cumsum_work(self):
        [work_ft, work_fn] = self.work_ftfn_many_steps()
        sum_ftwork = np.cumsum(work_ft, axis=0)
        sum_fnwork = np.cumsum(work_fn, axis=0)
        return [sum_ftwork, sum_fnwork]


def sum_results_all_contacts(
    func,
    f_read, id_i, idj_or_idw_list, step1, step2, error_tolerence, method,
    ):

    [number_contact_total, contact_id_collection_no_dup] = number_contact_total_id_collection(f_read, id_i, step1, step2)

    sum = 0
    for id in contact_id_collection_no_dup:
        sum = sum + func(f_read, id_i, id, step1, step2, error_tolerence, method)

    return sum


def max_id_step(f_read, step1, step2):

    df = pd.read_hdf(f_read, 'df')
    df_step = extract_dataframe(df, step1, step2)
    
    ids = df_step[['id']].values
    diff_next = (ids[:-1] != ids[1:])
    index_diff_next = np.nonzero(diff_next)[0]
    index_diff_last = index_diff_next + 1
    step_id_difflast = np.empty( (len(index_diff_last), 2), dtype=int)
    step_id_difflast[:, 0] = index_diff_last + step1
    step_id_difflast[:, 1:2] = ids[index_diff_last]
    step_id_initial = np.array([[step1, ids[0,0]]], dtype=int)
    step_id = np.concatenate((step_id_initial, step_id_difflast), axis=0)

    return step_id

def find_lost_id(f_read, step1, step2):

    df = pd.read_hdf(f_read, 'df')
    df = df[['step', 'id']].astype('int64')
    df_step = extract_dataframe(df, step1, step2)
    df_id_lastexiststep = df_step.groupby(['id']).max()
    df_id_lastexiststep = df_id_lastexiststep.loc[df_id_lastexiststep['step'] != int(step2-1)]

    return df_id_lastexiststep
