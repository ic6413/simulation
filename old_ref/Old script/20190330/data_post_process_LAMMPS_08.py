# data_post_process_LAMMPS
# import
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
import LAMMPS_attributes as La

# ====================================== import variable
# atom radius
dp0 = La.dp0
rp0 = dp0/2
density = La.density
# intersection point
z_zplane = La.z_zplane
r_in = La.r_in
r_out = La.r_out
# parameter
mu = La.mu
kn = La.kn 
kt = La.kt 
gamma_n = La.gamma_n
gamma_t = La.gamma_t

type_radius_list = La.type_radius_list
#type_radius_list = [[1,1]]
# zplane1 zplane2....
zplane_list = La.zplane_list
# zcylinder1 zcylinder2....
zcylinder_list = La.zcylinder_list
# timestep
ts = La.ts
# ====================================== end import variable

# define function for get radius of each atom
def radius_by_type(r0, type, type_radius_list):
    rp = 0
    for type_radius in type_radius_list:
        rp = rp + r0*type_radius[1]*(type_radius[0] == type)
    return rp
def mass_by_type(r0, type, type_radius_list, density):
    return 4/3*np.pi*radius_by_type(r0, type, type_radius_list)**3*density
# the closest point on zplane for atom
def point_zplane(position, z_zplane):
	return position*[1, 1, 0] + [0, 0, z_zplane]
# the closest point on zcylinder for atom
def point_zcylinder(position, r):
    R = np.sqrt(np.sum(position**2*[1, 1, 0], axis=-1))
    return position*np.stack([R]*3, axis=-1)
# np.square np.inner cross
def length(vector):
	return np.sqrt(np.sum(vector**2, axis=-1))
def unit(vector):
	return vector/np.stack([length(vector)]*3,axis=-1)
def projection_scalar(v, n):
	return np.sum(v*unit(n), axis=-1)
def projection_vector(v, n):
	return unit(n)*np.stack([projection_scalar(v, n)]*3,axis=-1)
def verticle_vector(v, n):
	return v - projection_vector(v, n)
def costheta(v, n):
    return projection_scalar(v, n)/length(v)

# r_p_wall_vector, wall_position is radius or plane height
def r_p_wall_vector(position, walltype, wall_position):
    if walltype == 'zp':
        return position - point_zplane(position, wall_position)
    elif walltype == 'zcy':
        return position - point_zcylinder(position, wall_position)

def fji_plus_check(typei, typej, xi, xj, vi, vj, fi, fj, previous_ktforce, total_length, total_displacement, omi, omj, tqi, tqj, xi_plus, xj_plus, fji_plus, tqji_plus):
    
    mi = mass_by_type(rp0, np.array([typei, typei, typei]), type_radius_list, density)
    mj = mass_by_type(rp0, np.array([typej, typej, typej]), type_radius_list, density)
    ri = radius_by_type(rp0, np.array([typei, typei, typei]), type_radius_list)
    rj = radius_by_type(rp0, np.array([typej, typej, typej]), type_radius_list)
    xij_plus = xi_plus - xj_plus
    overlapij_vector_plus = unit(-xij_plus)*(ri + rj - length(xij_plus))*(ri + rj - length(xij_plus) >= 0)

    fjin_plus = projection_vector(fji_plus, xij_plus)
    fjit_plus = verticle_vector(fji_plus, xij_plus)

    meff = mi*mj/(mi + mj)
    ai = fi/mi
    aj = fj/mj
    vi_half = vi + ai*ts/2
    vj_half = vj + aj*ts/2
    vij_half = vi_half - vj_half
    xij_cal = xi - xj
    xi_plus_cal = xi + vi_half*ts
    xj_plus_cal = xj + vj_half*ts
    xij_plus_cal = xi_plus_cal - xj_plus_cal
    overlapij_vector_plus_cal = unit(-xij_plus_cal)*(ri + rj - length(xij_plus_cal))*(ri + rj - length(xij_plus) >= 0)
    if (overlapij_vector_plus_cal == 0).all():
        sys.exit('no overlap')
    xi_plus_error = (xi_plus_cal - xi_plus)/xi_plus
    xj_plus_error = (xj_plus_cal - xj_plus)/xj_plus
    xij_plus_error = (xij_plus_cal - xij_plus)/xij_plus
    overlapij_vector_plus_error = (overlapij_vector_plus_cal - overlapij_vector_plus)/overlapij_vector_plus
    vijn_half = projection_vector(vij_half, xij_plus_cal)
    vijt_half_no_om = verticle_vector(vij_half, xij_plus_cal)
    I_i = 2/5*mi*ri**2
    I_j = 2/5*mj*rj**2
    omi_half = omi + tqi/I_i*ts/2
    omj_half = omj + tqj/I_j*ts/2
    vijt_half = vijt_half_no_om + ri*np.cross(omi_half, unit(-xij_plus_cal)) - rj*np.cross(omj_half, unit(xij_plus_cal))
    fnk = -kn*unit(xij_plus_cal)*(length(xij_plus_cal)-ri-rj)
    fngamma = - meff*gamma_n*vijn_half
    ftk = -kt*(vijt_half*ts)
    ftgamma = -meff*gamma_t*vijt_half
    fjin_plus_cal = fnk + fngamma
    displacement_singlestep = vijt_half*ts

    ftk_include_his = ftk + history_force(previous_ktforce, total_length, total_displacement, vijt_half, xij_cal, xij_plus_cal)
    fjit_plus_cal = ftgamma + ftk_include_his
    if length(fjit_plus_cal) > mu*length(fjin_plus_cal):
        fjit_plus_cal = mu*length(fjin_plus_cal)*unit(fjit_plus_cal)
    fji_plus_cal = fjin_plus_cal + fjit_plus_cal
    tqji_plus_cal = np.cross(ri*unit(-xij_plus_cal),fjit_plus_cal)
    fjin_error = (fjin_plus_cal - fjin_plus)/fjin_plus
    fjit_error = (fjit_plus_cal - fjit_plus)/fjit_plus
    fji_error = (fji_plus_cal - fji_plus)/fji_plus 
    tqji_error = (tqji_plus_cal - tqji_plus)/tqji_plus

    printdic = {}
    printdic['xj_plus'] = xj_plus
    printdic['xj_plus_cal'] = xj_plus_cal
    printdic["xi_plus_error"] = xi_plus_error
    printdic["xj_plus_error"] = xj_plus_error
    printdic["xij_plus_error"] = xij_plus_error
    printdic["overlapij_vector_plus_error"] = overlapij_vector_plus_error
    printdic["vijt_half_no_om"] = vijt_half_no_om
    printdic["vijt_half"] = vijt_half
    printdic["fjin_error"] = fjin_error
    printdic["fjit_error"] = fjit_error
    printdic["fji_error"] = fji_error
    printdic["tqji_error"] = tqji_error
    printdic["fji_plus_cal"] = fji_plus_cal
    printdic["length(fjit_plus_cal)/length(fjin_plus_cal)"] = length(fjit_plus_cal)/length(fjin_plus_cal)
    printdic["overlapij_vector_plus_cal/xij_plus_cal"] = overlapij_vector_plus_cal/xij_plus_cal
    printdic["fjit_plus - fjit_plus_cal"] = fjit_plus - fjit_plus_cal
    printdic["costheta(fjit_plus, xij_plus)"] = costheta(fjit_plus, xij_plus)
    printdic["costheta(fjit_plus_cal, xij_plus_cal)"] = costheta(fjit_plus_cal, xij_plus_cal)
    printdic["costheta(fjit_plus_cal, fjit_plus)"] = costheta(fjit_plus_cal, fjit_plus)
    printdic["length(fji_plus_cal)/length(fji_plus)"] = length(fji_plus_cal)/length(fji_plus)
    pprint.pprint (printdic)
    pprint.pprint (printdic, open("printdic",'a'))

    return [ftk_include_his, displacement_singlestep, fjit_error]

def history_force(previous_ktforce, total_length, total_displacement, vijt_half, xij_cal, xij_plus_cal):
    p_scalar = projection_scalar(total_displacement, vijt_half)
    n1 = unit(xij_cal)
    n2 = unit(xij_plus_cal)
    dn = n2 - n1
    dn_t1 = verticle_vector(dn, n1)
    dn_t2 = verticle_vector(dn, n2)
    history_ktforce_para = projection_scalar(previous_ktforce, dn_t1)*unit(dn_t2)
    history_ktforce_verti = verticle_vector(previous_ktforce, dn_t1)
    
    # wrong when direction change#===============
    #return -kt*total_length*unit(vijt_half)
    #return -kt*total_length*unit(vijt_half)*p_scalar/abs(p_scalar)
    #return -kt*length(total_displacement)*unit(vijt_half)
    return -kt*total_length*unit(vijt_half)*p_scalar/abs(p_scalar)
    
    # wrong definitely#==========================
    #return -kt*length(total_displacement)
    # ===========================================


    # correct================================================================================
    #return -kt*projection_scalar(total_displacement, vijt_half)*unit(vijt_half)
    #return -kt*length(total_displacement)*unit(vijt_half)*p_scalar/abs(p_scalar)
    #return -kt*verticle_vector(total_displacement, xij_plus_cal)
    
    # covariant  =============================================
    #return history_ktforce_para + history_ktforce_verti

def fji_plus_check_multistep(f_read, id_i, id_j, step1, step2, error_tolerence):
    df = pd.read_hdf(f_read, 'df')
    ftk_include_his = 0
    displacement_singlestep = 0
    total_length = 0
    total_displacement = 0
    total_project_length = 0
    errorsteplist = []
    for step in range(step1, step2):
        print(step)
        typei = df.loc[(df['step'] == step) & (df['id_i'] == id_i), ['type_i']].values[0,0]
        xi = df.loc[(df['step'] == step) & (df['id_i'] == id_i), ['x_i','y_i','z_i']].values[0]
        vi = df.loc[(df['step'] == step) & (df['id_i'] == id_i), ['vx_i','vy_i','vz_i']].values[0]
        fi = df.loc[(df['step'] == step) & (df['id_i'] == id_i), ['fx_i','fy_i','fz_i']].values[0]
        omi = df.loc[(df['step'] == step) & (df['id_i'] == id_i), ['omegax_i','omegay_i','omegaz_i']].values[0]
        tqi = df.loc[(df['step'] == step) & (df['id_i'] == id_i), ['tqx_i','tqy_i','tqz_i']].values[0]
        
        typej = df.loc[(df['step'] == step) & (df['id_i'] == id_j), ['type_i']].values[0,0]
        xj = df.loc[(df['step'] == step) & (df['id_i'] == id_j), ['x_i','y_i','z_i']].values[0]
        vj = df.loc[(df['step'] == step) & (df['id_i'] == id_j), ['vx_i','vy_i','vz_i']].values[0]
        fj = df.loc[(df['step'] == step) & (df['id_i'] == id_j), ['fx_i','fy_i','fz_i']].values[0]
        omj = df.loc[(df['step'] == step) & (df['id_i'] == id_j), ['omegax_i','omegay_i','omegaz_i']].values[0]
        tqj = df.loc[(df['step'] == step) & (df['id_i'] == id_j), ['tqx_i','tqy_i','tqz_i']].values[0]
        
        fji_plus = df.loc[(df['step'] == step+1) & (df['id_i'] == id_i), ['f_force_pair[1]_i','f_force_pair[2]_i','f_force_pair[3]_i']].values[0]
        tqji_plus = df.loc[(df['step'] == step+1) & (df['id_i'] == id_i), ['tqx_i','tqy_i','tqz_i']].values[0]
        xi_plus = df.loc[(df['step'] == step+1) & (df['id_i'] == id_i), ['x_i','y_i','z_i']].values[0]
        xj_plus = df.loc[(df['step'] == step+1) & (df['id_i'] == id_j), ['x_i','y_i','z_i']].values[0]

        previous_ktforce = ftk_include_his
        total_displacement = total_displacement + displacement_singlestep
        total_length = total_length + length(displacement_singlestep)
        total_project_length = total_project_length + length(displacement_singlestep)
        [ftk_include_his, displacement_singlestep, fjit_error] = fji_plus_check(typei, typej, xi, xj, vi, vj, fi, fj,previous_ktforce, total_length, total_displacement, omi, omj, tqi, tqj, xi_plus, xj_plus, fji_plus, tqji_plus)
        if (abs(fjit_error) > error_tolerence).any():
            errorsteplist.append([step, fjit_error])
    return errorsteplist