#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataSciece.changeDirOnImportExport setting
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'input'))
	print(os.getcwd())
except:
	pass

#%%
import re
import time
from itertools import chain
from itertools import repeat
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D

#%%
topdirectory = '/home/ic6413/mylammps/myfiles/run/'
runfolder = '0103_mu05_nofreeze_traceid15556/'
dumppath = topdirectory + runfolder + 'output/dump/'
pythonpath = topdirectory + runfolder + 'output/python/'

#%%
## build wall contact profile from dump custom15556
## read trace hdf5 file
# timer
start = time.process_time()
# read
filename = 'trace15556.h5'
df = pd.read_hdf( pythonpath + filename, 'df')
# error test
error_tolerence = 10**-5
# atom radius
dp0 = 0.00083
rp0 = dp0/2
# intersection point
z_zplane = 0*dp0
r_in = 37*dp0
r_out = (37+16)*dp0
# parameter
mu = 0.5
kn = 8837.32631448687 
kt = 2524.95037556768 
gamma_n = 34445.5603308471 
gamma_t = 17222.7801654236

def point_zp(position):
	return position*[1, 1, 0] + [0, 0, z_zplane]

def point_zc_in(position):
	r = r_in
	R = np.sqrt(np.sum(position**2*[1, 1, 0], axis=-1))
	return position*np.stack([R]*3, axis=-1)

def point_zc_out(position):
	r = r_out
	R = np.sqrt(np.sum(position**2*[1, 1, 0], axis=-1))
	return position*np.stack([R]*3, axis=-1)

wallpoint = point_zc_in #point_zp, point_zc_in, point_zc_out

# np.square np.inner cross
def length(vector):
	return np.sqrt(np.sum(vector**2, axis=-1))
def unit(vector):
	return vector/np.stack([length(vector)]*3,axis=-1)
def projection_length(v, n):
	return np.sum(v*unit(n), axis=-1)
def projection_vector(v, n):
	return unit(n)*np.stack([projection_length(v, n)]*3,axis=-1)
def verticle_vector(v, n):
	return v - projection_vector(v, n)

# r_p_wall_vector
def r_p_wall_vector(position):
	return position - wallpoint(position)
# calculate radius from type
def radius_type(type):
	return rp0*(0.9*(type == 1) + 1.0*(type == 2) + 1.1*(type == 3))
def mass_type(type):
    return 4/3*np.pi*radius_type(type)**3*2500
# ifoverlap
def ifoverlap(position, type):
	r_p_wall = length(position - wallpoint(position))
	return r_p_wall <= radius_type(type)

#header = ['id_i', 'step', 'type_i', 'v_EKP_atom_i', 'v_KEtr_atom_i', 'c_KEt_atom_i', 'c_KEr_atom_i',
# 'x_i', 'y_i', 'z_i', 'vx_i', 'vy_i', 'vz_i', 'fx_i', 'fy_i', 'fz_i',
# 'omegax_i', 'omegay_i', 'omegaz_i', 'tqx_i', 'tqy_i', 'tqz_i',
# 'f_force_all[1]_i', 'f_force_all[2]_i', 'f_force_all[3]_i', 
# 'f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i', 
# 'c_n_contact_i', 'v_check_nb_all_i',
# 'index', 'id_j', 'type_j', 
# 'distance', 
# 'fjinx', 'fjiny', 'fjinz', 'fjin', 
# 'fjitx', 'fjity', 'fjitz', 'fjit', 
# 'vijnx', 'vijny', 'vijnz',
# 'vijtx', 'vijty', 'vijtz', 
# 'c_KEt_atom_j', 'c_KEr_atom_j', 'x_j', 'y_j', 'z_j', 'vx_j', 'vy_j', 'vz_j', 'omegax_j', 'omegay_j', 'omegaz_j']

begin = 0
end = 300000

r_i = radius_type(df['type_i'].values)
r_j = radius_type(df['type_j'].values)

m_i = mass_type(df['type_i'].values)
m_j = mass_type(df['type_j'].values)

ts = 10**-6

meff = m_i*m_j/(m_i + m_j)

xyz_vector_i = df[['x_i', 'y_i', 'z_i']].values
xyz_vector_j = df[['x_j', 'y_j', 'z_j']].values
xyz_vector_ij = xyz_vector_i - xyz_vector_j

xyz_ij = length(xyz_vector_ij)
xyz_ij_unit = unit(xyz_vector_ij)
xyz_vector_ji = -xyz_vector_ij
xyz_ji_unit = -xyz_ij_unit
xyz_ci_vector = xyz_vector_ji*np.stack([r_i/(r_i + r_j)]*3, axis=-1)
xyz_cj_vector = xyz_vector_ij*np.stack([r_j/(r_i + r_j)]*3, axis=-1)
overlapij = xyz_ij - (r_i + r_j)  # negative means contact
overlapij[overlapij >= -2*10**-9] = np.nan 
# check overlap from atom, distance from pair consistant
distance = df['distance'].values
overlapij_bypair = distance - (r_i + r_j)
distance[overlapij_bypair >= -2*10**-9] = np.nan 
overlapij_bypair[overlapij_bypair >= -2*10**-9] = np.nan 
check_distance_error = np.absolute((xyz_ij - distance)/overlapij)  #use overlab as denominator
#if (check_distance_error > error_tolerence).any():
#    sys.exit('distance not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_distance_error), np.nanmax(check_distance_error))))


v_vector_i = df[['vx_i', 'vy_i', 'vz_i']].values
vn_i = projection_length(v_vector_i, xyz_ij_unit)
vn_vector_i = projection_vector(v_vector_i, xyz_ij_unit)
vt_vector_i = verticle_vector(v_vector_i, xyz_ij_unit)

v_vector_j = df[['vx_j', 'vy_j', 'vz_j']].values
vn_j = projection_length(v_vector_j, xyz_ji_unit)  # j-->i direction positive
vn_vector_j = projection_vector(v_vector_j, xyz_ji_unit)
vt_vector_j = verticle_vector(v_vector_j, xyz_ji_unit)

v_vector_ij = v_vector_i - v_vector_j
vn_vector_ij = projection_vector(v_vector_ij, xyz_vector_ij)   
vt_vector_ij = verticle_vector(v_vector_ij, xyz_vector_ij)

#check v consistency 
vt_vector_ij_from_pair = df[['vijtx', 'vijty', 'vijtz']].values
vn_vector_ij_from_pair = df[['vijnx', 'vijny', 'vijnz']].values
v_vector_ij_from_pair = vt_vector_ij_from_pair + vn_vector_ij_from_pair
check_vij_error = np.absolute((v_vector_ij_from_pair - v_vector_ij)/v_vector_ij)
check_vij_error[v_vector_ij_from_pair < 10**-7] = np.nan

#if (check_vnij_error > error_tolerence).any():
#	sys.exit('vijn not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_vnij_error, axis=0), np.nanmax(check_vnij_error, axis=0))))

if (check_vij_error > error_tolerence).any():
	sys.exit('vij not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_vij_error, axis=0), np.nanmax(check_vij_error, axis=0))))
#
## check vijt from pair verticle to xyzij
#check_vijt_xyzij_costheta = np.absolute(projection_length(xyz_ij_unit, vt_vector_ij_from_pair))
#if (check_vijt_xyzij_costheta > error_tolerence).any():
#	sys.exit(
#			'vijt xyzij not orthogonal, max |costheta| = '
#			+repr((np.nanargmax(check_vijt_xyzij_costheta), np.nanmax(check_vijt_xyzij_costheta)))
#			)
## check vijn from pair parallel to xyzij
#check_vijn_xyzij_costheta = np.absolute(projection_length(xyz_ij_unit, vn_vector_ij_from_pair))
#if ((1 - check_vijn_xyzij_costheta) > error_tolerence).any():
#	sys.exit(
#			'vijn xyzij not parallel, (idrow, min|costheta|) = ' 
#			+repr((
#                np.nanargmin(np.absolute(check_vijn_xyzij_costheta)),
#                np.nanmin(np.absolute(check_vijn_xyzij_costheta))
#                ))
#			)

omega_vector_i = df[['omegax_i', 'omegay_i', 'omegaz_i']].values
omega_vector_j = df[['omegax_j', 'omegay_j', 'omegaz_j']].values

omegar_c_vector_i = np.cross(omega_vector_i, xyz_ci_vector)
omegar_c_vector_j = np.cross(omega_vector_j, xyz_cj_vector)

vc_vector_i = omegar_c_vector_i + v_vector_i
vc_vector_j = omegar_c_vector_j + v_vector_j
vct_vector_i = omegar_c_vector_i + vt_vector_i
vct_vector_j = omegar_c_vector_j + vt_vector_j

vct_vector_ij = vct_vector_i - vct_vector_j

fn_vector_ji = df[['fjinx', 'fjiny', 'fjinz']].values
ft_vector_ji = df[['fjitx', 'fjity', 'fjitz']].values
f_vector_ji = fn_vector_ji + ft_vector_ji

# f direction check
check_fn_xyzij_costheta = np.absolute(projection_length(xyz_ij_unit, fn_vector_ji))
check_ft_xyzij_costheta = np.absolute(projection_length(xyz_ij_unit, ft_vector_ji))
if ((1 - check_fn_xyzij_costheta) > error_tolerence).any():
	sys.exit(
			'fn xyzij not parallel, (idrow, min|costheta|) = ' 
			+repr((
                np.nanargmin(np.absolute(check_fn_xyzij_costheta)),
                np.nanmin(np.absolute(check_fn_xyzij_costheta))
                ))
			)
if (check_ft_xyzij_costheta > error_tolerence).any():
	sys.exit(
			'ft xyzij not orthogonal, max |costheta| = '
			+repr((np.nanargmax(check_ft_xyzij_costheta), np.nanmax(check_ft_xyzij_costheta)))
			)


# check fall = f
fall_i_vector = df[['f_force_all[1]_i', 'f_force_all[2]_i', 'f_force_all[3]_i']].values
f_vector_i = df[['fx_i', 'fy_i', 'fz_i']].values
f_vector_j = df[['fx_j', 'fy_j', 'fz_j']].values
a_vector_i = f_vector_i/np.stack([m_i]*3,axis=-1)
a_vector_j = f_vector_j/np.stack([m_j]*3,axis=-1)
check_fall_f_error = np.absolute((fall_i_vector - f_vector_i)/fall_i_vector)
if (check_fall_f_error > error_tolerence).any():
	sys.exit(
			'fall f not equal, (idrow, maxerror) = '
			+repr((np.nanargmax(check_fall_f_error, axis=0), np.nanmax(check_fall_f_error, axis=0)))
			)

# check fpair = sum(fji) 
fpair_vector = df.groupby(['step'])[['f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i']].agg('first').values
sumfn_vector_ji = df.groupby(['step'])[['fjinx', 'fjiny', 'fjinz']].agg('sum').values
sumft_vector_ji = df.groupby(['step'])[['fjitx', 'fjity', 'fjitz']].agg('sum').values
check_fpair_sumfji_error = (sumfn_vector_ji + sumft_vector_ji - fpair_vector)/fpair_vector
#if ~(check_fpair_sumfji_error <= error_tolerence).all():
#	sys.exit(
#			'fpair sumfji not equal, (idrow, maxerror) = '
#			+repr((np.nanargmax(check_fpair_sumfji_error, axis=0), np.nanmax(check_fpair_sumfji_error, axis=0)))
#			)

# check x v a consistant
xyzva_each_steps = df.groupby(['step'])[['type_i', 'x_i', 'y_i', 'z_i', 'vx_i', 'vy_i', 'vz_i', 'fx_i', 'fy_i', 'fz_i',]].agg('first')
xyz = xyzva_each_steps[['x_i', 'y_i', 'z_i',]].values
vxyz = xyzva_each_steps[['vx_i', 'vy_i', 'vz_i',]].values
fxyz = xyzva_each_steps[['fx_i', 'fy_i', 'fz_i',]].values
mass_each_steps = 4/3*np.pi*radius_type(xyzva_each_steps['type_i'].values)**3*2500
axyz = fxyz/np.stack([mass_each_steps]*3,axis=-1)
dxyz = np.diff(xyz, axis=0)
ddxyz = xyz[2:]+xyz[:-2]-2*xyz[1:-1]
ddxyz[ddxyz==0] = np.nan
dvxyz = np.diff(vxyz, axis=0)
axyz_average_two_steps = (axyz[0:-1] + axyz[1:])/2
# check dv = average (a1+a2) * timestep
check_dv_a = np.absolute((dvxyz - axyz_average_two_steps*ts)/dvxyz)
if (check_dv_a > error_tolerence).any():
    sys.exit('dv a not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_dv_a, axis=0), np.nanmax(check_dv_a, axis=0))))

# check ddx = a*ts**2
#check_ddx_a = (ddxyz - (axyz[1:-1]*ts**2))/ddxyz
#if (check_ddx_a > error_tolerence).any():
#    sys.exit('ddx a not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_ddx_a, axis=0), np.nanmax(check_ddx_a, axis=0))))

# check dx = v*t + 1/2*a*ts**2
#check_dx_va = (dxyz - (vxyz[0:-1]*ts + 0.5*axyz[0:-1]*ts**2))/dxyz
#if (check_dx_va > error_tolerence).any():
#    sys.exit('dx v a not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_dx_va, axis=0), np.nanmax(check_dx_va, axis=0))))

# check vij error 1 for atom 15583
dftrace = df.iloc[94027]
traceid_checkv = dftrace['id_j']
stepstrace = np.arange(dftrace['step']-1, dftrace['step']+3)
xyzijdifftrace = np.diff(xyz_vector_ij[(df['id_j']==traceid_checkv) & (df['step'].isin(stepstrace))],axis=0)
v_vector_ij_from_pairtrace = v_vector_ij_from_pair[(df['id_j']==traceid_checkv) & (df['step'].isin(stepstrace))]

v_vector_ijtrace = v_vector_ij[(df['id_j']==traceid_checkv) & (df['step'].isin(stepstrace))]
xyz_ij_unit_trace = xyz_ij_unit[(df['id_j']==traceid_checkv) & (df['step'].isin(stepstrace))]
print (df.loc[df['step'].isin(stepstrace)]['step'])
print (v_vector_ijtrace*ts)
print (v_vector_ij_from_pairtrace*ts)
print (xyzijdifftrace)

a_vector_ij_trace = (a_vector_i - a_vector_j)[(df['id_j']==traceid_checkv) & (df['step'].isin(stepstrace))]
check_dx_va_trace_vpair = (xyzijdifftrace - (v_vector_ij_from_pairtrace[0:-1]*ts + 0.5*a_vector_ij_trace[0:-1]*ts**2))/xyzijdifftrace
#if (check_dx_va_trace_vpair > error_tolerence).any():
#    sys.exit('dx v a not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_dx_va_trace_vpair, axis=0), np.nanmax(check_dx_va_trace_vpair, axis=0))))

check_dx_va_trace = (xyzijdifftrace - (v_vector_ijtrace[0:-1]*ts + 0.5*a_vector_ij_trace[0:-1]*ts**2))/xyzijdifftrace
if (check_dx_va_trace > error_tolerence).any():
    sys.exit('dx v a not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_dx_va_trace, axis=0), np.nanmax(check_dx_va_trace, axis=0))))
print(np.nanmax(check_dx_va_trace, axis=0))

# check 

LHS = v_vector_ij_from_pairtrace[1:] - v_vector_ijtrace[1:]

RHS = 0.5*ts*gamma_n*(v_vector_ij_from_pairtrace[:-1] + ts*a_vector_ij_trace[:-1] - v_vector_ijtrace[1:])




print(LHS)
print(RHS)
print(projection_vector(LHS, xyz_ij_unit_trace[:-1]))
print(projection_vector(RHS, xyz_ij_unit_trace[:-1]))

print ((v_vector_ijtrace-v_vector_ij_from_pairtrace))
print (a_vector_ij_trace*0.5*ts)
# torque ji
tq_vector_ji = np.cross(xyz_ci_vector, ft_vector_ji)
# wall force zplane inwall outwall
# torque wall i
# check torque  'tqx_i', 'tqy_i', 'tqz_i', = sum tq_ji tq_wall_i

# timer
end = time.process_time()
print(end-start)


#%%
df.iloc[94027,:]['id_j']
tracej = 15583
tracenearstep = 76383
pickrow = (df['id_j']==tracej) & (df['step'].isin(list(range(tracenearstep-3, tracenearstep+4))))
print(np.where(pickrow))
print('v_vector_i')
print(v_vector_i[pickrow])
print('v_vector_j')
print(v_vector_j[pickrow])
print('v_vector_ij')
print(v_vector_ij[pickrow])
print('v_vector_ij_frompair')
print(v_vector_ij_from_pair[pickrow])
print('v_vector_ij_frompair - v_vector_ij')
print(v_vector_ij_from_pair[pickrow] - v_vector_ij[pickrow])
print('(a_vector_i-a_vector_j)*0.5ts')
print((a_vector_i[pickrow]- a_vector_j[pickrow])*0.5*ts)

#%%
ran = range(233249,233256)
#dtest = df.iloc[233249:233256,:]
#print(dtest[['step','id_j','vx_i','vy_i','vz_i','vx_j','vy_j','vz_j','vijtx', 'vijty', 'vijtz','vijnx', 'vijny', 'vijnz']])
print (vn_vector_i[ran])
an_vector_i = projection_vector(a_vector_i, xyz_ij_unit)
print (an_vector_i[ran])

#%%
print (np.diff(v_vector_i[ran],axis=0))
print (10**-6*a_vector_i[ran])
print(projection_vector(10**-6*a_vector_i, xyz_vector_ij)[ran])
print(projection_vector(np.diff(v_vector_i,axis=0), xyz_vector_ij[1:])[ran])
#%%
# check fpair = sum(fji) 
fpair_vector = df.groupby(['step'])[['f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i']].agg('first').values
sumfn_vector_ji = df.groupby(['step'])[['fjinx', 'fjiny', 'fjinz']].agg('sum').values
sumft_vector_ji = df.groupby(['step'])[['fjitx', 'fjity', 'fjitz']].agg('sum').values
check_fpair_sumfji_error = (sumfn_vector_ji + sumft_vector_ji - fpair_vector)/fpair_vector
if ~(check_fpair_sumfji_error <= error_tolerence).all():
	sys.exit(
			'fpair sumfji not equal, (idrow, maxerror) = '
			+repr((np.nanargmax(check_fpair_sumfji_error, axis=0), np.nanmax(check_fpair_sumfji_error, axis=0)))
			)
#%%
# check vpair = 

v_particular_idj = df.loc[df['id_j'] == 15583]

vi_particular_idj = v_particular_idj[['vx_i', 'vy_i', 'vz_i']].values

vj_particular_idj = v_particular_idj[['vx_j', 'vy_j', 'vz_j']].values
vijn_particular_idj = v_particular_idj[['vijnx', 'vijny', 'vijnz']].values
vijt_particular_idj = v_particular_idj[['vijtx', 'vijty', 'vijtz']].values


check_v15583_error = (vi_particular_idj - vj_particular_idj - (vijn_particular_idj + vijt_particular_idj))/(vijn_particular_idj + vijt_particular_idj)
print(v_particular_idj.iloc[34098:34101,:])
print(vi_particular_idj[34098:34101])
print(vj_particular_idj[34098:34101])
print(vijn_particular_idj[34098:34101])
print(vijt_particular_idj[34098:34101])
print(vijn_particular_idj[34098:34101] + vijt_particular_idj[34098:34101])
print(vi_particular_idj[34098:34101] - vj_particular_idj[34098:34101])
if ~(check_v15583_error <= error_tolerence).all():
	sys.exit(
			'fpair sumfji not equal, (idrow, maxerror) = '
			+repr((np.nanargmax(check_v15583_error, axis=0), np.nanmax(check_v15583_error, axis=0)))
			)  
#%%
print(fpair_vector[60559:60563])
print(sumfn_vector_ji[60559:60563])
print(sumft_vector_ji[60559:60563])
print(sumft_vector_ji[60559:60563]+sumfn_vector_ji[60559:60563])
#%%
dtestsingle = dtest.iloc[0::2,:]
dtestsinglevi = dtestsingle[['vx_i','vy_i','vz_i']].values
dtestsinglevj = dtestsingle[['vx_j','vy_j','vz_j']].values
dtestsinglemi =  mass_type(dtestsingle['type_i'].values)
dtestsinglemj =  mass_type(dtestsingle['type_j'].values)
dtestsingleai = dtestsingle[['fx_i','fy_i','fz_i']].values/np.stack([dtestsinglemi]*3, axis=-1)
dtestsingleaj = dtestsingle[['fx_j','fy_j','fz_j']].values/np.stack([dtestsinglemj]*3, axis=-1)
dtestsingleaij = dtestsingleai - dtestsingleaj
dtestsinglevijc = dtestsinglevi - dtestsinglevj
dtestsinglevijcp = dtestsinglevijc + dtestsingleaij*0.5*10**-6
dtestsinglevijcm = dtestsinglevijc - dtestsingleaij*0.5*10**-6
dtestsinglevij = dtestsingle[['vijnx','vijny','vijnz']].values + dtestsingle[['vijtx','vijty','vijtz']].values
print (dtestsingle[['step','id_j']])
print (dtestsinglevijc)
print (dtestsinglevijcp)
print (dtestsinglevijcm)
print (dtestsinglevij)
print (dtestsinglevijcp - dtestsinglevijcm)
print (np.diff(dtestsinglevij,axis=0))

dtestsingle = dtest.iloc[1::2,:]
dtestsinglevi = dtestsingle[['vx_i','vy_i','vz_i']].values
dtestsinglevj = dtestsingle[['vx_j','vy_j','vz_j']].values
dtestsinglemi =  mass_type(dtestsingle['type_i'].values)
dtestsinglemj =  mass_type(dtestsingle['type_j'].values)
dtestsingleai = dtestsingle[['fx_i','fy_i','fz_i']].values/np.stack([dtestsinglemi]*3, axis=-1)
dtestsingleaj = dtestsingle[['fx_j','fy_j','fz_j']].values/np.stack([dtestsinglemj]*3, axis=-1)
dtestsingleaij = dtestsingleai - dtestsingleaj
dtestsinglevijc = dtestsinglevi - dtestsinglevj
dtestsinglevijcp = dtestsinglevijc + dtestsingleaij*0.5*10**-6
dtestsinglevijcm = dtestsinglevijc - dtestsingleaij*0.5*10**-6
dtestsinglevij = dtestsingle[['vijnx','vijny','vijnz']].values + dtestsingle[['vijtx','vijty','vijtz']].values
print (dtestsingle[['step','id_j']])
print (dtestsinglevijc)
print (dtestsinglevijcp)
print (dtestsinglevijcm)
print (dtestsinglevij)
print (dtestsinglevijcp - dtestsinglevijcm)
print (np.diff(dtestsinglevij,axis=0))

#%%
# check fnjipair vijpair overlapijatom consistant 
fn_vector_ji_calculate_bypair = -kn*np.stack([overlapij]*3, axis=-1)*xyz_ij_unit - gamma_n*np.stack([meff]*3, axis=-1)*vn_vector_ij_from_pair
check_fji_bypair = np.absolute((fn_vector_ji_calculate_bypair - fn_vector_ji)/fn_vector_ji)
# if (check_fji_bypair > error_tolerence).any():
#    sys.exit('fjipair calculate not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_fji_bypair, axis=0), np.nanmax(check_fji_bypair, axis=0))))

# check verlet with damp
vij



# check if only one contact , in n direction , fpair=fjiatom from vijatom overlapijatom consistant
df_1contact = df.loc[df['c_n_contact_i'] == 1]
fpair_vector_1c = df_1contact.groupby(['step'])[['f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i']].agg('first').values
xyz_vector_ij_1c = xyz_vector_ij[df['c_n_contact_i'] == 1]
fpair_vector_1c_n = projection_vector(fpair_vector_1c, xyz_vector_ij_1c)
fpair_vector_n_calculate = -kn*np.stack([overlapij]*3, axis=-1)*xyz_ij_unit - gamma_n*np.stack([meff]*3, axis=-1)*vn_vector_ij
fpair_vector_1c_n_calculate = fpair_vector_n_calculate[df['c_n_contact_i'] == 1]
check_fjiatom_1c = np.absolute((fpair_vector_1c_n_calculate - fpair_vector_1c_n)/fpair_vector_1c_n)
if (check_fjiatom_1c > error_tolerence).any():
	sys.exit('fjipair calculate not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_fjiatom_1c, axis=0), np.nanmax(check_fjiatom_1c, axis=0))))

#%%
df.loc[df['c_n_contact_i'] == 1]['step'].values[281]
df.loc[df['step'].isin(np.arange(43089,43092))][['x_i','y_i','z_i', 'fx_i', 'fy_i', 'fz_i','fjinx','fjitx','fjiny','fjity','fjinz','fjitz', 'f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i','f_force_all[3]_i']]
#fpair_vector_1c43090 = 
#%%
(0.025565**2+0.017682**2)**0.5/0.00083-37
#%%   
# conclusion not to use overlap from pair
# now only use pair
checkrow1 = 29985
print(checkrow1)
print(check_fji_bypair[checkrow1])
print(check_fji_byatom[checkrow1])
print(check_fji_atom_pair[checkrow1])
print(fn_vector_ji[checkrow1])
print(fn_vector_ji_calculate_bypair[checkrow1])
print(fn_vector_ji_calculate_byatom[checkrow1])
print((-kn*np.stack([overlapij_bypair]*3, axis=-1)*xyz_ij_unit)[checkrow1])
print((-kn*np.stack([overlapij]*3, axis=-1)*xyz_ij_unit)[checkrow1])
print((- gamma_n*np.stack([meff]*3, axis=-1)*vn_vector_ij_from_pair)[checkrow1])
print((- gamma_n*np.stack([meff]*3, axis=-1)*vn_vector_ij)[checkrow1])

checkrow1 = 195989
print(checkrow1)
print(check_fji_bypair[checkrow1])
print(check_fji_byatom[checkrow1])
print(check_fji_atom_pair[checkrow1])
print(fn_vector_ji[checkrow1])
print(fn_vector_ji_calculate_bypair[checkrow1])
print(fn_vector_ji_calculate_byatom[checkrow1])
print((-kn*np.stack([overlapij_bypair]*3, axis=-1)*xyz_ij_unit)[checkrow1])
print((-kn*np.stack([overlapij]*3, axis=-1)*xyz_ij_unit)[checkrow1])
print((- gamma_n*np.stack([meff]*3, axis=-1)*vn_vector_ij_from_pair)[checkrow1])
print((- gamma_n*np.stack([meff]*3, axis=-1)*vn_vector_ij)[checkrow1])


#%%
print(fn_vector_ji[checkrow1])
print(overlapij_bypair[checkrow1])
print(overlapij[checkrow1])
print(vn_vector_ij_from_pair[checkrow1])
print(vn_vector_ij[checkrow1])
#%%
print(dxyz[84182])
print(vxyz[84182])
print((0.5*axyz[0:-1]*ts**2)[84182])
print(check_dx_va[84182])
#%%
print(fpair_vector[44509])
print(sumfn_vector_ji[44509])
print(sumft_vector_ji[44509])
#%%
aa = df.loc[(df['step'].isin([76381,76382,76383,76384])) & (df['id_j']==15583)]
vv = aa.values
insert = (vv[0:-1] + vv[1:])/2
dinsert = pd.DataFrame(data = insert, columns = list(aa), dtype = 'float64')
fillhalf = aa.append(dinsert).sort_values(by=['step'])
massi = 4/3*np.pi*radius_type(aa['type_i'].values)**3*2500
axi = aa['fx_i'].values/massi
massj = 4/3*np.pi*radius_type(aa['type_j'].values)**3*2500
axj = aa['fx_j'].values/massj
vxihalfp = aa['vx_i'].values + 0.5*10**-6*axi
vxjhalfp = aa['vx_j'].values + 0.5*10**-6*axj
vxihalfm = aa['vx_i'].values - 0.5*10**-6*axi
vxjhalfm = aa['vx_j'].values - 0.5*10**-6*axj
print(aa['vx_i'].values - aa['vx_j'].values)
print(aa['vijnx'].values + aa['vijtx'].values)
print(vxihalfp)
print(vxjhalfp)
print(vxihalfm)
print(vxjhalfm)
C1 = vxihalfp-vxjhalfp
C3 = vxihalfm-vxjhalfm
C2 = aa['vijnx'].values + aa['vijtx'].values
print(C1)
print(C2)
print(C1[0:-1]-C2[1:])
print(C1-C2)
print(C3)
print(C2)
print(C3[0:-1]-C2[1:])
print(C3-C2)

#%%
aa['c_n_contact_i']

#%%
insert2 = (vv[0:-2] + vv[2:])/2
dinsert2 = pd.DataFrame(data = insert2, columns = list(df9402594030), dtype = 'float64')
insert2ij = dinsert2[['vx_i','vy_i','vz_i']].values - dinsert2[['vx_j','vy_j','vz_j',]].values
insert2ijpair = dinsert2[['vijnx', 'vijny', 'vijnz']].values + dinsert2[['vijtx', 'vijty', 'vijtz']].values
print(fillhalfij[2]-insert2ijpair)
print(fillhalfijpair[2]-insert2ij)
#%%
list(fillhalf)

#%%
aa['fx_i']

#%%
xijdiff = np.diff(aa['x_i'].values-aa['x_j'].values)
print(xijdiff)

vijtttt = (aa['vijnx'].values + aa['vijtx'].values)
print(xijdiff/vijtttt[0:3])
print(xijdiff/vijtttt[1:4])
print(xijdiff/(vijtttt[1:4]+vijtttt[0:3])*2)
#%%
xidiff = np.diff(aa['x_i'].values)
print(xidiff)

vitttt = aa['vx_i'].values 
print(xidiff/vitttt[0:3])
print(xidiff/vitttt[1:4])
print(xidiff/(vitttt[1:4]+vitttt[0:3])*2)

#%%
aa[['x_i','y_i','z_i','x_j', 'y_j', 'z_j']]
aa[['f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i']]
aa[['f_force_all[1]_i', 'f_force_all[2]_i', 'f_force_all[3]_i']]
aa[['fx_i', 'fy_i', 'fz_i']]
bb = df.loc[df['step'].isin([76381,76382,76383,76384])]
bb[['fjinx','fjitx']]
#%%
ri = np.sqrt(np.sum(aa[['x_i','y_i']].values**2, axis=-1))/0.00083 - 37
rj = np.sqrt(np.sum(aa[['x_j','y_j']].values**2, axis=-1))/0.00083 - 37
print(ri)
print(rj)

#%%
check_distance_error

#%%
np.nanmax(check_distance_error)

#%%
v15556 = np.array([1.686196085711881*10**-6  ,-8.852552269706013*10**-6 ,-1.826084896779219*10**-5])
f15556 = np.array([-1.626872276392055*10**-7 ,-1.548920962408988*10**-8 ,-2.703051493712974*10**-7])
v15583 = np.array([1.685666610907108*10**-6  ,-2.88589931340842*10**-5  ,-1.915450565985885*10**-5]) 
f15583 = np.array([-5.947568645198572*10**-7 ,-6.283281379894822*10**-7 ,-1.922792346965301*10**-8])
vn1558315556 = np.array([-8.2148*10**-6      ,-4.56281*10**-6 			,	-2.71612*10**-6 ])
vt1558315556 = np.array([8.06519*10**-6      ,-1.54578*10**-5 			,	1.57477*10**-6])
vnerror = projection_vector((v15583-v15556),vn1558315556)-vn1558315556
print (vnerror)
vterror = projection_vector((v15583-v15556),vt1558315556)-vt1558315556
print (vterror)
print (v15583-v15556-projection_vector((v15583-v15556),vn1558315556)-projection_vector((v15583-v15556),vt1558315556))
#%%
v_h_15556 = v15556 - f15556/mass_type(1)*ts/2
#%%
(0.02560872364422996**2 + 0.01761750612522943**2)**0.5/0.00083


#%%
