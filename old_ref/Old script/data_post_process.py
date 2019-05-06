#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataSciece.changeDirOnImportExport setting
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'input'))
	print(os.getcwd())
except:
	pass

#%%
# import
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
import data_post_process_LAMMPS as dp
import LAMMPS_attributes as La
# ======================================
## setting
# select id, number or 'all'
select_id_list = [15555, 15556]
begin = 0
end = 300000

f_read_select_cipcj = (
	'/home/ic6413/mylammps/myfiles/run/'
	+ '0103_mu05_nofreeze_traceid15556/'
	+ 'output/python/'
	+ 'data_select_cipcj'
	+ '_'
	+ ''.join(str(e) for e in select_id_list)
	+ '.h5'
)

df = pd.read_hdf(f_read_select_cipcj, 'df')

# error test
error_tolerence = 10**-3

## read LAMMPS attribute
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
# zplane1 zplane2....
zplane_list = La.zplane_list
# zcylinder1 zcylinder2....
zcylinder_list = La.zcylinder_list
# timestep
ts = La.ts

##======== process data from df =============================================================
xyz_vector_i = df[['x_i', 'y_i', 'z_i']].values
xyz_vector_j = df[['x_j', 'y_j', 'z_j']].values
xyz_vector_ij = xyz_vector_i - xyz_vector_j

r_i = dp.radius_by_type(rp0, df[['type_i']].values, type_radius_list)
r_j = dp.radius_by_type(rp0, df[['type_j']].values, type_radius_list)
m_i = dp.mass_by_type(rp0, df[['type_i']].values, type_radius_list, density)
m_j = dp.mass_by_type(rp0, df[['type_j']].values, type_radius_list, density)
meff = m_i*m_j/(m_i + m_j)

xyz_ij = dp.length(xyz_vector_ij)
xyz_ij_unit = dp.unit(xyz_vector_ij)
xyz_vector_ji = -xyz_vector_ij
xyz_ji_unit = -xyz_ij_unit
xyz_ci_vector = xyz_vector_ji*np.stack([r_i/(r_i + r_j)]*3, axis=-1)
xyz_cj_vector = xyz_vector_ij*np.stack([r_j/(r_i + r_j)]*3, axis=-1)
overlapij = xyz_ij - (r_i + r_j)  # negative means contact
overlapij[overlapij >= -2*10**-9] = np.nan 
# check the "distance" in Pair file is consistant with the overlap calculate from Atom files
distance = df[['distance']].values
overlapij_bypair = distance - (r_i + r_j)
distance[overlapij_bypair >= -2*10**-9] = np.nan 
overlapij_bypair[overlapij_bypair >= -2*10**-9] = np.nan 
check_distance_error = np.absolute((xyz_ij - distance)/overlapij)  #use overlab as denominator
#if (check_distance_error > error_tolerence).any():
#    sys.exit('distance not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_distance_error), np.nanmax(check_distance_error))))


v_vector_i = df[['vx_i', 'vy_i', 'vz_i']].values
vn_i = dp.projection_scalar(v_vector_i, xyz_ij_unit)
vn_vector_i = dp.projection_vector(v_vector_i, xyz_ij_unit)
vt_vector_i = dp.verticle_vector(v_vector_i, xyz_ij_unit)

v_vector_j = df[['vx_j', 'vy_j', 'vz_j']].values
vn_j = dp.projection_scalar(v_vector_j, xyz_ji_unit)  # j-->i direction positive
vn_vector_j = dp.projection_vector(v_vector_j, xyz_ji_unit)
vt_vector_j = dp.verticle_vector(v_vector_j, xyz_ji_unit)

v_vector_ij = v_vector_i - v_vector_j
vn_vector_ij = dp.projection_vector(v_vector_ij, xyz_vector_ij)   
vt_vector_ij = dp.verticle_vector(v_vector_ij, xyz_vector_ij)

# check vn consistency 
vt_vector_ij_from_pair = df[['vijtx', 'vijty', 'vijtz']].values
vn_vector_ij_from_pair = df[['vijnx', 'vijny', 'vijnz']].values
v_vector_ij_from_pair = vt_vector_ij_from_pair + vn_vector_ij_from_pair
check_vnij_error = np.absolute((vn_vector_ij_from_pair - vn_vector_ij)/vn_vector_ij)
#check_vnij_error[vn_vector_ij_from_pair < 10**-7] = np.nan

#if (check_vnij_error > error_tolerence).any():
#	sys.exit('vijn not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_vnij_error, axis=0), np.nanmax(check_vnij_error, axis=0))))

# debug function
def debug_print_vn_incon(
						row1,
						row2,
						):
	for row in range(row1,row2):
		print(' ')
		print('(row, error) = ' + repr((row, check_vnij_error[row])))
		print('costheta Pair')
		print(dp.costheta(vn_vector_ij_from_pair[row], v_vector_ij_from_pair[row]))
		print('costheta Atom')
		print(dp.costheta(vn_vector_ij[row], v_vector_ij[row]))
		print('vnij Pair')
		print(vn_vector_ij_from_pair[row])
		print('vnij Atom calcu')
		print(vn_vector_ij[row])
		print('vni Atom calcu')
		print(vn_vector_i[row])
		print('vnj Atom calcu')
		print(vn_vector_j[row])
		print('distant error')
		print(check_distance_error[row])
		print('overlapij')
		print(overlapij[row])
		print('overlap by pair')
		print(overlapij_bypair[row])

# debug find the first row not consistent
for row, check_vnij_error_one_row in enumerate(check_vnij_error):
	if (check_vnij_error_one_row > error_tolerence).any():
		debug_print_vn_incon(
						row-1,
						row+2,
						)
		sys.exit('vij not consistent, (first row, error) = ' + repr((row, check_vnij_error_one_row)))
	
# debug 
if (check_vnij_error > error_tolerence).any():
	for row in np.nanargmax(check_vnij_error, axis=0):
		print('row = ' + str(row))
		print('costheta')
		print(dp.length(vn_vector_ij_from_pair[row])/dp.length(v_vector_ij_from_pair[row]))
		print('vnij')
		print(vn_vector_ij_from_pair[row])
		print('vnij Atom calcu')
		print(vn_vector_ij[row])
		print('vni Atom calcu')
		print(vn_vector_i[row])
		print('vnj Atom calcu')
		print(vn_vector_j[row])
		row = row -1
		print('row = ' + str(row))
		print('costheta')
		print(dp.length(vn_vector_ij_from_pair[row])/dp.length(v_vector_ij_from_pair[row]))
		print('vnij')
		print(vn_vector_ij_from_pair[row])
		print('vnij Atom calcu')
		print(vn_vector_ij[row])
		print('vni Atom calcu')
		print(vn_vector_i[row])
		print('vnj Atom calcu')
		print(vn_vector_j[row])

	sys.exit('vij not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_vnij_error, axis=0), np.nanmax(check_vnij_error, axis=0))))

if (check_vnij_error > error_tolerence).any():
	for row in np.nanargmax(check_vnij_error, axis=0):
		print('row = ' + str(row))
		print('costheta')
		print(dp.length(vn_vector_ij_from_pair[row])/dp.length(v_vector_ij_from_pair[row]))
		print('vnij')
		print(vn_vector_ij_from_pair[row])
		print('vnij Atom calcu')
		print(vn_vector_ij[row])
		print('vni Atom calcu')
		print(vn_vector_i[row])
		print('vnj Atom calcu')
		print(vn_vector_j[row])
		row = row -1
		print('row = ' + str(row))
		print('costheta')
		print(dp.length(vn_vector_ij_from_pair[row])/dp.length(v_vector_ij_from_pair[row]))
		print('vnij')
		print(vn_vector_ij_from_pair[row])
		print('vnij Atom calcu')
		print(vn_vector_ij[row])
		print('vni Atom calcu')
		print(vn_vector_i[row])
		print('vnj Atom calcu')
		print(vn_vector_j[row])

	sys.exit('vij not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_vnij_error, axis=0), np.nanmax(check_vnij_error, axis=0))))
#
## check vijt from pair verticle to xyzij
#check_vijt_xyzij_costheta = np.absolute(projection_scalar(xyz_ij_unit, vt_vector_ij_from_pair))
#if (check_vijt_xyzij_costheta > error_tolerence).any():
#	sys.exit(
#			'vijt xyzij not orthogonal, max |costheta| = '
#			+repr((np.nanargmax(check_vijt_xyzij_costheta), np.nanmax(check_vijt_xyzij_costheta)))
#			)
## check vijn from pair parallel to xyzij
#check_vijn_xyzij_costheta = np.absolute(projection_scalar(xyz_ij_unit, vn_vector_ij_from_pair))
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
check_fn_xyzij_costheta = np.absolute(dp.projection_scalar(xyz_ij_unit, fn_vector_ji))
check_ft_xyzij_costheta = np.absolute(dp.projection_scalar(xyz_ij_unit, ft_vector_ji))
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
mass_each_steps = np.vectorize(dp.mass_by_type)(rp0, xyzva_each_steps['type_i'], type_radius_list, density)
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
print(dp.projection_vector(LHS, xyz_ij_unit_trace[:-1]))
print(dp.projection_vector(RHS, xyz_ij_unit_trace[:-1]))

print ((v_vector_ijtrace-v_vector_ij_from_pairtrace))
print (a_vector_ij_trace*0.5*ts)
# torque ji
tq_vector_ji = np.cross(xyz_ci_vector, ft_vector_ji)
# wall force zplane inwall outwall
# torque wall i
# check torque  'tqx_i', 'tqy_i', 'tqz_i', = sum tq_ji tq_wall_i

#%%
len(r_i)

len(df['type_i'].values)
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
an_vector_i = dp.projection_vector(a_vector_i, xyz_ij_unit)
print (an_vector_i[ran])

#%%
print (np.diff(v_vector_i[ran],axis=0))
print (10**-6*a_vector_i[ran])
print(dp.projection_vector(10**-6*a_vector_i, xyz_vector_ij)[ran])
print(dp.projection_vector(np.diff(v_vector_i,axis=0), xyz_vector_ij[1:])[ran])
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
dtestsinglemi =  dp.mass_type(dtestsingle['type_i'].values)
dtestsinglemj =  dp.mass_type(dtestsingle['type_j'].values)
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


# check if only one contact , in n direction , fpair=fjiatom from vijatom overlapijatom consistant
df_1contact = df.loc[df['c_n_contact_i'] == 1]
fpair_vector_1c = df_1contact.groupby(['step'])[['f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i']].agg('first').values
xyz_vector_ij_1c = xyz_vector_ij[df['c_n_contact_i'] == 1]
fpair_vector_1c_n = dp.projection_vector(fpair_vector_1c, xyz_vector_ij_1c)
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
massi = 4/3*np.pi*dp.radius_by_type(r0=La.dp0/2, type=aa['type_i'].values, type_radius_list=La.type_radius_list)**3*2500
axi = aa['fx_i'].values/massi
massj = 4/3*np.pi*dp.radius_by_type(r0=La.dp0/2, type=aa['type_j'].values, type_radius_list=La.type_radius_list)**3*2500
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
vnerror = ds.projection_vector((v15583-v15556),vn1558315556)-vn1558315556
print (vnerror)
vterror = ds.projection_vector((v15583-v15556),vt1558315556)-vt1558315556
print (vterror)
print (v15583-v15556-ds.projection_vector((v15583-v15556),vn1558315556)-ds.projection_vector((v15583-v15556),vt1558315556))
#%%
v_h_15556 = v15556 - f15556/mass_type(1)*ts/2
#%%
(0.02560872364422996**2 + 0.01761750612522943**2)**0.5/0.00083




#%%
