#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataSciece.changeDirOnImportExport setting
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'inputvariable'))
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
# define radius function
def radius_by_type(type_atom):
	r0 =0.000415
	rp = r0*(
		0.9*(type_atom == 1)
		+1.0*(type_atom == 2)
		+1.1*(type_atom == 3)
		)
	return rp
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
dp = 0.00083
rp = dp/2
# intersection point
z_zplane = 0*dp
r_in = 37*dp
r_out = (37+16)*dp

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
	return n*np.stack([projection_length(v, n)]*3,axis=-1)
def verticle_vector(v, n):
	return v - projection_vector(v, n)

# r_p_wall_vector
def r_p_wall_vector(position):
	return position - wallpoint(position)
# calculate radius from type
def radius_type(type):
	return rp*(0.9*(type == 1) + 1.0*(type == 2) + 1.1*(type == 3))
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

m_i = 4/3*np.pi*r_i**3*2500
m_j = 4/3*np.pi*r_j**3*2500

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

# check distance from pair consistant
distance = df['distance'].values
check_distance_error = np.absolute((xyz_ij - distance)/xyz_ij)
if (check_distance_error > error_tolerence).any():
    sys.exit('distance not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_distance_error), np.nanmax(check_distance_error))))


v_vector_i = df[['vx_i', 'vy_i', 'vz_i']].values
vn_i = projection_length(v_vector_i, xyz_ij_unit)
vn_vector_i = projection_vector(v_vector_i, xyz_ij_unit)
vt_vector_i = verticle_vector(v_vector_i, xyz_ij_unit)

v_vector_j = df[['vx_j', 'vy_j', 'vz_j']].values
vn_j = projection_length(v_vector_j, xyz_ji_unit)  # j-->i direction positive
vn_vector_j = projection_vector(v_vector_j, xyz_ji_unit)
vt_vector_j = verticle_vector(v_vector_j, xyz_ji_unit)

v_vector_ij = v_vector_i - v_vector_j
vn_vector_ij = vn_vector_i - vn_vector_j     
vt_vector_ij = vt_vector_i - vt_vector_j

#check v consistency 
vt_vector_ij_from_pair = df[['vijtx', 'vijty', 'vijtz']].values
vn_vector_ij_from_pair = df[['vijnx', 'vijny', 'vijnz']].values
check_vij_error = np.absolute((vn_vector_ij_from_pair + vt_vector_ij_from_pair - v_vector_ij)/v_vector_ij)
if (check_vij_error > error_tolerence).any():
	sys.exit('vij not consistent, (idrow, maxerror) = ' + repr((np.nanargmax(check_vij_error, axis=0), np.nanmax(check_vij_error, axis=0))))
# check vijt from pair verticle to xyzij
check_vijt_xyzij_costheta = np.absolute(projection_length(xyz_ij_unit, vt_vector_ij_from_pair))
if (check_vijt_xyzij_costheta > error_tolerence).any():
	sys.exit(
			'vijt xyzij not orthogonal, max |costheta| = '
			+repr((np.nanargmax(check_vijt_xyzij_costheta), np.nanmax(check_vijt_xyzij_costheta)))
			)
# check vijn from pair parallel to xyzij
check_vijn_xyzij_costheta = np.absolute(projection_length(xyz_ij_unit, vn_vector_ij_from_pair))
if ((1 - check_vijn_xyzij_costheta) > error_tolerence).any():
	sys.exit(
			'vijn xyzij not parallel, (idrow, min|costheta|) = ' 
			+repr((
                np.nanargmin(np.absolute(check_vijn_xyzij_costheta)),
                np.nanmin(np.absolute(check_vijn_xyzij_costheta))
                ))
			)

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
f_i_vector = df[['fx_i', 'fy_i', 'fz_i']].values
check_fall_f_error = np.absolute((fall_i_vector - f_i_vector)/fall_i_vector)
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
if ~(check_fpair_sumfji_error <= error_tolerence).all():
	sys.exit(
			'fpair sumfji not equal, (idrow, maxerror) = '
			+repr((np.nanargmax(check_fpair_sumfji_error, axis=0), np.nanmax(check_fpair_sumfji_error, axis=0)))
			)


# torque ji
tq_vector_ji = np.cross(xyz_ci_vector, ft_vector_ji)
# wall force zplane inwall outwall
# torque wall i
# check torque  'tqx_i', 'tqy_i', 'tqz_i', = sum tq_ji tq_wall_i

# timer
end = time.process_time()
print(end-start)


#%%
# check fpair = sum(fji) 
fpair_vector = df.groupby(['step'])[['f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i']].agg('first').values
sumfn_vector_ji = df.groupby(['step'])[['fjinx', 'fjiny', 'fjinz']].agg('sum').values
sumft_vector_ji = df.groupby(['step'])[['fjitx', 'fjity', 'fjitz']].agg('sum').values
check_fpair_sumfji_error = (sumfn_vector_ji + sumft_vector_ji - fpair_vector)/fpair_vector
#%%
a=sumfn_vector_ji[44509]
#%%
b=sumft_vector_ji[44509]
#%%
c=fpair_vector[44509]

#%%
a+b-c
#%%
(a+b-c)/c
#%%
e = df.loc[df['step']==44509]
e1 = df.loc[df['step']==44508]
#%%
e[['id_j','f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i',
'fjinx', 'fjiny', 'fjinz',
'fjitx', 'fjity', 'fjitz' ]]
#%%
e1[['id_j','f_force_pair[1]_i', 'f_force_pair[2]_i', 'f_force_pair[3]_i',
'fjinx', 'fjiny', 'fjinz',
'fjitx', 'fjity', 'fjitz' ]]
#%%
e['fjiny'].values+e['fjity'].values

#%%
e['fjity'].values

#%%
e1['f_force_pair[2]_i']
#%%
vvv = df.loc[107630,['step','id_i', 'id_j', 'vx_i', 'vy_i', 'vz_i', 'vx_j', 'vy_j', 'vz_j', 'vijtx', 'vijty', 'vijtz']]
vvv
#%%
list(df)
#%%

#check v consistency 
vt_vector_ij_from_pair = df[['vijtx', 'vijty', 'vijtz']].values

#%%
vt_vector_i[107630]

#%%
vt_vector_j[107630]

#%%
vt_vector_i[107630] - vt_vector_j[107630]

#%%
vt_vector_ij[107630]
#%%
#check v consistency 
vt_vector_ij_from_pair [107630]
#%%
(vt_vector_ij_from_pair[107630] - vt_vector_ij[107630])
#%%
a=np.array([[1,2],[3,4]])
b=np.array([1,2])
a*b/b

#%%
#type(xyz_ci_vector) 
type(xyz_vector_ji)
#type(r_i)
#type(r_i + r_j)

#%%
check_vijt_error

#%%
np.where(np.isnan(check_vijt))[0]

#%%
check_vijt[10000]

#%%
np.where(np.isnan(check_vijt))

#%%
X = np.array([1,2,3,4])
N = 7
np.vstack([X]*N)

#%%
N*[X]

#%%
np.hstack((1,1))

#%%
np.nanmax(check_vijt_error)

#%%
np.nan!=1

#%%
a=np.array([1,2,3])
b=-a
a[0]=10
b

#%%
