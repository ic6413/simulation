import csv
import json

dict = {}

dict['ts'] = 10**-6
# atom radius
dict['dp0'] = 0.00083
dict['density'] = 2500
dict['r_in'] = 37*dict['dp0']
dict['r_out'] = (37+16)*dict['dp0']
dict['mu'] = 0.5
dict['kn'] = 8837.32631448687 
dict['kt'] = 2524.95037556768 
dict['gamma_n'] = 34445.5603308471
dict['gamma_t'] = 17222.7801654236
dict['type_radius_list'] = [
    [1, 0.9],
    [2, 1.0],
    [3, 1.1],
]
dict['zplane_list'] = [0]*dict['dp0']
dict['zcylinder_list'] = [dict['r_in'] , dict['r_out']]

with open('attribute.json', 'w') as fp:
	json.dump(dict, fp)
