# import python module
import csv
import json
import os.path

# import module in simulation folder
import osmanage as om
import datapath as dp

dict = {}

dict['ts'] = 4e-7
# atom radius
dict['dp0'] = 0.00083
dict['density'] = 2500
dict['r_in'] = 37*dict['dp0']
dict['r_out'] = (37+16)*dict['dp0']
dict['g'] = [0, 0, -9.8]
dict['mu'] = 0.5
dict['kn'] = 8837.32631448687 
dict['kt'] = 2524.95037556768 
dict['gamma_n'] = 34445.5603308471
dict['gamma_t'] = 17222.7801654236
dict['type_radius_list'] = [
    [1, 0.9*dict['dp0']/2],
    [2, 1.0*dict['dp0']/2],
    [3, 1.1*dict['dp0']/2],
]

dict['walls_p'] = [
    [
        'p',
        [0,0,0],
        [0,0,1],
    ],
]
dict['walls_cy'] = [
    [
        'cy',
        [0,0,0],
        [0,0,1],
        dict['r_in'],
    ],
    [
        'cy',
        [0,0,0],
        [0,0,1],
        dict['r_out'],
    ]
]

if not os.path.exists(dp.f_attribute):

    om.create_directory(dp.post_process_path)
    om.create_directory(dp.attribute_json_path)

    with open(dp.f_attribute, 'w') as f:
        json.dump(dict, f, indent=4)

else:
    print ("attribute file already exist")