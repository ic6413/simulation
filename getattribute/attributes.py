print('create attribute')

startstep = 0
ts = 4e-7
# atom radius
dp0 = 0.00083
density = 2500
r_in = 37*dp0
r_out = (37+16)*dp0
g = [0, 0, -9.8]
mu = 0.5
kn = 8837.32631448687 
kt = 2524.95037556768 
gamma_n = 34445.5603308471
gamma_t = 17222.7801654236
type_radius_list = [
    [1, 0.9*dp0/2],
    [2, 1.0*dp0/2],
    [3, 1.1*dp0/2],
]

walls_p = [
    [
        'p',
        [0,0,0],
        [0,0,1],
    ],
]
walls_cy = [
    [
        'cy',
        [0,0,0],
        [0,0,1],
        r_in,
    ],
    [
        'cy',
        [0,0,0],
        [0,0,1],
        r_out,
    ]
]

print ("finish creating attribute")