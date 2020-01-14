##### calculate hertz contact model parameter

#### import package
import numpy as np

####

def k_n(young_modules, poisson_ratio):
    return 2*young_modules/3/(1-poisson_ratio**2)

def k_0(density, d):
    r = d/2
    m = 4/3*np.pi*r**3*density
    return m*9.8/d

def gamma_n(young_modules, poisson_ratio, meff, shear_visocus, bulk_viscous):
    A = (
        1/3
        *(3*bulk_viscous-shear_visocus)**2
        /(3*bulk_viscous+2*shear_visocus)
        *(1-poisson_ratio**2)
        *(1-2*poisson_ratio)
        /(young_modules*poisson_ratio**2)
        )
    return 3/2*A/meff*k_n(young_modules, poisson_ratio)

def k_t(young_modules, poisson_ratio):
    return 2*(1-poisson_ratio)/(2-poisson_ratio)*k_n(young_modules, poisson_ratio)

def gamma_t(young_modules, poisson_ratio, meff, shear_visocus, bulk_viscous):
    A = (
        1/3
        *(3*bulk_viscous-shear_visocus)**2
        /(3*bulk_viscous+2*shear_visocus)
        *(1-poisson_ratio**2)
        *(1-2*poisson_ratio)
        /(young_modules*poisson_ratio**2)
        )
    return 3/2*A/meff*k_t(young_modules, poisson_ratio)

def collision_time(young_modules, poisson_ratio, meff, initial_v_n):
    return 2.94*(meff/k_n(young_modules, poisson_ratio))**(2/5)*(initial_v_n)**(-1/5)

def max_velocity_allowed(young_modules, poisson_ratio, meff, min_number_step_in_a_collision, timestep):
    min_collison_t = timestep*min_number_step_in_a_collision
    max_velocity_allowed = (
                            min_collison_t
                            /(
                            2.94*(meff/k_n(young_modules, poisson_ratio))**(2/5)
                            )
                            )**(-5)
    return max_velocity_allowed

def estitution_coefficients(young_modules, poisson_ratio, meff, ):


def print_all_contact_parameter(young_modules=8.36*10**10, poisson_ratio=0.22,
                                diameter=10**-4,density=2500,
                                shear_visocus=10**-15, bulk_viscous=0,
                                initial_v_n=0.001):
    r = diameter/2
    density = 2500
    m = 4/3*np.pi*r**3*density
    meff = m/2

    print(
        "kn "
        + '{:.3E}'.format(k_n(young_modules, poisson_ratio))
    )
    
    print(
        "gamma_n "
        + '{:.3E}'.format(gamma_n(young_modules, poisson_ratio, meff, shear_visocus, bulk_viscous))
    )

    print(
        "kt "
        + '{:.3E}'.format(k_t(young_modules, poisson_ratio))
    )

    print(
        "gamma_t "
        + '{:.3E}'.format(gamma_t(young_modules, poisson_ratio, meff, shear_visocus, bulk_viscous))
    )

    print(
        "collision_time "
        + '{:.3E}'.format(collision_time(young_modules, poisson_ratio, meff, initial_v_n))
    )

    print(
        "m "
        + '{:.3E}'.format(m)
    )