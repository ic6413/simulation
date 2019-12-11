##### calculate hertz contact model parameter

#### import package
import numpy as np

####
def k_n(young_modules, poisson_ratio):
    return 2*young_modules/3/(1-poisson_ratio**2)

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

