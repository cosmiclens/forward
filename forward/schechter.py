import numpy as np
from collections import namedtuple

SchechterParameters = namedtuple('SchechterParameters', ['a_phi', 'b_phi', 'a_m', 'b_m', 'alpha'])

def dv_domega_dz(z, cosmology):
    d_h = cosmology.hubble_distance
    d_m = cosmology.comoving_transverse_distance(z)
    e_fac = np.sqrt(cosmology.inv_efunc(z))
    return d_h * d_m * d_m * e_fac

def schechter(m, z, parameters):
    phi_star = parameters.b_phi * np.exp(parameters.a_phi*z)
    m_star = parameters.a_m * z + parameters.b_m
    out = 0.4 * np.log(10) * phi_star
    out = out * np.power(10, 0.4*(m_star-m)*(parameters.alpha+1))
    out = out * np.exp(-np.power(10, 0.4*(m_star-m)))
    return out
