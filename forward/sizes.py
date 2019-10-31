import numpy as np

# This code generates the half-light physical radius and the half-light physical angular radius
# Based on https://arxiv.org/pdf/1705.05386.pdf

# The physical radius of a galaxy lognormally distributed 
# The function's arguments: magnitude M, mean parameters amu, bmu, and standard deviation sigma
# The mean is given by equation [3.14]
# Free parameters: amu, bmu and sigma

def rphys50(M, amu = 1, bmu = 0, sigma = 1, size = None):
    muphys = amu * M + bmu
    r = np.random.lognormal(muphys, sigma, size)
    return r

# The next function transforms the half-light physical radius into the angular radius eq [3.15]
# Arguments: physical radius r, redshift z, cosmology c as given by astropy.cosmology
# The angular diameter distance dA is computed by using from e.g. astropy.cosmology import WMAP9 as cosmo

def rang50(r = 1, z = 0, c):
    dA = c.angular_diameter_distance(z)
    ang = r / dA
    return ang
