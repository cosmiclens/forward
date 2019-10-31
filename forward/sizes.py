import numpy as np

# This code generates the half-light physical radius and the half-light physical angular radius
# Based on https://arxiv.org/pdf/1705.05386.pdf

# The physical radius of a galaxy lognormally distributed 
# The function's arguments: magnitude M, mean parameters amu, bmu, and standard deviation sigma
# The mean is given by equation [3.14]
# Free parameters: amu, bmu and sigma

def rphys50(M, amu = 1, bmu = 0, sigma = 1, size = 100):
    muphys = amu * M + bmu
    r = np.random.lognormal(muphys, sigma, size)
    return r

