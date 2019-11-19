import numpy as np

# This code is based on https://arxiv.org/abs/astro-ph/0401052v1
# It computes the density laws given in table 3 in units of stars/pc^3

def disc_density(a = 0. , age = 0.):
    ''' This function computes the disc density
        Argunents: a, age
        a ? distance in pc
        age in Gyears
        rho_0 = 1.03E-3 stars/pc^3 local density
        d0 = ??? normalisation constant
        Constants:
        age_cutoff = 0.15 Gyr
        hr_plus = 5000. pc (if age <= 0.15 Gyr) 2530. pc (otherwise) scale length
        hr_minus = 3000. pc (if age <= 0.15 Gyr) 1320. pc (otherwise) hole scale length
    '''
    rho_0 = 1.03E-3
    d0 = ???
    age_cutoff = 0.15
    if age <= age_cutoff:
        hr_plus = 5000.0
        hr_minus = 3000.0
        discd = (rho_0/d0)*np.exp(- (a/hr_plus)**2 )  - (rho_0/d0)*np.exp(- (a/hr_minus)**2 )
    elif age > age_cutoff:
        hr_plus = 2530.0
        hr_minus = 1320.0
        discd = (rho_0/d0)*np.exp(- np.sqrt(0.25 + (a/hr_plus)**2 ))  - (rho_0/d0)*np.exp(- np.sqrt(0.25 + (a/hr_minus)**2 ))
    else:
        print('Age must be a positive number')
    return discd

def thickdisc_density(R,z):
    ''' This function computes the thick disc density
        Arguments: R, z
        R radius in parsec
        z euclidian height
        Constants:
        xl = 400 pc some length??
        hz = 800. pc scale length
        hr = 450. pc scale height
        rho_0 = 1.03E-3 stars/pc^3 local density
        d0 = ??? normalisation constant
        Rc = 2540 pc cutoff radius
    '''
    xl = 400.
    hz = 800.
    hr = 2500.
    rho_0 = 1.03E-3
    d0 = ???
    Rc = 2540
    if np.absolute(z) <= xl:
        thickd = (rho_0/d0) * np.exp((Rc - R)/hr) * (1. - z**2 / ((hz*xl)(2. +xl))  )
    else:
        thickd = rho_0 * np.exp((Rc - R)/hr) * np.exp( (xl - np.absolute(z))/ hz )/(1. + xl/(2. * hz))
    return thickd


def spheroid_density(a = 0., rho_0 = 1.03E-3,d0, Rc = 2540):
    ''' This function computes the spheroidal density
        Arguments: a, ac, rho_0, d0, Rc
            a   ? distance in pc
            ac = 500 pc cutoff distance
            rho_0 = 1.03E-3 stars/pc^3 local density
            d0 = ??? normalisation constant
            Rc = 2540 pc cutoff radius
    '''
    ac = 500.
    if a <= ac:
        sphed = (rho_0/d0) * (a/Rc)**(-2.44)
    else:
        sphed = rho_0 * (a/Rc)**(-2.44)
    return sphed


def bulge_density(x = 0.1, y = 0.1, z = 0.1, x0 = 1.59, y0 = 0.424, z0 = 0.424, N = 13.70, Rc = 2.54):
    ''' This function computes the bulge density
        Arguments: x,y,z, x0, y0, z0, N, Rc

            x,y,z euclidian dimensions in kpc
            x0 = 1.59 kpc scale length
            y0 = 0.424 kpc mayor axis
            z0 = 0.424 kpc minor axis
            N = 13.70 stars/pc^3 normalisation or star density at the center
            Rc = 2.54 kpc cutoff radius
    '''
    r = np.sqrt(x**2 + y**2)
    rn2 = (x/x0)**2 + (y/y0)**2
    rs2 = np.sqrt( rn2**2 + (z/z0)**4)
    if r < Rc:
        bulged = N * np.exp(- 0.5 * rs2)
    elif r > Rc:
        bulged = N * np.exp(- 0.5 * rs2 - 2 * ( r - Rc )**2)
    else:
        print('The bulge dimensions x and y cannot be null')

    return bulged

def ISM_density(R = 0., z = 0., rho_0 = 1.03E-3, Rmw = 8000., hz = 140., hr = 450.):
    ''' This function computes the ISM density
        Arguments: R, z, rho_0, Rmw, hz, hr

            R is the Galactocentric distance
            z is height above the Galactic plane
            rho_0 = 1.03E-3 stars/pc^3 local density
            Rmw = 8000. pc ? Galactocentric distance of the Milky Way
            hz = 140. pc scale length
            hr = 450. pc scale height
    '''
    ismd = rho_0 * np.exp( - ( R - Rmw )/ hr - np.absolute( z )/hz  )
    return ismd

def dmhalo_density(R = 0., z = 0. , epsilon = 1., rho_c =0.1079, Rc = 2697. ):
    ''' This function computes the dark halo density
        Arguments: R, z, epsilon, rho_c, Rc

            R is the Galactocentric distance
            z is height above the Galactic plane
            epsilon is the axis ratio
            rho_c =0.1079 critical density in units of ??
            Rc = 2697. Galactocentric? distance in pc

    '''
    a2 = R**2 + (z/epsilon)**2
    halodm = rho_c / ( 1. + a2 / R**2   )
    return halodm
