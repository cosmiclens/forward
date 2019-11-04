import numpy as np
import itertools as it


def abflux(z, le, fe, lx, Rx, extinction=None):
  '''AB flux for a single SED and filter at multiple redshifts.

  Does not include the cosmological scaling with luminosity distance.
  '''

  le = np.expand_dims(le, -1)
  fe = np.expand_dims(fe, -1)
  lo = le*(1 + np.atleast_1d(z))
  Ro = np.interp(lo, lx, Rx, left=0, right=0)

  result = np.trapz(lo*fe*Ro, le, axis=-2)/np.trapz(lo*Ro, lo, axis=-2)

  return result


def abmag(z, le, fe, lx, Rx, extinction=None):
  '''Compute magnitude at source redshifts from SEDs and filters.

  Does not include the cosmological scaling with (10pc/d_L)^2.

  Arguments:
  z (number or array): redshift values for observed SED
  le (array or list of arrays): emitted wavelength λ_e of SED
  fe (array or list of arrays): emitted flux density f_e(λ_e) of SED
  lx (array or list of arrays): filter wavelength λ_x
  Rx (array or list of arrays): filter response values R_x(λ_x) in counts/photon

  Returns:
  array of shape (#redshifts, #SEDs, #bands) with singular dimensions squeezed
  '''

  if extinction is not None:
    raise NotImplementedError('cannot handle extinction')

  # check input dimensions
  if np.ndim(z) > 1:
    raise ValueError('redshift z must be 1d array or scalar')

  if np.ndim(le) < 1 or np.ndim(fe) < 1 or \
      (np.ndim(le[0]) > 0 and np.ndim(le[0][0]) > 0) or \
      (np.ndim(fe[0]) > 0 and np.ndim(fe[0][0]) > 0):
    raise ValueError('SED le, fe must be 1d array or list of 1d arrays')

  if np.ndim(lx) < 1 or np.ndim(Rx) < 1 or \
      (np.ndim(lx[0]) > 0 and np.ndim(lx[0][0]) > 0) or \
      (np.ndim(Rx[0]) > 0 and np.ndim(Rx[0][0]) > 0):
    raise ValueError('filter lx, Rx must be 1d array or list of 1d arrays')

  # expand SEDs and filters to be 2d of the same length
  if np.ndim(le[0]) == 0:
    le = it.repeat(le)
  if np.ndim(fe[0]) == 0:
    fe = np.expand_dims(fe, 0)
  if np.ndim(lx[0]) == 0:
    lx = it.repeat(lx)
  if np.ndim(Rx[0]) == 0:
    Rx = np.expand_dims(Rx, 0)

  # compute AB fluxes for all bands and SEDs
  fluxes = [
    [
      abflux(z, le1, fe1, lx1, Rx1, extinction)
      for le1, fe1 in zip(le, fe)
    ]
    for lx1, Rx1 in zip(lx, Rx)
  ]

  # transpose to get shape (len(z), len(fe), len(Rx)) for easier processing
  fluxes = np.transpose(fluxes)

  # convert fluxes to magnitudes
  mag = -2.5*np.log10(fluxes)

  return mag

