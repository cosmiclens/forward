import numpy as np

def integrate_sed(z, le, fe, lx, Rx, extinction=None):
  '''Integrate redshifted SED against filter band.

  Arguments:
  z (number or array-like): redshift values for observed SED
  le (array-like): emitted wavelength λ_e of SED
  fe (array-like): emitted function values f_e(λ_e) of SED
  lx (array-like): filter wavelength λ_x
  Rx (array-like): filter response values R_x(λ_x)

  Returns:
  array of shape (#redshifts, #SEDs, #bands)
  '''

  if extinction is not None:
    raise NotImplementedError('cannot handle extinction')

  if np.ndim(z) > 1:
    raise ValueError('redshift z must be 1d array or scalar')

  # make sure Rx is a 2d array
  # then use broadcasting to make lx same shape
  Rx = np.atleast_2d(Rx)
  lx = np.ones(Rx.shape)*np.atleast_2d(lx)

  lo = np.expand_dims(le, -1)*np.add(z, 1)
  fo = np.expand_dims(fe, -1)
  Ro = [[np.interp(lo, lx1, Rx1, left=0, right=0)] for lx1, Rx1 in zip(lx, Rx)]

  result = np.trapz(lo*fo*Ro, lo, axis=-2).transpose()

  return result

