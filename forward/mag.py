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
  integral of (redshifted SED)*(filter pass) over observed wavelength
  '''

  if extinction is not None:
    raise NotImplementedError('cannot handle extinction')

  if np.ndim(z) > 1:
    raise ValueError('redshift z must be 1d array or scalar')

  if np.ndim(lx) != 1 or np.ndim(Rx) != 1:
      raise ValueError('filter lx, Rx must be 1d arrays')

  lo = np.expand_dims(le, -1)*np.add(z, 1)
  fo = np.expand_dims(fe, -1)
  Ro = np.interp(lo, lx, Rx, left=0., right=0.)

  result = np.trapz(lo*fo*Ro, lo, axis=-2)

  return result

