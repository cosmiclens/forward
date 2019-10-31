import numpy as np

def integrate_sed(z, le, fe, lx, Rx, extinction=None):
  '''Integrate one or many redshifted SEDs against filter bands.

  Arguments:
  z (float or array-like): redshift values for observed SEDs
  le (array or list of arrays): emitted wavelength λ_e of SEDs
  fe (array or list of arrays): emitted function values f_e(λ_e) of SEDs
  lx (array of float): filter wavelength λ_x
  Rx (array of float): filter response values R_x(λ_x)

  Returns:
  array with shape (number of SEDs, number of redshifts)
  '''

  if extinction is not None:
    raise NotImplementedError('cannot handle extinction')

  if len(np.shape(z)) > 1:
    raise ValueError('redshift z must be 1d array or scalar')

  if len(np.shape(lx)) != 1 or len(np.shape(Rx)) != 1:
      raise ValueError('filter lx, Rx must be 1d arrays')

  lo = np.expand_dims(le, -1)*np.add(z, 1)
  fo = np.expand_dims(fe, -1)
  Ro = np.apply_along_axis(lambda l: np.interp(l, lx, Rx, left=0., right=0.), -1, lo)

  result = np.trapz(lo*fo*Ro, lo, axis=-2)

  return result

