import numpy as np

def integrate_sed(z, le, fe, lx, Rx, extinction=None):
  if extinction is not None:
    raise NotImplementedError('cannot handle extinction')

  if len(np.shape(z)) > 1:
    raise ValueError('redshift z must be 1d array or scalar')

  if len(np.shape(lx)) != 1 or len(np.shape(Rx)) != 1:
      raise ValueError('filter lx, Rx must be 1d arrays')

  lo = np.expand_dims(le, -1)*np.add(z, 1)
  fo = np.expand_dims(fe, -1)
  Ro = np.apply_along_axis(lambda l: np.interp(l, lx, Rx), -1, lo)

  result = np.trapz(lo*fo*Ro, lo, axis=-2)

  return result

