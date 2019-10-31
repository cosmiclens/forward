import numpy as np

def integrate_sed(z, λe, fe, λx, Rx, extinction=None):
  if extinction is not None:
    raise NotImplementedError('cannot handle extinction')

  if len(np.shape(z)) > 1:
    raise ValueError('redshift z must be 1d array or scalar')

  if len(np.shape(λx)) != 1 or len(np.shape(Rx)) != 1:
      raise ValueError('filter λx, Rx must be 1d arrays')

  λo = np.expand_dims(λe, -1)*np.add(z, 1)
  fo = np.expand_dims(fe, -1)
  Ro = np.apply_along_axis(lambda λ: np.interp(λ, λx, Rx), -1, λo)

  return np.trapz(λo*fo*Ro, λo, axis=-2)

