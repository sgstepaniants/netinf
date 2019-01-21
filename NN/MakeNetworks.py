import numpy as np
from scipy.sparse import diags

# create connectivity networks

def tridiag(n, corners=True):
  diagonals = [np.ones(n-1), np.ones(n-1)]
  out = diags(diagonals, [-1, 1]).toarray()
  
  if (corners):
    out[0, n-1] = 1
    out[n-1, 0] = 1
  
  return out


def erdos_reyni(n, p, zero_diag=True):
  mat = np.random.rand(n, n)
  out = np.zeros((n, n))

  # if >= 1-p, add connection
  out[mat >= 1 - p] = 1
  if zero_diag:
      np.fill_diagonal(out, 0)
  return out


def symmetric_erdos_reyni(n, p, zero_diag=True):
  out = erdos_reyni(n, p, zero_diag)
  i_lower = np.tril_indices(n, -1)
  out[i_lower] = out.T[i_lower]
  return out
