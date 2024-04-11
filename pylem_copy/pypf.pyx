import numpy as np
import copy
cimport pypf
cimport numpy as np
from libc.stdlib cimport malloc, free

def flood(np.ndarray[double, ndim = 2, mode = 'c'] dem not None):

  m, n = dem.shape[0], dem.shape[1]

  pypfc(&dem[0,0], m, n);
