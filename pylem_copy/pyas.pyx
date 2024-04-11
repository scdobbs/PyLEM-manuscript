import numpy as np
import copy
cimport pyas
cimport numpy as np
from libc.stdlib cimport malloc, free

def area_dinf(np.ndarray[double, ndim = 2, mode = 'c'] dem not None, float dx):

  m, n = dem.shape[0], dem.shape[1]
  cdef np.ndarray[double, ndim = 2, mode = 'c'] a = np.zeros((m,n), dtype = float)
  cdef np.ndarray[double, ndim = 2, mode = 'c'] s = np.zeros((m,n), dtype = float)

  pyasc_dinf(&dem[0,0], dx, &a[0,0], &s[0,0], m, n)

  return a, s

def area(np.ndarray[double, ndim = 2, mode = 'c'] dem not None, float dx):

  m, n = dem.shape[0], dem.shape[1]
  cdef np.ndarray[double, ndim = 2, mode = 'c'] a = np.zeros((m,n), dtype = float)
  cdef np.ndarray[double, ndim = 2, mode = 'c'] s = np.zeros((m,n), dtype = float)

  pyasc(&dem[0,0], dx, &a[0,0], &s[0,0], m, n)

  return a, s

def length(np.ndarray[double, ndim = 2, mode = 'c'] dem not None, float dx):
  m, n = dem.shape[0], dem.shape[1]
  cdef np.ndarray[double, ndim = 2, mode = 'c'] l = np.zeros((m,n), dtype = float)

  pylc(&dem[0,0], dx, &l[0,0], m, n)

  return l
