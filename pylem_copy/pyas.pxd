# pyas.pxd

cdef extern from "pyasc.h" nogil:
  ctypedef signed int int32_t;
  void pyasc(double *dem, double dx, double *a, double *s, int32_t m, int32_t n);
  void pyasc_dinf(double *dem, double dx, double *a, double *s, int32_t m, int32_t n);
  void pylc(double *dem, double dx, double *l, int32_t m, int32_t n);
