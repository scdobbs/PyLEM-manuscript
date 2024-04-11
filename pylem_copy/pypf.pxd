# pypf.pxd

cdef extern from "pypfc.h" nogil:
    ctypedef signed int int32_t;
    void pypfc(double *dem, int32_t m, int32_t n)
