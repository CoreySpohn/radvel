# cimport the Cython declarations for numpy
from __future__ import absolute_import
cimport numpy as cnp
import numpy as np
import cython
from libc.math cimport sin
from libc.math cimport cos


# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
cnp.import_array()

# Wrapping kepler(M,e) a simple function that takes two doubles as
# arguments and returns a double
cdef extern from "kepler.h":
    double kepler(double M, double e)
    double rv_drive(double t, double per, double tp, double e, double cosom, double sinom, double k )
    void eccanom_orvara(double E[], double sinE[], double cosE[], const double M[], const double e, const int n);

DTYPE = np.float64
ctypedef cnp.float64_t DTYPE_t

# create the wrapper code, with numpy type annotations
@cython.boundscheck(False)
def kepler_array(double [:,] M, double e):
    cdef int size, i

    size = M.shape[0]
    cdef cnp.ndarray[double, ndim=1] E = \
        cnp.ndarray(shape=(size,), dtype=DTYPE) 

    for i in range(size):
        E[i] = kepler(M[i], e)

    return E 

# create the wrapper code, with numpy type annotations
@cython.boundscheck(False)
def rv_drive_array(cnp.ndarray[DTYPE_t, ndim=1] t, double per, double tp, 
                   double e, double om, double k):
    cdef int size, i 
    size = t.shape[0]

    cdef cnp.ndarray[DTYPE_t, ndim=1] rv = t.copy()
    cdef double cosom = cos(om)
    cdef double sinom = sin(om)
    for i in range(size):
        rv[i] = rv_drive(t[i], per, tp, e, cosom, sinom, k)

    return rv
