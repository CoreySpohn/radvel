# cimport the Cython declarations for numpy
from __future__ import absolute_import
cimport numpy as np
import numpy as np
import cython
from libc.math cimport sin
from libc.math cimport cos

# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# Wrapping kepler(M,e) a simple function that takes two doubles as
# arguments and returns a double
cdef extern from "rv.h":
    void rv_calc(double rv[], double t[], double params[], char basis[], int ntimes, int nplan, int specific_planet)
    double timetrans_to_timeperi_c(double tc, double per, double ecc, double omega)
    void model_call(double mod[], const double t[], const double params[],
                    const char basis[], const double time_base, const int jit_ind,
                    const int gamma_ind, const int ntimes, const int nplan,
                    const int specific_planet)
    void residuals(double res[], const double t[], const double rv[], 
                     const double rv_err[], const double params[], const char basis[], 
                     const double time_base, const int jit_ind, const int gamma_ind, 
                     const int ntimes, const int nplan)
    double logprob(const double t[], const double rv[], const double rv_err[],
                   const double params[], const char basis[],
                   const double time_base, const int jit_ind, const int gamma_ind,
                   const int ntimes, const int nplan)

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


# @cython.boundscheck(False)
# def calc_rv_from_params_old(np.ndarray[DTYPE_t, ndim=1] t, np.ndarray[DTYPE_t, ndim=2] params, char basis[], int nplan):

#     cdef int ntimes = len(t)
#     cdef np.ndarray[DTYPE_t, ndim=1] rv = np.zeros(ntimes)
#     rv_calc(<double*> rv.data, <double*> t.data, <double*> params.data, <char*> basis, ntimes, nplan)

#     return rv

@cython.boundscheck(False)
def calc_rv_forward_model(np.ndarray[DTYPE_t, ndim=1] t, param_obj, vector, planet_num=None):

    cdef np.ndarray[DTYPE_t, ndim=2] params = vector.vector
    basis_bytes = param_obj.basis.name.encode("utf-8")
    cdef char* basis = basis_bytes
    cdef int nplan = param_obj.num_planets
    cdef int ntimes = len(t)
    cdef int specific_planet = 0
    cdef np.ndarray[DTYPE_t, ndim=1] rv = np.zeros(ntimes)

    rv_calc(<double*> rv.data, <double*> t.data, <double*> params.data, <char*> basis, ntimes, nplan, specific_planet)

    return rv

@cython.boundscheck(False)
def timetrans_to_timeperi(double tc, double per, double ecc, double omega):

    tp = timetrans_to_timeperi_c(tc, per, ecc, omega)
    return tp

@cython.boundscheck(False)
def model_call_c(np.ndarray[DTYPE_t, ndim=1] t, param_obj, vector, double time_base, int jit_index, int gamma_index, *args, **kwargs):

    cdef np.ndarray[DTYPE_t, ndim=2] params = vector.vector
    basis = param_obj.basis.name.encode("utf-8")
    cdef int nplan = param_obj.num_planets
    cdef int ntimes = len(t)
    cdef np.ndarray[DTYPE_t, ndim=1] model_values = np.zeros(ntimes)
    if "planet_num" in kwargs.keys():
        pnum = kwargs['planet_num']
    else:
        pnum = 0
    # const int gamma_ind = vector.indices[self.gamma_param]
    # const int jit_ind = vector.indices[self.jit_param]

    model_call(<double*> model_values.data, <double*> t.data,
               <double*> params.data, <char*> basis, <double> time_base, 
               <int> jit_index, <int> gamma_index, <int> ntimes, <int> nplan, 
               <int> pnum)

    return model_values

@cython.boundscheck(False)
def residuals_c(np.ndarray[DTYPE_t, ndim=1] t, np.ndarray[DTYPE_t, ndim=1] rv, np.ndarray[DTYPE_t, ndim=1] rv_err, 
        param_obj, vector, double time_base, int jit_index, int gamma_index):

    cdef np.ndarray[DTYPE_t, ndim=2] params = vector.vector
    # basis = param_obj.basis.name.encode("utf-8")
    basis_bytes = param_obj.basis.name.encode("utf-8")
    cdef char* basis = basis_bytes
    cdef int nplan = param_obj.num_planets
    cdef int ntimes = len(t)
    cdef int specific_planet = 0
    cdef np.ndarray[DTYPE_t, ndim=1] res = np.zeros(ntimes)
    cdef np.ndarray[DTYPE_t, ndim=1] rv_err2 = np.zeros(ntimes)
    cdef np.ndarray[DTYPE_t, ndim=1] sum_sig_quad = np.zeros(ntimes)

    residuals(<double*> res.data,<double*> <double*> t.data, <double*> rv.data, <double*> rv_err.data,
                       <double*> params.data, <char*> basis, <double> time_base, 
                       <int> jit_index, <int> gamma_index, <int> ntimes, <int> nplan)

    return res

@cython.boundscheck(False)
def logprob_c(np.ndarray[DTYPE_t, ndim=1] t, np.ndarray[DTYPE_t, ndim=1] rv, np.ndarray[DTYPE_t, ndim=1] rv_err, 
        param_obj, vector, double time_base, int jit_index, int gamma_index):

    cdef np.ndarray[DTYPE_t, ndim=2] params = vector.vector
    # basis = param_obj.basis.name.encode("utf-8")
    basis_bytes = param_obj.basis.name.encode("utf-8")
    cdef char* basis = basis_bytes
    cdef int nplan = param_obj.num_planets
    cdef int ntimes = len(t)

    _logprob = logprob(<double*> t.data, <double*> rv.data, <double*> rv_err.data,
                       <double*> params.data, <char*> basis, <double> time_base, 
                       <int> jit_index, <int> gamma_index, <int> ntimes, <int> nplan)

    return _logprob
